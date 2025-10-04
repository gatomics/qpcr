
"""
Gatomis VCF + Phenotype Analyzer (RUO)
Core utilities:
- VCF parsing (VCF/VCF.GZ), supports SnpEff ANN and VEP CSQ
- Lightweight phenotype scoring (HPO → genes), panel filtering
"""

from typing import Dict, List, Tuple, Optional
import io, gzip, re, json
import pandas as pd

# ------------ VCF parsing ------------

def _open_text(file_obj):
    if isinstance(file_obj, (str, bytes)):
        if str(file_obj).endswith(".gz"):
            return io.TextIOWrapper(gzip.open(file_obj, "rb"))
        else:
            return open(file_obj, "r", encoding="utf-8", errors="ignore")
    head = file_obj.read(2)
    file_obj.seek(0)
    if head == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.GzipFile(fileobj=file_obj, mode="rb"))
    else:
        return io.TextIOWrapper(file_obj, encoding="utf-8", errors="ignore")

def parse_vcf(file_obj, max_variants: int = 500000):
    fh = _open_text(file_obj)
    csq_fields = None
    ann_fields = ["Allele","Annotation","Impact","Gene_Name","Gene_ID","Feature_Type","Feature_ID",
                  "Transcript_BioType","Rank/Total","HGVS.c","HGVS.p","cDNA.pos/cDNA.length",
                  "CDS.pos/CDS.length","AA.pos/AA.length","Distance","ERRORS/WARNINGS/INFO"]
    samples = []
    records = []
    stats = {"n_total":0,"n_pass":0,"snps":0,"indels":0,"by_chrom":{},"filters":{},"ti":0,"tv":0}

    for line in fh:
        if line.startswith("##"):
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format:\s*([^\">]+)", line)
                if m:
                    csq_fields = [f.strip() for f in m.group(1).split("|")]
            continue
        if line.startswith("#CHROM"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 8:
                samples = parts[9:]
            break

    for line in fh:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
        fmt = parts[8] if len(parts) > 8 else None
        sample_vals = parts[9:] if len(parts) > 9 else []
        gt = None
        ad = None
        dp = None
        if fmt and sample_vals:
            fmt_keys = fmt.split(":")
            sv = sample_vals[0].split(":")
            fm = dict(zip(fmt_keys, sv))
            gt = fm.get("GT", None)
            dp = _to_float(fm.get("DP", None))
            ad = fm.get("AD", None)

        stats["n_total"] += 1
        stats["filters"][flt] = stats["filters"].get(flt, 0) + 1
        if flt in ("PASS","."):
            stats["n_pass"] += 1
        alt_alleles = alt.split(",")
        for a in alt_alleles:
            if len(ref)==1 and len(a)==1:
                stats["snps"] += 1
                if (ref, a) in [("A","G"),("G","A"),("C","T"),("T","C")]:
                    stats["ti"] += 1
                else:
                    stats["tv"] += 1
            else:
                stats["indels"] += 1
        stats["by_chrom"][chrom]=stats["by_chrom"].get(chrom,0)+1

        info_dict = {}
        for kv in info.split(";"):
            if not kv: 
                continue
            if "=" in kv:
                k,v = kv.split("=",1)
                info_dict[k]=v
            else:
                info_dict[kv]=True

        af = first_float(info_dict, ["AF","AF_POPMAX","gnomAD_AF","VAF"])
        mq = first_float(info_dict, ["MQ"])

        gene = None; conseq=None; impact=None; hgvsc=None; hgvsp=None
        if "ANN" in info_dict:
            ann_entries = info_dict["ANN"].split(",")
            rank = {"HIGH":3,"MODERATE":2,"LOW":1,"MODIFIER":0}
            best=None; best_rank=-1
            for e in ann_entries:
                fields = e.split("|")
                fields += [""]*(len(ann_fields)-len(fields))
                ann, imp, sym, hgc, hgp = fields[1], fields[2], fields[3], fields[9], fields[10]
                r = rank.get(imp,0)
                if r>best_rank:
                    best_rank=r; best=(sym, ann, imp, hgc, hgp)
            if best:
                gene, conseq, impact, hgvsc, hgvsp = best
        elif "CSQ" in info_dict:
            csq_entries = info_dict["CSQ"].split(",")
            fields = csq_fields or []
            rank = {"HIGH":3,"MODERATE":2,"LOW":1,"MODIFIER":0}
            best=None; best_rank=-1
            for e in csq_entries:
                vals = e.split("|")
                d = { fields[i] if i<len(fields) else f"F{i}": (vals[i] if i<len(vals) else "") for i in range(max(len(vals), len(fields))) }
                imp = d.get("IMPACT","")
                sym = d.get("SYMBOL","")
                cons = d.get("Consequence","") or d.get("CONSEQUENCE","")
                hgc = d.get("HGVSc","")
                hgp = d.get("HGVSp","")
                r = rank.get(imp,0)
                if r>best_rank:
                    best_rank=r; best=(sym, cons, imp, hgc, hgp)
            if best:
                gene, conseq, impact, hgvsc, hgvsp = best

        rec = {
            "CHROM": chrom, "POS": int(pos),
            "REF": ref, "ALT": alt,
            "QUAL": _to_float(qual), "FILTER": flt,
            "AF": af, "DP": dp, "MQ": mq, "AD": ad,
            "GENE": gene, "CONSEQUENCE": conseq, "IMPACT": impact,
            "HGVSc": hgvsc, "HGVSp": hgvsp, "GT": gt
        }
        records.append(rec)
        if len(records) >= max_variants:
            break

    return pd.DataFrame(records), {"samples": samples, "stats": stats}

def _to_float(x):
    try:
        return float(x)
    except Exception:
        return None

def first_float(info_dict, keys):
    for k in keys:
        if k in info_dict:
            try:
                return float(str(info_dict[k]).split(",")[0])
            except:
                pass
    return None

# ------------ Phenotype utilities ------------

def load_hpo_map(df_or_file) -> pd.DataFrame:
    """
    Accepts CSV/TSV with columns: HPO_ID, GeneSymbol
    Returns DataFrame with normalized columns.
    """
    if isinstance(df_or_file, pd.DataFrame):
        df = df_or_file.copy()
    else:
        if str(df_or_file).lower().endswith(".tsv"):
            df = pd.read_csv(df_or_file, sep="\t")
        else:
            df = pd.read_csv(df_or_file)
    cols = {c.lower(): c for c in df.columns}
    # normalize expected names
    hpo_col = next((cols[k] for k in cols if k in ("hpo_id","hpo","term_id")), None)
    gene_col = next((cols[k] for k in cols if k in ("genesymbol","gene","symbol")), None)
    if not hpo_col or not gene_col:
        raise ValueError("HPO map must contain columns: HPO_ID, GeneSymbol (or equivalents).")
    out = df.rename(columns={hpo_col:"HPO_ID", gene_col:"GeneSymbol"})[["HPO_ID","GeneSymbol"]].dropna()
    out["GeneSymbol"] = out["GeneSymbol"].astype(str).str.upper().str.strip()
    out["HPO_ID"] = out["HPO_ID"].astype(str).str.strip()
    return out

def phenotype_score(variants_df: pd.DataFrame, hpo_terms: List[str], hpo_map: Optional[pd.DataFrame]=None, panel_genes: Optional[List[str]]=None):
    """
    Adds columns:
    - PHENO_MATCH (bool)
    - PHENO_SCORE (int): count of HPO terms mapping to the gene
    - PANEL_MATCH (bool): if gene is in provided panel
    """
    v = variants_df.copy()
    v["PHENO_MATCH"] = False
    v["PHENO_SCORE"] = 0
    v["PANEL_MATCH"] = False

    if hpo_map is not None and len(hpo_terms)>0:
        hpo_terms = [t.strip() for t in hpo_terms if t.strip()]
        gene_counts = (hpo_map[hpo_map["HPO_ID"].isin(hpo_terms)]
                       .groupby("GeneSymbol").size().to_dict())
        v["GENE_UP"] = v["GENE"].fillna("").str.upper()
        v["PHENO_SCORE"] = v["GENE_UP"].map(gene_counts).fillna(0).astype(int)
        v["PHENO_MATCH"] = v["PHENO_SCORE"] > 0
        v.drop(columns=["GENE_UP"], inplace=True)

    if panel_genes is not None and len(panel_genes)>0:
        panel_set = set([g.strip().upper() for g in panel_genes if g and isinstance(g, str)])
        v["PANEL_MATCH"] = v["GENE"].fillna("").str.upper().isin(panel_set)

    return v

def prioritize(df: pd.DataFrame):
    """
    Simple prioritization score combining impact + phenotype + AF (rarer is higher).
    Score = w1*(impact_rank) + w2*(pheno_score) + w3*(rarity)
    """
    imp_rank = {"HIGH":3, "MODERATE":2, "LOW":1, "MODIFIER":0}
    f = df.copy()
    f["IMPACT_RANK"] = f["IMPACT"].map(imp_rank).fillna(0)
    # rarity: 1 - AF (unknown AF treated as 0.5)
    rarity = 1 - f["AF"].fillna(0.5).clip(0,1)
    f["PRIORITY_SCORE"] = (2.0*f["IMPACT_RANK"]) + (1.5*f["PHENO_SCORE"]) + (1.0*rarity)
    return f.sort_values(["PANEL_MATCH","PHENO_MATCH","PRIORITY_SCORE"], ascending=[False, False, False]).reset_index(drop=True)
