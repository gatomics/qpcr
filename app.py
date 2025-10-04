
# Streamlit App: Gatomis VCF + Phenotype Analyzer (RUO)

import io, os
import pandas as pd
import streamlit as st
from datetime import datetime
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm

from vcf_pheno_core import parse_vcf, load_hpo_map, phenotype_score, prioritize

st.set_page_config(page_title="Gatomis VCF + Phenotype Analyzer (RUO)", layout="wide")

# Header with logo
col1, col2 = st.columns([1,3])
with col1:
    if os.path.exists("logo.png"):
        st.image("logo.png", width=160)
with col2:
    st.title("GATomics — VCF + Phenotype Analyzer (RUO)")
    st.caption("Upload a VCF (.vcf/.vcf.gz), provide HPO terms and optional panels, then review prioritized variants and export a branded PDF/CSV.")

with st.sidebar:
    st.header("Patient & Report")
    patient_id = st.text_input("Patient ID", value="PAT-001")
    sex = st.selectbox("Sex", ["Unknown","Male","Female"], index=0)
    age = st.text_input("Age", value="")
    report_title = st.text_input("Report Title", value="GATOMICS — VCF + Phenotype Report (RUO)")
    lab_name = st.text_input("Laboratory / Client", value="Gatomis AB")
    prepared_by = st.text_input("Prepared by", value="Analyst")
    reviewed_by = st.text_input("Reviewed by (optional)", value="")
    st.divider()
    st.header("Filters")
    pass_only = st.checkbox("PASS only", value=True)
    min_af = st.number_input("Min AF", 0.0, 1.0, 0.0, 0.01)
    min_dp = st.number_input("Min DP", 0.0, 100000.0, 0.0, 1.0)
    st.divider()
    st.header("Phenotype (HPO)")
    hpo_text = st.text_area("HPO terms (comma or newline separated, e.g., HP:0001250, HP:0000252)", height=100)
    hpo_file = st.file_uploader("HPO→Gene map (CSV/TSV with columns: HPO_ID, GeneSymbol)", type=["csv","tsv"])
    st.caption("Tip: you can export HPO-gene lists from external tools or provide your own curated mapping.")
    st.divider()
    panel_upload = st.file_uploader("Gene Panel (one symbol per line or CSV column)", type=["txt","csv"])
    st.caption("Optional: limit/prioritize to your gene list (e.g., cancer/rare disease panel).")
    st.divider()
    st.markdown("**RUO** · For research use only.")

uploaded = st.file_uploader("Upload VCF (.vcf/.vcf.gz)", type=["vcf","gz"])

def parse_panel(file) -> list:
    if file is None:
        return None
    try:
        if file.name.lower().endswith(".csv"):
            df = pd.read_csv(file)
            # try any column
            first_col = df.columns[0]
            return [str(x).strip() for x in df[first_col].dropna().unique().tolist()]
        else:
            txt = file.read().decode("utf-8", errors="ignore")
            return [l.strip() for l in txt.splitlines() if l.strip()]
    except Exception as e:
        st.warning(f"Could not parse panel: {e}")
        return None

def pdf_report(summary_lines, table: pd.DataFrame):
    buf = io.BytesIO()
    c = canvas.Canvas(buf, pagesize=A4)
    W,H=A4; margin=18*mm

    # header with logo
    if os.path.exists("logo.png"):
        try:
            c.drawImage("logo.png", margin, H - margin - 18*mm, width=38*mm, height=16*mm, mask='auto')
        except Exception:
            pass
    c.setFont("Helvetica-Bold", 14)
    c.drawRightString(W - margin, H - margin - 6*mm, report_title)
    c.setFont("Helvetica", 9)
    c.drawRightString(W - margin, H - margin - 12*mm, datetime.now().strftime("%Y-%m-%d %H:%M"))
    c.setFont("Helvetica", 9)
    c.drawString(margin, H - margin - 23*mm, f"Laboratory: {lab_name}")
    c.line(margin, H - margin - 25*mm, W - margin, H - margin - 25*mm)

    y = H - margin - 32*mm
    c.setFont("Helvetica-Bold", 11)
    c.drawString(margin, y, "Patient")
    y -= 6*mm
    c.setFont("Helvetica", 9)
    c.drawString(margin, y, f"Patient ID: {patient_id}   Sex: {sex}   Age: {age}")
    y -= 8*mm
    c.setFont("Helvetica-Bold", 11)
    c.drawString(margin, y, "Summary")
    y -= 6*mm
    c.setFont("Helvetica", 9)
    for line in summary_lines:
        if y < margin + 25*mm:
            c.showPage(); y = H - margin - 20*mm; c.setFont("Helvetica", 9)
        c.drawString(margin, y, line)
        y -= 5*mm

    # Table
    cols = ["CHROM","POS","REF","ALT","GENE","CONSEQUENCE","IMPACT","AF","DP","PHENO_SCORE","PANEL_MATCH","PRIORITY_SCORE"]
    if y < margin + 35*mm:
        c.showPage(); y = H - margin - 20*mm
    c.setFont("Helvetica-Bold", 11)
    c.drawString(margin, y, "Top Prioritized Variants")
    y -= 7*mm; c.setFont("Helvetica", 8)
    for _, row in table[cols].head(40).iterrows():
        line = f"{row['CHROM']}:{int(row['POS'])} {row['REF']}>{row['ALT']} | {str(row['GENE'])} | {str(row['CONSEQUENCE'])} | {str(row['IMPACT'])} | AF={row['AF']} | DP={row['DP']} | PHENO={row['PHENO_SCORE']} | PANEL={row['PANEL_MATCH']} | Score={round(row['PRIORITY_SCORE'],3)}"
        if y < margin + 20*mm:
            c.showPage(); y = H - margin - 20*mm; c.setFont("Helvetica", 8)
        c.drawString(margin, y, line[:150])
        y -= 4.5*mm

    # Signature
    if y < margin + 25*mm:
        c.showPage(); y = H - margin - 20*mm
    c.setFont("Helvetica-Bold", 10)
    c.drawString(margin, y, "Authorization")
    y -= 8*mm; c.setFont("Helvetica", 9)
    c.drawString(margin, y, f"Prepared by: {prepared_by}     Reviewed by: {reviewed_by}     Date: {datetime.now().strftime('%Y-%m-%d')}")
    c.line(margin, y-1*mm, margin+55*mm, y-1*mm)
    c.line(margin+75*mm, y-1*mm, margin+130*mm, y-1*mm)

    c.save()
    return buf.getvalue()

if uploaded is not None:
    with st.spinner("Parsing VCF..."):
        df, meta = parse_vcf(uploaded)
    st.success(f"Parsed {len(df)} variants. Samples: {', '.join(meta['samples']) if meta['samples'] else 'N/A'}")

    # Apply basic filters
    f = df.copy()
    if pass_only:
        f = f[f["FILTER"].isin(["PASS","."])]
    if min_af > 0:
        f = f[f["AF"].fillna(-1) >= min_af]
    if min_dp > 0:
        f = f[f["DP"].fillna(-1) >= min_dp]

    # Phenotype mapping
    panel_genes = parse_panel(panel_upload)
    hpo_terms = [t.strip() for t in (hpo_text.replace("\n",",").split(",")) if t.strip()]
    hpo_map = None
    if hpo_file is not None:
        try:
            import pandas as pd
            if hpo_file.name.lower().endswith(".tsv"):
                hpo_map = pd.read_csv(hpo_file, sep="\t")
            else:
                hpo_map = pd.read_csv(hpo_file)
            hpo_map = load_hpo_map(hpo_map)
            st.success(f"HPO map loaded with {len(hpo_map)} entries.")
        except Exception as e:
            st.error(f"Failed to load HPO map: {e}")

    scored = phenotype_score(f, hpo_terms, hpo_map=hpo_map, panel_genes=panel_genes)
    prioritized = prioritize(scored)

    # Summary lines
    s = meta["stats"]
    summary = [
        f"Total parsed: {len(df)} | After filters: {len(prioritized)}",
        f"PASS: {s['n_pass']} / Total: {s['n_total']} | SNPs: {s['snps']} | INDELs: {s['indels']} | Ti/Tv: {(s['ti']/s['tv']) if s['tv'] else 'NA'}",
        f"HPO terms: {', '.join(hpo_terms) if hpo_terms else 'None'}",
        f"Panel genes: {len(panel_genes) if panel_genes else 0}"
    ]

    st.subheader("📊 Summary")
    for line in summary: st.write("- " + line)

    st.subheader("🧬 Prioritized Variants")
    st.dataframe(prioritized, use_container_width=True)

    st.download_button("⬇️ Download prioritized CSV", data=prioritized.to_csv(index=False).encode("utf-8"), file_name=f"{patient_id}_prioritized.csv", mime="text/csv")
    pdf_bytes = pdf_report(summary, prioritized)
    st.download_button("⬇️ Download PDF report", data=pdf_bytes, file_name=f"{patient_id}_vcf_report.pdf", mime="application/pdf")

st.divider()
st.markdown("**Inputs:** VCF/VCF.GZ; optional HPO→Gene CSV/TSV; optional gene panel list/CSV.  \n**RUO:** For research use only.")
