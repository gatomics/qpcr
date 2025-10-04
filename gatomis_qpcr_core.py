
"""
Gatomis qPCR Analyzer (RUO) - Core Logic
Author: Ali + ChatGPT

This module provides reusable functions to parse Ct tables (CSV/XLSX),
apply QC rules, compute per-marker calls, and produce a final sample call.
"""

from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple
import pandas as pd
import numpy as np

# ----------- Config (versioned) -----------
CONFIG = {
    "version": "0.1.0",
    "qc": {
        "ref_ct_max": 32.0,        # reference gene must be <= this Ct
        "replicate_sd_max": 0.75,  # max allowed SD across replicates (per marker)
    },
    "marker_call": {
        "marker_ct_max": 41.0,
        "delta_ct_max": 10.5,      # ΔCt = Ct(marker) - Ct(ref)
    },
    "sample_call": {
        "min_positive_markers": 2  # >=2 positive markers => Positive
    }
}

@dataclass
class MarkerResult:
    sample_id: str
    marker: str
    n_reps: int
    ct_mean: float
    ct_sd: float
    ref_ct: float
    delta_ct: float
    qc_pass: bool
    marker_positive: bool

@dataclass
class SampleResult:
    sample_id: str
    ref_marker: str
    ref_ct: float
    ref_qc_pass: bool
    n_markers_positive: int
    final_call: str
    marker_results: List[MarkerResult]
    notes: List[str]

def _clean_df(df: pd.DataFrame) -> pd.DataFrame:
    # Normalize columns
    cols = {c.lower().strip(): c for c in df.columns}
    # Expect essential columns
    mapping = {}
    for key in ["sample_id", "marker", "ct"]:
        if key in cols:
            mapping[cols[key]] = key
        else:
            # try fuzzy
            for c in df.columns:
                lc = c.lower().strip()
                if key == "sample_id" and lc in ["sample", "id", "sample id", "sample_id"]:
                    mapping[c] = "sample_id"
                elif key == "marker" and lc in ["gene", "target", "assay", "marker"]:
                    mapping[c] = "marker"
                elif key == "ct" and lc in ["cq", "ct", "c_t", "ct value"]:
                    mapping[c] = "ct"
    if set(mapping.values()) != {"sample_id", "marker", "ct"}:
        raise ValueError("Input must have columns: sample_id, marker, ct (or equivalents).")

    df2 = df.rename(columns=mapping)[["sample_id", "marker", "ct"]].copy()
    # Coerce types
    df2["sample_id"] = df2["sample_id"].astype(str).str.strip()
    df2["marker"] = df2["marker"].astype(str).str.strip()
    # Ensure numeric ct
    df2["ct"] = pd.to_numeric(df2["ct"], errors="coerce")
    df2 = df2.dropna(subset=["sample_id", "marker", "ct"])
    return df2

def analyze_ct_table(df: pd.DataFrame, reference_marker: str = "COL2A1") -> Tuple[pd.DataFrame, List[SampleResult]]:
    """
    Input df with columns: sample_id, marker, ct (replicates allowed as multiple rows)
    Returns: (wide_table_df, sample_results_list)
    """
    df2 = _clean_df(df)

    # Group to compute replicate stats
    agg = df2.groupby(["sample_id", "marker"]).agg(
        ct_mean=("ct", "mean"),
        ct_sd=("ct", "std"),
        n_reps=("ct", "count")
    ).reset_index()

    # Fill NaN SD with 0 for single replicate
    agg["ct_sd"] = agg["ct_sd"].fillna(0.0)

    sample_results: List[SampleResult] = []
    rows = []

    for sample_id, g in agg.groupby("sample_id"):
        notes = []
        # Reference
        ref_rows = g[g["marker"].str.upper() == reference_marker.upper()]
        if ref_rows.empty:
            ref_ct = np.nan
            ref_qc_pass = False
            notes.append(f"Missing reference marker {reference_marker}.")
        else:
            ref_ct = float(ref_rows.iloc[0]["ct_mean"])
            ref_sd = float(ref_rows.iloc[0]["ct_sd"])
            # QC for reference: Ct threshold + replicate SD
            ref_qc_pass = (ref_ct <= CONFIG["qc"]["ref_ct_max"]) and (ref_sd <= CONFIG["qc"]["replicate_sd_max"])

        marker_results: List[MarkerResult] = []
        n_pos = 0

        for _, r in g.iterrows():
            marker = str(r["marker"])
            ct_mean = float(r["ct_mean"])
            ct_sd = float(r["ct_sd"])
            n_reps = int(r["n_reps"])

            # compute ΔCt using ref_ct if available
            delta_ct = ct_mean - ref_ct if pd.notna(ref_ct) else np.nan

            # marker QC: replicate SD
            m_qc = (ct_sd <= CONFIG["qc"]["replicate_sd_max"])

            # marker positive logic
            marker_positive = False
            if ref_qc_pass and m_qc and pd.notna(ref_ct):
                if (ct_mean <= CONFIG["marker_call"]["marker_ct_max"]) and (delta_ct <= CONFIG["marker_call"]["delta_ct_max"]):
                    marker_positive = True

            if marker_positive and marker.upper() != reference_marker.upper():
                n_pos += 1

            mr = MarkerResult(
                sample_id=sample_id,
                marker=marker,
                n_reps=n_reps,
                ct_mean=ct_mean,
                ct_sd=ct_sd,
                ref_ct=ref_ct if pd.notna(ref_ct) else np.nan,
                delta_ct=delta_ct if pd.notna(delta_ct) else np.nan,
                qc_pass=bool(m_qc),
                marker_positive=bool(marker_positive if marker.upper()!=reference_marker.upper() else False)
            )
            marker_results.append(mr)

            rows.append({
                "sample_id": sample_id,
                "marker": marker,
                "n_reps": n_reps,
                "ct_mean": round(ct_mean, 3),
                "ct_sd": round(ct_sd, 3),
                "ref_ct": None if pd.isna(ref_ct) else round(ref_ct, 3),
                "delta_ct": None if pd.isna(delta_ct) else round(delta_ct, 3),
                "qc_pass": m_qc,
                "marker_positive": (marker_positive if marker.upper()!=reference_marker.upper() else False)
            })

        # Final call
        if not ref_qc_pass:
            final_call = "Invalid (Reference QC Fail)"
        else:
            min_pos = CONFIG["sample_call"]["min_positive_markers"]
            if n_pos >= min_pos:
                final_call = "Positive"
            elif n_pos == 1:
                final_call = "Indeterminate"
            else:
                final_call = "Negative"

        sres = SampleResult(
            sample_id=sample_id,
            ref_marker=reference_marker,
            ref_ct=float(ref_ct) if pd.notna(ref_ct) else np.nan,
            ref_qc_pass=bool(ref_qc_pass),
            n_markers_positive=int(n_pos),
            final_call=final_call,
            marker_results=marker_results,
            notes=notes
        )
        sample_results.append(sres)

    wide = pd.DataFrame(rows)
    return wide, sample_results

def sample_results_to_df(sample_results: List[SampleResult]) -> pd.DataFrame:
    rows = []
    for s in sample_results:
        rows.append({
            "sample_id": s.sample_id,
            "ref_marker": s.ref_marker,
            "ref_ct": s.ref_ct,
            "ref_qc_pass": s.ref_qc_pass,
            "n_markers_positive": s.n_markers_positive,
            "final_call": s.final_call,
            "notes": "; ".join(s.notes) if s.notes else ""
        })
    return pd.DataFrame(rows)
