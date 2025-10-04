
# Gatomis qPCR Analyzer (RUO)

A minimal, vendor-agnostic web app to analyze methylation qPCR Ct tables and generate a RUO PDF report.

## Features
- Upload CSV/XLSX with columns: `sample_id, marker, ct` (replicates as multiple rows)
- QC rules (configurable): reference Ct threshold, replicate SD threshold
- Marker positivity by Ct + ΔCt, final sample call (Positive/Indeterminate/Negative)
- PDF report export (summary + per-marker metrics)
- Streamlit UI; reusable Python core

## Quick Start
1. Create a virtual environment (Python 3.9+ recommended)
2. Install requirements
3. Run Streamlit

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Input Schema
CSV or Excel with columns:
```
sample_id, marker, ct
P001, SDC2, 37.2
P001, VIM, 39.0
P001, COL2A1, 30.1
P001, SFRP2, 38.5
```

## Thresholds (default v0.1.0)
- Reference Ct ≤ 32.0
- Replicate SD ≤ 0.75
- Marker Ct ≤ 41.0
- ΔCt ≤ 10.5 (ΔCt = Ct(marker) - Ct(reference))
- Final call:
  - ≥2 positive markers → Positive
  - 1 positive → Indeterminate
  - 0 positive → Negative

You can change thresholds in `gatomis_qpcr_core.py` (CONFIG dict).

## Files
- `app.py` — Streamlit UI + PDF export
- `gatomis_qpcr_core.py` — reusable analysis logic
- `sample_ct.csv` — example input
- `requirements.txt` — Python deps

## RUO
For Research Use Only. Not for diagnostic use.
