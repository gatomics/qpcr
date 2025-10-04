
# GATOMICS — VCF + Phenotype Analyzer (RUO)

Web app to parse VCF/VCF.GZ, combine with phenotype (HPO) and panels, prioritize variants, and export branded PDF/CSV.

## Inputs
- **VCF/VCF.GZ** (supports SnpEff `ANN` and VEP `CSQ` columns)
- **HPO→Gene map**: CSV/TSV with columns `HPO_ID, GeneSymbol`
- **Gene panel** (optional): TXT (one gene per line) or CSV (first column)

## Run locally
```bash
python -m venv .venv
source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy (Streamlit Cloud)
Push to GitHub, then in Streamlit Cloud choose this repo and `app.py`.

## RUO
For Research Use Only. Not for diagnostic use.
