
# Gatomis qPCR Analyzer (RUO) - Streamlit App
# Run: streamlit run app.py

import io
import pandas as pd
import streamlit as st
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from datetime import datetime

from gatomis_qpcr_core import analyze_ct_table, sample_results_to_df, CONFIG

st.set_page_config(page_title="Gatomis qPCR Analyzer (RUO)", layout="wide")

st.title("🧪 Gatomis qPCR Analyzer (RUO)")
st.caption("Upload Ct table (CSV or Excel), set reference gene, and generate a RUO report.")

with st.sidebar:
    st.header("⚙️ Settings")
    ref_marker = st.text_input("Reference marker (default COL2A1)", value="COL2A1")
    report_title = st.text_input("Report Title", value="Gatomis qPCR Analyzer (RUO)")
    lab_name = st.text_input("Laboratory / Client Name", value="Gatomis AB")
    prepared_by = st.text_input("Prepared by", value="Analyst")
    reviewed_by = st.text_input("Reviewed by (optional)", value="")
    language = st.selectbox("Report Language", ["English", "Arabic", "Bilingual"])
    st.write("**Thresholds (v%s)**" % CONFIG["version"])
    st.write(f"- Ref Ct ≤ **{CONFIG['qc']['ref_ct_max']}**")
    st.write(f"- Replicate SD ≤ **{CONFIG['qc']['replicate_sd_max']}**")
    st.write(f"- Marker Ct ≤ **{CONFIG['marker_call']['marker_ct_max']}**")
    st.write(f"- ΔCt ≤ **{CONFIG['marker_call']['delta_ct_max']}**")
    st.write(f"- ≥ **{CONFIG['sample_call']['min_positive_markers']}** positive markers → **Positive**")
    st.divider()
    st.markdown("**RUO** · For Research Use Only · Not for diagnostic use.")

uploaded = st.file_uploader("Upload Ct table (CSV or XLSX)", type=["csv", "xlsx", "xls"])

def load_df(file) -> pd.DataFrame:
    if file.name.lower().endswith(".csv"):
        return pd.read_csv(file)
    else:
        return pd.read_excel(file)


def _tr(text_en, text_ar, language):
    if language == "English":
        return text_en
    elif language == "Arabic":
        return text_ar
    else:
        return f"{text_en} / {text_ar}"

def generate_pdf(sample_df: pd.DataFrame, detailed_df: pd.DataFrame, ref_marker: str,
                 report_title: str, lab_name: str, prepared_by: str, reviewed_by: str, language: str) -> bytes:
    buffer = io.BytesIO()
    c = canvas.Canvas(buffer, pagesize=A4)
    W, H = A4
    margin = 18*mm

    logo_path = "logo.png"
    def header():
        if os.path.exists(logo_path):
            try:
                c.drawImage(logo_path, margin, H - margin - 18*mm, width=38*mm, height=16*mm, mask='auto')
            except Exception:
                pass
        c.setFont("Helvetica-Bold", 14)
        c.drawRightString(W - margin, H - margin - 6*mm, report_title)
        c.setFont("Helvetica", 9)
        c.drawRightString(W - margin, H - margin - 12*mm, datetime.now().strftime("%Y-%m-%d %H:%M"))
        c.setFont("Helvetica", 9)
        c.drawString(margin, H - margin - 23*mm, _tr("Laboratory", "المختبر", language) + f": {lab_name}")
        c.line(margin, H - margin - 25*mm, W - margin, H - margin - 25*mm)

    def footer(page_num):
        c.setFont("Helvetica", 8)
        c.drawCentredString(W/2, margin/2, f"RUO — For Research Use Only | صفحة/Page {page_num}")

    page_num = 1
    header()
    y = H - margin - 32*mm

    c.setFont("Helvetica-Bold", 11)
    c.drawString(margin, y, _tr(f"Summary by Sample (Reference: {ref_marker})", f'ملخص حسب العينة (الجين المرجعي: {ref_marker})', language))
    y -= 7*mm

    c.setFont("Helvetica", 9)
    for _, row in sample_df.iterrows():
        ref_qc = "PASS" if row['ref_qc_pass'] else "FAIL"
        call = row['final_call']
        if language != "English":
            call = {"Positive":"موجب","Negative":"سلبي","Indeterminate":"غير حاسم","Invalid (Reference QC Fail)":"غير صالح (فشل ضبط الجودة)"}.get(call, call)
            ref_qc = {"PASS":"نجح","FAIL":"فشل"}.get(ref_qc, ref_qc)
        line = f"{row['sample_id']}  |  ref_ct={row['ref_ct']:.2f}  |  QC={ref_qc}  |  pos_markers={int(row['n_markers_positive'])}  |  call={call}"
        if y < margin + 28*mm:
            footer(page_num); c.showPage(); page_num += 1; header(); y = H - margin - 32*mm; c.setFont("Helvetica", 9)
        c.drawString(margin, y, line)
        y -= 6*mm

    if y < margin + 35*mm:
        footer(page_num); c.showPage(); page_num += 1; header(); y = H - margin - 32*mm

    c.setFont("Helvetica-Bold", 11)
    c.drawString(margin, y, _tr("Detailed Marker Metrics", "تفاصيل المؤشرات الجينية", language))
    y -= 7*mm
    c.setFont("Helvetica", 8)

    for _, row in detailed_df.iterrows():
        qc_text = "PASS" if row['qc_pass'] else "FAIL"
        if language != "English":
            qc_text = {"PASS":"نجح","FAIL":"فشل"}.get(qc_text, qc_text)
        line = f"{row['sample_id']} | {row['marker']}: n={int(row['n_reps'])}, ct_mean={row['ct_mean']}, sd={row['ct_sd']}, ref_ct={row['ref_ct']}, dCt={row['delta_ct']}, QC={qc_text}, positive={'YES' if row['marker_positive'] else 'NO'}"
        if y < margin + 35*mm:
            footer(page_num); c.showPage(); page_num += 1; header(); y = H - margin - 32*mm; c.setFont("Helvetica", 8)
        c.drawString(margin, y, line)
        y -= 5*mm

    if y < margin + 25*mm:
        footer(page_num); c.showPage(); page_num += 1; header(); y = H - margin - 32*mm

    c.setFont("Helvetica-Bold", 10)
    c.drawString(margin, y, _tr("Authorization", "التوقيع والاعتماد", language))
    y -= 10*mm
    c.setFont("Helvetica", 9)

    c.drawString(margin, y, _tr("Prepared by:", "أعدّه:", language) + f" {prepared_by}")
    c.line(margin+35*mm, y-1*mm, margin+90*mm, y-1*mm)

    y -= 10*mm
    c.drawString(margin, y, _tr("Reviewed by:", "راجعه:", language) + f" {reviewed_by}")
    c.line(margin+35*mm, y-1*mm, margin+90*mm, y-1*mm)

    y -= 10*mm
    c.drawString(margin, y, _tr("Date:", "التاريخ:", language) + f" {datetime.now().strftime('%Y-%m-%d')}")
    c.line(margin+35*mm, y-1*mm, margin+90*mm, y-1*mm)

    footer(page_num)
    c.save()
    pdf_bytes = buffer.getvalue()
    buffer.close()
    return pdf_bytes


if uploaded is not None:
    try:
        df = load_df(uploaded)
        st.success(f"Loaded file with {len(df)} rows.")
        detailed_df, sample_results = analyze_ct_table(df, reference_marker=ref_marker)
        sample_df = sample_results_to_df(sample_results)

        st.subheader("📊 Sample Summary")
        st.dataframe(sample_df, use_container_width=True)

        st.subheader("🔬 Per‑Marker Details")
        st.dataframe(detailed_df, use_container_width=True)

        # PDF
        pdf_bytes = generate_pdf(sample_df, detailed_df, ref_marker, report_title, lab_name, prepared_by, reviewed_by, language)
        st.download_button("⬇️ Download PDF Report", data=pdf_bytes, file_name="Gatomis_qPCR_Report.pdf", mime="application/pdf")

    except Exception as e:
        st.error(f"Error: {e}")

st.divider()
st.markdown("**Schema:** CSV/XLSX with columns: `sample_id, marker, ct` (replicates = multiple rows).")
st.markdown("_Example markers: SDC2, VIM, SFRP2, NPY; Reference: COL2A1._")
