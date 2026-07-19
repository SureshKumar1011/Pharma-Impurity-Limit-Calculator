import streamlit as st
import pandas as pd

# --- Helper Functions ---
def format_percent(value: float, capped: bool = False, substance: bool = False) -> str:
    """
    Format percentages according to rules:
    - Drug substance fixed thresholds: show 2 decimals (0.05%, 0.10%, 0.15%)
    - Drug product fixed thresholds: show 1 decimal (1.0%, 0.5%, 0.2%, 0.1%)
    - Capped/calculated values: show up to 3 decimals
    """
    if capped:
        return f"{value:.3f}%"
    else:
        if substance:
            return f"{value:.2f}%"
        else:
            return f"{value:.1f}%"

def calc_threshold(dose_mg, percent, cap_value, cap_unit="mg", substance=False):
    cap_mg = cap_value if cap_unit == "mg" else cap_value / 1000.0
    percent_value_mg = dose_mg * (percent / 100)

    if percent_value_mg <= cap_mg:
        return format_percent(percent, substance=substance)
    else:
        applied_percent = (cap_mg / dose_mg) * 100
        return f"{format_percent(applied_percent, capped=True)} (capped at {cap_value} {cap_unit}/day)"

# --- Drug Substance (ICH Q3A) ---
def ich_q3a_limits(dose_mg: float) -> dict:
    dose_g = dose_mg / 1000.0
    if dose_g <= 2:  # ≤ 2 g/day
        # Identification threshold: 0.10% or 1.0 mg/day
        id_threshold_mg = 0.001 * dose_mg
        if id_threshold_mg > 1.0:
            id_percent = (1.0 / dose_mg) * 100
            id_value = f"{id_percent:.2f}% (capped at 1.0 mg/day)"
        else:
            id_value = "0.10%"

        # Qualification threshold: 0.15% or 1.0 mg/day
        qt_threshold_mg = 0.0015 * dose_mg
        if qt_threshold_mg > 1.0:
            qt_percent = (1.0 / dose_mg) * 100
            qt_value = f"{qt_percent:.2f}% (capped at 1.0 mg/day)"
        else:
            qt_value = "0.15%"

        return {
            "Reporting threshold": "0.05%",
            "Identification threshold": id_value,
            "Qualification threshold": qt_value
        }
    else:  # > 2 g/day
        return {
            "Reporting threshold": "0.03%",
            "Identification threshold": "0.05%",
            "Qualification threshold": "0.05%"
        }


# --- Drug Product (ICH Q3B) ---
def calc_threshold(dose_mg, percent, cap_value, cap_unit):
    """Calculate threshold with % and cap logic, formatting per ICH style."""
    threshold_mg = (percent / 100.0) * dose_mg
    if cap_unit == "µg":
        cap_mg = cap_value / 1000.0
    else:
        cap_mg = cap_value

    if threshold_mg > cap_mg:
        capped_percent = (cap_mg / dose_mg) * 100
        return f"{capped_percent:.2f}% (capped at {cap_value} {cap_unit}/day)"
    else:
        # No capping → display as ICH % (no extra digits)
        return f"{percent:.1f}%"

def ich_q3b_limits(dose_mg):
    """Return Reporting, Identification, and Qualification thresholds for Drug Product (ICH Q3B)."""
    limits = {}

    # Reporting threshold
    if dose_mg <= 1000:  # ≤ 1 g
        limits["Reporting threshold"] = "0.1%"
    else:  # > 1 g
        limits["Reporting threshold"] = "0.05%"

    # Identification threshold
    if dose_mg < 1:
        limits["Identification threshold"] = calc_threshold(dose_mg, 1.0, 5, "µg")
    elif dose_mg <= 10:  # 1–10 mg band
        limits["Identification threshold"] = calc_threshold(dose_mg, 0.5, 20, "µg")
    elif dose_mg <= 2000:  # >10 mg – 2 g band
        limits["Identification threshold"] = calc_threshold(dose_mg, 0.2, 2.0, "mg")
    else:  # > 2 g
        limits["Identification threshold"] = "0.10%"

    # Qualification threshold
    if dose_mg < 10:
        limits["Qualification threshold"] = calc_threshold(dose_mg, 1.0, 50, "µg")
    elif dose_mg <= 100:  # 10–100 mg band
        limits["Qualification threshold"] = calc_threshold(dose_mg, 0.5, 200, "µg")
    elif dose_mg <= 2000:  # >100 mg – 2 g band
        limits["Qualification threshold"] = calc_threshold(dose_mg, 0.2, 3.0, "mg")
    else:  # > 2 g
        limits["Qualification threshold"] = "0.15%"

    return limits

# --- Mutagenic Impurities (ICH M7 TTC) ---
def ich_m7_ttc(dose_mg: float, duration_value: int, duration_unit: str, qual_threshold: str) -> str:
    if duration_value == 0:
        duration_months = 9999
    elif duration_unit == "Days":
        duration_months = duration_value / 30.0
    elif duration_unit == "Months":
        duration_months = duration_value
    else:
        duration_months = duration_value * 12

    if duration_months <= 1:
        ttc_ug = 120
        label = "≤ 1 month"
    elif duration_months <= 12:
        ttc_ug = 20
        label = "1–12 months"
    elif duration_months <= 120:
        ttc_ug = 10
        label = "1–10 years"
    else:
        ttc_ug = 1.5
        label = ">10 years to lifetime"

    ttc_mg = ttc_ug / 1000.0
    ppm = (ttc_mg / dose_mg) * 1e6
    percent = (ttc_mg / dose_mg) * 100

    result = f"{ttc_ug} µg/day ({label}) → {ppm:.3f} ppm"
    if ppm >= 1000:
        result += f" ({percent:.3f}%)"

    # Cap against ICH Q3 qualification threshold
    try:
        qual_percent = float(qual_threshold.split("%")[0])
        if percent > qual_percent:
            result += f" → Capped at Qualification threshold = {qual_threshold}"
    except:
        pass

    return result

# --- Nitrosamines ---
def nitrosamine_limits(dose_mg: float, market: str, qual_threshold: str) -> dict:
    ai_values = {
        "Category 1 (US)": 26.5,
        "Category 1 (Global)": 18,
        "Category 2": 100,
        "Category 3": 400,
        "Category 4": 1500,
        "Category 5": 1500
    }
    categories = {}
    for cat, ai in ai_values.items():
        if "Global" in cat and market == "US":
            continue
        if "US" in cat and market == "Global":
            continue

        ppm = ai / dose_mg
        result = f"{ppm:.3f} ppm ({ai} ng/day)"
        if ppm >= 1000:
            percent = ppm / 10000
            result += f" ({percent:.3f}%)"

        # Cap against ICH Q3 qualification threshold
        try:
            qual_percent = float(qual_threshold.split("%")[0])
            qual_ppm = qual_percent * 10000
            if ppm > qual_ppm:
                result += f" → Capped at Qualification threshold = {qual_threshold}"
        except:
            pass

        categories[cat] = result
    return categories

# --- Streamlit UI ---
st.title("Pharma Impurity Limit Calculator")

st.markdown("### Enter Drug Information")

# First row: Drug name, Category, Market
col1, col2, col3 = st.columns(3)
with col1:
    drug_name = st.text_input("Drug name:")
with col2:
    category = st.selectbox("Category:", ["Drug substance", "Drug product"])
with col3:
    market = st.selectbox("Market:", ["US", "Global"])

# Second row: Dose, Duration value, Duration unit
col1, col2, col3 = st.columns(3)
with col1:
    dose = st.number_input("Max daily dose (mg):", min_value=0.0, format="%.6g")
with col2:
    duration_value = st.number_input("Treatment duration:", min_value=0, step=1)
with col3:
    duration_unit = st.selectbox("Duration unit:", ["Days", "Months", "Years"])

# Results
if dose > 0 and drug_name:
    if category == "Drug substance":
        limits = ich_q3a_limits(dose)
    else:
        limits = ich_q3b_limits(dose)

    qual_threshold = limits["Qualification threshold"]

    m7_limit = ich_m7_ttc(dose, duration_value, duration_unit, qual_threshold)
    nitrosamine_ai = nitrosamine_limits(dose, market, qual_threshold)

    data = [
        ["Reporting threshold", limits["Reporting threshold"]],
        ["Identification threshold", limits["Identification threshold"]],
        ["Qualification threshold", limits["Qualification threshold"]],
        ["Mutagenic alert impurity (ICH M7)", m7_limit]
    ]
    for cat, val in nitrosamine_ai.items():
        data.append([f"Nitrosamine {cat}", val])

    df = pd.DataFrame(data, columns=["Threshold / Category", "Limit"])
    st.markdown("### Results")
    st.table(df)
else:
    st.info("Please enter drug name and dose to see impurity limits.")

# --- Footer (email only, hidden when printing) ---
st.markdown("---")
st.markdown("**Contact: dr.suresh_11@outlook.com**")

# CSS to hide footer when printing
hide_footer_css = """
@media print {
    .stMarkdown {
        display: none;
    }
}
"""
st.markdown(f"<style>{hide_footer_css}</style>", unsafe_allow_html=True)
