"""
streamlit_app.py — ADME Analysis Dashboard
Run with:  streamlit run streamlit_app.py
Requires adme_analysis.py and literature_intelligence.py in the same folder.
"""

import io, time
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adme_analysis import (
    analyze_smiles,
    build_adme_panels,
    build_adme_figure,
    generate_optimization_advice,
)
from literature_intelligence import get_literature_intelligence

st.set_page_config(page_title="LeadRefine · ADME Analyzer", page_icon="🧬", layout="wide")

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Syne:wght@400;500;600;700;800&family=Plus+Jakarta+Sans:ital,wght@0,300;0,400;0,500;0,600;0,700;1,400&family=JetBrains+Mono:wght@400;500&display=swap');

/* ── Design Tokens ──────────────────────────────────────────── */
:root {
    --bg:          #05070F;
    --surface:     #0C1220;
    --card:        #0F1825;
    --card-hi:     #131F30;
    --border:      #1A2840;
    --border-hi:   #253B5C;
    --accent:      #00C8E8;
    --accent2:     #F5A623;
    --accent3:     #9B71F0;
    --accent-glow: rgba(0, 200, 232, 0.18);
    --text:        #EDF4FF;
    --text-muted:  #7A8FA8;
    --text-dim:    #3D5068;
    --success:     #22C983;
    --danger:      #F05252;
    --warning:     #F5A623;
    --font-head:   'Syne', sans-serif;
    --font-body:   'Plus Jakarta Sans', sans-serif;
    --font-mono:   'JetBrains Mono', monospace;
    --radius:      12px;
    --radius-sm:   8px;
    --radius-pill: 30px;
}

/* ── Base ───────────────────────────────────────────────────── */
html, body, .stApp, [data-testid="stAppViewContainer"] {
    background-color: var(--bg) !important;
    font-family: var(--font-body) !important;
    color: var(--text) !important;
}
.main .block-container {
    padding: 2.2rem 3rem 5rem !important;
    max-width: 1380px !important;
}

/* ── Typography ─────────────────────────────────────────────── */
h1 {
    font-family: var(--font-head) !important;
    font-size: 2.55rem !important;
    font-weight: 800 !important;
    background: linear-gradient(125deg, #00C8E8 0%, #9B71F0 55%, #F5A623 100%);
    -webkit-background-clip: text !important;
    -webkit-text-fill-color: transparent !important;
    background-clip: text !important;
    letter-spacing: -0.6px !important;
    line-height: 1.15 !important;
    margin-bottom: 0.4rem !important;
}
h2, h3 {
    font-family: var(--font-head) !important;
    color: var(--text) !important;
}
p, li {
    font-family: var(--font-body) !important;
    color: var(--text) !important;
    font-size: 0.92rem !important;
    line-height: 1.68 !important;
}
code {
    font-family: var(--font-mono) !important;
    background: var(--surface) !important;
    color: var(--accent) !important;
    border: 1px solid var(--border) !important;
    border-radius: 5px !important;
    padding: 2px 8px !important;
    font-size: 0.82rem !important;
}
hr {
    border: none !important;
    border-top: 1px solid var(--border) !important;
    margin: 1.8rem 0 !important;
}

/* ── Sidebar ────────────────────────────────────────────────── */
[data-testid="stSidebar"] {
    background: linear-gradient(180deg, #060A14 0%, #08111F 100%) !important;
    border-right: 1px solid var(--border) !important;
}
[data-testid="stSidebar"] * {
    font-family: var(--font-body) !important;
}
[data-testid="stSidebar"] .stMarkdown h2 {
    font-family: var(--font-head) !important;
    font-size: 1.1rem !important;
    font-weight: 800 !important;
    background: linear-gradient(90deg, var(--accent), var(--accent3));
    -webkit-background-clip: text !important;
    -webkit-text-fill-color: transparent !important;
    background-clip: text !important;
    letter-spacing: 0.3px !important;
    margin-bottom: 0 !important;
}
[data-testid="stSidebar"] .stRadio label span {
    font-size: 0.88rem !important;
    color: var(--text-muted) !important;
    font-weight: 500 !important;
}
[data-testid="stSidebar"] .stSelectbox label {
    font-size: 0.74rem !important;
    font-weight: 700 !important;
    text-transform: uppercase !important;
    letter-spacing: 0.8px !important;
    color: var(--text-dim) !important;
}
[data-testid="stSidebar"] .stCaption p {
    font-size: 0.74rem !important;
    color: var(--text-dim) !important;
    line-height: 1.5 !important;
}

/* ── Inputs ─────────────────────────────────────────────────── */
.stTextInput > div > div > input,
.stTextArea  > div > div > textarea {
    background-color: var(--surface) !important;
    color: var(--text) !important;
    border: 1.5px solid var(--border) !important;
    border-radius: var(--radius) !important;
    font-family: var(--font-mono) !important;
    font-size: 0.86rem !important;
    padding: 12px 16px !important;
    transition: border-color 0.2s, box-shadow 0.2s !important;
}
.stTextInput > div > div > input:focus,
.stTextArea  > div > div > textarea:focus {
    border-color: var(--accent) !important;
    box-shadow: 0 0 0 3px var(--accent-glow) !important;
    outline: none !important;
}
.stTextInput > label, .stTextArea > label {
    font-family: var(--font-body) !important;
    font-size: 0.76rem !important;
    font-weight: 700 !important;
    color: var(--text-dim) !important;
    letter-spacing: 0.8px !important;
    text-transform: uppercase !important;
    margin-bottom: 6px !important;
}
.stTextArea > div > div > textarea::placeholder,
.stTextInput > div > div > input::placeholder {
    color: var(--text-dim) !important;
    font-family: var(--font-mono) !important;
}

/* ── Buttons ────────────────────────────────────────────────── */
.stButton > button {
    font-family: var(--font-head) !important;
    font-weight: 700 !important;
    font-size: 0.88rem !important;
    letter-spacing: 0.4px !important;
    border-radius: var(--radius) !important;
    padding: 10px 22px !important;
    border: 1.5px solid var(--border) !important;
    background: var(--surface) !important;
    color: var(--text-muted) !important;
    transition: all 0.2s ease !important;
    cursor: pointer !important;
}
.stButton > button:hover {
    background: var(--card-hi) !important;
    border-color: var(--accent) !important;
    color: var(--accent) !important;
    box-shadow: 0 0 18px var(--accent-glow) !important;
    transform: translateY(-1px) !important;
}
[data-testid="baseButton-primary"] {
    background: linear-gradient(130deg, #0090AA 0%, #6B46D9 100%) !important;
    border: none !important;
    color: #fff !important;
    box-shadow: 0 4px 22px rgba(0, 200, 232, 0.28) !important;
    font-size: 0.92rem !important;
    padding: 12px 28px !important;
}
[data-testid="baseButton-primary"]:hover {
    background: linear-gradient(130deg, #00B8D4 0%, #8B5CF6 100%) !important;
    box-shadow: 0 6px 32px rgba(0, 200, 232, 0.42) !important;
    transform: translateY(-2px) !important;
}

/* ── Download Button ────────────────────────────────────────── */
.stDownloadButton > button {
    font-family: var(--font-body) !important;
    font-size: 0.8rem !important;
    font-weight: 600 !important;
    background: transparent !important;
    border: 1px solid var(--border) !important;
    color: var(--text-muted) !important;
    border-radius: var(--radius-sm) !important;
    padding: 8px 16px !important;
    transition: all 0.2s !important;
}
.stDownloadButton > button:hover {
    border-color: var(--accent2) !important;
    color: var(--accent2) !important;
    background: rgba(245, 166, 35, 0.06) !important;
}

/* ── Toggle ─────────────────────────────────────────────────── */
.stToggle label p {
    font-family: var(--font-body) !important;
    font-size: 0.9rem !important;
    font-weight: 600 !important;
    color: var(--text) !important;
}

/* ── Select ─────────────────────────────────────────────────── */
.stSelectbox > div > div {
    background: var(--surface) !important;
    border: 1.5px solid var(--border) !important;
    border-radius: var(--radius) !important;
    color: var(--text) !important;
    font-family: var(--font-body) !important;
    font-size: 0.88rem !important;
}

/* ── Metrics ────────────────────────────────────────────────── */
[data-testid="stMetric"] {
    background: var(--card) !important;
    border: 1px solid var(--border) !important;
    border-radius: 14px !important;
    padding: 18px 20px !important;
    position: relative !important;
    overflow: hidden !important;
    transition: border-color 0.25s, transform 0.2s !important;
}
[data-testid="stMetric"]:hover {
    border-color: var(--accent) !important;
    transform: translateY(-2px) !important;
}
[data-testid="stMetric"]::after {
    content: '';
    position: absolute;
    top: 0; left: 0; right: 0;
    height: 2px;
    background: linear-gradient(90deg, var(--accent), var(--accent3));
    border-radius: 14px 14px 0 0;
}
[data-testid="stMetricLabel"] p {
    font-family: var(--font-body) !important;
    font-size: 0.72rem !important;
    font-weight: 700 !important;
    text-transform: uppercase !important;
    letter-spacing: 1px !important;
    color: var(--text-muted) !important;
}
[data-testid="stMetricValue"] {
    font-family: var(--font-head) !important;
    font-size: 1.55rem !important;
    font-weight: 700 !important;
    color: var(--text) !important;
}

/* ── Expanders ──────────────────────────────────────────────── */
[data-testid="stExpander"] {
    background: var(--card) !important;
    border: 1px solid var(--border) !important;
    border-radius: var(--radius) !important;
    margin-bottom: 10px !important;
    overflow: hidden !important;
    transition: border-color 0.2s !important;
}
[data-testid="stExpander"]:hover {
    border-color: var(--border-hi) !important;
}
[data-testid="stExpander"] summary {
    font-family: var(--font-body) !important;
    font-weight: 600 !important;
    font-size: 0.92rem !important;
    color: var(--text) !important;
    padding: 14px 18px !important;
}
[data-testid="stExpanderDetails"] {
    padding: 4px 18px 18px !important;
    border-top: 1px solid var(--border) !important;
}

/* ── Progress bar ───────────────────────────────────────────── */
.stProgress > div > div {
    background: var(--border) !important;
    border-radius: 4px !important;
    height: 5px !important;
}
.stProgress > div > div > div {
    background: linear-gradient(90deg, var(--accent), var(--accent3)) !important;
    border-radius: 4px !important;
    transition: width 0.3s ease !important;
}

/* ── Alerts ─────────────────────────────────────────────────── */
[data-testid="stAlert"] {
    border-radius: var(--radius) !important;
    font-family: var(--font-body) !important;
    font-size: 0.88rem !important;
}

/* ── Caption ────────────────────────────────────────────────── */
.stCaption p {
    color: var(--text-dim) !important;
    font-size: 0.76rem !important;
    font-family: var(--font-body) !important;
}

/* ── DataFrames ─────────────────────────────────────────────── */
[data-testid="stDataFrame"] {
    border: 1px solid var(--border) !important;
    border-radius: var(--radius) !important;
    overflow: hidden !important;
}

/* ── Spinner ────────────────────────────────────────────────── */
[data-testid="stSpinner"] > div {
    border-top-color: var(--accent) !important;
}

/* ── ──────────────────────────────────────────────────────── */
/* CUSTOM HTML COMPONENT CLASSES                               */
/* ─────────────────────────────────────────────────────────── */

/* Section headings */
.sec-head {
    font-family: var(--font-head);
    font-size: 1.08rem;
    font-weight: 700;
    color: var(--accent);
    display: flex;
    align-items: center;
    gap: 7px;
    padding: 8px 0 10px;
    margin: 1.6rem 0 0.7rem;
    border-bottom: 1px solid var(--border);
    letter-spacing: 0.2px;
}
.sub-head {
    font-family: var(--font-head);
    font-size: 0.78rem;
    font-weight: 700;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 1px;
    margin: 1rem 0 0.35rem;
}

/* Glossary */
.gloss-term {
    font-family: var(--font-head);
    font-weight: 700;
    font-size: 0.87rem;
    color: var(--accent);
    display: block;
    margin-top: 0.7rem;
}
.gloss-def {
    color: var(--text-muted);
    font-family: var(--font-body);
    font-size: 0.81rem;
    line-height: 1.65;
}

/* PubChem block */
.pubchem-name {
    font-family: var(--font-head);
    font-size: 1.65rem;
    font-weight: 800;
    color: var(--text);
    margin-bottom: 6px;
    letter-spacing: -0.4px;
    line-height: 1.2;
}
.pubchem-meta {
    font-family: var(--font-body);
    font-size: 0.86rem;
    color: var(--text-muted);
    line-height: 1.85;
}
.info-pill {
    display: inline-block;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: var(--radius-pill);
    padding: 3px 13px;
    margin: 3px 3px;
    font-family: var(--font-body);
    font-size: 0.76rem;
    color: var(--text-muted);
    transition: border-color 0.15s, color 0.15s;
    cursor: default;
}
.info-pill:hover { border-color: var(--accent); color: var(--text); }

/* ── Literature Intelligence ────────────────────────────────── */
.lit-header {
    font-family: var(--font-head);
    font-size: 1.2rem;
    font-weight: 800;
    color: var(--text);
    border-bottom: 1.5px solid var(--border);
    padding-bottom: 9px;
    margin: 1.8rem 0 0.4rem;
    letter-spacing: -0.2px;
}
.lit-sub {
    font-family: var(--font-body);
    font-size: 0.81rem;
    color: var(--text-muted);
    margin-bottom: 1rem;
    line-height: 1.6;
}
.lit-source {
    font-family: var(--font-body);
    font-size: 0.74rem;
    color: var(--success);
    font-weight: 600;
}
.lit-cat {
    font-family: var(--font-head);
    font-size: 0.92rem;
    font-weight: 700;
    padding-left: 9px;
    border-left: 3px solid var(--ac);
    color: var(--ac);
    margin: 1.3rem 0 0.4rem;
}
.lit-card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 12px 16px;
    margin-bottom: 8px;
    transition: border-color 0.2s, transform 0.15s;
}
.lit-card:hover {
    border-color: var(--border-hi);
    transform: translateX(3px);
}
.lit-title {
    font-family: var(--font-body);
    font-size: 0.89rem;
    font-weight: 600;
    color: var(--text);
    margin-bottom: 3px;
    line-height: 1.4;
}
.lit-meta {
    font-family: var(--font-body);
    font-size: 0.74rem;
    color: var(--text-muted);
}
.lit-abs {
    font-family: var(--font-body);
    font-size: 0.74rem;
    color: var(--text-muted);
    font-style: italic;
    margin-top: 5px;
    line-height: 1.5;
}
.lit-pill {
    display: inline-block;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 3px 12px;
    margin: 3px;
    font-family: var(--font-body);
    font-size: 0.74rem;
    color: var(--text-muted);
}
.lit-empty {
    background: var(--card);
    border: 1px dashed var(--border);
    border-radius: 10px;
    padding: 26px;
    text-align: center;
    font-family: var(--font-body);
    color: var(--text-muted);
    font-size: 0.88rem;
}

/* ── Scope / Warning Banners ────────────────────────────────── */
.scope-banner {
    background: linear-gradient(135deg,
        rgba(0,200,232,0.04) 0%,
        rgba(155,113,240,0.04) 100%);
    border: 1px solid var(--border);
    border-left: 4px solid var(--accent);
    border-radius: var(--radius);
    padding: 14px 18px;
    margin-bottom: 16px;
}
.scope-title {
    font-family: var(--font-head);
    font-size: 0.86rem;
    font-weight: 700;
    color: var(--accent);
    margin-bottom: 8px;
    letter-spacing: 0.3px;
}
.scope-body {
    font-family: var(--font-body);
    font-size: 0.81rem;
    color: var(--text-muted);
    line-height: 1.75;
}
.oos-banner {
    background: rgba(245, 166, 35, 0.04);
    border: 1px solid rgba(245, 166, 35, 0.5);
    border-left: 4px solid var(--warning);
    border-radius: var(--radius);
    padding: 14px 18px;
    margin: 14px 0;
}
.oos-title {
    font-family: var(--font-head);
    font-size: 0.86rem;
    font-weight: 700;
    color: var(--warning);
    margin-bottom: 8px;
}
.oos-item {
    font-family: var(--font-body);
    font-size: 0.81rem;
    color: var(--text);
    margin-bottom: 5px;
    line-height: 1.55;
}
.oos-footer {
    font-family: var(--font-body);
    font-size: 0.76rem;
    color: var(--text-muted);
    margin-top: 10px;
    line-height: 1.5;
}

/* ── Brand badge ─────────────────────────────────────────────── */
.brand-badge {
    display: inline-flex;
    align-items: center;
    gap: 8px;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: var(--radius-pill);
    padding: 4px 14px 4px 10px;
    margin-bottom: 14px;
    font-family: var(--font-body);
    font-size: 0.76rem;
    font-weight: 600;
    color: var(--text-muted);
    letter-spacing: 0.4px;
}
.brand-dot {
    width: 7px; height: 7px;
    border-radius: 50%;
    background: var(--success);
    box-shadow: 0 0 7px var(--success);
    display: inline-block;
    animation: blink 2.4s ease-in-out infinite;
}
@keyframes blink {
    0%, 100% { opacity: 1; }
    50%       { opacity: 0.3; }
}

/* ── Sidebar custom brand ────────────────────────────────────── */
.sb-brand {
    display: flex;
    align-items: center;
    gap: 10px;
    margin-bottom: 4px;
}
.sb-icon {
    width: 34px; height: 34px;
    border-radius: 9px;
    background: linear-gradient(135deg, #0891B2, #7C3AED);
    display: flex; align-items: center; justify-content: center;
    font-size: 1rem;
    flex-shrink: 0;
}
.sb-name {
    font-family: var(--font-head);
    font-size: 1.05rem;
    font-weight: 800;
    background: linear-gradient(90deg, var(--accent), var(--accent3));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    letter-spacing: 0.2px;
}
.sb-tagline {
    font-family: var(--font-body);
    font-size: 0.72rem;
    color: var(--text-dim);
    margin-top: 2px;
}
.sb-divider {
    border: none;
    border-top: 1px solid var(--border);
    margin: 14px 0;
}

</style>
""", unsafe_allow_html=True)

# ── Helpers ───────────────────────────────────────────────────────────────────
def fig_to_bytes(fig, dpi=300):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    buf.seek(0)
    return buf.read()

def risk_icon(val):
    return {
        "HIGH":     "🟠", "CRITICAL": "🔴", "MODERATE": "🟡", "LOW":      "🟢",
        "CLEAN":    "🟢", "POSITIVE": "🔴", "NEGATIVE": "🟢", "FLAGGED":  "🟠"
    }.get(str(val).upper(), "⚪")

def sh(t):  return f'<p class="sec-head">{t}</p>'
def sub(t): return f'<p class="sub-head">{t}</p>'

@st.cache_data(show_spinner=False)
def cached_analyze(smiles, fetch_pubchem):
    return analyze_smiles(smiles, fetch_pubchem=fetch_pubchem)


# ── GLOSSARY ──────────────────────────────────────────────────────────────────
GLOSSARY = {
    "Drug-likeness Score (0–100)": (
        "A composite heuristic score estimating how suitable a molecule is as an oral drug. "
        "≥ 70 = ACCEPT, 50–69 = borderline, < 50 = REJECT."
    ),
    "Lipinski Rule-of-Five": (
        "Oral bioavailability filter: MW ≤ 500 Da, LogP ≤ 5, H-bond donors ≤ 5, "
        "H-bond acceptors ≤ 10. Failing ≥ 2 rules predicts poor absorption."
    ),
    "LogP": (
        "Octanol–water partition coefficient — measure of lipophilicity. "
        "High LogP (> 5) → poor solubility, high metabolism. Low LogP (< 0) → poor membrane permeation."
    ),
    "TPSA (Topological Polar Surface Area)": (
        "Sum of polar atom surface areas in Å². Predicts gut absorption and BBB penetration. "
        "≤ 90 Å² → good oral absorption; ≤ 60 Å² → likely CNS penetration; > 140 Å² → poor oral bioavailability."
    ),
    "CYP450 (Cytochrome P450)": (
        "Family of liver enzymes responsible for metabolising ~75% of drugs. "
        "Key isoforms: CYP3A4, CYP2D6, CYP2C9, CYP2C19, CYP1A2. "
        "Structural alerts predict inhibition risk and potential drug–drug interactions (DDIs)."
    ),
    "hERG Liability": (
        "hERG is a cardiac ion channel. Inhibition prolongs the QT interval, "
        "risking fatal arrhythmia (torsades de pointes). High LogP + basic nitrogen → HIGH risk."
    ),
    "Ames Mutagenicity": (
        "Bacterial reverse-mutation assay predicting genotoxicity. "
        "POSITIVE indicates structural alerts (nitroaromatics, alkylating agents) associated with DNA damage."
    ),
    "Microsomal Stability / t½": (
        "Predicted in vitro half-life in liver microsomes. "
        "HIGH (t½ > 60 min) = metabolically stable. LOW (t½ < 30 min) = rapid clearance → short in vivo half-life."
    ),
    "Reactive Metabolites": (
        "Some functional groups (furans, thiophenes, catechols) are bioactivated to electrophilic species "
        "that covalently modify proteins → idiosyncratic toxicity. CRITICAL = strong concern."
    ),
    "GSH Trapping": (
        "Glutathione (GSH) trapping assay detects reactive intermediates in vitro. "
        "POSITIVE = molecule forms adducts with GSH, indicating bioactivation liability."
    ),
    "MBI (Mechanism-Based Inactivation)": (
        "Irreversible CYP inhibition caused by a reactive metabolite binding covalently to the enzyme. "
        "Results in non-linear DDIs and requires time-dependent inhibition (TDI) follow-up assays."
    ),
    "PAINS (Pan-Assay Interference Compounds)": (
        "Structural motifs (e.g. rhodanines, catechols, quinones) that produce false positives in HTS assays "
        "through non-specific mechanisms (aggregation, redox cycling, fluorescence). FLAGGED = assay artefact risk."
    ),
    "Brenk Alerts": (
        "Substructure filter from Brenk et al. (2008) identifying fragments with unfavourable "
        "physicochemical properties or known toxicity, used for lead-hopping and library design."
    ),
    "SureChEMBL / ICH S2(R1)": (
        "Genotoxic impurity alerts based on structural classes from regulatory guidelines (ICH S2(R1), REACH). "
        "HIGH = regulatory concern requiring dedicated genotoxicity battery before clinical trials."
    ),
    "PPB (Plasma Protein Binding)": (
        "Fraction of drug bound to plasma proteins (albumin, AGP). "
        "Only unbound (free) fraction is pharmacologically active. "
        "> 95% binding → reduced free drug; may affect efficacy and duration."
    ),
    "Caco-2 Permeability (Papp)": (
        "Human colon carcinoma cell monolayer assay modelling intestinal absorption. "
        "Papp > 20 nm/s = HIGH (good oral absorption); < 5 nm/s = LOW (poor absorption)."
    ),
    "P-gp (P-glycoprotein) Efflux": (
        "ABC transporter (MDR1/ABCB1) that pumps substrates out of cells — intestinal epithelium, BBB, tumour cells. "
        "HIGH P-gp risk → reduced oral absorption and CNS penetration; potential multidrug resistance."
    ),
    "Metabolic Soft Spots": (
        "Structural positions predicted to undergo preferential oxidative metabolism (CYP-mediated). "
        "Blocking these sites (fluorination, deuteration, steric shielding) can improve metabolic stability."
    ),
}


def render_glossary():
    st.markdown(sh("📖 Term Glossary"), unsafe_allow_html=True)
    st.caption("Every heading used in this report is explained here. Expand to read.")
    with st.expander("📖 Click to read the full glossary", expanded=False):
        terms = list(GLOSSARY.items())
        mid   = (len(terms) + 1) // 2
        c1, c2 = st.columns(2)
        for col, chunk in [(c1, terms[:mid]), (c2, terms[mid:])]:
            with col:
                for term, defn in chunk:
                    st.markdown(
                        f'<span class="gloss-term">{term}</span><br>'
                        f'<span class="gloss-def">{defn}</span>',
                        unsafe_allow_html=True,
                    )
                    st.markdown("---")


# ── PubChem block ─────────────────────────────────────────────────────────────
def render_pubchem(pc):
    if not pc.get("found"):
        st.warning(f"⚠️ {pc.get('error', 'Compound not found in PubChem database.')}")
        return

    cid  = pc["cid"]
    name = pc.get("common_name") or pc.get("iupac_name") or f"CID {cid}"

    img_col, id_col = st.columns([1, 3])
    with img_col:
        st.image(pc["links"]["structure_image"], use_container_width=True,
                 caption=f"CID {cid}")
    with id_col:
        st.markdown(f'<div class="pubchem-name">{name}</div>', unsafe_allow_html=True)
        meta = (
            f'Formula: <b>{pc.get("formula","")}</b> &nbsp;|&nbsp; '
            f'MW: <b>{pc.get("molecular_weight","")}&thinsp;Da</b> &nbsp;|&nbsp; '
            f'Charge: <b>{pc.get("charge","")}</b>'
        )
        if pc.get("xlogp") is not None:
            meta += f' &nbsp;|&nbsp; XLogP: <b>{pc["xlogp"]}</b>'
        if pc.get("exact_mass"):
            meta += f' &nbsp;|&nbsp; Exact mass: <b>{pc["exact_mass"]}&thinsp;Da</b>'
        st.markdown(f'<div class="pubchem-meta">{meta}</div>', unsafe_allow_html=True)
        st.markdown(
            f'🔗 **[Open on PubChem]({pc["links"]["compound_page"]})**'
            f' &nbsp;|&nbsp; [SDF]({pc["links"]["sdf_download"]})'
            f' &nbsp;|&nbsp; [JSON]({pc["links"]["json_api"]})'
        )
        if pc.get("synonyms"):
            pills = "".join(f'<span class="info-pill">{s}</span>'
                            for s in pc["synonyms"][:8])
            st.markdown(f"**Known names:** {pills}", unsafe_allow_html=True)

    st.markdown("---")

    with st.expander("🔑 Chemical Identifiers", expanded=True):
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"**CID:** `{cid}`")
            st.markdown("**IUPAC name:**"); st.code(pc.get("iupac_name",""), language=None)
            st.markdown("**Canonical SMILES:**"); st.code(pc.get("canonical_smiles",""), language=None)
        with c2:
            st.markdown(f"**InChIKey:** `{pc.get('inchikey','')}`")
            if pc.get("inchi"):
                st.markdown("**InChI:**"); st.code(pc["inchi"], language=None)


# ── Detail sections ───────────────────────────────────────────────────────────
def render_detail(r):
    props = r["properties"]

    with st.expander("⚗️ Physicochemical Properties", expanded=False):
        ca, cb = st.columns(2)
        with ca:
            st.markdown(sub("Lipinski Rule-of-Five"), unsafe_allow_html=True)
            for pname, pval, plim in [
                ("Molecular Weight", props["molecular_weight"], 500),
                ("LogP",             props["logp"],             5),
                ("H-bond Donors",    props["hbd"],              5),
                ("H-bond Acceptors", props["hba"],              10),
            ]:
                st.markdown(f"{'✅' if pval<=plim else '❌'} **{pname}:** `{pval:.2f}` *(≤ {plim})*")
        with cb:
            st.markdown(sub("Extended"), unsafe_allow_html=True)
            st.markdown(f"{'✅' if props['tpsa']<=140 else '❌'} **TPSA:** `{props['tpsa']:.1f} Å²` *(≤ 140)*")
            st.markdown(f"{'✅' if props['rotatable_bonds']<=10 else '❌'} **Rot. Bonds:** `{props['rotatable_bonds']}` *(≤ 10)*")
        if r["violations"]:
            st.warning(f"⚠️ **{len(r['violations'])} violation(s)**")
            for v in r["violations"]:
                st.markdown(f"- **{v['property']}** = {v['value']:.2f} (limit {v['limit']}) — {v['issue']}")

    with st.expander("🔬 Metabolic Stability", expanded=False):
        ms, rm, gs, mb, ss = (
            r["microsomal_stability"], r["reactive_metabolites"],
            r["gsh_trapping"],        r["mbi_risk"],
            r["metabolic_soft_spots"],
        )
        mc1, mc2 = st.columns(2)
        with mc1:
            cls  = ms.get("stability_class","N/A")
            icon = {"HIGH":"🟢","MODERATE":"🟡","LOW":"🔴"}.get(cls,"⚪")
            st.markdown(sub("Microsomal Stability"), unsafe_allow_html=True)
            st.markdown(f"{icon} **{cls}** — Score `{ms.get('microsomal_stability_score')}/100`, "
                        f"t½ `{ms.get('predicted_t_half')}`")
            st.markdown(sub("Metabolic Soft Spots"), unsafe_allow_html=True)
            if ss["soft_spots"]:
                for s in ss["soft_spots"]: st.markdown(f"- {s}")
                st.caption(f"Clearance: **{ss.get('clearance_prediction')}**")
            else:
                st.success("None detected.")
        with mc2:
            rm_risk = rm.get("reactive_metabolite_risk","N/A")
            st.markdown(sub("Reactive Metabolites"), unsafe_allow_html=True)
            st.markdown(f"{risk_icon(rm_risk)} **{rm_risk}** — {rm.get('alert_count',0)} alert(s)")
            for a in rm.get("reactive_metabolite_alerts",[]): st.markdown(f"- {a}")
            for a in rm.get("recommended_assays",[]): st.markdown(f"  ↳ {a}")
            st.markdown(sub("GSH Trapping"), unsafe_allow_html=True)
            gsh_r = gs.get("gsh_trapping_risk","N/A")
            st.markdown(f"{risk_icon(gsh_r)} **{gsh_r}**")
            for a in gs.get("gsh_trap_alerts",[]): st.markdown(f"- {a}")
            st.markdown(sub("MBI Risk"), unsafe_allow_html=True)
            mbi_r = mb.get("mbi_risk","N/A")
            st.markdown(f"{risk_icon(mbi_r)} **{mbi_r}**")
            for a in mb.get("mbi_alerts",[]): st.markdown(f"- {a}")

    with st.expander("☠️ Toxicophore Alerts", expanded=False):
        summary = r["toxicophore_summary"]
        overall = summary.get("overall_toxicophore_risk","N/A")
        oi = {"CRITICAL":"🔴","HIGH":"🟠","MODERATE":"🟡","CLEAN":"🟢"}.get(overall,"⚪")
        st.markdown(f"### {oi} Overall: **{overall}**")
        st.info(summary.get("filter_recommendation",""))
        bd = summary.get("breakdown",{})
        tc1, tc2, tc3, tc4, tc5 = st.columns(5)
        for col, (lbl, key) in zip(
            [tc1,tc2,tc3,tc4,tc5],
            [("PAINS","pains"),("Brenk","brenk"),("SureChEMBL","surechembl"),
             ("Ames","ames"),("Reactive","reactive")]
        ):
            cnt = bd.get(key,0)
            col.metric(lbl, f"{'🔴' if cnt>=3 else '🟠' if cnt>=1 else '🟢'} {cnt}")
        tc_a, tc_b, tc_c = st.columns(3)
        for col, key, lbl in [(tc_a,"pains","PAINS"),(tc_b,"brenk","Brenk"),(tc_c,"surechembl","SureChEMBL")]:
            with col:
                alerts = r[key].get(f"{key}_alerts",[])
                flag   = r[key].get(f"{key}_flag", r[key].get("regulatory_risk",""))
                st.markdown(f"**{lbl} — {flag}**")
                if alerts:
                    for a in alerts: st.markdown(f"- {a}")
                else:
                    st.success(f"No {lbl} alerts.")

    with st.expander("🧬 ADME Properties", expanded=False):
        ac1, ac2 = st.columns(2)
        with ac1:
            cyp  = r["cyp450"]
            herg = r["herg"]
            ames = r["mutagenicity"]
            st.markdown(sub("CYP450"), unsafe_allow_html=True)
            st.markdown(f"{risk_icon(cyp['cyp_risk'])} **{cyp['cyp_risk']}** ({cyp['isoforms_flagged_count']} isoform(s))")
            for iso, alerts in cyp.get("flagged_isoforms",{}).items():
                st.markdown(f"  *{iso}:* {len(alerts)} alert(s)")
            st.markdown(sub("hERG Cardiac Liability"), unsafe_allow_html=True)
            st.markdown(f"{risk_icon(herg['herg_risk'])} **{herg['herg_risk']}** (score {herg['herg_risk_score']})")
            st.markdown(sub("Ames Mutagenicity"), unsafe_allow_html=True)
            st.markdown(f"{risk_icon(ames['mutagenicity_flag'])} **{ames['mutagenicity_flag']}** — {ames['alert_count']} alert(s)")
            for a in ames.get("alerts",[]): st.markdown(f"- {a}")
        with ac2:
            ppb   = r["plasma_protein_binding"]
            caco2 = r["caco2_permeability"]
            pgp   = r["pgp_efflux"]
            st.markdown(sub("Plasma Protein Binding (PPB)"), unsafe_allow_html=True)
            st.markdown(f"{risk_icon(ppb['ppb_class'])} **{ppb['ppb_class']}** — ~{ppb['ppb_estimate_pct']}% bound")
            st.caption(ppb.get("ppb_note",""))
            st.markdown(sub("Caco-2 Permeability"), unsafe_allow_html=True)
            ci = {"HIGH":"🟢","MEDIUM":"🟡","LOW":"🔴"}.get(caco2["caco2_class"],"⚪")
            st.markdown(f"{ci} **{caco2['caco2_class']}** — {caco2['caco2_papp_nm_s']} nm/s")
            st.caption(caco2.get("caco2_note",""))
            st.markdown(sub("P-glycoprotein Efflux"), unsafe_allow_html=True)
            st.markdown(f"{risk_icon(pgp['pgp_risk'])} **{pgp['pgp_risk']}** (score {pgp['pgp_risk_score']})")
            for reason in pgp.get("pgp_risk_reasons",[]): st.markdown(f"- {reason}")

    with st.expander("💡 Optimisation Advice", expanded=False):
        adv = generate_optimization_advice(r)
        if not adv:
            st.success("✅ No optimisation needed — all filters passed.")
        else:
            for a in sorted(adv):
                st.markdown(f"{'⚠️' if a.startswith('PRIORITY') else '▸'} {a}")


# ── Panel display helper ──────────────────────────────────────────────────────
PANEL_META = {
    "overview":    ("📊 Overview",                  "Drug-likeness score and physicochemical radar chart"),
    "risk_matrix": ("🚦 ADME Risk Matrix",           "Full 14-row risk rating across all ADME endpoints"),
    "alerts":      ("🔎 Alerts & Stability",         "Structural alert counts + microsomal stability breakdown"),
    "cyp_tox":     ("⚠️ CYP450 & Toxicophores",      "CYP isoform liability and toxicophore alert distribution"),
    "properties":  ("📋 Property Table",             "Numerical properties vs. drug-likeness thresholds"),
    "advice":      ("💡 Optimisation Recommendations","Priority chemistry modifications suggested"),
}

def render_panels(result, key_prefix=""):
    panels = build_adme_panels(result)
    if not panels:
        st.error("Could not build charts.")
        return

    panel_order = ["overview", "risk_matrix", "alerts", "cyp_tox", "properties", "advice"]

    fig_combined = build_adme_figure(result)
    if fig_combined:
        combined_bytes = fig_to_bytes(fig_combined, dpi=300)
        plt.close(fig_combined)
        st.download_button(
            label="⬇️ Download combined report (300 DPI)",
            data=combined_bytes,
            file_name=f"adme_combined_{key_prefix or 'report'}.png",
            mime="image/png",
            key=f"dl_combined_{key_prefix}",
        )

    st.markdown("")
    for key in panel_order:
        fig = panels.get(key)
        if fig is None:
            continue
        title, caption = PANEL_META[key]
        st.markdown(sh(title), unsafe_allow_html=True)
        st.caption(caption)
        st.pyplot(fig, use_container_width=True)
        png_bytes = fig_to_bytes(fig, dpi=300)
        plt.close(fig)
        st.download_button(
            label=f"⬇️ Download — {title} (300 DPI)",
            data=png_bytes,
            file_name=f"adme_{key}_{key_prefix or 'report'}.png",
            mime="image/png",
            key=f"dl_{key}_{key_prefix}",
        )
        st.markdown("")


# ══════════════════════════════════════════════════════════════════════════════
# LITERATURE INTELLIGENCE
# ══════════════════════════════════════════════════════════════════════════════

_LIT_CATS = {
    "Metabolism":   ("⚗️",  "#00C8E8"),
    "Toxicity":     ("☠️",  "#F05252"),
    "Resistance":   ("🛡️", "#F5A623"),
    "Pharmacology": ("💊",  "#22C983"),
    "Clinical":     ("🏥",  "#9B71F0"),
    "Synthesis":    ("🔬",  "#F97316"),
    "General":      ("📄",  "#7A8FA8"),
}
_LIT_ORDER = ["Metabolism","Toxicity","Resistance","Pharmacology","Clinical","Synthesis","General"]


def _lit_paper_card(paper: dict, accent: str) -> None:
    parts = []
    if paper.get("year"):    parts.append(f"📅 {paper['year']}")
    if paper.get("journal"): parts.append(f"📰 {paper['journal']}")
    if paper.get("doi"):     parts.append(f"DOI: {paper['doi']}")
    meta = "  ·  ".join(parts)
    url  = paper.get("url", "")
    link = (f'<a href="{url}" target="_blank" style="color:{accent};'
            f'font-size:0.77rem;text-decoration:none;font-weight:600;">🔗 View</a>') if url else ""
    abs_html = (f'<div class="lit-abs">"{paper["abstract"]}"</div>'
                if paper.get("abstract") else "")
    st.markdown(
        f'<div class="lit-card">'
        f'<div class="lit-title">{paper.get("title","Untitled")}</div>'
        f'<div class="lit-meta">{meta} &nbsp; {link}</div>'
        + abs_html + "</div>",
        unsafe_allow_html=True,
    )


def render_literature_intelligence(pubchem_metadata: dict, key_prefix: str = "") -> None:
    st.markdown(
        '<div class="lit-header">📚 Literature Intelligence</div>'
        '<div class="lit-sub">Recent papers retrieved automatically '
        '(PubMed → Semantic Scholar → OpenAlex → CrossRef) '
        '— grouped by research area.</div>',
        unsafe_allow_html=True,
    )

    pc_name = ""
    if pubchem_metadata and pubchem_metadata.get("found"):
        pc_name = (pubchem_metadata.get("common_name") or
                   pubchem_metadata.get("iupac_name") or "").strip()

    c1, c2 = st.columns([4, 1])
    with c1:
        query = st.text_input(
            "lit_input", value=pc_name,
            placeholder="e.g. lidocaine, aspirin, ibuprofen",
            key=f"lit_q_{key_prefix}", label_visibility="collapsed",
        )
    with c2:
        clicked = st.button("🔍 Search", key=f"lit_btn_{key_prefix}",
                            type="primary", use_container_width=True)

    if not pc_name and not clicked:
        st.caption("💡 Enable Fetch PubChem in sidebar for auto-fill, or type a name and press Search.")

    auto_k = f"lit_done_{key_prefix}"
    if clicked:
        st.session_state[auto_k] = False
        should = True
    elif pc_name and not st.session_state.get(auto_k, False):
        should = True
    else:
        should = False

    if not should or not query.strip():
        return

    if not clicked:
        st.session_state[auto_k] = True

    with st.spinner(f"Searching for **{query.strip()}** across PubMed, Semantic Scholar, OpenAlex…"):
        res = get_literature_intelligence(query.strip())

    if res.get("error"):
        st.error(f"⚠️ {res['error']}")
        return

    if res["total"] == 0:
        reason = res.get("_empty_reason", "")
        st.markdown(
            '<div class="lit-empty">📭 No papers found for this compound.'
            + (f"<br><small style='color:#7A8FA8'>{reason}</small>" if reason else "")
            + "</div>",
            unsafe_allow_html=True,
        )
        return

    cats  = res["categories"]
    total = res["total"]
    src   = res.get("source", "")

    pills = f'<span class="lit-pill">🗂️ {total} paper{"s" if total!=1 else ""}</span>'
    for cat, papers in cats.items():
        em, _ = _LIT_CATS.get(cat, ("📄", "#7A8FA8"))
        pills += f'<span class="lit-pill">{em} {cat}: {len(papers)}</span>'
    if src and src != "none":
        pills += f'<span class="lit-source" style="margin-left:8px;">via {src}</span>'
    st.markdown(pills, unsafe_allow_html=True)
    st.markdown("")

    ordered = [c for c in _LIT_ORDER if c in cats] + [c for c in cats if c not in _LIT_ORDER]
    for cat in ordered:
        papers = cats.get(cat, [])
        if not papers:
            continue
        em, ac = _LIT_CATS.get(cat, ("📄", "#7A8FA8"))
        n = len(papers)
        st.markdown(
            f'<div class="lit-cat" style="--ac:{ac};">{em} {cat} '
            f'<span style="font-size:0.78rem;font-weight:400;color:#7A8FA8;">'
            f'({n} paper{"s" if n!=1 else ""})</span></div>',
            unsafe_allow_html=True,
        )
        for p in papers:
            _lit_paper_card(p, ac)


# ══════════════════════════════════════════════════════════════════════════════
# SCOPE BANNER + OUT-OF-SCOPE VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

def _render_scope_banner() -> None:
    st.markdown("""
<div class="scope-banner">
  <div class="scope-title">🎯 Designed for drug-like organic small molecules</div>
  <div class="scope-body">
    <b style="color:#EDF4FF;">Works well with:</b> oral drug candidates, investigational compounds,
    natural product derivatives, synthetic intermediates with MW 150–800 Da.<br>
    <b style="color:#EDF4FF;">Limited accuracy for:</b> heavy metals, inorganic salts, simple
    industrial solvents, polymers, biologics, or any molecule with fewer than 5 heavy atoms.<br>
    <b style="color:#EDF4FF;">All outputs are rule-based structural predictions</b> —
    not a substitute for experimental ADME data.
  </div>
</div>
""", unsafe_allow_html=True)


def _out_of_scope_warning(result: dict) -> None:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Descriptors

    mol = result.get("mol")
    if mol is None:
        return

    props    = result.get("properties", {})
    mw       = props.get("molecular_weight", 0)
    warnings = []

    HEAVY_METALS = {82, 80, 33, 48, 26, 29, 30, 28, 27, 24, 25, 34, 81}
    metal_symbols = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in HEAVY_METALS:
            metal_symbols.append(atom.GetSymbol())
    if metal_symbols:
        unique = list(dict.fromkeys(metal_symbols))
        warnings.append(
            f"**Heavy metal atom(s) detected ({', '.join(unique)})** — "
            "ADME models are designed for organic molecules. CYP, hERG, PPB and "
            "microsomal stability outputs are not meaningful for metal-containing compounds."
        )

    if mw < 150:
        warnings.append(
            f"**Molecular weight is very low ({mw:.0f} Da)** — "
            "molecules below 150 Da are typically solvents, reagents or fragments, "
            "not drug candidates. Lipinski and permeability scores are unreliable at this size."
        )

    if mw > 900:
        warnings.append(
            f"**Molecular weight is very high ({mw:.0f} Da)** — "
            "oral small-molecule drug space is typically 150–600 Da. "
            "Above 900 Da, Lipinski rules and Caco-2 predictions lose validity."
        )

    carbon_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    if carbon_count == 0:
        warnings.append(
            "**No carbon atoms detected** — this appears to be an inorganic compound. "
            "All ADME modules assume organic chemistry. Results are not meaningful."
        )

    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count < 5:
        warnings.append(
            f"**Very few heavy atoms ({heavy_atom_count})** — "
            "this molecule is too simple to produce meaningful drug-likeness predictions."
        )

    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count == 0 and mw < 200:
        warnings.append(
            "**No ring systems and low molecular weight** — "
            "this looks like a simple acyclic compound (solvent, reagent or industrial chemical). "
            "ADME predictions are designed for more complex drug-like scaffolds."
        )

    if not warnings:
        return

    items_html = "".join(
        f'<div class="oos-item">▸ {w}</div>' for w in warnings
    )
    st.markdown(
        f'<div class="oos-banner">'
        f'<div class="oos-title">⚠️ Out-of-scope molecule detected — interpret results with caution</div>'
        + items_html +
        f'<div class="oos-footer">'
        f'The full analysis is still shown below — but the values above are not pharmacologically '
        f'meaningful for this molecule type. LeadRefine is designed for drug-like organic '
        f'small molecules (MW 150–900 Da, carbon-containing, structurally complex).'
        f'</div></div>',
        unsafe_allow_html=True,
    )


# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("""
<div class="sb-brand">
  <div class="sb-icon">🧬</div>
  <div>
    <div class="sb-name">LeadRefine</div>
    <div class="sb-tagline">ADME Analysis Platform</div>
  </div>
</div>
<hr class="sb-divider">
""", unsafe_allow_html=True)

    mode = st.radio("Mode", ["Single molecule", "Batch (multiple SMILES)"])

    st.markdown('<hr class="sb-divider">', unsafe_allow_html=True)

    fetch_pubchem = st.toggle(
        "🔗 Fetch PubChem data",
        value=True,
        help="Queries PubChem for compound name, structure, identifiers and synonyms (~3–5 s). "
             "Turn OFF for instant results.",
    )
    if fetch_pubchem:
        st.caption("~3–5 s per molecule (name + structure only).")
    else:
        st.info("⚡ Fast mode — PubChem section skipped.")

    st.markdown('<hr class="sb-divider">', unsafe_allow_html=True)
    st.markdown("**Demo molecules**")
    demos = {
        "Aspirin":      "CC(=O)Oc1ccccc1C(=O)O",
        "CCl4":         "ClC(Cl)(Cl)Cl",
        "Ibuprofen":    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        "Caffeine":     "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        "Morphine":     "OC1=CC=C2CC3N(CCC34CCc5c4cc(O)c(c5)OC)C2=C1",
    }
    selected_demo = st.selectbox("Load a demo", ["— select —"] + list(demos.keys()))
    st.markdown('<hr class="sb-divider">', unsafe_allow_html=True)
    st.caption("Rule-based screening.\nNot a substitute for experimental data.")


# ── Title + Glossary (always shown) ──────────────────────────────────────────
st.markdown(
    '<div class="brand-badge"><span class="brand-dot"></span>Rule-based · Structural Prediction</div>',
    unsafe_allow_html=True,
)
st.title("🧬 ADME Analysis Dashboard")
st.markdown(
    "Screen molecules for drug-likeness, metabolic stability, toxicophore alerts, "
    "and ADME properties — enriched with live PubChem data."
)
render_glossary()

# ══════════════════════════════════════════════════════════════════
# SINGLE MOLECULE
# ══════════════════════════════════════════════════════════════════
if mode == "Single molecule":

    _render_scope_banner()
    default_smiles = demos[selected_demo] if selected_demo != "— select —" else ""
    smiles_input = st.text_input(
        "SMILES string",
        value=default_smiles,
        placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O",
    )
    if fetch_pubchem:
        st.caption("⏱️ With PubChem: ~5–10 s  |  Toggle OFF in sidebar for < 1 s")

    if st.button("Analyse", type="primary", use_container_width=True) and smiles_input.strip():
        _t0 = time.time()
        _st = st.empty()
        _st.info("🔬 Running ADME analysis…")
        _res = cached_analyze(smiles_input.strip(), fetch_pubchem)
        _st.success(f"✅ Done in {time.time()-_t0:.1f} s")
        if not _res["valid"]:
            st.error("❌ Invalid SMILES string.")
            st.stop()
        st.session_state["single_result"] = _res
        st.session_state["lit_done_single"] = False

    # ── Render results ───────────────────────────────────────────────────────
    result = st.session_state.get("single_result")
    if result and result.get("valid"):

        _out_of_scope_warning(result)

        st.markdown("---")
        k1, k2, k3, k4, k5 = st.columns(5)
        d_icon = "🟢" if result["decision"]=="ACCEPT" else "🔴"
        k1.metric("Score",       f"{result['score']} / 100")
        k2.metric("Status",      result["physchem_status"])
        k3.metric("Decision",    f"{d_icon} {result['decision']}")
        k4.metric("Lipinski",    "✅ PASS" if result["lipinski"]["passes"] else "❌ FAIL")
        k5.metric("Tox Alerts",  result["toxicophore_summary"]["total_toxicophore_alerts"])

        st.markdown(sh("📊 ADME Charts"), unsafe_allow_html=True)
        render_panels(result, key_prefix="single")

        st.markdown(sh("🔗 PubChem Information"), unsafe_allow_html=True)
        render_pubchem(result.get("pubchem_metadata", {}))

        st.markdown(sh("🔬 Detailed Analysis"), unsafe_allow_html=True)
        st.caption("Click any section to expand.")
        render_detail(result)

        st.markdown("---")
        render_literature_intelligence(
            pubchem_metadata=result.get("pubchem_metadata", {}),
            key_prefix="single",
        )


# ══════════════════════════════════════════════════════════════════
# BATCH MODE
# ══════════════════════════════════════════════════════════════════
else:
    _render_scope_banner()
    st.markdown("Enter one SMILES per line:")
    batch_input = st.text_area(
        "SMILES list", height=180,
        placeholder="CC(=O)Oc1ccccc1C(=O)O\nClC(Cl)(Cl)Cl\nCC(C)Cc1ccc(cc1)C(C)C(=O)O",
    )
    if fetch_pubchem:
        st.caption("⏱️ PubChem ON: ~5–10 s per molecule. Toggle OFF in sidebar for instant results.")
    else:
        st.caption("⚡ PubChem OFF: < 1 s per molecule.")

    if st.button("Analyse all", type="primary", use_container_width=True) and batch_input.strip():
        smiles_list = [s.strip() for s in batch_input.strip().splitlines() if s.strip()]
        n = len(smiles_list)

        est = n * (7 if fetch_pubchem else 1)
        est_str = f"{est//60} min {est%60} s" if est >= 60 else f"~{est} s"
        st.info(f"⏱️ Analysing **{n} molecule(s)** — estimated time: **{est_str}**")

        progress_bar = st.progress(0)
        status_text  = st.empty()
        results      = []
        t_start      = time.time()

        for idx, smi in enumerate(smiles_list):
            status_text.markdown(f"🔬 Molecule **{idx+1}/{n}**: `{smi[:50]}`")
            r = analyze_smiles(smi.strip(), fetch_pubchem=fetch_pubchem)
            results.append(r)
            elapsed = time.time() - t_start
            if idx > 0:
                rem = (elapsed/(idx+1)) * (n-idx-1)
                rem_str = f"{int(rem//60)} min {int(rem%60)} s" if rem>=60 else f"{int(rem)} s"
                status_text.markdown(f"✅ **{idx+1}/{n}** done  |  ⏳ ~{rem_str} remaining")
            progress_bar.progress((idx+1)/n)

        total = time.time() - t_start
        progress_bar.progress(1.0)
        status_text.success(f"✅ All {n} done in {total:.1f} s ({total/n:.1f} s/molecule avg)")

        valid   = [r for r in results if r["valid"]]
        invalid = [r for r in results if not r["valid"]]
        st.markdown(f"**{len(valid)} valid** | **{len(invalid)} invalid**")
        if invalid:
            with st.expander("❌ Invalid SMILES"):
                for r in invalid: st.markdown(f"- `{r['smiles']}`")
        if not valid: st.stop()

        st.markdown("---")
        st.markdown(sh("📋 Results Summary"), unsafe_allow_html=True)
        import pandas as pd
        rows = []
        for r in valid:
            pc      = r.get("pubchem_metadata",{})
            pc_name = (pc.get("common_name") or pc.get("iupac_name","")) if pc.get("found") else ""
            rows.append({
                "Name":            pc_name[:30] if pc_name else r["smiles"][:30],
                "Score":           r["score"],
                "Decision":        r["decision"],
                "Lipinski":        "PASS" if r["lipinski"]["passes"] else "FAIL",
                "Reactive Met.":   r["reactive_metabolites"]["reactive_metabolite_risk"],
                "Tox Alerts":      r["toxicophore_summary"]["total_toxicophore_alerts"],
                "hERG":            r["herg"]["herg_risk"],
                "CYP Risk":        r["cyp450"]["cyp_risk"],
                "Microsomal":      r["microsomal_stability"]["stability_class"],
            })
        st.dataframe(pd.DataFrame(rows), use_container_width=True)

        st.markdown("---")
        st.markdown(sh("📊 Individual Reports"), unsafe_allow_html=True)
        for i, r in enumerate(valid):
            pc   = r.get("pubchem_metadata",{})
            name = (pc.get("common_name") or pc.get("iupac_name","")) if pc.get("found") else ""
            label = f"Molecule {i+1}"
            if name: label += f" — {name}"
            label += f":  `{r['smiles'][:50]}`"

            with st.expander(label, expanded=False):
                _out_of_scope_warning(r)
                d_icon = "🟢" if r["decision"]=="ACCEPT" else "🔴"
                sc1, sc2, sc3, sc4 = st.columns(4)
                sc1.metric("Score",     f"{r['score']} / 100")
                sc2.metric("Decision",  f"{d_icon} {r['decision']}")
                sc3.metric("Lipinski",  "✅ PASS" if r["lipinski"]["passes"] else "❌ FAIL")
                sc4.metric("Tox Alerts",r["toxicophore_summary"]["total_toxicophore_alerts"])
                st.markdown("---")

                st.markdown(sub("📊 ADME Charts"), unsafe_allow_html=True)
                render_panels(r, key_prefix=f"mol{i}")

                st.markdown("---")
                st.markdown(sub("🔗 PubChem Information"), unsafe_allow_html=True)
                render_pubchem(pc)

                st.markdown("---")
                st.markdown(sub("🔬 Detailed Analysis"), unsafe_allow_html=True)
                render_detail(r)

                st.markdown("---")
                render_literature_intelligence(
                    pubchem_metadata=r.get("pubchem_metadata", {}),
                    key_prefix=f"mol{i}",
                )