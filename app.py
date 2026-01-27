# app.py

import streamlit as st
from rdkit.Chem import Draw
from analysis import (
    analyze_multiple_smiles,
    rank_molecules,
    dataset_decision,
    generate_optimization_advice,
    explain_imperfection
)

# -----------------------------
# Page config
# -----------------------------
st.set_page_config(
    page_title="Drug-Likeness Analyzer",
    layout="wide"
)

st.title("🧪 Drug-Likeness & Physicochemical Analysis")
st.write(
    "Enter **multiple SMILES codes** (one per line). "
    "Molecules are evaluated using physicochemical drug-likeness rules "
    "(not ADMET prediction)."
)

# -----------------------------
# Input
# -----------------------------
smiles_input = st.text_area(
    "Enter SMILES codes (one per line)",
    height=220,
    placeholder="CC(=O)Oc1ccccc1C(=O)O\nCC(C)Cc1ccc(cc1)C(=O)O"
)

analyze_btn = st.button("Analyze molecules")

# -----------------------------
# Processing
# -----------------------------
if analyze_btn:
    smiles_list = [s.strip() for s in smiles_input.splitlines() if s.strip()]

    if not smiles_list:
        st.warning("Please enter at least one SMILES string.")
    else:
        results = rank_molecules(analyze_multiple_smiles(smiles_list))

        # -----------------------------
        # Dataset decision
        # -----------------------------
        st.subheader("📊 Dataset Summary")
        st.info(dataset_decision(results))

        if not results:
            st.stop()

        top_score = results[0]["score"]

        # -----------------------------
        # Per-molecule analysis
        # -----------------------------
        st.subheader("🧬 Molecule-wise Analysis")

        for i, r in enumerate(results, 1):
            with st.expander(f"Rank {i} — {r['smiles']}"):

                if not r["valid"]:
                    st.error("Invalid SMILES string")
                    continue

                # ---- layout ----
                col_img, col_info = st.columns([1, 2])

                # ---- molecule image ----
                with col_img:
                    st.image(
                        Draw.MolToImage(r["mol"], size=(260, 260)),
                        caption="Molecular structure"
                    )

                # ---- summary ----
                with col_info:
                    st.markdown("### Decision Summary")
                    st.write(f"**Decision:** {r['decision']}")
                    st.write(f"**Physicochemical status:** {r['physchem_status']}")
                    st.write(f"**Drug-likeness score:** {r['score']}")

                    st.markdown("### Properties")
                    for k, v in r["properties"].items():
                        st.write(f"- **{k}**: {v}")

                # -----------------------------
                # Violations
                # -----------------------------
                if r["violations"]:
                    st.markdown("### ⚠ Rule Violations")
                    for v in r["violations"]:
                        st.write(
                            f"- **{v['property']}** "
                            f"(value: {v['value']} > limit: {v['limit']})"
                        )
                        st.write(f"  - Issue: {v['issue']}")
                else:
                    st.success("No physicochemical rule violations detected.")

                # -----------------------------
                # Rejection reason
                # -----------------------------
                if r["decision"] == "REJECT":
                    st.markdown("### ❌ Rejection Reason")
                    st.error(r["reject_reason"])
                    continue

                # -----------------------------
                # Rank-aware explanation
                # -----------------------------
                is_top = (i == 1)
                tied = (r["score"] == top_score and not is_top)

                if is_top:
                    st.markdown("### 🏆 Best Candidate Explanation")
                elif tied:
                    st.markdown("### 🧪 Equally Optimal Candidate")
                else:
                    st.markdown("### 🔬 Why this ranks below the top")

                for reason in explain_imperfection(
                    r,
                    is_top=is_top,
                    tied=tied
                ):
                    st.write(f"- {reason}")

                # -----------------------------
                # Optimization guidance
                # -----------------------------
                advice = generate_optimization_advice(r)
                if advice:
                    st.markdown("### 🔧 Optimization Suggestions")
                    for a in advice:
                        st.write(f"- {a}")

        st.success("Analysis complete ✅")
