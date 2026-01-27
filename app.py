import streamlit as st
from rdkit.Chem import Draw
from analysis import (
    analyze_multiple_smiles,
    rank_molecules,
    dataset_decision,
    generate_optimization_advice
)

st.set_page_config(page_title="LeadRefine", layout="wide")

st.title("🧪 LeadRefine — Drug-Likeness Screening Tool")

smiles_input = st.text_area(
    "Enter SMILES (one per line):",
    height=150,
    placeholder="CC(=O)Oc1ccccc1C(=O)O\nCC(C)Cc1ccc(cc1)C(=O)O"
)

if st.button("Run Analysis"):
    smiles_list = [s.strip() for s in smiles_input.splitlines() if s.strip()]
    results = analyze_multiple_smiles(smiles_list)
    ranked = rank_molecules(results)
    decision = dataset_decision(ranked)

    st.subheader(f"📊 Dataset Decision: {decision}")

    for i, r in enumerate(ranked, 1):
        st.markdown("---")
        st.subheader(f"Rank {i}: {r['smiles']}")
        st.write(f"**Physchem status:** {r['physchem_status']}")
        st.write(f"**Decision:** {r['decision']}")
        st.write(f"**Score:** {r['score']}")

        st.image(Draw.MolToImage(r["mol"], size=(250, 250)))

        st.json(r["properties"])

        if r["violations"]:
            st.warning("Rule Violations")
            for v in r["violations"]:
                st.write(f"- {v['issue']}")
        else:
            st.success("No rule violations")

        if r["physchem_status"] == "NEAR_OPTIMAL":
            advice = generate_optimization_advice(r)
            if advice:
                st.info("Optimization Suggestions:")
                for a in advice:
                    st.write(f"- {a}")

        if r["decision"] == "REJECT":
            st.error(f"Rejected: {r['reject_reason']}")
