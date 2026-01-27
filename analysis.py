from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# ------------------ RULES ------------------

RULES = {
    "molecular_weight": 500,
    "logp": 5,
    "hbd": 5,
    "hba": 10,
    "tpsa": 140,
    "rotatable_bonds": 10
}

VIOLATION_MEANINGS = {
    "molecular_weight": {"category": "size", "issue": "Molecule may be too large for oral absorption"},
    "logp": {"category": "lipophilicity", "issue": "Excessive hydrophobicity may reduce solubility"},
    "hbd": {"category": "polarity", "issue": "Too many hydrogen bond donors"},
    "hba": {"category": "polarity", "issue": "Too many hydrogen bond acceptors"},
    "tpsa": {"category": "polarity", "issue": "High polar surface area may hinder permeability"},
    "rotatable_bonds": {"category": "flexibility", "issue": "High molecular flexibility may reduce binding efficiency"}
}

SUGGESTION_MAP = {
    "size": ["Remove bulky substituents", "Simplify ring systems"],
    "lipophilicity": ["Add polar functional groups", "Reduce hydrophobic substituents"],
    "polarity": ["Reduce hydrogen bond donors/acceptors", "Use bioisosteres"],
    "flexibility": ["Reduce rotatable bonds", "Introduce ring systems"]
}

# ------------------ CORE ANALYSIS ------------------

def analyze_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "valid": False, "error": "Invalid SMILES"}

    properties = {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": Descriptors.TPSA(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol)
    }

    violations = []
    for prop, limit in RULES.items():
        if properties[prop] > limit:
            meaning = VIOLATION_MEANINGS[prop]
            violations.append({
                "property": prop,
                "value": properties[prop],
                "limit": limit,
                "category": meaning["category"],
                "issue": meaning["issue"],
                "suggestions": SUGGESTION_MAP.get(meaning["category"], [])
            })

    score = max(0, 100 - 15 * len(violations))

    severe_flags = 0
    if properties["logp"] > 6: severe_flags += 1
    if properties["rotatable_bonds"] > 12: severe_flags += 1
    if properties["tpsa"] < 10: severe_flags += 1

    reject_reason = None
    if properties["molecular_weight"] < 150:
        reject_reason = "Too small to be a viable drug scaffold"
    if properties["hbd"] == 0 and properties["hba"] == 0:
        reject_reason = "No functional groups for target binding"

    if len(violations) == 0 and score >= 85:
        physchem_status = "OPTIMAL"
    elif len(violations) <= 1 and score >= 75 and severe_flags == 0:
        physchem_status = "NEAR_OPTIMAL"
    elif score >= 50:
        physchem_status = "POOR"
    else:
        physchem_status = "WEAK"

    decision = "REJECT" if reject_reason else "ACCEPT"

    return {
        "smiles": smiles,
        "valid": True,
        "mol": mol,
        "properties": properties,
        "violations": violations,
        "score": score,
        "physchem_status": physchem_status,
        "decision": decision,
        "reject_reason": reject_reason
    }

def analyze_multiple_smiles(smiles_list):
    return [analyze_smiles(s) for s in smiles_list]

def rank_molecules(results):
    return sorted(
        [r for r in results if r["valid"]],
        key=lambda x: x["score"],
        reverse=True
    )

def dataset_decision(results):
    decisions = [r["decision"] for r in results]
    if all(d == "REJECT" for d in decisions):
        return "ALL_REJECTED"
    if all(d == "ACCEPT" for d in decisions):
        return "ALL_ACCEPTED"
    return "MIXED"

def generate_optimization_advice(molecule):
    if not molecule["violations"]:
        return None
    advice = []
    for v in molecule["violations"]:
        advice.extend(v["suggestions"])
    return advice
