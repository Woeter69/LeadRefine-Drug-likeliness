from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Draw

# -----------------------------
# Rules
# -----------------------------
RULES = {
    "molecular_weight": 500,
    "logp": 5,
    "hbd": 5,
    "hba": 10,
    "tpsa": 140,
    "rotatable_bonds": 10
}

VIOLATION_MEANINGS = {
    "molecular_weight": {"category": "size", "issue": "Too large for oral absorption"},
    "logp": {"category": "lipophilicity", "issue": "Poor solubility risk"},
    "hbd": {"category": "polarity", "issue": "Too many H-bond donors"},
    "hba": {"category": "polarity", "issue": "Too many H-bond acceptors"},
    "tpsa": {"category": "polarity", "issue": "Poor membrane permeability"},
    "rotatable_bonds": {"category": "flexibility", "issue": "Excessive flexibility"}
}

SUGGESTION_MAP = {
    "size": ["Remove bulky groups", "Simplify scaffold"],
    "lipophilicity": ["Add polar groups", "Reduce alkyl chains"],
    "polarity": ["Reduce donors/acceptors", "Use bioisosteres"],
    "flexibility": ["Reduce rotatable bonds", "Add rings"]
}

# -----------------------------
# Core analysis
# -----------------------------
def analyze_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "valid": False}

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
            v = VIOLATION_MEANINGS[prop]
            violations.append({
                "property": prop,
                "value": properties[prop],
                "limit": limit,
                "issue": v["issue"],
                "suggestions": SUGGESTION_MAP[v["category"]]
            })

    # -----------------------------
    # Scoring
    # -----------------------------
    score = 100
    score -= 15 * len(violations)

    for prop, limit in RULES.items():
        score -= max(0, properties[prop] - limit) * 0.1

    # Hard rejection rules
    reject_reason = None
    if properties["molecular_weight"] < 150:
        reject_reason = "Too small to be a viable drug scaffold"
    if properties["hbd"] == 0 and properties["hba"] == 0:
        reject_reason = "No functional groups for target binding"

    # 🔥 CRITICAL FIX
    if reject_reason:
        score -= 40

    score = max(score, 0)

    # Status
    if score >= 85 and not violations:
        status = "OPTIMAL"
    elif score >= 70:
        status = "NEAR_OPTIMAL"
    elif score >= 50:
        status = "WEAK"
    else:
        status = "POOR"

    decision = "REJECT" if reject_reason else "ACCEPT"

    return {
        "smiles": smiles,
        "valid": True,
        "mol": mol,
        "properties": properties,
        "violations": violations,
        "score": round(score, 2),
        "physchem_status": status,
        "decision": decision,
        "reject_reason": reject_reason
    }


def analyze_multiple_smiles(smiles_list):
    return [analyze_smiles(s) for s in smiles_list]


def rank_molecules(results):
    valid = [r for r in results if r["valid"]]
    return sorted(valid, key=lambda x: x["score"], reverse=True)


def dataset_decision(results):
    decisions = [r["decision"] for r in results]
    if all(d == "REJECT" for d in decisions):
        return "ALL_REJECTED"
    if all(d == "ACCEPT" for d in decisions):
        return "ALL_ACCEPTED"
    return "MIXED"


# -----------------------------
# Explanation helpers
# -----------------------------
def explain_imperfection(mol, is_top=False, tied=False):
    reasons = []

    if mol["violations"]:
        reasons.append("Violates one or more physicochemical drug-likeness rules")

    if mol["properties"]["logp"] > 4.5:
        reasons.append("Lipophilicity is close to the upper acceptable limit")

    if mol["properties"]["rotatable_bonds"] > 8:
        reasons.append("Molecular flexibility may reduce binding specificity")

    if mol["properties"]["tpsa"] > 120:
        reasons.append("High polarity may reduce membrane permeability")

    # Clean, honest messaging
    if not reasons:
        if tied:
            reasons.append(
                "This molecule is equally optimal by all evaluated drug-likeness criteria"
            )
        elif is_top:
            reasons.append(
                "This molecule satisfies all evaluated physicochemical drug-likeness rules "
                "and is selected as the top representative candidate in this dataset"
            )
        else:
            reasons.append(
                "This molecule is optimal by current rules but ranks below the top due to "
                "relative differences within this dataset"
            )

    return reasons



def generate_optimization_advice(mol):
    advice = []
    for v in mol["violations"]:
        advice.extend(v["suggestions"])
    return list(set(advice))
