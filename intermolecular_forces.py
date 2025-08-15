import re
from compound import Compound
from molecular_structure import analyze_molecular_properties # CORRECTED IMPORT

# A dictionary to rank IMF strength
IMF_STRENGTH = {
    "London Dispersion Forces": 1,
    "Dipole-Dipole": 2,
    "Hydrogen Bonding": 3
}

def get_predominant_imf(formula: str) -> tuple[str, float]:
    """
    Determines the predominant intermolecular force for a given molecule by parsing
    the analysis string.
    Returns the IMF type and a secondary ranking metric (molar mass or ΔEN).
    """
    analysis = analyze_molecular_properties(formula)
    c = Compound(formula)

    # --- Determine Polarity from the analysis string ---
    is_polar = "Polarity: Polar" in analysis
    
    # --- Check for Hydrogen Bonding ---
    # Condition: H is present and bonded to N, O, or F.
    has_h = "H" in c.composition
    has_nof = any(x in ["N", "O", "F"] for x in c.composition)
    # A simple check for molecules like NH3, H2O, HF
    can_h_bond = has_h and has_nof and is_polar

    # --- Extract ΔEN for ranking polar molecules ---
    max_en_diff = 0.0
    match = re.search(r"Max Bond ΔEN: (\d+\.\d+)", analysis)
    if match:
        max_en_diff = float(match.group(1))

    # --- Assign the predominant IMF ---
    if can_h_bond:
        return "Hydrogen Bonding", max_en_diff
    elif is_polar:
        return "Dipole-Dipole", max_en_diff
    else:
        # For nonpolar molecules, the ranking metric is molar mass
        return "London Dispersion Forces", c.molar_mass

def rank_by_boiling_point(formulas: list[str]) -> str:
    """
    Ranks a list of molecules by increasing boiling point based on their IMFs.
    """
    ranked_molecules = []
    for f in formulas:
        imf_type, metric = get_predominant_imf(f)
        ranked_molecules.append({
            "formula": f,
            "imf": imf_type,
            "strength": IMF_STRENGTH[imf_type],
            "metric": metric
        })
    
    # Sort first by IMF strength, then by the secondary metric
    ranked_molecules.sort(key=lambda x: (x["strength"], x["metric"]))
    
    ranked_list = [m["formula"] for m in ranked_molecules]
    return " < ".join(ranked_list)