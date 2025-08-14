# intermolecular_forces.py

from compound import Compound
from molecular_structure import predict_molecular_properties # Assuming this is in the same directory

# A dictionary to rank IMF strength
IMF_STRENGTH = {
    "London Dispersion Forces": 1,
    "Dipole-Dipole": 2,
    "Hydrogen Bonding": 3
}

def get_predominant_imf(formula: str) -> tuple[str, float]:
    """
    Determines the predominant intermolecular force for a given molecule.
    Returns the IMF type and a secondary ranking metric (molar mass or ΔEN).
    """
    # Use the detailed analysis function we built
    analysis = predict_molecular_properties(formula)
    
    # Check for Hydrogen Bonding condition
    c = Compound(formula)
    has_h = "H" in c.composition
    has_nof = any(x in c.composition for x in ["N", "O", "F"])
    # This is a simplification; we'd need to know the actual bonds for a perfect check
    # But for simple molecules like NH3, H2O, HF it works.
    
    is_polar = "Polarity: Polar" in analysis
    
    if has_h and has_nof and is_polar:
        # For H-bonders, the secondary metric is electronegativity difference
        return "Hydrogen Bonding", analysis["max_en_diff"]
    elif is_polar:
        # For other polar molecules, the secondary metric is also ΔEN
        return "Dipole-Dipole", analysis["max_en_diff"]
    else:
        # For nonpolar molecules, the secondary metric is molar mass
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
    
    # Sort first by IMF strength, then by the secondary metric (molar mass or ΔEN)
    ranked_molecules.sort(key=lambda x: (x["strength"], x["metric"]))
    
    ranked_list = [m["formula"] for m in ranked_molecules]
    return " < ".join(ranked_list)