# solutions.py

from compound import Compound
import re

def analyze_van_hoff_factor(formula: str) -> str:
    """
    Calculates the ideal van't Hoff factor for an ionic compound and
    explains the concept of the real van't Hoff factor.
    """
    try:
        # --- Calculate Ideal Factor ---
        # A simple way to estimate the number of ions for simple salts
        # by counting uppercase letters and polyatomic groups.
        
        # Count polyatomic ions first
        polyatomic_ions = ['NH4', 'OH', 'NO3', 'SO4', 'CO3', 'PO4']
        ion_count = 0
        temp_formula = formula
        for ion in polyatomic_ions:
            if ion in temp_formula:
                ion_count += 1
                temp_formula = temp_formula.replace(ion, "")
        
        # Count remaining elements (assumed to be monatomic ions)
        ion_count += len(re.findall(r'[A-Z]', temp_formula))

        ideal_factor = ion_count
        
        # --- Generate Explanation ---
        summary = (
            f"--- Van't Hoff Factor Analysis for {formula} ---\n"
            f"  - Ideal Factor (i_ideal): {ideal_factor}\n"
            f"    (Calculated based on the dissociation into {ideal_factor} ions.)\n\n"
            f"  - Real Factor (i_real): The real van't Hoff factor is typically slightly less than the ideal value.\n"
            f"    This is due to 'ion pairing,' where some dissolved ions are attracted to each other\n"
            f"    and briefly act as a single particle, reducing the total effective number of particles.\n"
            f"    For {formula}, the real factor would be expected to be just under {ideal_factor}."
        )
        return summary

    except Exception as e:
        return f"Could not analyze {formula}: {e}"
    
def calculate_total_ion_concentration(formula: str, molarity: float) -> float:
    """
    Calculates the total ion concentration for a given ionic solution,
    assuming ideal (complete) dissociation.
    """
    # --- Logic to find the ideal van't Hoff factor (i) ---
    # Count polyatomic ions first
    polyatomic_ions = ['NH4', 'OH', 'NO3', 'SO4', 'CO3', 'PO4']
    ion_count = 0
    temp_formula = formula
    for ion in polyatomic_ions:
        if ion in temp_formula:
            # This simple logic assumes the polyatomic ion appears once.
            # A more advanced version would parse coefficients.
            ion_count += 1
            temp_formula = temp_formula.replace(ion, "")
    
    # Count the number of uppercase letters to find monatomic ions
    # and then parse any numbers that follow them.
    for match in re.finditer(r'([A-Z][a-z]*)(\d*)', temp_formula):
        element, count = match.groups()
        ion_count += int(count) if count else 1
            
    ideal_factor = ion_count
    
    # --- Calculate and return the total concentration ---
    return molarity * ideal_factor

def calculate_specific_ion_concentration(
    compound_formula: str,
    ion_formula: str,
    molarity: float
) -> float:
    """
    Calculates the concentration of a specific ion in a solution,
    assuming ideal (complete) dissociation.
    """
    # Parse the main compound to find the count of the target ion
    compound = Compound(compound_formula)
    ion = Compound(ion_formula)
    
    # Check how many times the ion's atoms fit into the compound's atoms
    # This handles simple cases like Cl in MgCl₂ (returns 2) or K in KCl (returns 1)
    
    # Get atom counts for the ion
    ion_atoms = ion.composition
    
    # Find how many of that ion are in the compound
    # For a simple ion like "Cl", this is straightforward
    if len(ion_atoms) == 1:
        ion_symbol = next(iter(ion_atoms))
        count = compound.composition.get(ion_symbol, 0)
    else: # A more complex polyatomic ion check would be needed for full generality
        count = 1 # Simplified assumption for now

    return molarity * count

def raoults_law_mole_fraction_solvent(p_solution: float, p_solvent_pure: float) -> float:
    """
    Calculates the mole fraction of the solvent in an ideal solution
    using Raoult's Law.
    
    Args:
        p_solution: The vapor pressure of the solution.
        p_solvent_pure: The vapor pressure of the pure solvent.
        (Pressures must be in the same units)
        
    Returns:
        The mole fraction of the solvent (X_solvent).
    """
    if p_solvent_pure == 0:
        raise ValueError("Vapor pressure of the pure solvent cannot be zero.")
    return p_solution / p_solvent_pure

def _get_ideal_van_hoff_factor(formula: str) -> int:
    """Helper function to calculate the ideal van't Hoff factor."""
    polyatomic_ions = ['NH4', 'OH', 'NO3', 'SO4', 'CO3', 'PO4']
    ion_count = 0
    temp_formula = formula
    for ion in polyatomic_ions:
        if ion in temp_formula:
            ion_count += 1
            temp_formula = temp_formula.replace(ion, "", 1)
    
    for match in re.finditer(r'([A-Z][a-z]*)(\d*)', temp_formula):
        _, count = match.groups()
        ion_count += int(count) if count else 1
            
    return ion_count

def calculate_freezing_point_depression(
    solute_formula: str,
    molality: float,
    kf_solvent: float
) -> float:
    """
    Calculates the freezing point depression of a solution.
    ΔTf = i * Kf * m
    """
    # 1. Determine the ideal van't Hoff factor
    i = _get_ideal_van_hoff_factor(solute_formula)
    
    # 2. Calculate the depression
    delta_tf = i * kf_solvent * molality
    return delta_tf

# --- Mass Conversion Helper ---
_MASS_TO_GRAM = {
    "g": 1.0,
    "kg": 1000.0,
    "mg": 1e-3,
    "ug": 1e-6, # microgram
}

def _to_gram(value: float, unit: str) -> float:
    """Converts a mass value to grams."""
    try:
        return value * _MASS_TO_GRAM[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown mass unit '{unit}'")

def calculate_ppb_from_mass(
    mass_solute: float,
    solute_unit: str,
    mass_solution: float,
    solution_unit: str
) -> float:
    """
    Calculates concentration in parts per billion (ppb) from mass with unit handling.
    """
    # Convert both masses to a consistent base unit (grams)
    solute_g = _to_gram(mass_solute, solute_unit)
    solution_g = _to_gram(mass_solution, solution_unit)
    
    if solution_g == 0:
        raise ValueError("Mass of solution cannot be zero.")
    
    return (solute_g / solution_g) * 1e9

def calculate_ppm_from_mass(
    mass_solute: float,
    solute_unit: str,
    mass_solution: float,
    solution_unit: str
) -> float:
    """
    Calculates concentration in parts per billion (ppb) from mass with unit handling.
    """
    # Convert both masses to a consistent base unit (grams)
    solute_g = _to_gram(mass_solute, solute_unit)
    solution_g = _to_gram(mass_solution, solution_unit)
    
    if solution_g == 0:
        raise ValueError("Mass of solution cannot be zero.")
    
    return (solute_g / solution_g) * 1e6
