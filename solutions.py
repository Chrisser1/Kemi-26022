# solutions.py

from typing import List, Tuple
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
    assuming ideal (complete) dissociation. (Corrected Version)
    """
    # Use the robust helper function to get the correct van't Hoff factor
    ideal_factor = _get_ideal_van_hoff_factor(formula)
    
    # Calculate and return the total concentration
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
    """
    Calculates the ideal van't Hoff factor (i) by correctly identifying
    polyatomic and monatomic ions in a formula. (Fully Corrected Version)
    """
    c = Compound(formula)
    composition = c.composition.copy()
    ion_count = 0
    
    # List of common polyatomic ions
    polyatomic_ions = ['NH4', 'OH', 'NO3', 'SO4', 'CO3', 'PO4']

    # 1. Handle polyatomic ions in parentheses first, e.g., Al₂(SO₄)₃
    poly_in_parens = re.findall(r'\((\w+)\)(\d*)', formula)
    for part, count_str in poly_in_parens:
        count = int(count_str) if count_str else 1
        ion_count += count
        
        poly_comp = Compound(part).composition
        for atom, num_in_poly in poly_comp.items():
            composition[atom] -= num_in_poly * count
            if composition[atom] == 0:
                del composition[atom]

    # 2. Handle common polyatomic ions not in parentheses, e.g., KNO₃
    for ion in polyatomic_ions:
        poly_comp = Compound(ion).composition
        # Check if the remaining atoms can form this polyatomic ion
        if all(composition.get(atom, 0) >= num for atom, num in poly_comp.items()):
            # Find how many times this ion can be formed
            num_ions = min(composition[atom] // num for atom, num in poly_comp.items())
            ion_count += num_ions
            # Subtract the atoms of the polyatomic ion from the total composition
            for atom, num in poly_comp.items():
                composition[atom] -= num * num_ions
                if composition[atom] == 0:
                    del composition[atom]

    # 3. Count the remaining atoms as monatomic ions
    for atom, count in composition.items():
        ion_count += count
        
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

def calculate_colligative_constant(delta_T: float, molality: float, i: int = 1) -> float:
    """
    Calculates a colligative constant (Kb or Kf) from experimental data.
    K = ΔT / (i * m)
    
    Args:
        delta_T: The measured change in temperature (°C).
        molality: The molality of the solution (m).
        i: The van't Hoff factor for the solute (defaults to 1).
        
    Returns:
        The colligative constant in °C/m.
    """
    if molality == 0 or i == 0:
        raise ValueError("Molality and van't Hoff factor cannot be zero.")
    return delta_T / (i * molality)

def raoults_law_moles_solute(
    solute_formula: str,
    n_solvent: float,
    p_solution: float,
    p_solvent_pure: float
) -> float:
    """
    Calculates the moles of an ionic solute in a solution using Raoult's Law.
    """
    if p_solvent_pure == 0:
        raise ValueError("Vapor pressure of pure solvent cannot be zero.")
    
    # Step 1: Calculate mole fraction of the solvent
    x_solvent = p_solution / p_solvent_pure
    
    # Step 2: Get the van't Hoff factor for the solute using the corrected helper
    i = _get_ideal_van_hoff_factor(solute_formula)
    
    # Step 3: Calculate total moles of all dissolved particles
    if x_solvent >= 1:
        return 0
    n_particles = (n_solvent / x_solvent) - n_solvent
    
    # Step 4: Calculate moles of the solute compound by dividing by i
    n_solute = n_particles / i
    return n_solute

def rank_by_colligative_effect(formulas: list[str]) -> str:
    """
    Ranks electrolytes by the magnitude of their colligative effect
    (e.g., lowest to highest boiling point) based on their van't Hoff factor.
    """
    ranked_molecules = []
    for f in formulas:
        i = _get_ideal_van_hoff_factor(f)
        ranked_molecules.append((f, i))
    
    ranked_molecules.sort(key=lambda x: x[1])
    
    ranked_list = [item[0] for item in ranked_molecules]
    return " < ".join(ranked_list)

def raoults_law_total_pressure(components: List[Tuple[float, float]]) -> float:
    """
    Calculates the total vapor pressure of an ideal solution of volatile liquids.
    
    Args:
        components: A list of tuples, where each tuple contains
                    (mole_fraction_X, pure_vapor_pressure_P_X).
                    
    Returns:
        The total vapor pressure of the solution.
    """
    total_pressure = 0
    for x, p_pure in components:
        partial_pressure = x * p_pure
        total_pressure += partial_pressure
    return total_pressure

def calculate_molar_mass_from_bp_elevation(
    mass_solute_g: float,
    mass_solvent_g: float,
    delta_Tb: float,
    kb_solvent: float,
    i: int = 1
) -> float:
    """
    Calculates the molar mass of a solute from boiling point elevation data.
    """
    if i == 0 or kb_solvent == 0:
        raise ValueError("van't Hoff factor and Kb cannot be zero.")
    
    # 1. Calculate molality from ΔTb
    molality = delta_Tb / (i * kb_solvent)
    
    # 2. Calculate moles of solute from molality and solvent mass
    solvent_kg = mass_solvent_g / 1000.0
    moles_solute = molality * solvent_kg
    
    # 3. Calculate molar mass
    if moles_solute == 0:
        raise ValueError("Calculated moles of solute is zero.")
    molar_mass = mass_solute_g / moles_solute
    return molar_mass

def raoults_law_solution_pressure(
    x_solvent: float,
    p_solvent_pure: float
) -> float:
    """
    Calculates the vapor pressure of an ideal solution using Raoult's Law.
    
    Args:
        x_solvent: The mole fraction of the solvent.
        p_solvent_pure: The vapor pressure of the pure solvent.
        
    Returns:
        The vapor pressure of the solution in the same units as p_solvent_pure.
    """
    return x_solvent * p_solvent_pure

def calculate_molality(
    moles_solute: float,
    mass_solvent: float,
    solvent_unit: str = 'g'
) -> float:
    """
    Calculates the molality of a solution.
    
    Args:
        moles_solute: The moles of the solute.
        mass_solvent: The mass of the solvent.
        solvent_unit: The unit of the solvent's mass (e.g., 'g', 'kg').
        
    Returns:
        The molality of the solution in mol/kg (m).
    """
    # Convert solvent mass to kg if it's in grams
    solvent_kg = mass_solvent / 1000.0 if solvent_unit.lower() == 'g' else mass_solvent
    
    if solvent_kg == 0:
        raise ValueError("Mass of solvent cannot be zero.")
        
    return moles_solute / solvent_kg

def calculate_dilution(M1=None, V1=None, M2=None, V2=None):
    """
    Solves for a single unknown variable in the dilution equation M1*V1 = M2*V2.
    Provide the three known values and 'None' for the one to be calculated.
    Units for concentration and volume must be consistent.
    """
    knowns = {'M1': M1, 'V1': V1, 'M2': M2, 'V2': V2}
    unknowns = [k for k, v in knowns.items() if v is None]

    if len(unknowns) != 1:
        raise ValueError("Exactly one variable must be None to solve for.")

    var_to_solve = unknowns[0]

    if var_to_solve == 'V1':
        result = (M2 * V2) / M1
        return {'V1': result}
    elif var_to_solve == 'M1':
        result = (M2 * V2) / V1
        return {'M1': result}
    elif var_to_solve == 'V2':
        result = (M1 * V1) / M2
        return {'V2': result}
    elif var_to_solve == 'M2':
        result = (M1 * V1) / V2
        return {'M2': result}