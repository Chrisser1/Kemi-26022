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
