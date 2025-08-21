import math
import re

from compound import Compound

def calculate_pOH(oh_concentration: float) -> float:
    """Calculates pOH from the hydroxide ion concentration."""
    return -math.log10(oh_concentration)

def convert_pOH_to_pH(pOH: float) -> float:
    """Converts a pOH value to a pH value (at 25°C)."""
    return 14.0 - pOH

def describe_solution_from_oh(oh_concentration: float) -> str:
    """Calculates pH from [OH⁻] and describes the solution."""
    pOH = calculate_pOH(oh_concentration)
    pH = convert_pOH_to_pH(pOH)
    
    description = f"The pOH is {pOH:.1f} and the pH is {pH:.1f}. "
    
    if pH < 3:
        description += "The solution is highly acidic."
    elif pH < 6:
        description += "The solution is weakly acidic."
    elif pH <= 8:
        description += "The solution is neutral."
    elif pH <= 11:
        description += "The solution is weakly basic."
    else:
        description += "The solution is highly basic."
        
    return description

def calculate_pH(h_concentration: float) -> float:
    """Calculates pH from the hydrogen ion concentration."""
    if h_concentration <= 0:
        raise ValueError("Hydrogen ion concentration must be positive.")
    return -math.log10(h_concentration)

# --- Data for Salt Acidity ---
STRONG_ACIDS = ["HCl", "HBr", "HI", "HNO3", "H2SO4", "HClO4", "HClO3"]
STRONG_BASES = ["LiOH", "NaOH", "KOH", "RbOH", "CsOH", "Ca(OH)2", "Sr(OH)2", "Ba(OH)2"]
WEAK_BASE_CATIONS = ["NH4"] # Cations that are conjugate acids of weak bases

def predict_salt_acidity(formula: str) -> str:
    """
    Predicts if a salt solution will be acidic, basic, or neutral.
    """
    # --- 1. Split the salt into cation and anion parts ---
    # This is a simplified parser for common salt types
    cation_part = ""
    anion_part = ""
    
    # Check for ammonium cation first
    if formula.startswith("NH4"):
        cation_part = "NH4"
        anion_part = formula[3:]
    else:
        # Find the first metal cation
        match = re.match(r'([A-Z][a-z]*)(\d*)', formula)
        if match:
            cation_part = match.group(1)
            anion_part = formula[len(match.group(0)):]

    if not cation_part or not anion_part:
        return f"Could not parse the salt '{formula}'."

    # --- 2. Analyze the Cation ---
    is_cation_acidic = cation_part in WEAK_BASE_CATIONS
    
    # --- 3. Analyze the Anion ---
    # Form the parent acid by adding H+
    # This is a simplification; H₂SO₄ from SO₄, etc., would need more logic
    parent_acid = "H" + anion_part
    is_anion_basic = parent_acid not in STRONG_ACIDS

    # --- 4. Determine the overall acidity ---
    if is_cation_acidic and not is_anion_basic:
        return "Acidic"
    elif not is_cation_acidic and is_anion_basic:
        return "Basic"
    elif not is_cation_acidic and not is_anion_basic:
        return "Neutral"
    else: # Both are from weak parents
        return "Depends on Ka and Kb values"
    
def rank_oxyacids_by_strength(formulas: list[str]) -> str:
    """
    Ranks a list of simple oxyacids (HXO) from strongest to weakest
    based on the electronegativity of the central atom X.
    """
    acid_strengths = []
    
    for formula in formulas:
        try:
            # Assumes formula is in HXO format, e.g., "HClO"
            central_atom = formula.replace("H", "").replace("O", "")
            
            # Get electronegativity from the Compound class's internal dictionary
            en_value = Compound._EN.get(central_atom)
            
            if en_value is None:
                raise ValueError(f"Electronegativity for '{central_atom}' not found.")
                
            acid_strengths.append((formula, en_value))
        except Exception as e:
            print(f"Could not process {formula}: {e}")
            
    # Sort the list in descending order of electronegativity (strongest first)
    acid_strengths.sort(key=lambda x: x[1], reverse=True)
    
    # Format the ranked list into a string
    ranked_list = [acid[0] for acid in acid_strengths]
    return " > ".join(ranked_list)

# Ion-product constant for water at 25°C
KW_25C = 1.0e-14

def calculate_hydroxide_from_hydronium(h_concentration: float) -> float:
    """Calculates [OH⁻] from [H⁺] using Kw."""
    if h_concentration <= 0:
        raise ValueError("Hydrogen ion concentration must be positive.")
    return KW_25C / h_concentration

def calculate_hydronium_from_hydroxide(oh_concentration: float) -> float:
    """Calculates [H⁺] from [OH⁻] using Kw."""
    if oh_concentration <= 0:
        raise ValueError("Hydroxide ion concentration must be positive.")
    return KW_25C / oh_concentration

def calculate_weak_acid_ph(initial_conc: float, Ka: float) -> float:
    """
    Calculates the pH of a weak acid solution using the small x approximation.
    """
    if initial_conc <= 0 or Ka <= 0:
        raise ValueError("Initial concentration and Ka must be positive.")
        
    # Ka ≈ x² / [HA]₀  =>  x = sqrt(Ka * [HA]₀)
    # where x is the equilibrium [H₃O⁺] concentration
    h3o_conc = math.sqrt(Ka * initial_conc)
    
    # pH = -log[H₃O⁺]
    return -math.log10(h3o_conc)

# Ion-product constant for water at 25°C
KW_25C = 1.0e-14

def convert_ka_kb(k_value: float) -> float:
    """
    Converts between Ka and Kb for a conjugate acid-base pair at 25°C.
    """
    if k_value <= 0:
        raise ValueError("K value must be positive.")
    return KW_25C / k_value

def rank_oxyacids_by_strength(formulas: list[str]) -> str:
    """
    Ranks a list of simple oxyacids (HXO) from strongest to weakest
    based on the electronegativity of the central atom X.
    """
    acid_strengths = []
    
    for formula in formulas:
        try:
            # Assumes formula is in HXO format, e.g., "HClO"
            # This logic finds the atom that is not H or O
            central_atom = formula.replace("H", "").replace("O", "")
            
            # Get electronegativity from the Compound class's internal dictionary
            en_value = Compound._EN.get(central_atom)
            
            if en_value is None:
                raise ValueError(f"Electronegativity for '{central_atom}' not found.")
                
            acid_strengths.append((formula, en_value))
        except Exception as e:
            print(f"Could not process {formula}: {e}")
            
    # Sort the list in descending order of electronegativity (strongest first)
    acid_strengths.sort(key=lambda x: x[1], reverse=True)
    
    # Format the ranked list into a string
    ranked_list = [acid[0] for acid in acid_strengths]
    return " > ".join(ranked_list)

def _get_acid_protons(acid_formula: str) -> int:
    """Helper to find the number of acidic protons (na)."""
    if acid_formula in ["H2SO4", "H2CO3"]: return 2
    if acid_formula == "H3PO4": return 3
    return 1

def _get_base_hydroxides(base_formula: str) -> int:
    """Helper to find the number of hydroxide ions (nb)."""
    if "(OH)2" in base_formula: return 2
    if "(OH)3" in base_formula: return 3
    return 1

def solve_titration(
    Ma: float = None, Va_mL: float = None, na: int = 1,
    Mb: float = None, Vb_mL: float = None, nb: int = 1
) -> dict[str, float]:
    """
    Solves for an unknown variable in an acid-base titration using the
    formula na*Ma*Va = nb*Mb*Vb.
    Use 'None' for the single variable you want to solve for.

    Args:
        Ma: Molarity of the acid.
        Va_mL: Volume of the acid in mL.
        na: Stoichiometric coefficient of the acid (protons it donates). Defaults to 1.
        Mb: Molarity of the base.
        Vb_mL: Volume of the base in mL.
        nb: Stoichiometric coefficient of the base (protons it accepts). Defaults to 1.

    Returns:
        A dictionary containing the value of the unknown variable.
    """
    params = {'Ma': Ma, 'Va_mL': Va_mL, 'Mb': Mb, 'Vb_mL': Vb_mL}
    unknown = [k for k, v in params.items() if v is None]

    if len(unknown) != 1:
        raise ValueError("Exactly one variable (Ma, Va_mL, Mb, or Vb_mL) must be None.")

    unknown_var = unknown[0]

    if unknown_var == 'Va_mL':
        result = (nb * Mb * Vb_mL) / (na * Ma)
        return {'Va_mL': result}
    elif unknown_var == 'Ma':
        result = (nb * Mb * Vb_mL) / (na * Va_mL)
        return {'Ma': result}
    elif unknown_var == 'Vb_mL':
        result = (na * Ma * Va_mL) / (nb * Mb)
        return {'Vb_mL': result}
    elif unknown_var == 'Mb':
        result = (na * Ma * Va_mL) / (nb * Vb_mL)
        return {'Mb': result}
    
def calculate_buffer_h3o_conc(Ka: float, conc_acid: float, conc_base: float) -> float:
    """
    Calculates the hydronium concentration of a buffer solution.
    [H₃O⁺] = Ka * ([HA] / [A⁻])
    """
    if conc_base == 0:
        raise ValueError("Concentration of the conjugate base cannot be zero.")
    
    return Ka * (conc_acid / conc_base)

def calculate_buffer_pH(Ka: float, conc_acid: float, conc_base: float) -> float:
    """
    Calculates the pH of a buffer solution using the Henderson-Hasselbalch equation.
    pH = pKa + log([A⁻]/[HA])
    """
    if conc_acid == 0 or conc_base == 0:
        raise ValueError("Concentrations cannot be zero.")
        
    pKa = -math.log10(Ka)
    return pKa + math.log10(conc_base / conc_acid)

def calculate_buffer_ph_from_mass(
    Ka: float,
    mass_acid_g: float,
    fm_acid: float,
    mass_base_g: float,
    fm_base: float
) -> float:
    """
    Calculates the pH of a buffer solution directly from the mass and
    formula mass (FM) of its components.

    Args:
        Ka: The acid dissociation constant of the weak acid.
        mass_acid_g: Mass of the weak acid in grams.
        fm_acid: Formula mass of the weak acid (g/mol).
        mass_base_g: Mass of the conjugate base in grams.
        fm_base: Formula mass of the conjugate base (g/mol).

    Returns:
        The pH of the buffer solution.
    """
    # 1. Calculate moles of acid and base
    moles_acid = mass_acid_g / fm_acid
    moles_base = mass_base_g / fm_base

    if moles_acid <= 0 or moles_base <= 0:
        raise ValueError("Calculated moles must be positive.")

    # 2. Calculate pKa
    pKa = -math.log10(Ka)

    # 3. Use Henderson-Hasselbalch with the mole ratio
    pH = pKa + math.log10(moles_base / moles_acid)
    
    return pH

def calculate_weak_base_pOH(initial_conc: float, Kb: float) -> float:
    """
    Calculates the pOH of a weak base solution using the small x approximation.
    """
    if initial_conc <= 0 or Kb <= 0:
        raise ValueError("Initial concentration and Kb must be positive.")
        
    # Kb ≈ x² / [B]₀  =>  x = sqrt(Kb * [B]₀)
    # where x is the equilibrium [OH⁻] concentration
    oh_conc = math.sqrt(Kb * initial_conc)
    
    # pOH = -log[OH⁻]
    return -math.log10(oh_conc)

def calculate_weak_acid_dissociation(
    initial_conc: float,
    Ka: float
) -> dict[str, float]:
    """
    Calculates equilibrium concentrations for a weak acid solution
    using the small x approximation.
    """
    if initial_conc <= 0 or Ka <= 0:
        raise ValueError("Initial concentration and Ka must be positive.")
        
    # x = [H₃O⁺] = [A⁻]
    x = math.sqrt(Ka * initial_conc)
    
    concentrations = {
        "H3O+": x,
        "A-": x, # A- represents the conjugate base, e.g., HS⁻
        "HA": initial_conc - x # HA represents the weak acid, e.g., H₂S
    }
    return concentrations


def rank_acids_by_inductive_effect(formulas: list[str]) -> str:
    """
    Ranks a list of haloacetic acids based on the inductive effect of the halogens.
    """
    acid_strengths = []
    halogens = ["F", "Cl", "Br", "I"]

    for formula in formulas:
        try:
            # Create a score based on the number and electronegativity of halogens
            strength_score = 0
            c = Compound(formula)
            for atom, count in c.composition.items():
                if atom in halogens:
                    en_value = Compound._EN.get(atom, 0)
                    strength_score += count * en_value
            
            acid_strengths.append((formula, strength_score))
        except Exception as e:
            print(f"Could not process {formula}: {e}")
    
    # Sort by the calculated strength score in ascending order (weakest to strongest)
    acid_strengths.sort(key=lambda x: x[1])
    
    ranked_list = [acid[0] for acid in acid_strengths]
    return " < ".join(ranked_list)

def calculate_hydronium_from_pH(pH: float) -> float:
    """Calculates [H⁺] concentration from pH."""
    return 10**(-pH)

def calculate_hydroxide_from_pOH(pOH: float) -> float:
    """Calculates [OH⁻] concentration from pOH."""
    return 10**(-pOH)

def calculate_buffer_ph_from_pka(
    pKa: float,
    conc_acid: float,
    conc_base: float
) -> float:
    """
    Calculates the pH of a buffer solution from pKa using the
    Henderson-Hasselbalch equation.
    """
    if conc_acid <= 0 or conc_base <= 0:
        raise ValueError("Concentrations must be positive.")
    
    return pKa + math.log10(conc_base / conc_acid)

def calculate_pka_from_buffer_ph(
    pH: float,
    conc_acid: float,
    conc_base: float
) -> float:
    """
    Calculates the pKa of a weak acid from the pH of a buffer solution.
    pKa = pH - log([A⁻]/[HA])
    """
    if conc_acid <= 0 or conc_base <= 0:
        raise ValueError("Concentrations must be positive.")
    
    return pH - math.log10(conc_base / conc_acid)

def calculate_buffer_after_addition(
    initial_vol_L: float,
    initial_conc_acid: float,
    initial_conc_base: float,
    added_vol_L: float,
    added_conc: float,
    added_is_acid: bool
) -> dict[str, float]:
    """
    Calculates the final concentrations of a buffer's components after
    a strong acid or base is added.

    Returns a dictionary with the final concentrations of the acid and base.
    """
    # 1. Calculate initial moles
    moles_acid = initial_vol_L * initial_conc_acid
    moles_base = initial_vol_L * initial_conc_base
    moles_added = added_vol_L * added_conc

    # 2. Perform stoichiometric reaction
    if added_is_acid:
        moles_base -= moles_added
        moles_acid += moles_added
    else: # Added a base
        moles_acid -= moles_added
        moles_base += moles_added

    # Check for buffer capacity exceeded
    if moles_acid < 0 or moles_base < 0:
        raise ValueError("Buffer capacity has been exceeded.")

    # 3. Calculate new total volume and final concentrations
    total_volume = initial_vol_L + added_vol_L
    final_conc_acid = moles_acid / total_volume
    final_conc_base = moles_base / total_volume

    return {
        "final_conc_acid": final_conc_acid,
        "final_conc_base": final_conc_base
    }

def compare_buffer_capacities(buffers: list[dict]) -> str:
    """
    Analyzes a list of buffers and determines which has the greatest capacity.

    Args:
        buffers: A list of dictionaries, where each dict has 'conc_acid'
                 and 'conc_base' keys.

    Returns:
        A string report explaining which buffer is best and why.
    """
    if not buffers:
        return "No buffers provided for comparison."

    best_buffer = None
    max_capacity_score = -1

    report = "--- Buffer Capacity Analysis ---\n\n"
    report += "| Buffer | Total Conc. (M) | Ratio [Base]/[Acid] |\n"
    report += "|:-------|:---------------:|:--------------------:|\n"

    for i, b in enumerate(buffers):
        total_conc = b['conc_acid'] + b['conc_base']
        ratio = b['conc_base'] / b['conc_acid']
        
        # A simple scoring: prioritize total concentration, but penalize skewed ratios.
        # A ratio far from 1 makes a buffer less effective.
        ratio_penalty = abs(1 - ratio)
        capacity_score = total_conc / (1 + ratio_penalty) # Higher is better

        report += f"|   #{i+1}   |      {total_conc:.3f}      |        {ratio:.2f}         |\n"

        if capacity_score > max_capacity_score:
            max_capacity_score = capacity_score
            best_buffer_index = i

    best = buffers[best_buffer_index]
    report += "\n--- Conclusion ---\n"
    report += (
        f"Buffer #{best_buffer_index + 1} ({best['conc_acid']} M Acid / {best['conc_base']} M Base) "
        "has the greatest buffering capacity.\n"
        "This is because it has the highest concentration of buffer components "
        "while maintaining a ratio close to 1."
    )
    return report

def titration_ph_at_half_equivalence(Ka_values: list[float], point_number: int) -> float:
    """
    Calculates the pH at a half-equivalence point for a titration.

    Args:
        Ka_values: A list of the acid dissociation constants (e.g., [Ka1, Ka2, ...]).
        point_number: The half-equivalence point of interest (e.g., 1 for the first).

    Returns:
        The pH value at that specific point.
    """
    if not 1 <= point_number <= len(Ka_values):
        raise IndexError("The requested half-equivalence point number is out of range for the given Ka values.")

    # At the nth half-equivalence point, pH = pKa_n
    Ka = Ka_values[point_number - 1]
    pKa = -math.log10(Ka)
    
    return pKa

def classify_lewis_acid_base(compound: Compound) -> str:
    """
    Classifies a chemical species as a Lewis acid, Lewis base, or neither,
    based on a provided Compound object.
    """
    try:
        formula = compound.input_formula or compound.formula()

        # Rule 1: Cations are Lewis acids
        if compound.charge is not None and compound.charge > 0:
            return "Lewis Acid (cation)"

        # Rule 2: Common molecules with incomplete octets are Lewis acids
        central_atom = min(compound.composition, key=compound.composition.get)
        if central_atom in ['B', 'Al']: # Group 13 elements
            return f"Lewis Acid (incomplete octet on {central_atom})"

        # Rule 3: Molecules with a central atom that can be an electron acceptor
        if formula in ["CO2", "SO2", "SO3"]:
            return "Lewis Acid (central atom can accept e- pair)"

        # Rule 4: Common molecules with lone pairs are Lewis bases
        if formula in ["H2O", "NH3", "OH"]:
            return "Lewis Base (has lone pair(s) to donate)"
        if compound.charge is not None and compound.charge < 0 and len(compound.composition) == 1:
            return "Lewis Base (anion with lone pair(s))" # e.g., F-, Cl-

        # Rule 5: Neutral, elemental atoms are typically neither
        if len(compound.composition) == 1 and compound.charge is None:
            return "Neither (neutral atom)"

        return "Classification uncertain based on simple rules."

    except Exception as e:
        return f"Could not analyze compound '{compound}': {e}"
    
def calculate_buffer_ph_after_addition_from_moles(
    initial_moles_acid: float,
    initial_moles_base: float,
    added_moles: float,
    added_is_acid: bool,
    Ka: float
) -> float:
    """
    Calculates the final pH of a buffer after a strong acid or base is added,
    starting from the initial moles of the buffer components.

    Args:
        initial_moles_acid: Initial moles of the weak acid.
        initial_moles_base: Initial moles of the conjugate base.
        added_moles: Moles of the strong acid or base being added.
        added_is_acid: True if adding a strong acid, False for a strong base.
        Ka: The acid dissociation constant of the weak acid.

    Returns:
        The final pH of the buffer solution.
    """
    # 1. Perform stoichiometric reaction
    if added_is_acid:
        final_moles_base = initial_moles_base - added_moles
        final_moles_acid = initial_moles_acid + added_moles
    else:  # Added a base
        final_moles_acid = initial_moles_acid - added_moles
        final_moles_base = initial_moles_base + added_moles

    # Check for buffer capacity exceeded
    if final_moles_acid <= 0 or final_moles_base <= 0:
        raise ValueError("Buffer capacity has been exceeded.")

    # 2. Calculate pKa
    pKa = -math.log10(Ka)

    # 3. Use Henderson-Hasselbalch with the final mole ratio
    pH = pKa + math.log10(final_moles_base / final_moles_acid)
    
    return pH

def calculate_titration_ph_at_equivalence(
    initial_conc_acid: float,
    initial_vol_mL_acid: float,
    Ka_acid: float,
    conc_base: float
) -> float:
    """
    Calculates the pH at the equivalence point for a weak acid-strong base titration.
    """
    # 1. Calculate initial moles of acid and volume of base needed
    initial_moles_acid = initial_conc_acid * (initial_vol_mL_acid / 1000)
    vol_L_base = initial_moles_acid / conc_base
    
    # 2. Calculate total volume and concentration of conjugate base
    total_volume_L = (initial_vol_mL_acid / 1000) + vol_L_base
    conc_conjugate_base = initial_moles_acid / total_volume_L
    
    # 3. Calculate Kb for the conjugate base
    Kb_base = convert_ka_kb(Ka_acid)
    
    # 4. Calculate pOH from the weak base equilibrium
    pOH = calculate_weak_base_pOH(conc_conjugate_base, Kb_base)
    
    # 5. Convert pOH to pH and return
    return convert_pOH_to_pH(pOH)