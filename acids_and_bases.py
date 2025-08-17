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
    acid_formula: str, Va_mL: float,
    base_formula: str, Vb_mL: float,
    known_conc: float, known_is_base: bool
) -> float:
    """
    Solves for the unknown molarity in an acid-base titration.
    na*Ma*Va = nb*Mb*Vb
    """
    na = _get_acid_protons(acid_formula)
    nb = _get_base_hydroxides(base_formula)
    
    if known_is_base:
        Mb = known_conc
        # Solve for Ma
        Ma = (nb * Mb * Vb_mL) / (na * Va_mL)
        return Ma
    else: # Known is acid
        Ma = known_conc
        # Solve for Mb
        Mb = (na * Ma * Va_mL) / (nb * Vb_mL)
        return Mb
    
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