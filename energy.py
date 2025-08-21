# energy.py

from collections import defaultdict
import math
import pandas as pd
import re
from typing import Dict, Sequence, Tuple
from compound import Compound
from equations import ChemicalEquation
from sympy import solve, symbols, Eq

# ─────────────────────────────────────────────────────────────────────────────
# 1) HEAT CAPACITY CALCULATIONS (q = m·c·ΔT or q = n·Cp·ΔT)
# ─────────────────────────────────────────────────────────────────────────────

def heat_mass(mass_g: float, specific_heat: float, delta_T: float) -> float:
    return mass_g * specific_heat * delta_T

def heat_moles(moles: float, molar_heat_cap: float, delta_T: float) -> float:
    return moles * molar_heat_cap * delta_T

def delta_T_from_heat_mass(q: float, mass_g: float, specific_heat: float) -> float:
    return q / (mass_g * specific_heat)

def delta_T_from_heat_moles(q: float, moles: float, molar_heat_cap: float) -> float:
    return q / (moles * molar_heat_cap)


# ─────────────────────────────────────────────────────────────────────────────
# 2) REACTION ENTHALPY (ΔH°rxn via standard formation enthalpies)
# ─────────────────────────────────────────────────────────────────────────────

def reaction_enthalpy(
    reaction: ChemicalEquation,
    formation_energies: Dict[Compound, float]
) -> float:
    ν_react, ν_prod = reaction.balance()
    ΔH_p = sum(ν * formation_energies[c]
               for c, ν in zip(reaction.products,  ν_prod))
    ΔH_r = sum(ν * formation_energies[c]
               for c, ν in zip(reaction.reactants, ν_react))
    return ΔH_p - ΔH_r


# ─────────────────────────────────────────────────────────────────────────────
# 3) REACTION ENTROPY & FREE ENERGY
# ─────────────────────────────────────────────────────────────────────────────

def reaction_entropy(
    reaction: ChemicalEquation,
    entropy_data: Dict[Compound, float]
) -> float:
    """
    ΔS°rxn = Σ ν·S°(products) − Σ ν·S°(reactants)
    Returns J/(mol·K).
    """
    ν_react, ν_prod = reaction.balance()
    ΔS_p = sum(ν * entropy_data.get(c, 0.0)
               for c, ν in zip(reaction.products,  ν_prod))
    ΔS_r = sum(ν * entropy_data.get(c, 0.0)
               for c, ν in zip(reaction.reactants, ν_react))
    return ΔS_p - ΔS_r

def heat_of_combustion(
    compound: Compound,
    mass_g: float,
    delta_h_rxn_per_mole: float
) -> float:
    """
    Calculates the total heat released for the combustion of a given mass of a substance.

    - compound: A Compound object representing the substance being burned.
    - mass_g: The mass of the substance in grams.
    - delta_h_rxn_per_mole: The standard enthalpy of reaction in kJ/mol for the substance.
    Returns the total heat released in kJ.
    """
    # 1. Set the mass of the compound to calculate the moles
    compound.set_mass(mass_g)
    
    # 2. Calculate the total heat released
    total_heat = compound.amount_mol * delta_h_rxn_per_mole
    
    return total_heat

def reaction_free_energy(
    reaction: ChemicalEquation,
    gibbs_data: Dict[Compound, float],
    temperature_K: float = 298.15
) -> float:
    """
    ΔG°rxn: try direct Σ ΔGf°; if missing, fall back to ΔH°–T·ΔS°.
    Returns kJ/mol.
    """
    ν_react, ν_prod = reaction.balance()
    # first attempt: direct ΔGf°
    try:
        ΔG_p = sum(ν * gibbs_data[c]
                   for c, ν in zip(reaction.products,  ν_prod))
        ΔG_r = sum(ν * gibbs_data[c]
                   for c, ν in zip(reaction.reactants, ν_react))
        return ΔG_p - ΔG_r
    except KeyError:
        # fallback: ΔH° – T·ΔS°  (convert ΔS J→kJ)
        ΔH = reaction_enthalpy(reaction, _Hf_map)
        ΔS = reaction_entropy(reaction, _S_map) / 1000.0
        return ΔH - temperature_K * ΔS


# ─────────────────────────────────────────────────────────────────────────────
# 4) CSV LOADER FOR ΔHf°, ΔGf°, S°
# ─────────────────────────────────────────────────────────────────────────────

# match “Substance” like "Na+(aq)" or "SO4^2-(aq)”
# capture both the core (formula+charge) and the phase
_SUBSTANCE_RE = re.compile(
    r'^\s*'
      r'(?P<core>.+?)'          # everything up to the "("
      r'\s*\(\s*'
      r'(?P<phase>aq|s|l|g)'    # Named phase group!
      r'\s*\)\s*$'
)
_ION_RE = re.compile(r'''
    ^\s*
    (?P<formula>[A-Za-z0-9\(\)]+?)   # e.g. Ag or SO4
    \^?                              # optional caret
    (?P<charge>(?:\d+[+\-]|[+\-]))   # "2+" or "+" or "-"
    \s*$
''', re.X)

def load_standard_thermo(csv_path: str
        ) -> Tuple[
            Dict[Compound,float],  # ΔHf° (kJ/mol)
            Dict[Compound,float],  # ΔGf° (kJ/mol)
            Dict[Compound,float]   # S° (J/(mol·K))
        ]:
    df = pd.read_csv(csv_path, skipinitialspace=True, quotechar='"')
    Hf_map: Dict[Compound,float] = {}
    Gf_map: Dict[Compound,float] = {}
    S_map:  Dict[Compound,float] = {}

    for _, row in df.iterrows():
        raw    = str(row['Substance']).strip()
        h_str  = str(row['ΔHf∘ (kJ/mol)']).strip()
        g_str  = str(row['ΔGf∘ (kJ/mol)']).strip()
        s_str  = str(row['S∘ (J/(mol·K))']).strip()
        # skip if all blank
        if pd.isna(h_str) and pd.isna(g_str) and pd.isna(s_str):
            continue

        # parse values (Unicode minus → ASCII minus)
        h_val = float(h_str.replace('−','-')) if h_str and not pd.isna(h_str) else None
        g_val = float(g_str.replace('−','-')) if g_str and not pd.isna(g_str) else None
        s_val = float(s_str.replace('−','-')) if s_str and not pd.isna(s_str) else None

        # split off phase
        m1 = _SUBSTANCE_RE.match(raw)
        if not m1:
            continue
        core, phase = m1.group('core'), m1.group('phase')

        # split off charge
        cm = _ION_RE.match(core)
        if cm:
            formula = cm.group('formula')
            charg_s = cm.group('charge')
            sign     = +1 if charg_s.endswith('+') else -1
            mag      = int(charg_s[:-1]) if len(charg_s)>1 else 1
            charge   = mag * sign
        else:
            formula, charge = core, None

        # build compound
        cmpd = Compound(formula)
        cmpd.phase  = phase
        cmpd.charge = charge

        if h_val is not None: Hf_map[cmpd] = h_val
        if g_val is not None: Gf_map[cmpd] = g_val
        if s_val is not None: S_map [cmpd] = s_val

    return Hf_map, Gf_map, S_map


# ─────────────────────────────────────────────────────────────────────────────
# 5) AUTO‐LOAD YOUR CSV ONCE
# ─────────────────────────────────────────────────────────────────────────────

_Hf_map, _Gf_map, _S_map = load_standard_thermo("enthalpy_data.csv")

# expose getters in case you want them:
def get_Hf_map() -> Dict[Compound,float]: return _Hf_map
def get_Gf_map() -> Dict[Compound,float]: return _Gf_map
def get_S_map()  -> Dict[Compound,float]: return _S_map

# Add near the top of energy.py

# Conversion factors to the base unit (Joule)
_ENERGY_TO_JOULE = {
    "j": 1.0,
    "kj": 1000.0,
    "cal": 4.184,
    "kcal": 4184.0,
    "btu": 1055.06,
}

def _to_joule(E: float, unit: str) -> float:
    try:
        return E * _ENERGY_TO_JOULE[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown energy unit '{unit}'")

def _from_joule(E_j: float, unit: str) -> float:
    try:
        return E_j / _ENERGY_TO_JOULE[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown energy unit '{unit}'")

def convert_energy(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert an energy value from one unit to another.
    """
    e_joule = _to_joule(value, from_unit)
    return _from_joule(e_joule, to_unit)


L_ATM_TO_JOULE = 101.325  # The conversion factor from the problem

def pressure_volume_work(delta_V: float, P: float) -> float:
    """
    Calculates pressure-volume work in Joules.
    Assumes P is in atm and ΔV is in Liters.
    w = -PΔV
    """
    work_L_atm = -P * delta_V
    return work_L_atm * L_ATM_TO_JOULE

def specific_heat_from_mass(q: float, mass_g: float, delta_T: float) -> float:
    """
    Calculates the specific heat capacity of a substance.
    c = q / (m * ΔT)
    """
    if mass_g == 0 or delta_T == 0:
        raise ValueError("Mass and delta_T cannot be zero.")
    return q / (mass_g * delta_T)

def parse_reaction_string(reaction_str: str) -> Dict[str, int]:
    """
    Parses a reaction string into a dictionary of species and their net coefficients.
    Handles integer and floating-point coefficients.
    """
    net_coeffs = defaultdict(int)
    
    # Split reaction into reactants and products
    try:
        reactants_str, products_str = reaction_str.split('->')
    except ValueError:
        raise ValueError(f"Invalid reaction format: '{reaction_str}'. Must contain '->'.")

    # --- Process reactants (negative coefficients) ---
    for term in reactants_str.split('+'):
        term = term.strip()
        parts = term.split()
        try:
            if len(parts) == 1:
                coeff, species = -1.0, parts[0]
            else:
                # Use float() to handle decimals like "0.5"
                coeff, species = -float(parts[0]), " ".join(parts[1:])
            net_coeffs[species] += coeff
        except (ValueError, IndexError):
            raise ValueError(f"Could not parse reactant term: '{term}'")

    # --- Process products (positive coefficients) ---
    for term in products_str.split('+'):
        term = term.strip()
        parts = term.split()
        try:
            if len(parts) == 1:
                coeff, species = 1.0, parts[0]
            else:
                # Use float() here as well
                coeff, species = float(parts[0]), " ".join(parts[1:])
            net_coeffs[species] += coeff
        except (ValueError, IndexError):
            raise ValueError(f"Could not parse product term: '{term}'")
            
    return net_coeffs

def parse_reaction_string_special(reaction_str: str) -> Dict[str, int]:
    """
    Parses a reaction string into a dictionary of species and their net coefficients.
    Handles integer and floating-point coefficients. (Corrected Version)
    """
    net_coeffs = defaultdict(float)
    
    # Regex to capture an optional coefficient and the species name that follows.
    # The species can include letters, numbers, and a state like (g) or (graphite).
    term_re = re.compile(r"(\d*\.?\d*)\s*([A-Za-z0-9()]+(?:\([a-z]+\))?)")

    try:
        reactants_str, products_str = reaction_str.split('->')
    except ValueError:
        raise ValueError(f"Invalid reaction format: '{reaction_str}'. Must contain '->'.")

    # --- Process reactants (negative coefficients) ---
    for term in reactants_str.split('+'):
        term = term.strip()
        match = term_re.match(term)
        if not match:
            raise ValueError(f"Could not parse reactant term: '{term}'")
        
        coeff_str, species = match.groups()
        coeff = float(coeff_str) if coeff_str else 1.0
        net_coeffs[species] -= coeff

    # --- Process products (positive coefficients) ---
    for term in products_str.split('+'):
        term = term.strip()
        match = term_re.match(term)
        if not match:
            raise ValueError(f"Could not parse product term: '{term}'")
            
        coeff_str, species = match.groups()
        coeff = float(coeff_str) if coeff_str else 1.0
        net_coeffs[species] += coeff
            
    return net_coeffs

def solve_hess_law(
    target_reaction_str: str,
    known_reaction_strs: list[str],
    known_enthalpies: list[float]
) -> tuple[float, dict[str, float]]:
    """
    Solves for the enthalpy of a target reaction using Hess's Law from reaction strings.
    """
    target_coeffs = parse_reaction_string_special(target_reaction_str)
    
    all_species = set(target_coeffs.keys())
    known_coeffs_list = []
    for r_str in known_reaction_strs:
        coeffs = parse_reaction_string_special(r_str)
        known_coeffs_list.append(coeffs)
        all_species.update(coeffs.keys())
        
    species_list = sorted(list(all_species))
    
    multipliers = symbols(f'x_:{len(known_reaction_strs)}')
    
    equations = []
    for species in species_list:
        lhs = sum(
            multipliers[i] * known_coeffs_list[i].get(species, 0)
            for i in range(len(known_reaction_strs))
        )
        rhs = target_coeffs.get(species, 0)
        equations.append(Eq(lhs, rhs))
        
    solution = solve(equations, multipliers)
    
    if not solution or not isinstance(solution, dict):
        raise ValueError("Could not solve for the reaction multipliers.")

    total_enthalpy = sum(solution[m] * h for m, h in zip(multipliers, known_enthalpies))
    
    multiplier_dict = {known_reaction_strs[i]: float(solution[multipliers[i]]) for i in range(len(known_reaction_strs))}
    
    return float(total_enthalpy), multiplier_dict

# Add this new function to your energy.py module
def mass_from_specific_heat(q: float, c: float, delta_T: float) -> float:
    """
    Calculates the mass of a substance from the heat exchanged.
    m = q / (c * ΔT)
    """
    if c == 0 or delta_T == 0:
        raise ValueError("Specific heat and delta_T cannot be zero.")
    return q / (c * delta_T)

def gibbs_free_energy(
    delta_H: float,
    delta_S: float,
    T: float,
    S_unit: str = 'J/mol*K',
    T_unit: str = 'C'
) -> tuple[float, str]:
    """
    Calculates Gibbs Free Energy (ΔG) and determines spontaneity.
    ΔG = ΔH - TΔS.
    """
    # 1. Convert units to kJ and K
    T_k = T + 273.15 if T_unit.lower() == 'c' else T
    
    # CORRECTED LINE: Checks if the unit string contains 'j/' (for Joules)
    delta_S_kJ = delta_S / 1000.0 if 'j/' in S_unit.lower() else delta_S

    # 2. Calculate ΔG
    delta_G = delta_H - (T_k * delta_S_kJ)

    # 3. Determine spontaneity
    if delta_G < 0:
        spontaneity = "ΔG is negative, so the reaction is spontaneous in the forward direction."
    elif delta_G > 0:
        spontaneity = "ΔG is positive, so the reaction is spontaneous in the reverse direction."
    else:
        spontaneity = "ΔG is zero, so the reaction is at equilibrium."

    return delta_G, spontaneity

def analyze_reaction_thermodynamics(
    reaction: ChemicalEquation,
    hf_data: Dict[Compound, float],
    s_data: Dict[Compound, float]
) -> str:
    """
    Calculates ΔH° and ΔS° for a reaction and provides a descriptive summary.

    Args:
        reaction: A ChemicalEquation object.
        hf_data: A dictionary mapping compounds to their ΔH_f° values.
        s_data: A dictionary mapping compounds to their S° values.

    Returns:
        A formatted string summarizing the thermodynamic properties.
    """
    try:
        # --- Calculate Enthalpy (ΔH°) ---
        delta_H = reaction_enthalpy(reaction, hf_data)
        
        if delta_H > 0:
            h_summary = f"Endothermic (ΔH° = {delta_H:.1f} kJ/mol)"
        else:
            h_summary = f"Exothermic (ΔH° = {delta_H:.1f} kJ/mol)"

        # --- Calculate Entropy (ΔS°) ---
        delta_S = reaction_entropy(reaction, s_data)

        if delta_S > 0:
            s_summary = f"Entropy Increases (ΔS° = {delta_S:.1f} J/mol·K)"
        else:
            s_summary = f"Entropy Decreases (ΔS° = {delta_S:.1f} J/mol·K)"
            
        return f"{reaction.as_molecular()}\n  - Enthalpy: {h_summary}\n  - Entropy:  {s_summary}"

    except KeyError as e:
        return f"Could not analyze {reaction.as_molecular()}: Missing data for {e}"
    
def entropy_of_phase_change(
    delta_H: float,
    T: float,
    H_unit: str = 'kJ/mol',
    T_unit: str = 'C'
) -> float:
    """
    Calculates the entropy change for a phase transition at constant T.
    ΔS = ΔH / T. Returns ΔS in J/(K*mol).
    """
    # 1. Convert units
    T_k = T + 273.15 if T_unit.lower() == 'c' else T
    delta_H_J = delta_H * 1000.0 if H_unit.lower() == 'kj/mol' else delta_H

    if T_k == 0:
        raise ValueError("Temperature cannot be zero Kelvin.")
    
    # 2. Calculate and return ΔS
    return delta_H_J / T_k

def clausius_clapeyron_pressure(
    P1: float,
    T1_C: float,
    T2_C: float,
    delta_H_vap_kJ: float
) -> float:
    """
    Calculates vapor pressure P2 at temperature T2 using the Clausius-Clapeyron equation.
    
    Args:
        P1: Known pressure (e.g., in mmHg or atm) at T1.
        T1_C: Known temperature in Celsius.
        T2_C: Target temperature in Celsius.
        delta_H_vap_kJ: Enthalpy of vaporization in kJ/mol.
        
    Returns:
        The new pressure P2 in the same units as P1.
    """
    R = 8.314  # J/(mol*K)
    
    # Convert units
    T1_K = T1_C + 273.15
    T2_K = T2_C + 273.15
    delta_H_vap_J = delta_H_vap_kJ * 1000
    
    # Clausius-Clapeyron equation
    ln_P2_over_P1 = (-delta_H_vap_J / R) * (1/T2_K - 1/T1_K)
    
    P2 = P1 * math.exp(ln_P2_over_P1)
    return P2

def clausius_clapeyron_temperature(
    P1: float,
    T1_C: float,
    P2: float,
    delta_H_vap_kJ: float
) -> float:
    """
    Calculates temperature T2 at which a substance boils at pressure P2.
    """
    R = 8.314  # J/(mol*K)
    
    # Convert units
    T1_K = T1_C + 273.15
    delta_H_vap_J = delta_H_vap_kJ * 1000
    
    # Rearranged Clausius-Clapeyron equation
    inv_T2 = (1/T1_K) - (R / delta_H_vap_J) * math.log(P2 / P1)
    
    T2_K = 1 / inv_T2
    
    # Convert back to Celsius and return
    return T2_K - 273.15