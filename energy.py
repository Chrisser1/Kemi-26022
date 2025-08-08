# energy.py

import pandas as pd
import re
from typing import Dict, Sequence, Tuple
from compound import Compound
from equations import ChemicalEquation

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
