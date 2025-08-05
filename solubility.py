# solubility.py

import re
from typing import Tuple
from compound import Compound

# ──────────────────────────────────────────────────────────────────────────────
# All the INSOLUBLE cation–anion pairs from Figure 4.12 (blue squares)
# ──────────────────────────────────────────────────────────────────────────────
INSOLUBLE_PAIRS = {
    ("Ag", "Cl"), ("Pb", "Cl"), ("Hg", "Cl"),
    ("Ag", "Br"), ("Pb", "Br"), ("Hg", "Br"),
    ("Ag", "I"),  ("Pb", "I"),  ("Hg", "I"),
    ("Ag", "S"),  ("Pb", "S"),  ("Hg", "S"), ("Mg", "S"),
    ("Ag", "SO4"),  ("Pb", "SO4"),  ("Hg", "SO4"), ("Ca", "SO4"), ("Ba", "SO4"),
    ("Ag", "OH"),  ("Pb", "OH"),  ("Hg", "OH"), ("Ca", "OH"), ("Mg", "OH"), ("Fe", "OH"), ("Zn", "OH"), ("Cu", "OH"), ("Al", "OH"),
    ("Ag", "SO3"),  ("Pb", "SO3"),  ("Hg", "SO3"), ("Ca", "SO3"), ("Mg", "SO3"), ("Ba", "SO3"), ("Fe", "SO3"), ("Zn", "SO3"), ("Cu", "SO3"), ("Al", "SO3"),
    ("Ag", "CO3"),  ("Pb", "CO3"),  ("Hg", "CO3"), ("Ca", "CO3"), ("Mg", "CO3"), ("Ba", "CO3"), ("Fe", "CO3"), ("Zn", "CO3"), ("Cu", "CO3"), ("Al", "CO3"),
    ("Ag", "PO4"),  ("Pb", "PO4"),  ("Hg", "PO4"), ("Ca", "PO4"), ("Mg", "PO4"), ("Ba", "PO4"), ("Fe", "PO4"), ("Zn", "PO4"), ("Cu", "PO4"), ("Al", "PO4"), ("Li", "PO4"),
}

# add sulfide, sulfite, carbonate, phosphate insoluble pairs
_special_anions = ["S", "SO3", "CO3", "PO4"]
_cations_for_poly = ["Ag","Mg","Ca","Sr","Ba","Fe","Zn","Cu","Pb","Hg","Al"]
for an in _special_anions:
    for cat in _cations_for_poly:
        INSOLUBLE_PAIRS.add((cat, an))


# ──────────────────────────────────────────────────────────────────────────────
# Splitting formula → (cation, anion)
# ──────────────────────────────────────────────────────────────────────────────
# match polyatomics first, longest first
_POLYATOMS = ["C2H3O2","ClO4","ClO3","NO3","SO4","SO3","CO3","PO4","OH"]
# we'll sort by length so "C2H3O2" beats "CO3", etc.
_POLYATOMS = sorted(_POLYATOMS, key=len, reverse=True)

_RX = re.compile(r"^([A-Z][a-z]?)(.+)$")

def _split_ion_pair(formula: str) -> Tuple[str,str]:
    """
    Return (cation, anion) for a simple ionic formula like "CaCO3" or "K2SO4".
    Polyatomics (CO3, PO4, SO4, etc.) are recognized as units.
    """
    # 1) try polyatomics first
    for poly in _POLYATOMS:
        if formula.endswith(poly):
            cat = re.sub(r"\d+", "", formula[:-len(poly)])
            return cat, poly

    # 2) otherwise single‐element anion
    m = _RX.match(formula)
    if not m:
        raise ValueError(f"Cannot split ionic formula '{formula}'")
    cat, rest = m.group(1), m.group(2)
    # strip any digits from the anion part:
    an = re.sub(r"\d+", "", rest)
    return cat, an

# ──────────────────────────────────────────────────────────────────────────────
# Fallback “high‐school” rules for anion/cation pairs NOT in INSOLUBLE_PAIRS
# ──────────────────────────────────────────────────────────────────────────────
def _default_rules(cat: str, an: str) -> bool:
    # large oxyanions always soluble
    if an in {"NO3","ClO4","ClO3","C2H3O2"}:
        return True
    # halides
    if an in {"Cl","Br","I"}:
        return cat not in {"Ag","Pb","Hg"}
    # sulfates
    if an == "SO4":
        return cat not in {"Ca","Sr","Ba","Pb","Hg"}
    # hydroxides
    if an == "OH":
        return cat in {"Li","Na","K","Rb","Cs","NH4"} or cat == "Ba"
    # carbonates & phosphates
    if an in {"CO3","PO4"}:
        return cat in {"Li","Na","K","Rb","Cs","NH4"}
    # assume soluble
    return True


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────
def is_soluble(c: Compound) -> bool:
    """
    True if the SOLID ionic compound `c` will dissolve in water,
    based first on the exhaustive Figure 4.12 lookup, then falling
    back to the default high-school rules when needed.
    """
    f = c.input_formula or c.formula()
    cat, an = _split_ion_pair(f)
    cat = re.sub(r"\d+", "", cat)   # drop any stoichiometric subscript on cat

    # explicit chart lookup
    if (cat, an) in INSOLUBLE_PAIRS:
        return False
    # if this anion appeared in the chart at all, it must be soluble
    if any(a == an for (_, a) in INSOLUBLE_PAIRS):
        return True
    # fallback
    return _default_rules(cat, an)
