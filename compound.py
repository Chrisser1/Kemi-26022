from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Union
import re


@dataclass(frozen=True)
class Element:
    symbol: str
    name: str
    atomic_number: int
    atomic_mass: float  # g mol⁻¹
    group: Optional[int] = None
    period: Optional[int] = None
    category: Optional[str] = None

    def __repr__(self) -> str:
        return f"<{self.symbol} ({self.atomic_number}): {self.atomic_mass} g/mol>"

PT: Dict[str, Element] = {}

def _load_periodic_table() -> None:
    # symbol, name, Z, atomic_mass, group, period
    elements = [
        ("H",  "Hydrogen",      1,   1.0079,    1, 1),
        ("He", "Helium",        2,   4.0026,   18, 1),
        ("Li", "Lithium",       3,   6.941,     1, 2),
        ("Be", "Beryllium",     4,   9.0122,    2, 2),
        ("B",  "Boron",         5,  10.81,     13, 2),
        ("C",  "Carbon",        6,  12.011,    14, 2),
        ("N",  "Nitrogen",      7,  14.007,    15, 2),
        ("O",  "Oxygen",        8,  15.999,    16, 2),
        ("F",  "Fluorine",      9,  18.998,    17, 2),
        ("Ne", "Neon",         10,  20.180,    18, 2),
        ("Na", "Sodium",       11,  22.990,     1, 3),
        ("Mg", "Magnesium",    12,  24.305,     2, 3),
        ("Al", "Aluminium",    13,  26.982,    13, 3),
        ("Si", "Silicon",      14,  28.086,    14, 3),
        ("P",  "Phosphorus",   15,  30.974,    15, 3),
        ("S",  "Sulfur",       16,  32.06,     16, 3),
        ("Cl", "Chlorine",     17,  35.45,     17, 3),
        ("Ar", "Argon",        18,  39.948,    18, 3),
        ("K",  "Potassium",    19,  39.098,     1, 4),
        ("Ca", "Calcium",      20,  40.078,     2, 4),
        ("Sc", "Scandium",     21,  44.956,     3, 4),
        ("Ti", "Titanium",     22,  47.867,     4, 4),
        ("V",  "Vanadium",     23,  50.942,     5, 4),
        ("Cr", "Chromium",     24,  51.996,     6, 4),
        ("Mn", "Manganese",    25,  54.938,     7, 4),
        ("Fe", "Iron",         26,  55.845,     8, 4),
        ("Co", "Cobalt",       27,  58.933,     9, 4),
        ("Ni", "Nickel",       28,  58.693,    10, 4),
        ("Cu", "Copper",       29,  63.546,    11, 4),
        ("Zn", "Zinc",         30,  65.38,     12, 4),
        ("Ga", "Gallium",      31,  69.723,    13, 4),
        ("Ge", "Germanium",    32,  72.63,     14, 4),
        ("As", "Arsenic",      33,  74.922,    15, 4),
        ("Se", "Selenium",     34,  78.971,    16, 4),
        ("Br", "Bromine",      35,  79.904,    17, 4),
        ("Kr", "Krypton",      36,  83.798,    18, 4),
        ("Rb", "Rubidium",     37,  85.468,     1, 5),
        ("Sr", "Strontium",    38,  87.62,      2, 5),
        ("Y",  "Yttrium",      39,  88.906,     3, 5),
        ("Zr", "Zirconium",    40,  91.224,     4, 5),
        ("Nb", "Niobium",      41,  92.906,     5, 5),
        ("Mo", "Molybdenum",   42,  95.95,      6, 5),
        ("Tc", "Technetium",   43,  98.0,       7, 5),
        ("Ru", "Ruthenium",    44, 101.07,      8, 5),
        ("Rh", "Rhodium",      45, 102.91,      9, 5),
        ("Pd", "Palladium",    46, 106.42,     10, 5),
        ("Ag", "Silver",       47, 107.868,    11, 5),
        ("Cd", "Cadmium",      48, 112.414,    12, 5),
        ("In", "Indium",       49, 114.818,    13, 5),
        ("Sn", "Tin",          50, 118.71,     14, 5),
        ("Sb", "Antimony",     51, 121.76,     15, 5),
        ("Te", "Tellurium",    52, 127.60,     16, 5),
        ("I",  "Iodine",       53, 126.904,    17, 5),
        ("Xe", "Xenon",        54, 131.293,    18, 5),
        ("Cs", "Caesium",      55, 132.905,     1, 6),
        ("Ba", "Barium",       56, 137.327,     2, 6),
        ("La", "Lanthanum",    57, 138.905,     3, 9),  # La–Lu in period 6 lanthanoid row
        ("Ce", "Cerium",       58, 140.116,     4, 9),
        ("Pr", "Praseodymium", 59, 140.908,     5, 9),
        ("Nd", "Neodymium",    60, 144.242,     6, 9),
        ("Pm", "Promethium",   61, 145.0,       7, 9),
        ("Sm", "Samarium",     62, 150.36,      8, 9),
        ("Eu", "Europium",     63, 151.964,     9, 9),
        ("Gd", "Gadolinium",   64, 157.25,     10, 9),
        ("Tb", "Terbium",      65, 158.925,    11, 9),
        ("Dy", "Dysprosium",   66, 162.500,    12, 9),
        ("Ho", "Holmium",      67, 164.930,    13, 9),
        ("Er", "Erbium",       68, 167.259,    14, 9),
        ("Tm", "Thulium",      69, 168.934,    15, 9),
        ("Yb", "Ytterbium",    70, 173.045,    16, 9),
        ("Lu", "Lutetium",     71, 174.967,    17, 9),
        ("Hf", "Hafnium",      72, 178.49,      4, 6),
        ("Ta", "Tantalum",     73, 180.948,     5, 6),
        ("W",  "Tungsten",     74, 183.84,      6, 6),
        ("Re", "Rhenium",      75, 186.207,     7, 6),
        ("Os", "Osmium",       76, 190.23,      8, 6),
        ("Ir", "Iridium",      77, 192.217,     9, 6),
        ("Pt", "Platinum",     78, 195.084,    10, 6),
        ("Au", "Gold",         79, 196.966,    11, 6),
        ("Hg", "Mercury",      80, 200.592,    12, 6),
        ("Tl", "Thallium",     81, 204.38,     13, 6),
        ("Pb", "Lead",         82, 207.2,      14, 6),
        ("Bi", "Bismuth",      83, 208.980,    15, 6),
        ("Po", "Polonium",     84, 209.0,      16, 6),
        ("At", "Astatine",     85, 210.0,      17, 6),
        ("Rn", "Radon",        86, 222.0,      18, 6),
        ("Fr", "Francium",     87, 223.0,       1, 7),
        ("Ra", "Radium",       88, 226.0,       2, 7),
        ("Ac", "Actinium",     89, 227.0,       3,10),  # Ac–Lr in period 7 actinoid row
        ("Th", "Thorium",      90, 232.038,     4,10),
        ("Pa", "Protactinium", 91, 231.036,     5,10),
        ("U",  "Uranium",      92, 238.029,     6,10),
        ("Np", "Neptunium",    93, 237.0,       7,10),
        ("Pu", "Plutonium",    94, 244.0,       8,10),
        ("Am", "Americium",    95, 243.0,       9,10),
        ("Cm", "Curium",       96, 247.0,      10,10),
        ("Bk", "Berkelium",    97, 247.0,      11,10),
        ("Cf", "Californium",  98, 251.0,      12,10),
        ("Es", "Einsteinium",  99, 252.0,      13,10),
        ("Fm", "Fermium",     100, 257.0,      14,10),
        ("Md", "Mendelevium",101, 258.0,       15,10),
        ("No", "Nobelium",   102, 259.0,       16,10),
        ("Lr", "Lawrencium",103, 266.0,        17,10),
        ("Rf", "Rutherfordium",104, 267.0,     4, 7),
        ("Db", "Dubnium",     105, 268.0,      5, 7),
        ("Sg", "Seaborgium",  106, 269.0,      6, 7),
        ("Bh", "Bohrium",     107, 270.0,      7, 7),
        ("Hs", "Hassium",     108, 270.0,      8, 7),
        ("Mt", "Meitnerium",  109, 278.0,      9, 7),
        ("Ds", "Darmstadtium",110, 281.0,     10, 7),
        ("Rg", "Roentgenium", 111, 282.0,     11, 7),
        ("Cn", "Copernicium", 112, 285.0,     12, 7),
        ("Nh", "Nihonium",    113, 286.0,     13, 7),
        ("Fl", "Flerovium",   114, 289.0,     14, 7),
        ("Mc", "Moscovium",   115, 290.0,     15, 7),
        ("Lv", "Livermorium", 116, 293.0,     16, 7),
        ("Ts", "Tennessine",  117, 294.0,     17, 7),
        ("Og", "Oganesson",   118, 294.0,     18, 7),
    ]

    for sym, name, z, mass, grp, prd in elements:
        PT[sym] = Element(sym, name, z, mass, group=grp, period=prd)

_load_periodic_table()

# ──────────────────────────────────────────────────────────────────────────────
# FORMULA PARSING
# ──────────────────────────────────────────────────────────────────────────────

_token_re = re.compile(r"([A-Z][a-z]?|\d+|\(|\))")

def _parse_formula_tokens(tokens: list[str]) -> Dict[str, int]:
    """
    Recursive descent helper used by :func:`parse_formula`.
    """
    comp: Dict[str, int] = {}

    i = 0
    last_el: Optional[str] = None
    while i < len(tokens):
        tok = tokens[i]
        if tok == '(':
            # Find matching ')' and recurse
            depth = 1
            j = i + 1
            while j < len(tokens) and depth:
                if tokens[j] == '(': depth += 1
                elif tokens[j] == ')': depth -= 1
                j += 1
            if depth != 0:
                raise ValueError("Unbalanced parentheses in formula")
            inner = _parse_formula_tokens(tokens[i + 1 : j - 1])
            multiplier = 1
            if j < len(tokens) and tokens[j].isdigit():
                multiplier = int(tokens[j]); j += 1
            for el, cnt in inner.items():
                comp[el] = comp.get(el, 0) + cnt * multiplier
            i = j
            last_el = None
            continue
        elif tok.isdigit():
            if last_el is None:
                raise ValueError("Unexpected number in formula")
            comp[last_el] += int(tok) - 1  # we already counted 1
            i += 1
            continue
        else:  # element symbol
            if tok not in PT:
                raise ValueError(f"Unknown element symbol '{tok}'")
            comp[tok] = comp.get(tok, 0) + 1
            last_el = tok
            i += 1
    return comp


def parse_formula(formula: str) -> Dict[str, int]:
    """Convert a chemical formula string (e.g. "C6H12O6") → composition dict."""
    tokens = _token_re.findall(formula)
    if not tokens:
        raise ValueError("Empty formula")
    return _parse_formula_tokens(tokens)


class Compound:
    """Represents a chemical compound with optional mass / mol information."""

    _COVALENT_NONMETALS = {'H','C','N','O','F','P','S','Cl','Se','Br','I'}

    # for bond‐type classification:
    _NONMETALS = _COVALENT_NONMETALS | {"He","Ne","Ar","Kr","Xe","Rn","Og"}

    _POLYATOMIC_IONS = {
    'NH4': 'ammonium',
    'OH': 'hydroxide',
    'NO3': 'nitrate',
    'SO4': 'sulfate',
    'CO3': 'carbonate',
    'PO4': 'phosphate',
    'C2H3O2': 'acetate',
    'CN': 'cyanide',
    'HCO3': 'bicarbonate',
    'ClO3': 'chlorate',
    'ClO4': 'perchlorate',
    'Cr2O7': 'dichromate',
    'MnO4': 'permanganate',
    'MgCO3': 'magnesium carbonate',
    'SO3': 'sulfite',
    'CrO4': 'chromate',
    }

    _KNOWN_COMPOUNDS = {
        'HCl':    'hydrochloric acid',
        'HNO3':   'nitric acid',
        'H2SO4':  'sulfuric acid',
        'H3PO4':  'phosphoric acid',
        'CH3COOH':'acetic acid',
        'NH3':    'ammonia',
        'H2O':    'water',
        'O3':     'ozone',
        'CO2':    'carbon dioxide',
        'NO3':    'nitrate',
        'NaCl':   'sodium chloride',
        'KCl':    'potassium chloride',
        'CaCO3':  'calcium carbonate',
        'MgSO4':  'magnesium sulfate',
        # …add any others you want to recognize…
    }

    _ANION_CHARGES = {
        'NH4': +1,
        'OH': -1,
        'NO3': -1,
        'SO4': -2,
        'CO3': -2,
        'PO4': -3,
        'C2H3O2': -1,
        'CN': -1,
        'HCO3': -1,
        'SO3': -2,
        'CrO4': -2,
        'ClO3': -1,
        'ClO4': -1,
        # …add more if you like
    }

    # Pauling electronegativity (add more as needed)
    _EN: dict[str,float] = {
        "H":2.20,"C":2.55,"N":3.04,"O":3.44,"F":3.98,
        "P":2.19,"S":2.58,"Cl":3.16,"Br":2.96,"I":2.66,
        "Na":0.93,"Mg":1.31,"Al":1.61,"Si":1.90,"K":0.82,
        "Ca":1.00,"Fe":1.83,"Cu":1.90,"Zn":1.65,"Ag":1.93,
        # …extend for any others you’ll encounter…
    }

    # regex to pull off a trailing ionic charge, e.g. "Fe2+", "Cl-", "SO4^2-"
    _ion_re = re.compile(r'^(?P<formula>[A-Za-z0-9()]+)'
                         r'(?P<charge>⁺⁻0-9\+\-]+)$')

    phase: Optional[str] = None  # 's', 'l', 'g', 'aq' (default: None)
    charge: Optional[int] = None  # e.g. +2, -1, etc.

    def __init__(self, formula_or_comp: Union[str, Dict[str, int]]):
        # accept either a formula string (possibly with charge) or a dict
        self.input_formula = None
        self.charge: Optional[int] = None

        # accept either a formula string or a pre-built composition dict
        if isinstance(formula_or_comp, str):
            inp = formula_or_comp.strip()
            m = Compound._ion_re.match(inp)
            if m:
                # split out the formula part from the charge part
                raw, ch = m.group("formula"), m.group("charge")
                # unify +/− signs, parse into integer
                ch_norm = ch.replace('⁺','+').replace('⁻','-')
                # handle e.g. "+", "2+" or "-","3-"
                if ch_norm in ('+','-'):
                    val = 1 if ch_norm=='+' else -1
                else:
                    # e.g. "2+", "3-"
                    num, sign = re.match(r'(\d+)([+\-])', ch_norm).groups()
                    val = int(num) * (1 if sign=='+' else -1)
                self.charge = val
                inp = raw
            self.input_formula = inp
            self.composition = parse_formula(inp)
        else:
            self.composition = dict(formula_or_comp)
        self.mass_g: Optional[float] = None
        self.amount_mol: Optional[float] = None

    @classmethod
    def from_formula(cls,
                     formula: str,
                     name: Optional[str] = None) -> "Compound":
        return cls(formula, name)

    @property
    def molar_mass(self) -> float:
        """Return M, the molar mass in g mol⁻¹."""
        return sum(PT[el].atomic_mass * n for el, n in self.composition.items())

    @property
    def bond_types(self) -> list[str]:
        """Guess bond types: metallic, ionic, covalent (polar/nonpolar), H-bonding."""
        elems = set(self.composition)
        nonm = { *self._COVALENT_NONMETALS, "He","Ne","Ar","Kr","Xe","Rn","Og" }

        types: list[str] = []

        # single‐element “compound” (O2, Fe, etc.)
        if len(elems)==1:
            el = next(iter(elems))
            if el not in nonm:
                types.append("Metallic")
            else:
                types.append("Nonpolar covalent")
            # still check H-bonding if it’s H2?
            if el=="H":
                types.append("Hydrogen-bond donor")
            return types

        # Ionic: metal + nonmetal
        if elems & nonm and elems - nonm:
            types.append("Ionic")
        # Metallic clusters (all metal)
        elif not (elems & nonm):
            types.append("Metallic")
        # Covalent (all nonmetal)
        else:
            # estimate polarity by max EN diff among bonded pairs
            maxdiff = 0.0
            for a in elems:
                for b in elems:
                    if a==b: continue
                    en_a = Compound._EN.get(a)
                    en_b = Compound._EN.get(b)
                    if en_a and en_b:
                        maxdiff = max(maxdiff, abs(en_a-en_b))
            if maxdiff < 0.5:
                types.append("Nonpolar covalent")
            else:
                types.append("Polar covalent")

        # Hydrogen‐bond capability?
        if "H" in elems and elems & {"N","O","F"}:
            types.append("Hydrogen-bond donor/acceptor")

        return types

    @staticmethod
    def _anion_name(el_name: str) -> str:
        """Convert a neutral element name into its anion form (–ide)."""
        exceptions = {
            "Oxygen": "oxide",
            "Sulfur": "sulfide",
            "Nitrogen": "nitride",
            "Phosphorus": "phosphide",
            "Carbon": "carbide",
            "Hydrogen": "hydride",
        }
        if el_name in exceptions:
            return exceptions[el_name]
        if el_name.endswith("ine"):
            return el_name[:-3] + "ide"
        return el_name + "ide"

    @staticmethod
    def _to_roman(num: int) -> str:
        """Convert an integer to a Roman numeral (supports 1–20)."""
        vals = [(10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')]
        res = ''
        for val, sym in vals:
            while num >= val:
                res += sym
                num -= val
        return res

    @staticmethod
    def _element_type(symbol: str) -> Optional[str]:
        el = PT[symbol]
        z, g = el.atomic_number, el.group

        # special cases
        if z == 1:
            return "nonmetal"

        # f-block
        if 57 <= z <= 71:
            return "lanthanoid"
        if 89 <= z <= 103:
            return "actinoid"

        # s-block
        if g == 1:
            return "alkali metal"
        if g == 2:
            return "alkaline earth metal"

        # d-block
        if 3 <= g <= 12:
            return "transition metal"

        # p-block metals
        post = {"Al","Ga","In","Sn","Tl","Pb","Bi","Nh","Fl","Mc","Lv"}
        if symbol in post:
            return "post-transition metal"

        # metalloids
        met = {"B","Si","Ge","As","Sb","Te","Po"}
        if symbol in met:
            return "metalloid"

        # remaining p-block nonmetals
        if symbol in {"C","N","O","P","S","Se"}:
            return "nonmetal"

        # group-based nonmetals
        if g == 17:
            return "halogen"
        if g == 18:
            return "noble gas"

        return None

    def set_mass(self, grams: float) -> None:
        self.mass_g = grams
        self.amount_mol = grams / self.molar_mass

    def set_moles(self, moles: float) -> None:
        self.amount_mol = moles
        self.mass_g = moles * self.molar_mass

    def element_atom_count(self) -> int:
        """Return total number of atoms in one formula unit."""
        return sum(self.composition.values())

    def total_formula_units(self) -> float:
        """Return number of formula units (mol × Avogadro's number)."""
        if self.amount_mol is None:
            raise ValueError("Set mass or moles first.")
        return self.amount_mol * 6.02214076e23

    def total_atoms(self) -> float:
        """Return total atoms in the current sample."""
        return self.total_formula_units() * self.element_atom_count()

    def element_moles(self) -> Dict[str, float]:
        """Return a dict of moles of each element in this compound."""
        if self.amount_mol is None:
            raise ValueError("You must set mass or moles before querying element amounts.")
        return {el: n * self.amount_mol for el, n in self.composition.items()}

    def element_mole_percent(self) -> Dict[str, float]:
        """Return a dict of mole percentages of each element in this compound."""
        if self.amount_mol is None:
            raise ValueError("You must set mass or moles before querying element amounts.")
        total = sum(self.composition.values()) * self.amount_mol
        # Alternatively, percent is simply based on mole fractions:
        return {
            el: (count * self.amount_mol / self.amount_mol) * 100 / sum(self.composition.values())
            for el, count in self.composition.items()
        }

    def element_molar_contributions(self) -> Dict[str, float]:
        """Return each element’s g/mol contribution plus the total molar mass of the compound."""
        contributions = {el: PT[el].atomic_mass * n
                        for el, n in self.composition.items()}
        contributions['total'] = self.molar_mass
        return contributions

    def element_mass_contributions(self) -> Dict[str, float]:
        """Return each element’s mass contribution in grams."""
        if self.mass_g is None:
            raise ValueError("You must set mass or moles before querying element contributions.")
        return {
            el: (PT[el].atomic_mass * count / self.molar_mass) * self.mass_g
            for el, count in self.composition.items()
        }

    def element_mass_percent(self) -> Dict[str, float]:
        """Return a dict of mass percentages of each element in this compound."""
        total_molar = self.molar_mass
        if total_molar == 0:
            raise ValueError("Compound has zero molar mass")
        return {
            el: (PT[el].atomic_mass * count) / total_molar * 100
            for el, count in self.composition.items()
        }


    def formula(self) -> str:
        """Return a Hill-system formatted formula string."""
        # Hill: C first, then H, then all others alphabetically
        elements = sorted(
            self.composition.keys(),
            key=lambda el: (0, el) if el == "C"
                        else (1, el) if el == "H"
                        else (2, el)
        )
        parts = [
            f"{el}{self.composition[el] if self.composition[el] > 1 else ''}"
            for el in elements
        ]
        return "".join(parts)


    def name(self) -> str:
        """Generate an IUPAC‐style name for covalent molecules and ionic salts."""
        comp = self.composition
        elems = list(comp.items())  # [(symbol, count), ...]

        if len(comp)==1 and self.charge is not None:
            el = next(iter(comp))
            base = PT[el].name.capitalize()
            z = abs(self.charge)
            # cation: e.g. Fe2+ → Iron(II)
            if self.charge>0:
                return f"{base}({self._to_roman(z)})"
            # anion: e.g. Cl− → Chloride
            else:
                root = self._anion_name(PT[el].name).capitalize()
                return root

        # Otherwise, if it’s a neutral element (no charge) and you still want types:
        if len(comp)==1 and self.charge is None:
            el = next(iter(comp))
            base = PT[el].name.capitalize()
            et = self._element_type(el)
            return f"{base} ({et})" if et else base

        if self.input_formula in self._KNOWN_COMPOUNDS:
            return self._KNOWN_COMPOUNDS[self.input_formula].capitalize()

        # Pure covalent: both atoms are in the covalent-nonmetal list
        if (len(comp) == 2
            and all(el in self._COVALENT_NONMETALS for el in comp)):
            elems = list(comp.items())
            prefixes = {
                1: "mono", 2: "di", 3: "tri", 4: "tetra",
                5: "penta", 6: "hexa", 7: "hepta", 8: "octa",
                9: "nona", 10: "deca"
            }
            parts = []
            for i, (el, cnt) in enumerate(elems):
                name_low = PT[el].name.lower()
                pre = prefixes.get(cnt, str(cnt))
                if i == 0:
                    parts.append(name_low if cnt == 1 else f"{pre}{name_low}")
                else:
                    root = self._anion_name(PT[el].name).lower()
                    parts.append(f"{pre}{root}")
            return " ".join(parts).capitalize()


        # Polyatomic ionic salts
        for formula, ion_name in self._POLYATOMIC_IONS.items():
            anion_comp = parse_formula(formula)
            if all(comp.get(e,0) >= n for e,n in anion_comp.items()):
                # leftover is the cation
                cat = {e: comp[e] for e in comp if e not in anion_comp}
                if len(cat) == 1:
                    cation, n_cat = next(iter(cat.items()))
                    n_an = sum(anion_comp.values())
                    charge_an = self._ANION_CHARGES.get(formula)
                    if charge_an is None:
                        raise ValueError(f"No charge for {formula}")
                    ox_cat = -(n_an*charge_an)//n_cat
                    return f"{PT[cation].name}({self._to_roman(ox_cat)}) {ion_name}"

        # Simple binary ionic
        if len(elems) == 2:
            (cat, n_cat), (an, n_an) = elems
            # guess anion charge from a small lookup or group‐rule
            ANION_OX = {'O':-2,'S':-2,'N':-3,'P':-3,'F':-1,'Cl':-1,'Br':-1,'I':-1}
            an_ox = ANION_OX.get(an)
            if an_ox is None:
                grp = PT[an].group
                if grp is None:
                    raise ValueError(f"No oxidation for {an}")
                an_ox = grp - 18
            ox_cat = -(n_an*an_ox)//n_cat
            return f"{PT[cat].name}({self._to_roman(ox_cat)}) {self._anion_name(PT[an].name)}"

        raise NotImplementedError("Cannot generate name for this compound.")


    def display(self) -> None:
        """Pretty-print all available info about this compound, including generated name and atom count."""
        print(f"Name:         {self.name()}")
        print(f"Bond type(s): {', '.join(self.bond_types)}")
        print(f"Formula:      {self.formula()}")
        print(f"Molar mass:   {self.molar_mass:.4f} g/mol")
        # show atom count per formula unit
        print(f"Atom count:   {self.element_atom_count()} atoms/unit")
        if self.amount_mol is not None:
            print(f"Total units:  {self.total_formula_units():.2e} formula units"
              f" ({self.total_atoms():.2e} atoms in total)")
        if self.mass_g is not None:
            print(f"Mass:         {self.mass_g:.4f} g")
        if self.amount_mol is not None:
            print(f"Amount:       {self.amount_mol:.6f} mol")


    def __repr__(self) -> str:
        extra = []
        if self.mass_g is not None:
            extra.append(f"{self.mass_g} g")
        if self.amount_mol is not None:
            extra.append(f"{self.amount_mol} mol")
        extras = ", ".join(extra)
        return f"<Compound {self.formula()} ({extras})>"
