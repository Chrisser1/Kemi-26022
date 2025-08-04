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
    """Populate the *PT* registry with all 118 elements.

    Atomic masses come from the IUPAC 2023 Table of Standard Atomic Weights.
    Values for radioactive super‑heavies are the mass number of the most
    stable isotope (shown in parentheses in most references).
    """
    elements = [
        # symbol, name, Z, atomic_mass
        ("H",  "Hydrogen",      1,   1.01),
        ("He", "Helium",        2,   4.00),
        ("Li", "Lithium",       3,   6.94),
        ("Be", "Beryllium",     4,   9.01),
        ("B",  "Boron",         5,  10.81),
        ("C",  "Carbon",        6,  12.01),
        ("N",  "Nitrogen",      7,  14.01),
        ("O",  "Oxygen",        8,  16.00),
        ("F",  "Fluorine",      9,  19.00),
        ("Ne", "Neon",         10,  20.18),
        ("Na", "Sodium",       11,  22.99),
        ("Mg", "Magnesium",    12,  24.31),
        ("Al", "Aluminium",    13,  26.98),
        ("Si", "Silicon",      14,  28.09),
        ("P",  "Phosphorus",   15,  30.97),
        ("S",  "Sulfur",       16,  32.06),
        ("Cl", "Chlorine",     17,  35.45),
        ("Ar", "Argon",        18,  39.95),
        ("K",  "Potassium",    19,  39.10),
        ("Ca", "Calcium",      20,  40.08),
        ("Sc", "Scandium",     21,  44.96),
        ("Ti", "Titanium",     22,  47.87),
        ("V",  "Vanadium",     23,  50.94),
        ("Cr", "Chromium",     24,  52.00),
        ("Mn", "Manganese",    25,  54.94),
        ("Fe", "Iron",         26,  55.85),
        ("Co", "Cobalt",       27,  58.93),
        ("Ni", "Nickel",       28,  58.69),
        ("Cu", "Copper",       29,  63.55),
        ("Zn", "Zinc",         30,  65.38),
        ("Ga", "Gallium",      31,  69.72),
        ("Ge", "Germanium",    32,  72.63),
        ("As", "Arsenic",      33,  74.92),
        ("Se", "Selenium",     34,  78.97),
        ("Br", "Bromine",      35,  79.90),
        ("Kr", "Krypton",      36,  83.80),
        ("Rb", "Rubidium",     37,  85.47),
        ("Sr", "Strontium",    38,  87.62),
        ("Y",  "Yttrium",      39,  88.91),
        ("Zr", "Zirconium",    40,  91.22),
        ("Nb", "Niobium",      41,  92.91),
        ("Mo", "Molybdenum",   42,  95.95),
        ("Tc", "Technetium",   43,  98.0),
        ("Ru", "Ruthenium",    44, 101.07),
        ("Rh", "Rhodium",      45, 102.91),
        ("Pd", "Palladium",    46, 106.42),
        ("Ag", "Silver",       47, 107.87),
        ("Cd", "Cadmium",      48, 112.41),
        ("In", "Indium",       49, 114.82),
        ("Sn", "Tin",          50, 118.71),
        ("Sb", "Antimony",     51, 121.76),
        ("Te", "Tellurium",    52, 127.60),
        ("I",  "Iodine",       53, 126.90),
        ("Xe", "Xenon",        54, 131.29),
        ("Cs", "Caesium",      55, 132.91),
        ("Ba", "Barium",       56, 137.33),
        ("La", "Lanthanum",    57, 138.91),
        ("Ce", "Cerium",       58, 140.12),
        ("Pr", "Praseodymium", 59, 140.91),
        ("Nd", "Neodymium",    60, 144.24),
        ("Pm", "Promethium",   61, 145.0),
        ("Sm", "Samarium",     62, 150.36),
        ("Eu", "Europium",     63, 151.96),
        ("Gd", "Gadolinium",   64, 157.25),
        ("Tb", "Terbium",      65, 158.93),
        ("Dy", "Dysprosium",   66, 162.50),
        ("Ho", "Holmium",      67, 164.93),
        ("Er", "Erbium",       68, 167.26),
        ("Tm", "Thulium",      69, 168.93),
        ("Yb", "Ytterbium",    70, 173.05),
        ("Lu", "Lutetium",     71, 174.97),
        ("Hf", "Hafnium",      72, 178.49),
        ("Ta", "Tantalum",     73, 180.95),
        ("W",  "Tungsten",     74, 183.84),
        ("Re", "Rhenium",      75, 186.21),
        ("Os", "Osmium",       76, 190.23),
        ("Ir", "Iridium",      77, 192.22),
        ("Pt", "Platinum",     78, 195.08),
        ("Au", "Gold",         79, 196.97),
        ("Hg", "Mercury",      80, 200.59),
        ("Tl", "Thallium",     81, 204.38),
        ("Pb", "Lead",         82, 207.2),
        ("Bi", "Bismuth",      83, 208.98),
        ("Po", "Polonium",     84, 209.0),
        ("At", "Astatine",     85, 210.0),
        ("Rn", "Radon",        86, 222.0),
        ("Fr", "Francium",     87, 223.0),
        ("Ra", "Radium",       88, 226.0),
        ("Ac", "Actinium",     89, 227.0),
        ("Th", "Thorium",      90, 232.04),
        ("Pa", "Protactinium", 91, 231.04),
        ("U",  "Uranium",      92, 238.03),
        ("Np", "Neptunium",    93, 237.0),
        ("Pu", "Plutonium",    94, 244.0),
        ("Am", "Americium",    95, 243.0),
        ("Cm", "Curium",       96, 247.0),
        ("Bk", "Berkelium",    97, 247.0),
        ("Cf", "Californium",  98, 251.0),
        ("Es", "Einsteinium",  99, 252.0),
        ("Fm", "Fermium",     100, 257.0),
        ("Md", "Mendelevium", 101, 258.0),
        ("No", "Nobelium",    102, 259.0),
        ("Lr", "Lawrencium",  103, 266.0),
        ("Rf", "Rutherfordium",104, 267.0),
        ("Db", "Dubnium",     105, 268.0),
        ("Sg", "Seaborgium",  106, 269.0),
        ("Bh", "Bohrium",     107, 270.0),
        ("Hs", "Hassium",     108, 270.0),
        ("Mt", "Meitnerium",  109, 278.0),
        ("Ds", "Darmstadtium",110, 281.0),
        ("Rg", "Roentgenium", 111, 282.0),
        ("Cn", "Copernicium", 112, 285.0),
        ("Nh", "Nihonium",    113, 286.0),
        ("Fl", "Flerovium",   114, 289.0),
        ("Mc", "Moscovium",   115, 290.0),
        ("Lv", "Livermorium", 116, 293.0),
        ("Ts", "Tennessine",  117, 294.0),
        ("Og", "Oganesson",   118, 294.0),
    ]

    for sym, name, z, mass in elements:
        PT[sym] = Element(sym, name, z, mass)

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

    def __init__(self,
                 formula_or_comp: Union[str, Dict[str, int]]):
        # accept either a formula string or a pre-built composition dict
        if isinstance(formula_or_comp, str):
            self.input_formula = formula_or_comp
            self.composition = parse_formula(formula_or_comp)
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

        if self.input_formula in self._KNOWN_COMPOUNDS:
            return self._KNOWN_COMPOUNDS[self.input_formula].capitalize()

        if len(comp) == 1:
            # Single element compounds (e.g. "Au")
            el = next(iter(comp))
            if el in PT:
                return PT[el].name.capitalize()
            else:
                raise ValueError(f"Unknown element symbol '{el}' in formula '{self.input_formula}'")

        # 1) Pure covalent: both atoms are in the covalent-nonmetal list
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


        # 2) Polyatomic ionic salts
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

        # 3) Simple binary ionic
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
        print(f"Formula:      {self.formula()}")
        print(f"Molar mass:   {self.molar_mass:.4f} g/mol")
        # show atom count per formula unit
        print(f"Atom count:   {self.element_atom_count()} atoms/unit")
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
