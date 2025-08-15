from itertools import product
from sympy import Matrix
from typing import Dict, List, Tuple
from math import gcd, inf
from functools import reduce
from collections import Counter

from compound import Compound, PT, parse_formula


class ChemicalEquation:
    """
    Represent and balance a chemical equation using sympy's nullspace.
    Reactants and products are lists of Compound instances (which may have a .phase attribute).

    Includes automatic classification into standard reaction types:
      - Synthesis (Combination)
      - Decomposition
      - Single Displacement
      - Double Displacement
      - Combustion
      - Oxidation-Reduction
    Also supports molecular, ionic, and net ionic representations.
    """


    def __init__(self,
                 reactants: List[Compound],
                 products:  List[Compound]):
        self.reactants = reactants
        self.products  = products

    def _element_list(self) -> List[str]:
        # collect all unique element symbols in order
        elems: List[str] = []
        for cmpd in self.reactants + self.products:
            for el in cmpd.composition:
                if el not in elems:
                    elems.append(el)
        return elems

    def _validate(self):
        # ensure no new elements appear in products
        elems_in  = {el for cmpd in self.reactants for el in cmpd.composition}
        elems_out = {el for cmpd in self.products  for el in cmpd.composition}
        extra = elems_out - elems_in
        if extra:
            raise ValueError(f"Cannot produce {extra} on product side when absent from reactants")

    @staticmethod
    def _split_ions(cmpd: Compound, coeff: int):
        """
        Splits an aqueous compound into its constituent ions with correct charges.
        """
        # If it's not an aqueous compound or already has a defined charge, do not split.
        if cmpd.phase != "aq" or cmpd.charge is not None:
            return [cmpd] * coeff

        comp = cmpd.composition.copy()
        ions: list[Compound] = []

        # 1. Handle known polyatomic ions first.
        # (This part of the logic remains the same)
        poly_found = False
        for poly in sorted(Compound._POLYATOMIC_IONS, key=len, reverse=True):
            needed = parse_formula(poly)
            if all(comp.get(el, 0) >= n for el, n in needed.items()):
                poly_found = True
                count = min(comp.get(el, 0) // n for el, n in needed.items())
                for _ in range(count):
                    ion = Compound(poly)
                    ion.charge = Compound._ANION_CHARGES.get(poly)
                    ion.phase = "aq"
                    ions.append(ion)
                    for el, n in needed.items():
                        comp[el] -= n
                        if comp[el] == 0:
                            del comp[el]
        
        # 2. Handle the remaining atoms (either the cation for a polyatomic salt or a full binary salt).
        if not poly_found:
            # Logic for simple binary ionic compounds (e.g., SnCl2, FeCl2)
            cation_part = {el: n for el, n in comp.items() if el not in Compound._NONMETALS}
            anion_part = {el: n for el, n in comp.items() if el in Compound._NONMETALS}

            # Proceed if we have one metal and one non-metal part
            if len(cation_part) == 1 and len(anion_part) == 1:
                cation_sym, cation_count = list(cation_part.items())[0]
                anion_sym, anion_count = list(anion_part.items())[0]

                # Determine anion charge from common rules (halogens are -1, group 16 are -2)
                anion_charge = -1 if PT[anion_sym].group == 17 else -2

                # Calculate total negative charge from the anion(s)
                total_negative_charge = anion_count * anion_charge

                # The total positive charge must balance this out
                cation_charge = abs(total_negative_charge) // cation_count

                # Add the correctly charged cation ions
                for _ in range(cation_count):
                    ion = Compound(cation_sym); ion.charge = cation_charge; ion.phase = "aq"
                    ions.append(ion)

                # Add the correctly charged anion ions
                for _ in range(anion_count):
                    ion = Compound(anion_sym); ion.charge = anion_charge; ion.phase = "aq"
                    ions.append(ion)
        
        # This logic handles what's left over after a polyatomic is removed
        elif comp: 
            total_anion_charge = sum(i.charge for i in ions if i.charge < 0)
            cation_atoms = sum(comp.values())
            cation_charge = abs(total_anion_charge) // cation_atoms if cation_atoms > 0 else 0
            for el, n in comp.items():
                for _ in range(n):
                    ion = Compound(el); ion.charge = cation_charge; ion.phase = "aq"
                    ions.append(ion)

        return ions * coeff

    def balance(self) -> Tuple[List[int], List[int]]:
        """
        Return minimal positive integer coefficients (reactants, products).
        """
        self._validate()
        elems = self._element_list()
        n_r = len(self.reactants)
        n_p = len(self.products)

        # build stochiometric matrix (+ for RHS, - for LHS)
        M = [[self.reactants[i].composition.get(el, 0) for i in range(n_r)] +
             [-self.products[j].composition.get(el, 0) for j in range(n_p)]
             for el in elems]

        mat = Matrix(M)
        nulls = mat.nullspace()
        if not nulls:
            raise ValueError("Cannot balance equation: no nullspace basis.")

        if len(nulls) == 1:
            vec = nulls[0]
        else:
            for a,b in product(range(1,10), repeat=2):
                v = a*nulls[0] + b*nulls[1]
                if all(x>0 for x in v):
                    vec = v
                    break

        # clear denominators
        dens = [v.q for v in vec]
        def lcm(a: int, b: int) -> int:
            return a * b // gcd(a, b)
        L = reduce(lcm, dens, 1)
        coeffs = [int(v * L) for v in vec]
        # force all positive
        if any(x < 0 for x in coeffs):
            coeffs = [-x for x in coeffs]

        return coeffs[:n_r], coeffs[n_r:]

    def reaction_type(self) -> str:
        """
        Classify the reaction into standard types, including detecting
        redox when the same element’s ionic charge changes.
        """
        n_r, n_p = len(self.reactants), len(self.products)

        # 1) Build ionic-list for the LHS (all reactants) with their balanced coefficients
        r_coeffs, _ = self.balance()
        ionic_lhs = []
        for cmpd, coeff in zip(self.reactants, r_coeffs):
            ionic_lhs += self._split_ions(cmpd, coeff)

        # 2) Acid–base (neutralization) check:
        #    look at input_formula (not formula()), so that polyatomic OH comes through as "OH"
        has_Hp     = any(i.input_formula=="H"  and i.charge==+1 for i in ionic_lhs)
        has_OHm    = any(i.input_formula=="OH" and i.charge==-1 for i in ionic_lhs)
        makes_water = any(p.phase=="l" and p.formula()=="H2O" for p in self.products)
        if has_Hp and has_OHm and makes_water:
            return "Syre‐base‐reaktion"

        # Precipitation: two aqueous reactants → at least one solid product
        if (n_r == 2
            and all(c.phase=="aq" for c in self.reactants)
            and any(p.phase=="s" for p in self.products)):
            return "Precipitation (Double Displacement)"   # or "Precipitation (Double Displacement)"

        # Check for simple redox: same element appears as an ion on both
        # sides with different charges.
        #    Build maps element→set of charges
        from collections import defaultdict
        charges_l = defaultdict(set)
        charges_r = defaultdict(set)
        for cmpd in self.reactants:
            if cmpd.phase == "aq" and cmpd.charge is not None:
                # assume single-element ion: e.g. Fe2+
                sym = next(iter(cmpd.composition))
                charges_l[sym].add(cmpd.charge)
        for cmpd in self.products:
            if cmpd.phase == "aq" and cmpd.charge is not None:
                sym = next(iter(cmpd.composition))
                charges_r[sym].add(cmpd.charge)

        for el in charges_l:
            if el in charges_r:
                # if the charge sets differ, it’s redox
                if charges_l[el] != charges_r[el]:
                    return "Oxidation-Reduction"

        # Synthesis
        if n_r > 1 and n_p == 1:
            return "Synthesis (Combination)"
        # Decomposition
        if n_r == 1 and n_p > 1:
            return "Decomposition"

        # Combustion
        forms_r = [c.formula() for c in self.reactants]
        forms_p = [c.formula() for c in self.products]
        if "O2" in forms_r and any(set(c.composition)<={"C","H"} for c in self.reactants):
            if "CO2" in forms_p and "H2O" in forms_p:
                return "Combustion"

        # Single vs Double Displacement
        if n_r == 2 and n_p == 2:
            def is_elem(c):
                return len(c.composition)==1 and next(iter(c.composition.values()))==1
            reac_e = sum(is_elem(c) for c in self.reactants)
            prod_e = sum(is_elem(c) for c in self.products)
            if reac_e==1 and prod_e==1:
                return "Single Displacement"
            if reac_e==0 and prod_e==0:
                return "Double Displacement"

        return "Other"


    def as_molecular(self) -> str:
        """Return the balanced molecular equation."""
        r_coeffs, p_coeffs = self.balance()
        return self._fmt_equation(self.reactants, self.products,
                                  r_coeffs, p_coeffs)

    def as_ionic(self) -> str:
        """
        Return the **complete ionic** equation: every (aq) compound
        split into its proper ions, with their true charges.
        """

        # 3) build the flat lists of all ions (with multiplicity)
        r_coeffs, p_coeffs = self.balance()
        flat_lhs = []
        flat_rhs = []
        for cmpd, c in zip(self.reactants, r_coeffs):
            flat_lhs += self._split_ions(cmpd, c)
        for cmpd, c in zip(self.products, p_coeffs):
            flat_rhs += self._split_ions(cmpd, c)

        # 4) tally up identical ions for their stoichiometric coefficients
        def key(i: Compound):
            return (i.formula(), i.charge, i.phase)

        lhs_cnt = Counter(key(i) for i in flat_lhs)
        rhs_cnt = Counter(key(i) for i in flat_rhs)

        # 5) re-assemble in order, with counts
        lhs_uniques, lhs_coefs = [], []
        seen = set()
        for ion in flat_lhs:
            k = key(ion)
            if k in seen: continue
            seen.add(k)
            lhs_uniques.append(ion)
            lhs_coefs.append(lhs_cnt[k])

        rhs_uniques, rhs_coefs = [], []
        seen = set()
        for ion in flat_rhs:
            k = key(ion)
            if k in seen: continue
            seen.add(k)
            rhs_uniques.append(ion)
            rhs_coefs.append(rhs_cnt[k])

        return self._fmt_equation(lhs_uniques, rhs_uniques,
                                  lhs_coefs,  rhs_coefs)

    def as_net_ionic(self) -> str:
        """
        Remove spectators: any ion appearing on both sides cancels,
        *then* group identical ions to show stoichiometric coefficients.
        """

        # 1) build the full ionic lists exactly as in as_ionic (but WITHOUT grouping)
        def split_all(cmpd, coeff):
            if cmpd.phase != "aq" or cmpd.charge is not None:
                return [cmpd] * coeff
            return self._split_ions(cmpd, coeff)

        r_coeffs, p_coeffs = self.balance()
        flat_lhs = sum((split_all(c,n) for c,n in zip(self.reactants, r_coeffs)), [])
        flat_rhs = sum((split_all(c,n) for c,n in zip(self.products,  p_coeffs)), [])

        # 2) tally them
        def key(i): return (i.formula(), i.charge, i.phase)
        lhs_cnt = Counter(key(i) for i in flat_lhs)
        rhs_cnt = Counter(key(i) for i in flat_rhs)

        # 3) cancel spectators
        for k in set(lhs_cnt) & set(rhs_cnt):
            m = min(lhs_cnt[k], rhs_cnt[k])
            lhs_cnt[k] -= m
            rhs_cnt[k] -= m

        # 4) rebuild uniques+counts
        lhs_uniques, lhs_coefs = [], []
        for k,count in lhs_cnt.items():
            if count>0:
                fm, ch, ph = k
                ion = Compound(fm)
                ion.charge = ch
                ion.phase  = ph
                lhs_uniques.append(ion)
                lhs_coefs.append(count)

        rhs_uniques, rhs_coefs = [], []
        for k,count in rhs_cnt.items():
            if count>0:
                fm, ch, ph = k
                ion = Compound(fm)
                ion.charge = ch
                ion.phase  = ph
                rhs_uniques.append(ion)
                rhs_coefs.append(count)

        # 5) format
        return self._fmt_equation(lhs_uniques, rhs_uniques,
                                  lhs_coefs,    rhs_coefs)

    # internal helpers
    def _species(self, side, r_coeffs, p_coeffs):
        """Expand each Compound by its coefficient."""
        return ([c for c, n in zip(self.reactants, r_coeffs) for _ in range(n)],
                [c for c, n in zip(self.products,  p_coeffs) for _ in range(n)])

    def _fmt_equation(self, lhs: List[Compound], rhs: List[Compound],
                      lcoeffs: List[int], rcoeffs: List[int]) -> str:
        """Generic formatter once you have expanded species & unit coefficients."""
        def part(sp, coeffs):
            out = []
            for c, n in zip(sp, coeffs):
                prefix = str(n) if n!=1 else ""
                fm = c.input_formula or c.formula()
                if getattr(c, 'phase', None):
                    fm += f"({c.phase})"
                # superscript the charge:
                if getattr(c, 'charge', None):
                    sign = '+' if c.charge>0 else '-'
                    mag  = abs(c.charge)
                    fm   = f"{fm}^{'%d'%mag if mag>1 else ''}{sign}"
                out.append(prefix+fm)
            return " + ".join(out)

        return f"{part(lhs, lcoeffs)} -> {part(rhs, rcoeffs)}"

    def theoretical_yield(self, given: Dict[Compound, float],
                      given_in_mass: bool = False
                     ) -> Tuple[Compound, Dict[Compound, Tuple[float, float]]]:
        """
        Returns:
        - limiting_reagent (Compound)
        - yields: dict mapping each product to (mol, grams) theoretical yield
        """

        r_coeffs, p_coeffs = self.balance()
        avail_mol = {}
        for cmpd, qty in given.items():
            if given_in_mass:
                cmpd.set_mass(qty)
                avail_mol[cmpd] = cmpd.amount_mol
            else:
                cmpd.set_moles(qty)
                avail_mol[cmpd] = qty

        rxn_units = {}
        for cmpd, coeff in zip(self.reactants, r_coeffs):
            if coeff == 0:
                continue
            available = avail_mol.get(cmpd, float('inf'))
            rxn_units[cmpd] = available / coeff

        limiting = min(rxn_units, key=rxn_units.get)
        max_rxn = rxn_units[limiting]

        yields = {}
        for prod, coeff in zip(self.products, p_coeffs):
            prod_mol = max_rxn * coeff
            prod_mass = prod_mol * prod.molar_mass
            yields[prod] = (prod_mol, prod_mass)

        return limiting, yields


    def __str__(self) -> str:
        """
        Return a multi-line string with:
         * molecular equation  [type]
         * ionic equation      [type]
         * net ionic equation  [type]
        """
        mol = self.as_molecular()
        ionic = self.as_ionic()
        net = self.as_net_ionic()
        return (
            f"Reaction type: {self.reaction_type()}\n"
            f"Molecular: {mol}\n"
            f"Ionic:     {ionic}\n"
            f"Net ionic: {net}"
        )
