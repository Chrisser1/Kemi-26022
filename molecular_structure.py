# molecular_structure.py

import re
from typing import Dict, Tuple

from compound import PT, Compound

def estimate_reaction_enthalpy(
    bonds_broken: Dict[str, int],
    bonds_formed: Dict[str, int],
    bond_energies: Dict[str, float]
) -> float:
    """
    Estimates the enthalpy of reaction from dictionaries of bonds and energies.
    
    Args:
        bonds_broken: A dictionary of {bond_type: count} for reactants.
        bonds_formed: A dictionary of {bond_type: count} for products.
        bond_energies: A dictionary of {bond_type: energy_in_kJ/mol}.
        
    Returns:
        The estimated ΔH_rxn in kJ/mol.
    """
    try:
        energy_in = sum(bond_energies[bond] * count for bond, count in bonds_broken.items())
        energy_out = sum(bond_energies[bond] * count for bond, count in bonds_formed.items())
        
        return energy_in - energy_out
    except KeyError as e:
        raise ValueError(f"Energy for bond '{e}' not found in the provided bond_energies dictionary.")

# --- Data for VSEPR Theory ---
VALENCE_ELECTRONS = {el.symbol: (el.group if el.group <= 2 else el.group - 10) for el in PT.values() if el.group is not None}
# Halogens contribute 1 electron for bonding
BONDING_ELECTRONS = {'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'H': 1, 'O': 2} 

VSEPR_GEOMETRIES = {
    2: {0: ("Linear", "Linear")},
    3: {0: ("Trigonal Planar", "Trigonal Planar"), 1: ("Trigonal Planar", "Bent")},
    4: {0: ("Tetrahedral", "Tetrahedral"), 1: ("Tetrahedral", "Trigonal Pyramidal"), 2: ("Tetrahedral", "Bent")},
    5: {0: ("Trigonal Bipyramidal", "Trigonal Bipyramidal"), 1: ("Trigonal Bipyramidal", "Seesaw"), 2: ("Trigonal Bipyramidal", "T-shaped"), 3: ("Trigonal Bipyramidal", "Linear")},
    6: {0: ("Octahedral", "Octahedral"), 1: ("Octahedral", "Square Pyramidal"), 2: ("Octahedral", "Square Planar")}
}

# --- VSEPR Function ---

# NEW: A dictionary for how many electrons surrounding atoms typically accept in bonds.
ELECTRONS_ACCEPTED_IN_BONDS = {'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'H': 1, 'O': 2}

# In molecular_structure.py, make sure you have this version of the function.
def predict_molecular_geometry(compound_formula: str) -> Tuple[str, str, int, int]:
    """
    Predicts the molecular geometry of a simple compound using VSEPR theory.
    (Updated to handle multiple bonds and ions)
    """
    c = Compound(compound_formula)
    
    if not c.composition:
        raise ValueError("Cannot parse compound formula.")

    # --- Improved Central Atom Logic ---
    if len(c.composition) > 1:
        non_h_atoms = {el: count for el, count in c.composition.items() if el != 'H'}
        if 1 in non_h_atoms.values():
            central_atom = [el for el, count in non_h_atoms.items() if count == 1][0]
        else:
            central_atom = min(non_h_atoms, key=lambda el: Compound._EN.get(el, 99))
    else:
        central_atom = list(c.composition.keys())[0]

    surrounding_atoms = {el: count for el, count in c.composition.items() if el != central_atom}
    
    # --- VSEPR Logic based on Electron Domains ---
    total_valence_e = VALENCE_ELECTRONS.get(central_atom, 0) + \
                      sum(VALENCE_ELECTRONS.get(el, 0) * count for el, count in surrounding_atoms.items())
    # Account for the overall charge of the ion
    if c.charge is not None:
        total_valence_e -= c.charge

    bonding_domains = sum(surrounding_atoms.values())
    electrons_in_bonds = bonding_domains * 2
    remaining_electrons = total_valence_e - electrons_in_bonds

    electrons_for_central_atom = remaining_electrons
    for atom, count in surrounding_atoms.items():
        if atom != 'H':
            electrons_needed = count * 6 
            placed = min(electrons_for_central_atom, electrons_needed)
            electrons_for_central_atom -= placed

    lone_pairs_on_central = electrons_for_central_atom // 2
    total_domains = bonding_domains + lone_pairs_on_central
    
    try:
        electron_geom, molecular_geom = VSEPR_GEOMETRIES[total_domains][lone_pairs_on_central]
        return electron_geom, molecular_geom, bonding_domains, lone_pairs_on_central
    except KeyError:
        raise NotImplementedError(f"Geometry for {total_domains} domains with {lone_pairs_on_central} lone pairs is not defined.")

def _is_symmetrical(molecular_geom: str, surrounding_atoms: Dict[str, int]) -> bool:
    """
    A helper function to determine if a geometry is symmetrical,
    considering if all surrounding atoms are the same.
    """
    symmetrical_geometries = [
        "Linear", "Trigonal Planar", "Tetrahedral",
        "Trigonal Bipyramidal", "Octahedral", "Square Planar"
    ]
    
    # If the shape is inherently asymmetrical (like 'Bent' or 'Seesaw'), it's not symmetrical.
    if molecular_geom not in symmetrical_geometries:
        return False
        
    # If the shape is symmetrical, we then check if all surrounding atoms are identical.
    # If there's more than one type of surrounding atom (e.g., CH₂Cl₂), it's not symmetrical.
    if len(surrounding_atoms) > 1:
        return False
        
    return True

def analyze_molecular_properties(formula: str) -> str:
    """
    Provides a full VSEPR and polarity analysis for a molecule, including
    quantitative bond polarity (ΔEN).
    """
    try:
        c = Compound(formula)
        atoms = list(c.composition.keys())
        
        # --- Handle Diatomic Molecules ---
        if sum(c.composition.values()) == 2:
            max_en_diff = 0.0
            bond_summary = "Nonpolar Covalent"
            polarity_summary = "Nonpolar"

            if len(atoms) == 2: # Heteronuclear (e.g., HBr)
                en_diff = abs(Compound._EN.get(atoms[0], 0) - Compound._EN.get(atoms[1], 0))
                max_en_diff = en_diff
                if en_diff > 0.4:
                    bond_summary = "Polar Covalent"
                    polarity_summary = "Polar"
            
            summary = (
                f"--- Analysis for {formula} ---\n"
                f"  - Note: This is a diatomic molecule.\n"
                f"  - Molecular Geometry: Linear\n"
                f"  - Bond Polarity: {bond_summary}\n"
                f"  - Polarity: {polarity_summary}\n"
                f"  - Max Bond ΔEN: {max_en_diff:.2f}"
            )
            return summary

        # --- Logic for Polyatomic Molecules ---
        electron_geom, molecular_geom, bp, lp = predict_molecular_geometry(formula)
        central_atom = min(c.composition, key=c.composition.get)
        surrounding_atoms = {el: count for el, count in c.composition.items() if el != central_atom}
        
        # --- Analyze Bond Polarity ---
        max_en_diff = 0.0
        for atom in surrounding_atoms.keys():
            en_diff = abs(Compound._EN.get(central_atom, 0) - Compound._EN.get(atom, 0))
            max_en_diff = max(max_en_diff, en_diff)
        
        bond_summary = "Polar Covalent" if max_en_diff > 0.4 else "Nonpolar Covalent"

        # --- Determine Overall Polarity ---
        symmetrical = _is_symmetrical(molecular_geom, surrounding_atoms)
        is_nonpolar = (bond_summary == "Nonpolar Covalent") or symmetrical
        polarity_summary = "Nonpolar" if is_nonpolar else "Polar"
        
        summary = (
            f"--- Analysis for {formula} ---\n"
            f"  - Central Atom: {central_atom}\n"
            f"  - Bonding Pairs: {bp} | Lone Pairs: {lp}\n"
            f"  - Electron Geometry: {electron_geom}\n"
            f"  - Molecular Geometry: {molecular_geom}\n"
            f"  - Symmetry: {'Symmetrical' if symmetrical else 'Asymmetrical'}\n"
            f"  - Bond Polarity: {bond_summary}\n"
            f"  - Polarity: {polarity_summary}\n"
            f"  - Max Bond ΔEN: {max_en_diff:.2f}" # ADDED THIS LINE
        )
        return summary
        
    except (ValueError, NotImplementedError, KeyError) as e:
        return f"Could not analyze {formula}: {e}"
    
def find_most_polar_molecule(formulas: list[str]) -> str:
    """
    Analyzes a list of molecules and identifies the one with the largest net dipole moment.
    (Upgraded to quantitatively compare polar molecules)
    """
    results = {}
    most_polar_molecule = None
    max_en_diff = -1.0

    for formula in formulas:
        summary = analyze_molecular_properties(formula)
        results[formula] = summary
        
        # Use regex to find the ΔEN value in the summary
        match = re.search(r"Max Bond ΔEN: (\d+\.\d+)", summary)
        if match:
            en_diff = float(match.group(1))
            if en_diff > max_en_diff:
                max_en_diff = en_diff
                most_polar_molecule = formula
    
    # --- Compile the final report ---
    report = "--- Polarity Analysis Report ---\n\n"
    for formula in formulas:
        report += f"{results[formula]}\n\n"
    
    report += "--- Conclusion ---\n"
    if most_polar_molecule and max_en_diff > 0.4:
         report += (f"Based on the largest electronegativity difference (ΔEN = {max_en_diff:.2f}),\n"
                   f"the molecule with the largest net dipole moment is {most_polar_molecule}.")
    else:
        report += "All molecules are nonpolar or have very low polarity."
        
    return report

def rank_by_ionic_character(formulas: list[str]) -> str:
    """
    Ranks a list of binary/diatomic compounds in order of increasing ionic character
    based on their electronegativity difference.
    """
    bond_characters = []
    for f in formulas:
        try:
            c = Compound(f)
            composition = c.composition
            delta_en = 0.0

            # Case 1: Homonuclear diatomic molecule (e.g., F2)
            if len(composition) == 1 and sum(composition.values()) == 2:
                delta_en = 0.0
            # Case 2: Heteronuclear binary compound (e.g., HCl)
            elif len(composition) == 2:
                elem1, elem2 = composition.keys()
                en1 = Compound._EN.get(elem1, 0)
                en2 = Compound._EN.get(elem2, 0)
                delta_en = abs(en1 - en2)
            else:
                print(f"Warning: '{f}' is not a diatomic/binary compound. Skipping.")
                continue
            
            bond_characters.append((f, delta_en))
        except Exception as e:
            print(f"Could not process {f}: {e}")

    # Sort the list based on the calculated ΔEN
    bond_characters.sort(key=lambda x: x[1])
    
    # Format the ranked list into a string
    ranked_list = [item[0] for item in bond_characters]
    return " < ".join(ranked_list)

# A pre-computed dictionary for the number of structural isomers of alkanes (CnH2n+2)
ALKANE_ISOMER_COUNT = {
    1: 1,  # Methane
    2: 1,  # Ethane
    3: 1,  # Propane
    4: 2,  # Butane
    5: 3,  # Pentane
    6: 5,  # Hexane
    7: 9,  # Heptane
    8: 18, # Octane
    9: 35, # Nonane
    10: 75, # Decane
}

def count_alkane_isomers(formula: str) -> int:
    """
    Counts the number of possible structural isomers for a saturated alkane
    using a lookup table for formulas up to C10H22.
    """
    # Use regex to extract the number of carbon atoms
    match = re.match(r'C(\d+)H(\d+)', formula)
    if not match:
        raise ValueError("Formula must be in the format CnH2n+2 (e.g., 'C6H14').")
    
    n_carbons = int(match.group(1))
    n_hydrogens = int(match.group(2))

    # Verify it's a saturated alkane
    if n_hydrogens != 2 * n_carbons + 2:
        raise ValueError("The formula does not match a saturated alkane (CnH2n+2).")

    # Look up the number of isomers from the pre-computed dictionary
    count = ALKANE_ISOMER_COUNT.get(n_carbons)
    if count is None:
        raise NotImplementedError(f"Isomer count for alkanes with {n_carbons} carbons is not available in the lookup table.")
        
    return count