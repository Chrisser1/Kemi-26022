# molecular_structure.py

from typing import Dict

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