# Add Faraday's constant at the top of the file
FARADAY_CONSTANT = 96485  # J/(V·mol) or C/mol

def find_strongest_agent(
    half_reactions: dict[str, float],
    agent_type: str
) -> tuple[str, float]:
    """
    Finds the strongest oxidizing or reducing agent from a dictionary of species
    and their standard reduction potentials.

    Args:
        half_reactions: A dictionary where keys are the ion formulas (e.g., "Au+")
                        and values are their E° values in Volts.
        agent_type: The type of agent to find, either "oxidizing" or "reducing".

    Returns:
        A tuple containing the formula of the strongest agent and its E° value.
    """
    if not half_reactions:
        raise ValueError("The half_reactions dictionary cannot be empty.")

    if agent_type.lower() == 'oxidizing':
        # Strongest oxidizing agent has the MOST POSITIVE E°
        strongest_agent = max(half_reactions, key=half_reactions.get)
    elif agent_type.lower() == 'reducing':
        # Strongest reducing agent has the MOST NEGATIVE E°
        strongest_agent = min(half_reactions, key=half_reactions.get)
    else:
        raise ValueError("agent_type must be 'oxidizing' or 'reducing'.")

    return strongest_agent, half_reactions[strongest_agent]

def calculate_gibbs_from_potential(E_cell_volts: float, n_electrons: int) -> float:
    """
    Calculates the standard free-energy change (ΔG°) from the standard
    cell potential (E°cell).

    Args:
        E_cell_volts: The standard cell potential in Volts.
        n_electrons: The number of moles of electrons transferred in the reaction.

    Returns:
        The standard free-energy change in kJ.
    """
    # ΔG° = -nFE°
    delta_g_joules = -n_electrons * FARADAY_CONSTANT * E_cell_volts
    
    # Convert from Joules to kilojoules
    delta_g_kj = delta_g_joules / 1000.0
    
    return delta_g_kj

def solve_cell_potential(E_cell: float = None, E_cathode: float = None, E_anode: float = None) -> dict[str, float]:
    """
    Solves for an unknown standard potential (E°cell, E°cathode, or E°anode)
    given the other two values. Provide 'None' for the value to be calculated.

    Args:
        E_cell: The standard cell potential in Volts.
        E_cathode: The standard reduction potential of the cathode in Volts.
        E_anode: The standard reduction potential of the anode in Volts.

    Returns:
        A dictionary containing the name and value of the calculated potential.
    """
    knowns = {'E_cell': E_cell, 'E_cathode': E_cathode, 'E_anode': E_anode}
    unknown = [k for k, v in knowns.items() if v is None]

    if len(unknown) != 1:
        raise ValueError("Exactly one potential (E_cell, E_cathode, or E_anode) must be None.")

    var_to_solve = unknown[0]

    if var_to_solve == 'E_anode':
        result = E_cathode - E_cell
        return {'E_anode': result}
    elif var_to_solve == 'E_cathode':
        result = E_cell + E_anode
        return {'E_cathode': result}
    else: # Solving for E_cell
        result = E_cathode - E_anode
        return {'E_cell': result}
    
    from chempy import balance_stoichiometry

def format_equation(reac_dict, prod_dict):
    """Helper function to format the balanced dictionaries into a string."""
    def format_side(d):
        parts = []
        for species, coeff in sorted(d.items()):
            # Don't show coefficient if it's 1, unless it's the only species
            if coeff == 1 and len(d) > 1 and species not in ('H2O', 'H+', 'OH-'):
                 parts.append(species)
            else:
                 parts.append(f"{coeff} {species}")
        return " + ".join(parts)

    return f"{format_side(reac_dict)} -> {format_side(prod_dict)}"

# Activity series from most to least reactive
ACTIVITY_SERIES = [
    'K', 'Ba', 'Ca', 'Na', 'Mg', 'Al', 'Zn', 'Fe', 'Ni', 'Sn', 'Pb', 'H', 'Cu', 'Ag', 'Au', 'Pt'
]

def will_reaction_occur(solid_metal: str, ion_metal: str) -> bool:
    """
    Uses the activity series to predict if a single-displacement reaction will occur.

    Args:
        solid_metal: The symbol of the solid metal reactant (e.g., "Ni").
        ion_metal: The symbol of the metal in the aqueous ion (e.g., "Sn").

    Returns:
        True if the reaction is spontaneous, False otherwise.
    """
    try:
        # Get the reactivity index (lower index is more reactive)
        solid_reactivity = ACTIVITY_SERIES.index(solid_metal)
        ion_reactivity = ACTIVITY_SERIES.index(ion_metal)
        
        # The reaction occurs if the solid metal is more reactive (has a lower index)
        return solid_reactivity < ion_reactivity
    except ValueError as e:
        raise ValueError(f"Metal not found in the activity series: {e}")