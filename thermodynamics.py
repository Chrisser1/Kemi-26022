def enthalpy_from_surroundings_entropy(delta_s_surr: float, T_kelvin: float) -> float:
    """
    Calculates the enthalpy change of a reaction from the entropy change
    of the surroundings.

    Args:
        delta_s_surr: The change in entropy for the surroundings in J/(mol·K).
        T_kelvin: The temperature in Kelvin.

    Returns:
        The enthalpy change for the reaction in kJ/mol.
    """
    # ΔH_rxn = -T * ΔS_surr
    delta_h_joules = -T_kelvin * delta_s_surr
    
    # Convert from J/mol to kJ/mol
    delta_h_kj = delta_h_joules / 1000.0
    
    return delta_h_kj

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
    
    # Checks if the unit string contains 'j/' (for Joules)
    delta_S_kJ = delta_S / 1000.0 if 'j/' in S_unit.lower() else delta_S

    # 2. Calculate ΔG
    delta_G = delta_H - (T_k * delta_S_kJ)

    # 3. Determine spontaneity
    if delta_G < 0:
        spontaneity = "ΔG is negative, so the reaction is spontaneous."
    elif delta_G > 0:
        spontaneity = "ΔG is positive, so the reaction is non-spontaneous."
    else:
        spontaneity = "ΔG is zero, so the reaction is at equilibrium."

    return delta_G, spontaneity

def entropy_of_universe(delta_s_sys: float, delta_h_rxn_kj: float, T_kelvin: float) -> dict[str, float]:
    """
    Calculates the entropy change for the surroundings and the universe.

    Args:
        delta_s_sys: The entropy change for the system in J/(mol·K).
        delta_h_rxn_kj: The enthalpy change for the reaction in kJ/mol.
        T_kelvin: The temperature in Kelvin.

    Returns:
        A dictionary containing ΔS_surr and ΔS_univ in J/(mol·K).
    """
    # 1. Convert enthalpy from kJ to J
    delta_h_rxn_j = delta_h_rxn_kj * 1000.0

    # 2. Calculate entropy of the surroundings
    delta_s_surr = -delta_h_rxn_j / T_kelvin

    # 3. Calculate entropy of the universe
    delta_s_univ = delta_s_sys + delta_s_surr

    return {
        "delta_s_surr": delta_s_surr,
        "delta_s_univ": delta_s_univ
    }

def entropy_of_surroundings(delta_h_rxn_kj: float, T_kelvin: float) -> float:
    """
    Calculates the entropy change for the surroundings from the reaction's
    enthalpy change.

    Args:
        delta_h_rxn_kj: The enthalpy change for the reaction in kJ/mol.
        T_kelvin: The temperature in Kelvin.

    Returns:
        The entropy change for the surroundings in J/(mol·K).
    """
    # 1. Convert enthalpy from kJ to J
    delta_h_rxn_j = delta_h_rxn_kj * 1000.0

    # 2. Calculate entropy of the surroundings
    delta_s_surr = -delta_h_rxn_j / T_kelvin
    
    return delta_s_surr