# universal R in L·bar·K⁻¹·mol⁻¹
R_BAR = 0.083145

# how many bar is 1 of these units?
_PRESSURE_TO_BAR = {
    "bar":   1.0,
    "atm":   1.01325,
    "pa":    1e-5,        # 1 Pa = 1e-5 bar
    "kpa":   1e-2,        # 1 kPa = 0.01 bar
    "torr":  0.00133322,  # 1 Torr = 133.322 Pa = 0.00133322 bar
    "mmhg":  0.00133322,
    "psi":   0.0689476,
}

# how many litres is 1 of these units?
_VOLUME_TO_LITRE = {
    "l":    1.0,
    "ml":   1e-3,
    "m3":   1000.0,
    "cm3":  1e-3,
}

def _to_bar(p: float, unit: str) -> float:
    try:
        return p * _PRESSURE_TO_BAR[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown pressure unit '{unit}'")

def _from_bar(p_bar: float, unit: str) -> float:
    try:
        return p_bar / _PRESSURE_TO_BAR[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown pressure unit '{unit}'")

def _to_litre(V: float, unit: str) -> float:
    try:
        return V * _VOLUME_TO_LITRE[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown volume unit '{unit}'")

def _from_litre(V_l: float, unit: str) -> float:
    try:
        return V_l / _VOLUME_TO_LITRE[unit.lower()]
    except KeyError:
        raise ValueError(f"Unknown volume unit '{unit}'")

def c_to_k(T_c: float) -> float:
    return T_c + 273.15

def ideal_gas_volume(n: float,
                     T: float,
                     p: float,
                     T_unit: str = "C",
                     P_unit: str = "bar",
                     V_unit: str = "L") -> float:
    """
    Compute V from PV = nRT.
      - n in mol
      - T in °C or K
      - p in bar, atm, Pa, …
    Returns volume in L, mL, m3, …
    """
    if T_unit.lower() == "c":
        T_k = c_to_k(T)
    elif T_unit.lower() == "k":
        T_k = T
    else:
        raise ValueError(f"Unknown temperature unit '{T_unit}'")

    p_bar = _to_bar(p, P_unit)
    V_l   = n * R_BAR * T_k / p_bar
    return _from_litre(V_l, V_unit)

def ideal_gas_pressure(n: float,
                       T: float,
                       V: float,
                       T_unit: str = "C",
                       V_unit: str = "L",
                       p_unit: str = "bar") -> float:
    """
    Compute p from PV = nRT.
      - V in L, mL, m3 …
    Returns pressure in bar, atm, Pa, …
    """
    if T_unit.lower() == "c":
        T_k = c_to_k(T)
    elif T_unit.lower() == "k":
        T_k = T
    else:
        raise ValueError(f"Unknown temperature unit '{T_unit}'")

    V_l = _to_litre(V, V_unit)
    p_bar = n * R_BAR * T_k / V_l
    return _from_bar(p_bar, p_unit)

def ideal_gas_amount(V: float,
                     T: float,
                     P: float,
                     V_unit: str = "L",
                     T_unit: str = "C",
                     p_unit: str = "bar") -> float:
    """
    Return n (mol) from PV = nRT.
    V in L/mL/m³, T in °C or K, P in bar/atm/Pa/…
    """
    # 1) temperature → K
    if T_unit.lower() == "c":
        T_k = T + 273.15
    elif T_unit.lower() == "k":
        T_k = T
    else:
        raise ValueError(f"Unknown T unit '{T_unit}'")
    # 2) pressure → bar
    P_bar = _to_bar(P, p_unit)
    # 3) volume → L
    V_L   = _to_litre(V, V_unit)
    # 4) n = PV / (R·T)
    return P_bar * V_L / (R_BAR * T_k)

def total_pressure_and_mole_fractions(
    partials: dict[str, float],
    p_unit: str = "atm"
) -> tuple[float, dict[str, float]]:
    """
    Given a dict of {species: partial_pressure}, returns
    (P_total, {species: mole_fraction}).
    """
    # 1) convert all partials into bar (or atm, your choice)
    p_bar = {sp: _to_bar(p, p_unit) for sp, p in partials.items()}
    # 2) total pressure
    P_tot = sum(p_bar.values())
    # 3) mole fractions
    x = {sp: p_val / P_tot for sp, p_val in p_bar.items()}
    # if you want P_tot back in the original unit:
    P_tot_in = _from_bar(P_tot, p_unit)
    return P_tot_in, x

def ideal_gas_temperature(n: float,
                          P: float,
                          V: float,
                          p_unit: str = "bar",
                          V_unit: str = "L",
                          T_unit: str = "C") -> float:
    """
    Solve T from PV = n R T.
      - n in mol
      - P in bar, atm, Pa, …
      - V in L, mL, m3, …
    Returns T in °C or K.
    """
    # 1) convert P → bar, V → L
    p_bar = _to_bar(P, p_unit)
    V_l   = _to_litre(V, V_unit)

    # 2) compute T in Kelvin
    T_k = p_bar * V_l / (n * R_BAR)

    # 3) return in requested unit
    if T_unit.lower() == "c":
        return T_k - 273.15
    elif T_unit.lower() == "k":
        return T_k
    else:
        raise ValueError(f"Unknown temperature unit '{T_unit}'")

def ideal_gas_moles(P: float,
                    V: float,
                    T: float,
                    p_unit: str = "bar",
                    V_unit: str = "L",
                    T_unit: str = "C") -> float:
    """
    Solve n from PV = nRT
      - P in bar, atm, Pa, …
      - V in L, mL, m3, …
      - T in °C or K
    Returns n in mol.
    """
    P_bar = _to_bar(P, p_unit)
    V_l   = _to_litre(V, V_unit)
    if T_unit.lower()=="c":
        T_k = c_to_k(T)
    elif T_unit.lower()=="k":
        T_k = T
    else:
        raise ValueError(f"Unknown T unit {T_unit!r}")

    # 2) n = P V / (R T)
    return P_bar * V_l / (R_BAR * T_k)

def convert_pressure(value: float,
                     from_unit: str,
                     to_unit: str) -> float:
    """
    Convert a pressure from one unit to another.
      - value: numeric pressure
      - from_unit: any of 'bar','atm','Pa','kPa','torr','mmHg','psi', etc.
      - to_unit:   same possible units
    Returns the converted pressure.
    """
    # 1) send it to bar
    p_bar = _to_bar(value, from_unit)
    # 2) back out to the target unit
    return _from_bar(p_bar, to_unit)


def convert_volume(value: float,
                   from_unit: str,
                   to_unit: str) -> float:
    """
    Convert a volume from one unit to another.
      - value: numeric volume
      - from_unit: 'L','mL','m3','cm3',…
      - to_unit:   same possible units
    Returns the converted volume.
    """
    # 1) send it to litres
    v_L = _to_litre(value, from_unit)
    # 2) back out to the target unit
    return _from_litre(v_L, to_unit)

def component_moles(total_moles: float, mole_fraction: float) -> float:
    """
    Calculate the moles of a component in a mixture from its mole fraction.
    
    - total_moles: The total number of moles in the gas mixture.
    - mole_fraction: The mole fraction of the component of interest.
    Returns the number of moles of the component.
    """
    return total_moles * mole_fraction

def mole_fraction(component_moles: float, total_moles: float) -> float:
    """
    Calculate the mole fraction of a component in a mixture.
    
    - component_moles: The number of moles of the component of interest.
    - total_moles: The total number of moles in the gas mixture.
    Returns the mole fraction of the component.
    """
    if total_moles == 0:
        return 0
    return component_moles / total_moles

def ideal_gas_molar_mass(density: float,
                         T: float,
                         P: float,
                         density_unit: str = "g/L",
                         T_unit: str = "K",
                         P_unit: str = "bar") -> float:
    """
    Calculate molar mass from density using the Ideal Gas Law.
    
    Returns molar mass in g/mol.
    """
    # Convert inputs to standard units (L, K, bar)
    V_unit_part = density_unit.split('/')[-1]
    V_for_density = _to_litre(1.0, V_unit_part)
    density_g_per_L = density / V_for_density
    
    if T_unit.lower() == "c":
        T_k = c_to_k(T)
    else:
        T_k = T
        
    p_bar = _to_bar(P, P_unit)

    # M = (rho * R * T) / P
    molar_mass = (density_g_per_L * R_BAR * T_k) / p_bar
    return molar_mass

def partial_pressure(component_moles: float, 
                     total_moles: float, 
                     total_pressure: float) -> float:
    """
    Calculates the partial pressure of a component in a gas mixture.
    
    - component_moles: Moles of the gas of interest.
    - total_moles: Total moles of gas in the mixture.
    - total_pressure: The total pressure of the mixture (in any unit).
    Returns the partial pressure in the same unit as the total pressure.
    """
    # This function assumes 'mole_fraction' is also in gases.py
    x = mole_fraction(component_moles, total_moles)
    return x * total_pressure

def ideal_gas_density(molar_mass: float,
                      T: float,
                      P: float,
                      T_unit: str = "K",
                      p_unit: str = "bar") -> float:
    """
    Calculate gas density from molar mass using the Ideal Gas Law.
    
    Returns density in g/L.
    """
    # Convert inputs to standard units (K, bar)
    if T_unit.lower() == "c":
        T_k = c_to_k(T)
    else:
        T_k = T
        
    p_bar = _to_bar(P, p_unit)

    # rho = (P * M) / (R * T)
    density = (p_bar * molar_mass) / (R_BAR * T_k)
    return density