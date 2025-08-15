import math

# --- Constants ---
AVOGADRO_NUMBER = 6.022e23

# --- Unit Conversion Helpers ---
_DENSITY_CONVERSION = { # to g/cm³
    "g/cm3": 1.0,
    "g/ml": 1.0,
    "kg/m3": 0.001
}
_LENGTH_CONVERSION = { # to cm
    "cm": 1.0,
    "m": 100.0,
    "pm": 1e-10,
    "nm": 1e-7
}

def _to_g_cm3(value: float, unit: str) -> float:
    return value * _DENSITY_CONVERSION.get(unit.lower(), 1.0)

def _from_cm(value_cm: float, to_unit: str) -> float:
    return value_cm / _LENGTH_CONVERSION.get(to_unit.lower(), 1.0)

# --- Data for Crystal Structures ---
CRYSTAL_STRUCTURES = {
    "BCC": {"Z": 2, "radius_formula": lambda a: (math.sqrt(3) * a) / 4},
    "FCC": {"Z": 4, "radius_formula": lambda a: (math.sqrt(2) * a) / 4},
    "SCC": {"Z": 1, "radius_formula": lambda a: a / 2}
}

# --- Main Calculation Function ---
def calculate_atomic_radius(
    density: float,
    molar_mass_g_mol: float,
    crystal_type: str,
    density_unit: str = "g/cm3",
    output_unit: str = "pm"
) -> float:
    """
    Calculates atomic radius from density and crystal structure with unit handling.
    """
    crystal_type = crystal_type.upper()
    if crystal_type not in CRYSTAL_STRUCTURES:
        raise ValueError(f"Crystal type '{crystal_type}' not supported.")
        
    structure_data = CRYSTAL_STRUCTURES[crystal_type]
    Z = structure_data["Z"]
    
    # 1. Convert density to g/cm³
    density_g_cm3 = _to_g_cm3(density, density_unit)
    
    # 2. Calculate the edge length 'a' in cm
    a_cubed = (Z * molar_mass_g_mol) / (density_g_cm3 * AVOGADRO_NUMBER)
    a_cm = a_cubed**(1/3)
    
    # 3. Calculate the radius 'r' in cm
    r_cm = structure_data["radius_formula"](a_cm)
    
    # 4. Convert radius to the desired output unit and return
    return _from_cm(r_cm, output_unit)