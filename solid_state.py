import math

from compound import PT

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

def _to_cm(value: float, from_unit: str) -> float:
    """Helper to convert various length units to centimeters."""
    return value * _LENGTH_CONVERSION.get(from_unit.lower(), 1.0)

def calculate_density_from_edge_length(
    molar_mass_g_mol: float,
    edge_length: float,
    crystal_type: str,
    edge_length_unit: str = "pm"
) -> float:
    """
    Calculates the density of a metal in g/cm³ from its unit cell edge length.
    """
    crystal_type = crystal_type.upper()
    if crystal_type not in CRYSTAL_STRUCTURES:
        raise ValueError(f"Crystal type '{crystal_type}' not supported.")
        
    Z = CRYSTAL_STRUCTURES[crystal_type]["Z"]
    
    # 1. Convert edge length to cm and calculate volume
    a_cm = _to_cm(edge_length, edge_length_unit)
    volume_cm3 = a_cm**3
    
    # 2. Calculate density using the standard formula
    density = (Z * molar_mass_g_mol) / (volume_cm3 * AVOGADRO_NUMBER)
    
    return density


def calculate_atoms_in_spherical_particle(
    element_symbol: str,
    density_g_cm3: float,
    radius: float,
    radius_unit: str = "nm"
) -> int:
    """
    Calculates the number of atoms in a spherical particle.
    """
    # 1. Get molar mass and Avogadro's number
    molar_mass = PT[element_symbol].atomic_mass
    
    # 2. Convert radius to cm and calculate volume
    radius_cm = _to_cm(radius, radius_unit)
    volume_cm3 = (4/3) * math.pi * (radius_cm**3)
    
    # 3. Calculate the mass of the particle
    mass_g = volume_cm3 * density_g_cm3
    
    # 4. Calculate moles and then the number of atoms
    moles = mass_g / molar_mass
    num_atoms = moles * AVOGADRO_NUMBER
    
    return int(round(num_atoms))
