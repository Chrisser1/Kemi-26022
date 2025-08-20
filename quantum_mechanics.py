import math
import re

from compound import PT

# --- Physical Constants ---
SPEED_OF_LIGHT_MS = 2.998e8  # Speed of light in m/s
PLANCK_CONSTANT_JS = 6.626e-34 # Planck's constant in J·s
PI = math.pi


# --- Core Functions ---
def calculate_wavelength(frequency: float) -> float:
    """
    Calculates wavelength in meters from frequency in Hz.
    λ = c / ν
    """
    if frequency == 0:
        raise ValueError("Frequency cannot be zero.")
    return SPEED_OF_LIGHT_MS / frequency

# --- Unit Conversions ---
def m_to_nm(meters: float) -> float:
    """Converts meters to nanometers."""
    return meters * 1e9

def nm_to_m(nanometers: float) -> float:
    """Converts nanometers to meters."""
    return nanometers * 1e-9

    return nanometers * 1e-9

# --- Wave-Particle Duality Functions ---
def de_broglie_wavelength(mass_kg: float, velocity_ms: float) -> float:
    """
    Calculates the de Broglie wavelength (in meters) of a particle.
    λ = h / (m * v)
    """
    if mass_kg == 0 or velocity_ms == 0:
        raise ValueError("Mass and velocity cannot be zero.")
    return PLANCK_CONSTANT_JS / (mass_kg * velocity_ms)

def de_broglie_velocity(mass_kg: float, wavelength_m: float) -> float:
    """
    Calculates the velocity of a particle from its de Broglie wavelength.
    v = h / (m * λ)
    """
    if mass_kg == 0 or wavelength_m == 0:
        raise ValueError("Mass and wavelength cannot be zero.")
    return PLANCK_CONSTANT_JS / (mass_kg * wavelength_m)

# --- Heisenberg Uncertainty Principle ---
def heisenberg_position_uncertainty(mass_kg: float, delta_v_ms: float) -> float:
    """
    Calculates the minimum uncertainty in position (Δx).
    Δx = h / (4 * π * m * Δv)
    """
    if mass_kg <= 0 or delta_v_ms <= 0:
        raise ValueError("Mass and velocity uncertainty must be positive.")
    
    delta_p = mass_kg * delta_v_ms  # Uncertainty in momentum
    return PLANCK_CONSTANT_JS / (4 * PI * delta_p)

def analyze_simple_mo_filling(element_symbol: str, num_atoms: int) -> dict:
    """
    Analyzes the electron filling of molecular orbitals formed from s-orbitals.
    
    Returns a dictionary with the number of electrons in bonding and antibonding orbitals.
    """
    if element_symbol not in PT:
        raise ValueError(f"Element '{element_symbol}' not found.")
        
    # Assumes s-orbital, so we look at group 1 or 2 for valence electrons
    valence_electrons_per_atom = PT[element_symbol].group
    total_electrons = num_atoms * valence_electrons_per_atom
    
    num_bonding_mos = num_atoms // 2
    bonding_capacity = num_bonding_mos * 2
    
    electrons_in_bonding = min(total_electrons, bonding_capacity)
    electrons_in_antibonding = total_electrons - electrons_in_bonding
    
    return {
        "total_valence_electrons": total_electrons,
        "bonding_mo_electrons": electrons_in_bonding,
        "antibonding_mo_electrons": electrons_in_antibonding
    }

# Simplified orbital filling order for elements up to Krypton (Z=36)
ORBITAL_ORDER = [
    ('1s', 2), ('2s', 2), ('2p', 6), ('3s', 2), ('3p', 6),
    ('4s', 2), ('3d', 10), ('4p', 6)
]

def get_ion_configuration(element_symbol: str, charge: int) -> str:
    """
    Determines the electron configuration of an ion for elements up to Z=36.
    """
    if element_symbol not in PT:
        raise ValueError(f"Element '{element_symbol}' not found in periodic table data.")

    atomic_number = PT[element_symbol].atomic_number
    if atomic_number > 36:
        raise NotImplementedError("This function only supports elements up to Krypton (Z=36).")

    # 1. Determine the configuration of the neutral atom
    neutral_electrons = atomic_number
    neutral_config = {}
    electrons_to_place = neutral_electrons
    
    for orbital, capacity in ORBITAL_ORDER:
        if electrons_to_place > 0:
            n_electrons = min(electrons_to_place, capacity)
            neutral_config[orbital] = n_electrons
            electrons_to_place -= n_electrons
        else:
            break
            
    # 2. Remove electrons to form the cation
    electrons_to_remove = charge
    ion_config = neutral_config.copy()
    
    # Sort orbitals by principal number (n) descending, then by l descending (p>s)
    # This ensures 4s is removed before 3d
    sorted_orbitals = sorted(
        ion_config.keys(), 
        key=lambda o: (int(o[0]), o[1]), 
        reverse=True
    )

    for orbital in sorted_orbitals:
        if electrons_to_remove > 0:
            electrons_in_orbital = ion_config[orbital]
            removed = min(electrons_to_remove, electrons_in_orbital)
            ion_config[orbital] -= removed
            electrons_to_remove -= removed
            if ion_config[orbital] == 0:
                del ion_config[orbital]

    # 3. Format the output string with noble gas shorthand
    noble_gases = {18: "Ar", 10: "Ne", 2: "He"}
    ion_electron_count = atomic_number - charge
    
    base_noble_gas = ""
    noble_gas_electrons = 0
    # Find the correct noble gas to use as a base
    for num, sym in sorted(noble_gases.items(), reverse=True):
        if atomic_number > num:
            base_noble_gas = f"[{sym}]"
            noble_gas_electrons = num
            break

    # Rebuild the final string part
    final_config_str = ""
    # Sort for conventional output order (e.g., 3d before 4s)
    final_sorted_orbitals = sorted(ion_config.keys(), key=lambda o: (int(o[0]), o[1]))
    
    electrons_accounted_for = noble_gas_electrons
    for orbital, capacity in ORBITAL_ORDER:
         # Only include orbitals that come after the noble gas
        if electrons_accounted_for < atomic_number:
            if orbital in ion_config:
                 final_config_str += f"{orbital}{ion_config[orbital]} "
        # Update accounted electrons based on the neutral atom's full shell
        electrons_accounted_for += capacity


    return (base_noble_gas + " " + final_config_str).strip()

def count_unpaired_electrons(element_symbol: str) -> int:
    """
    Counts the number of unpaired electrons for a neutral atom in its ground state.
    """
    config_str = get_ion_configuration(element_symbol, charge=0)
    
    # Data for subshell properties
    subshell_data = {'s': 1, 'p': 3, 'd': 5, 'f': 7}
    
    # Find the last (valence) subshell
    last_part = config_str.split()[-1]
    match = re.match(r'\d+([spdf])(\d+)', last_part)
    
    if not match:
         # Handle cases like He (1s2) without a noble gas shorthand
        match_no_shorthand = re.match(r'(\d+)([spdf])(\d+)', config_str)
        if not match_no_shorthand or element_symbol in ["He", "Ne", "Ar", "Kr", "Xe", "Rn"]:
            return 0 # Noble gases and fully-filled s-shells have 0 unpaired electrons
        orbital_type = match_no_shorthand.group(2)
        electron_count = int(match_no_shorthand.group(3))
    else:
        orbital_type = match.group(1)
        electron_count = int(match.group(2))

    num_orbitals = subshell_data[orbital_type]
    
    # Apply Hund's rule
    if electron_count <= num_orbitals:
        # All electrons are unpaired
        return electron_count
    else:
        # Some electrons are paired
        return (2 * num_orbitals) - electron_count
    
AVOGADRO_CONSTANT = 6.022e23

def calculate_photoelectron_kinetic_energy(
    wavelength_nm: float,
    work_function: float,
    work_function_unit: str = 'kJ/mol'
) -> float:
    """
    Calculates the kinetic energy of a photoelectron in Joules.
    """
    # 1. Calculate photon energy in Joules
    wavelength_m = nm_to_m(wavelength_nm)
    energy_photon_J = (PLANCK_CONSTANT_JS * SPEED_OF_LIGHT_MS) / wavelength_m

    # 2. Convert work function to Joules per electron
    work_function_J_per_electron = 0
    if work_function_unit.lower() == 'kj/mol':
        work_function_J_per_electron = (work_function * 1000) / AVOGADRO_CONSTANT
    elif work_function_unit.lower() == 'ev':
        work_function_J_per_electron = work_function * 1.602e-19
    elif work_function_unit.lower() == 'j':
         work_function_J_per_electron = work_function
    else:
        raise ValueError(f"Unsupported work function unit: {work_function_unit}")

    # 3. Calculate kinetic energy
    kinetic_energy = energy_photon_J - work_function_J_per_electron
    
    # Kinetic energy cannot be negative
    return max(0, kinetic_energy)

def calculate_wavelength_from_energy(energy_J: float) -> float:
    """
    Calculates wavelength in meters from energy in Joules.
    λ = hc / E
    """
    if energy_J == 0:
        raise ValueError("Energy cannot be zero.")
    return (PLANCK_CONSTANT_JS * SPEED_OF_LIGHT_MS) / energy_J

def classify_radiation_from_wavelength(wavelength_nm: float) -> str:
    """
    Classifies electromagnetic radiation based on its wavelength in nanometers.
    """
    if wavelength_nm < 0.01:
        return "Gamma ray"
    elif wavelength_nm < 10:
        return "X-ray"
    elif wavelength_nm < 400:
        return "Ultraviolet"
    elif wavelength_nm < 700:
        return "Visible light"
    elif wavelength_nm < 1e6: # 1 mm
        return "Infrared"
    elif wavelength_nm < 1e9: # 1 m
        return "Microwave"
    else:
        return "Radio wave"