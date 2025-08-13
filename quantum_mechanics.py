import math

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