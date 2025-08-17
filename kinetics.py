import math


def calculate_relative_rate(
    rate_known: float,
    coeff_known: int,
    coeff_target: int
) -> float:
    """
    Calculates the rate of change of a target species based on the rate of a known species.
    Assumes the sign convention (negative for reactants, positive for products) is handled
    by the input and the desired output.
    
    The relationship is: rate_target / coeff_target = rate_known / coeff_known
    """
    if coeff_known == 0:
        raise ValueError("Coefficient of the known species cannot be zero.")
    
    # We use the absolute rate for calculation, the signs are a matter of convention.
    # Rate of reaction = abs(rate_known) / coeff_known
    # Rate of target = Rate of reaction * coeff_target
    return (abs(rate_known) / coeff_known) * coeff_target

def determine_reaction_order(
    rate1: float,
    conc1: float,
    rate2: float,
    conc2: float
) -> float:
    """
    Determines the reaction order for a single reactant using two experiments.
    Solves for n in: (rate2/rate1) = (conc2/conc1)^n
    """
    if any(x == 0 for x in [rate1, conc1, rate2, conc2]):
        raise ValueError("Rates and concentrations cannot be zero.")
        
    rate_ratio = rate2 / rate1
    conc_ratio = conc2 / conc1
    
    # n = log(rate_ratio) / log(conc_ratio)
    order = math.log(rate_ratio) / math.log(conc_ratio)
    
    # Round to the nearest integer or common half-integer
    return round(order * 2) / 2

# --- Constants ---
LN_2 = math.log(2)

# --- Core Functions ---
def calculate_decay_constant(half_life: float) -> float:
    """Calculates the decay constant (λ) from the half-life."""
    if half_life <= 0:
        raise ValueError("Half-life must be a positive number.")
    return LN_2 / half_life

def calculate_initial_amount(
    amount_remaining: float,
    time_elapsed: float,
    half_life: float
) -> float:
    """
    Calculates the initial amount of a substance before radioactive decay.
    """
    # First, get the decay constant
    decay_constant = calculate_decay_constant(half_life)
    
    # Then, calculate the initial amount using N₀ = N(t) * e^(λt)
    initial_amount = amount_remaining * math.exp(decay_constant * time_elapsed)
    return initial_amount

# Add to your kinetics.py file

def calculate_average_reaction_rate(
    conc_initial: float,
    conc_final: float,
    t_initial: float,
    t_final: float,
    coefficient: int,
    is_reactant: bool = True
) -> float:
    """
    Calculates the average rate of a reaction from the change in
    concentration of one species.
    """
    delta_conc = conc_final - conc_initial
    delta_t = t_final - t_initial
    
    if delta_t == 0:
        raise ValueError("Time interval cannot be zero.")
    if coefficient == 0:
        raise ValueError("Stoichiometric coefficient cannot be zero.")
        
    rate_of_species = delta_conc / delta_t
    
    # Reaction rate is positive by convention
    if is_reactant:
        return -rate_of_species / coefficient
    else:
        return rate_of_species / coefficient

def solve_second_order_time(
    k: float,
    conc_initial: float,
    conc_final: float
) -> float:
    """
    Calculates the time elapsed for a second-order reaction
    to reach a final concentration.
    t = (1/k) * (1/[A]t - 1/[A]₀)
    """
    if k <= 0:
        raise ValueError("Rate constant k must be positive.")
    if conc_initial <= 0 or conc_final <= 0:
        raise ValueError("Concentrations must be positive.")
    
    time = (1 / k) * ((1 / conc_final) - (1 / conc_initial))
    return time

def calculate_first_order_half_life(k: float) -> float:
    """
    Calculates the half-life of a first-order reaction from the rate constant.
    t_1/2 = ln(2) / k
    """
    if k <= 0:
        raise ValueError("Rate constant k must be positive.")
    return math.log(2) / k

def calculate_half_life(
    initial_amount: float,
    final_amount: float,
    time_elapsed: float
) -> float:
    """
    Calculates the half-life from initial and final amounts over a time interval.
    """
    if final_amount > initial_amount:
        raise ValueError("Final amount cannot be greater than initial amount.")

    # 1. Calculate decay constant λ from N(t) = N₀ * e^(-λt)
    decay_constant = -(1 / time_elapsed) * math.log(final_amount / initial_amount)
    
    # 2. Calculate half-life from λ = ln(2) / t_½
    half_life = math.log(2) / decay_constant
    return half_life

def arrhenius_solve_for_temperature(
    k1: float,
    T1_C: float,
    k2: float,
    Ea_kJ_mol: float
) -> float:
    """
    Calculates the temperature (T2) in Kelvin at which a reaction
    has a rate constant k2, using the two-point Arrhenius equation.
    """
    R = 8.3145  # J/(mol*K)
    
    # Convert units
    T1_K = T1_C + 273.15
    Ea_J_mol = Ea_kJ_mol * 1000
    
    # Rearranged Arrhenius equation
    inv_T2 = (1 / T1_K) - (R / Ea_J_mol) * math.log(k2 / k1)
    
    T2_K = 1 / inv_T2
    return T2_K

def calculate_second_order_half_life(k: float, conc_initial: float) -> float:
    """
    Calculates the half-life of a second-order reaction.
    t_1/2 = 1 / (k * [A]₀)
    """
    if k <= 0 or conc_initial <= 0:
        raise ValueError("Rate constant and initial concentration must be positive.")
    return 1 / (k * conc_initial)

def calculate_second_order_rate_constant(half_life: float, conc_initial: float) -> float:
    """
    Calculates the rate constant (k) of a second-order reaction from its half-life.
    k = 1 / (t_1/2 * [A]₀)
    """
    if half_life <= 0 or conc_initial <= 0:
        raise ValueError("Half-life and initial concentration must be positive.")
    return 1 / (half_life * conc_initial)

def solve_second_order_final_conc(
    k: float,
    conc_initial: float,
    time: float
) -> float:
    """
    Calculates the final concentration for a second-order reaction.
    1/[A]t = kt + 1/[A]₀
    """
    if k <= 0 or conc_initial <= 0 or time < 0:
        raise ValueError("Rate constant, initial concentration, and time must be positive.")
    
    inv_conc_final = (k * time) + (1 / conc_initial)
    return 1 / inv_conc_final

def solve_first_order_initial_conc(
    k: float,
    conc_final: float,
    time: float
) -> float:
    """
    Calculates the initial concentration [A]₀ for a first-order reaction.
    [A]₀ = e^(ln[A]t + kt)
    """
    if k <= 0 or conc_final <= 0 or time < 0:
        raise ValueError("Rate constant, final concentration, and time must be positive.")
    
    exponent = math.log(conc_final) + (k * time)
    return math.exp(exponent)

def calculate_first_order_rate_constant(half_life: float) -> float:
    """
    Calculates the rate constant (k) of a first-order reaction from its half-life.
    k = ln(2) / t_1/2
    """
    if half_life <= 0:
        raise ValueError("Half-life must be a positive number.")
    return math.log(2) / half_life

def calculate_rate_of_change(
    conc_initial_known: float,
    conc_final_known: float,
    t_initial: float,
    t_final: float,
    coeff_known: int,
    coeff_target: int,
    known_is_reactant: bool
) -> float:
    """
    Calculates the rate of change of a target species given data for a known species.
    """
    # 1. Calculate the rate of change of the known species
    delta_conc = conc_final_known - conc_initial_known
    delta_t = t_final - t_initial
    if delta_t == 0:
        raise ValueError("Time interval cannot be zero.")
    rate_of_change_known = delta_conc / delta_t
    
    # 2. Use stoichiometry to find the rate of the target species
    # General rate = (-1/coeff_known) * rate_change_known
    # Rate change target = coeff_target * General rate
    sign = -1 if known_is_reactant else 1
    
    rate_of_change_target = (rate_of_change_known / (sign * coeff_known)) * coeff_target
    
    return rate_of_change_target

