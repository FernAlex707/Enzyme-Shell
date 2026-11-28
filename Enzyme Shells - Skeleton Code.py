# Overall reaction: Ethanol + NAD+ <--> Acetaldehyde + NADH + H+
# S    = Substrate (ethanol)
# P    = Product (acetaldehyde)
# NAD  = Oxidized cofactor
# NADH = Reduced cofactor

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

# SECTION 1: INITIAL CONDITIONS AND MODEL PARAMETERS
# --------------------------------------------------------------

R_gas = 8.314             # Universal gas constant [J/(mol*K)]
N_A = 6.02214076e23       # Avogadro's constant [molecules/mol]
T = 298.15                # Operating temperature [K], assumed
T_ref = 298.15            # Reference temperature [K], provisional
P = 101325                # Operating pressure [Pa]
rho = 997.0               # Water density [kg/m^3]
mu = 8.90e-4              # Water dynamic viscosity [Pa*s]

# Pressure drop is not included in the baseline model because the system is a stirred batch reactor, not a packed bed.

V_reactor = 1.0e-3        # Bulk liquid volume [m^3], assumed
t_batch = 3600.0          # Batch duration [s], assumed
N_cycles = 10             # Number of operating cycles
N_shells = 1.0e15         # Total protein shells, assumed - determines total enzyme inventory catalytic capacity for reactor
N_RPM = 300               # Impeller speed [rpm], assumed - affects external mass-transfer resistance

# Part 1 assumes sufficient cofactor availability and neglects reversibility and product inhibition.
# These assumptions allow the reaction to be represented using single-substrate Michaelis-Menten kinetics.

# Reactor/Reaction Initial Conditions
C_S0 = 10.0               # Initial bulk substrate [mol/m^3]
C_P0 = 0.0                # Initial bulk product [mol/m^3]
C_NAD0 = 10.0             # Initial bulk NAD+ [mol/m^3]
C_NADH0 = 0.0             # Initial bulk NADH [mol/m^3]

C_S_in0 = 0.0             # Initial internal substrate
C_P_in0 = 0.0             # Initial internal product
C_NAD_in0 = 0.0           # Initial internal NAD+
C_NADH_in0 = 0.0          # Initial internal NADH

X_0 = 0.0                 # Initial substrate conversion [-]
activity_0 = 1.0          # Initial relative enzyme activity [-]

#Enzyme Shell Geometry
R_shell = 29.0e-9         # External shell radius [m]
L_shell = 4.0e-9          # Shell thickness [m], assumed
R_inner = R_shell - L_shell  # Internal lumen radius [m]

r_pore = 1.25e-9          # Individual pore radius [m]
N_pores = 60              # Number of pores, assumed
tau = 1.5                 # Shell tortuosity [-], assumed - accounts for the increased diffusion distance caused by non-straight transport pathways.

# Enzymatic properties
N_enzyme = 204            # Enzyme molecules per shell
kcat = 0.46               # Turnover frequency [1/s], provisional
Km = 4.7                  # Michaelis constant [mol/m^3], provisional
r_enzyme = 3.0e-9         # Enzyme hydrodynamic radius [m], assumed

# Each enzyme molecule is assumed to contribute one active catalytic site.

# Molecular diffusion properties
D_S = 1.0e-9              # Substrate diffusivity [m^2/s]
D_P = 1.0e-9              # Product diffusivity [m^2/s]
D_NAD = 4.0e-10           # NAD+ diffusivity [m^2/s]
D_NADH = 4.0e-10          # NADH diffusivity [m^2/s]

r_S = 0.23e-9             # Substrate hydrodynamic radius [m]
r_P = 0.23e-9             # Product hydrodynamic radius [m]
r_NAD = 0.70e-9           # NAD+ hydrodynamic radius [m]
r_NADH = 0.70e-9          # NADH hydrodynamic radius [m]

D_lumen = D_S             # Substrate diffusivity inside the lumen [m^2/s]

# Shell permeability
P_shell_min = 1.0e-9      # Minimum permeability [m/s]
P_shell_max = 1.0e-1      # Maximum permeability [m/s]
P_shell = 5.0e-9          # Independent test permeability [m/s]

# Temperature and enzyme deactivation
Ea = 40.0e3               # Reaction activation energy [J/mol], assumed
Ed = 100.0e3              # Deactivation activation energy [J/mol], assumed
kd = 0.0                  # Baseline deactivation constant [1/s]
kd_max = 9.0e-5           # Upper deactivation investigation bound [1/s]

# Enzyme leaching and pore blockage
k_leach = 0.0             # Fraction of enzyme lost per wash [-]
kb = 0.0                  # First-order blockage constant [1/s]
f_open_0 = 1.0            # Initial fraction of open pores [-]
f_enzyme_0 = 1.0          # Initial fraction of retained enzyme [-]
porosity_loss_max = 0.50  # Maximum porosity loss scenario [-]

# Numerical parameters
N_cells = 100             # Number of radial finite-volume cells
N_time = 1000             # Number of time evaluation points
N_sweep = 200             # Number of parameter sweep points

rtol = 1.0e-8             # Relative solver tolerance
atol = 1.0e-10            # Absolute solver tolerance

t_span = (0.0, t_batch)
t_eval = np.linspace(0.0, t_batch, N_time)

solver_method = "BDF"     # Implicit solver for stiff differential equations

# Verify that the initial parameters are physically meaningful.
assert R_shell > L_shell > 0
assert r_pore > 0
assert N_pores >= 0
assert tau >= 1
assert Km > 0 and kcat > 0
assert N_enzyme > 0 and N_shells > 0
assert V_reactor > 0 and C_S0 > 0
assert D_S > 0
assert P_shell_max > P_shell_min > 0

# SECTION 2: DERIVED PARAMETERS
# --------------------------------------------------------------

# The surface areas and volumes of the protein shell.

A_shell = 4 * np.pi * R_shell**2       # External surface area [m^2]
A_inner = 4 * np.pi * R_inner**2       # Internal surface area [m^2]

V_outer = (4/3) * np.pi * R_shell**3   # Total capsule volume [m^3]
V_inner = (4/3) * np.pi * R_inner**3   # Internal lumen volume [m^3]

V_shell_material = V_outer - V_inner  # Geometric shell-wall volume [m^3]

# The fraction of the external shell surface occupied by open pores.
def porosity(pore_radius, pore_number, shell_radius):

    if np.any(np.asarray(pore_radius) < 0):
        raise ValueError("Pore radius cannot be negative.")

    if pore_number < 0 or shell_radius <= 0:
        raise ValueError("Invalid pore geometry.")

    A_pore = np.pi * pore_radius**2
    A_pores_total = pore_number * A_pore
    A_surface = 4 * np.pi * shell_radius**2

    epsilon = A_pores_total / A_surface

    if np.any(epsilon > 1):
        raise ValueError("Pore area exceeds shell surface area.")

    return epsilon


A_pore = np.pi * r_pore**2             # Individual pore area [m^2]
A_pores_total = N_pores * A_pore       # Total open pore area [m^2]

epsilon = porosity(r_pore, N_pores, R_shell)

# The enzyme inventory and concentration inside one protein shell.
n_enzyme_shell = N_enzyme / N_A        # Enzyme per shell [mol]
C_E_shell = n_enzyme_shell / V_inner   # Internal enzyme concentration [mol/m^3]

Vmax_shell = kcat * C_E_shell          # Internal maximum reaction rate [mol/(m^3*s)]

n_enzyme_total = N_shells * n_enzyme_shell  # Total enzyme inventory [mol]
V_lumen_total = N_shells * V_inner          # Total lumen volume [m^3]
V_capsules_total = N_shells * V_outer       # Total capsule volume [m^3]
A_shell_total = N_shells * A_shell          # Total external area [m^2]

f_lumen = V_lumen_total / V_reactor         # Lumen-to-bulk volume ratio [-]

C_E_free = n_enzyme_total / V_reactor       # Free-enzyme concentration [mol/m^3]
Vmax_free = kcat * C_E_free                 # Free-enzyme maximum rate [mol/(m^3*s)]

# Both systems must have the same total maximum catalytic capacity.
assert np.isclose(
    Vmax_shell * V_lumen_total,
    Vmax_free * V_reactor,
    rtol=1e-10,
    atol=0)

# SECTION 3: FREE-ENZYME MICHAELIS-MENTEN MODEL
# --------------------------------------------------------------

# Michaelis-Menten equation: v = Vmax*C/(Km+C)
def mm_rate(C, Vmax, Km):

    if Vmax < 0 or Km <= 0:
        raise ValueError("Invalid kinetic parameters.")

    C = np.maximum(C, 0)  # Prevent negative rates from numerical errors

    rate = Vmax * C / (Km + C)

    return rate

# Establish the batch reactor mole balances using 1:1 reaction stoichiometry.
def free_batch(t, y, Vmax, Km):

    C_S = y[0]
    C_P = y[1]

    rate = mm_rate(C_S, Vmax, Km)

    dC_S_dt = -rate
    dC_P_dt = rate

    return [dC_S_dt, dC_P_dt]


# Solve the free-enzyme batch reactor using the initial substrate and product concentrations.
def solve_free_batch(C_S0, C_P0, Vmax, Km, times):

    initial = [C_S0, C_P0]

    solution = solve_ivp(
        free_batch,
        (times[0], times[-1]),
        initial,
        args=(Vmax, Km),
        method=solver_method,
        t_eval=times,
        rtol=rtol,
        atol=atol
    )

    if not solution.success:
        raise RuntimeError(solution.message)

    return solution

# Conversion represents the fraction of initial substrate consumed.
def conversion(C_initial, C_current):

    X = (C_initial - C_current) / C_initial

    return X

# Product yield represents the amount of product formed relative to the initial substrate.
def product_yield(C_S_initial, C_P_initial, C_P_current):

    Y = (C_P_current - C_P_initial) / C_S_initial

    return Y

# SECTION 4: FREE-ENZYME MODEL VERIFICATION
# --------------------------------------------------------------

# Integrated Michaelis-Menten equation.
def analytical_mm(t, C_S0, Vmax, Km):

    if t == 0 or Vmax == 0:
        return C_S0

    def residual(z):
        return -C_S0 * np.expm1(z) - Km*z - Vmax*t

    z_lower = -(Vmax*t/Km) - 1
    z_upper = 0

    z = brentq(residual, z_lower, z_upper, xtol=1e-14)

    C_S = C_S0 * np.exp(z)

    return C_S

# The absolute and relative errors between the numerical and analytical solutions.
def calculate_errors(numerical, analytical):

    absolute_error = np.abs(numerical - analytical)

    relative_error = absolute_error / np.maximum(
        np.abs(analytical), 1e-15
    )

    return absolute_error, relative_error

# The substrate-product mole balance.
def mass_balance_error(C_S, C_P, C_S0, C_P0):

    initial_total = C_S0 + C_P0
    current_total = C_S + C_P

    return np.abs(current_total - initial_total)

# The first-order and zero-order limiting cases of Michaelis-Menten kinetics.
def verify_kinetic_limits(Vmax, Km):

    C_low = Km * 1e-6

    rate_low = mm_rate(C_low, Vmax, Km)
    rate_first_order = (Vmax/Km) * C_low

    C_high = Km * 1e6

    rate_high = mm_rate(C_high, Vmax, Km)
    rate_zero_order = Vmax

    assert np.isclose(
        rate_low, rate_first_order, rtol=1e-5, atol=0
    )

    assert np.isclose(
        rate_high, rate_zero_order, rtol=1e-5, atol=0
    )

    print("Michaelis-Menten limiting cases: PASSED")

# SECTION 5: PORE GEOMETRY AND HINDERED DIFFUSION
# --------------------------------------------------------------

# Renkin equation accounts for steric exclusion and hydrodynamic interactions within narrow pores.
def renkin(solute_radius, pore_radius):

    if np.any(np.asarray(solute_radius) < 0):
        raise ValueError("Solute radius cannot be negative.")

    if np.any(np.asarray(pore_radius) <= 0):
        raise ValueError("Pore radius must be positive.")

    lam = solute_radius / pore_radius
    lam = np.clip(lam, 0, 1)

    H = (1 - lam)**2 * (
        1 - 2.1044*lam + 2.089*lam**3 - 0.948*lam**5
    )

    return np.clip(H, 0, 1)

# Effective diffusivity accounts for porosity, molecular hindrance and tortuosity.
def effective_diffusivity(D, epsilon, H, tau):

    if D < 0 or tau < 1:
        raise ValueError("Invalid diffusion parameters.")

    if np.any((np.asarray(epsilon) < 0) | (np.asarray(epsilon) > 1)):
        raise ValueError("Invalid porosity.")

    if np.any((np.asarray(H) < 0) | (np.asarray(H) > 1)):
        raise ValueError("Invalid hindrance factor.")

    D_eff = D * epsilon * H / tau

    return D_eff

# Shell permeability relates the concentration difference across the shell to molecular flux.
def shell_permeability(
    D, solute_radius, pore_radius, pore_number,
    shell_radius, shell_thickness, tau):

    if shell_thickness <= 0:
        raise ValueError("Shell thickness must be positive.")

    epsilon_local = porosity(
        pore_radius, pore_number, shell_radius
    )

    H = renkin(solute_radius, pore_radius)

    D_eff = effective_diffusivity(D, epsilon_local, H, tau)

    permeability = D_eff / shell_thickness

    return permeability

# The hindrance factor for each molecular species.
H_S = renkin(r_S, r_pore)
H_P = renkin(r_P, r_pore)
H_NAD = renkin(r_NAD, r_pore)
H_NADH = renkin(r_NADH, r_pore)

# Species-specific effective diffusivities [m^2/s].
D_eff_S = effective_diffusivity(D_S, epsilon, H_S, tau)
D_eff_P = effective_diffusivity(D_P, epsilon, H_P, tau)
D_eff_NAD = effective_diffusivity(D_NAD, epsilon, H_NAD, tau)
D_eff_NADH = effective_diffusivity(D_NADH, epsilon, H_NADH, tau)

# Species-specific shell permeabilities [m/s].
P_shell_S = D_eff_S / L_shell
P_shell_P = D_eff_P / L_shell
P_shell_NAD = D_eff_NAD / L_shell
P_shell_NADH = D_eff_NADH / L_shell

# SECTION 6: WELL-MIXED SHELL INTERIOR MODEL
# --------------------------------------------------------------

# Steady-state mole balance: P_shell*A_shell*(C_out-C_in) = V_inner*Vmax*C_in/(Km+C_in)
def shell_balance(C_in, C_out, permeability, area, volume, Vmax, Km):

    transport = permeability * area * (C_out - C_in)
    reaction = volume * mm_rate(C_in, Vmax, Km)

    return transport - reaction

# The nonlinear mole balance for the internal substrate concentration.
def internal_concentration(
    C_out, permeability, area, volume, Vmax, Km):

    if C_out < 0 or permeability < 0:
        raise ValueError("Invalid concentration or permeability.")

    if area <= 0 or volume <= 0 or Km <= 0 or Vmax < 0:
        raise ValueError("Invalid geometry or kinetic parameters.")

    if permeability == 0 and Vmax == 0:
        raise ValueError("Internal concentration is undetermined.")

    if C_out == 0 or permeability == 0:
        return 0.0

    if Vmax == 0 or np.isinf(permeability):
        return C_out

    Da = Vmax * volume / (Km * permeability * area)
    s_out = C_out / Km

    def residual(x):

        return 1 - x - Da*x/(1 + s_out*x)

    x = brentq(residual, 0, 1, xtol=1e-15)

    C_in = x * C_out

    return C_in


# The effectiveness factor compares the actual reaction rate to the unrestricted reaction rate.
def effectiveness_factor(C_in, C_out, Vmax, Km):

    actual_rate = mm_rate(C_in, Vmax, Km)
    unrestricted_rate = mm_rate(C_out, Vmax, Km)

    if unrestricted_rate == 0:
        return np.nan

    eta = actual_rate / unrestricted_rate

    return eta


# The internal concentration, reaction rate, transport rate and effectiveness factor.
def shell_state(C_out, permeability, area, volume, Vmax, Km):

    C_in = internal_concentration(
        C_out, permeability, area, volume, Vmax, Km)

    actual_rate = mm_rate(C_in, Vmax, Km)
    unrestricted_rate = mm_rate(C_out, Vmax, Km)

    eta = effectiveness_factor(C_in, C_out, Vmax, Km)

    reaction_total = actual_rate * volume

    if np.isinf(permeability):
        transport_total = reaction_total
    else:
        transport_total = permeability * area * (C_out - C_in)

    return {
        "C_in": C_in,
        "C_out": C_out,
        "actual_rate": actual_rate,
        "unrestricted_rate": unrestricted_rate,
        "eta": eta,
        "reaction_total": reaction_total,
        "transport_total": transport_total}


# Coupling the well-mixed shell model
def encapsulated_batch(
    t, y, permeability, area, volume,
    Vmax, Km, shell_number, reactor_volume):

    C_bulk = max(y[0], 0)

    state = shell_state(
        C_bulk, permeability, area, volume, Vmax, Km)

    rate_per_shell = state["reaction_total"]

    rate_bulk = shell_number * rate_per_shell / reactor_volume

    dC_S_dt = -rate_bulk
    dC_P_dt = rate_bulk

    return [dC_S_dt, dC_P_dt]

# The quasi-steady encapsulated batch reactor.
def solve_encapsulated_batch(permeability):

    initial = [C_S0, C_P0]

    solution = solve_ivp(
        encapsulated_batch,
        t_span,
        initial,
        args=(
            permeability, A_shell, V_inner,
            Vmax_shell, Km, N_shells, V_reactor),
        method=solver_method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol)

    if not solution.success:
        raise RuntimeError(solution.message)

    return solution

# SECTION 7: DIMENSIONLESS ANALYSIS
# --------------------------------------------------------------

# The Thiele modulus compares the reaction rate to the internal diffusion rate.
def thiele_modulus(radius, Vmax, Km, diffusivity):

    phi = radius * np.sqrt(Vmax / (Km * diffusivity))

    return phi


# The shell Damkohler number compares the reaction capacity to the shell transport capacity.
# Da << 1 indicates reaction-controlled behavior.
# Da >> 1 indicates potentially significant shell transport resistance.

def damkohler_number(Vmax, Km, volume, permeability, area):

    if permeability == 0:
        return np.inf

    Da = Vmax * volume / (Km * permeability * area)

    return Da


# The shell Biot number compares transport through the shell to diffusion inside the lumen.
def biot_number(
    permeability, radius, diffusivity, area_outer, area_inner
):

    P_inner = permeability * area_outer / area_inner

    Bi = P_inner * radius / diffusivity

    return Bi


# Analytical internal effectiveness factor for a sphere with first-order kinetics.
def spherical_effectiveness(phi):

    if phi < 0:
        raise ValueError("Thiele modulus cannot be negative.")

    if phi < 1e-3:

        # Series approximation prevents numerical cancellation at small phi.

        eta = 1 - phi**2/15 + 2*phi**4/315

    else:

        eta = (3/phi**2) * (phi/np.tanh(phi) - 1)

    return eta


# Analytical shell effectiveness factor in the first-order kinetic limit.
def first_order_shell_effectiveness(Da):

    return 1 / (1 + Da)


# Calculate all dimensionless parameters using a common set of inputs.
def dimensionless(
    permeability, radius, diffusivity, Vmax,
    Km, volume, area_outer, area_inner):

    phi = thiele_modulus(radius, Vmax, Km, diffusivity)

    Da = damkohler_number(
        Vmax, Km, volume, permeability, area_outer)

    Bi = biot_number(
        permeability, radius, diffusivity, area_outer, area_inner)

    eta_internal = spherical_effectiveness(phi)

    return {
        "phi": phi,
        "Da": Da,
        "Bi": Bi,
        "eta_internal": eta_internal}

P_critical = Vmax_shell * V_inner / (Km * A_shell)


# Sweep shell permeability across the full uncertain range.
def permeability_sweep(C_out):

    P_values = np.logspace(
        np.log10(P_shell_min),
        np.log10(P_shell_max),
        N_sweep)

    C_in_values = np.zeros(N_sweep)
    eta_values = np.zeros(N_sweep)
    Da_values = np.zeros(N_sweep)

    for i, permeability in enumerate(P_values):

        state = shell_state(
            C_out, permeability, A_shell, V_inner, Vmax_shell, Km)

        C_in_values[i] = state["C_in"]
        eta_values[i] = state["eta"]

        Da_values[i] = damkohler_number(
            Vmax_shell, Km, V_inner, permeability, A_shell)

    return P_values, C_in_values, eta_values, Da_values


# Verify the physical limiting cases of the shell model.
def verify_shell_model():

    # Zero permeability must prevent steady-state substrate consumption.
    state_zero = shell_state(
        C_S0, 0, A_shell, V_inner, Vmax_shell, Km)

    assert np.isclose(state_zero["eta"], 0)

    # Infinite permeability must eliminate shell transport resistance.
    state_high = shell_state(
        C_S0, np.inf, A_shell, V_inner, Vmax_shell, Km)

    assert np.isclose(state_high["eta"], 1)

    # Molecules larger than the pore must be excluded.
    assert np.isclose(renkin(r_pore, r_pore), 0)
    assert np.isclose(renkin(2*r_pore, r_pore), 0)
    assert np.isclose(renkin(0, r_pore), 1)

    # Verify that substrate transport equals reaction consumption.
    state = shell_state(
        C_S0, P_shell, A_shell, V_inner, Vmax_shell, Km)

    assert np.isclose(
        state["reaction_total"],
        state["transport_total"],
        rtol=1e-8,
        atol=0)

    # Verify the first-order shell effectiveness factor.
    C_low = Km * 1e-6

    state_low = shell_state(
        C_low, P_shell, A_shell, V_inner, Vmax_shell, Km)

    Da = damkohler_number(
        Vmax_shell, Km, V_inner, P_shell, A_shell)

    assert np.isclose(
        state_low["eta"],
        first_order_shell_effectiveness(Da),
        rtol=1e-5,
        atol=0)

    # Verify the limiting cases of internal diffusion.
    assert np.isclose(spherical_effectiveness(0), 1)

    assert np.isclose(
        spherical_effectiveness(1000),
        3/1000,
        rtol=0.002,
        atol=0)

    # Verify the relationship between the dimensionless parameters.
    groups = dimensionless(
        P_shell, R_inner, D_lumen, Vmax_shell,
        Km, V_inner, A_shell, A_inner)

    assert np.isclose(
        groups["Da"],
        groups["phi"]**2 / (3*groups["Bi"]),
        rtol=1e-10,
        atol=0)

    print("Shell model verification: PASSED")

def run_part1():

    # Solve the free-enzyme reference reactor.
    free_solution = solve_free_batch(
        C_S0, C_P0, Vmax_free, Km, t_eval)

    time = free_solution.t

    C_S_free = free_solution.y[0]
    C_P_free = free_solution.y[1]

    X_free = conversion(C_S0, C_S_free)
    Y_free = product_yield(C_S0, C_P0, C_P_free)

    # Compare the numerical solution against the integrated Michaelis-Menten equation.
    C_S_analytical = np.array([
        analytical_mm(t, C_S0, Vmax_free, Km)
        for t in time])

    absolute_error, relative_error = calculate_errors(
        C_S_free, C_S_analytical
    )

    balance_error = mass_balance_error(
        C_S_free, C_P_free, C_S0, C_P0
    )

    assert np.allclose(
        C_S_free, C_S_analytical, rtol=1e-6, atol=1e-8
    )

    assert np.max(balance_error) < 1e-8
    assert np.min(free_solution.y) >= -atol

    print("Free-enzyme analytical verification: PASSED")
    print("Free-enzyme mass balance: PASSED")

    verify_kinetic_limits(Vmax_free, Km)
    verify_shell_model()

    # Evaluate the shell using independent and geometry-derived permeabilities.
    state_test = shell_state(
        C_S0, P_shell, A_shell, V_inner, Vmax_shell, Km
    )

    state_geometry = shell_state(
        C_S0, P_shell_S, A_shell, V_inner, Vmax_shell, Km
    )

    # Solve the encapsulated batch reactor.
    encapsulated_solution = solve_encapsulated_batch(P_shell)

    C_S_encapsulated = encapsulated_solution.y[0]
    C_P_encapsulated = encapsulated_solution.y[1]

    X_encapsulated = conversion(C_S0, C_S_encapsulated)

    Y_encapsulated = product_yield(
        C_S0, C_P0, C_P_encapsulated
    )

    # Verify that encapsulation does not increase the intrinsic catalytic capacity.
    assert np.all(C_S_encapsulated >= C_S_free - 1e-7)

    assert np.max(
        mass_balance_error(
            C_S_encapsulated, C_P_encapsulated, C_S0, C_P0
        )
    ) < 1e-8

    assert np.min(encapsulated_solution.y) >= -atol

    # Calculate dimensionless parameters.
    groups = dimensionless(
        P_shell, R_inner, D_lumen, Vmax_shell,
        Km, V_inner, A_shell, A_inner
    )

    phi = groups["phi"]
    Da = groups["Da"]
    Bi = groups["Bi"]
    eta_internal = groups["eta_internal"]

    # Perform the permeability sweep and verify that effectiveness increases with permeability.
    P_values, C_in_values, eta_values, Da_values = permeability_sweep(C_S0)

    assert np.all(np.diff(eta_values) >= -1e-12)
    assert np.all((eta_values >= 0) & (eta_values <= 1 + 1e-12))

    # Print shell geometry and enzyme properties.
    print("\nSHELL GEOMETRY")
    print(f"Outer radius: {R_shell*1e9:.3f} nm")
    print(f"Inner radius: {R_inner*1e9:.3f} nm")
    print(f"Pore radius: {r_pore*1e9:.3f} nm")
    print(f"Number of pores: {N_pores}")
    print(f"Shell porosity: {epsilon*100:.4f}%")
    print(f"Internal volume: {V_inner:.6e} m^3")

    print("\nENZYME PROPERTIES")
    print(f"Enzymes per shell: {N_enzyme}")
    print(f"Internal enzyme concentration: {C_E_shell:.6e} mM")
    print(f"Total enzyme inventory: {n_enzyme_total:.6e} mol")
    print(f"Internal Vmax: {Vmax_shell:.6e} mol/(m^3*s)")
    print(f"Free-enzyme Vmax: {Vmax_free:.6e} mol/(m^3*s)")

    print("\nSHELL DIFFUSION")
    print(f"Substrate hindrance: {H_S:.6f}")
    print(f"Product hindrance: {H_P:.6f}")
    print(f"NAD+ hindrance: {H_NAD:.6f}")
    print(f"NADH hindrance: {H_NADH:.6f}")
    print(f"Substrate permeability: {P_shell_S:.6e} m/s")
    print(f"Product permeability: {P_shell_P:.6e} m/s")
    print(f"NAD+ permeability: {P_shell_NAD:.6e} m/s")
    print(f"NADH permeability: {P_shell_NADH:.6e} m/s")

    print("\nSTEADY-STATE SHELL MODEL")
    print(f"Bulk substrate concentration: {C_S0:.6f} mM")
    print(f"Internal substrate concentration: {state_test['C_in']:.6f} mM")
    print(f"Test permeability: {P_shell:.6e} m/s")
    print(f"Test effectiveness: {state_test['eta']:.6f}")
    print(f"Geometry-derived effectiveness: {state_geometry['eta']:.6f}")

    print("\nDIMENSIONLESS ANALYSIS")
    print(f"Thiele modulus: {phi:.6e}")
    print(f"Damkohler number: {Da:.6e}")
    print(f"Biot number: {Bi:.6e}")
    print(f"Internal effectiveness: {eta_internal:.9f}")
    print(f"Critical permeability: {P_critical:.6e} m/s")

    print("\nBATCH REACTOR PERFORMANCE")
    print(f"Free-enzyme conversion: {X_free[-1]*100:.4f}%")
    print(f"Free-enzyme yield: {Y_free[-1]*100:.4f}%")
    print(f"Encapsulated conversion: {X_encapsulated[-1]*100:.4f}%")
    print(f"Encapsulated yield: {Y_encapsulated[-1]*100:.4f}%")

    print("\nNUMERICAL VERIFICATION")
    print(f"Maximum absolute error: {np.max(absolute_error):.6e}")
    print(f"Maximum relative error: {np.max(relative_error):.6e}")
    print(f"Maximum mass balance error: {np.max(balance_error):.6e}")

    # Figure 1: Free-enzyme substrate and product concentrations versus time.
    plt.figure(figsize=(8, 5))
    plt.plot(time/60, C_S_free, label="Substrate")
    plt.plot(time/60, C_P_free, label="Product")
    plt.xlabel("Time (min)")
    plt.ylabel("Concentration (mM)")
    plt.title("Free-Enzyme Batch Reactor")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 2: Analytical and numerical substrate concentration comparison.
    plt.figure(figsize=(8, 5))
    plt.plot(time/60, C_S_free, label="Numerical")
    plt.plot(time/60, C_S_analytical, "--", label="Analytical")
    plt.xlabel("Time (min)")
    plt.ylabel("Substrate Concentration (mM)")
    plt.title("Analytical and Numerical Verification")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 3: Absolute numerical error versus time.
    plt.figure(figsize=(8, 5))
    plt.plot(time/60, absolute_error)
    plt.xlabel("Time (min)")
    plt.ylabel("Absolute Error (mM)")
    plt.title("Numerical Solution Error")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 4: Renkin hindrance factor versus pore radius.
    pore_values = np.linspace(0.1e-9, 5e-9, N_sweep)
    plt.figure(figsize=(8, 5))
    plt.plot(pore_values*1e9, renkin(r_S, pore_values), label="Substrate")
    plt.plot(pore_values*1e9, renkin(r_P, pore_values), "--", label="Product")
    plt.plot(pore_values*1e9, renkin(r_NAD, pore_values), label="NAD+")
    plt.plot(pore_values*1e9, renkin(r_NADH, pore_values), "--", label="NADH")
    plt.xlabel("Pore Radius (nm)")
    plt.ylabel("Renkin Hindrance Factor")
    plt.title("Molecular Hindrance Through Protein-Shell Pores")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 5: Shell permeability versus pore radius for different pore numbers.
    plt.figure(figsize=(8, 5))
    for pore_number in [20, 60, 100]:

        epsilon_sweep = pore_number * np.pi * pore_values**2 / A_shell
        valid_radii = pore_values[epsilon_sweep <= 1]

        P_geometry = shell_permeability(
            D_S, r_S, valid_radii, pore_number,
            R_shell, L_shell, tau)

        plt.plot(
            valid_radii*1e9, P_geometry,
            label=f"{pore_number} pores")

    plt.xlabel("Pore Radius (nm)")
    plt.ylabel("Shell Permeability (m/s)")
    plt.title("Effect of Pore Geometry on Shell Permeability")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 6: Effectiveness factor versus shell permeability.
    plt.figure(figsize=(8, 5))
    plt.semilogx(P_values, eta_values, label="Michaelis-Menten")
    plt.axvline(P_critical, linestyle="--", label="Da = 1")
    plt.axvline(P_shell, linestyle=":", label="Test permeability")
    plt.xlabel("Shell Permeability (m/s)")
    plt.ylabel("Effectiveness Factor")
    plt.title("Effectiveness Factor vs. Shell Permeability")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 7: Effectiveness factor versus Damkohler number.
    order = np.argsort(Da_values)

    plt.figure(figsize=(8, 5))
    plt.semilogx(
        Da_values[order], eta_values[order],
        label="Michaelis-Menten")
    plt.semilogx(
        Da_values[order],
        first_order_shell_effectiveness(Da_values[order]),
        "--", label="First-order limit")
    
    plt.axvline(1, linestyle=":", label="Da = 1")
    plt.xlabel("Shell Damkohler Number")
    plt.ylabel("Effectiveness Factor")
    plt.title("Effectiveness Factor vs. Damkohler Number")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 8: Internal substrate concentration versus shell permeability.
    plt.figure(figsize=(8, 5))
    plt.semilogx(P_values, C_in_values, label="Internal substrate")
    plt.axhline(C_S0, linestyle="--", label="Bulk substrate")
    plt.xlabel("Shell Permeability (m/s)")
    plt.ylabel("Substrate Concentration (mM)")
    plt.title("Internal Substrate Concentration")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 9: Internal effectiveness factor versus Thiele modulus.
    phi_values = np.logspace(-5, 3, N_sweep)

    eta_phi = np.array([
        spherical_effectiveness(value)
        for value in phi_values])

    plt.figure(figsize=(8, 5))
    plt.semilogx(phi_values, eta_phi, label="First-order spherical solution")
    plt.axvline(phi, linestyle="--", label="Baseline Thiele modulus")
    plt.xlabel("Thiele Modulus")
    plt.ylabel("Internal Effectiveness Factor")
    plt.title("Internal Diffusion Effectiveness")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Figure 10: Free-enzyme versus encapsulated enzyme conversion.
    plt.figure(figsize=(8, 5))
    plt.plot(time/60, X_free*100, label="Free enzyme")
    plt.plot(time/60, X_encapsulated*100, "--", label="Encapsulated enzyme")
    plt.xlabel("Time (min)")
    plt.ylabel("Substrate Conversion (%)")
    plt.title("Free vs. Encapsulated Enzyme Performance")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return {
        "free_solution": free_solution,
        "encapsulated_solution": encapsulated_solution,
        "shell_state": state_test,
        "geometry_state": state_geometry,
        "dimensionless": groups,
        "permeability_values": P_values,
        "effectiveness_values": eta_values,
        "damkohler_values": Da_values}


# PART 2
# --------------------------------------------------------------

# SECTION 8: RADIAL FINITE-VOLUME DISCRETIZATION
# --------------------------------------------------------------

# Divide the spherical lumen into radial control volumes.
# Calculate radial cell centers, cell boundaries, surface areas and volumes.
# Establish the spherical reaction-diffusion equation.
# Calculate diffusive fluxes between neighboring cells.
# Apply Michaelis-Menten kinetics within each control volume.
# Apply the symmetry boundary condition at r = 0.
# Apply the permeable-shell boundary condition at r = R_inner.
# Account for diffusion between the outermost cell center and the internal shell surface.
# Establish the initial radial concentration profile.

# SECTION 9: TRANSIENT REACTION-DIFFUSION MODEL
# --------------------------------------------------------------

# Construct the system of finite-volume differential equations.
# Couple the internal concentration profile to the bulk batch reactor.
# Include substrate accumulation inside the lumen.
# Solve the system using solve_ivp with BDF or Radau.
# Calculate the volume-integrated reaction rate.
# Calculate the transient effectiveness factor.
# Calculate substrate conversion and product yield.
# Generate radial concentration profiles at different times.
# Compare the steady-state results against the well-mixed model in Section 6.
# Verify conservation across the bulk liquid and internal lumen.

# SECTION 10: EXTERNAL MASS TRANSFER
# --------------------------------------------------------------

# Calculate the characteristic particle diameter.
# Establish a justified relationship between mixing rate and particle-fluid relative velocity.
# Calculate the Reynolds number: Re = rho*u*d_particle/mu.
# Calculate the Schmidt number: Sc = mu/(rho*D).
# Select an appropriate Sherwood number correlation.
# Calculate the external mass-transfer coefficient: kL = Sh*D/d_particle.
# Combine external film and shell transport resistances using the same reference area.
# Modify the external boundary condition.
# Investigate the effect of mixing rate on effectiveness.
# Verify the stagnant-sphere limiting case: Sh = 2.

# SECTION 11: COFACTOR AND PRODUCT INHIBITION
# --------------------------------------------------------------

# Introduce NAD+ as a second substrate and NADH as a reaction product.
# Establish the complete reaction stoichiometry.
# Select a literature-supported two-substrate kinetic expression.
# Incorporate product inhibition and reaction reversibility where supported by data.
# Reuse the species-specific permeabilities calculated in Section 5.
# Establish coupled radial and bulk species balances.
# Solve the multi-species reaction-diffusion equations.
# Verify substrate-product and NAD+-NADH conservation.

# SECTION 12: TEMPERATURE EFFECTS
# --------------------------------------------------------------

# Verify the reference temperature of the kinetic parameters.
# Introduce Arrhenius temperature dependence for kcat.
# Introduce temperature-dependent enzyme deactivation.
# Establish the active-enzyme balance: dE/dt = -kd(T)*E.
# Update Vmax as active enzyme concentration changes.
# Incorporate temperature-dependent molecular diffusivities.
# Investigate the trade-off between catalytic activity and enzyme stability.
# Verify that the reference temperature reproduces the baseline model.

# PART 3
# --------------------------------------------------------------
# SECTION 13: ENZYME LEACHING
# --------------------------------------------------------------

# Compare the enzyme hydrodynamic radius against the pore radius.
# Establish a size-dependent enzyme leaching model.
# Account for possible enzyme attachment to scaffold proteins.
# Calculate enzyme loss during the washing step.
# Update the retained and active enzyme inventories separately.
# Recalculate internal enzyme concentration and Vmax.

# SECTION 14: PORE BLOCKAGE
# --------------------------------------------------------------

# Establish a physically justified pore blockage mechanism.
# Define the foulant concentration or exposure variable.
# Calculate the remaining fraction of open pores.
# Update the effective shell porosity.
# Recalculate species-specific shell permeabilities.
# Prevent negative porosity and permeability.
# Investigate the effect of pore blockage on catalytic activity.

# SECTION 15: REPEATED-BATCH REACTOR
# --------------------------------------------------------------

# Establish arrays for storing the results of each batch cycle.
# Create a loop over the total number of cycles.
# Introduce fresh bulk substrate and cofactor at the beginning of each cycle.
# Solve the complete reaction-diffusion model.
# Calculate final conversion, yield and product formation.
# Update active enzyme concentration.
# Perform the washing step and calculate enzyme leaching.
# Update pore blockage and shell permeability.
# Carry retained internal species into the next cycle.
# Calculate cumulative product formation.
# Generate activity and conversion versus cycle number.

# SECTION 16: NUMERICAL VERIFICATION
# --------------------------------------------------------------

# Perform radial mesh convergence using 50, 100 and 200 control volumes.
# Perform solver tolerance convergence.
# Verify nonnegative species concentrations.
# Verify mass balance closure throughout the reaction and washing steps.
# Verify zero-deactivation, zero-leaching and zero-blockage limiting cases.
# Compare the results against the Part 1 and Part 2 reference models.

# SECTION 17: PARAMETER ESTIMATION AND VALIDATION
# --------------------------------------------------------------

# Collect experimental data and kinetic parameters from literature.
# Establish parameter uncertainty ranges.
# Define experimental and predicted responses.
# Calculate the residuals between experimental and predicted results.
# Perform parameter estimation using scipy.optimize.least_squares.
# Evaluate parameter identifiability.
# Compare free and encapsulated enzyme kinetics.
# Compare different P22 shell morphologies.
# Validate the model using independent experimental data.
# Treat parameters without sufficient experimental data as uncertain scenario variables.

# SECTION 18: SENSITIVITY ANALYSIS
# --------------------------------------------------------------

# Establish lower and upper bounds for uncertain parameters.
# Investigate sensitivity to pore radius, pore number and shell radius.
# Investigate sensitivity to enzyme loading and molecular diffusivity.
# Investigate sensitivity to shell permeability and tortuosity.
# Investigate sensitivity to mixing rate and temperature.
# Investigate sensitivity to enzyme deactivation, leaching and pore blockage.
# Perform one-at-a-time sensitivity analysis.
# Perform global Sobol sensitivity analysis if feasible.
# Rank parameters according to their influence on model outputs.
# Avoid independently varying parameters that are mathematically dependent unless justified.

# SECTION 19: PORE DESIGN OPTIMIZATION
# --------------------------------------------------------------

# Establish a range of candidate pore radii and pore numbers.
# Construct a two-dimensional parameter grid.
# Reject geometries that violate pore-area or structural constraints.
# Simulate repeated-batch operation for each feasible pore design.
# Calculate cumulative product formation and enzyme retention.
# Identify the optimal feasible pore design.
# Generate a heatmap of cumulative product formation versus pore radius and pore number.

# SECTION 20: FINAL FIGURES AND RESULTS
# --------------------------------------------------------------

# Compile free-enzyme and encapsulated enzyme kinetic results.
# Compile analytical and numerical verification results.
# Generate radial concentration profiles.
# Generate effectiveness factor versus Thiele modulus and Damkohler number.
# Generate temperature-dependent performance curves.
# Generate activity, conversion and yield versus cycle number.
# Generate cumulative product formation versus cycle number.
# Generate parameter sensitivity rankings.
# Generate the optimal pore design heatmap.
# Export numerical results and figures.
# Document model assumptions, parameter sources and limitations.

if __name__ == "__main__":

    results_part1 = run_part1()
    # results_part2 = run_part2()
    # results_part3 = run_part3()
