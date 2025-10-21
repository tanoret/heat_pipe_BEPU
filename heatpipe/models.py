from math import pi, sin, sqrt, log
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple, Dict

GRAV = 980.0  # cm/s^2
PI = pi

class Regime(Enum):
    LAMINAR = "laminar"
    TURBULENT = "turbulent"

@dataclass
class SectionLengths:
    L_e: float  # evaporator [cm]
    L_a: float  # adiabatic [cm]
    L_c: float  # condenser [cm]

@dataclass
class Geometry:
    radius_cm: float # Inner radius
    t_a_cm: float # Annulus thickness
    t_w_cm: float # Wick thickness
    effective_pore_radius_cm: float = 0.0
    wavelength_cm: Optional[float] = None # Characteristic length for entrainment
    passages: int = 1 

@dataclass
class FlowFlags:
    vapor_regime_evap: Regime = Regime.LAMINAR
    vapor_regime_adiab: Regime = Regime.LAMINAR
    vapor_regime_cond: Regime = Regime.LAMINAR

def cross_section_area(radius_cm: float) -> float:
    return PI * radius_cm**2

def vapor_cross_section_area(radius_cm: float, t_a_cm: float, t_w_cm: float) -> float:
    return PI * (radius_cm-t_a_cm-t_w_cm)**2

def wick_cross_section_area(radius_cm: float, t_a_cm: float, t_w_cm: float) -> float:
    return PI * (radius_cm**2 - (radius_cm-t_a_cm)**2) # Assume flow in the annulus

def reynolds(m_dot_g_s: float, mu_poise: float, rho_g_cm3: float, area_cm2: float, d_h_cm: float) -> float:
    V = m_dot_g_s / (rho_g_cm3 * area_cm2)
    return rho_g_cm3 * V * d_h_cm / mu_poise

def fanning_friction_factor(Re: float, regime: Regime) -> float:
    if Re <= 0:
        return 0.0
    if regime == Regime.LAMINAR:
        return 16.0 / Re
    return 0.0791 * (Re ** -0.25)

def velocity_profile_factor(regime: Regime) -> float:
    return 1.234 if regime == Regime.LAMINAR else 2.22

def dp_friction_fanning(f: float, rho: float, V: float, L_cm: float, D_cm: float) -> float:
    return 2.0 * f * rho * (V**2) * (L_cm / D_cm)

def dp_inertial(A_prof: float, rho_v: float, V_exit: float) -> float:
    return A_prof * rho_v * (V_exit**2)

def dp_hydrostatic(rho_l: float, delta_z_cm: float, theta_deg: float) -> float:
    return rho_l * GRAV * delta_z_cm * sin(theta_deg * PI / 180.0)

def capillary_pressure_drop(sigma: float, r_eff_cm: float) -> float:
    return 0.0 if r_eff_cm <= 0 else 2.0 * sigma / r_eff_cm

# --- Fanno relations (compressible adiabatic with friction) ---
def fanno_F(M: float, gamma: float) -> float:
    if M <= 0:
        return float("inf")
    term1 = (1.0 - M**2) / (gamma * M**2)
    term2 = ((gamma + 1.0) / (2.0 * gamma)) * log(((gamma + 1.0) * M**2) / (2.0 + (gamma - 1.0) * M**2))
    return term1 + term2

def fmax_to_sonic(M1: float, gamma: float) -> float:
    return max(fanno_F(M1, gamma), 0.0)

def pressr_Fanno(M: float, gamma: float) -> float:
    from math import sqrt
    num = (gamma + 1.0) / 2.0
    den = 1.0 + (gamma - 1.0) * 0.5 * M**2
    return (1.0 / M) * sqrt(num / den)

def fanno_downstream_mach(M1: float, gamma: float, four_f_L_over_D: float) -> Tuple[bool, float]:
    Fmax = fmax_to_sonic(M1, gamma)
    if four_f_L_over_D >= Fmax - 1e-9:
        return True, 1.0
    target = fanno_F(M1, gamma) - four_f_L_over_D
    a, b = 1e-6, 0.999999
    for _ in range(60):
        mid = 0.5 * (a + b)
        Fmid = fanno_F(mid, gamma)
        if Fmid > target:
            a = mid
        else:
            b = mid
    return False, 0.5 * (a + b)

def adiabatic_section_pressure_outlet(P_in: float, gamma: float, f: float, L_cm: float, D_cm: float, M_in: float) -> Tuple[float, bool, float]:
    four_f_L_over_D = 4.0 * f * (L_cm / D_cm)
    choked, M_out = fanno_downstream_mach(M_in, gamma, four_f_L_over_D)
    if choked:
        M_out = 1.0
    PR = pressr_Fanno(M_in, gamma) / pressr_Fanno(M_out, gamma)
    return P_in / PR, choked, M_out

# --- Entrainment limits ---
def entrainment_velocity_weber(sigma: float, rho_v: float, wavelength_cm: float) -> float:
    from math import pi, sqrt
    return sqrt(2.0 * pi * sigma / (rho_v * wavelength_cm))

def flooding_coefficient_Bo(gamma_bond: float) -> float:
    from math import tanh
    return tanh(0.5 * gamma_bond ** 0.25)

def entrainment_limit_wickless_smooth(rho_v: float, rho_l: float, sigma: float, g: float = GRAV) -> float:
    from math import sqrt
    Bo = g * (rho_l - rho_v) / sigma
    Cw = flooding_coefficient_Bo(Bo)
    return Cw * sqrt(sigma / max(rho_v, 1e-9))

@dataclass
class PressureBreakdown:
    dp_fric_evap: float
    dp_inert_evap: float
    dp_fric_cond: float
    dp_adiabatic: float
    dp_hydro: float
    dp_capillary: float
    choked_adiab: bool
    M_out_adiab: float

def pressure_breakdown_cgs(
    lengths: SectionLengths,
    geom: Geometry,
    fluids: Dict[str, float],
    flags: FlowFlags,
    m_dot_g_s: float,
    theta_deg: float = 0.0,
) -> PressureBreakdown:
    from math import sqrt
    R = geom.radius_cm
    D = 2.0 * R
    AV = PI * R**2
    V_exit = m_dot_g_s / (fluids["rho_v"] * AV)
    Re_e = reynolds(m_dot_g_s, fluids["mu_v"], fluids["rho_v"], AV, D)
    f_e = fanning_friction_factor(Re_e, flags.vapor_regime_evap)
    Re_c = reynolds(m_dot_g_s, fluids["mu_v"], fluids["rho_v"], AV, D)
    f_c = fanning_friction_factor(Re_c, flags.vapor_regime_cond)
    dp_fric_evap = dp_friction_fanning(f_e, fluids["rho_v"], V_exit, lengths.L_e, D)
    dp_fric_cond = dp_friction_fanning(f_c, fluids["rho_v"], V_exit, lengths.L_c, D)
    A_prof = velocity_profile_factor(flags.vapor_regime_evap)
    dp_inert = dp_inertial(A_prof, fluids["rho_v"], V_exit)

    delta_z = lengths.L_e + lengths.L_a + lengths.L_c
    dp_hyd = dp_hydrostatic(fluids["rho_l"], delta_z, theta_deg)

    dp_adiab = 0.0
    choked = False
    M_out = 0.0
    if lengths.L_a > 0.0:
        a_in = sqrt(fluids["gamma"] * fluids["p_sat"] / max(fluids["rho_v"], 1e-9))
        M_in = V_exit / max(a_in, 1e-9)
        Re_a = reynolds(m_dot_g_s, fluids["mu_v"], fluids["rho_v"], AV, D)
        f_a = fanning_friction_factor(Re_a, flags.vapor_regime_adiab)
        P2, choked, M_out = adiabatic_section_pressure_outlet(fluids["p_sat"], fluids["gamma"], f_a, lengths.L_a, D, max(min(M_in, 0.99), 1e-6))
        dp_adiab = max(fluids["p_sat"] - P2, 0.0)

    dp_cap = capillary_pressure_drop(fluids["sigma"], geom.effective_pore_radius_cm)

    return PressureBreakdown(
        dp_fric_evap=dp_fric_evap,
        dp_inert_evap=dp_inert,
        dp_fric_cond=dp_fric_cond,
        dp_adiabatic=dp_adiab,
        dp_hydro=dp_hyd,
        dp_capillary=dp_cap,
        choked_adiab=choked,
        M_out_adiab=M_out
    )

@dataclass
class Limits:
    q_capillary_W: Optional[float]
    q_sonic_W: Optional[float]
    q_entrainment_W: Optional[float]
    q_viscous_W: Optional[float]

def estimate_limits_cgs(
    lengths: SectionLengths,
    geom: Geometry,
    fluids: Dict[str, float],
    flags: FlowFlags,
) -> Limits:
    from math import sqrt
    R = geom.radius_cm
    AV_cm2 = PI * R**2
    AV_m2 = AV_cm2 * 1e-4
    h_fg_J_kg = fluids["h_fg"] * 1e3
    L_eff = ((lengths.L_e+lengths.L_c)/2 + lengths.L_a)/100 # effective length [m]

    # Viscous
    q_visc = AV_m2**2*h_fg_J_kg*fluids["rho_v"]*1000*(fluids["p_sat"]*0.1)/16/PI/(fluids["mu_v"]*0.1)/L_eff

    # Sonic (characteristic)
    a = sqrt(fluids["gamma"] * fluids["p_sat"] / max(fluids["rho_v"], 1e-9))
    G = fluids["rho_v"] * a  # g/(cm^2 s)
    q_sonic = G * h_fg_J_kg * AV_m2 * 1e-3

    # Capillary (very rough, use available dp to infer a characteristic V)
    dp_avail = capillary_pressure_drop(fluids["sigma"], geom.effective_pore_radius_cm)
    if dp_avail <= 0:
        q_cap = None
    else:
        D_cm = 2.0 * R
        f_nom = 0.02
        L_char = max(1.0, lengths.L_e + lengths.L_a + lengths.L_c)
        Vchar = sqrt(dp_avail * D_cm / (2.0 * f_nom * max(fluids["rho_v"], 1e-9) * L_char))
        Gcap = fluids["rho_v"] * Vchar
        q_cap = Gcap * h_fg_J_kg * AV_m2 * 1e-3

    # Entrainment (if wavelength known)
    q_ent = None
    if geom.wavelength_cm:
        Vcrit = sqrt(2.0 * pi * fluids["sigma"] / (max(fluids["rho_v"], 1e-9) * geom.wavelength_cm))
        Gent = fluids["rho_v"] * Vcrit
        q_ent = Gent * h_fg_J_kg * AV_m2 * 1e-3

    return Limits(q_capillary_W=q_cap, q_sonic_W=q_sonic, q_entrainment_W=q_ent, q_viscous_W=q_visc)
