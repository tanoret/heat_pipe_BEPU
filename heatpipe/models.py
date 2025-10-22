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
    elif Re < 2000:
        return 16.0 / Re
    elif Re < 20000:
        return 0.0791 / (Re**0.25)
    else:
        return 0.046 / (Re**0.2)

def velocity_profile_factor(regime: Regime) -> float:
    return 1.234 if regime == Regime.LAMINAR else 2.22

def dp_friction_fanning(f: float, rho: float, V: float, L_cm: float, D_cm: float) -> float:
    return 2.0 * f * rho * (V**2) * (L_cm / D_cm)

def dp_inertial(A_prof: float, rho_v: float, V_exit: float) -> float:
    return A_prof * rho_v * (V_exit**2)

def dp_hydrostatic(rho_l: float, delta_z_cm: float, theta_deg: float) -> float:
    return rho_l * GRAV * delta_z_cm * sin(theta_deg * PI / 180.0)

def dp_liquid(mdot: float, t_a:float, rho: float, mu: float, L_eff: float, A_l: float) -> float:
    K = t_a**2 / 12
    return mu*mdot/(rho*t_a**2*A_l)

def max_capillary_pressure(sigma: float, r_eff_cm: float) -> float:
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
    dp_tot: float

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

    dp_cap = max_capillary_pressure(fluids["sigma"], geom.effective_pore_radius_cm)
    dp_tot = dp_inert+dp_fric_evap+dp_adiab+dp_fric_cond+dp_hyd
    return PressureBreakdown(
        dp_fric_evap=dp_fric_evap,
        dp_inert_evap=dp_inert,
        dp_fric_cond=dp_fric_cond,
        dp_adiabatic=dp_adiab,
        dp_hydro=dp_hyd,
        dp_capillary=dp_cap,
        choked_adiab=choked,
        M_out_adiab=M_out,
        dp_tot = dp_tot,
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
    theta_deg: float = 0.0,
    q_max_search_W: float = 5.0e4,
    tol_rel: float = 1e-3,
    ) -> Limits:

    """Return classical transport limits for the supplied configuration.

    The routine evaluates the sonic, entrainment and capillary limits by
    solving simplified balance relations that are common in the heat-pipe
    literature (e.g. Faghri, 1995; Reay & Kew, 2006).  All calculations are
    performed in SI internally to avoid unit drift before the results are
    reported in Watts.
    """

    from math import sqrt

    # --- common geometric conversions ---
    R_cm = geom.radius_cm
    area_cm2 = PI * R_cm**2
    area_m2 = area_cm2 * 1.0e-4
    length_total_cm = lengths.L_e + lengths.L_a + lengths.L_c
    length_total_m = length_total_cm * 1.0e-2

    # --- fluid property conversions ---
    rho_v_kg_m3 = fluids["rho_v"] * 1.0e3
    rho_l_kg_m3 = fluids["rho_l"] * 1.0e3
    mu_l_Pa_s = fluids["mu_l"] * 0.1
    p_sat_Pa = fluids["p_sat"] * 0.1
    sigma_N_m = fluids["sigma"] * 1.0e-3
    h_fg_J_kg = fluids["h_fg"] * 1.0e3
    h_fg_J_g = fluids["h_fg"]  # (kJ/kg == J/g)

    # --- sonic limit ---
    if rho_v_kg_m3 <= 0 or p_sat_Pa <= 0:
        q_sonic = None
    else:
        a = sqrt(max(fluids["gamma"], 1e-9) * p_sat_Pa / max(rho_v_kg_m3, 1e-12))
        mass_flux_kg_m2_s = rho_v_kg_m3 * a
        q_sonic = mass_flux_kg_m2_s * h_fg_J_kg * area_m2

    # --- viscous (liquid) limit ---
    # Treat the wick as a simple capillary annulus with permeability ~ r_p^2 / 8
    # (see e.g. Faghri 1995, Eq. 4.53).  Use Darcy's law with gravity head.
    if geom.effective_pore_radius_cm > 0:
        r_eff_m = geom.effective_pore_radius_cm * 1.0e-2
        permeability_m2 = r_eff_m**2 / 8.0
        R_outer_m = R_cm * 1.0e-2
        R_inner_m = max(R_cm - geom.t_a_cm, 0.0) * 1.0e-2
        wick_area_m2 = PI * max(R_outer_m**2 - R_inner_m**2, 0.0)
        dp_available_Pa = max_capillary_pressure(sigma_N_m * 1.0e3, geom.effective_pore_radius_cm) * 0.1
        dp_gravity_Pa = (rho_l_kg_m3 - rho_v_kg_m3) * 9.80 * length_total_m * abs(
            sin(theta_deg * PI / 180.0)
        )
        dp_driving_Pa = max(dp_available_Pa - dp_gravity_Pa, 0.0)
        if dp_driving_Pa <= 0:
            q_visc = None
        else:
            volumetric_flow_m3_s = permeability_m2 * wick_area_m2 * dp_driving_Pa / (
                mu_l_Pa_s * max(length_total_m, 1e-9)
            )
            mass_flow_kg_s = volumetric_flow_m3_s * rho_l_kg_m3
            q_visc = mass_flow_kg_s * h_fg_J_kg
    else:
        q_visc = None

    # --- capillary limit (solve m_dot so vapor Δp <= capillary Δp) ---
    dp_cap_cgs = max_capillary_pressure(fluids["sigma"], geom.effective_pore_radius_cm)
    if dp_cap_cgs <= 0:
        q_cap = None
    else:
        def dp_error(q_watts: float) -> float:
            if q_watts <= 0:
                return -dp_cap_cgs
            m_dot_g_s = q_watts / max(h_fg_J_g, 1e-9)
            pb = pressure_breakdown_cgs(
                lengths, geom, fluids, flags, m_dot_g_s=m_dot_g_s, theta_deg=theta_deg
            )
            return pb.dp_tot - pb.dp_capillary

        lo, hi = 0.0, q_max_search_W
        err_lo = dp_error(lo)
        err_hi = dp_error(hi)
        if err_hi <= 0:
            q_cap = hi
        else:
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                err_mid = dp_error(mid)
                if abs(err_mid) <= tol_rel * dp_cap_cgs:
                    q_cap = mid
                    break
                if err_mid > 0:
                    hi = mid
                    err_hi = err_mid
                else:
                    lo = mid
                    err_lo = err_mid
            else:
                q_cap = 0.5 * (lo + hi)

    # --- entrainment limit ---
    if geom.wavelength_cm:
        wavelength_m = geom.wavelength_cm * 1.0e-2
        if wavelength_m <= 0 or rho_v_kg_m3 <= 0:
            q_ent = None
        else:
            Vcrit = sqrt(2.0 * PI * sigma_N_m / (rho_v_kg_m3 * wavelength_m))
            mass_flux_kg_m2_s = rho_v_kg_m3 * Vcrit
            q_ent = mass_flux_kg_m2_s * h_fg_J_kg * area_m2
    else:
        q_ent = None

    return Limits(q_capillary_W=q_cap, q_sonic_W=q_sonic, q_entrainment_W=q_ent, q_viscous_W=q_visc)
