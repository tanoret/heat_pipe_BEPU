from math import pi, sin, sqrt, log
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Tuple, Dict
from heatpipe.io import props_for
from heatpipe.fluids import Fluid, tsat_from_p_cgs

GRAV = 980.0  # Gravitational acceleration [cm/s^2]
PI = pi 
RBAR = 8.3144621 # Universal gas constant [J/mol-K]

@dataclass
class SectionLengths:
    L_e: float  # evaporator length [cm]
    L_a: float  # adiabatic length [cm]
    L_c: float  # condenser length [cm]

@dataclass
class Geometry:
    radius_cm: float # Inner radius [cm]
    t_a_cm: float # Annulus thickness [cm]
    t_w_cm: float # Wick thickness [cm]
    effective_pore_radius_cm: float = 0.0 # Effective pore radius [cm]
    wavelength_cm: Optional[float] = None # Characteristic length for entrainment [cm]

def wick_cross_section_area(radius_cm: float, t_a_cm: float, t_w_cm: float) -> float: 
    # Will be adding different wick designs, currently annular wick only
    return PI * (radius_cm**2 - (radius_cm-t_a_cm)**2) # Assume flow in the annulus

def reynolds(m_dot_g_s: float, mu_poise: float, rho_g_cm3: float, area_cm2: float, d_h_cm: float) -> float: 
    # Axial Reynolds number 
    V = m_dot_g_s / (rho_g_cm3 * area_cm2)
    return rho_g_cm3 * V * d_h_cm / mu_poise

def fanning_friction_factor(Re: float) -> float:
    # Fanning friction factor for viscous pressure drops
    if Re <= 0:
        return 0.0
    elif Re < 2000: # Laminar
        return 16.0 / Re
    elif Re < 20000: # Transition
        return 0.0791 * (Re**-0.25)
    else: # Turbulent
        return 0.046 * (Re**-0.2)

def dp_friction_fanning(f: float, rho: float, V: float, L: float, D: float) -> float: # Friction pressure drop 
    return 2.0 * f * rho * (V**2) * (L / D) 

def dp_evap(mdot: float, r_v: float, rho_v: float, mu_v: float, L_e: float) -> Tuple[float, float]:
    # Evaporator inertial pressure drop from Busse
    Re = 2 * mdot / PI / mu_v / r_v # Axial Reynolds number
    A_v = PI*r_v*r_v 
    Re_r = mdot/2/PI/L_e/mu_v # Radial Reynolds number
    psi = 0.61*Re_r + 0.61*Re_r / (3.6 + Re_r) # Velocity correction factor for mass addition effects
    phi = 16/Re*L_e/r_v/2 
    dpve = (mdot/A_v)**2 / rho_v * phi
    dpie = (mdot/A_v)**2 / rho_v * psi * phi
    return dpie, dpve

def dp_cond_inertial(fluid_name: str, mdot: float, r_v: float, rho_v: float, mu_v: float, lengths: SectionLengths) -> float:
    # Condenser inertial pressure drop calculations
    # For abs(Re_r) < 2.25 --> Busse
    # For abs(Re_r) > 2.25 --> Kemme
    # Properties based on condenser inlet 
    AV = PI*r_v*r_v
    Re_r = -1*mdot/2/PI/lengths.L_c/mu_v # Radial Reynolds number
    v = mdot/rho_v/(PI*r_v*r_v)
    if Re_r > -2.25:
        temp = 5+18/Re_r
        beta = 15/22*(temp + sqrt(temp**2 - 44/5))
        coeff = Re_r*(7/9 - 8/27*beta + 23/405*beta**2)
        dpic = coeff*4*mu_v*v*lengths.L_c/r_v/r_v
    else:
        temp = (2*lengths.L_e+4*lengths.L_a)/lengths.L_c
        recov = (Re_r + 2) / (1.23 * Re_r - temp)
        dpic = -1*recov*rho_v*v**2
    return dpic

def dp_hydrostatic(rho_l: float, L_t: float, theta_deg: float) -> float:
    # Hydrostatic pressure drop calculations
    return rho_l * GRAV * L_t * sin(theta_deg * PI / 180.0)

def dp_liquid(mdot: float, r_v:float, t_a:float, rho: float, mu: float, L_eff: float, A_l: float) -> float:
    # Liquid viscous pressure drop calculations
    # Will be adding different wick designs, currently annular wick only
    K = t_a**2 / 12
    dpl = 6 * mu * mdot * L_eff / (pi*rho*r_v*t_a**3)
    return dpl

def max_capillary_pressure(sigma: float, r_eff_cm: float) -> float:
    # Maximum capillary pressure provided by the wick
    return 0.0 if r_eff_cm <= 0 else 2.0 * sigma / r_eff_cm

# --- Fanno flow relations (compressible adiabatic with friction) ---
def speed_of_sound(gamma:float, rho: float, p:float) -> float:
    # Speed of sound calculation
    return sqrt(gamma * p / max(rho, 1e-9))

def fanno_F(M: float, gamma: float) -> float:
    # Fanno flow momentum equation with L_a as the choking length
    if M <= 0:
        return float("inf")
    term1 = (1.0 - M**2) / (gamma * M**2)
    term2 = ((gamma + 1.0) / (2.0 * gamma)) * log(((gamma + 1.0) * M**2) / (2.0 + (gamma - 1.0) * M**2))
    return term1 + term2

def pressr_Fanno(M: float, gamma: float) -> float:
    # This function returns the pressure ratio P/P^*
    from math import sqrt
    num = (gamma + 1.0) / 2.0
    den = 1.0 + (gamma - 1.0) * 0.5 * M**2
    return (1.0 / M) * sqrt(num / den)

def fanno_downstream_mach(M1: float, gamma: float, four_f_L_over_D: float) -> Tuple[bool, float]:
    # This function finds the Mach number at the adiabatic section outlet
    # from the Fanno flow relations
    Fmax = max(fanno_F(M1, gamma), 0.0)
    if four_f_L_over_D >= Fmax - 1e-9:
        return True, 1.0
    target = fanno_F(M1, gamma) - four_f_L_over_D
    a, b = 1e-6, 1.0
    for i in range(100):
        mid = 0.5 * (a + b)
        Fmid = fanno_F(mid, gamma)
        if Fmid > target:
            a = mid
        else:
            b = mid
        err = Fmid - target
        if abs(err) < 1e-4*target:
            # print("M_1 = ", M1, "M_2 = ", mid)
            break
        if i >= 100:
            raise RuntimeError(f"Maximum number of iterations (100) reached for Mach number calculation.")
    return False, 0.5 * (a + b)

def adiabatic_section_pressure_outlet(P_in: float, gamma: float, f: float, L_cm: float, D_cm: float, M_in: float) -> Tuple[float, bool, float]:
    # This function finds the pressure at the adiabatic section outlet from the Mach number
    four_f_L_over_D = 4.0 * f * (L_cm / D_cm)
    choked, M_out = fanno_downstream_mach(M_in, gamma, four_f_L_over_D)
    if choked:
        M_out = 1.0
    PR = pressr_Fanno(M_in, gamma) / pressr_Fanno(M_out, gamma) # Ratios of pressures at the adiabatic outlet over inlet (P2/P1)
    if (PR < 1)  or (PR > 2.08):
        PR = 2.08
    return P_in / PR, choked, M_out

# --- Sonic Limit
def sonic_limit(fluids: Dict[str, float], geom: Geometry, lengths: SectionLengths) -> float:
    # Finds sonic limit based on the Mach number at the adiabatic exit being 1
    mu_v = fluids["mu_v"]
    hfg = fluids["h_fg"]
    rhov = fluids["rho_v"]
    Pv = fluids["p_sat"]
    rk = fluids["gamma"]
    RV = geom.radius_cm-geom.t_a_cm-geom.t_w_cm
    AV = PI * RV**2
    qtest = 5000 # Initial guess
    mdot = qtest / hfg

    # Iterate to find the evaporator exit Mach number
    for _ in range(100):
        Re = 2 * mdot / (PI * RV * mu_v)
        f = fanning_friction_factor(Re)
        f_fanno = 4 * f * lengths.L_a / (2 * RV)
        target = f_fanno
        a, b = 1e-6, 1.0
        for _ in range(60):
            mid = 0.5 * (a + b)
            Fmid = fanno_F(mid, rk)
            if Fmid > target:
                a = mid
            else:
                b = mid
            err = Fmid - target
            if abs(err) < 1e-4*target:
                break
        rmach = 0.5 * (a + b)
        
        # Use evaporator exit Mach number to calculate the next limit guess
        Qsonic = rmach * AV * hfg * sqrt(rhov * Pv) 
        #Pchoke, choked, M_out = adiabatic_section_pressure_outlet(Pv, rk, f, lengths.L_a, DV, rmach)
        
        # Fixed point iteration for Qsonic
        dQ = abs(Qsonic - mdot * hfg)
        if dQ <= 1:
            break
        mdot = Qsonic / hfg

    return Qsonic

# --- Entrainment limits ---
def entrainment_velocity_weber(sigma: float, rho_v: float, wavelength_cm: float) -> float:
    # Vapor velocity resulting in We' = 1
    from math import pi, sqrt
    return sqrt(2.0 * pi * sigma / (rho_v * wavelength_cm))

def flooding_coefficient_Bo(gamma_bond: float) -> float:
    # Flooding coefficient for wickless
    from math import tanh
    return tanh(0.5 * gamma_bond ** 0.25)

def entrainment_limit_wickless_smooth(rho_v: float, rho_l: float, sigma: float, g: float = GRAV) -> float:
    # Flooding limit for wickless
    from math import sqrt
    Bo = g * (rho_l - rho_v) / sigma # Bond number
    Cw = flooding_coefficient_Bo(Bo)
    return Cw * sqrt(sigma / max(rho_v, 1e-9))

def entrainment_limit(fluid_name: str, fluids: Dict[str, float], geom: Geometry, lengths: SectionLengths) -> float:
    # Calculate entrainment limit based on the critical velocity occoring at the adiabatic exit
    RV = geom.radius_cm-geom.t_a_cm-geom.t_w_cm
    AV = PI * RV**2

    P2 = fluids["p_sat"] # initial guess for adiabatic exit pressure

    for _ in range(100):
        Tbc = tsat_from_p_cgs(Fluid(fluid_name), P2)
        fluids_cond_inlet = props_for(fluid_name, Tbc) # Adiabatic section exit fluid properties
    
        V_crit = sqrt(2.0 * pi * fluids_cond_inlet["sigma"] / (fluids_cond_inlet["rho_v"]*geom.wavelength_cm)) # type: ignore
        rmach = V_crit / speed_of_sound(fluids_cond_inlet["gamma"], fluids_cond_inlet["rho_v"], P2)
        if rmach > 1.0:
            rmach = 1.0
    
        Re = 2 * RV * fluids_cond_inlet["rho_v"] * V_crit / fluids_cond_inlet["mu_v"]
        f = fanning_friction_factor(Re)
        f_fanno = 4 * f * lengths.L_a / (2 * RV)
        fanno_F(rmach, fluids_cond_inlet["gamma"])
        choked, M_in = fanno_downstream_mach(rmach, fluids_cond_inlet["gamma"], -1*f_fanno) # find evaporator exit (upstream) Mach
        P2 = fluids["p_sat"] / pressr_Fanno(M_in, fluids["gamma"]) * pressr_Fanno(rmach, fluids_cond_inlet["gamma"])
        Tbc2 = tsat_from_p_cgs(Fluid(fluid_name), P2)
        # print("Tbc = ", Tbc2)
        if abs(Tbc-Tbc2) < 1:
            break
        Tbc = Tbc2
    Q_entrn = sqrt(2 * PI * fluids_cond_inlet["rho_v"] * fluids_cond_inlet["sigma"] / geom.wavelength_cm) * fluids_cond_inlet["h_fg"] * AV  # type: ignore
    return Q_entrn

def capillary_limit(fluid_name: str, fluids: Dict[str, float], geom: Geometry, lengths: SectionLengths, theta_deg: float) -> float:
    q_max_search_W = 1.0e5 # Max search power for bisection
    tol_rel = 1.0e-4 # Relative tolerance for convergence
    dp_cap_cgs = max_capillary_pressure(fluids["sigma"], geom.effective_pore_radius_cm)
    if dp_cap_cgs <= 0:
        q_cap = None
    else:
        def dp_error(q_watts: float) -> float:
            if q_watts <= 0:
                return -dp_cap_cgs
            m_dot_g_s = q_watts / max(fluids["h_fg"], 1e-9)
            pb = pressure_breakdown_cgs(fluid_name, lengths, geom, fluids, m_dot_g_s, theta_deg)
            return pb.dp_tot - pb.dp_capillary

        lo, hi = 0.0, q_max_search_W
        err_lo = dp_error(lo)
        err_hi = dp_error(hi)
        if err_hi <= 0:
            q_cap = hi
        else:
            for _ in range(100):
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
    return q_cap # type: ignore

@dataclass
class PressureBreakdown:
    dp_fric_evap: float
    dp_inert_evap: float
    dp_inert_cond: float
    dp_fric_cond: float
    dp_a: float
    dp_hydro: float
    dp_capillary: float
    choked_adiab: bool
    M_out_adiab: float
    dp_tot: float
    dp_liq: float
    dp_e: float
    dp_c: float
    T_bc: float
    T_be: float

def pressure_breakdown_cgs(
    fluid_name: str,
    lengths: SectionLengths,
    geom: Geometry,
    fluids: Dict[str, float],
    m_dot_g_s: float,
    theta_deg: float,
) -> PressureBreakdown:
    # This function calculates the pressure drops for a given heat pipe and mdot
    from math import sqrt
    # Calculate geometric parameters
    R = geom.radius_cm
    RV = geom.radius_cm-geom.t_a_cm-geom.t_w_cm
    A = PI * R * R
    DV = 2.0 * RV
    AV = PI * RV**2
    AL = A - PI*(R-geom.t_a_cm)**2 
    L_eff = lengths.L_e/2 + lengths.L_a + lengths.L_c/2
    L_t = lengths.L_e + lengths.L_a + lengths.L_c

    # Maximum capillary pressure 
    dp_cap = max_capillary_pressure(fluids["sigma"], geom.effective_pore_radius_cm) 

    V_exit = m_dot_g_s / (fluids["rho_v"] * AV) # Evaporator exit velocity

    # Evaporator viscous and inertial pressure drops
    Re_e = reynolds(m_dot_g_s, fluids["mu_v"], fluids["rho_v"], AV, DV)
    f_e = fanning_friction_factor(Re_e)
    dpei, dp_fric_evap = dp_evap(m_dot_g_s, RV, fluids["rho_v"], fluids["mu_v"], lengths.L_e)
    # dp_fric_evap = dp_friction_fanning(f_e, fluids["rho_v"], V_exit, lengths.L_e/2, DV)
    dp_e = dp_fric_evap + dpei
    Tbe = tsat_from_p_cgs(Fluid(fluid_name), fluids["p_sat"]+dp_e)

    # Adiabatic viscous pressure drop
    dp_adiab = 0.0
    choked = False
    M_out = 0.0
    if lengths.L_a > 0.0:
        # a_in = speed_of_sound(fluids["gamma"], fluids["rho_v"], fluids["p_sat"]) # Use this normally
        a_in = sqrt(RBAR*fluids["T_sat"]/(fluids["m_w"]/1000))*100 # This is not correct, it is missing the gamma but this is what HTPIPE has, use this for validation only 
        M_in = V_exit / max(a_in, 1e-9)
        Re_a = reynolds(m_dot_g_s, fluids["mu_v"], fluids["rho_v"], AV, DV)
        f_a = fanning_friction_factor(Re_a)
        if M_in > 0.2:
            P2, choked, M_out = adiabatic_section_pressure_outlet(fluids["p_sat"], fluids["gamma"], f_a, lengths.L_a, DV, M_in)
        if M_out > 0.3:
            dp_adiab = max(fluids["p_sat"] - P2, 0.0)/2 # HTPIPE has an 1/2 factor here but not sure why
            # print("M_in = ", M_in, "M_out = ", M_out)
        else:
            dp_adiab = dp_friction_fanning(f_a, fluids["rho_v"], V_exit, lengths.L_a, DV)
            P2 = fluids["p_sat"] - dp_adiab
        Tbc = tsat_from_p_cgs(Fluid(fluid_name), P2)
        # print("Tbc = ", Tbc)
        fluids_cond_inlet = props_for(fluid_name, Tbc) # Adiabatic section exit fluid properties
    else:
        fluids_cond_inlet = fluids # If there is no adiabatic section use the evaporator exit fluid properties 
    
    # Condenser viscous and inertial pressure drops
    Re_c = reynolds(m_dot_g_s, fluids_cond_inlet["mu_v"], fluids_cond_inlet["rho_v"], AV, DV)
    f_c = fanning_friction_factor(Re_c)
    dp_fric_cond = dp_friction_fanning(f_c, fluids_cond_inlet["rho_v"], V_exit, lengths.L_c/2, DV)
    dpci = dp_cond_inertial(fluid_name, m_dot_g_s, RV, fluids_cond_inlet["rho_v"], fluids_cond_inlet["mu_v"], lengths)
    dp_c = dp_fric_cond + dpci

    # Liquid viscous pressure drop
    dp_liq = dp_liquid(m_dot_g_s, RV, geom.t_a_cm, fluids["rho_l"], fluids["mu_l"], L_eff, AL)
    # dp_liq_e = dp_liq * (lengths.L_e/2/L_eff)
    # dp_liq_a = dp_liq * (lengths.L_a/L_eff)
    dp_liq_c = dp_liq * (lengths.L_e/2/L_eff)

    # Hydrostatic pressure drop
    dp_hyd = dp_hydrostatic(fluids["rho_l"], L_t, theta_deg)

    # Calculate total pressure drops based on the wet point location
    if dp_hyd < 0: # if gravity-assist, wet point at evaporator endcap
        dp_tot = dpei+dpci+dp_fric_evap+dp_adiab+dp_fric_cond+dp_liq+dp_hyd
    else: # if horizontal or against gravity
        dplcz = dp_liq_c + dp_hyd*lengths.L_c/L_t
        dplvc = dp_c + dplcz
        if dplvc < 0: # if there is net pressure gain in the condenser, wet point at condenser inlet
            dp_hyd = dp_hyd * (lengths.L_a + lengths.L_e)/L_t
            dp_liq = dp_liq - dp_liq_c
            dp_c = 0
            dp_tot = dp_e + dp_adiab + dp_liq + dp_hyd
        else: # else wet point is at condenser endcap
            dp_tot = dpei+dpci+dp_fric_evap+dp_adiab+dp_fric_cond+dp_liq+dp_hyd

    return PressureBreakdown(
        dp_fric_evap=dp_fric_evap,
        dp_inert_evap=dpei,
        dp_inert_cond=dpci,
        dp_fric_cond=dp_fric_cond,
        dp_a=dp_adiab,
        dp_hydro=dp_hyd,
        dp_capillary=dp_cap,
        choked_adiab=choked,
        M_out_adiab=M_out,
        dp_liq = dp_liq,
        dp_e = dp_e,
        dp_c = dp_c,
        dp_tot = dp_tot,
        T_bc = Tbc,
        T_be = Tbe,
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
    fluid_name: str,
    theta_deg: float,
    ) -> Limits:
    # This function finds the power limits at the given geometry and fluid properties
    from math import sqrt

    L_eff = lengths.L_e/2 + lengths.L_a + lengths.L_c/2
    AV = PI * (geom.radius_cm-geom.t_a_cm-geom.t_w_cm)**2

    # --- Sonic Limit ---
    if fluids["rho_v"] <= 0 or fluids["p_sat"] <= 0:
        q_sonic = None
    else:
        q_sonic = sonic_limit(fluids, geom, lengths)

    # --- Viscous Limit ---
    q_visc = AV**2*fluids["h_fg"]*fluids["rho_v"]*fluids["p_sat"]/16/PI/fluids["mu_v"]/L_eff

    # --- Capillary Limit ---
    q_cap = capillary_limit(fluid_name, fluids, geom, lengths, theta_deg)


    # --- Entrainment Limit ---
    if geom.wavelength_cm:
        if geom.wavelength_cm <= 0 or fluids["rho_v"] <= 0:
            q_ent = None
        else:
            q_ent = entrainment_limit(fluid_name, fluids, geom, lengths)
    else:
        q_ent = None

    return Limits(q_capillary_W=q_cap, q_sonic_W=q_sonic, q_entrainment_W=q_ent, q_viscous_W=q_visc)
