
from math import exp, log, log10
import numpy as np
from typing import cast
# from numpy.polynomial import polynomial
from dataclasses import dataclass
from enum import Enum

# ---- Constants (cgs) ----
R_BAR_CGS = 8.314e7  # erg/(mol·K)
TORR_TO_DYNCM2 = 1333.0  # 1 torr = 1333 dyn/cm^2
DYNCM2_TO_PA = 0.1  # 1 dyn/cm^2 = 0.1 Pa
DYN_PER_CM_TO_N_PER_M = 1.0e-3  # 1 dyn/cm = 0.001 N/m
POISE_TO_PA_S = 0.1  # 1 P = 0.1 Pa·s
GCML3_TO_KGM3 = 1000.0  # 1 g/cm^3 = 1000 kg/m^3

class Fluid(Enum):
    LITHIUM = "lithium"
    SODIUM = "sodium"
    POTASSIUM = "potassium"
    MERCURY = "mercury"
    WATER = "water"

@dataclass
class PropertiesCGS:
    pv: float          # saturation vapor pressure [dyn/cm^2]
    mw: float          # molecular weight [g/mol]
    rhol: float        # liquid density [g/cm^3]
    muv: float         # vapor dynamic viscosity [poise]
    mul: float         # liquid dynamic viscosity [poise]
    hfg: float         # latent heat of vaporization [kJ/kg]
    sigma: float       # surface tension [dyn/cm]
    gamma: float       # ratio of specific heats [-]
    rhov: float        # vapor density via ideal gas at saturation [g/cm^3]

def _pv_cgs(fluid: Fluid, T: float) -> float:
    """Saturation vapor pressure [dyn/cm^2] at temperature T [K]."""
    if fluid == Fluid.LITHIUM:
        return 10 ** (7.67 - 7740.0 / T) * TORR_TO_DYNCM2
    elif fluid == Fluid.SODIUM:
        return 3.83e10 / exp(12160.0 / T)
    elif fluid == Fluid.POTASSIUM:
        # return 2.197e10 / exp(10223.0 / T)   
        result = np.polyval([1.3299e-14, -8.6184e-11, 2.231e-07, -0.00029381, 0.20636, -51.574], T)
        return exp(result)*10
        # return val
    elif fluid == Fluid.MERCURY:
        return 1.332e3 * exp(17.85 - 7059.5 / T)
    elif fluid == Fluid.WATER:
        return 3.975e11 * exp(-4872.0 / T)
    else:
        raise ValueError("Unknown fluid")

def _mw(fluid: Fluid) -> float:
    return {
        Fluid.LITHIUM: 6.94,
        Fluid.SODIUM: 23.0,
        Fluid.POTASSIUM: 39.0983,
        Fluid.MERCURY: 200.59,
        Fluid.WATER: 18.0,
    }[fluid]

def _rhol_cgs(fluid: Fluid, T: float) -> float:
    if fluid == Fluid.LITHIUM:
        return 0.555 - 0.934e-4 * T
    elif fluid == Fluid.SODIUM:
        return 1.018 - 2.34e-4 * T
    elif fluid == Fluid.POTASSIUM:
        # return 0.909 - 2.41e-4 * T
        result = np.polyval([-3.4521e-14, 1.9884e-10, -4.2981e-07, 0.00042662, -0.42954, 940.43], T)
        return float(result)/1000
    elif fluid == Fluid.MERCURY:
        return 12.75 - 2.50e-3 * T
    elif fluid == Fluid.WATER:
        return 1.49 - 1.40e-3 * T
    else:
        raise ValueError("Unknown fluid")

def _muv_cgs(fluid: Fluid, T: float) -> float:
    if fluid == Fluid.LITHIUM:
        return 1.2e-7 * T - 6.0e-6
    elif fluid == Fluid.SODIUM:
        #return 1.6e-7 * T - 5.0e-6
        return 1.131e-7*T + 9.414e-5
    elif fluid == Fluid.POTASSIUM:
        # return 1.46e-7 * T - 5.0e-6
        return float(np.polyval([-3.5469e-21, 1.9365e-17, -3.9354e-14, 3.0725e-11, 6.9661e-09, 5.4558e-06], T))/0.1
    elif fluid == Fluid.MERCURY:
        return 1.033e-6 * T - 2.0e-5
    elif fluid == Fluid.WATER:
        return 6.91e-5 * exp(4.67e-6 * T * T)
    else:
        raise ValueError("Unknown fluid")

def _mul_cgs(fluid: Fluid, T: float) -> float:
    if fluid == Fluid.LITHIUM:
        return 0.00142 * exp(5.48e10 / (R_BAR_CGS * T))
    elif fluid == Fluid.SODIUM:
        # The OCR shows "10**(-3.0494 + 30.9/T)"
        #return 10 ** (-3.0494 + 30.9 / T)
        return exp(-6.4406 - 0.39580*log(T) + 556.8/T)/0.1
    elif fluid == Fluid.POTASSIUM:
        # return 0.75 * 10 ** (-2.9995 + 245.0 / T)
        return float(np.polyval([-4.5263e-19, 2.8447e-15, -7.0683e-12, 8.7656e-09, -5.5722e-06, 0.0016193], T))/0.1
    elif fluid == Fluid.MERCURY:
        return 5.138e-3 * exp(364.3 / T)
    elif fluid == Fluid.WATER:
        return 6.22e-5 * exp(1.478e3 / T)
    else:
        raise ValueError("Unknown fluid")

def _hfg_kj_per_kg(fluid: Fluid, T: float) -> float:
    if fluid == Fluid.LITHIUM:
        # HFG = 0.2412E5 + T*(-0.0952 + T*(-0.2282E-2 + T*0.6261E-6))
        return (0.2412e5 + T * (-0.0952 + T * (-0.2282e-2 + T * 0.6261e-6)))
    elif fluid == Fluid.SODIUM:
        # 5.226E3 + T*(-1.474 + T*(3.292E-4 - 5.462E-8*T))
        #return 5.226e3 + T * (-1.474 + T * (3.292e-4 - 5.462e-8 * T))
        return -1*0.9311*T + 4.973e3
    elif fluid == Fluid.POTASSIUM:
        # 2.92E3 + T*(-1.104 - T*(0.001323 - 4.123E-7*T))
        #return 2.92e3 + T * (-1.104 - T * (0.001323 - 4.123e-7 * T))
        # return -1*0.511*T + 2.449e3
        return float(np.polyval([-2.0228e-11, 4.2749e-08, 0.00015934, -0.60104, 135.36, 2226400], T))/1000
    elif fluid == Fluid.MERCURY:
        return 355.0 * exp(-3.45e-4 * T)
    elif fluid == Fluid.WATER:
        return 3800.0 - 4.333 * T
    else:
        raise ValueError("Unknown fluid")

def _sigma_cgs(fluid: Fluid, T: float) -> float:
    if fluid == Fluid.LITHIUM:
        return 453.0 - 0.148 * T
    elif fluid == Fluid.SODIUM:
        return 220.0 - 0.091 * T
    elif fluid == Fluid.POTASSIUM:
        # return 136.0 - 0.0645 * T
        return float(np.polyval([-5.0487e-27, 1.7371e-23, -5.7598e-20, 2.7756e-17, -6.6e-05, 0.13127], T))*1000 
        # return float(np.polyval([12.12e-9, -32.65e-6, -45.72e-3, 130.48], T))# Russian book
    elif fluid == Fluid.MERCURY:
        return 562.4 - 0.308 * T
    elif fluid == Fluid.WATER:
        return 133.5 - 0.205 * T
    else:
        raise ValueError("Unknown fluid")

def _gamma(fluid: Fluid, T: float) -> float:
    if fluid == Fluid.LITHIUM:
        return 1.7997 - 0.0001479 * T
    elif fluid == Fluid.SODIUM:
        # Not legible in the OCR; assume monatomic ideal-gas value (~5/3).
        # return 1.667
        return 1.7590 - 0.0001720 * T
    elif fluid == Fluid.POTASSIUM:
        return 1.7402 - 0.0001230 * T
    elif fluid == Fluid.MERCURY:
        return 1.667
    elif fluid == Fluid.WATER:
        return 1.324
    else:
        raise ValueError("Unknown fluid")
    
def _rho_v(fluid: Fluid, T: float) -> float: # Future development, use for HTPIPE validation with potassium
    mw = _mw(fluid)
    if fluid == Fluid.LITHIUM:
        return (mw * _pv_cgs(fluid, T) / (R_BAR_CGS * T))*GCML3_TO_KGM3
    elif fluid == Fluid.SODIUM:
        # Not legible in the OCR; assume monatomic ideal-gas value (~5/3).
        # return 1.667
        return (mw * _pv_cgs(fluid, T) / (R_BAR_CGS * T))*GCML3_TO_KGM3
    elif fluid == Fluid.POTASSIUM:
        return exp(np.polyval([1.268e-14, -8.2191e-11, 2.1319e-07, -0.00028146, 0.19765, -60.872],T))/1000
    elif fluid == Fluid.MERCURY:
        return (mw * _pv_cgs(fluid, T) / (R_BAR_CGS * T))*GCML3_TO_KGM3
    elif fluid == Fluid.WATER:
        return (mw * _pv_cgs(fluid, T) / (R_BAR_CGS * T))*GCML3_TO_KGM3
    else:
        raise ValueError("Unknown fluid")

def tsat_from_p_cgs(fluid: Fluid, P: float) -> float:
    """Saturation temperature [K] from pressure P [dyn/cm^2]."""
    if fluid == Fluid.LITHIUM:
        return 7740.0 / (7.67 - log10(P / TORR_TO_DYNCM2))
    elif fluid == Fluid.SODIUM:
        return 12180.0 / log(3.83e10 / P)
    elif fluid == Fluid.POTASSIUM:
        # return 10223.0 / log(2.197e10 / P)
        # return exp(np.polyval([0.032364, -0.127580, -0.112436, 1.112872, -1.557274, 8.915520],T*1e-3))/1000
        return 62.71*(P*0.1)**(0.1923)+459.001 # From fit of HTPIPE properties
    elif fluid == Fluid.MERCURY:
        return 7059.5 / (17.85 + log(1.332e3 / P))
    elif fluid == Fluid.WATER:
        return 4872.0 / log(3.975e11 / P)
    else:
        raise ValueError("Unknown fluid")

def properties_cgs(fluid: Fluid, T: float):
    """Return a bundle of properties at temperature T [K] in cgs units."""
    pv = _pv_cgs(fluid, T)
    mw = _mw(fluid)
    rhol = _rhol_cgs(fluid, T)
    muv = _muv_cgs(fluid, T)
    mul = _mul_cgs(fluid, T)
    hfg = _hfg_kj_per_kg(fluid, T)
    sigma = _sigma_cgs(fluid, T)
    gamma = _gamma(fluid, T)
    rhov = mw * pv / (R_BAR_CGS * T)  # ideal gas, g/cm^3, use otherwise
    if fluid == Fluid.POTASSIUM:
        rhov = _rho_v(fluid, T) # Use for HTPIPE validation with potassium
    return PropertiesCGS(pv, mw, rhol, muv, mul, hfg, sigma, gamma, rhov)

# SI wrapper if needed (not used in this demo notebook)
@dataclass
class PropertiesSI:
    pv: float       # Pa
    mw: float          # g/mol (kept for convenience)
    rhol: float       # kg/m^3
    muv: float        # Pa·s
    mul: float        # Pa·s
    hfg: float        # J/kg
    sigma: float       # N/m
    gamma: float       # [-]
    rhov: float       # kg/m^3

def properties_SI(fluid: Fluid, T: float):
    pv = _pv_cgs(fluid, T)*POISE_TO_PA_S
    mw = _mw(fluid)
    rhol = _rhol_cgs(fluid, T)*GCML3_TO_KGM3    
    muv = _muv_cgs(fluid, T)*POISE_TO_PA_S
    mul = _mul_cgs(fluid, T)*POISE_TO_PA_S
    hfg = _hfg_kj_per_kg(fluid, T)*1000
    sigma = _sigma_cgs(fluid, T)*DYN_PER_CM_TO_N_PER_M
    gamma = _gamma(fluid, T)
    rhov = (mw * _pv_cgs(fluid, T) / (R_BAR_CGS * T))*GCML3_TO_KGM3
    return PropertiesSI(pv, mw, rhol, muv, mul, hfg, sigma, gamma, rhov)
