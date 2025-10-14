# Module: `heatpipe.models`

Hydrodynamic model primitives: Reynolds/friction, pressure components, compressible adiabatic (Fanno) section, entrainment criteria, and convenience wrappers to compute a pressure budget and first-cut performance limits. All physics uses **cgs**.

## Data types

```python
from enum import Enum
from dataclasses import dataclass
from typing import Optional

class Regime(Enum):
    LAMINAR = "laminar"
    TURBULENT = "turbulent"

@dataclass
class SectionLengths:   # [cm]
    L_e: float          # evaporator length
    L_a: float          # adiabatic length
    L_c: float          # condenser length

@dataclass
class Geometry:
    radius_cm: float
    effective_pore_radius_cm: float = 0.0     # capillary characteristic radius
    wavelength_cm: Optional[float] = None     # interfacial wavelength (entrainment)
    passages: int = 1

@dataclass
class FlowFlags:
    vapor_regime_evap: Regime = Regime.LAMINAR
    vapor_regime_adiab: Regime = Regime.LAMINAR
    vapor_regime_cond: Regime = Regime.LAMINAR

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

@dataclass
class Limits:
    q_capillary_W: Optional[float]
    q_sonic_W: Optional[float]
    q_entrainment_W: Optional[float]
```

`GRAV = 980.0` cm/s² is used for hydrostatics. Circular tube: $D_h = 2r$.

---

## Utilities

### `cross_section_area(radius_cm) -> cm²`
$A = \pi r^2$

### `reynolds(m_dot, \mu, \rho, A, D) -> Re`
Flow speed $V = \dot m / (\rho A)$, so
$$
\mathrm{Re} = \frac{\rho V D}{\mu}
$$

### `fanning_friction_factor(Re, regime) -> f`
- Laminar: $f = 16/\mathrm{Re}$
- Turbulent: $f = 0.0791\,\mathrm{Re}^{-1/4}$

### `velocity_profile_factor(regime) -> A_\text{prof}`
$A = 1.234$ (laminar), $A = 2.22$ (turbulent)

---

## Pressure terms (all $\Delta p$ in dyn/cm²)

### `dp_friction_fanning(f, \rho, V, L, D)`
$$
\Delta p_f = 2\,f\,\rho\,V^2\,\frac{L}{D}
$$

### `dp_inertial(A, \rho_v, V_\text{exit})`
$$
\Delta p_i = A\,\rho_v\,V_{\text{exit}}^2
$$

### `dp_hydrostatic(\rho_l, \Delta z, \theta_\deg)`
$$
\Delta p_z = \rho_l\,g\,\Delta z\,\sin\theta
$$

### `capillary_pressure_drop(\sigma, R_\text{eff})`
$$
\Delta p_\text{cap} = \frac{2\,\sigma}{R_\text{eff}}
$$

---

## Compressible adiabatic section (Fanno model)

**Goal:** given inlet state $(P_1, M_1, \gamma)$ and a duct with constant area, wall friction $f$, length $L$, diameter $D$, compute if choking occurs and the outlet pressure $P_2$.

### Fanno function $F(M)$
$$
F(M) = \frac{1 - M^2}{\gamma M^2}
      + \frac{\gamma + 1}{2\gamma}\,
        \ln\!\left(\frac{(\gamma+1)M^2}{2 + (\gamma-1)M^2}\right)
$$

### Distance to sonic
Maximum $4fL/D$ to reach the sonic condition from $M_1$:
$$
\left(\frac{4fL}{D}\right)_\max = F(M_1)
$$

### Pressure ratio function (heatpipe `PRESSR` analog)
$$
\frac{P}{P^*} = \frac{1}{M}
\sqrt{\frac{(\gamma+1)/2}{1 + (\gamma-1)M^2/2}}
$$

### Implementation
- `fanno_downstream_mach(M1, γ, four_f_L_over_D)`  
  Returns `(choked, M2)` on the subsonic branch by solving $F(M_2) = F(M_1) - 4fL/D$ (or `choked=True` if the length exceeds the distance to sonic).
- `adiabatic_section_pressure_outlet(P_in, γ, f, L, D, M_in)`  
  Uses `PRESSR(M)` to compute $P_2$: $P_1/P_2 = \text{PRESSR}(M_1)/\text{PRESSR}(M_2)$.

---

## Entrainment limits

### Screened/wicked (Weber criterion)
Critical vapor velocity when interfacial waves entrain liquid:
$$
V_\text{crit} = \sqrt{\frac{2\pi\,\sigma}{\rho_v\,\lambda}}
$$
Mass flux $G = \rho_v V_\text{crit}$; a power limit follows as $Q \approx (G\,h_{fg})A$ (convert g→kg).

Function: `entrainment_velocity_weber(σ, ρ_v, λ)`

### Smooth, wickless (flooding correlation)
Bond number $Bo = g(\rho_l-\rho_v)/\sigma$, empirical
$$
C_w = \tanh\!\left(0.5\,Bo^{1/4}\right),
\qquad
V_\text{scale} \sim C_w\sqrt{\frac{\sigma}{\rho_v}}
$$

Function: `entrainment_limit_wickless_smooth(ρ_v, ρ_l, σ, g=GRAV) -> V_scale`

> These are **first-cut** estimates used for sweeps/plots. Detailed wick/artery liquid-side losses and wet-point location will further constrain the capillary limit in full designs.

---

## Orchestration helpers

### `pressure_breakdown_cgs(lengths, geom, fluids, flags, m_dot_g_s, theta_deg=0.0) -> PressureBreakdown`
Inputs:
- `lengths`: `SectionLengths` (cm)
- `geom`: `Geometry` ($r$, $R_\text{eff}$, optional $\lambda$)
- `fluids`: dict with keys  
  `p_sat, rho_v, rho_l, mu_v, mu_l, sigma, gamma, h_fg` (from `heatpipe.io.props_for`)
- `flags`: `FlowFlags` (laminar/turbulent selections)
- `m_dot_g_s`: mass flow, g/s
- `theta_deg`: inclination angle, degrees

Computes:
- $V$ from $\dot m, \rho_v, A$  
- frictional drops in evaporator & condenser  
- inertial term at the evaporator exit  
- hydrostatic head over total axial length and $\theta$  
- adiabatic compressible drop via Fanno (and whether it **chokes**)  
- available capillary head $2\sigma/R_\text{eff}$

Returns all components (dyn/cm²), plus `choked_adiab` and `M_out_adiab`.

### `estimate_limits_cgs(lengths, geom, fluids, flags) -> Limits`
First-cut performance limits (W) using area $A=\pi r^2$ and $h_{fg}$:

- **Sonic**: speed of sound $a=\sqrt{\gamma p/\rho}$ at the evaporator exit; mass flux $G=\rho_v a$, then $Q \approx G\,h_{fg}\,A$ (g→kg).
- **Capillary**: use available $\Delta p_{\text{cap}}=2\sigma/R_\text{eff}$ and a nominal friction balance to get a characteristic velocity
  $V \approx \sqrt{\Delta p\,D / (2 f_\text{nom}\,\rho_v L_\text{char})}$ with $f_\text{nom}\approx 0.02$, then $Q \approx \rho_v V h_{fg} A$.
- **Entrainment**: if `wavelength_cm` is provided, use the Weber criterion to compute $V_\text{crit}$ and convert to $Q$.

> **Note**: these estimators are intentionally conservative and meant for parametric studies/plots.

---

## Example

```python
from heatpipe.io import props_for
from heatpipe.models import *

# State / fluid (cgs dict)
T = 850.0
f = props_for("sodium", T)

# Geometry & lengths
geom = Geometry(radius_cm=0.5, effective_pore_radius_cm=0.02, wavelength_cm=0.05)
L = SectionLengths(L_e=30.0, L_a=20.0, L_c=30.0)
flags = FlowFlags(vapor_regime_evap=Regime.LAMINAR,
                  vapor_regime_adiab=Regime.TURBULENT,
                  vapor_regime_cond=Regime.LAMINAR)

# Pressure budget at trial mass flow
pb = pressure_breakdown_cgs(L, geom, f, flags, m_dot_g_s=2.0, theta_deg=0.0)

# First-cut limits
lims = estimate_limits_cgs(L, geom, f, flags)

print(pb)
print(lims)
```

---

## Limitations & next steps
- The capillary limit currently uses a single-path viscous balance with a nominal $f$; full wick/artery **porous** models and wet-point location will refine this.
- Fanno treatment assumes constant area and properties across the adiabatic section (consistent with heatpipe’s 1-D approach).
