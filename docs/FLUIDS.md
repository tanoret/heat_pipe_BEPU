# Module: `heatpipe.fluids`

## Purpose
Property correlations and helpers for five working fluids:
`lithium`, `sodium`, `potassium`, `mercury`, `water`.

The module uses **cgs** units internally to match legacy heatpipe practice.

## Constants (cgs)
- $\bar R = 8.314\times 10^{7}\ \text{erg/(mol·K)}$
- $1\ \text{torr} = 1333\ \text{dyn/cm}^2$

## Public API

### `class Fluid(Enum)`
`LITHIUM`, `SODIUM`, `POTASSIUM`, `MERCURY`, `WATER`

### `@dataclass PropertiesCGS`
All fields in **cgs** unless noted.
- `pv`: saturation pressure, dyn/cm²  
- `mw`: molecular weight, g/mol  
- `rhol`: liquid density, g/cm³  
- `muv`: vapor dynamic viscosity, P  
- `mul`: liquid dynamic viscosity, P  
- `hfg`: latent heat, **kJ/kg**  
- `sigma`: surface tension, dyn/cm  
- `gamma`: ratio of specific heats, –  
- `rhov`: vapor density at saturation (ideal gas), g/cm³

### `properties_cgs(fluid: Fluid, T: float) -> PropertiesCGS`
Returns a complete property bundle at temperature $T$ [K].

### `tsat_from_p_cgs(fluid: Fluid, P: float) -> float`
Returns saturation temperature $T$ [K] from pressure $P$ [dyn/cm²].

---

## Correlations

Let $T$ be Kelvin.

### Saturation pressure $p_v(T)$
- **Li**: $p_v = 10^{\,7.67 - \frac{7740}{T}}\times 1333$
- **Na**: $p_v = 3.83\times 10^{10}\,\exp\!\left(-\frac{12160}{T}\right)$
- **K**: $p_v = 2.197\times 10^{10}\,\exp\!\left(-\frac{10223}{T}\right)$
- **Hg**: $p_v = 1.332\times 10^3 \exp\!\left(17.85 - \frac{7059.5}{T}\right)$
- **H₂O**: $p_v = 3.975\times 10^{11}\,\exp\!\left(-\frac{4872}{T}\right)$

### Saturation temperature $T(p)$ (inverse forms)
- **Li**: $T = \dfrac{7740}{7.67 - \log_{10}(p/1333)}$
- **Na**: $T = \dfrac{12180}{\ln\!\left(3.83\times 10^{10}/p\right)}$
- **K**: $T = \dfrac{10223}{\ln\!\left(2.197\times 10^{10}/p\right)}$
- **Hg**: $T = \dfrac{7059.5}{17.85 + \ln\!\left(1.332\times 10^3/p\right)}$
- **H₂O**: $T = \dfrac{4872}{\ln\!\left(3.975\times 10^{11}/p\right)}$

### Liquid density $\rho_l(T)$ [g/cm³]
- **Li**: $\rho_l = 0.555 - 0.934\times 10^{-4} T$
- **Na**: $\rho_l = 1.018 - 2.34\times 10^{-4} T$
- **K**: $\rho_l = 0.909 - 2.41\times 10^{-4} T$
- **Hg**: $\rho_l = 12.75 - 2.50\times 10^{-3} T$
- **H₂O**: $\rho_l = 1.49 - 1.40\times 10^{-3} T$

### Vapor viscosity $\mu_v(T)$ [P]
- **Li**: $\mu_v = 1.2\times 10^{-7}\,T - 6.0\times 10^{-6}$
- **Na**: $\mu_v = 1.6\times 10^{-7}\,T - 5.0\times 10^{-6}$
- **K**: $\mu_v = 1.46\times 10^{-7}\,T - 5.0\times 10^{-6}$
- **Hg**: $\mu_v = 1.033\times 10^{-6}\,T - 2.0\times 10^{-5}$
- **H₂O**: $\mu_v = 6.91\times 10^{-5}\,\exp\!\left(4.67\times 10^{-6} T^2\right)$

### Liquid viscosity $\mu_l(T)$ [P]
- **Li**: $\mu_l = 1.42\times 10^{-3}\,\exp\!\left(\dfrac{5.48\times 10^{10}}{\bar R\,T}\right)$
- **Na**: $\mu_l = 10^{\,-3.0494 + 30.9/T}$
- **K**: $\mu_l = 0.75\times 10^{\,-2.9995 + 245/T}$
- **Hg**: $\mu_l = 5.138\times 10^{-3}\,\exp\!\left(\dfrac{364.3}{T}\right)$
- **H₂O**: $\mu_l = 6.22\times 10^{-5}\,\exp\!\left(\dfrac{1.478\times 10^{3}}{T}\right)$

### Latent heat $h_{fg}(T)$ [kJ/kg]
- **Li**: $h_{fg} = 0.2412\times 10^{5} + T\left(-0.0952 + T\left(-0.2282\times 10^{-2} + 0.6261\times 10^{-6} T\right)\right)$
- **Na**: $h_{fg} = 5.226\times 10^{3} + T\left(-1.474 + T\left(3.292\times 10^{-4} - 5.462\times 10^{-8} T\right)\right)$
- **K**: $h_{fg} = 2.92\times 10^{3} + T\left(-1.104 - T\left(1.323\times 10^{-3} - 4.123\times 10^{-7} T\right)\right)$
- **Hg**: $h_{fg} = 355.0\,\exp\!\left(-3.45\times 10^{-4} T\right)$
- **H₂O**: $h_{fg} = 3800.0 - 4.333\,T$

### Surface tension $\sigma(T)$ [dyn/cm]
- **Li**: $\sigma = 453.0 - 0.148\,T$
- **Na**: $\sigma = 220.0 - 0.091\,T$
- **K**: $\sigma = 136.0 - 0.0645\,T$
- **Hg**: $\sigma = 562.4 - 0.308\,T$
- **H₂O**: $\sigma = 133.5 - 0.205\,T$

### Specific heat ratio $\gamma(T)$ [$-$]
- **Li**: $\gamma = 1.7997 - 1.479\times 10^{-4} T$
- **Na**: $\gamma \approx 1.667$  *(placeholder in code)*
- **K**: $\gamma = 1.7402 - 1.230\times 10^{-4} T$
- **Hg**: $\gamma \approx 1.667$
- **H₂O**: $\gamma \approx 1.324$

### Vapor density at saturation (ideal gas)
$$
\rho_v = \frac{MW \cdot p_v}{\bar R\, T}\quad \text{[g/cm}^3\text{]}
$$

---

## Examples

```python
from heatpipe.fluids import Fluid, properties_cgs, tsat_from_p_cgs

P = properties_cgs(Fluid.SODIUM, 900.0)
T = tsat_from_p_cgs(Fluid.SODIUM, P.pv)   # ≈ 900 K

print(P.pv, P.rhol, P.muv, P.mul, P.hfg, P.sigma, P.gamma, P.rhov)
```

> **Notes**
> - Validity: correlations are intended for typical heat-pipe temperature ranges used in heatpipe. The module does not enforce bounds; avoid extreme extrapolation.
> - Units: all inputs/outputs are **cgs** except `hfg` (kJ/kg).
