
# heat_pipe_BEPU — Hydrodynamic Heat Pipe Modeling with BEPU

heatpipe is a unit-consistent Python package for **hydrodynamic heat-pipe** analysis inspired by heatpipe.
It provides:
- Core models for pressure components, Fanno (compressible adiabatic) section, and limit estimators (capillary, sonic, entrainment).
- A **UQ/BEPU layer** with Monte Carlo + **Wilks**, **PCE** + Sobol indices, and a simple **GP** surrogate.
- Plotting helpers and clean JSON/YAML I/O.
- Markdown documentation and a demo notebook.

> Units: internal physics uses **cgs** (dyn, cm, g, s) to match legacy heatpipe.

## Installation
```bash
pip install -e .
# Optional: YAML config support
pip install 'heatpipe[yaml]'
```

## Quick start
```python
from heatpipe.io import props_for
from heatpipe.models import (Regime, SectionLengths, Geometry, FlowFlags,
                           pressure_breakdown_cgs, estimate_limits_cgs)
T = 900.0
f = props_for("sodium", T)
geom = Geometry(radius_cm=0.5, effective_pore_radius_cm=0.02, wavelength_cm=0.05)
L = SectionLengths(30.0, 20.0, 30.0)
flags = FlowFlags(Regime.LAMINAR, Regime.TURBULENT, Regime.LAMINAR)
pb = pressure_breakdown_cgs(L, geom, f, flags, m_dot_g_s=2.0, theta_deg=0.0)
lims = estimate_limits_cgs(L, geom, f, flags)
print(pb); print(lims)
```

### Uncertainty (BEPU)
```python
from heatpipe.uq import UncertainInput, UQSpec, uq_monte_carlo
spec = UQSpec(inputs=[
    UncertainInput("radius_cm", "normal", {"mean":0.5, "std":0.02}, lower=0.3, upper=0.8),
    UncertainInput("effective_pore_radius_cm", "uniform", {"low":0.01, "high":0.04}),
    UncertainInput("wavelength_cm", "uniform", {"low":0.04, "high":0.06}),
    UncertainInput("L_e", "deterministic", {"value":30.0}),
    UncertainInput("L_a", "deterministic", {"value":20.0}),
    UncertainInput("L_c", "deterministic", {"value":30.0}),
    UncertainInput("T_K", "normal", {"mean":850.0, "std":30.0})
])
X, outs = uq_monte_carlo(spec, n=300, seed=7)
print(outs["Q_lim_W"].min(), outs["Q_lim_W"].mean())
```

## Documentation
See `documentation/` for FLUIDS, MODELS, UQ, and supporting guides.
