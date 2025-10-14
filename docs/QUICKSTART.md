# Quickstart

## Setup
```bash
python -m venv .venv && source .venv/bin/activate   # on Windows: .venv\Scripts\activate
pip install -U pip matplotlib
# Optional for YAML configs:
pip install pyyaml
```

## Add the package
Unzip `heatpipe_demo.zip`, then either add its parent to `PYTHONPATH` or install editable:
```bash
pip install -e .
```

## Run the demo notebook
Open `heatpipe_demo_notebook.ipynb` and run all cells. It showcases fluids, pressure breakdowns, limit sweeps, geometry sweeps, and config I/O.

## Minimal scripting example
```python
from heatpipe.io import props_for
from heatpipe.models import (Regime, SectionLengths, Geometry, FlowFlags,
                           pressure_breakdown_cgs, estimate_limits_cgs)

T = 900.0
f = props_for("sodium", T)  # cgs keys
geom = Geometry(radius_cm=0.5, effective_pore_radius_cm=0.02, wavelength_cm=0.05)
L = SectionLengths(30.0, 20.0, 30.0)
flags = FlowFlags(Regime.LAMINAR, Regime.TURBULENT, Regime.LAMINAR)

pb = pressure_breakdown_cgs(L, geom, f, flags, m_dot_g_s=2.0, theta_deg=0.0)
lims = estimate_limits_cgs(L, geom, f, flags)
print(pb); print(lims)
```
