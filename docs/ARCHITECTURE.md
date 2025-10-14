# Architecture

```
heatpipe/
  fluids.py     # fluid correlations & property bundles (cgs; optional SI wrapper)
  models.py     # Δp components, Fanno/sonic, capillary & entrainment limits
  io.py         # dataclass configs, JSON/YAML load/save, property bridge
  plotting.py   # matplotlib helpers (single-plot functions)
```

**Principles**
- Keep physics in **cgs** internally; only convert at edges.
- Pure/stateless functions where possible for easy tests.
- Small dataclasses for geometry/sections/flags.
- Modules independent: change fluids or geometry without touching I/O or plotting.
