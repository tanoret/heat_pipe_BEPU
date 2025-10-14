# Module: `heatpipe.io`

## Dataclasses & enums
```python
class Regime(str, Enum): laminar, turbulent
class Option(int, Enum): design_passages=1, limits_vs_temperature=2, axial_profiles=3, design_dimension=4

@dataclass class SectionLengths: L_e:float; L_a:float; L_c:float
@dataclass class Geometry: radius_cm:float; effective_pore_radius_cm:float=0.0; wavelength_cm:float|None=None; passages:int=1
@dataclass class FlowFlags: vapor_regime_evap:Regime="laminar"; vapor_regime_adiab:Regime="laminar"; vapor_regime_cond:Regime="laminar"
@dataclass class PipeConfig: fluid:str; theta_deg:float; geometry:Geometry; lengths:SectionLengths; flags:FlowFlags
@dataclass class SweepSpec: start_K:float; stop_K:float; step_K:float
@dataclass class RunSpec: option:Option; sweep:SweepSpec|None=None
```

## Serialization
- `dump_config_file(path, cfg, run)` — JSON (.json) or YAML (.yaml/.yml if PyYAML installed).
- `load_config(path) -> (PipeConfig, RunSpec)`.

## Property bridge
- `props_for(fluid_name: str, T_K: float) -> dict` returns the model‑ready cgs dictionary.

## Example JSON
```json
{
  "pipe": {
    "fluid": "sodium",
    "theta_deg": 0.0,
    "geometry": {"radius_cm": 0.5, "effective_pore_radius_cm": 0.02, "wavelength_cm": 0.05, "passages": 6},
    "lengths": {"L_e": 30.0, "L_a": 20.0, "L_c": 30.0},
    "flags": {"vapor_regime_evap": "laminar", "vapor_regime_adiab": "turbulent", "vapor_regime_cond": "laminar"}
  },
  "run": {"option": 2, "sweep": {"start_K": 700.0, "stop_K": 950.0, "step_K": 25.0}}
}
```
