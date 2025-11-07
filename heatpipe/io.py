
from __future__ import annotations
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Optional, Dict, Any, Tuple

try:
    import yaml  # type: ignore
    _HAVE_YAML = True
except Exception:
    yaml = None
    _HAVE_YAML = False

from .fluids import Fluid, properties_cgs  # type: ignore

class Option(int, Enum):
    design_passages = 1
    limits_vs_temperature = 2
    axial_profiles = 3
    design_dimension = 4

@dataclass
class SectionLengths:
    L_e: float
    L_a: float
    L_c: float

@dataclass
class Geometry:
    radius_cm: float
    effective_pore_radius_cm: float = 0.0
    wavelength_cm: Optional[float] = None
    passages: int = 1

@dataclass
class PipeConfig:
    fluid: str
    theta_deg: float
    geometry: Geometry
    lengths: SectionLengths

@dataclass
class SweepSpec:
    start_K: float
    stop_K: float
    step_K: float

@dataclass
class RunSpec:
    option: Option
    sweep: Optional[SweepSpec] = None

def _to_dict(obj):
    d = asdict(obj)
    def _fix(v):
        if isinstance(v, Enum):
            return v.value
        if isinstance(v, dict):
            return {k:_fix(x) for k,x in v.items()}
        if isinstance(v, (list, tuple)):
            return [_fix(x) for x in v]
        return v
    return _fix(d)

def dumps_config(cfg: PipeConfig, run: RunSpec, fmt: str = "json") -> str:
    payload = {"pipe": _to_dict(cfg), "run": _to_dict(run)}
    fmt = fmt.lower()
    if fmt == "json":
        import json
        return json.dumps(payload, indent=2)
    elif fmt in ("yaml", "yml"):
        if not _HAVE_YAML:
            raise RuntimeError("PyYAML is not available")
        return yaml.safe_dump(payload, sort_keys=False)  # type: ignore
    else:
        raise ValueError("fmt must be 'json' or 'yaml'")

def dump_config_file(path: str, cfg: PipeConfig, run: RunSpec) -> None:
    ext = path.split(".")[-1].lower()
    s = dumps_config(cfg, run, fmt="yaml" if ext in ("yaml", "yml") else "json")
    with open(path, "w") as f:
        f.write(s)

def _enum_val(E, v):
    if isinstance(v, E):
        return v
    return E(v)

def _load_payload(d: Dict[str, Any]) -> Tuple[PipeConfig, RunSpec]:
    p = d["pipe"]
    g = p["geometry"]
    l = p["lengths"]
    f = p["flags"]
    geom = Geometry(
        radius_cm=float(g["radius_cm"]),
        effective_pore_radius_cm=float(g.get("effective_pore_radius_cm", 0.0)),
        wavelength_cm=(None if g.get("wavelength_cm") in (None, "null") else float(g.get("wavelength_cm", 0.0))),
        passages=int(g.get("passages", 1)),
    )
    cfg = PipeConfig(
        fluid=str(p["fluid"]),
        theta_deg=float(p.get("theta_deg", 0.0)),
        geometry=geom,
        lengths=SectionLengths(L_e=float(l["L_e"]), L_a=float(l["L_a"]), L_c=float(l["L_c"])),
    )
    r = d.get("run", {"option": 2})
    opt = _enum_val(Option, r.get("option", 2))
    sweep = None
    if r.get("sweep"):
        s = r["sweep"]
        sweep = SweepSpec(start_K=float(s["start_K"]), stop_K=float(s["stop_K"]), step_K=float(s["step_K"]))
    run = RunSpec(option=opt, sweep=sweep)
    return cfg, run

def load_config(path: str):
    ext = path.split(".")[-1].lower()
    with open(path, "r") as f:
        raw = f.read()
    if ext in ("yaml", "yml") and _HAVE_YAML:
        d = yaml.safe_load(raw)  # type: ignore
    else:
        import json
        d = json.loads(raw)
    return _load_payload(d)

def props_for(fluid_name: str, T_K: float):
    from .fluids import Fluid, properties_cgs
    f = Fluid(fluid_name) if not isinstance(fluid_name, Fluid) else fluid_name
    P = properties_cgs(f, T_K)
    return {"T_sat": T_K, "p_sat": P.pv, "rho_v": P.rhov, "rho_l": P.rhol, "mu_v": P.muv, "mu_l": P.mul, "sigma": P.sigma, "gamma": P.gamma, "h_fg": P.hfg, "m_w": P.mw}
