
from __future__ import annotations

import matplotlib.pyplot as plt
from typing import Optional, Dict, Sequence

def _finalize(fig, title: str, xlabel: str, ylabel: str, savepath: Optional[str] = None):
    ax = plt.gca()
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
    return fig

def plot_limits_vs_temperature(Ts, q_cap=None, q_sonic=None, q_ent=None, savepath: Optional[str] = None):
    fig = plt.figure()
    if Ts:
        if q_cap:   plt.plot(Ts, [None if v is None else v for v in q_cap], label="Capillary limit")
        if q_sonic: plt.plot(Ts, [None if v is None else v for v in q_sonic], label="Sonic limit")
        if q_ent:   plt.plot(Ts, [None if v is None else v for v in q_ent], label="Entrainment limit")
        plt.legend()
    return _finalize(fig, "Heat-pipe performance limits vs temperature", "Evaporator exit temperature [K]", "Q limit [W]", savepath)

def plot_axial_temperature(x_cm, T_K, savepath: Optional[str] = None):
    fig = plt.figure()
    plt.plot(x_cm, T_K, marker="o")
    return _finalize(fig, "Axial temperature profile", "Axial position x [cm]", "Temperature [K]", savepath)

def plot_cumulative_pressure(x_cm, dp_dyn_cm2, savepath: Optional[str] = None):
    fig = plt.figure()
    plt.plot(x_cm, dp_dyn_cm2, marker="o")
    return _finalize(fig, "Cumulative pressure drop", "Axial position x [cm]", "Δp cumulative [dyn/cm²]", savepath)

def plot_pressure_breakdown(components: Dict[str, float], savepath: Optional[str] = None):
    fig = plt.figure()
    keys = list(components.keys())
    vals = [components[k] for k in keys]
    plt.bar(keys, vals)
    plt.xticks(rotation=30, ha="right")
    return _finalize(fig, "Pressure-drop breakdown", "Component", "Δp [dyn/cm²]", savepath)

def plot_Q_vs_geometry(x_vals: Sequence[float], q_vals_W: Sequence[float], xlabel: str = "Characteristic dimension [cm]", savepath: Optional[str] = None):
    fig = plt.figure()
    plt.plot(x_vals, q_vals_W, marker="o")
    return _finalize(fig, "Q vs geometry", xlabel, "Q [W]", savepath)
