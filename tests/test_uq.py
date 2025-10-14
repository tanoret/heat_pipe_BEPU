
from heatpipe.uq import UncertainInput, UQSpec, uq_monte_carlo

def test_uq_mc_runs():
    spec = UQSpec(inputs=[
        UncertainInput("radius_cm", "normal", {"mean":0.5, "std":0.02}, lower=0.3, upper=0.8),
        UncertainInput("effective_pore_radius_cm", "uniform", {"low":0.01, "high":0.04}),
        UncertainInput("wavelength_cm", "uniform", {"low":0.04, "high":0.06}),
        UncertainInput("L_e", "deterministic", {"value":30.0}),
        UncertainInput("L_a", "deterministic", {"value":20.0}),
        UncertainInput("L_c", "deterministic", {"value":30.0}),
        UncertainInput("T_K", "normal", {"mean":850.0, "std":30.0})
    ])
    X, outs = uq_monte_carlo(spec, n=20, seed=1)
    assert "Q_lim_W" in outs and outs["Q_lim_W"].shape[0] == 20
