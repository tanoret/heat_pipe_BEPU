# HTPIPE UQ Layer (`htpipe.uq`)

This module adds Best Estimate Plus Uncertainty (BEPU) capabilities on top of the deterministic hydrodynamic models. It includes Monte Carlo with Wilks tolerance logic, non-intrusive Hermite PCE with Sobol indices, and a simple Gaussian Process surrogate.

## APIs

```python
# Spec
UncertainInput(name, dist, params, lower=None, upper=None)
UQSpec(inputs: List[UncertainInput], correlation: np.ndarray|None=None)

# Sampling
sample(spec, n, method="lhs"|"iid", seed=None)

# Wilks
wilks_min_samples_one_sided(beta, gamma, side="upper"|"lower")
wilks_order_index_two_sided(N, beta, gamma)

# Workflows
uq_monte_carlo(spec, n, run_fn=default_run_function, method="lhs", seed=None)
uq_pce(spec, n_train, order=2, run_fn=..., target="Q_lim_W", method="lhs", seed=None)
uq_gp(spec, n_train, run_fn=..., target="Q_lim_W", method="lhs", seed=None)
```

### Default run function
`default_run_function(row: Dict[str,float]) -> Dict[str,float]` expects keys like `radius_cm`, `effective_pore_radius_cm`, `wavelength_cm`, `L_e`, `L_a`, `L_c`, `T_K`, and `fluid`. It computes `Q_cap_W`, `Q_sonic_W`, `Q_entr_W`, and `Q_lim_W` via the existing `htpipe` models.

## Notes
- Correlation is applied using a Gaussian copula with Cholesky factorization.
- PCE uses a rank-based isoprobabilistic transform to $Z\sim \mathcal N(0,1)$, fits an orthonormal Hermite basis (total order), then derives Sobol indices analytically.
- GP uses an RBF kernel with heuristic hyperparameters for quick emulation.
