from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from numpy.random import default_rng
from numpy.linalg import lstsq, cholesky, LinAlgError
from numpy import sqrt

# --- Helper: inverse standard normal CDF (Acklam) and Phi using math.erf ---
import math as _math
import numpy as _np

def _norm_ppf(u: _np.ndarray) -> _np.ndarray:
    # Return z = Phi^{-1}(u) for u in (0,1) using Acklam's rational approximation.
    u = _np.asarray(u, dtype=float)
    u = _np.clip(u, 1e-12, 1-1e-12)

    a = [-3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
          1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00]
    b = [-5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
          6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00]
    d = [ 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
          3.754408661907416e+00]

    plow = 0.02425
    phigh = 1 - plow

    z = _np.empty_like(u)

    mask_low = u < plow
    if _np.any(mask_low):
        q = _np.sqrt(-2 * _np.log(u[mask_low]))
        z[mask_low] = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                       ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)

    mask_c = (u >= plow) & (u <= phigh)
    if _np.any(mask_c):
        q = u[mask_c] - 0.5
        r = q*q
        z[mask_c] = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
                     (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)

    mask_high = u > phigh
    if _np.any(mask_high):
        q = _np.sqrt(-2 * _np.log(1 - u[mask_high]))
        z[mask_high] = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
                         ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    return z

def _phi_from_z(z: _np.ndarray) -> _np.ndarray:
    # Standard normal CDF Phi(z) using math.erf (vectorized)
    return 0.5 * (1.0 + _np.vectorize(_math.erf)(z / _math.sqrt(2.0)))


# -----------------------------
# Distributions and transforms
# -----------------------------

@dataclass
class UncertainInput:
    name: str
    dist: str  # "normal" | "uniform" | "lognormal" | "triangular" | "deterministic"
    params: Dict[str, float]  # e.g., {"mean":..., "std":...} or {"low":...,"high":...}
    lower: Optional[float] = None
    upper: Optional[float] = None

    def ppf(self, u: np.ndarray) -> np.ndarray:
        """Inverse CDF transform U ~ U(0,1) -> X in physical units."""
        d = self.dist.lower()
        u = np.clip(u, 1e-12, 1-1e-12)
        if d == "deterministic":
            val = float(self.params.get("value", 0.0))
            return np.full_like(u, val, dtype=float)
        if d == "uniform":
            a = float(self.params["low"]); b = float(self.params["high"])
            return a + (b - a) * u
        if d == "normal":
            m = float(self.params["mean"]); s = float(self.params["std"])
            z = _norm_ppf(u)
            return m + s * z
        if d == "lognormal":
            m = float(self.params["mean"]); s = float(self.params["std"])
            z = _norm_ppf(u)
            return np.exp(m + s * z)
        if d == "triangular":
            a = float(self.params["low"]); b = float(self.params["mode"]); c = float(self.params["high"])
            Fm = (b - a) / (c - a)
            x = np.empty_like(u, dtype=float)
            mask = (u < Fm)
            x[mask] = a + np.sqrt((b - a) * (c - a) * u[mask])
            x[~mask] = c - np.sqrt((c - b) * (c - a) * (1.0 - u[~mask]))
            return x
        raise ValueError(f"Unknown dist: {self.dist}")

    def clip(self, x: np.ndarray) -> np.ndarray:
        lo = self.lower if self.lower is not None else -np.inf
        hi = self.upper if self.upper is not None else np.inf
        return np.clip(x, lo, hi)


@dataclass
class UQSpec:
    inputs: List[UncertainInput]
    correlation: Optional[np.ndarray] = None  # Gaussian copula correlation (n x n)

    def names(self) -> List[str]:
        return [u.name for u in self.inputs]

    def index(self, name: str) -> int:
        return self.names().index(name)


# -----------------------------
# Sampling engines
# -----------------------------

def _lhs_u01(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    """Latin Hypercube in [0,1]^d."""
    U = np.zeros((n, d), dtype=float)
    for j in range(d):
        bins = (np.arange(n) + rng.random(n)) / n
        rng.shuffle(bins)
        U[:, j] = bins
    return U

def sample(spec: UQSpec, n: int, method: str = "lhs", seed: Optional[int] = None) -> np.ndarray:
    """
    Draw samples in physical space X for the UQSpec.
    method: 'lhs' | 'iid'
    Applies optional Gaussian copula correlation if spec.correlation is provided.
    Returns X [n,d] with columns in the order of spec.inputs.
    """
    rng = default_rng(seed)
    d = len(spec.inputs)
    # base uniforms
    if method == "lhs":
        U = _lhs_u01(n, d, rng)
    elif method == "iid":
        U = rng.random((n, d))
    else:
        raise ValueError("method must be 'lhs' or 'iid'")

    # apply Gaussian copula correlation if provided
    if spec.correlation is not None:
        R = np.array(spec.correlation, dtype=float)
        if R.shape != (d, d):
            raise ValueError("correlation must be (d,d)")
        # transform U to Z ~ N(0,1)^d, correlate via Cholesky, back to U via Phi
        Z = _norm_ppf(U)
        try:
            L = cholesky(R + 1e-12 * np.eye(d))
        except LinAlgError:
            raise ValueError("Correlation not positive definite")
        Zc = Z @ L.T
        U = _phi_from_z(Zc)

    # inverse transforms per variable
    X = np.zeros((n, d), dtype=float)
    for j, ui in enumerate(spec.inputs):
        xj = ui.ppf(U[:, j])
        X[:, j] = ui.clip(xj)
    return X


# -----------------------------
# Wilks tolerance limits
# -----------------------------

def wilks_min_samples_one_sided(beta: float, gamma: float, side: str = "upper") -> int:
    """
    Return minimal N for Wilks 1st-order one-sided tolerance limits
    with confidence 'beta' and coverage 'gamma'.
    side: 'upper' uses the maximum order statistic; 'lower' uses the minimum.
    Formulas:
      upper:  1 - gamma^N >= beta  ->  N >= ln(1 - beta)/ln(gamma)
      lower:  1 - (1 - gamma)^N >= beta  ->  N >= ln(1 - beta)/ln(1 - gamma)
    """
    if not (0 < beta < 1 and 0 < gamma < 1):
        raise ValueError("beta, gamma must be in (0,1)")
    if side == "upper":
        N = math.log(1.0 - beta) / math.log(gamma)
    elif side == "lower":
        N = math.log(1.0 - beta) / math.log(1.0 - gamma)
    else:
        raise ValueError("side must be 'upper' or 'lower'")
    return int(math.ceil(N))

def wilks_order_index_two_sided(N: int, beta: float, gamma: float) -> int:
    """
    Return the Wilks order 'k' (symmetric two-sided with orders k and N-k+1) by
    solving sum_{i=0}^{k-1} C(N,i) gamma^i (1-gamma)^{N-i} <= 1 - beta.
    """
    from math import comb
    for k in range(1, N//2 + 1):
        cdf = sum(comb(N, i) * (gamma ** i) * ((1 - gamma) ** (N - i)) for i in range(k))
        if cdf <= 1.0 - beta:
            return k
    return max(1, N // 20)


# -----------------------------
# Polynomial Chaos Expansion (Hermite PCE on Gaussian space)
# -----------------------------

def _hermite_normed(z: np.ndarray, order: int) -> np.ndarray:
    """
    Return [psi_0(z), ..., psi_order(z)] where psi_n are orthonormal
    probabilists' Hermite polynomials on N(0,1): psi_n = He_n / sqrt(n!).
    """
    z = np.asarray(z, dtype=float)
    Ps = np.zeros((z.shape[0], order + 1), dtype=float)
    Ps[:, 0] = 1.0
    if order >= 1:
        Ps[:, 1] = z
    for n in range(1, order):
        Ps[:, n + 1] = z * Ps[:, n] - n * Ps[:, n - 1]
    # normalize by sqrt(n!)
    fact = 1.0
    for n in range(1, order + 1):
        fact *= n
        Ps[:, n] /= math.sqrt(fact)
    return Ps  # (N, order+1)

def _multiindex_total_degree(d: int, p: int) -> List[Tuple[int, ...]]:
    """All multi-indices a in N^d with total degree <= p, in lexicographic order."""
    idx: List[Tuple[int, ...]] = []
    def rec(prefix, remaining_degree, i):
        if i == d:
            idx.append(tuple(prefix))
            return
        for k in range(remaining_degree + 1):
            prefix.append(k)
            rec(prefix, remaining_degree - k, i + 1)
            prefix.pop()
    for tot in range(0, p + 1):
        rec([], tot, 0)
    return idx

@dataclass
class PCEModel:
    coeffs: np.ndarray          # shape (M,)
    multi_index: List[Tuple[int,...]]  # length M
    order: int
    input_names: List[str]

    @property
    def mean(self) -> float:
        return float(self.coeffs[0])

    @property
    def variance(self) -> float:
        c = self.coeffs
        return float(np.sum(c[1:]**2))

    def predict(self, Z: np.ndarray) -> np.ndarray:
        """Evaluate surrogate on Z ~ N(0,1)^d, Z shape (N,d)"""
        N, d = Z.shape
        # precompute univariate psi up to order for each dim
        Ps_list = [_hermite_normed(Z[:, j], self.order) for j in range(d)]  # list of (N,order+1)
        A = np.ones((N, len(self.multi_index)), dtype=float)
        for m, alpha in enumerate(self.multi_index):
            for j, a in enumerate(alpha):
                A[:, m] *= Ps_list[j][:, a]
        return A @ self.coeffs

    def sobol_first_total(self) -> Tuple[np.ndarray, np.ndarray]:
        """First-order and total Sobol indices from orthonormal Hermite PCE."""
        V = self.variance
        if V <= 0: 
            return np.zeros(len(self.input_names)), np.zeros(len(self.input_names))
        first = np.zeros(len(self.input_names))
        total = np.zeros(len(self.input_names))
        for m, alpha in enumerate(self.multi_index):
            if m == 0:  # constant term
                continue
            active = [i for i,a in enumerate(alpha) if a > 0]
            cm2 = self.coeffs[m]**2
            for i in active:
                total[i] += cm2
            if len(active) == 1:
                first[active[0]] += cm2
        return first / V, total / V


def fit_pce(
    X: np.ndarray, y: np.ndarray, spec: UQSpec, order: int = 2
) -> Tuple[PCEModel, Dict[str, np.ndarray]]:
    """
    Fit a Hermite PCE on Gaussian space via isoprobabilistic transform:
    1) For each physical X_ij, compute u_ij via ranks in [0,1].
    2) Map U -> Z via Z = Phi^{-1}(U) (Acklam approx).
    3) Build tensor-product Hermite basis up to 'order' and solve least squares A c = y.
    Returns the PCE model and the Z used during fitting.
    """
    n, d = X.shape
    # Empirical CDF via ranks to get U in [0,1] per input
    U = np.zeros_like(X, dtype=float)
    for j, ui in enumerate(spec.inputs):
        xj = X[:, j]
        ranks = np.argsort(np.argsort(xj))
        U[:, j] = (ranks + 0.5) / n
    # Map to standard normal
    Z = _norm_ppf(U)

    # Build basis
    multi = _multiindex_total_degree(d, order)  # length M
    # Precompute univariate psi for each dim
    Ps_list = [_hermite_normed(Z[:, j], order) for j in range(d)]  # list of (N,order+1)
    A = np.ones((n, len(multi)), dtype=float)
    for m, alpha in enumerate(multi):
        for j, a in enumerate(alpha):
            A[:, m] *= Ps_list[j][:, a]
    # Solve least squares
    coeffs, *_ = lstsq(A, y, rcond=None)
    model = PCEModel(coeffs=coeffs, multi_index=multi, order=order, input_names=spec.names())
    return model, {"Z": Z}


# -----------------------------
# Gaussian Process (simple RBF GP)
# -----------------------------

@dataclass
class GPModel:
    X_train: np.ndarray  # (N,d), inputs normalized to [0,1]
    y_train: np.ndarray  # (N,)
    lengthscale: np.ndarray  # (d,)
    sigma_f2: float
    sigma_n2: float
    L: np.ndarray        # Cholesky factor of K (N,N)
    _min: np.ndarray
    _max: np.ndarray

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Predict mean and variance at X (in original units)."""
        Xt = _normalize(X, self._min, self._max)
        K_star = _rbf_kernel(self.X_train, Xt, self.lengthscale, self.sigma_f2)
        alpha = np.linalg.solve(self.L.T, np.linalg.solve(self.L, self.y_train))
        mu = K_star.T @ alpha
        # variance
        v = np.linalg.solve(self.L, K_star)
        k_ss = _rbf_self(Xt, self.lengthscale, self.sigma_f2) + 1e-12
        var = k_ss - np.sum(v * v, axis=0)
        var = np.maximum(var, 1e-12)
        return mu, var

def _normalize(X: np.ndarray, xmin: np.ndarray, xmax: np.ndarray) -> np.ndarray:
    return (X - xmin) / np.maximum(xmax - xmin, 1e-12)

def _rbf_kernel(X1: np.ndarray, X2: np.ndarray, l: np.ndarray, sigma_f2: float) -> np.ndarray:
    # squared distance with ARD
    d = (X1[:, None, :] - X2[None, :, :]) / l[None, None, :]
    r2 = np.sum(d**2, axis=2)
    return sigma_f2 * np.exp(-0.5 * r2)

def _rbf_self(X: np.ndarray, l: np.ndarray, sigma_f2: float) -> np.ndarray:
    return np.full(X.shape[0], sigma_f2)

def fit_gp(X: np.ndarray, y: np.ndarray, jitter: float = 1e-8) -> GPModel:
    """
    Fit a basic GP with RBF kernel.
    Heuristics:
      - Normalize inputs to [0,1] per feature.
      - Lengthscale l_j = 0.2
      - sigma_f^2 = Var(y) (or 1 if degenerate)
      - sigma_n^2 = jitter * Var(y) + 1e-12
    """
    X = np.asarray(X, float); y = np.asarray(y, float).ravel()
    xmin = X.min(axis=0); xmax = X.max(axis=0)
    Xn = _normalize(X, xmin, xmax)
    l = 0.2 * np.ones(X.shape[1])
    sf2 = float(np.var(y)) if np.var(y) > 0 else 1.0
    sn2 = float(jitter * (np.var(y) + 1.0))
    K = _rbf_kernel(Xn, Xn, l, sf2) + (sn2 + 1e-12) * np.eye(X.shape[0])
    L = cholesky(K + 1e-12 * np.eye(K.shape[0]))
    return GPModel(X_train=Xn, y_train=y, lengthscale=l, sigma_f2=sf2, sigma_n2=sn2, L=L, _min=xmin, _max=xmax)


# -----------------------------
# Wrapper: run model and collect outputs
# -----------------------------

def default_run_function(x_row: Dict[str, float]) -> Dict[str, float]:
    """
    Example run function that expects (subset of) the following keys:
    - radius_cm, effective_pore_radius_cm, wavelength_cm, L_e, L_a, L_c, T_K, fluid
    Returns dict with 'Q_lim_W' and components.
    Users may override this with a project-specific function.
    """
    from .io import props_for
    from .models import Geometry, SectionLengths, FlowFlags, Regime, estimate_limits_cgs
    T = float(x_row.get("T_K", 850.0))
    geom = Geometry(
        radius_cm=float(x_row.get("radius_cm", 0.5)),
        t_a_cm=float(x_row.get("t_a_cm", 0.05)),
        t_w_cm=float(x_row.get("t_w_cm", 0.05)),
        effective_pore_radius_cm=float(x_row.get("effective_pore_radius_cm", 0.02)),
        wavelength_cm=float(x_row.get("wavelength_cm", 0.05)),
        passages=int(x_row.get("passages", 1)),
    )
    L = SectionLengths(L_e=float(x_row.get("L_e", 30.0)),
                       L_a=float(x_row.get("L_a", 20.0)),
                       L_c=float(x_row.get("L_c", 30.0)))
    flags = FlowFlags(vapor_regime_evap=Regime.LAMINAR, vapor_regime_adiab=Regime.TURBULENT, vapor_regime_cond=Regime.LAMINAR)
    f = props_for(str(x_row.get("fluid", "sodium")), T)
    lims = estimate_limits_cgs(L, geom, f, flags)
    qvals = [q for q in [lims.q_capillary_W, lims.q_sonic_W, lims.q_entrainment_W] if q is not None]
    Qlim = float(min(qvals)) if qvals else float("nan")
    return {
        "Q_cap_W": float(lims.q_capillary_W) if lims.q_capillary_W is not None else float("nan"),
        "Q_sonic_W": float(lims.q_sonic_W) if lims.q_sonic_W is not None else float("nan"),
        "Q_entr_W": float(lims.q_entrainment_W) if lims.q_entrainment_W is not None else float("nan"),
        "Q_lim_W": Qlim,
    }

def run_batch(spec: UQSpec, X: np.ndarray, run_fn: Callable[[Dict[str,float]], Dict[str,float]]) -> Dict[str, np.ndarray]:
    """
    Evaluate the wrapped model on samples X and return a dict of outputs stacked as arrays.
    """
    names = spec.names()
    outs: Dict[str, List[float]] = {}
    for i in range(X.shape[0]):
        row = {names[j]: float(X[i, j]) for j in range(X.shape[1])}
        out = run_fn(row)
        for k, v in out.items():
            outs.setdefault(k, []).append(float(v))
    return {k: np.asarray(v, dtype=float) for k, v in outs.items()}


# -----------------------------
# High-level UQ workflows
# -----------------------------

def uq_monte_carlo(spec: UQSpec, n: int, run_fn: Callable[[Dict[str,float]], Dict[str,float]] = default_run_function,
                   method: str = "lhs", seed: Optional[int] = None) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    X = sample(spec, n, method=method, seed=seed)
    outs = run_batch(spec, X, run_fn)
    return X, outs

def uq_pce(spec: UQSpec, n_train: int, order: int = 2,
           run_fn: Callable[[Dict[str,float]], Dict[str,float]] = default_run_function,
           target: str = "Q_lim_W", method: str = "lhs", seed: Optional[int] = None) -> Tuple[PCEModel, Dict[str, np.ndarray]]:
    X = sample(spec, n_train, method=method, seed=seed)
    y = run_batch(spec, X, run_fn)[target]
    model, extras = fit_pce(X, y, spec, order=order)
    return model, {"X": X, "y": y, **extras}

def uq_gp(spec: UQSpec, n_train: int,
          run_fn: Callable[[Dict[str,float]], Dict[str,float]] = default_run_function,
          target: str = "Q_lim_W", method: str = "lhs", seed: Optional[int] = None) -> Tuple[GPModel, Dict[str, np.ndarray]]:
    X = sample(spec, n_train, method=method, seed=seed)
    y = run_batch(spec, X, run_fn)[target]
    gp = fit_gp(X, y)
    return gp, {"X": X, "y": y}
