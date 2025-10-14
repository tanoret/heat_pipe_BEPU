# Validation & QA

## Built‑in checks
- **TSAT inversion**: verify `tsat_from_p_cgs(fluid, properties_cgs(fluid,T).pv) ≈ T` over a temperature grid.
- **Dimensional sanity**: ensure Δp terms scale with \(L/D\), \(V^2\), and \(\sin\theta\) as expected.
- **Fanno identities**: confirm \(F(M)\to 0\) as \(M\to 1^-\); monotonicity on subsonic branch; pressure ratio continuity.

## Golden tests (to add)
- Reproduce the heatpipe sample problems: tables and plots with agreed tolerances per variable.
- Add regressions for limit curves vs temperature and Q vs geometry sweeps.

## Acceptance criteria (first pass)
1) Options-style limit sweeps produce monotone, physically reasonable curves vs \(T\).  
2) Pressure budgets are non-negative and consistent under unit changes.  
3) `io.props_for` keys align with what `models` expects.
