# Module: `heatpipe.plotting`

Single-plot functions (no styles/colors set). Each returns a `matplotlib.figure.Figure` and optionally saves a PNG when `savepath` is provided.

- `plot_limits_vs_temperature(Ts, q_cap=None, q_sonic=None, q_ent=None, savepath=None)`  
- `plot_axial_temperature(x_cm, T_K, savepath=None)`  
- `plot_cumulative_pressure(x_cm, dp_dyn_cm2, savepath=None)`  
- `plot_pressure_breakdown(components: dict, savepath=None)`  
- `plot_Q_vs_geometry(x_vals, q_vals_W, xlabel="Characteristic dimension [cm]", savepath=None)`  
