
from heatpipe.io import props_for
from heatpipe.models import Geometry, SectionLengths, FlowFlags, Regime, pressure_breakdown_cgs

def test_pressure_breakdown_runs():
    T = 850.0
    f = props_for("sodium", T)
    geom = Geometry(radius_cm=1.57, effective_pore_radius_cm=32.0e-4, wavelength_cm=0.05, t_a_cm=0.07, t_w_cm=0.07)
    L = SectionLengths(30.0, 20.0, 30.0)
    flags = FlowFlags(Regime.LAMINAR, Regime.TURBULENT, Regime.LAMINAR)
    pb = pressure_breakdown_cgs(L, geom, f, flags, m_dot_g_s=2.0, theta_deg=0.0)
    assert pb.dp_capillary > 0
