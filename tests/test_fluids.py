
from heatpipe.fluids import Fluid, properties_cgs, tsat_from_p_cgs

def test_tsat_inversion_sodium():
    for T in [700.0, 800.0, 900.0]:
        P = properties_cgs(Fluid.SODIUM, T).pv
        T_back = tsat_from_p_cgs(Fluid.SODIUM, P)
        assert abs(T_back - T) / T < 1e-2
