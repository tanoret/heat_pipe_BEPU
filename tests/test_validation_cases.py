import pytest

from heatpipe.io import props_for
from heatpipe.models import pressure_breakdown_cgs, estimate_limits_cgs
from heatpipe.validation import ReferenceDatum, validation_cases


CASES = validation_cases()


def _assert_with_reference(key: str, datum: ReferenceDatum, actual: float) -> None:
    if datum.abs_tol is not None:
        diff = abs(actual - datum.value)
        assert (
            diff <= datum.abs_tol
        ), f"{key}: abs err {diff:.3e} exceeds {datum.abs_tol:.3e}"
    if datum.rel_tol is not None:
        if datum.value == 0.0:
            raise AssertionError(f"{key}: rel_tol provided for zero reference value")
        rel_err = abs(actual - datum.value) / abs(datum.value)
        assert (
            rel_err <= datum.rel_tol
        ), f"{key}: rel err {rel_err:.3e} exceeds {datum.rel_tol:.3e}"


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_pressure_breakdown_regressions(case):
    if not case.pressure_refs:
        pytest.skip(f"{case.name} has no pressure_refs")

    fluids = props_for(case.fluid, case.temperature_K)
    limits = estimate_limits_cgs(
        case.lengths,
        case.geometry,
        fluids,
        case.fluid,
        theta_deg=case.theta_deg,
    )
    m_dot_g_s = limits.q_capillary_W / fluids["h_fg"] # type: ignore
    breakdown = pressure_breakdown_cgs(
        case.fluid,
        case.lengths,
        case.geometry,
        fluids,
        m_dot_g_s=m_dot_g_s,
        theta_deg=case.theta_deg,
    )

    for key, datum in case.pressure_refs.items():
        actual = getattr(breakdown, key)
        _assert_with_reference(key, datum, actual)


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.name)
def test_limit_regressions(case):
    fluids = props_for(case.fluid, case.temperature_K)
    limits = estimate_limits_cgs(
        case.lengths,
        case.geometry,
        fluids,
        case.fluid,
        theta_deg=case.theta_deg,
    )

    for key, datum in case.limit_refs_W.items():
        actual = getattr(limits, key)
        assert actual is not None, f"{key}: expected finite value"
        _assert_with_reference(key, datum, actual)