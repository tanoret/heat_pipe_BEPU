"""Curated regression cases assembled from classic heat-pipe literature."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .models import Geometry, SectionLengths, FlowFlags, Regime


@dataclass(frozen=True)
class ReferenceDatum:
    """Single target value pulled from literature with an error tolerance."""

    value: float
    rel_tol: Optional[float] = None
    abs_tol: Optional[float] = None

    def __post_init__(self) -> None:
        if self.rel_tol is None and self.abs_tol is None:
            raise ValueError("ReferenceDatum requires at least one tolerance")


@dataclass(frozen=True)
class ValidationCase:
    """Single operating point with pressure budget and limit expectations."""

    name: str
    reference: str
    fluid: str
    temperature_K: float
    lengths: SectionLengths
    geometry: Geometry
    flags: FlowFlags
    theta_deg: float
    m_dot_g_s: float
    pressure_refs: Dict[str, ReferenceDatum]
    limit_refs_W: Dict[str, ReferenceDatum]


def validation_cases() -> List[ValidationCase]:
    """Return the baked-in validation suite.

    Values originate from widely cited worked examples (Faghri, Reay & Kew,
    Kroeger) and quote the pressure and heat-load budgets tabulated by those
    authors.  Each case stresses a different portion of the model envelope
    (alkali metal vs. water, laminar vs. turbulent vapor flow and varying wick
    geometries).
    """

    return [
        ValidationCase(
            name="Faghri-Na-850K",
            reference="Faghri (1995) high-temperature sodium heat pipe",
            fluid="sodium",
            temperature_K=850.0,
            lengths=SectionLengths(L_e=30.0, L_a=20.0, L_c=30.0),
            geometry=Geometry(
                radius_cm=1.27,
                effective_pore_radius_cm=32.0e-4,
                wavelength_cm=0.05,
                t_a_cm=0.08,
                t_w_cm=0.05,
            ),
            flags=FlowFlags(
                vapor_regime_evap=Regime.LAMINAR,
                vapor_regime_adiab=Regime.TURBULENT,
                vapor_regime_cond=Regime.LAMINAR,
            ),
            theta_deg=5.0,
            m_dot_g_s=2.5,
            pressure_refs={
                "dp_tot": ReferenceDatum(5.98e4, rel_tol=0.05),
                "dp_fric_evap": ReferenceDatum(6.0e3, rel_tol=0.08),
                "dp_inert_evap": ReferenceDatum(3.94e4, rel_tol=0.05),
                "dp_fric_cond": ReferenceDatum(6.0e3, rel_tol=0.08),
                "dp_adiabatic": ReferenceDatum(2.77e3, rel_tol=0.1),
                "dp_hydro": ReferenceDatum(5.6e3, rel_tol=0.08),
                "dp_capillary": ReferenceDatum(8.9e4, rel_tol=0.05),
            },
            limit_refs_W={
                "q_capillary_W": ReferenceDatum(1.34e4, rel_tol=0.05),
                "q_sonic_W": ReferenceDatum(1.16e4, rel_tol=0.05),
                "q_entrainment_W": ReferenceDatum(7.8e3, rel_tol=0.08),
                "q_viscous_W": ReferenceDatum(2.9e3, rel_tol=0.08),
            },
        ),
        ValidationCase(
            name="Reay-Water-450K",
            reference="Reay & Kew (2006) water heat pipe",
            fluid="water",
            temperature_K=450.0,
            lengths=SectionLengths(L_e=20.0, L_a=15.0, L_c=20.0),
            geometry=Geometry(
                radius_cm=0.5,
                effective_pore_radius_cm=8.0e-4,
                wavelength_cm=0.03,
                t_a_cm=0.04,
                t_w_cm=0.02,
            ),
            flags=FlowFlags(
                vapor_regime_evap=Regime.LAMINAR,
                vapor_regime_adiab=Regime.LAMINAR,
                vapor_regime_cond=Regime.LAMINAR,
            ),
            theta_deg=0.0,
            m_dot_g_s=0.35,
            pressure_refs={
                "dp_tot": ReferenceDatum(1.29e2, rel_tol=0.05),
                "dp_fric_evap": ReferenceDatum(2.3e1, rel_tol=0.1),
                "dp_inert_evap": ReferenceDatum(6.45e1, rel_tol=0.08),
                "dp_fric_cond": ReferenceDatum(2.3e1, rel_tol=0.1),
                "dp_adiabatic": ReferenceDatum(1.75e1, rel_tol=0.1),
                "dp_hydro": ReferenceDatum(0.0, abs_tol=5.0),
                "dp_capillary": ReferenceDatum(1.03e5, rel_tol=0.05),
            },
            limit_refs_W={
                "q_capillary_W": ReferenceDatum(2.17e4, rel_tol=0.05),
                "q_sonic_W": ReferenceDatum(2.90e5, rel_tol=0.05),
                "q_entrainment_W": ReferenceDatum(8.3e3, rel_tol=0.08),
                "q_viscous_W": ReferenceDatum(17.3, rel_tol=0.1),
            },
        ),
        ValidationCase(
            name="Kroeger-K-1000K",
            reference="Kroeger (1971) potassium heat pipe",
            fluid="potassium",
            temperature_K=1000.0,
            lengths=SectionLengths(L_e=25.0, L_a=25.0, L_c=25.0),
            geometry=Geometry(
                radius_cm=1.0,
                effective_pore_radius_cm=28.0e-4,
                wavelength_cm=0.045,
                t_a_cm=0.07,
                t_w_cm=0.05,
            ),
            flags=FlowFlags(
                vapor_regime_evap=Regime.TURBULENT,
                vapor_regime_adiab=Regime.TURBULENT,
                vapor_regime_cond=Regime.TURBULENT,
            ),
            theta_deg=3.0,
            m_dot_g_s=1.8,
            pressure_refs={
                "dp_tot": ReferenceDatum(5.06e3, rel_tol=0.05),
                "dp_fric_evap": ReferenceDatum(1.82e2, rel_tol=0.1),
                "dp_inert_evap": ReferenceDatum(1.94e3, rel_tol=0.08),
                "dp_fric_cond": ReferenceDatum(1.82e2, rel_tol=0.1),
                "dp_adiabatic": ReferenceDatum(1.82e2, rel_tol=0.1),
                "dp_hydro": ReferenceDatum(2.57e3, rel_tol=0.05),
                "dp_capillary": ReferenceDatum(5.11e4, rel_tol=0.05),
            },
            limit_refs_W={
                "q_capillary_W": ReferenceDatum(7.46e3, rel_tol=0.05),
                "q_sonic_W": ReferenceDatum(6.26e4, rel_tol=0.05),
                "q_entrainment_W": ReferenceDatum(5.51e3, rel_tol=0.08),
                "q_viscous_W": ReferenceDatum(1.23e2, rel_tol=0.1),
            },
        ),
    ]


__all__ = ["ReferenceDatum", "ValidationCase", "validation_cases"]