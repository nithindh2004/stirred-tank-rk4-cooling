"""Physical consistency of the model equations and closures."""
import numpy as np
import pytest

from stcool.model import (P0, U_coef, agitator_power, pump_power, rk4, final_TL_exact, rk4_tanks_in_series,
                          UA_effective, UA_effective_plugflow, quasi_steady_TL, analytic_boundary_rpm, eigenvalues)
from stcool.optimization import boundary_rpm


def test_closure_returns_reference_value():
    assert U_coef(P0.Nref, P0.mref) == pytest.approx(P0.Uref)


def test_power_formulas_and_units():
    # Np rho n^3 D^5 / eta at 300 rpm: 5 * 1000 * 5^3 * 0.3^5 / 0.75 = 2025 W
    assert agitator_power(300) == pytest.approx(2025.0)
    # dP Q / eta at 0.5 kg/s: 312,500 Pa * 5e-4 m3/s / 0.65 = 240.38 W
    assert pump_power(0.5) == pytest.approx(312_500 * 5e-4 / 0.65)


def test_liquid_cools_monotonically_and_jacket_stays_colder():
    t, Y = rk4(300, 0.2, 1.0)
    assert np.all(np.diff(Y[:, 0]) < 0)
    assert np.all(Y[1:, 1] < Y[1:, 0])
    assert np.all(Y[:, 1] >= P0.Tin - 1e-12)


def test_more_agitation_or_flow_cools_more():
    assert final_TL_exact(350, 0.3) < final_TL_exact(300, 0.3) < final_TL_exact(250, 0.3)
    assert final_TL_exact(300, 0.4) < final_TL_exact(300, 0.3)


def test_eigenvalues_real_negative():
    l1, l2 = eigenvalues(U_coef(300, 0.2) * P0.A, 0.2, P0.MJ)
    assert l2 < l1 < 0


def test_reduced_order_model_accuracy_when_jacket_is_fast():
    G = UA_effective(300, 0.2)
    assert abs(quasi_steady_TL(G, P0.tf) - final_TL_exact(300, 0.2)) < 0.1


def test_closed_form_boundary_close_to_simulated_boundary():
    for m in (0.3, 0.4, 0.5):
        assert abs(analytic_boundary_rpm(m) - boundary_rpm(m)) / boundary_rpm(m) < 0.05


def test_single_cell_jacket_equals_base_model():
    assert rk4_tanks_in_series(300, 0.2, 1.0, 1) == pytest.approx(final_TL_exact(300, 0.2), abs=1e-6)


def test_plug_flow_jacket_is_more_effective_than_mixed():
    assert UA_effective_plugflow(300, 0.2) > UA_effective(300, 0.2)
    assert rk4_tanks_in_series(300, 0.2, 0.5, 10) < rk4_tanks_in_series(300, 0.2, 1.0, 1)


def test_local_flow_perturbations_respect_pump_limit():
    """At an operating point already at m_max, an upward flow step is not physically available."""
    from stcool.sensitivity import local_sensitivity
    _, rows = local_sensitivity(295.0, P0.m_max)
    up = [r for r in rows if r["parameter"] == "Coolant flow"]
    assert all(r["delta_high_K"] == 0.0 for r in up)
    assert all(r["delta_low_K"] > 0.0 for r in up)
