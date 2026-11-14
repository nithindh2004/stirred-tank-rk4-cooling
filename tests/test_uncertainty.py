"""Actuator model and statistics."""
import numpy as np
import pytest

from stcool.model import P0
from stcool.uncertainty import (actual_flow, dataset, wilson_interval, wilson_lower_one_sided, required_pass_fraction,
                                pass_fraction)


def test_actual_flow_never_exceeds_pump_limit():
    S = dataset("VAL")
    for m in (0.3, 0.45, 0.5):
        md = actual_flow(m, S["z_flow"])
        assert md.max() <= P0.m_max and md.min() >= 0


def test_point_mass_at_limit_for_full_command():
    md = actual_flow(0.5, dataset("VAL")["z_flow"])
    assert np.mean(md == P0.m_max) == pytest.approx(0.5, abs=0.005)
    assert md.mean() == pytest.approx(0.5 * (1 - 0.05 / np.sqrt(2 * np.pi)), abs=2e-4)


def test_low_commands_unaffected_by_saturation():
    z = dataset("VAL")["z_flow"]
    assert np.allclose(actual_flow(0.3, z), 0.3 * (1 + 0.05 * z))


def test_wilson_bounds():
    lo, hi = wilson_interval(9500, 10000)
    assert lo < 0.95 < hi
    lb = wilson_lower_one_sided(9500, 10000)
    assert lo < lb < 0.95


def test_required_pass_fraction_meets_lower_bound():
    n = 20_000
    p_req = required_pass_fraction(0.95, n)
    assert wilson_lower_one_sided(p_req * n, n) == pytest.approx(0.95, abs=1e-9)
    assert 0.952 < p_req < 0.953


def test_pass_fraction_monotone_in_speed():
    S = dataset("OPT")
    pf = [pass_fraction(r, 0.5, S) for r in (150, 200, 250, 300, 350)]
    assert np.all(np.diff(pf) >= 0)
