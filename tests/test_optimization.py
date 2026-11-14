"""Optimisation logic and the validation protocol."""
import numpy as np
import pytest

from stcool.model import P0, final_TL_exact
from stcool.optimization import (grid_search, nominal_optimum, chance_constrained_optimum, practical_setpoint,
                                 validation_protocol, min_rpm_meeting)
from stcool.uncertainty import dataset, pass_fraction, required_pass_fraction


def test_grid_reproduces_v1_point():
    g = grid_search()
    assert (g["rpm"], g["m"]) == (125.0, 0.49)


def test_nominal_optimum_is_on_constraint_and_flow_bound():
    nom = nominal_optimum()
    assert nom["T"] == pytest.approx(P0.Ttarget, abs=1e-6)
    assert nom["at_flow_bound"]
    assert nom["E"] <= grid_search()["E"]


def test_chance_constrained_optimum_meets_requirement_minimally():
    S = dataset("OPT")
    best, _ = chance_constrained_optimum(0.95, S)
    p_req = required_pass_fraction(0.95, len(S["T0"]))
    assert pass_fraction(best["rpm"], best["m"], S) >= p_req
    assert pass_fraction(best["rpm"] - 1.0, best["m"], S) < p_req


def test_rounding_uses_only_optimisation_data_and_rounds_up():
    S = dataset("OPT")
    best, _ = chance_constrained_optimum(0.95, S)
    sp = practical_setpoint(best, 0.95, S)
    assert sp["rpm"] >= best["rpm"] and sp["rpm"] % 5 == 0


def test_validation_protocol_uses_val2_only_after_failure():
    VAL, VAL2 = dataset("VAL"), dataset("VAL2")
    ok, trail = validation_protocol(dict(rpm=400.0, m=0.5), VAL, VAL2)
    assert len(trail) == 1 and trail[0]["accepted"]
    bad, trail = validation_protocol(dict(rpm=150.0, m=0.5), VAL, VAL2)
    assert len(trail) == 2 and trail[1]["dataset"] == "VAL2" and bad["rpm"] == 155.0
