"""RK4 implementation: accuracy, order and stability."""
import numpy as np
import pytest

from stcool.model import P0, rk4, trajectory_exact, original_model_TL
from stcool.simulation import convergence_study, rk4_stability_limit, original_model_rk4_error


def test_rk4_matches_exact_solution():
    t, Y = rk4(300, 0.2, 1.0)
    TL, TJ = trajectory_exact(300, 0.2, t)
    assert np.abs(Y[:, 0] - TL).max() < 1e-7 and np.abs(Y[:, 1] - TJ).max() < 1e-7


def test_fourth_order_convergence():
    _, order = convergence_study(300, 0.2)
    assert 3.8 < order < 4.6


def test_rk4_unstable_beyond_stability_limit():
    h_max = rk4_stability_limit(300, 0.5)
    _, Y_ok = rk4(300, 0.5, 0.8 * h_max, t_end=P0.tf)
    assert np.all(np.isfinite(Y_ok)) and Y_ok[:, 1].max() < 400
    _, Y_bad = rk4(300, 0.5, 1.3 * h_max, t_end=1800 - 1800 % (1.3 * h_max))
    assert np.abs(Y_bad[:, 1]).max() > 1e3          # amplification > 1 for the fast mode


def test_original_model_reproduced():
    assert original_model_rk4_error() < 1e-10
    assert original_model_TL(1800) == pytest.approx(354.061486, abs=1e-6)
