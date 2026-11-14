"""Soft sensor and adaptive policy behaviour."""
import numpy as np
import pytest

from stcool.model import P0, U_coef
from stcool.uncertainty import dataset, actual_flow
from stcool.adaptive_policy import run_policy


def _subset(n=2000):
    return {k: v[:n] if isinstance(v, np.ndarray) else v for k, v in dataset("CAL").items()}


def test_soft_sensor_is_nearly_unbiased_without_noise():
    r = run_policy(_subset(), 200, 0.5, 600, 0.2, sigma=0.0)
    err = r["UA_hat"] / r["UA_true"] - 1
    assert abs(err.mean()) < 0.01 and err.std() < 0.01


def test_speed_respects_bounds_and_cap_flag():
    r = run_policy(_subset(), 200, 0.5, 600, 0.2, noise_seed=2)
    assert r["rpm2"].min() >= P0.rpm_min and r["rpm2"].max() <= P0.rpm_max
    assert np.all(r["rpm2"][r["capped"]] == P0.rpm_max)


def test_larger_margin_never_lowers_reliability():
    S = _subset()
    a = run_policy(S, 200, 0.5, 600, 0.1, noise_seed=4)["pass_frac"]
    b = run_policy(S, 200, 0.5, 600, 0.4, noise_seed=4)["pass_frac"]
    assert b >= a


def test_practical_margin_rounds_up_to_declared_step():
    from stcool.adaptive_policy import practical_margin, DELTA_INCREMENT
    assert DELTA_INCREMENT == 0.05
    assert practical_margin(0.1862) == pytest.approx(0.20)
    assert practical_margin(0.20) == pytest.approx(0.20)       # already on the grid: unchanged
    assert practical_margin(0.2001) == pytest.approx(0.25)


def test_adaptive_validation_uses_val2_only_after_failure():
    from stcool.adaptive_policy import validation_protocol
    S = _subset(3000)
    S_val2 = {k: v[:3000] if isinstance(v, np.ndarray) else v for k, v in dataset("STRESS").items()}
    seeds = {"VAL": 1, "VAL2": 2}
    cfg = dict(rpm1=200.0, t_obs=600.0, delta=1.0, m_cmd=0.5)          # generous margin: passes first time
    final, _, s, trail = validation_protocol(cfg, S, S_val2, seeds)
    assert len(trail) == 1 and final["delta"] == 1.0 and s["accepted"]
    cfg = dict(rpm1=200.0, t_obs=600.0, delta=-1.0, m_cmd=0.5)         # negative margin: must fail VAL
    final, _, s, trail = validation_protocol(cfg, S, S_val2, seeds)
    assert len(trail) == 2 and trail[1]["dataset"] == "STRESS" and final["delta"] == pytest.approx(-0.95)


def test_matlab_uses_the_python_selected_configuration():
    """The Octave replication hard-codes the adaptive configuration; it must equal the Python CAL selection."""
    import re
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).resolve().parents[1]
    scan = root / "data" / "adaptive_policy" / "configuration_scan_CAL.csv"
    if not scan.exists():
        pytest.skip("run `python run_all.py` first")
    sel = pd.read_csv(scan).query("selected").iloc[0]
    m = (root / "src" / "matlab" / "run_all_matlab.m").read_text()
    assert float(re.search(r"rpm1 = (\d+);", m).group(1)) == sel.rpm1
    assert float(re.search(r"tobs = (\d+);", m).group(1)) == sel.t_obs_s
