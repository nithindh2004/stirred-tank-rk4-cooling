"""First-law closure and energy accounting."""
import numpy as np
import pytest

from stcool.model import P0, batch_energy_kWh, agitator_power, pump_power
from stcool.simulation import energy_balance
from stcool.uncertainty import dataset, mc_energy_kWh, actual_flow
from stcool.adaptive_policy import run_policy


@pytest.mark.parametrize("rpm,m", [(300, 0.2), (295, 0.5), (120, 0.5)])
def test_energy_balance_closes(rpm, m):
    eb = energy_balance(rpm, m)
    assert eb["relative_residual_simpson"] < 1e-8


def test_batch_energy_is_power_times_time():
    assert batch_energy_kWh(300, 0.2) == pytest.approx((2025 + pump_power(0.2)) * 1800 / 3.6e6)


def test_monte_carlo_pump_energy_uses_actual_flow():
    S = dataset("OPT")
    E = mc_energy_kWh(295, 0.5, S)
    md = actual_flow(0.5, S["z_flow"])
    assert np.allclose(E, (agitator_power(295) + pump_power(md)) * 1800 / 3.6e6)
    assert E.mean() < batch_energy_kWh(295, 0.5)    # saturation can only lower the flow


def test_adaptive_energy_is_integral_of_piecewise_power():
    S = {k: v[:500] if isinstance(v, np.ndarray) else v for k, v in dataset("CAL").items()}
    r = run_policy(S, 200, 0.5, 600, 0.2, noise_seed=1)
    md = actual_flow(0.5, S["z_flow"])
    expect = (agitator_power(200) * 600 + agitator_power(r["rpm2"]) * 1200 + pump_power(md) * 1800) / 3.6e6
    assert np.allclose(r["E"], expect)
