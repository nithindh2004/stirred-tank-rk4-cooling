"""Uncertainty model, coolant actuator model, Monte Carlo evaluation and statistics.

Coolant flow: commanded vs actual
---------------------------------
The operator commands a flow m_cmd in [m_min, m_max]. A flow-control loop
(flow meter + control valve) delivers it with a proportional error, 5 % (1 sd)
of the command. The pump cannot deliver more than m_max, so

    m_actual = min(m_max, max(0, m_cmd (1 + 0.05 z))),   z ~ N(0, 1).

This is a *saturated actuator*: if the loop asks for more than the pump can
give, the valve sits fully open and the flow equals m_max. The flow
distribution therefore has a point mass at m_max. For m_cmd = m_max this mass
is P(z >= 0) = 50 %, the remaining 50 % lies below m_max, and
E[m_actual] = m_max (1 - 0.05 / sqrt(2 pi)) = 0.980 m_max.

A truncated distribution was considered and rejected: truncation rescales the
error density so that the loop would appear *more accurate* near the limit,
which has no physical basis. Saturation keeps the controller-error model and
adds the hard equipment limit. The lower clip at zero never binds for
m_cmd >= 0.05 kg/s (it would need z < -20).

Other uncertain inputs (independent): U_ref ~ N(150, 15^2) W/(m2 K),
M_J ~ N(10, 1^2) kg, T_0 ~ N(370, 1^2) K. The guards max(1, .) on U_ref and
M_J are > 9 sd away and never bind.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

from .model import P0, Params, U_coef, agitator_power, pump_power, propagate_exact

FLOW_ERROR_REL_SD = 0.05

UNCERTAINTY_TABLE = [
    # name, distribution, mean, sd, unit, basis
    ("U_ref", "Normal", 150.0, 15.0, "W/(m2 K)", "assumed: 10 % (1 sd) uncertainty in the reference coefficient"),
    ("M_J", "Normal", 10.0, 1.0, "kg", "assumed: 10 % (1 sd) uncertainty in effective jacket inventory"),
    ("flow error z", "Normal", 0.0, 1.0, "-", "assumed: actual = min(m_max, m_cmd (1 + 0.05 z)) (flow-loop error, pump limit)"),
    ("T_0", "Normal", 370.0, 1.0, "K", "assumed: 1 K (1 sd) batch-to-batch variation of the initial temperature"),
]

# Fixed, documented datasets. Each has one purpose only. VAL and VAL2 are
# final-release seeds that were never examined during development: the
# development runs used seeds 20261110, 20261111 and 20261115 for validation,
# and those seeds are retired (see CHANGELOG.md).
DATASETS = {
    "OPT":    dict(seed=20261109, n=20_000,  purpose="optimisation of fixed settings (chance-constrained search, rounding check)"),
    "CAL":    dict(seed=20261112, n=20_000,  purpose="calibration and configuration selection of the adaptive policy"),
    "VAL":    dict(seed=20261120, n=200_000, purpose="independent validation of every selected design and policy (used once)"),
    "VAL2":   dict(seed=20261121, n=200_000, purpose="fresh validation, used only if a design fails VAL and is modified by the pre-declared rule"),
    "STRESS": dict(seed=20261114, n=50_000,  purpose="model-form and robustness stress tests (evaluation only, no selection)"),
}
NOISE_SEEDS = {"CAL": 31, "VAL": 40, "VAL2": 41, "STRESS": 34}
SOBOL_SEED = 20261113


def sample_inputs(n: int, seed: int, p: Params = P0) -> dict:
    """Draw the four uncertain inputs. Draw order: U_ref, M_J, z_flow, T_0."""
    g = np.random.default_rng(seed)
    return dict(
        Uref=np.maximum(1.0, p.Uref + 0.10 * p.Uref * g.standard_normal(n)),
        MJ=np.maximum(1.0, p.MJ + 0.10 * p.MJ * g.standard_normal(n)),
        z_flow=g.standard_normal(n),
        T0=p.T0 + 1.0 * g.standard_normal(n),
    )


def dataset(name: str, p: Params = P0) -> dict:
    d = DATASETS[name]
    S = sample_inputs(d["n"], d["seed"], p)
    S["name"] = name
    return S


def actual_flow(m_cmd, z_flow, p: Params = P0):
    """Saturated actuator: proportional loop error, clipped to [0, m_max]."""
    return np.clip(np.asarray(m_cmd, float) * (1.0 + FLOW_ERROR_REL_SD * np.asarray(z_flow, float)), 0.0, p.m_max)


def mc_final_TL(rpm, m_cmd, S: dict, p: Params = P0):
    """T_L(t_f) for every sample (exact solution; identical to RK4 within 1e-7 K)."""
    md = actual_flow(m_cmd, S["z_flow"], p)
    UA = U_coef(rpm, md, p, S["Uref"]) * p.A
    return propagate_exact(S["T0"], p.Tin, UA, md, S["MJ"], p.tf, p)[0]


def mc_energy_kWh(rpm, m_cmd, S: dict, p: Params = P0):
    """Per-batch auxiliary electricity with the pump at the *actual* flow."""
    md = actual_flow(m_cmd, S["z_flow"], p)
    return (agitator_power(rpm, p) + pump_power(md, p)) * p.tf / 3.6e6


def expected_energy_kWh(rpm, m_cmd, S: dict, p: Params = P0):
    return float(np.mean(mc_energy_kWh(rpm, m_cmd, S, p)))


def pass_fraction(rpm, m_cmd, S: dict, p: Params = P0) -> float:
    return float(np.mean(mc_final_TL(rpm, m_cmd, S, p) <= p.Ttarget))


# ---------------------------------------------------------------- statistics
def wilson_interval(k, n, conf=0.95):
    """Two-sided Wilson score interval."""
    z = norm.ppf(0.5 + conf / 2)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def wilson_lower_one_sided(k, n, conf=0.95):
    """One-sided lower Wilson bound (z = 1.645 for 95 %)."""
    z = norm.ppf(conf)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return c - h


def required_pass_fraction(alpha, n, conf=0.95):
    """Smallest sample pass fraction whose one-sided Wilson lower bound reaches alpha."""
    if alpha <= 0.0:
        return 0.0
    f = lambda ph: wilson_lower_one_sided(ph * n, n, conf) - alpha
    return brentq(f, alpha, 1.0 - 1e-12)


def summarise(x, target, energy=None):
    k = int(np.sum(x <= target)); n = len(x)
    lo, hi = wilson_interval(k, n)
    out = dict(n=n, pass_frac=k / n, pass_ci_lo=lo, pass_ci_hi=hi, pass_lb95_one_sided=wilson_lower_one_sided(k, n),
               mean_T=float(np.mean(x)), sd_T=float(np.std(x, ddof=1)), p05_T=float(np.quantile(x, .05)),
               median_T=float(np.median(x)), p95_T=float(np.quantile(x, .95)), p99_T=float(np.quantile(x, .99)))
    if energy is not None:
        out.update(mean_E=float(np.mean(energy)), median_E=float(np.median(energy)),
                   p95_E=float(np.quantile(energy, .95)), p99_E=float(np.quantile(energy, .99)), max_E=float(np.max(energy)))
    return out
