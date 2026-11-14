"""Simulation-based adaptive "measure-then-set" operating policy.

Phase 1 (0 <= t < t_obs): run at speed N1 with commanded flow m_cmd; log T_L and
T_J every `dt_meas` seconds with Gaussian sensor noise (sd `sigma`, assumed).

Estimation (soft sensor, quasi-steady, using data after `t_fit_start` once the
fast jacket transient has decayed):
    Q_hat     = -M_L Cp_L * slope(T_L)
    UA_hat    = Q_hat / (T_L - T_J)
    mCp_hat   = (Q_hat - M_J,nominal Cp dT_J/dt) / (T_J - T_in)
Both temperatures are fitted with straight lines over the window.

Phase 2 (t_obs <= t <= t_f): choose the lowest speed N2 in [rpm_min, rpm_max]
whose exact-model prediction from the estimated state reaches T_target - delta
by t_f, assuming U scales as (N2/N1)^eN_plan. If no speed suffices, run at
rpm_max (the batch is "speed-capped").

Energy is the time integral of instantaneous electrical power, which is
piecewise constant: P_ag(N1) t_obs + P_ag(N2) (t_f - t_obs) + P_pump(m_actual) t_f.
"""
from __future__ import annotations

import math

import numpy as np

from .model import P0, Params, U_coef, agitator_power, pump_power, propagate_exact
from .uncertainty import actual_flow, required_pass_fraction, summarise

DELTA_INCREMENT = 0.05   # K; practical rounding step of the safety margin, declared before validation
SENSOR_SIGMA_K = 0.10    # K; assumed temperature-sensor noise (1 sd), both T_L and T_J


def _observe_chunk(S, idx, rpm1, m_cmd, t_obs, p, sigma, dt_meas, t_fit_start, rng, truth_p):
    """Phase 1 and the soft-sensor estimate for one chunk of samples. Nothing here
    depends on the margin delta, so calibration can reuse it for every trial delta."""
    Uref, MJ, T0 = S["Uref"][idx], S["MJ"][idx], S["T0"][idx]
    md = actual_flow(m_cmd, S["z_flow"][idx], p)
    UA1 = U_coef(rpm1, md, truth_p, Uref) * p.A
    tm = np.arange(0.0, t_obs + 1e-9, dt_meas)
    TL, TJ = propagate_exact(T0[:, None], p.Tin, UA1[:, None], md[:, None], MJ[:, None], tm[None, :], p)
    yL = TL + sigma * rng.standard_normal(TL.shape)
    yJ = TJ + sigma * rng.standard_normal(TJ.shape)
    w = tm >= t_fit_start
    tw = tm[w]; tc = tw - tw.mean(); sxx = np.sum(tc ** 2)
    mL = yL[:, w].mean(1); sL = (yL[:, w] @ tc) / sxx
    mJ = yJ[:, w].mean(1); sJ = (yJ[:, w] @ tc) / sxx
    Q = -p.ML * p.CpL * sL
    return dict(Uref=Uref, MJ=MJ, md=md, UA1=UA1, TL_end=TL[:, -1], TJ_end=TJ[:, -1],
                UA_hat=Q / (mL - mJ), mCp_hat=(Q - p.MJ * p.CpC * sJ) / (mJ - p.Tin),
                TL_obs=mL + sL * (t_obs - tw.mean()), TJ_obs=mJ + sJ * (t_obs - tw.mean()))


def _decide_chunk(ob, rpm1, t_obs, delta, p, eN_plan, truth_p):
    """Phase 2 for one chunk: choose N2 from the estimate and simulate the true plant."""
    eN = p.eN if eN_plan is None else eN_plan
    goal = p.Ttarget - delta
    UA_hat, mCp_hat, TL_obs, TJ_obs = ob["UA_hat"], ob["mCp_hat"], ob["TL_obs"], ob["TJ_obs"]

    def predict(r):
        return propagate_exact(TL_obs, TJ_obs, UA_hat * (r / rpm1) ** eN, mCp_hat / p.CpC, p.MJ, p.tf - t_obs, p)[0]

    n = len(UA_hat)
    lo = np.full(n, p.rpm_min); hi = np.full(n, p.rpm_max)
    ok_min = predict(lo) <= goal
    ok_max = predict(hi) <= goal
    for _ in range(25):                           # vectorised bisection, prediction monotone in speed
        mid = 0.5 * (lo + hi)
        good = predict(mid) <= goal
        hi = np.where(good, mid, hi); lo = np.where(good, lo, mid)
    rpm2 = np.where(ok_min, p.rpm_min, np.where(ok_max, hi, p.rpm_max))
    md = ob["md"]
    UA2 = U_coef(rpm2, md, truth_p, ob["Uref"]) * p.A
    TLf = propagate_exact(ob["TL_end"], ob["TJ_end"], UA2, md, ob["MJ"], p.tf - t_obs, p)[0]
    E = (agitator_power(rpm1, p) * t_obs + agitator_power(rpm2, p) * (p.tf - t_obs) + pump_power(md, p) * p.tf) / 3.6e6
    return dict(TLf=TLf, E=E, rpm2=rpm2, capped=~ok_max, UA_true=ob["UA1"], UA_hat=UA_hat,
                mCp_true=md * p.CpC, mCp_hat=mCp_hat)


def observe(S, rpm1, m_cmd, t_obs, p: Params = P0, sigma=SENSOR_SIGMA_K, dt_meas=10.0, t_fit_start=120.0,
            noise_seed=0, truth_p: Params | None = None, chunk=50_000):
    """Phase-1 observations and estimates for every sample of S (list of chunks).
    The noise generator is consumed chunk by chunk in order, exactly as in run_policy."""
    truth_p = p if truth_p is None else truth_p
    rng = np.random.default_rng(noise_seed)
    n = len(S["T0"])
    return [_observe_chunk(S, np.arange(i, min(i + chunk, n)), rpm1, m_cmd, t_obs, p, sigma, dt_meas, t_fit_start,
                           rng, truth_p) for i in range(0, n, chunk)]


def run_policy(S, rpm1, m_cmd, t_obs, delta, p: Params = P0, sigma=SENSOR_SIGMA_K, dt_meas=10.0, t_fit_start=120.0,
               noise_seed=0, eN_plan=None, truth_p: Params | None = None, chunk=50_000, observed=None):
    """Simulate the policy on every sample of S. `truth_p` lets the true plant
    differ from the planner's model (model-form stress tests). `observed` (from
    `observe`, same arguments) skips recomputing phase 1 when only delta changes."""
    truth_p = p if truth_p is None else truth_p
    if observed is None:
        observed = observe(S, rpm1, m_cmd, t_obs, p, sigma, dt_meas, t_fit_start, noise_seed, truth_p, chunk)
    parts = [_decide_chunk(ob, rpm1, t_obs, delta, p, eN_plan, truth_p) for ob in observed]
    out = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    out["pass_frac"] = float(np.mean(out["TLf"] <= p.Ttarget))
    return out


def calibrate_delta(S, alpha, rpm1, m_cmd, t_obs, p: Params = P0, noise_seed=0, lo=0.0, hi=3.0, iters=18,
                    sigma=SENSOR_SIGMA_K, eN_plan=None, truth_p: Params | None = None):
    """Smallest margin delta whose pass fraction on S meets the one-sided 95 % lower-bound rule.
    Phase 1 is simulated once and reused for every trial delta (it does not depend on delta)."""
    p_req = required_pass_fraction(alpha, len(S["T0"]))
    ob = observe(S, rpm1, m_cmd, t_obs, p, sigma, noise_seed=noise_seed, truth_p=truth_p)
    f = lambda d: run_policy(S, rpm1, m_cmd, t_obs, d, p, sigma=sigma, eN_plan=eN_plan, truth_p=truth_p,
                             observed=ob)["pass_frac"]
    if f(hi) < p_req:
        return np.nan
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if f(mid) >= p_req:
            hi = mid
        else:
            lo = mid
    return hi


def policy_summary(res, p: Params = P0):
    s = summarise(res["TLf"], p.Ttarget, res["E"])
    s.update(frac_speed_capped=float(np.mean(res["capped"])), median_rpm2=float(np.median(res["rpm2"])),
             p95_rpm2=float(np.quantile(res["rpm2"], 0.95)),
             UA_rel_err_sd=float(np.std(res["UA_hat"] / res["UA_true"] - 1)),
             UA_rel_err_mean=float(np.mean(res["UA_hat"] / res["UA_true"] - 1)),
             mCp_rel_err_sd=float(np.std(res["mCp_hat"] / res["mCp_true"] - 1)))
    return s


def practical_margin(delta_cal, increment=DELTA_INCREMENT):
    """Pre-declared rounding rule (applied on calibration data only): round the
    calibrated margin UP to the next `increment`. A larger margin can only
    raise the pass fraction, so rounding up never weakens the calibration."""
    return math.ceil(delta_cal / increment - 1e-9) * increment


def validation_protocol(cfg, S_val, S_val2, noise_seeds, p: Params = P0, alpha=0.95, sigma=SENSOR_SIGMA_K,
                        increment=DELTA_INCREMENT):
    """Validate the policy once on S_val. If (and only if) it fails the one-sided
    lower-bound rule, apply the pre-declared modification (margin + increment)
    and validate the modified policy on the fresh S_val2. Returns the final
    configuration, the simulation result used for reporting and the audit trail."""
    def evaluate(S, delta, seed):
        res = run_policy(S, cfg["rpm1"], cfg["m_cmd"], cfg["t_obs"], delta, p, sigma=sigma, noise_seed=seed)
        s = policy_summary(res, p)
        s["accepted"] = bool(s["pass_lb95_one_sided"] >= alpha)
        return res, s

    res, s = evaluate(S_val, cfg["delta"], noise_seeds["VAL"])
    trail = [dict(stage="VAL", dataset=S_val["name"], delta=cfg["delta"], **s)]
    final = dict(cfg)
    if not s["accepted"]:
        final["delta"] = cfg["delta"] + increment
        res, s = evaluate(S_val2, final["delta"], noise_seeds["VAL2"])
        trail.append(dict(stage="VAL2 (after pre-declared modification)", dataset=S_val2["name"], delta=final["delta"], **s))
    return final, res, s, trail
