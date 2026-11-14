"""Deterministic and chance-constrained optimisation of the fixed operating point.

Decision variables : agitator speed N (rpm, continuous) and commanded coolant
                     flow m_cmd (kg/s).
Bounds             : 100 <= N <= 600 rpm; 0.05 <= m_cmd <= m_max = 0.50 kg/s.
Objective          : auxiliary electricity per batch (agitator + pump), kWh.
                     Deterministic problem: actual flow = commanded flow.
                     Chance-constrained problem: expected energy over the
                     uncertain actual flow.
Constraint         : T_L(1800 s) <= 355 K (deterministic), or
                     P[T_L(1800 s) <= 355 K] >= alpha (chance constraint),
                     enforced on the optimisation dataset through the
                     one-sided 95 % Wilson lower bound (sample-average
                     approximation with a confidence margin).
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar

from .model import P0, Params, final_TL_exact, batch_energy_kWh
from .uncertainty import (mc_final_TL, pass_fraction, expected_energy_kWh, required_pass_fraction,
                          summarise, mc_energy_kWh)


# ---------------------------------------------------------------- deterministic
def grid_search(p: Params = P0, rpm_step=5.0, m_step=0.005):
    rv = np.arange(p.rpm_min, p.rpm_max + 1e-9, rpm_step)
    mv = np.round(np.arange(p.m_min, p.m_max + 1e-9, m_step), 6)
    R, M = np.meshgrid(rv, mv, indexing="ij")
    T = final_TL_exact(R, M, p)
    E = batch_energy_kWh(R, M, p)
    Ef = np.where(T <= p.Ttarget, E, np.inf)
    i = np.unravel_index(np.argmin(Ef), Ef.shape)
    return dict(rpm=float(R[i]), m=float(M[i]), T=float(T[i]), E=float(E[i]), R=R, M=M, T_grid=T, E_grid=E)


def boundary_rpm(m, p: Params = P0):
    """Speed at which T_L(t_f) = T_target for actual flow m (nan if unreachable)."""
    g = lambda r: final_TL_exact(r, m, p) - p.Ttarget
    if g(p.rpm_max) > 0:
        return np.nan
    if g(p.rpm_min) <= 0:
        return p.rpm_min
    return brentq(g, p.rpm_min, p.rpm_max, xtol=1e-10)


def nominal_optimum(p: Params = P0):
    """Continuous deterministic optimum: 1-D minimisation of energy along the active constraint."""
    def E_of(m):
        r = boundary_rpm(m, p)
        return 1e9 if np.isnan(r) else float(batch_energy_kWh(r, m, p))
    ms = np.linspace(p.m_min, p.m_max, 181)
    Es = np.array([E_of(m) for m in ms])
    j = int(np.argmin(Es))
    lo, hi = ms[max(j - 1, 0)], ms[min(j + 1, len(ms) - 1)]
    res = minimize_scalar(E_of, bounds=(lo, hi), method="bounded", options={"xatol": 1e-10})
    m = float(res.x) if res.fun <= Es[j] else float(ms[j])
    r = boundary_rpm(m, p)
    # sensitivity of optimal energy to the flow bound (finite difference along the boundary)
    h = 1e-4
    dE_dm = (E_of(m) - E_of(m - h)) / h
    return dict(rpm=float(r), m=m, E=float(batch_energy_kWh(r, m, p)), T=float(final_TL_exact(r, m, p)),
                at_flow_bound=bool(abs(m - p.m_max) < 1e-6), dE_dm_kWh_per_kgps=float(dE_dm))


# ---------------------------------------------------------------- chance-constrained
def min_rpm_meeting(m_cmd, p_req, S, p: Params = P0, tol=1e-3):
    """Lowest speed whose sample pass fraction >= p_req (pass fraction is monotone in speed)."""
    if pass_fraction(p.rpm_max, m_cmd, S, p) < p_req:
        return np.nan
    if pass_fraction(p.rpm_min, m_cmd, S, p) >= p_req:
        return p.rpm_min
    lo, hi = p.rpm_min, p.rpm_max
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if pass_fraction(mid, m_cmd, S, p) >= p_req:
            hi = mid
        else:
            lo = mid
    return hi


def chance_constrained_optimum(alpha, S, p: Params = P0, m_step=0.005, m_max=None, use_lower_bound=True):
    """Minimise expected energy s.t. the pass-rate requirement on dataset S.

    With use_lower_bound=True the requirement is that the one-sided 95 % Wilson
    lower bound of the sample pass fraction is >= alpha.
    """
    m_max = p.m_max if m_max is None else m_max
    n = len(S["T0"])
    p_req = required_pass_fraction(alpha, n) if use_lower_bound else alpha
    best, rows = None, []
    for m in np.round(np.arange(p.m_min, m_max + 1e-9, m_step), 6):
        r = min_rpm_meeting(m, p_req, S, p)
        if np.isnan(r):
            continue
        E = expected_energy_kWh(r, m, S, p)
        rows.append((m, r, E))
        if best is None or E < best["E"]:
            best = dict(rpm=float(r), m=float(m), E=E)
    if best is not None:
        best["p_req"] = p_req
        best["in_sample_pass"] = pass_fraction(best["rpm"], best["m"], S, p)
    return best, np.array(rows)


def practical_setpoint(design, alpha, S_opt, p: Params = P0, rpm_increment=5.0):
    """Pre-declared rounding rule, applied using the optimisation data only:
    round the speed UP to the next `rpm_increment`; keep the commanded flow.
    Rounding up can only raise the pass fraction (monotonicity), which is
    re-checked on S_opt."""
    rpm = math.ceil(design["rpm"] / rpm_increment - 1e-9) * rpm_increment
    p_req = required_pass_fraction(alpha, len(S_opt["T0"]))
    while pass_fraction(rpm, design["m"], S_opt, p) < p_req:        # never triggers in practice
        rpm += rpm_increment
    return dict(rpm=float(rpm), m=design["m"])


def validate(rpm, m_cmd, S, p: Params = P0, alpha=0.95):
    x = mc_final_TL(rpm, m_cmd, S, p)
    E = mc_energy_kWh(rpm, m_cmd, S, p)
    s = summarise(x, p.Ttarget, E)
    s["accepted"] = bool(s["pass_lb95_one_sided"] >= alpha)
    s["dataset"] = S["name"]
    return s


def validation_protocol(setpoint, S_val, S_val2, p: Params = P0, alpha=0.95, rpm_increment=5.0):
    """Validate once on VAL. If (and only if) the design fails acceptance, apply the
    pre-declared modification (+rpm_increment) and validate the modified design on
    the fresh VAL2 set. Returns the full audit trail."""
    trail = [dict(stage="VAL", rpm=setpoint["rpm"], m=setpoint["m"], **validate(setpoint["rpm"], setpoint["m"], S_val, p, alpha))]
    final = dict(setpoint)
    if not trail[0]["accepted"]:
        final = dict(rpm=setpoint["rpm"] + rpm_increment, m=setpoint["m"])
        trail.append(dict(stage="VAL2 (after modification)", rpm=final["rpm"], m=final["m"],
                          **validate(final["rpm"], final["m"], S_val2, p, alpha)))
    return final, trail
