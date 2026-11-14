"""Deterministic simulation studies: trajectories, verification, energy balance,
reduced-order model and the jacket flow-pattern (tanks-in-series) check."""
from __future__ import annotations

import numpy as np
from scipy.integrate import simpson

from .model import (P0, Params, rk4, trajectory_exact, original_model_TL, time_constants, U_coef,
                    UA_effective, UA_effective_plugflow, quasi_steady_TL, rk4_tanks_in_series, eigenvalues)


def time_to_target(t, TL, target):
    """Linear interpolation of the first crossing of `target` (nan if never reached)."""
    if TL[-1] > target:
        return np.nan
    i = int(np.argmax(TL <= target))
    if i == 0:
        return 0.0
    return float(t[i - 1] + (target - TL[i - 1]) * (t[i] - t[i - 1]) / (TL[i] - TL[i - 1]))


def convergence_study(rpm, m, p: Params = P0, dts=(60, 30, 20, 10, 5, 2, 1, 0.5, 0.1)):
    rows = []
    for dt in dts:
        t, Y = rk4(rpm, m, dt, p)
        eL, eJ = trajectory_exact(rpm, m, t, p)
        rows.append(dict(dt_s=dt, T_liquid_tf_K=Y[-1, 0], T_jacket_tf_K=Y[-1, 1],
                         max_traj_error_K=float(np.abs(np.c_[Y[:, 0] - eL, Y[:, 1] - eJ]).max()),
                         final_TL_error_K=float(abs(Y[-1, 0] - eL[-1]))))
    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        b["observed_order"] = float(np.log(a["max_traj_error_K"] / b["max_traj_error_K"]) / np.log(a["dt_s"] / b["dt_s"]))
    rows[0]["observed_order"] = np.nan
    sel = [r for r in rows if 2 <= r["dt_s"] <= 30]
    order = float(np.polyfit(np.log([r["dt_s"] for r in sel]), np.log([r["max_traj_error_K"] for r in sel]), 1)[0])
    return rows, order


def rk4_stability_limit(rpm, m, p: Params = P0):
    """Largest stable RK4 step for the fast (jacket) mode: |h lambda| <= 2.785 on the real axis."""
    _, l_fast = eigenvalues(U_coef(rpm, m, p) * p.A, m, p.MJ, p)
    return float(2.785 / abs(l_fast))


def original_model_rk4_error(p: Params = P0, dt=10.0):
    k = p.Uref * p.A / (p.ML * p.CpL)
    t = np.arange(0.0, p.tf + 1e-9, dt); T = np.empty_like(t); T[0] = p.T0
    f = lambda T: -k * (T - p.Tin)
    for i in range(len(t) - 1):
        k1 = f(T[i]); k2 = f(T[i] + dt / 2 * k1); k3 = f(T[i] + dt / 2 * k2); k4 = f(T[i] + dt * k3)
        T[i + 1] = T[i] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return float(np.abs(T - original_model_TL(t, p)).max())


def energy_balance(rpm, m, p: Params = P0, dt=1.0):
    """Liquid enthalpy loss vs heat carried out by the coolant + heat stored in the jacket."""
    t, Y = rk4(rpm, m, dt, p)
    Q_liq = p.ML * p.CpL * (p.T0 - Y[-1, 0])
    q_out = m * p.CpC * (Y[:, 1] - p.Tin)
    Q_trap = float(np.trapezoid(q_out, t)); Q_simp = float(simpson(q_out, x=t))
    Q_store = p.MJ * p.CpC * (Y[-1, 1] - p.Tin)
    return dict(rpm=rpm, mdot=m, Q_liquid_MJ=Q_liq / 1e6, Q_coolant_out_MJ=Q_simp / 1e6, Q_jacket_stored_MJ=Q_store / 1e6,
                residual_trapezoid_J=Q_liq - Q_trap - Q_store, residual_simpson_J=Q_liq - Q_simp - Q_store,
                relative_residual_simpson=abs(Q_liq - Q_simp - Q_store) / Q_liq)


def reduced_order_comparison(points, p: Params = P0):
    rows = []
    for name, rpm, m in points:
        t, Y = rk4(rpm, m, 1.0, p)
        G = float(UA_effective(rpm, m, p))
        Tqs = quasi_steady_TL(G, t, p)
        tl, tj = time_constants(rpm, m, p)
        rows.append(dict(case=name, rpm=rpm, mdot=m, UA_W_K=float(U_coef(rpm, m, p) * p.A), mCp_W_K=m * p.CpC,
                         UA_eff_W_K=G, TL_full_K=Y[-1, 0], TL_reduced_K=float(Tqs[-1]),
                         error_K=float(Tqs[-1] - Y[-1, 0]), max_traj_error_K=float(np.abs(Tqs - Y[:, 0]).max()),
                         tau_jacket_over_tf=tj / p.tf))
    return rows


def jacket_flow_pattern(points, p: Params = P0, cells=(1, 2, 5, 20), dt=1.0):
    """Base well-mixed jacket (1 cell) vs tanks-in-series approaching plug flow."""
    rows = []
    for name, rpm, m in points:
        for n in cells:
            rows.append(dict(case=name, rpm=rpm, mdot=m, n_cells=n,
                             TL_tf_K=float(rk4_tanks_in_series(rpm, m, dt, n, p)),
                             UA_eff_mixed=float(UA_effective(rpm, m, p)), UA_eff_plug=float(UA_effective_plugflow(rpm, m, p))))
    return rows
