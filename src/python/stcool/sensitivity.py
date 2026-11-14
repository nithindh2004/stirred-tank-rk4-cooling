"""Local (one-at-a-time) and global (Sobol) sensitivity of T_L(t_f)."""
from __future__ import annotations

import numpy as np

from .model import P0, Params, final_TL_exact
from .uncertainty import actual_flow, SOBOL_SEED


def _flow(m, p):
    """Perturbed flow, saturated at the pump limit (an increase above m_max is not physically available)."""
    return min(p.m_max, m)

SOBOL_INPUTS = ["U_ref", "Jacket inventory M_J", "Coolant-flow error", "Initial temperature T_0"]


def local_sensitivity(rpm, m, p: Params = P0):
    """Two perturbation sets around an operating point (actual flow = m).

    'v1 steps' reproduces the original study (+/-20 % or +/-1 K).
    'one sd' moves each uncertain input by one standard deviation of its
    assumed distribution, so swings are comparable as importance measures.
    Flow perturbations are saturated at the pump limit, so at a point that
    already runs at m_max the upward flow step has no effect.
    """
    base = float(final_TL_exact(rpm, m, p))
    rows = []
    def add(scale, name, step, lo, hi):
        rows.append(dict(scale=scale, parameter=name, perturbation=step, T_low_K=float(lo), T_high_K=float(hi),
                         delta_low_K=float(lo - base), delta_high_K=float(hi - base), swing_K=float(abs(hi - lo))))
    add("v1 steps", "Overall UA scale (U_ref)", "+/-20 %", final_TL_exact(rpm, m, p, Uref=0.8 * p.Uref), final_TL_exact(rpm, m, p, Uref=1.2 * p.Uref))
    add("v1 steps", "Coolant flow", "+/-20 %", final_TL_exact(rpm, 0.8 * m, p), final_TL_exact(rpm, _flow(1.2 * m, p), p))
    add("v1 steps", "Agitator speed", "+/-20 %", final_TL_exact(0.8 * rpm, m, p), final_TL_exact(1.2 * rpm, m, p))
    add("v1 steps", "Jacket inventory M_J", "+/-20 %", final_TL_exact(rpm, m, p, MJ=0.8 * p.MJ), final_TL_exact(rpm, m, p, MJ=1.2 * p.MJ))
    add("v1 steps", "Initial temperature T_0", "+/-1 K", final_TL_exact(rpm, m, p, T0=p.T0 - 1), final_TL_exact(rpm, m, p, T0=p.T0 + 1))
    add("one sd", "U_ref", "+/-15 W/(m2 K)", final_TL_exact(rpm, m, p, Uref=p.Uref - 15), final_TL_exact(rpm, m, p, Uref=p.Uref + 15))
    add("one sd", "Coolant flow", "+/-5 %", final_TL_exact(rpm, 0.95 * m, p), final_TL_exact(rpm, _flow(1.05 * m, p), p))
    add("one sd", "Jacket inventory M_J", "+/-1 kg", final_TL_exact(rpm, m, p, MJ=p.MJ - 1), final_TL_exact(rpm, m, p, MJ=p.MJ + 1))
    add("one sd", "Initial temperature T_0", "+/-1 K", final_TL_exact(rpm, m, p, T0=p.T0 - 1), final_TL_exact(rpm, m, p, T0=p.T0 + 1))
    return base, rows


def _model_from_unit_normals(X, rpm, m_cmd, p):
    Uref = np.maximum(1.0, p.Uref + 0.10 * p.Uref * X[:, 0])
    MJ = np.maximum(1.0, p.MJ + 0.10 * p.MJ * X[:, 1])
    md = actual_flow(m_cmd, X[:, 2], p)
    T0 = p.T0 + X[:, 3]
    return final_TL_exact(rpm, md, p, Uref=Uref, MJ=MJ, T0=T0)


def sobol_indices(rpm, m_cmd, p: Params = P0, N=2 ** 15, seed=SOBOL_SEED, n_boot=500):
    """First-order (Saltelli 2010) and total-effect (Jansen 1999) indices.

    Inputs are the four independent standard-normal drivers of the uncertainty
    model, so the indices are conditional on those assumed distributions.
    Outputs are centred before estimation for numerical stability.
    """
    g = np.random.default_rng(seed)
    k = 4
    A = g.standard_normal((N, k)); B = g.standard_normal((N, k))
    fA = _model_from_unit_normals(A, rpm, m_cmd, p)
    fB = _model_from_unit_normals(B, rpm, m_cmd, p)
    fAB = []
    for i in range(k):
        ABi = A.copy(); ABi[:, i] = B[:, i]
        fAB.append(_model_from_unit_normals(ABi, rpm, m_cmd, p))
    fAB = np.array(fAB)
    f0 = np.mean(np.r_[fA, fB])
    fA, fB, fAB = fA - f0, fB - f0, fAB - f0

    def estimate(idx):
        a, b, ab = fA[idx], fB[idx], fAB[:, idx]
        V = np.var(np.r_[a, b], ddof=1)
        S1 = np.array([np.mean(b * (ab[i] - a)) / V for i in range(k)])
        ST = np.array([0.5 * np.mean((a - ab[i]) ** 2) / V for i in range(k)])
        return S1, ST

    S1, ST = estimate(np.arange(N))
    boots = [estimate(g.integers(0, N, N)) for _ in range(n_boot)]
    S1b = np.array([b[0] for b in boots]); STb = np.array([b[1] for b in boots])
    return dict(S1=S1, ST=ST, S1_lo=np.quantile(S1b, .025, 0), S1_hi=np.quantile(S1b, .975, 0),
                ST_lo=np.quantile(STb, .025, 0), ST_hi=np.quantile(STb, .975, 0), N=N, n_runs=N * (k + 2))
