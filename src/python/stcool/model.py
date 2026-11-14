"""Physical model of the jacketed stirred-tank cooling process.

States: T_L (bulk liquid) and T_J (well-mixed jacket coolant), both in K.

    M_L Cp_L dT_L/dt = -U A (T_L - T_J)                                   (liquid)
    M_J Cp_C dT_J/dt =  U A (T_L - T_J) + mdot Cp_C (T_in - T_J)          (jacket)
    U = U_ref (N / N_ref)^eN (mdot / mdot_ref)^eM                         (closure)

`mdot` is always the *actual* coolant flow. How a commanded flow becomes an
actual flow (actuator limit, control error) lives in `uncertainty.py`.

Two solvers are provided:
  * classical RK4 -- the principal method for every deterministic trajectory;
  * the exact matrix-exponential solution -- an independent verification
    reference, also used to evaluate T_L(t_f) cheaply for large Monte Carlo
    sample sets (the system is linear with constant coefficients over any
    interval in which N and mdot are constant).
"""
from __future__ import annotations

from dataclasses import dataclass, replace, asdict
import numpy as np


@dataclass(frozen=True)
class Params:
    # liquid charge (carried over from the original ESC113M model)
    rhoL: float = 1000.0      # kg/m3
    V: float = 0.50           # m3
    CpL: float = 4180.0       # J/(kg K)
    T0: float = 370.0         # K
    A: float = 2.0            # m2
    Uref: float = 150.0       # W/(m2 K) at the reference operating point
    Tin: float = 300.0        # K, coolant inlet (the original fixed coolant temperature)
    tf: float = 1800.0        # s, batch horizon
    # coolant and jacket (case-study assumptions)
    CpC: float = 4180.0       # J/(kg K)
    rhoC: float = 1000.0      # kg/m3
    MJ: float = 10.0          # kg, effective well-mixed jacket inventory
    # heat-transfer closure (assumed exponents)
    Nref: float = 300.0       # rpm
    mref: float = 0.20        # kg/s
    eN: float = 0.25
    eM: float = 0.33
    # design target (illustrative)
    Ttarget: float = 355.0    # K at t_f
    # auxiliary electrical power model (assumed screening values)
    D: float = 0.30           # impeller diameter, m
    Np: float = 5.0           # impeller power number (turbulent, Rushton-type)
    eta_ag: float = 0.75      # agitator drive efficiency
    dp_ref: float = 50_000.0  # Pa, coolant-loop pressure drop at mref
    eta_p: float = 0.65       # pump efficiency
    # operating bounds
    rpm_min: float = 100.0
    rpm_max: float = 600.0
    m_min: float = 0.05       # kg/s, lowest commanded flow considered
    m_max: float = 0.50       # kg/s, physical maximum deliverable coolant flow

    @property
    def ML(self) -> float:
        return self.rhoL * self.V

    def with_(self, **kw) -> "Params":
        return replace(self, **kw)

    def as_dict(self) -> dict:
        d = asdict(self)
        d["ML"] = self.ML
        return d


P0 = Params()

# Film-coefficient exponents used only to bound the closure exponents in the
# model-form stress test (textbook correlations: agitated-side Nu ~ Re^(2/3),
# turbulent channel flow Nu ~ Re^0.8). If U is a series of inside-film,
# coolant-film and wall/fouling resistances with shares f_in, f_co, f_w, then
# eN = (2/3) f_in and eM = 0.8 f_co locally.
FILM_EXP_INSIDE, FILM_EXP_COOLANT = 2.0 / 3.0, 0.8


# ---------------------------------------------------------------- closures
def U_coef(rpm, mdot, p: Params = P0, Uref=None):
    """Overall heat-transfer coefficient, W/(m2 K), at actual flow `mdot`."""
    Uref = p.Uref if Uref is None else Uref
    return Uref * (np.asarray(rpm, float) / p.Nref) ** p.eN * (np.asarray(mdot, float) / p.mref) ** p.eM


def agitator_power(rpm, p: Params = P0):
    """Agitator electrical power, W:  N_P rho n^3 D^5 / eta (n in rev/s)."""
    n = np.asarray(rpm, float) / 60.0
    return p.Np * p.rhoL * n ** 3 * p.D ** 5 / p.eta_ag


def pump_power(mdot, p: Params = P0):
    """Pump electrical power, W, at actual flow:  dP * Q / eta, dP = dP_ref (m/m_ref)^2, Q = m/rho."""
    mdot = np.asarray(mdot, float)
    dp = p.dp_ref * (mdot / p.mref) ** 2          # Pa
    return dp * (mdot / p.rhoC) / p.eta_p         # Pa * m3/s = W


def loop_pressure_drop(mdot, p: Params = P0):
    return p.dp_ref * (np.asarray(mdot, float) / p.mref) ** 2


def batch_energy_kWh(rpm, mdot, p: Params = P0, duration=None):
    """Auxiliary electricity (agitator + pump) for constant speed and actual flow."""
    duration = p.tf if duration is None else duration
    return (agitator_power(rpm, p) + pump_power(mdot, p)) * duration / 3.6e6


# ---------------------------------------------------------------- exact solution
def _coefficients(UA, mdot, MJ, p: Params):
    a = UA / (p.ML * p.CpL)
    b = UA / (MJ * p.CpC)
    f = mdot / MJ
    return a, b, f


def eigenvalues(UA, mdot, MJ, p: Params = P0):
    """Eigenvalues (slow l1, fast l2) of M = [[-a, a], [b, -(b+f)]]; always real, negative, distinct."""
    a, b, f = _coefficients(UA, mdot, MJ, p)
    tr = -(a + b + f)
    disc = np.sqrt(tr * tr - 4.0 * a * f)
    return 0.5 * (tr + disc), 0.5 * (tr - disc)


def propagate_exact(TL0, TJ0, UA, mdot, MJ, t, p: Params = P0):
    """Exact state after time t from (TL0, TJ0) for constant UA, mdot, MJ.

    With x = T - T_in and x' = M x,
        e^{Mt} = [e^{l1 t}(M - l2 I) - e^{l2 t}(M - l1 I)] / (l1 - l2).
    All arguments broadcast (NumPy).
    """
    a, b, f = _coefficients(UA, mdot, MJ, p)
    l1, l2 = eigenvalues(UA, mdot, MJ, p)
    x1, x2 = TL0 - p.Tin, TJ0 - p.Tin
    e1, e2 = np.exp(l1 * t), np.exp(l2 * t)
    d = l1 - l2
    TL = p.Tin + (e1 * ((-a - l2) * x1 + a * x2) - e2 * ((-a - l1) * x1 + a * x2)) / d
    TJ = p.Tin + (e1 * (b * x1 + (-(b + f) - l2) * x2) - e2 * (b * x1 + (-(b + f) - l1) * x2)) / d
    return TL, TJ


def final_TL_exact(rpm, mdot, p: Params = P0, Uref=None, MJ=None, T0=None):
    """Exact T_L(t_f) at constant speed and actual flow; uncertain inputs may be arrays."""
    UA = U_coef(rpm, mdot, p, Uref) * p.A
    MJ = p.MJ if MJ is None else MJ
    T0 = p.T0 if T0 is None else T0
    return propagate_exact(T0, p.Tin, UA, mdot, MJ, p.tf, p)[0]


def trajectory_exact(rpm, mdot, t, p: Params = P0):
    UA = U_coef(rpm, mdot, p) * p.A
    return propagate_exact(p.T0, p.Tin, UA, mdot, p.MJ, np.asarray(t, float), p)


def time_constants(rpm, mdot, p: Params = P0):
    """(tank, jacket) time constants -1/l, s."""
    l1, l2 = eigenvalues(U_coef(rpm, mdot, p) * p.A, mdot, p.MJ, p)
    return -1.0 / l1, -1.0 / l2


# ---------------------------------------------------------------- RK4
def rhs(y, UA, mdot, MJ, p: Params = P0):
    TL, TJ = y[0], y[1]
    Q = UA * (TL - TJ)                                   # W, liquid -> jacket
    return np.array([-Q / (p.ML * p.CpL),
                     (Q + mdot * p.CpC * (p.Tin - TJ)) / (MJ * p.CpC)])


def rk4_step(y, h, UA, mdot, MJ, p: Params = P0):
    k1 = rhs(y, UA, mdot, MJ, p)
    k2 = rhs(y + 0.5 * h * k1, UA, mdot, MJ, p)
    k3 = rhs(y + 0.5 * h * k2, UA, mdot, MJ, p)
    k4 = rhs(y + h * k3, UA, mdot, MJ, p)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def rk4(rpm, mdot, dt, p: Params = P0, y0=None, t_end=None):
    """Classical RK4 trajectory at constant speed and actual flow. Returns (t, Y[:, (T_L, T_J)])."""
    t_end = p.tf if t_end is None else t_end
    n = int(round(t_end / dt))
    h = t_end / n
    UA = float(U_coef(rpm, mdot, p) * p.A)
    Y = np.empty((n + 1, 2))
    Y[0] = (p.T0, p.Tin) if y0 is None else y0
    for k in range(n):
        Y[k + 1] = rk4_step(Y[k], h, UA, mdot, p.MJ, p)
    return np.arange(n + 1) * h, Y


def rk4_tanks_in_series(rpm, mdot, dt, n_cells, p: Params = P0, Uref=None, MJ=None, T0=None):
    """Vectorised RK4 for a jacket split into `n_cells` equal well-mixed cells in series.

    n_cells = 1 recovers the base model; large n_cells approaches a plug-flow
    jacket. The wall area, U and jacket inventory are divided equally among the
    cells; the liquid sees the sum of the cell heat flows. Returns T_L(t_f).
    """
    Uref = p.Uref if Uref is None else Uref
    MJ = p.MJ if MJ is None else MJ
    T0 = p.T0 if T0 is None else T0
    UA = np.asarray(U_coef(rpm, mdot, p, Uref) * p.A, float)
    shape = np.broadcast(UA, np.asarray(MJ), np.asarray(T0), np.asarray(mdot)).shape
    UAc = np.broadcast_to(UA, shape) / n_cells
    Cc = np.broadcast_to(np.asarray(MJ, float), shape) / n_cells * p.CpC
    mC = np.broadcast_to(np.asarray(mdot, float), shape) * p.CpC
    TL = np.broadcast_to(np.asarray(T0, float), shape).copy()
    TJ = np.full(shape + (n_cells,), p.Tin)

    def f(TL, TJ):
        q = UAc[..., None] * (TL[..., None] - TJ)                  # W per cell
        upstream = np.concatenate([np.full(shape + (1,), p.Tin), TJ[..., :-1]], axis=-1)
        dTJ = (q + mC[..., None] * (upstream - TJ)) / Cc[..., None]
        return -q.sum(-1) / (p.ML * p.CpL), dTJ

    n = int(round(p.tf / dt)); h = p.tf / n
    for _ in range(n):
        a1, b1 = f(TL, TJ)
        a2, b2 = f(TL + 0.5 * h * a1, TJ + 0.5 * h * b1)
        a3, b3 = f(TL + 0.5 * h * a2, TJ + 0.5 * h * b2)
        a4, b4 = f(TL + h * a3, TJ + h * b3)
        TL = TL + (h / 6) * (a1 + 2 * a2 + 2 * a3 + a4)
        TJ = TJ + (h / 6) * (b1 + 2 * b2 + 2 * b3 + b4)
    return TL


# ---------------------------------------------------------------- reduced-order model
def UA_effective(rpm, mdot, p: Params = P0):
    """Quasi-steady series-resistance conductance (1/UA + 1/(mdot Cp))^-1, W/K."""
    UA = U_coef(rpm, mdot, p) * p.A
    C = np.asarray(mdot, float) * p.CpC
    return UA * C / (UA + C)


def UA_effective_plugflow(rpm, mdot, p: Params = P0):
    """Quasi-steady conductance for a plug-flow jacket: C (1 - exp(-UA/C))."""
    UA = U_coef(rpm, mdot, p) * p.A
    C = np.asarray(mdot, float) * p.CpC
    return C * (1.0 - np.exp(-UA / C))


def quasi_steady_TL(UAeff, t, p: Params = P0):
    return p.Tin + (p.T0 - p.Tin) * np.exp(-UAeff * np.asarray(t, float) / (p.ML * p.CpL))


def UA_eff_required(p: Params = P0):
    """Effective conductance needed to reach T_target at t_f (quasi-steady)."""
    return p.ML * p.CpL * np.log((p.T0 - p.Tin) / (p.Ttarget - p.Tin)) / p.tf


def analytic_boundary_rpm(mdot, p: Params = P0):
    """Closed-form feasibility boundary N(mdot) from 1/UA + 1/(mdot Cp) = 1/UA_eff,req."""
    G = UA_eff_required(p)
    C = np.asarray(mdot, float) * p.CpC
    with np.errstate(divide="ignore", invalid="ignore"):
        UA = np.where(C > G, G * C / (C - G), np.nan)
        return p.Nref * (UA / p.A / (p.Uref * (np.asarray(mdot) / p.mref) ** p.eM)) ** (1.0 / p.eN)


def original_model_TL(t, p: Params = P0):
    """Original ESC113M single-state model (coolant fixed at T_in), analytical solution."""
    k = p.Uref * p.A / (p.ML * p.CpL)
    return p.Tin + (p.T0 - p.Tin) * np.exp(-k * np.asarray(t, float))
