"""End-to-end V3 workflow: every dataset, figure, output table and the LaTeX
value file used by both reports is produced here, from fixed seeds."""
from __future__ import annotations

import json, math, platform, time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy, matplotlib

from .model import (P0, FILM_EXP_INSIDE, FILM_EXP_COOLANT, U_coef, agitator_power, pump_power, loop_pressure_drop,
                    batch_energy_kWh, final_TL_exact, trajectory_exact, rk4, time_constants, UA_effective, UA_effective_plugflow,
                    UA_eff_required, original_model_TL, rk4_tanks_in_series)
from .uncertainty import (DATASETS, NOISE_SEEDS, SOBOL_SEED, UNCERTAINTY_TABLE, FLOW_ERROR_REL_SD, dataset, actual_flow,
                          mc_final_TL, mc_energy_kWh, pass_fraction, summarise, wilson_interval, required_pass_fraction)
from . import simulation as sim
from . import optimization as opt
from . import sensitivity as sens
from . import adaptive_policy as ap
from .plotting import apply_style
from .figures import make_figures
from .report_values import write_latex

ALPHA = 0.95
RPM_INCREMENT = 5.0
FRONTIER_ALPHAS = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.975, 0.99]
ADAPTIVE_GRID = dict(rpm1=[120.0, 160.0, 200.0, 240.0], t_obs=[300.0, 600.0, 900.0])


class Run:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.data = self.root / "data"; self.fig = self.root / "figures"; self.out = self.root / "outputs"
        self.gen = self.root / "reports" / "build" / "generated"
        for d in (self.data, self.fig, self.out, self.gen):
            d.mkdir(parents=True, exist_ok=True)
        self.macros: dict[str, str] = {}
        self.summary: list[dict] = []
        self.validations: list[dict] = []
        self.t0 = time.time()

    def csv(self, rows, rel):
        df = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
        path = self.data / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False, float_format="%.10g")
        return df

    def m(self, name, value, fmt="{:.2f}"):
        assert name.isalpha(), name
        self.macros[name] = fmt.format(value) if not isinstance(value, str) else value

    def key(self, metric, value, unit, source):
        self.summary.append(dict(metric=metric, value=value, unit=unit, source=source))

    def log(self, msg):
        print(f"[{time.time() - self.t0:6.1f}s] {msg}", flush=True)


def pct(x, d=1):
    return f"{100 * x:.{d}f}"


def main(root: Path):
    R = Run(root)
    p = P0
    apply_style()

    # ================================================================ 0. inputs
    R.log("0. inputs")
    params = [
        ("V", p.V, "m3", "liquid volume", "carried over from ESC113M"),
        ("rho_L", p.rhoL, "kg/m3", "liquid density", "carried over from ESC113M"),
        ("Cp_L", p.CpL, "J/(kg K)", "liquid heat capacity", "carried over from ESC113M"),
        ("T_0", p.T0, "K", "initial liquid temperature", "carried over from ESC113M"),
        ("T_in", p.Tin, "K", "coolant inlet temperature (original fixed coolant)", "carried over from ESC113M"),
        ("A", p.A, "m2", "heat-transfer area", "carried over from ESC113M"),
        ("U_ref", p.Uref, "W/(m2 K)", "overall coefficient at the reference point", "carried over; needs calibration"),
        ("t_f", p.tf, "s", "batch horizon", "carried over from ESC113M"),
        ("Cp_C", p.CpC, "J/(kg K)", "coolant heat capacity (water)", "literature property value"),
        ("rho_C", p.rhoC, "kg/m3", "coolant density (water)", "literature property value"),
        ("N_P", p.Np, "-", "impeller power number (turbulent Rushton-type turbine)", "literature-supported; vendor calibration"),
        ("M_J", p.MJ, "kg", "effective well-mixed jacket inventory", "case-study assumption; needs geometry"),
        ("N_ref", p.Nref, "rpm", "reference agitator speed", "case-study assumption"),
        ("mdot_ref", p.mref, "kg/s", "reference coolant flow", "case-study assumption"),
        ("e_N", p.eN, "-", "U exponent on speed", "case-study assumption; needs calibration"),
        ("e_M", p.eM, "-", "U exponent on coolant flow", "case-study assumption; needs calibration"),
        ("D", p.D, "m", "impeller diameter", "case-study assumption"),
        ("eta_ag", p.eta_ag, "-", "agitator drive efficiency", "case-study assumption; vendor calibration"),
        ("dP_ref", p.dp_ref, "Pa", "coolant-loop pressure drop at mdot_ref", "case-study assumption; vendor calibration"),
        ("eta_p", p.eta_p, "-", "pump efficiency", "case-study assumption; vendor calibration"),
        ("m_max", p.m_max, "kg/s", "maximum deliverable coolant flow (pump limit)", "case-study assumption (equipment limit)"),
        ("N range", f"{p.rpm_min:.0f}-{p.rpm_max:.0f}", "rpm", "agitator speed bounds", "case-study assumption"),
        ("T_target", p.Ttarget, "K", "liquid temperature required at t_f", "illustrative design constraint"),
    ]
    R.csv([dict(symbol=a, value=b, unit=c, description=d, classification=e) for a, b, c, d, e in params], "inputs/parameters.csv")
    R.csv([dict(input=a, distribution=b, mean=c, sd=d, unit=e, basis=f) for a, b, c, d, e, f in UNCERTAINTY_TABLE], "inputs/uncertainty_model.csv")
    R.csv([dict(dataset=k, seed=v["seed"], n=v["n"], noise_seed=NOISE_SEEDS.get(k, ""), purpose=v["purpose"]) for k, v in DATASETS.items()]
          + [dict(dataset="SOBOL", seed=SOBOL_SEED, n=2 ** 15 * 6, noise_seed="", purpose="Saltelli/Jansen Sobol estimation")], "inputs/datasets.csv")

    # ================================================================ 1. verification
    R.log("1. verification")
    conv, order = sim.convergence_study(300.0, 0.20, p)
    R.csv(conv, "verification/timestep_convergence.csv")
    tB, YB = rk4(300.0, 0.20, 1.0, p)
    eL, eJ = trajectory_exact(300.0, 0.20, tB, p)
    err_dt1 = float(np.abs(np.c_[YB[:, 0] - eL, YB[:, 1] - eJ]).max())
    orig_err = sim.original_model_rk4_error(p)
    R.m("ConvOrder", order); R.m("ErrDtOne", f"{err_dt1:.1e}".replace("e-0", "\\times10^{-").replace("e-", "\\times10^{-") + "}")
    R.m("OrigErr", f"{orig_err:.1e}".replace("e-", "\\times10^{-") + "}")
    stab_base = sim.rk4_stability_limit(300.0, 0.20, p)
    R.m("StabBase", stab_base, "{:.0f}")
    R.key("Original model T_L(1800 s) at 300 rpm, 0.20 kg/s", float(original_model_TL(p.tf, p)), "K", "deterministic/case_comparison.csv")
    R.key("Coupled model T_L(1800 s) at 300 rpm, 0.20 kg/s (RK4)", float(YB[-1, 0]), "K", "deterministic/case_comparison.csv")
    R.key("RK4 observed order (dt 2-30 s)", order, "-", "verification/timestep_convergence.csv")
    R.key("RK4 max trajectory error at dt = 1 s", err_dt1, "K", "verification/timestep_convergence.csv")

    # ================================================================ 2. baseline physics
    R.log("2. baseline and reduced-order model")
    T_orig = original_model_TL(tB, p)
    tl_b, tj_b = time_constants(300.0, 0.20, p)
    R.m("TOrig", T_orig[-1]); R.m("TBase", YB[-1, 0]); R.m("TJBase", YB[-1, 1]); R.m("TBaseGap", YB[-1, 0] - T_orig[-1])
    R.m("PeakTJBase", YB[:, 1].max(), "{:.1f}"); R.m("PeakTJBaseMin", tB[YB[:, 1].argmax()] / 60, "{:.1f}")
    R.m("TauTankBase", tl_b, "{:,.0f}"); R.m("TauJacketBase", tj_b, "{:.1f}")
    G_b = float(UA_effective(300.0, 0.20, p)); G_req = float(UA_eff_required(p))
    R.m("UAEffBase", G_b, "{:.1f}"); R.m("UAEffReq", G_req, "{:.1f}"); R.m("MinFlow", G_req / p.CpC, "{:.3f}")
    R.m("UAEffDrop", 100 * (1 - G_b / 300.0), "{:.0f}")

    # ================================================================ 3. deterministic optimum
    R.log("3. deterministic optimisation")
    grid = opt.grid_search(p)
    nom = opt.nominal_optimum(p)
    NOM = (nom["rpm"], nom["m"])
    R.m("NomRpm", nom["rpm"], "{:.1f}"); R.m("NomE", nom["E"], "{:.3f}"); R.m("GridRpm", grid["rpm"], "{:.0f}")
    R.m("GridM", grid["m"], "{:.2f}"); R.m("GridE", grid["E"], "{:.4f}"); R.m("NomEFour", nom["E"], "{:.4f}")
    R.csv([dict(method="grid search (5 rpm x 0.005 kg/s)", rpm=grid["rpm"], mdot=grid["m"], T_L_tf_K=grid["T"], energy_kWh=grid["E"]),
           dict(method="continuous (Brent root + bounded scalar minimisation)", rpm=nom["rpm"], mdot=nom["m"], T_L_tf_K=nom["T"],
                energy_kWh=nom["E"], at_flow_bound=nom["at_flow_bound"], dE_dmdot_kWh_per_kg_s=nom["dE_dm_kWh_per_kgps"])],
          "optimization/nominal_optimum.csv")
    R.csv(pd.DataFrame(dict(rpm=grid["R"].ravel(), mdot=grid["M"].ravel(), T_L_tf_K=grid["T_grid"].ravel(), energy_kWh=grid["E_grid"].ravel())),
          "optimization/operating_grid.csv")

    # ================================================================ 4. datasets + actuator
    R.log("4. datasets and actuator model")
    OPT, VAL, VAL2, CAL, STRESS = (dataset(k, p) for k in ("OPT", "VAL", "VAL2", "CAL", "STRESS"))
    R.m("NOpt", DATASETS["OPT"]["n"], "{:,}"); R.m("NVal", DATASETS["VAL"]["n"], "{:,}"); R.m("NCal", DATASETS["CAL"]["n"], "{:,}")
    R.m("NStress", DATASETS["STRESS"]["n"], "{:,}")
    for k in ("OPT", "VAL", "VAL2", "CAL", "STRESS"):
        R.m("Seed" + k.title().replace("2", "Two").replace("3", "Three"), DATASETS[k]["seed"], "{}")
    p_req = required_pass_fraction(ALPHA, DATASETS["OPT"]["n"])
    R.m("PReqOpt", 100 * p_req, "{:.2f}")
    act_rows = []
    for mc in (0.40, 0.45, 0.50):
        md = actual_flow(mc, VAL["z_flow"], p)
        unsat = mc * (1 + FLOW_ERROR_REL_SD * VAL["z_flow"])
        act_rows.append(dict(m_cmd=mc, frac_at_limit=float(np.mean(md >= p.m_max - 1e-12)), mean_actual=float(md.mean()),
                             sd_actual=float(md.std()), p05=float(np.quantile(md, .05)), median=float(np.median(md)),
                             max_actual=float(md.max()), v2_frac_above_limit=float(np.mean(unsat > p.m_max)), v2_max=float(unsat.max()),
                             mean_pump_W=float(pump_power(md, p).mean()), pump_W_at_command=float(pump_power(mc, p))))
    act = R.csv(act_rows, "monte_carlo/flow_actuator_distribution.csv")
    a50 = act[act.m_cmd == 0.50].iloc[0]
    R.m("FracAtCap", 100 * a50.frac_at_limit, "{:.1f}"); R.m("MeanActual", a50.mean_actual, "{:.4f}")
    R.m("VTwoAbove", 100 * a50.v2_frac_above_limit, "{:.1f}"); R.m("VTwoMax", a50.v2_max, "{:.3f}")
    R.m("PumpWCmd", a50.pump_W_at_command, "{:.1f}"); R.m("PumpWMean", a50.mean_pump_W, "{:.1f}")

    # ================================================================ 5. chance-constrained optimum + validation
    R.log("5. chance-constrained optimisation (OPT) and validation (VAL)")
    best, scan = opt.chance_constrained_optimum(ALPHA, OPT, p)
    R.csv(pd.DataFrame(scan, columns=["m_cmd", "min_rpm_meeting_requirement", "expected_energy_kWh"]), "optimization/chance_constrained_flow_scan.csv")
    sp = opt.practical_setpoint(best, ALPHA, OPT, p, RPM_INCREMENT)
    final, trail = opt.validation_protocol(sp, VAL, VAL2, p, ALPHA, RPM_INCREMENT)
    REC = (final["rpm"], final["m"])
    audit = [dict(stage="OPT: continuous chance-constrained optimum", rpm=best["rpm"], m=best["m"], dataset="OPT",
                  pass_frac=best["in_sample_pass"], expected_energy_kWh=best["E"], requirement=f"sample pass >= {p_req:.5f}"),
             dict(stage="OPT: practical setpoint (speed rounded up to 5 rpm)", rpm=sp["rpm"], m=sp["m"], dataset="OPT",
                  pass_frac=pass_fraction(sp["rpm"], sp["m"], OPT, p), expected_energy_kWh=float(mc_energy_kWh(sp["rpm"], sp["m"], OPT, p).mean()),
                  requirement=f"sample pass >= {p_req:.5f}")]
    for t in trail:
        audit.append(dict(stage=t["stage"], rpm=t["rpm"], m=t["m"], dataset=t["dataset"], pass_frac=t["pass_frac"],
                          pass_ci_lo=t["pass_ci_lo"], pass_ci_hi=t["pass_ci_hi"], pass_lb95_one_sided=t["pass_lb95_one_sided"],
                          expected_energy_kWh=t["mean_E"], accepted=t["accepted"], requirement="one-sided 95 % lower bound >= 0.95"))
    R.csv(audit, "optimization/recommended_setpoint_audit_trail.csv")
    vrec = trail[-1]
    R.m("CCRpm", best["rpm"], "{:.1f}"); R.m("CCE", best["E"], "{:.3f}"); R.m("CCPassOpt", 100 * best["in_sample_pass"], "{:.2f}")
    R.m("RecRpm", REC[0], "{:.0f}"); R.m("RecM", REC[1], "{:.2f}")
    R.m("RecE", vrec["mean_E"], "{:.3f}"); R.m("RecPass", 100 * vrec["pass_frac"], "{:.2f}")
    R.m("RecLo", 100 * vrec["pass_ci_lo"], "{:.2f}"); R.m("RecHi", 100 * vrec["pass_ci_hi"], "{:.2f}")
    R.m("RecLB", 100 * vrec["pass_lb95_one_sided"], "{:.2f}"); R.m("RecPNF", vrec["p95_T"], "{:.2f}")
    R.m("RecMeanT", vrec["mean_T"], "{:.2f}"); R.m("RecSdT", vrec["sd_T"], "{:.2f}")
    R.m("RecValStage", "accepted on VAL; VAL2 not needed" if len(trail) == 1 else "modified and validated on VAL2")
    R.m("RecDetE", float(batch_energy_kWh(*REC, p)), "{:.3f}")
    for t in trail:
        R.validations.append(dict(design="recommended fixed setting", stage=t["stage"], rpm=t["rpm"], m_cmd=t["m"], dataset=t["dataset"],
                                  n=t["n"], pass_frac=t["pass_frac"], ci_lo=t["pass_ci_lo"], ci_hi=t["pass_ci_hi"],
                                  lb95_one_sided=t["pass_lb95_one_sided"], mean_energy_kWh=t["mean_E"], accepted=t["accepted"]))
    R.key("Nominal (deterministic) optimum speed", nom["rpm"], "rpm", "optimization/nominal_optimum.csv")
    R.key("Nominal optimum energy", nom["E"], "kWh/batch", "optimization/nominal_optimum.csv")
    R.key("Recommended agitator speed", REC[0], "rpm", "optimization/recommended_setpoint_audit_trail.csv")
    R.key("Recommended commanded coolant flow", REC[1], "kg/s", "optimization/recommended_setpoint_audit_trail.csv")
    R.key("Recommended setting expected energy (VAL)", vrec["mean_E"], "kWh/batch", "optimization/recommended_setpoint_audit_trail.csv")
    R.key("Recommended setting P(pass) (VAL)", vrec["pass_frac"], "-", "optimization/recommended_setpoint_audit_trail.csv")
    R.key("Recommended setting P(pass) one-sided 95% lower bound (VAL)", vrec["pass_lb95_one_sided"], "-", "optimization/recommended_setpoint_audit_trail.csv")

    # legacy designs re-evaluated under the corrected model
    legacy = {"v1 nominal grid optimum": (125.0, 0.49), "v2/v3 nominal optimum": NOM, "v1 robust candidate": (300.0, 0.50),
              "v2 recommended": (285.0, 0.50), "v3 recommended": REC}
    leg_rows = []
    for name, (r, mc) in legacy.items():
        s = opt.validate(r, mc, VAL, p)
        leg_rows.append(dict(design=name, rpm=r, m_cmd=mc, pass_frac=s["pass_frac"], ci_lo=s["pass_ci_lo"], ci_hi=s["pass_ci_hi"],
                             lb95_one_sided=s["pass_lb95_one_sided"], p95_T=s["p95_T"], mean_energy_kWh=s["mean_E"],
                             deterministic_energy_kWh=float(batch_energy_kWh(r, mc, p)), meets_rule=s["accepted"]))
    leg = R.csv(leg_rows, "monte_carlo/designs_under_v3_model.csv")
    v2r = leg[leg.design == "v2 recommended"].iloc[0]; v1r = leg[leg.design == "v1 robust candidate"].iloc[0]
    nomv = leg[leg.design == "v2/v3 nominal optimum"].iloc[0]
    R.m("VTwoRecPass", 100 * v2r.pass_frac, "{:.2f}"); R.m("VTwoRecLB", 100 * v2r.lb95_one_sided, "{:.2f}")
    R.m("VOneRobPass", 100 * v1r.pass_frac, "{:.2f}"); R.m("VOneRobE", v1r.mean_energy_kWh, "{:.3f}")
    R.m("NomPass", 100 * nomv.pass_frac, "{:.1f}")
    R.m("RecVsVOne", 100 * (1 - vrec["mean_E"] / v1r.mean_energy_kWh), "{:.1f}")

    # robust screen around the recommendation (OPT and VAL)
    scr = []
    for r in np.arange(REC[0] - 15, REC[0] + 15.1, 5):
        so = summarise(mc_final_TL(r, REC[1], OPT, p), p.Ttarget); sv = summarise(mc_final_TL(r, REC[1], VAL, p), p.Ttarget)
        scr.append(dict(rpm=r, m_cmd=REC[1], expected_energy_kWh=float(mc_energy_kWh(r, REC[1], VAL, p).mean()),
                        opt_pass=so["pass_frac"], opt_meets_requirement=so["pass_frac"] >= p_req, val_pass=sv["pass_frac"],
                        val_ci_lo=sv["pass_ci_lo"], val_ci_hi=sv["pass_ci_hi"], val_lb95=sv["pass_lb95_one_sided"]))
    R.csv(scr, "optimization/speed_screen_near_recommendation.csv")

    # ================================================================ 6. frontier
    R.log("6. reliability-energy frontier")
    fr = []
    for a in FRONTIER_ALPHAS:
        b, _ = opt.chance_constrained_optimum(a, OPT, p)
        s = opt.validate(b["rpm"], b["m"], VAL, p, alpha=a)
        fr.append(dict(alpha=a, rpm=b["rpm"], m_cmd=b["m"], expected_energy_kWh=s["mean_E"], opt_requirement=b["p_req"],
                       opt_pass=b["in_sample_pass"], val_pass=s["pass_frac"], val_ci_lo=s["pass_ci_lo"], val_ci_hi=s["pass_ci_hi"],
                       val_lb95=s["pass_lb95_one_sided"], val_meets_alpha_lb=s["accepted"]))
        R.validations.append(dict(design=f"frontier alpha={a}", stage="VAL", rpm=b["rpm"], m_cmd=b["m"], dataset="VAL", n=s["n"],
                                  pass_frac=s["pass_frac"], ci_lo=s["pass_ci_lo"], ci_hi=s["pass_ci_hi"], lb95_one_sided=s["pass_lb95_one_sided"],
                                  mean_energy_kWh=s["mean_E"], accepted=s["accepted"]))
    fr = R.csv(fr, "optimization/reliability_energy_frontier.csv")
    f50 = fr[fr.alpha == 0.5].iloc[0]; f90 = fr[fr.alpha == 0.9].iloc[0]; f99 = fr[fr.alpha == 0.99].iloc[0]
    R.m("FrFiftyE", f50.expected_energy_kWh, "{:.3f}"); R.m("FrNinetyE", f90.expected_energy_kWh, "{:.3f}")
    R.m("FrNinetyNineE", f99.expected_energy_kWh, "{:.2f}"); R.m("FrNinetyNineRpm", f99.rpm, "{:.0f}")
    R.m("FrRatioA", f90.expected_energy_kWh / f50.expected_energy_kWh, "{:.1f}"); R.m("FrRatioB", f99.expected_energy_kWh / f90.expected_energy_kWh, "{:.1f}")
    R.m("FrRatioAll", f99.expected_energy_kWh / f50.expected_energy_kWh, "{:.0f}")
    R.m("FrAllAccepted", "all" if fr.val_meets_alpha_lb.all() else f"{int(fr.val_meets_alpha_lb.sum())} of {len(fr)}")

    # hypothetical pump-capacity relaxation
    relax = []
    for mm in (0.50, 0.55, 0.60, 0.70):
        pp = p.with_(m_max=mm)
        b, _ = opt.chance_constrained_optimum(ALPHA, OPT, pp, m_step=0.01)
        relax.append(dict(flow_limit_kg_s=mm, rpm=b["rpm"], m_cmd=b["m"], expected_energy_kWh=b["E"],
                          pump_W_at_command=float(pump_power(b["m"], pp)), loop_dp_bar_at_command=float(loop_pressure_drop(b["m"], pp) / 1e5)))
    relax = R.csv(relax, "optimization/pump_limit_relaxation_hypothetical.csv")
    r55 = relax[relax.flow_limit_kg_s == 0.55].iloc[0]
    R.m("RelaxFiftyFiveE", r55.expected_energy_kWh, "{:.3f}"); R.m("RelaxFiftyFiveCut", 100 * (1 - r55.expected_energy_kWh / relax.iloc[0].expected_energy_kWh), "{:.0f}")

    # ================================================================ 7. case comparison + profiles + energy balance
    R.log("7. deterministic case comparison")
    cases = [("Coupled baseline", 300.0, 0.20), ("Nominal optimum", *NOM), ("v1 robust candidate", 300.0, 0.50),
             ("V3 recommended", *REC)]
    crow, prof = [], []
    for name, r, mc in cases:
        t, Y = rk4(r, mc, 1.0, p)
        tl, tj = time_constants(r, mc, p)
        crow.append(dict(case=name, rpm=r, mdot=mc, U_W_m2K=float(U_coef(r, mc, p)), UA_W_K=float(U_coef(r, mc, p) * p.A),
                         UA_eff_W_K=float(UA_effective(r, mc, p)), T_L_tf_K=Y[-1, 0], T_J_tf_K=Y[-1, 1], peak_T_J_K=Y[:, 1].max(),
                         time_to_target_s=sim.time_to_target(t, Y[:, 0], p.Ttarget), tau_tank_s=tl, tau_jacket_s=tj,
                         rk4_stable_dt_max_s=sim.rk4_stability_limit(r, mc, p), agitator_W=float(agitator_power(r, p)),
                         pump_W=float(pump_power(mc, p)), energy_kWh=float(batch_energy_kWh(r, mc, p))))
        for k in range(0, len(t), 30):
            prof.append(dict(case=name, time_s=t[k], T_L_K=Y[k, 0], T_J_K=Y[k, 1]))
    cdf = R.csv(crow, "deterministic/case_comparison.csv")
    R.csv(prof, "deterministic/temperature_profiles.csv")
    rc = cdf[cdf.case == "V3 recommended"].iloc[0]
    R.m("RecTL", rc.T_L_tf_K, "{:.2f}"); R.m("RecTtt", rc.time_to_target_s, "{:.0f}"); R.m("RecTttMin", rc.time_to_target_s / 60, "{:.1f}")
    R.m("RecTauJ", rc.tau_jacket_s, "{:.1f}"); R.m("StabRec", rc.rk4_stable_dt_max_s, "{:.0f}"); R.m("RecPeakTJ", rc.peak_T_J_K, "{:.1f}")
    R.m("RecAgW", rc.agitator_W, "{:,.0f}"); R.m("RecPumpW", rc.pump_W, "{:.0f}"); R.m("RecTJ", rc.T_J_tf_K, "{:.1f}")
    eb = [sim.energy_balance(300.0, 0.20, p), sim.energy_balance(*REC, p)]
    R.csv(eb, "deterministic/energy_balance.csv")
    R.m("EBQliq", eb[0]["Q_liquid_MJ"], "{:.4f}"); R.m("EBQout", eb[0]["Q_coolant_out_MJ"], "{:.4f}"); R.m("EBQst", eb[0]["Q_jacket_stored_MJ"], "{:.4f}")
    R.m("EBResTrap", eb[0]["residual_trapezoid_J"], "{:.1f}"); R.m("EBResSimp", eb[0]["residual_simpson_J"], "{:.4f}")
    R.m("EBRel", f"{eb[0]['relative_residual_simpson']:.0e}".replace("e-", "\\times10^{-").replace("{-0", "{-") + "}")
    R.m("EBResSimpRec", eb[1]["residual_simpson_J"], "{:.4f}")
    rom = sim.reduced_order_comparison([("Coupled baseline", 300.0, 0.20), ("Low flow", 300.0, 0.08), ("Nominal optimum", *NOM),
                                        ("V3 recommended", *REC), ("Large jacket (M_J = 100 kg)", 300.0, 0.20)], p)
    # the last row needs a different inventory
    rom[-1] = sim.reduced_order_comparison([("Large jacket (M_J = 100 kg)", 300.0, 0.20)], p.with_(MJ=100.0))[0]
    rom = R.csv(rom, "deterministic/reduced_order_model.csv")
    R.m("RomErrBase", rom.iloc[0].error_K, "{:.2f}"); R.m("RomErrMaxAbs", rom.error_K.abs().max(), "{:.2f}")
    R.m("RomErrRec", rom[rom.case == "V3 recommended"].iloc[0].error_K, "{:.2f}")

    # ================================================================ 8. sensitivity
    R.log("8. sensitivity")
    base_T, loc = sens.local_sensitivity(300.0, 0.20, p)
    _, loc_rec = sens.local_sensitivity(*REC, p)
    for r_ in loc_rec:
        r_["scale"] = "one sd (at V3 recommended)" if r_["scale"] == "one sd" else "v1 steps (at V3 recommended)"
    R.csv(loc + loc_rec, "sensitivity/local_sensitivity.csv")
    sob = {}
    srows = []
    for name, (r, mc) in [("Nominal optimum", NOM), ("V3 recommended", REC)]:
        s = sens.sobol_indices(r, mc, p)
        sob[name] = s
        for i, inp in enumerate(sens.SOBOL_INPUTS):
            srows.append(dict(design=name, input=inp, S1=s["S1"][i], S1_ci_lo=s["S1_lo"][i], S1_ci_hi=s["S1_hi"][i],
                              ST=s["ST"][i], ST_ci_lo=s["ST_lo"][i], ST_ci_hi=s["ST_hi"][i], N_base=s["N"], model_runs=s["n_runs"]))
    R.csv(srows, "sensitivity/sobol_indices.csv")
    sr = sob["V3 recommended"]
    for i, nm in enumerate(["U", "MJ", "Flow", "TZero"]):
        R.m("SobOne" + nm, 100 * max(sr["S1"][i], 0), "{:.0f}"); R.m("SobTot" + nm, 100 * max(sr["ST"][i], 0), "{:.0f}")
    R.m("SobolRuns", sr["n_runs"], "{:,}")

    # ================================================================ 9. adaptive policy
    R.log("9. adaptive policy")
    # Selection and calibration use CAL only. For every configuration the margin is
    # calibrated on CAL and then rounded up by the pre-declared practical rule; the
    # configuration with the lowest CAL mean energy at that rounded margin is selected.
    scan_rows = []
    for r1 in ADAPTIVE_GRID["rpm1"]:
        for to in ADAPTIVE_GRID["t_obs"]:
            d = ap.calibrate_delta(CAL, ALPHA, r1, REC[1], to, p, noise_seed=NOISE_SEEDS["CAL"], sigma=ap.SENSOR_SIGMA_K)
            row = dict(rpm1=r1, t_obs_s=to, delta_calibrated_K=d)
            if not np.isnan(d):
                dr = ap.practical_margin(d)
                rc_ = ap.run_policy(CAL, r1, REC[1], to, dr, p, sigma=ap.SENSOR_SIGMA_K, noise_seed=NOISE_SEEDS["CAL"])
                row.update(delta_K=dr, cal_pass=rc_["pass_frac"], cal_mean_energy_kWh=float(rc_["E"].mean()))
            scan_rows.append(row)
    ascan = pd.DataFrame(scan_rows)
    feas = ascan.dropna(subset=["delta_calibrated_K"])
    pick = feas.loc[feas.cal_mean_energy_kWh.idxmin()]
    ascan["selected"] = (ascan.rpm1 == pick.rpm1) & (ascan.t_obs_s == pick.t_obs_s)
    R.csv(ascan, "adaptive_policy/configuration_scan_CAL.csv")
    cfg = dict(rpm1=float(pick.rpm1), t_obs=float(pick.t_obs_s), delta=float(pick.delta_K), m_cmd=REC[1])
    AD, res, asum, ad_trail = ap.validation_protocol(cfg, VAL, VAL2, NOISE_SEEDS, p, ALPHA, ap.SENSOR_SIGMA_K)
    AD["delta_cal"] = float(pick.delta_calibrated_K)
    stage = "accepted on VAL; VAL2 not needed" if len(ad_trail) == 1 else "modified by the pre-declared rule and validated on VAL2"
    rows_ad = []
    for t in ad_trail:
        rows_ad.append(dict(stage=t["stage"], dataset=t["dataset"], rpm1=AD["rpm1"], t_obs=AD["t_obs"], delta_calibrated=AD["delta_cal"],
                            delta=t["delta"], sigma_K=ap.SENSOR_SIGMA_K, m_cmd=REC[1], rpm_max=p.rpm_max,
                            **{k: v for k, v in t.items() if k not in ("stage", "dataset", "delta")}))
        R.validations.append(dict(design="adaptive policy", stage=t["stage"], rpm=f"{AD['rpm1']:.0f} then N2", m_cmd=REC[1], dataset=t["dataset"],
                                  n=t["n"], pass_frac=t["pass_frac"], ci_lo=t["pass_ci_lo"], ci_hi=t["pass_ci_hi"],
                                  lb95_one_sided=t["pass_lb95_one_sided"], mean_energy_kWh=t["mean_E"], accepted=t["accepted"]))
    R.csv(rows_ad, "adaptive_policy/validation_summary.csv")
    R.m("AdDeltaCal", AD["delta_cal"], "{:.3f}")
    qs = np.linspace(0, 1, 201)
    R.csv(pd.DataFrame(dict(quantile=qs, energy_kWh=np.quantile(res["E"], qs), rpm2=np.quantile(res["rpm2"], qs),
                            T_L_tf_K=np.quantile(res["TLf"], qs))), "adaptive_policy/validation_quantiles.csv")
    R.m("AdRpmOne", AD["rpm1"], "{:.0f}"); R.m("AdTobs", AD["t_obs"], "{:.0f}"); R.m("AdTobsMin", AD["t_obs"] / 60, "{:.0f}")
    R.m("AdDelta", AD["delta"], "{:.2f}"); R.m("AdPass", 100 * asum["pass_frac"], "{:.2f}"); R.m("AdLo", 100 * asum["pass_ci_lo"], "{:.2f}")
    R.m("AdHi", 100 * asum["pass_ci_hi"], "{:.2f}"); R.m("AdLB", 100 * asum["pass_lb95_one_sided"], "{:.2f}")
    R.m("AdMeanE", asum["mean_E"], "{:.3f}"); R.m("AdMedE", asum["median_E"], "{:.3f}"); R.m("AdPNFE", asum["p95_E"], "{:.2f}")
    R.m("AdPNNE", asum["p99_E"], "{:.2f}"); R.m("AdMaxE", asum["max_E"], "{:.2f}"); R.m("AdCap", 100 * asum["frac_speed_capped"], "{:.2f}")
    R.m("AdUAErr", 100 * asum["UA_rel_err_sd"], "{:.1f}"); R.m("AdMedRpm", asum["median_rpm2"], "{:.0f}"); R.m("AdPNFRpm", asum["p95_rpm2"], "{:.0f}")
    R.m("AdSave", 100 * (1 - asum["mean_E"] / vrec["mean_E"]), "{:.0f}"); R.m("AdStage", stage)
    # V2's figure labelled every batch above 3.2 kWh as speed-capped; report both quantities separately
    R.m("AdFracAboveThree", 100 * float(np.mean(res["E"] > 3.2)), "{:.2f}")
    R.m("AdCapEnergy", float((agitator_power(AD["rpm1"], p) * AD["t_obs"] + agitator_power(p.rpm_max, p) * (p.tf - AD["t_obs"])
                              + pump_power(actual_flow(REC[1], 0.0, p), p) * p.tf) / 3.6e6), "{:.2f}")
    R.key("Nominal optimum P(pass) (VAL)", float(nomv.pass_frac), "-", "monte_carlo/designs_under_v3_model.csv")
    R.key("Adaptive policy margin delta", AD["delta"], "K", "adaptive_policy/validation_summary.csv")
    R.key("Adaptive policy P(pass) (independent)", asum["pass_frac"], "-", "adaptive_policy/validation_summary.csv")
    R.key("Adaptive policy mean energy (independent)", asum["mean_E"], "kWh/batch", "adaptive_policy/validation_summary.csv")
    R.key("Adaptive policy speed-capped batches", asum["frac_speed_capped"], "-", "adaptive_policy/validation_summary.csv")
    R.key("Adaptive policy 95th-percentile energy", asum["p95_E"], "kWh/batch", "adaptive_policy/validation_summary.csv")

    rob = []
    for sg in (0.05, 0.10, 0.20, 0.50):
        d = ap.practical_margin(ap.calibrate_delta(CAL, ALPHA, AD["rpm1"], REC[1], AD["t_obs"], p, noise_seed=NOISE_SEEDS["CAL"], sigma=sg, hi=4.0))
        rr = ap.run_policy(STRESS, AD["rpm1"], REC[1], AD["t_obs"], d, p, sigma=sg, noise_seed=NOISE_SEEDS["STRESS"])
        rob.append(dict(test="sensor noise", setting=f"sd = {sg} K", delta_K=d, pass_frac=rr["pass_frac"], mean_energy_kWh=float(rr["E"].mean()),
                        frac_speed_capped=float(np.mean(rr["capped"]))))
    for e in (0.10, 0.40):
        rr = ap.run_policy(STRESS, AD["rpm1"], REC[1], AD["t_obs"], AD["delta"], p, sigma=ap.SENSOR_SIGMA_K, noise_seed=NOISE_SEEDS["STRESS"], eN_plan=e)
        rob.append(dict(test="planner speed exponent wrong (truth 0.25)", setting=f"planner eN = {e}", delta_K=AD["delta"], pass_frac=rr["pass_frac"],
                        mean_energy_kWh=float(rr["E"].mean()), frac_speed_capped=float(np.mean(rr["capped"]))))
    rob = R.csv(rob, "adaptive_policy/robustness_STRESS.csv")
    R.m("RobHalfE", rob[rob.setting == "sd = 0.5 K"].iloc[0].mean_energy_kWh, "{:.2f}")
    R.m("RobLowEN", 100 * rob[rob.setting == "planner eN = 0.1"].iloc[0].pass_frac, "{:.0f}")
    R.m("RobHighEN", 100 * rob[rob.setting == "planner eN = 0.4"].iloc[0].pass_frac, "{:.0f}")

    # ================================================================ 10. model-form stress test
    R.log("10. model-form stress test")
    f_in, f_co = p.eN / FILM_EXP_INSIDE, p.eM / FILM_EXP_COOLANT
    f_w_nom = 1 - f_in - f_co
    R.m("FIn", f_in, "{:.3f}"); R.m("FCo", f_co, "{:.3f}"); R.m("FWall", f_w_nom, "{:.3f}")
    scen = {"Weak scaling (wall/fouling share 0.60)": 0.60, "Nominal closure": f_w_nom, "Strong scaling (no wall/fouling resistance)": 0.0}
    sc_rows = []
    for name, fw in scen.items():
        k = (1 - fw) / (f_in + f_co)
        eN, eM = round(p.eN * k, 3), round(p.eM * k, 3)
        pp = p.with_(eN=eN, eM=eM)
        b, _ = opt.chance_constrained_optimum(ALPHA, OPT, pp, m_step=0.01)
        fixed = pass_fraction(*REC, STRESS, pp)
        rr = ap.run_policy(STRESS, AD["rpm1"], REC[1], AD["t_obs"], AD["delta"], p, sigma=ap.SENSOR_SIGMA_K, noise_seed=NOISE_SEEDS["STRESS"], truth_p=pp)
        sc_rows.append(dict(scenario=name, wall_fouling_share=fw, eN=eN, eM=eM, T_L_tf_recommended_K=float(final_TL_exact(*REC, pp)),
                            fixed_recommended_pass=fixed, adaptive_pass=rr["pass_frac"], adaptive_mean_energy_kWh=float(rr["E"].mean()),
                            reoptimised_rpm=b["rpm"] if b else np.nan, reoptimised_m=b["m"] if b else np.nan,
                            reoptimised_expected_energy_kWh=b["E"] if b else np.nan))
    sc = R.csv(sc_rows, "model_form/exponent_scenarios.csv")
    w = sc.iloc[0]; s_ = sc.iloc[2]
    R.m("WeakEN", w.eN, "{:.3f}"); R.m("WeakEM", w.eM, "{:.3f}"); R.m("StrongEN", s_.eN, "{:.3f}"); R.m("StrongEM", s_.eM, "{:.3f}")
    R.m("WeakFixed", 100 * w.fixed_recommended_pass, "{:.0f}"); R.m("StrongFixed", 100 * s_.fixed_recommended_pass, "{:.1f}")
    R.m("WeakAd", 100 * w.adaptive_pass, "{:.1f}"); R.m("StrongAd", 100 * s_.adaptive_pass, "{:.1f}")
    R.m("WeakReE", "infeasible" if np.isnan(w.reoptimised_expected_energy_kWh) else f"{w.reoptimised_expected_energy_kWh:.2f}")
    R.m("StrongReE", s_.reoptimised_expected_energy_kWh, "{:.3f}")
    eNs = np.round(np.arange(0.05, 0.701, 0.05), 3); eMs = np.round(np.arange(0.05, 0.801, 0.05), 3)
    grid_rows = []
    S_small = {k: v[:20_000] if isinstance(v, np.ndarray) else v for k, v in STRESS.items()}
    for eN in eNs:
        for eM in eMs:
            pp = p.with_(eN=float(eN), eM=float(eM))
            grid_rows.append(dict(eN=eN, eM=eM, admissible=bool(eN / FILM_EXP_INSIDE + eM / FILM_EXP_COOLANT <= 1 + 1e-9),
                                  fixed_recommended_pass=pass_fraction(*REC, S_small, pp)))
    mgrid = R.csv(grid_rows, "model_form/exponent_grid_recommended_pass.csv")
    adm = mgrid[mgrid.admissible]
    R.m("GridMinPass", 100 * adm.fixed_recommended_pass.min(), "{:.0f}"); R.m("GridFracAdmFail", 100 * float(np.mean(adm.fixed_recommended_pass < ALPHA)), "{:.0f}")

    # jacket flow pattern (tanks in series): deterministic and reliability
    jf = []
    for name, r, mc in [("Coupled baseline", 300.0, 0.20), ("V3 recommended", *REC)]:
        for ncell in (1, 2, 5, 20):
            dt = 1.0 if ncell <= 5 else 0.25
            jf.append(dict(case=name, rpm=r, mdot=mc, n_cells=ncell, T_L_tf_K=float(rk4_tanks_in_series(r, mc, dt, ncell, p)),
                           UA_eff_mixed=float(UA_effective(r, mc, p)), UA_eff_plugflow=float(UA_effective_plugflow(r, mc, p))))
    jf = pd.DataFrame(jf)
    S_j = {k: v[:10_000] if isinstance(v, np.ndarray) else v for k, v in STRESS.items()}
    md = actual_flow(REC[1], S_j["z_flow"], p)
    for ncell, dt in ((1, 2.0), (10, 1.0)):
        TLs = rk4_tanks_in_series(REC[0], md, dt, ncell, p, Uref=S_j["Uref"], MJ=S_j["MJ"], T0=S_j["T0"])
        jf.loc[(jf.case == "V3 recommended") & (jf.n_cells == ncell), "pass_frac_STRESS10k"] = float(np.mean(TLs <= p.Ttarget))
        if ncell == 10:
            jf = pd.concat([jf, pd.DataFrame([dict(case="V3 recommended", rpm=REC[0], mdot=REC[1], n_cells=10,
                                                   T_L_tf_K=float(rk4_tanks_in_series(REC[0], REC[1], 0.5, 10, p)),
                                                   pass_frac_STRESS10k=float(np.mean(TLs <= p.Ttarget)))])], ignore_index=True)
    jf = R.csv(jf, "model_form/jacket_flow_pattern.csv")
    j1 = jf[(jf.case == "V3 recommended") & (jf.n_cells == 1)].iloc[0]; j20 = jf[(jf.case == "V3 recommended") & (jf.n_cells == 20)].iloc[0]
    jb1 = jf[(jf.case == "Coupled baseline") & (jf.n_cells == 1)].iloc[0]; jb20 = jf[(jf.case == "Coupled baseline") & (jf.n_cells == 20)].iloc[0]
    j10 = jf[(jf.case == "V3 recommended") & (jf.n_cells == 10)].dropna(subset=["pass_frac_STRESS10k"]).iloc[-1]
    R.m("JetRecOne", j1.T_L_tf_K, "{:.2f}"); R.m("JetRecTwenty", j20.T_L_tf_K, "{:.2f}")
    R.m("JetBaseOne", jb1.T_L_tf_K, "{:.2f}"); R.m("JetBaseTwenty", jb20.T_L_tf_K, "{:.2f}")
    R.m("JetPassOne", 100 * j1.pass_frac_STRESS10k, "{:.1f}"); R.m("JetPassTen", 100 * j10.pass_frac_STRESS10k, "{:.1f}")
    R.m("JetGainBase", 100 * (jb1.UA_eff_plugflow / jb1.UA_eff_mixed - 1), "{:.0f}"); R.m("JetGainRec", 100 * (j1.UA_eff_plugflow / j1.UA_eff_mixed - 1), "{:.0f}")

    # isolated thermal-inventory sensitivity
    mj = []
    for MJ in (5.0, 10.0, 25.0, 50.0, 100.0):
        pp = p.with_(MJ=MJ)
        mj.append(dict(MJ_kg=MJ, tau_jacket_s=time_constants(*REC, pp)[1], T_L_tf_recommended_K=float(final_TL_exact(*REC, pp)),
                       T_L_tf_baseline_K=float(final_TL_exact(300.0, 0.20, pp))))
    mj = R.csv(mj, "model_form/jacket_inventory_isolated.csv")
    R.m("MJSpread", mj.T_L_tf_recommended_K.max() - mj.T_L_tf_recommended_K.min(), "{:.2f}")
    R.m("MJTauHundred", mj[mj.MJ_kg == 100].iloc[0].tau_jacket_s, "{:.0f}")

    # ================================================================ 11. equipment and process checks
    R.log("11. equipment checks")
    Dt = (4 * p.V / math.pi) ** (1 / 3)
    mu = 0.80e-3   # Pa s, water near 305 K (approximate)
    perim = math.pi * Dt
    eq = [
        ("Tank diameter for H = T (0.50 m3)", Dt, "m", "geometry implied by V; no tank drawing exists"),
        ("Wetted wall + base area at H = T", math.pi * Dt * Dt + math.pi * Dt ** 2 / 4, "m2", "A = 2.0 m2 is of the same order"),
        ("Impeller-to-tank diameter ratio D/T", p.D / Dt, "-", "within the commonly used 0.3-0.5 range"),
        ("Impeller Reynolds number at 120 rpm (water-like)", p.rhoL * (120 / 60) * p.D ** 2 / mu, "-", "turbulent: constant N_P applicable"),
        ("Shaft power per volume at recommended speed", agitator_power(REC[0], p) * p.eta_ag / p.V / 1e3, "kW/m3", "high-intensity duty; check with vendor"),
        ("Loop pressure drop at 0.50 kg/s", loop_pressure_drop(0.5, p) / 1e5, "bar", "extrapolated from assumed 0.5 bar at 0.2 kg/s"),
        ("Pump flow turndown modelled (0.05-0.50 kg/s)", 10.0, ":1", "model assumes variable-speed pump, constant efficiency"),
        ("Jacket gap implied by M_J over A", p.MJ / p.rhoC / p.A * 1e3, "mm", "thin; see isolated-inventory study"),
        ("Plain annular jacket Re at 0.50 kg/s, Re = 2 m/(pi T mu)", 2 * 0.5 / (perim * mu), "-", "laminar: assumed U needs a baffled/spiral jacket"),
        ("Coolant outlet maximum (baseline peak T_J)", float(YB[:, 1].max()), "K", "about 45 C; confirm against site cooling-water limits"),
        ("Liquid density range 300-370 K (water)", 960.0, "kg/m3", "about 996 to 960 kg/m3; constant 1000 used (<= 4 %)"),
        ("Liquid Cp range 300-370 K (water)", 4.21, "kJ/(kg K)", "about 4.18 to 4.21; constant 4.18 used (< 1 %)"),
    ]
    R.csv([dict(check=a, value=b, unit=c, comment=d) for a, b, c, d in eq], "deterministic/equipment_and_process_checks.csv")
    R.m("TankD", Dt, "{:.2f}"); R.m("DT", p.D / Dt, "{:.2f}"); R.m("PV", agitator_power(REC[0], p) * p.eta_ag / p.V / 1e3, "{:.1f}")
    R.m("DPMax", loop_pressure_drop(0.5, p) / 1e5, "{:.2f}"); R.m("Gap", p.MJ / p.rhoC / p.A * 1e3, "{:.0f}")
    R.m("ReJacket", 2 * 0.5 / (perim * mu), "{:.0f}"); R.m("ReImp", p.rhoL * 2 * p.D ** 2 / mu, "{:,.0f}")

    # ================================================================ 12. Monte Carlo data for figures
    R.log("12. Monte Carlo summaries")
    q = np.linspace(0, 1, 501)
    ecdf = {"V3 recommended": mc_final_TL(*REC, VAL, p), "Nominal optimum": mc_final_TL(*NOM, VAL, p),
            "v1 robust candidate": mc_final_TL(300.0, 0.5, VAL, p), "v2 recommended": mc_final_TL(285.0, 0.5, VAL, p),
            "Adaptive policy": res["TLf"]}
    R.csv(pd.DataFrame({"quantile": q, **{k: np.quantile(v, q) for k, v in ecdf.items()}}), "monte_carlo/final_temperature_quantiles_VAL.csv")
    xr = ecdf["V3 recommended"] <= p.Ttarget
    nn = np.unique(np.logspace(2, np.log10(len(xr)), 70).astype(int))
    cs = np.cumsum(xr)
    run = [dict(n=n, pass_frac=cs[n - 1] / n, ci_lo=wilson_interval(cs[n - 1], n)[0], ci_hi=wilson_interval(cs[n - 1], n)[1]) for n in nn]
    R.csv(run, "monte_carlo/convergence_recommended_VAL.csv")

    # ================================================================ 13. figures
    R.log("13. figures")
    make_figures(R, p, dict(tB=tB, YB=YB, T_orig=T_orig, conv=pd.DataFrame(conv), order=order, grid=grid, NOM=NOM, REC=REC,
                            OPT=OPT, VAL=VAL, fr=fr, leg=leg, loc=pd.DataFrame(loc), sobol=pd.DataFrame(srows), res=res, asum=asum,
                            AD=AD, vrec=vrec, ecdf=ecdf, run=pd.DataFrame(run), mgrid=mgrid, sc=sc, jf=jf, p_req=p_req))

    # ================================================================ 14. outputs
    R.log("14. outputs and report values")
    pd.DataFrame(R.summary).to_csv(R.out / "summary_results.csv", index=False, float_format="%.10g")
    pd.DataFrame(R.validations).to_csv(R.out / "validation_results.csv", index=False, float_format="%.10g")
    meta = dict(project="jacketed-stirred-tank-cooling-optimization", version="3.0.0",
                python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__,
                matplotlib=matplotlib.__version__, platform=platform.system(),
                datasets={k: dict(seed=v["seed"], n=v["n"], purpose=v["purpose"]) for k, v in DATASETS.items()},
                noise_seeds=NOISE_SEEDS, sobol_seed=SOBOL_SEED, reliability_target=ALPHA,
                acceptance_rule="one-sided 95% Wilson lower bound of P(pass) >= 0.95",
                runtime_s=round(time.time() - R.t0, 1))
    (R.out / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    from .crosscheck import build as build_crosscheck
    xc = build_crosscheck(R.root, dict(T_orig=float(T_orig[-1]), T_base=float(YB[-1, 0]), TJ_base=float(YB[-1, 1]), err_dt1=err_dt1,
                                       order=order, eb_simp=eb[0]["residual_simpson_J"], nom_rpm=nom["rpm"], nom_E=nom["E"],
                                       cc_rpm=best["rpm"], rec_rpm=REC[0], rec_pass=vrec["pass_frac"], rec_E=vrec["mean_E"],
                                       v2_pass=float(v2r.pass_frac), ad_delta_cal=AD["delta_cal"], ad_delta=AD["delta"],
                                       ad_pass=asum["pass_frac"], ad_lb=asum["pass_lb95_one_sided"], ad_E=asum["mean_E"]))
    if xc is not None:
        R.csv(xc, "verification/python_vs_octave.csv")
        R.m("OctCCRpm", xc[xc.quantity == "Chance-constrained speed on OPT [rpm]"].iloc[0].octave, "{:.1f}")
        R.m("OctRecPass", 100 * xc[xc.quantity == "Recommended P(pass) on VAL"].iloc[0].octave, "{:.2f}")
        R.m("OctAdPass", 100 * xc[xc.quantity == "Adaptive P(pass) on VAL"].iloc[0].octave, "{:.2f}")
        R.m("OctAdLB", 100 * xc[xc.quantity == "Adaptive one-sided lower bound on VAL"].iloc[0].octave, "{:.2f}")
        R.m("OctRuntime", xc.runtime.iloc[0])
    write_latex(R, p, dict(xc=xc, fr=fr, leg=leg, cdf=cdf, rom=rom, sc=sc, ascan=ascan, ad_val=pd.DataFrame(rows_ad), rob=rob, conv=pd.DataFrame(conv),
                           loc=pd.DataFrame(loc), sobol=pd.DataFrame(srows), audit=pd.DataFrame(audit), eb=pd.DataFrame(eb),
                           relax=relax, act=act, jf=jf))
    R.log("done")
    return R


