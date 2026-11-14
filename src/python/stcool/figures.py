"""All report and repository figures, drawn from results computed in pipeline.py."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .model import FILM_EXP_INSIDE, FILM_EXP_COOLANT, rk4, UA_effective, quasi_steady_TL, analytic_boundary_rpm
from .uncertainty import FLOW_ERROR_REL_SD, actual_flow
from . import optimization as opt
from .adaptive_policy import SENSOR_SIGMA_K
from .plotting import C, ENTITY, save


def _math(label):
    """Typeset plain-text symbols from the CSV labels for figure tick labels."""
    return (label.replace("U_ref", "$U_{ref}$").replace("M_J", "$M_J$").replace("T_0", "$T_0$")
            .replace("+/-", "\u00b1").replace("W/(m2 K)", "W/(m\u00b2 K)"))


def make_figures(R, p, d):
    """R: the pipeline Run (paths); p: Params; d: results computed in pipeline.main."""
    F = R.fig
    tB, YB = d["tB"], d["YB"]
    REC, NOM = d["REC"], d["NOM"]

    # model comparison
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    tm = tB / 60
    ax.plot(tm, d["T_orig"], "--", color=ENTITY["original"], label="Original model (coolant fixed at 300 K)")
    ax.plot(tm, YB[:, 0], color=ENTITY["baseline"], label="Coupled model, RK4 (liquid)")
    ax.plot(tm, quasi_steady_TL(UA_effective(300, 0.2, p), tB, p), ":", color=C["orange"], lw=2.2, label="Reduced-order model, $UA_{eff}$")
    ax.axhline(p.Ttarget, color=ENTITY["target"], lw=1.1, ls="-.")
    ax.text(0.3, p.Ttarget + 0.35, "target 355 K", color=ENTITY["target"], fontsize=8.5)
    for y, col, lab, dy in [(d["T_orig"][-1], ENTITY["original"], f"{d['T_orig'][-1]:.2f} K original (passes)", -1.7),
                            (YB[-1, 0], ENTITY["baseline"], f"{YB[-1, 0]:.2f} K coupled (fails)", 0.9)]:
        ax.plot(30, y, "o", color=col, ms=6, zorder=5)
        ax.annotate(lab, xy=(30, y), xytext=(30.8, y - 0.25), fontsize=8.5, color=C["ink"])
    ax.set_xlabel("Time (min)"); ax.set_ylabel("Liquid temperature (K)"); ax.set_ylim(350, 371); ax.set_xlim(-1, 41); ax.set_xticks(range(0, 31, 5))
    ax.set_title("Original vs coupled model at 300 rpm, 0.20 kg/s")
    ax.legend(loc="upper right", fontsize=8)
    save(fig, F / "model_comparison" / "original_vs_coupled.png")

    # temperature profiles
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    for lab, ent, (r, mc) in [("Coupled baseline (300 rpm, 0.20 kg/s)", "baseline", (300.0, 0.2)),
                              (f"Nominal optimum ({NOM[0]:.0f} rpm, 0.50 kg/s)", "nominal", NOM),
                              (f"V3 recommended ({REC[0]:.0f} rpm, {REC[1]:.2f} kg/s)", "recommended", REC)]:
        t, Y = rk4(r, mc, 1.0, p)
        ax.plot(t / 60, Y[:, 0], color=ENTITY[ent], label=lab)
        ax.plot(t / 60, Y[:, 1], color=ENTITY[ent], ls=":", lw=1.5)
    ax.axhline(p.Ttarget, color=ENTITY["target"], lw=1.1, ls="-.")
    ax.text(0.3, p.Ttarget + 1.2, "target 355 K", color=ENTITY["target"], fontsize=8.5)
    ax.text(12, 321, "dotted: jacket coolant", color=C["ink2"], fontsize=8.5)
    ax.set_xlabel("Time (min)"); ax.set_ylabel("Temperature (K)"); ax.set_title("Liquid (solid) and jacket (dotted) temperatures, RK4")
    ax.legend(loc="center right", fontsize=8)
    save(fig, F / "temperature_profiles" / "temperature_profiles.png")

    # verification
    conv = d["conv"]
    fig, ax = plt.subplots(figsize=(5.8, 3.7))
    ax.loglog(conv.dt_s, conv.max_traj_error_K, "o-", color=C["blue"], ms=5, label="RK4 max trajectory error vs exact solution")
    xx = np.array([2, 60]); ref = conv.max_traj_error_K[conv.dt_s == 10].item() * (xx / 10) ** 4
    ax.loglog(xx, ref, "--", color=C["muted"], lw=1.2, label="slope 4 reference")
    ax.set_xlabel("RK4 time step, dt (s)"); ax.set_ylabel("Max |T$_{RK4}$ - T$_{exact}$| over 0-1800 s (K)")
    ax.set_title("Fourth-order convergence of RK4")
    ax.text(0.04, 0.84, f"fitted order = {d['order']:.2f} (dt = 2-30 s)", transform=ax.transAxes, color=C["ink2"], fontsize=8.5)
    ax.legend(loc="lower right", fontsize=8)
    save(fig, F / "verification" / "rk4_convergence.png")

    # actuator distribution
    VAL = d["VAL"]
    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.4), gridspec_kw=dict(width_ratios=[1.25, 1]))
    ax = axs[0]
    unsat = np.sort(0.5 * (1 + FLOW_ERROR_REL_SD * VAL["z_flow"])); sat = np.sort(actual_flow(0.5, VAL["z_flow"], p))
    cdf = np.arange(1, len(sat) + 1) / len(sat)
    ax.plot(unsat, cdf, color=C["muted"], lw=1.8, ls="--", label="V2: unbounded")
    ax.plot(sat, cdf, color=C["aqua"], lw=2.2, label="V3: saturated")
    mass = float(np.mean(sat >= p.m_max - 1e-12))
    ax.axvline(p.m_max, color=C["red"], lw=1.1, ls="-.")
    ax.annotate(f"jump of {100 * mass:.1f} % at the\n0.50 kg/s pump limit", xy=(p.m_max, 0.75), xytext=(0.515, 0.42),
                fontsize=8, color=C["ink2"], arrowprops=dict(arrowstyle="->", color=C["muted"]))
    ax.text(0.515, 0.08, f"V2 model: {100 * np.mean(unsat > p.m_max):.1f} % of flows\nabove the pump limit", fontsize=8, color=C["ink2"])
    ax.set_xlim(0.42, 0.60); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Actual coolant flow at a 0.50 kg/s command (kg/s)"); ax.set_ylabel("Cumulative probability (VAL)")
    ax.set_title("(a) Saturated flow actuator"); ax.legend(fontsize=7.6, loc="upper left")
    ax = axs[1]
    leg = d["leg"]
    order = ["v1 nominal grid optimum", "v2/v3 nominal optimum", "v2 recommended", "v1 robust candidate", "v3 recommended"]
    lab = ["V1 nominal\n125/0.49", "nominal\n120/0.50", "V2 rec.\n285/0.50", "V1 robust\n300/0.50", f"V3 rec.\n{REC[0]:.0f}/{REC[1]:.2f}"]
    vals = [leg[leg.design == o].iloc[0] for o in order]
    cols = [ENTITY["nominal"], ENTITY["nominal"], ENTITY["v2 recommended"], ENTITY["v1 robust"], ENTITY["recommended"]]
    xs = np.arange(len(order))
    ax.bar(xs, [100 * v.pass_frac for v in vals], color=cols, width=0.62)
    ax.errorbar(xs, [100 * v.pass_frac for v in vals], yerr=[[100 * (v.pass_frac - v.ci_lo) for v in vals], [100 * (v.ci_hi - v.pass_frac) for v in vals]],
                fmt="none", ecolor=C["ink"], capsize=3, lw=1)
    for x, v in zip(xs, vals):
        ax.text(x, 100 * v.pass_frac + 1.2, f"{100 * v.pass_frac:.1f}", ha="center", fontsize=8)
    ax.axhline(95, color=C["red"], ls="-.", lw=1.1)
    ax.set_ylim(40, 101); ax.set_xticks(xs, lab, fontsize=7.6); ax.set_ylabel("P(T$_L$ <= 355 K), % (VAL, n = 200k)")
    ax.set_title("(b) Designs under the corrected model")
    save(fig, F / "monte_carlo" / "flow_actuator_and_designs.png")

    # operating window
    g = d["grid"]
    fig, ax = plt.subplots(figsize=(6.8, 4.3))
    cf = ax.contourf(g["M"], g["R"], g["E_grid"], levels=np.arange(0, 8.6, 0.4), cmap="Blues")
    cb = fig.colorbar(cf, ax=ax); cb.set_label("Auxiliary electricity per batch (kWh)")
    ax.contour(g["M"], g["R"], g["T_grid"], levels=[p.Ttarget], colors=[C["ink"]], linewidths=2)
    mm = np.linspace(0.07, 0.5, 400)
    ax.plot(mm, analytic_boundary_rpm(mm, p), color=C["yellow"], ls="--", lw=1.6)
    flows = np.round(np.arange(0.30, 0.5001, 0.01), 3)
    b95 = [opt.min_rpm_meeting(m_, d["p_req"], d["OPT"], p, tol=0.05) for m_ in flows]
    ax.plot(flows, b95, color=C["aqua"], lw=2.2)
    for (r, mc), col, mk in [(NOM, ENTITY["nominal"], "*"), ((300.0, 0.5), ENTITY["v1 robust"], "D"), (REC, ENTITY["recommended"], "o")]:
        ax.plot(mc, r, mk, color=col, mec="white", mew=1.1, ms=12 if mk == "*" else 8, clip_on=False, zorder=6)
    hd = [Line2D([], [], color=C["ink"], lw=2, label="T$_L$(1800 s) = 355 K (nominal model)"),
          Line2D([], [], color=C["yellow"], ls="--", lw=1.6, label="closed-form quasi-steady boundary"),
          Line2D([], [], color=C["aqua"], lw=2.2, label="P(pass) lower bound = 95 % (OPT data)"),
          Line2D([], [], marker="*", ls="", color=ENTITY["nominal"], mec="white", ms=11, label=f"nominal optimum {NOM[0]:.1f} rpm / 0.50"),
          Line2D([], [], marker="D", ls="", color=ENTITY["v1 robust"], mec="white", ms=7, label="V1 robust 300 rpm / 0.50"),
          Line2D([], [], marker="o", ls="", color=ENTITY["recommended"], mec="white", ms=7, label=f"V3 recommended {REC[0]:.0f} rpm / {REC[1]:.2f}")]
    ax.legend(handles=hd, loc="upper left", fontsize=7.2, frameon=True, facecolor="white", framealpha=0.92)
    ax.set_xlim(0.05, 0.505); ax.set_ylim(100, 600)
    ax.set_xlabel("Commanded coolant flow (kg/s); 0.50 = pump limit"); ax.set_ylabel("Agitator speed (rpm)")
    ax.set_title("Operating window")
    save(fig, F / "optimization" / "operating_window.png")

    # frontier
    fr, leg = d["fr"], d["leg"]
    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    ax.errorbar(100 * fr.val_pass, fr.expected_energy_kWh, xerr=[100 * (fr.val_pass - fr.val_ci_lo), 100 * (fr.val_ci_hi - fr.val_pass)],
                fmt="o-", color=C["blue"], ms=5, capsize=2, label="minimum-energy fixed setting (frontier)")
    for _, rw in fr.iterrows():
        if rw.alpha in (0.5, 0.9, 0.99):
            ax.annotate(f"{rw.rpm:.0f} rpm", (100 * rw.val_pass, rw.expected_energy_kWh), xytext=(6, 8) if rw.alpha == 0.5 else ((-92, -4) if rw.alpha == 0.99 else (-44, 6)),
                        textcoords="offset points", fontsize=8, color=C["ink2"])
    v1r = leg[leg.design == "v1 robust candidate"].iloc[0]; v2r = leg[leg.design == "v2 recommended"].iloc[0]
    ax.plot(100 * v2r.pass_frac, v2r.mean_energy_kWh, "s", color=ENTITY["v2 recommended"], ms=8, mec="white", label="V2 recommendation 285 rpm (now < 95 %)")
    ax.plot(100 * v1r.pass_frac, v1r.mean_energy_kWh, "D", color=ENTITY["v1 robust"], ms=8, mec="white", label="V1 robust 300 rpm")
    ax.plot(100 * d["vrec"]["pass_frac"], d["vrec"]["mean_E"], "o", color=ENTITY["recommended"], ms=9, mec="white", label=f"V3 recommended {REC[0]:.0f} rpm")
    ax.plot(100 * d["asum"]["pass_frac"], d["asum"]["mean_E"], "P", color=ENTITY["adaptive"], ms=10, mec="white", label="adaptive policy (mean)")
    ax.axvline(95, color=C["red"], ls="-.", lw=1)
    ax.set_xlabel("P(T$_L$(1800 s) <= 355 K), % (independent VAL set, n = 200k)"); ax.set_ylabel("Expected auxiliary energy (kWh/batch)")
    ax.set_title("Reliability-energy frontier"); ax.legend(fontsize=7.4, loc="upper left")
    save(fig, F / "optimization" / "reliability_energy_frontier.png")

    # ECDF + convergence
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.6), gridspec_kw=dict(width_ratios=[1.3, 1]))
    ax = axs[0]
    for name, ent in [("Nominal optimum", "nominal"), ("v2 recommended", "v2 recommended"), ("V3 recommended", "recommended"),
                      ("v1 robust candidate", "v1 robust"), ("Adaptive policy", "adaptive")]:
        x = np.sort(d["ecdf"][name])
        shown = {"v2 recommended": "V2 recommended (285 rpm)", "v1 robust candidate": "V1 robust (300 rpm)",
                 "V3 recommended": f"V3 recommended ({REC[0]:.0f} rpm)"}.get(name, name)
        ax.plot(x, np.arange(1, len(x) + 1) / len(x), color=ENTITY[ent], lw=1.8, label=f"{shown}: {100 * np.mean(x <= p.Ttarget):.1f} %")
    ax.axvline(p.Ttarget, color=C["red"], ls="-.", lw=1); ax.axhline(0.95, color=C["muted"], ls=":", lw=1)
    ax.set_xlim(346, 361); ax.set_xlabel("T$_L$ at 1800 s (K)"); ax.set_ylabel("Cumulative probability")
    ax.set_title("(a) Final-temperature distributions"); ax.legend(fontsize=7.2, loc="lower right")
    ax = axs[1]; run = d["run"]
    ax.fill_between(run.n, 100 * run.ci_lo, 100 * run.ci_hi, color=C["aqua"], alpha=0.25, lw=0, label="95 % Wilson interval")
    ax.plot(run.n, 100 * run.pass_frac, color=C["aqua"], lw=1.8, label="running estimate")
    ax.axhline(95, color=C["red"], ls="-.", lw=1)
    ax.set_xscale("log"); ax.set_ylim(85, 100); ax.set_xlabel("Number of Monte Carlo batches"); ax.set_ylabel("P(pass) estimate (%)")
    ax.set_title(f"(b) Sampling error at {REC[0]:.0f} rpm"); ax.legend(fontsize=7.6, loc="lower right")
    save(fig, F / "monte_carlo" / "final_temperature_and_convergence.png")

    # sensitivity
    loc, sob = d["loc"], d["sobol"]
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.8), gridspec_kw=dict(width_ratios=[1.05, 1]))
    ax = axs[0]
    sl = loc[loc.scale == "one sd"].sort_values("swing_K")
    yy = np.arange(len(sl))
    for y, (_, rw) in zip(yy, sl.iterrows()):
        ax.plot([rw.delta_low_K, rw.delta_high_K], [y, y], color=C["grid"], lw=6, solid_capstyle="butt", zorder=1)
    ax.scatter(sl.delta_low_K, yy, color=C["blue"], s=40, zorder=3, label="input at mean - 1 sd")
    ax.scatter(sl.delta_high_K, yy, color=C["orange"], s=40, zorder=3, label="input at mean + 1 sd")
    ax.set_yticks(yy, [f"{_math(a)}\n({_math(b)})" for a, b in zip(sl.parameter, sl.perturbation)], fontsize=8)
    ax.axvline(0, color=C["ink2"], lw=1); ax.set_xlabel("Change in T$_L$(1800 s) at the baseline (K)")
    ax.set_title("(a) Local, one standard deviation"); ax.legend(fontsize=7.6, loc="lower right")
    ax = axs[1]
    sd_ = sob[sob.design == "V3 recommended"].iloc[::-1]
    yy = np.arange(len(sd_)); h = 0.36
    ax.barh(yy + h / 2, sd_.S1.clip(lower=0), h, color=C["blue"], label="first-order S$_1$",
            xerr=[(sd_.S1 - sd_.S1_ci_lo).clip(lower=0), (sd_.S1_ci_hi - sd_.S1).clip(lower=0)], ecolor=C["ink2"], capsize=2)
    ax.barh(yy - h / 2, sd_.ST.clip(lower=0), h, color=C["orange"], label="total effect S$_T$",
            xerr=[(sd_.ST - sd_.ST_ci_lo).clip(lower=0), (sd_.ST_ci_hi - sd_.ST).clip(lower=0)], ecolor=C["ink2"], capsize=2)
    for y, s1 in zip(yy, sd_.S1):
        ax.text(max(s1, 0) + 0.03, y + h / 2, f"{100 * max(s1, 0):.0f}%", va="center", fontsize=8)
    ax.set_yticks(yy, [_math(x) for x in sd_.input], fontsize=8); ax.set_xlim(0, 1); ax.set_xlabel("Share of variance of T$_L$(1800 s)")
    ax.set_title(f"(b) Sobol indices at {REC[0]:.0f} rpm, {REC[1]:.2f} kg/s"); ax.legend(fontsize=7.6, loc="lower right")
    save(fig, F / "sensitivity" / "local_and_sobol.png")

    # adaptive policy
    res, asum, AD = d["res"], d["asum"], d["AD"]
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.8))
    ax = axs[0]
    idx = np.random.default_rng(3).choice(len(res["UA_true"]), 3000, replace=False)
    ax.scatter(res["UA_true"][idx], res["UA_hat"][idx], s=5, alpha=0.35, color=C["blue"], edgecolor="none")
    lim = [res["UA_true"][idx].min(), res["UA_true"][idx].max()]
    ax.plot(lim, lim, color=C["ink2"], lw=1.1, label="perfect estimate")
    ax.set_xlabel("True UA in phase 1 (W/K)"); ax.set_ylabel("Soft-sensor estimate of UA (W/K)")
    ax.set_title(f"(a) UA from {AD['t_obs'] / 60:.0f} min of noisy T$_L$, T$_J$ data")
    ax.text(0.04, 0.86, f"relative error sd = {100 * asum['UA_rel_err_sd']:.1f} %\nsensor noise sd = {SENSOR_SIGMA_K} K", transform=ax.transAxes, fontsize=8.5, color=C["ink2"])
    ax.legend(loc="lower right", fontsize=8)
    ax = axs[1]
    E = res["E"]; capped = res["capped"]
    bins = np.linspace(0, float(np.ceil(E.max() * 2) / 2), 70)
    ax.hist(E[~capped], bins=bins, color=ENTITY["adaptive"], alpha=0.85, label="batches below the speed cap")
    ax.hist(E[capped], bins=bins, color=C["red"], alpha=0.9, label=f"speed-capped at 600 rpm ({100 * capped.mean():.2f} %)")
    ax.axvline(d["vrec"]["mean_E"], color=ENTITY["recommended"], lw=2, label=f"fixed {REC[0]:.0f} rpm mean ({d['vrec']['mean_E']:.2f} kWh)")
    ax.axvline(E.mean(), color=C["ink"], lw=1.3, ls=":", label=f"adaptive mean ({E.mean():.2f} kWh)")
    ax.set_yscale("log")
    ax.set_xlabel("Auxiliary electricity per batch (kWh)"); ax.set_ylabel("Batches (of 200,000), log scale")
    ax.set_title(f"(b) Energy per batch, P(pass) = {100 * asum['pass_frac']:.1f} %"); ax.legend(fontsize=7.4, loc="upper right")
    save(fig, F / "adaptive_policy" / "adaptive_policy.png")

    # model form
    mg, sc = d["mgrid"], d["sc"]
    fig, axs = plt.subplots(1, 2, figsize=(9.8, 3.9), gridspec_kw=dict(width_ratios=[1.15, 1]))
    ax = axs[0]
    piv = mg.pivot(index="eN", columns="eM", values="fixed_recommended_pass")
    ext = [piv.columns.min() - 0.025, piv.columns.max() + 0.025, piv.index.min() - 0.025, piv.index.max() + 0.025]
    im = ax.imshow(100 * piv.values, origin="lower", extent=ext, aspect="auto", cmap="RdBu", vmin=50, vmax=100)
    cb = fig.colorbar(im, ax=ax); cb.set_label(f"P(pass) of {REC[0]:.0f} rpm setting (%)")
    ax.contour(piv.columns, piv.index, 100 * piv.values, levels=[95], colors=[C["ink"]], linewidths=1.4)
    em = np.linspace(0, 0.8, 50); ax.plot(em, FILM_EXP_INSIDE * (1 - em / FILM_EXP_COOLANT), color=C["ink2"], ls="--", lw=1.4)
    ax.text(0.55, 0.62, "inadmissible\n(resistance shares > 1)", fontsize=7.5, color=C["ink"], ha="center")
    for _, rw in sc.iterrows():
        ax.plot(rw.eM, rw.eN, "o", color="white", mec=C["ink"], ms=7)
    ax.annotate("weak", (sc.iloc[0].eM, sc.iloc[0].eN), xytext=(6, -10), textcoords="offset points", fontsize=8)
    ax.annotate("nominal", (sc.iloc[1].eM, sc.iloc[1].eN), xytext=(6, 4), textcoords="offset points", fontsize=8)
    ax.annotate("strong", (sc.iloc[2].eM, sc.iloc[2].eN), xytext=(6, 4), textcoords="offset points", fontsize=8)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3]); ax.grid(False)
    ax.set_xlabel("coolant-flow exponent e$_M$"); ax.set_ylabel("agitation exponent e$_N$")
    ax.set_title("(a) Fixed setting under other U closures")
    ax = axs[1]
    xs = np.arange(3); w_ = 0.36
    ax.bar(xs - w_ / 2, 100 * sc.fixed_recommended_pass, w_, color=ENTITY["recommended"], label=f"fixed {REC[0]:.0f} rpm")
    ax.bar(xs + w_ / 2, 100 * sc.adaptive_pass, w_, color=ENTITY["adaptive"], label="adaptive policy (planner assumes nominal)")
    for x, a, b in zip(xs, sc.fixed_recommended_pass, sc.adaptive_pass):
        ax.text(x - w_ / 2, 100 * a + 1, f"{100 * a:.0f}", ha="center", fontsize=8); ax.text(x + w_ / 2, 100 * b + 1, f"{100 * b:.0f}", ha="center", fontsize=8)
    ax.axhline(95, color=C["red"], ls="-.", lw=1)
    ax.set_xticks(xs, ["weak\nscaling", "nominal", "strong\nscaling"]); ax.set_ylim(40, 103); ax.set_ylabel("P(pass), % (STRESS set)")
    ax.set_title("(b) Named closure scenarios"); ax.legend(fontsize=7.4, loc="lower right")
    save(fig, F / "model_form" / "closure_stress_test.png")

