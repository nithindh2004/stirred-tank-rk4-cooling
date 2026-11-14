"""Compare the Python results with the independent MATLAB/Octave run (outputs/octave_run_log.txt).

Deterministic quantities must agree to the printed precision. Monte Carlo
quantities are compared as two independent estimates: MATLAB/Octave's
Mersenne Twister and NumPy's PCG64 draw different samples from the same
distributions, so they agree only to within sampling error.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


def build(root: Path, py: dict):
    log_path = root / "outputs" / "octave_run_log.txt"
    if not log_path.exists():
        return None
    log = log_path.read_text()
    runtime = re.search(r"Runtime: ([^\n]*)", log)

    def num(pattern, group=1):
        m = re.search(pattern, log)
        return float(m.group(group)) if m else float("nan")

    rows = [
        ("Original model T_L(1800 s) [K]", py["T_orig"], num(r"Original fixed-coolant T_L\(1800\) = ([\d.]+)"), "deterministic"),
        ("Coupled baseline T_L(1800 s) [K]", py["T_base"], num(r"Coupled baseline T_L\(1800\) = ([\d.]+)"), "deterministic"),
        ("Coupled baseline T_J(1800 s) [K]", py["TJ_base"], num(r"T_J = ([\d.]+)"), "deterministic"),
        ("Max |RK4 - exact|, dt = 1 s [K]", py["err_dt1"], num(r"Max \|RK4\(dt=1\) - exact\| = ([\d.e+-]+)"), "deterministic"),
        ("Observed RK4 order", py["order"], num(r"order \(dt 2-30 s\) = ([\d.]+)"), "deterministic"),
        ("Energy residual, Simpson [J]", py["eb_simp"], num(r"Energy residual \(Simpson\) = ([-\d.]+)"), "deterministic"),
        ("Nominal optimum speed [rpm]", py["nom_rpm"], num(r"Nominal optimum: ([\d.]+) rpm"), "deterministic"),
        ("Nominal optimum energy [kWh]", py["nom_E"], num(r"Nominal optimum: .* E = ([\d.]+) kWh"), "deterministic"),
        ("Chance-constrained speed on OPT [rpm]", py["cc_rpm"], num(r"Chance-constrained optimum \(OPT\): ([\d.]+) rpm"), "Monte Carlo (independent samples)"),
        ("Recommended speed after rounding [rpm]", py["rec_rpm"], num(r"Recommended setting: (\d+) rpm"), "Monte Carlo (independent samples)"),
        ("Recommended P(pass) on VAL", py["rec_pass"], num(r"VAL P\(pass\) = ([\d.]+) \[[\d.]+, [\d.]+\], one-sided"), "Monte Carlo (independent samples)"),
        ("Recommended expected energy [kWh]", py["rec_E"], num(r"one-sided LB = [\d.]+, E = ([\d.]+) kWh"), "Monte Carlo (independent samples)"),
        ("V2 285 rpm P(pass) under V3 model", py["v2_pass"], num(r"285 rpm under V3 model: P\(pass\) = ([\d.]+)"), "Monte Carlo (independent samples)"),
        ("Adaptive calibrated margin on CAL [K]", py["ad_delta_cal"], num(r"Adaptive policy: calibrated delta ([\d.]+) K"), "Monte Carlo (independent samples)"),
        ("Adaptive margin after pre-declared rounding [K]", py["ad_delta"], num(r"rounded to ([\d.]+) K"), "Monte Carlo (independent samples)"),
        ("Adaptive P(pass) on VAL", py["ad_pass"], num(r"Adaptive policy: .*?VAL P\(pass\) = ([\d.]+)", ), "Monte Carlo (independent samples)"),
        ("Adaptive one-sided lower bound on VAL", py["ad_lb"], num(r"Adaptive policy: .*?LB ([\d.]+)"), "Monte Carlo (independent samples)"),
        ("Adaptive mean energy [kWh]", py["ad_E"], num(r"Adaptive policy: .*?mean E ([\d.]+) kWh"), "Monte Carlo (independent samples)"),
    ]
    df = pd.DataFrame(rows, columns=["quantity", "python", "octave", "type"])
    df["abs_diff"] = (df.python - df.octave).abs()
    df["runtime"] = runtime.group(1).strip() if runtime else "unknown"
    return df
