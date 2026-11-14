"""Seeds, dataset separation and agreement between stored outputs and the code."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from stcool.uncertainty import DATASETS, dataset, NOISE_SEEDS
from stcool.optimization import chance_constrained_optimum, practical_setpoint, validate

ROOT = Path(__file__).resolve().parents[1]


def test_datasets_are_reproducible():
    a, b = dataset("OPT"), dataset("OPT")
    assert all(np.array_equal(a[k], b[k]) for k in ("Uref", "MJ", "z_flow", "T0"))


def test_datasets_use_distinct_seeds_and_do_not_overlap():
    seeds = [d["seed"] for d in DATASETS.values()]
    assert len(set(seeds)) == len(seeds)
    assert len(set(NOISE_SEEDS.values())) == len(NOISE_SEEDS)
    opt, val = dataset("OPT"), dataset("VAL")
    assert len(np.intersect1d(opt["T0"], val["T0"])) == 0


@pytest.mark.skipif(not (ROOT / "outputs" / "summary_results.csv").exists(), reason="run `python run_all.py` first")
def test_stored_recommendation_matches_code():
    s = pd.read_csv(ROOT / "outputs" / "summary_results.csv").set_index("metric")["value"]
    OPT, VAL = dataset("OPT"), dataset("VAL")
    best, _ = chance_constrained_optimum(0.95, OPT)
    sp = practical_setpoint(best, 0.95, OPT)
    assert sp["rpm"] == pytest.approx(float(s["Recommended agitator speed"]))
    v = validate(sp["rpm"], sp["m"], VAL)
    assert v["pass_frac"] == pytest.approx(float(s["Recommended setting P(pass) (VAL)"]), abs=1e-12)
    meta = json.loads((ROOT / "outputs" / "run_metadata.json").read_text())
    assert meta["datasets"]["VAL"]["seed"] == DATASETS["VAL"]["seed"]


@pytest.mark.skipif(not (ROOT / "outputs" / "summary_results.csv").exists(), reason="run `python run_all.py` first")
def test_readme_headline_numbers_match_outputs():
    """The hand-written README must quote the numbers that the final run produced."""
    s = pd.read_csv(ROOT / "outputs" / "summary_results.csv").set_index("metric")["value"].astype(float)
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    expected = [
        f"{s['Recommended agitator speed']:.0f} rpm at a commanded {s['Recommended commanded coolant flow']:.2f} kg/s",
        f"{s['Recommended setting expected energy (VAL)']:.3f} kWh per batch",
        f"P(pass) = {100 * s['Recommended setting P(pass) (VAL)']:.2f} %",
        f"one-sided lower bound {100 * s['Recommended setting P(pass) one-sided 95% lower bound (VAL)']:.2f} %",
        f"**{100 * s['Nominal optimum P(pass) (VAL)']:.1f} %**",
        f"{100 * s['Adaptive policy P(pass) (independent)']:.2f} % reliability",
        f"{s['Adaptive policy mean energy (independent)']:.3f} kWh mean",
        f"{100 * s['Adaptive policy speed-capped batches']:.2f} % of batches hit the 600 rpm cap",
        f"{s['Nominal (deterministic) optimum speed']:.1f} rpm",
    ]
    missing = [e for e in expected if e not in readme]
    assert not missing, missing
