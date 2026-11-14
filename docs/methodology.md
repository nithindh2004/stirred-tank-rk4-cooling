# Methodology

## Numerical methods
* **RK4** (classical, Δt = 1 s) integrates every deterministic trajectory. The largest stable step is 2.785 τ_J (47 s at the recommended point). Convergence is measured against the exact solution (observed order 4.22).
* The **exact solution** is used for verification and to evaluate `T_L(1800 s)` for Monte Carlo samples (up to about 10⁸ evaluations); the difference from RK4 is 3.2 × 10⁻⁸ K.
* The **tanks-in-series jacket** (model-form check) is integrated with a vectorised RK4.
* The **first law** is checked with Simpson's rule (residual 0.0044 J out of 25.5 MJ).

## Monte Carlo design and data separation
| Set | Seed | n | Single use |
|---|---|---|---|
| OPT | 20261109 | 20,000 | choosing fixed settings (common random numbers) |
| CAL | 20261112 | 20,000 | calibrating the margin and selecting the adaptive-policy configuration |
| VAL | 20261120 | 200,000 | validating each selected design and the adaptive policy, once |
| VAL2 | 20261121 | 200,000 | reserved for a design that fails VAL and is modified by the pre-declared rule (not needed) |
| STRESS | 20261114 | 50,000 | robustness and model-form tests (evaluation only) |

VAL and VAL2 are final-release seeds. The validation seeds examined during development (20261110, 20261111, 20261115) were retired, so no reported validation result was seen before the final run. Sensor-noise seeds: CAL 31, VAL 40, VAL2 41, STRESS 34. Sobol sampling uses seed 20261113.

Acceptance rule: the one-sided 95 % Wilson lower bound of P(pass) must be at least 0.95. On OPT and CAL (n = 20,000) this means a sample pass fraction of at least 95.25 %.

## Fixed-setting protocol
1. Chance-constrained optimum on OPT: bisection on speed at each commanded flow (the pass fraction is monotone in speed under common random numbers), then minimisation of expected energy over flow. Result: 293.1 rpm at 0.50 kg/s.
2. Practical setpoint by a pre-declared rule: speed rounded **up** to the next 5 rpm (295 rpm), re-checked on OPT only.
3. One validation on VAL: P(pass) = 95.39 % (95 % CI 95.29–95.48 %), lower bound 95.31 %, accepted. Had it failed, the pre-declared modification (+5 rpm) would have been validated on VAL2.

## Adaptive policy protocol
1. For each of 12 (N1, t_obs) configurations, the margin δ is calibrated on CAL as the smallest value meeting the same lower-bound rule, then rounded **up** to the next 0.05 K by a rule fixed in the code before any validation.
2. The configuration with the lowest CAL mean energy at its rounded margin is selected: 200 rpm for 600 s, δ calibrated at 0.186 K and rounded to 0.20 K.
3. One validation on VAL: P(pass) = 95.97 % (95 % CI 95.88–96.05 %), lower bound 95.89 %, accepted. Had it failed, δ + 0.05 K would have been validated on VAL2.
4. The independent Octave implementation calibrates its own margin on its own samples (0.183 K), applies the same rounding (0.20 K) and validates at 95.93 %.

## Sensitivity
* Local: v1 steps (±20 %, ±1 K) and ±1 standard deviation of each uncertain input.
* Global: Sobol first-order (Saltelli 2010) and total (Jansen 1999) indices, N = 32,768, 500 bootstrap resamples.

## Model form
Named closure scenarios obtained by varying the wall-plus-fouling resistance share (0.60 / 0.21 / 0), an (eN, eM) grid with the admissibility bound, a jacket split into 1–20 cells, and an isolated jacket-inventory study. No probability distributions are assigned to structural assumptions.
