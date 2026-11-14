# Changelog

## [3.0.0] — V3

### Corrected
- **Coolant flow above the pump limit.** Flow uncertainty now passes through a saturated actuator, `ṁ = min(ṁ_max, ṁ_cmd(1 + 0.05z))`. Half of V2's simulated flows at a 0.50 kg/s command had exceeded the pump capacity (up to 0.62 kg/s).
- **Pump energy** is computed at the delivered flow. The chance-constrained objective is the expected energy, written as the time integral of electrical power.
- **Validation independence.** Separate optimisation (OPT), calibration (CAL), validation (VAL), reserve (VAL2) and stress (STRESS) sets, each with one purpose. One-sided 95 % lower-bound acceptance rule. Every rounding rule is written into the code before validation. V2's rounding loop that consulted the validation set was removed.
- **No post-validation modification.** A development run of V3 rounded the adaptive-policy margin up after an Octave replication fell short of the acceptance rule, then re-validated on a third set. The final release declares that rounding in advance, validates once, and replaces all validation seeds examined during development (20261110, 20261111, 20261115 → 20261120, 20261121).
- **Adaptive-policy figure.** Speed-capped batches (1.27 %) are shown separately from high-energy batches (2.34 % above 3.2 kWh). V2 had labelled all batches above 3.2 kWh as speed-capped.
- **Local sensitivity above the pump limit.** Upward flow perturbations are now saturated at 0.50 kg/s.

### Changed results
- Recommended fixed setting: **295 rpm** at a commanded 0.50 kg/s (V2: 285 rpm, which reaches only 94.73 % under the corrected model). Expected energy 1.076 kWh; P(pass) on VAL 95.39 % [95.29, 95.48], one-sided lower bound 95.31 %.
- Adaptive policy: 200 rpm for 600 s, margin calibrated at 0.186 K and rounded to 0.20 K by the pre-declared rule. P(pass) on VAL 95.97 % [95.88, 96.05]; mean energy 0.496 kWh (54 % below the fixed setting); 1.27 % of batches speed-capped. Simulation only.

### Added
- Model-form stress test: admissible exponent region from film-coefficient scaling, weak/nominal/strong closures, an exponent grid, the adaptive policy under each closure.
- Jacket flow-pattern check (tanks-in-series RK4): the well-mixed jacket is the conservative limit.
- Parameter classification and equipment/process feasibility checks.
- RK4 stability limit and test; applicability table for the reduced-order model.
- Repository restructure (`src/python/stcool`, `src/matlab/supporting_functions`, `data/`, `figures/`, `outputs/`, `reports/`); 41 tests; a generated-values pipeline that feeds both reports; a test that keeps the README headline numbers in step with the outputs.
- Detailed report and a 3-page two-column paper, each as PDF, DOCX and LaTeX source. The detailed report's layout follows the author's CHE453 (IIT Kanpur) course-report template; this is a formatting choice only. The project itself is the ESC113M extension, and no CHE453 material is included.

### Release preparation (first public release; no change to any scientific result)
- **Licensing.** MIT for the code, CC BY 4.0 for reports, figures, documentation and generated data (`LICENSE`, `LICENSE-CC-BY-4.0.txt`). The original ESC113M Group 03 materials in `docs/original_esc113m/` are excluded from both licences and redistributed with the group's permission (`LICENSING.md`). `CITATION.cff` updated with affiliation, licences and the original project as a reference.
- **Runtime.** Phase 1 of the adaptive policy (simulated observation and soft-sensor estimate) does not depend on the safety margin, so `calibrate_delta` now simulates it once and reuses it for every trial margin instead of recomputing it 19 times. Same seeds, sample sizes, data sets and random-number order; every CSV, figure and report table is byte-identical. `run_all.py` went from 51 s to 36 s on the reference machine.
- **Report metadata.** Title, author, subject, keywords and institution set in both PDFs (`\hypersetup`) and both DOCX files (core properties). Rebuilt PDFs are pixel-identical to the previous ones at 150 dpi, and their extracted text is identical.
- **CHE453 reference clarified.** The detailed report's layout follows the author's CHE453 course-report template, at the author's request; the wording now states that this is a formatting choice only.
- **Figure rendering pinned.** `apply_style()` now sets `figure.dpi` and `text.hinting` to the values used for the published figures, so a local `matplotlibrc` cannot change the rendering. The committed figures are unchanged; a fresh environment now reproduces them pixel for pixel.
- **README** restructured; it states the simulation-only scope, the assumed distributions and closure, the auxiliary-energy scope and the measured runtime.
- **Review corrections.** `CITATION.cff` no longer lists two licences, because the Citation File Format reads a licence list as alternatives (OR); the scope of each licence stays in `LICENSING.md` and the CFF abstract. Detailed report: "real batches" changed to "simulated batches". Paper conclusion and README: the one-command claim now covers only the Python simulation results and figures; the Octave cross-check and report builds are separate commands.
- **Tests.** Seven release checks added (`tests/test_release.py`); no existing test changed. `pypdf` added as a test dependency.

## [2.0.0] — V2
- Exact matrix-exponential solution, series-resistance reduced model, continuous nominal optimum, chance-constrained optimisation, reliability–energy frontier, Sobol indices, adaptive policy, model-form grid, MATLAB code executed in Octave, automated tests.

## [1.0.0] — V1
- Coupled liquid–jacket model with operating-dependent U, power model, grid optimum and Monte Carlo study, built on the ESC113M single-state RK4 model.
