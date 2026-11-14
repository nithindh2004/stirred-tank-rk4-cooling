# Licensing and attribution

This repository combines material from two sources with different authors. Each part is licensed separately.

| Part | Authors | Terms |
|---|---|---|
| Source code: `run_all.py`, `src/`, `tests/`, `reports/build/*.py` | Nithin D H (V3 extension) | MIT License — see [`LICENSE`](LICENSE) |
| Reports, figures, documentation and simulation-generated data: `reports/` (PDF, DOCX, LaTeX sources), `figures/`, `data/`, `outputs/`, `docs/` (except `docs/original_esc113m/`), `README.md`, `CHANGELOG.md` | Nithin D H (V3 extension) | Creative Commons Attribution 4.0 International (CC BY 4.0) — see [`LICENSE-CC-BY-4.0.txt`](LICENSE-CC-BY-4.0.txt) |
| Original ESC113M course project: everything in `docs/original_esc113m/` except its `README.md` | ESC113M Group 03, IIT Kanpur | **Not covered by either licence above.** Copyright remains with the members of Group 03. It is included and redistributed here with their permission, unchanged, as part of this repository. Any other use requires their permission. |

## What came from where

* **Original group project (ESC113M Group 03).** The single-state energy balance for a perfectly mixed liquid cooled by a jacket at a constant 300 K, its RK4 solution with Δt = 10 s, the comparison with the analytical solution, the base parameter set (V = 0.5 m³, ρ = 1000 kg/m³, Cp = 4180 J/(kg K), U = 150 W/(m² K), A = 2 m², T₀ = 370 K, T_c = 300 K, 1800 s horizon), the original MATLAB script and the original report.
* **V3 extension (Nithin D H).** The coupled liquid–jacket model, the operating-dependent heat-transfer closure, the actuator and electricity models, the exact-solution verification, the sensitivity, optimisation, reliability and adaptive-policy analyses, the model-form stress test, the Python and MATLAB/Octave implementations of the extension, the tests, the figures, the generated data and both V3 reports. The V3 code re-implements the original model independently in order to reproduce its result; it does not copy the original script.

## Attribution when reusing

When reusing the code or reports, credit "Nithin D H, *Jacketed Stirred-Tank Cooling: Dynamic Modelling, Reliability-Based Optimisation and Adaptive Operation*" (see `CITATION.cff`). When referring to the base model, also credit ESC113M Group 03, IIT Kanpur.

The original report identifies its authors only as "Group – 03". The individual members' names are not recorded in the original files and are therefore not listed here.

Third-party references cited in the reports are the property of their respective authors and publishers and are cited, not reproduced.
