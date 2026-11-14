# Reproducibility

## One-command workflow
```bash
pip install -r requirements.txt
python run_all.py          # writes data/, figures/, outputs/ and reports/build/generated/ (about 36 s)
python -m pytest -q        # 48 tests: 41 scientific tests + 7 release checks
```
`run_all.py` is deterministic: every random draw comes from a fixed seed (`data/inputs/datasets.csv`; Sobol seed 20261113). Re-running from a clean checkout reproduces every CSV byte-for-byte on the same software versions (`outputs/run_metadata.json`).

## Independent implementation
`src/matlab/run_all_matlab.m` (with `supporting_functions/`) re-implements the model, verification, chance-constrained optimisation and adaptive-policy validation. It runs unchanged in MATLAB R2019b+ and GNU Octave 8+. **It has been executed in GNU Octave 8.4.0 only.**
```bash
octave --no-gui --eval "addpath('src/matlab'); run_all_matlab" > outputs/octave_run_log.txt
python run_all.py                 # rebuilds data/verification/python_vs_octave.csv from the log
```
The MATLAB code writes its own CSV files to `data/verification/octave/`. The adaptive-policy configuration (200 rpm, 600 s) is the one selected by the Python calibration scan; a test checks that the two agree.
Deterministic results agree to the printed precision. Monte Carlo results agree within sampling error, because MATLAB/Octave (Mersenne Twister) and NumPy (PCG64) draw different samples even with the same seed numbers.

## Reports
```bash
python reports/build/build_reports.py              # both reports
python reports/build/build_reports.py scientific_paper
```
Requirements: TeX Live with `pdflatex`, IEEEtran, booktabs, fancyhdr, enumitem, caption; pandoc ≥ 3; python-docx. `pypdf` is used only by the metadata test. Every number in both reports is a macro from `reports/build/generated/values.tex` and every table is generated from `data/`, so the reports cannot drift from the code. The DOCX files are converted from the same LaTeX sources; their equation, table and figure numbers are taken from the LaTeX `.aux` file.

## Software used for the published results
Python 3.13, NumPy 2.5, SciPy 1.18, pandas 3.0, Matplotlib 3.11, pytest (see `outputs/run_metadata.json` for exact versions); GNU Octave 8.4.0; TeX Live (pdflatex, IEEEtran); pandoc 3.1.3; python-docx 1.2.

## Runtime
Measured from a clean copy with every generated file deleted, on a 2-vCPU Intel Xeon (2.1 GHz) Linux machine:

| Stage | Time |
|---|---|
| Octave script | 11 s |
| `python run_all.py` | 36 s |
| `python -m pytest -q` | 4 s |
| `python reports/build/build_reports.py` | 7.5 s |
| Complete chain | 58 s |

Inside `run_all.py` the largest stages are the adaptive-policy calibration and validation (11 s), the model-form stress test (11 s) and the reliability–energy frontier (7.5 s); `run_all.py` prints the elapsed time at the start of each stage. Before the release preparation `run_all.py` took 51 s: the margin calibration re-simulated the margin-independent observation phase for every trial margin. It now simulates that phase once, with the same random numbers in the same order, and all outputs are unchanged.
