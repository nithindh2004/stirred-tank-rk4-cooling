# MATLAB / GNU Octave implementation

`run_all_matlab.m` is an independent implementation of the V3 workflow, using the functions in `supporting_functions/`. It needs no toolboxes and runs unchanged in MATLAB R2019b+ and GNU Octave 8+.

```
cd src/matlab
octave --no-gui --eval run_all_matlab > ../../outputs/octave_run_log.txt     % or: run_all_matlab  (MATLAB)
```

It writes CSV files to `data/verification/octave/` and prints the values that `src/python/stcool/crosscheck.py` compares with the Python reference. Monte Carlo values differ from Python by sampling error only, because the two languages use different random generators.

Tested with: GNU Octave 8.4.0. It has not been run in MATLAB.
