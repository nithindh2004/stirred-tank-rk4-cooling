"""stcool: jacketed stirred-tank cooling -- model, verification, optimisation and reliability analysis.

Modules
-------
model            energy balances, closures, RK4, exact solution, reduced-order model
uncertainty      actuator model, uncertain inputs, data sets, Monte Carlo evaluation, Wilson bounds
simulation       deterministic studies: convergence, stability, energy balance, reduced model
optimization     deterministic and chance-constrained optimisation, validation protocol
sensitivity      local and Sobol sensitivity
adaptive_policy  simulated measure-then-set policy with soft-sensor estimation
plotting         shared figure style
crosscheck       Python vs MATLAB/Octave comparison
pipeline         end-to-end workflow used by run_all.py
"""
__version__ = "3.0.0"
