# Original ESC113M model

The starting point of this project: a single-state energy balance for a perfectly mixed liquid cooled by a jacket held at a constant 300 K,

    M_L Cp_L dT/dt = -U A (T - T_c),

solved with RK4 (Δt = 10 s) and compared with its analytical solution. Files:

* `ESC113M_Project_Code.m`: the original MATLAB script.
* `Original_ESC113M_Project_Report.pdf`: the original course report (group submission, Group 03).

These files are kept unchanged for reference.

**Authorship and licence.** Both files are the work of ESC113M Group 03, IIT Kanpur, and are redistributed here with the group members' permission. They are **not** covered by the repository's MIT or CC BY 4.0 licences; copyright remains with their authors (see `../../LICENSING.md`). The V3 extension built on them is described in the main `README.md`.

The extended model in `src/` reproduces the original result exactly (T(1800 s) = 354.06 K; RK4 vs analytical error about 5 × 10⁻¹³ K, see `tests/test_rk4.py`).
