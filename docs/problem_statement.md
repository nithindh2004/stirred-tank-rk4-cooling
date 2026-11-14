# Problem statement and objectives

## Design question
Which agitator speed `N` and commanded coolant flow `ṁ_cmd` cool 500 kg of a water-like liquid in a jacketed, well-stirred batch tank from 370 K to **355 K or below within 1800 s** at the lowest auxiliary electricity? And in what fraction of simulated batches does it succeed when the heat-transfer coefficient, jacket inventory, delivered coolant flow and initial temperature are uncertain?

The 355 K / 1800 s target is an **illustrative case-study constraint**, chosen close to the 30-minute result of the original model. It is not a plant specification. The required reliability is 95 %, demonstrated on data that played no part in choosing the design.

## Objectives
1. Extend the original single-state model to coupled liquid and jacket-coolant energy balances with a finite, continuously supplied coolant inventory and operating-condition-dependent heat transfer.
2. Keep RK4 as the principal solver and verify it against the exact matrix-exponential solution and the first law.
3. Find the minimum-energy operating point, first deterministically and then with a 95 % reliability requirement under parameter uncertainty and a physical pump limit (0.50 kg/s).
4. Quantify which uncertainties matter (local and Sobol sensitivity).
5. Assess a simulated adaptive measure-then-set policy that estimates the heat-transfer conductance from early temperature data.
6. Test how far the conclusions depend on the model structure (heat-transfer closure, jacket mixing).

## Scope boundary
Inside: the tank liquid, the jacket coolant, the agitator, the coolant pump and flow loop. Outside: the chiller or cooling tower, water treatment, control hardware, and any reaction or phase change in the liquid.
