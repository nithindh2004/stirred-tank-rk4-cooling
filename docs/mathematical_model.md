# Mathematical model

## Energy balances
```
M_L Cp_L dT_L/dt = -U A (T_L - T_J)                                  (1)
M_J Cp_C dT_J/dt =  U A (T_L - T_J) + ṁ Cp_C (T_in - T_J)            (2)
T_L(0) = 370 K,  T_J(0) = T_in = 300 K
```
Adding (1) and (2) and integrating gives the first-law check used in the verification:
`M_L Cp_L [T_0 - T_L(t)] = ∫ ṁ Cp_C (T_J - T_in) dt + M_J Cp_C [T_J(t) - T_in]`.
Because `T_J < T_L` at all times, `T_L` decreases monotonically, so "reach 355 K within 1800 s" is equivalent to `T_L(1800 s) <= 355 K`.

## Heat-transfer closure (assumed)
`U = U_ref (N/N_ref)^eN (ṁ/ṁ_ref)^eM`, with U_ref = 150 W/(m² K) at 300 rpm and 0.20 kg/s, eN = 0.25 and eM = 0.33.
Physical bound: writing `1/U = 1/h_i + R_w + 1/h_o` with `h_i ∝ N^(2/3)` and `h_o ∝ ṁ^b`, `b <= 0.8`, the local exponents satisfy `eN/(2/3) + eM/0.8 <= 1`.

## Coolant actuator
`ṁ = min[ṁ_max, max(0, ṁ_cmd (1 + 0.05 z))]`, `z ~ N(0, 1)`, `ṁ_max = 0.50 kg/s`.
This is a saturated flow loop. At `ṁ_cmd = ṁ_max`, 50 % of batches receive exactly the limit and the mean delivered flow is 0.98 `ṁ_max`.

## Electricity (auxiliary only)
```
P_ag   = N_P ρ_L n^3 D^5 / η_ag            (n = N/60 rev/s)
ΔP     = ΔP_ref (ṁ/ṁ_ref)^2
P_pump = ΔP (ṁ/ρ_C) / η_p                   (delivered flow)
E      = (1/3.6e6) ∫ [P_ag(t) + P_pump(t)] dt   kWh
```

## Exact solution
With `x = (T_L - T_in, T_J - T_in)`, `x' = M x`, `M = [[-a, a], [b, -(b+f)]]`, `a = UA/(M_L Cp_L)`, `b = UA/(M_J Cp_C)`, `f = ṁ/M_J`.
The eigenvalues are real, negative and distinct, and
`e^{Mt} = [e^{λ1 t}(M - λ2 I) - e^{λ2 t}(M - λ1 I)] / (λ1 - λ2)`.

## Reduced-order model
Quasi-steady jacket: `1/UA_eff = 1/UA + 1/(ṁ Cp_C)` and `T_L(t) = T_in + (T_0 - T_in) exp(-UA_eff t / (M_L Cp_L))`.
Required conductance: `UA_eff,req = M_L Cp_L ln[(T_0 - T_in)/(T_target - T_in)] / t_f = 280.0 W/K`.
Plug-flow jacket limit: `UA_eff = ṁ Cp_C (1 - exp(-UA/(ṁ Cp_C)))`, which is always at least the mixed value.

## Optimisation problems
Deterministic: `min E(N, ṁ_cmd)` s.t. `T_L(t_f) <= 355 K`, `100 <= N <= 600 rpm`, `0.05 <= ṁ_cmd <= 0.50 kg/s`.
Chance-constrained: `min E[E(N, ṁ)]` s.t. `P[T_L(t_f) <= 355 K] >= 0.95`, enforced on the optimisation sample through the one-sided 95 % Wilson lower bound.
