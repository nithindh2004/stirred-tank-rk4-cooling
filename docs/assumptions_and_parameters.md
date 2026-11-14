# Assumptions and parameters

Every parameter is classified by origin. **None of the equipment values is an industrial specification**; the values marked as assumptions are screening values for this case study.

Classification key:

* **carried over from ESC113M**: the original course model;
* **literature property value / literature-supported**: standard property or correlation values;
* **case-study assumption**: introduced here to make the problem well-posed;
* **needs calibration**: should be replaced by experimental or vendor data before any real use.

| Symbol | Value | Unit | Description | Classification |
|---|---|---|---|---|
| V | 0.5 | m3 | liquid volume | carried over from ESC113M |
| rho_L | 1000.0 | kg/m3 | liquid density | carried over from ESC113M |
| Cp_L | 4180.0 | J/(kg K) | liquid heat capacity | carried over from ESC113M |
| T_0 | 370.0 | K | initial liquid temperature | carried over from ESC113M |
| T_in | 300.0 | K | coolant inlet temperature (original fixed coolant) | carried over from ESC113M |
| A | 2.0 | m2 | heat-transfer area | carried over from ESC113M |
| U_ref | 150.0 | W/(m2 K) | overall coefficient at the reference point | carried over; needs calibration |
| t_f | 1800.0 | s | batch horizon | carried over from ESC113M |
| Cp_C | 4180.0 | J/(kg K) | coolant heat capacity (water) | literature property value |
| rho_C | 1000.0 | kg/m3 | coolant density (water) | literature property value |
| N_P | 5.0 | - | impeller power number (turbulent Rushton-type turbine) | literature-supported; vendor calibration |
| M_J | 10.0 | kg | effective well-mixed jacket inventory | case-study assumption; needs geometry |
| N_ref | 300.0 | rpm | reference agitator speed | case-study assumption |
| mdot_ref | 0.2 | kg/s | reference coolant flow | case-study assumption |
| e_N | 0.25 | - | U exponent on speed | case-study assumption; needs calibration |
| e_M | 0.33 | - | U exponent on coolant flow | case-study assumption; needs calibration |
| D | 0.3 | m | impeller diameter | case-study assumption |
| eta_ag | 0.75 | - | agitator drive efficiency | case-study assumption; vendor calibration |
| dP_ref | 50000.0 | Pa | coolant-loop pressure drop at mdot_ref | case-study assumption; vendor calibration |
| eta_p | 0.65 | - | pump efficiency | case-study assumption; vendor calibration |
| m_max | 0.5 | kg/s | maximum deliverable coolant flow (pump limit) | case-study assumption (equipment limit) |
| N range | 100-600 | rpm | agitator speed bounds | case-study assumption |
| T_target | 355.0 | K | liquid temperature required at t_f | illustrative design constraint |

## Modelling assumptions

1. Perfectly mixed liquid; closed batch; no reaction, evaporation or heat loss.
2. Jacket coolant is one well-mixed inventory. This is the conservative limit of jacket mixing; a plug-flow jacket would raise the effective conductance by 9–14 %.
3. Constant properties. For water over 300–370 K the density changes by up to 4 % and the heat capacity by less than 1 %.
4. One overall coefficient U on a fixed area; wall heat capacity neglected.
5. Power-law closure for U with assumed exponents (stress-tested as scenarios).
6. Saturated flow actuator: proportional 5 % loop error, hard limit 0.50 kg/s.
7. Auxiliary electricity only (agitator + pump). The pump is modelled as variable-speed with constant efficiency; a throttled fixed-speed pump would use more energy at low flow.
8. Uncertain inputs are independent and normal with **assumed** spreads (`data/inputs/uncertainty_model.csv`).

## Equipment and process checks (indicative)

| Check | Value | Unit | Comment |
|---|---|---|---|
| Tank diameter for H = T (0.50 m3) | 0.8603 | m | geometry implied by V; no tank drawing exists |
| Wetted wall + base area at H = T | 2.906 | m2 | A = 2.0 m2 is of the same order |
| Impeller-to-tank diameter ratio D/T | 0.3487 | - | within the commonly used 0.3-0.5 range |
| Impeller Reynolds number at 120 rpm (water-like) | 2.25e+05 | - | turbulent: constant N_P applicable |
| Shaft power per volume at recommended speed | 2.888 | kW/m3 | high-intensity duty; check with vendor |
| Loop pressure drop at 0.50 kg/s | 3.125 | bar | extrapolated from assumed 0.5 bar at 0.2 kg/s |
| Pump flow turndown modelled (0.05-0.50 kg/s) | 10 | :1 | model assumes variable-speed pump, constant efficiency |
| Jacket gap implied by M_J over A | 5 | mm | thin; see isolated-inventory study |
| Plain annular jacket Re at 0.50 kg/s, Re = 2 m/(pi T mu) | 462.5 | - | laminar: assumed U needs a baffled/spiral jacket |
| Coolant outlet maximum (baseline peak T_J) | 318.1 | K | about 45 C; confirm against site cooling-water limits |
| Liquid density range 300-370 K (water) | 960 | kg/m3 | about 996 to 960 kg/m3; constant 1000 used (<= 4 %) |
| Liquid Cp range 300-370 K (water) | 4.21 | kJ/(kg K) | about 4.18 to 4.21; constant 4.18 used (< 1 %) |

## Jacket geometry
The 10 kg inventory implies a 5 mm gap. The isolated inventory study (`data/model_form/jacket_inventory_isolated.csv`) shows that the thermal inventory alone hardly matters, because the jacket time constant is short compared with the batch. That study does **not** show that jacket geometry is unimportant: a different gap or flow path also changes the coolant velocity, the pressure drop and the coolant-side coefficient, and therefore U, e_M and the pump energy. For a plain annular jacket with axial flow, Re = 2ṁ/(πTμ) is about 460 at 0.50 kg/s (laminar), so the assumed coolant-side heat transfer would in practice need a baffled or spiral jacket.
