function TL = final_TL(rpm, mdot, p, Uref, MJ, T0)
%FINAL_TL  Exact T_L(t_f) for constant speed and actual flow (arrays allowed).
TL = propagate_exact(T0, p.Tin, U_coef(rpm, mdot, p, Uref)*p.A, mdot, MJ, p.tf, p);
end
