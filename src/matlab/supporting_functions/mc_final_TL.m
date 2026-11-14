function x = mc_final_TL(rpm, mcmd, S, p)
%MC_FINAL_TL  T_L(t_f) for every sample, using the saturated actual flow.
md = actual_flow(mcmd, S.z, p);
x = final_TL(rpm, md, p, S.Uref, S.MJ, S.T0);
end
