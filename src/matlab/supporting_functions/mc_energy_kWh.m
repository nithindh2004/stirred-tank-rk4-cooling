function E = mc_energy_kWh(rpm, mcmd, S, p)
%MC_ENERGY_KWH  Per-batch electricity with the pump at the actual flow.
E = batch_energy_kWh(rpm, actual_flow(mcmd, S.z, p), p);
end
