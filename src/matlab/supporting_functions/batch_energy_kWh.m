function E = batch_energy_kWh(rpm, mdot, p)
%BATCH_ENERGY_KWH  Agitator + pump electricity over the batch (constant speed, actual flow).
E = (agitator_power(rpm, p) + pump_power(mdot, p))*p.tf/3.6e6;
end
