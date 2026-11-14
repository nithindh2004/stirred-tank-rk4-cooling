function P = pump_power(mdot, p)
%PUMP_POWER  Electrical pump power (W) at actual flow: dP*Q/eta, dP = dPref (m/mref)^2.
P = p.dpRef*(mdot/p.mref).^2 .* (mdot/p.rhoC)/p.etaP;
end
