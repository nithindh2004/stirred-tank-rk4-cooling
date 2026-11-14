function P = agitator_power(rpm, p)
%AGITATOR_POWER  Electrical agitator power (W): Np rho n^3 D^5 / eta, n in rev/s.
P = p.Np*p.rhoL*(rpm/60).^3*p.D^5/p.etaAg;
end
