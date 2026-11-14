function U = U_coef(rpm, mdot, p, Uref)
%U_COEF  Overall heat-transfer coefficient (W/m2/K) at actual coolant flow mdot.
U = Uref .* (rpm/p.Nref).^p.eN .* (mdot/p.mref).^p.eM;
end
