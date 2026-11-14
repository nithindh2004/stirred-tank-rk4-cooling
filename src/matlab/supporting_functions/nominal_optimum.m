function [r, m, E] = nominal_optimum(p)
%NOMINAL_OPTIMUM  Deterministic minimum-energy point on the T_L(t_f) = T_target boundary.
ms = linspace(p.mMin, p.mMax, 901); E = inf; r = NaN; m = NaN;
for q = 1:numel(ms)
    rb = boundary_rpm(ms(q), p);
    if ~isnan(rb)
        Eq = batch_energy_kWh(rb, ms(q), p);
        if Eq < E, E = Eq; r = rb; m = ms(q); end
    end
end
end

function r = boundary_rpm(m, p)
g = @(r) final_TL(r, m, p, p.Uref, p.MJ, p.T0) - p.Ttarget;
if g(p.rpmMax) > 0, r = NaN; return; end
if g(p.rpmMin) <= 0, r = p.rpmMin; return; end
lo = p.rpmMin; hi = p.rpmMax;
for k = 1:60, mid = 0.5*(lo+hi); if g(mid) <= 0, hi = mid; else, lo = mid; end, end
r = hi;
end
