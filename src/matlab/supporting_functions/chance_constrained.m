function [rb, mb, Eb] = chance_constrained(alpha, S, p)
%CHANCE_CONSTRAINED  Minimise expected energy s.t. the lower-bound pass rule on dataset S.
preq = required_pass_fraction(alpha, numel(S.T0));
pf = @(r, m) mean(mc_final_TL(r, m, S, p) <= p.Ttarget);
Eb = inf; rb = NaN; mb = NaN;
for m = round((p.mMin:0.005:p.mMax)*1e6)/1e6
    if pf(p.rpmMax, m) < preq, continue; end
    if pf(p.rpmMin, m) >= preq
        r = p.rpmMin;
    else
        lo = p.rpmMin; hi = p.rpmMax;
        while hi - lo > 1e-3
            mid = 0.5*(lo+hi);
            if pf(mid, m) >= preq, hi = mid; else, lo = mid; end
        end
        r = hi;
    end
    E = mean(mc_energy_kWh(r, m, S, p));
    if E < Eb, Eb = E; rb = r; mb = m; end
end
end
