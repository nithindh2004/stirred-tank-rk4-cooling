function [pf, E, capped, UAerr] = adaptive_policy(S, rpm1, mcmd, tobs, delta, p, noiseSeed, sigma)
%ADAPTIVE_POLICY  Measure-then-set policy (see src/python/stcool/adaptive_policy.py).
rng(noiseSeed, 'twister');
n = numel(S.Uref); md = actual_flow(mcmd, S.z, p);
UA1 = U_coef(rpm1, md, p, S.Uref)*p.A;
tm = 0:10:tobs;
[TL, TJ] = propagate_exact(S.T0, p.Tin, UA1, md, S.MJ, tm, p);
yL = TL + sigma*randn(size(TL)); yJ = TJ + sigma*randn(size(TJ));
w = tm >= 120; tw = tm(w); tc = tw - mean(tw);
mL = mean(yL(:,w), 2); sL = (yL(:,w)*tc')/sum(tc.^2);
mJ = mean(yJ(:,w), 2); sJ = (yJ(:,w)*tc')/sum(tc.^2);
Q = -p.ML*p.CpL*sL; UAh = Q./(mL - mJ); mCph = (Q - p.MJ*p.CpC*sJ)./(mJ - p.Tin);
TLo = mL + sL*(tobs - mean(tw)); TJo = mJ + sJ*(tobs - mean(tw));
goal = p.Ttarget - delta;
pred = @(r) propagate_exact(TLo, TJo, UAh.*(r/rpm1).^p.eN, mCph/p.CpC, p.MJ, p.tf - tobs, p);
lo = p.rpmMin*ones(n,1); hi = p.rpmMax*ones(n,1);
okMin = pred(lo) <= goal; okMax = pred(hi) <= goal;
for k = 1:25
    mid = 0.5*(lo+hi); good = pred(mid) <= goal;
    hi(good) = mid(good); lo(~good) = mid(~good);
end
rpm2 = hi; rpm2(okMin) = p.rpmMin; rpm2(~okMax) = p.rpmMax;
TLf = propagate_exact(TL(:,end), TJ(:,end), U_coef(rpm2, md, p, S.Uref)*p.A, md, S.MJ, p.tf - tobs, p);
E = (agitator_power(rpm1, p)*tobs + agitator_power(rpm2, p)*(p.tf - tobs) + pump_power(md, p)*p.tf)/3.6e6;
pf = mean(TLf <= p.Ttarget); capped = ~okMax; UAerr = UAh./UA1 - 1;
end
