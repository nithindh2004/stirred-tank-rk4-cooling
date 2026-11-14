function [TL, TJ] = propagate_exact(TL0, TJ0, UA, mdot, MJ, t, p)
%PROPAGATE_EXACT  Exact solution of the linear two-state model (implicit expansion).
a = UA/(p.ML*p.CpL); b = UA./(MJ*p.CpC); f = mdot./MJ;
tr = -(a+b+f); disc = sqrt(tr.^2 - 4*a.*f);
l1 = 0.5*(tr+disc); l2 = 0.5*(tr-disc);
x1 = TL0 - p.Tin; x2 = TJ0 - p.Tin; e1 = exp(l1.*t); e2 = exp(l2.*t); d = l1 - l2;
TL = p.Tin + (e1.*((-a-l2).*x1 + a.*x2) - e2.*((-a-l1).*x1 + a.*x2))./d;
TJ = p.Tin + (e1.*(b.*x1 + (-(b+f)-l2).*x2) - e2.*(b.*x1 + (-(b+f)-l1).*x2))./d;
end
