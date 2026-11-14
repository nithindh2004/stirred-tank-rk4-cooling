function [lo, hi, lb1] = wilson(k, n)
%WILSON  Two-sided 95 % Wilson interval and one-sided 95 % lower bound.
z = 1.959963984540054; ph = k/n; den = 1 + z^2/n;
c = (ph + z^2/(2*n))/den; h = z*sqrt(ph*(1-ph)/n + z^2/(4*n^2))/den; lo = c - h; hi = c + h;
z1 = 1.644853626951472; den1 = 1 + z1^2/n;
lb1 = (ph + z1^2/(2*n))/den1 - z1*sqrt(ph*(1-ph)/n + z1^2/(4*n^2))/den1;
end
