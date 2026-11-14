function I = simpson_rule(y, h)
%SIMPSON_RULE  Composite Simpson's rule on an even number of equal intervals.
n = numel(y) - 1;
if mod(n, 2) == 1, error('simpson_rule needs an even number of intervals'); end
I = h/3*(y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)));
end
