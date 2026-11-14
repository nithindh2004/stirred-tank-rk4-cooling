function preq = required_pass_fraction(alpha, n)
%REQUIRED_PASS_FRACTION  Smallest sample pass fraction whose one-sided lower bound >= alpha.
lo = alpha; hi = 1 - 1e-12;
for k = 1:100
    mid = 0.5*(lo+hi); [~, ~, lb] = wilson(mid*n, n);
    if lb >= alpha, hi = mid; else, lo = mid; end
end
preq = hi;
end
