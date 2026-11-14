function S = sample_inputs(n, seed, p)
%SAMPLE_INPUTS  Uncertain inputs (MATLAB/Octave Mersenne Twister; differs from NumPy's PCG64).
rng(seed, 'twister');
S.Uref = max(1, p.Uref + 0.10*p.Uref*randn(n,1));
S.MJ   = max(1, p.MJ + 0.10*p.MJ*randn(n,1));
S.z    = randn(n,1);
S.T0   = p.T0 + randn(n,1);
end
