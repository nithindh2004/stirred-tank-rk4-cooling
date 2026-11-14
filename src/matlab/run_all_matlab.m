function run_all_matlab(outDir)
%RUN_ALL_MATLAB  Independent MATLAB/Octave implementation of the V3 workflow.
%
%   run_all_matlab                 % writes to ../../data/verification/octave
%   run_all_matlab('my_out_dir')
%
% Runs unchanged in MATLAB (R2019b or newer) and GNU Octave 8.x; no toolboxes.
% It reproduces the deterministic results of the Python reference
% implementation exactly and the Monte Carlo results to within sampling error:
% MATLAB/Octave's Mersenne Twister generator produces different samples from
% NumPy's PCG64 even with the same seed numbers.
%
% Coolant flow is commanded; the delivered flow is the saturated actuator
% min(mMax, mcmd*(1 + 0.05 z)). Pump energy uses the delivered flow.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'supporting_functions'));
if nargin < 1, outDir = fullfile(here, '..', '..', 'data', 'verification', 'octave'); end
if ~exist(outDir, 'dir'), mkdir(outDir); end
p = default_params(); alpha = 0.95; t0 = tic;
if exist('OCTAVE_VERSION', 'builtin'), fprintf('Runtime: GNU Octave %s\n', OCTAVE_VERSION);
else, fprintf('Runtime: MATLAB %s\n', version); end

%% 1. Deterministic model and verification ---------------------------------
[tB, YB] = rk4_sim(300, 0.20, p, 1.0);
[eL, eJ] = propagate_exact(p.T0, p.Tin, U_coef(300, 0.20, p, p.Uref)*p.A, 0.20, p.MJ, tB, p);
Told = p.Tin + (p.T0 - p.Tin)*exp(-p.Uref*p.A/(p.ML*p.CpL)*p.tf);
fprintf('Original fixed-coolant T_L(1800) = %.6f K\n', Told);
fprintf('Coupled baseline T_L(1800) = %.6f K, T_J = %.6f K\n', YB(end,1), YB(end,2));
fprintf('Max |RK4(dt=1) - exact| = %.3e K\n', max(max(abs([YB(:,1)-eL, YB(:,2)-eJ]))));
dts = [60 30 20 10 5 2 1 0.5 0.1]; conv = zeros(numel(dts), 3);
for q = 1:numel(dts)
    [t, Y] = rk4_sim(300, 0.20, p, dts(q));
    [aL, aJ] = propagate_exact(p.T0, p.Tin, U_coef(300, 0.20, p, p.Uref)*p.A, 0.20, p.MJ, t, p);
    conv(q,:) = [dts(q), Y(end,1), max(max(abs([Y(:,1)-aL, Y(:,2)-aJ])))];
end
write_csv(fullfile(outDir, 'timestep_convergence.csv'), {'dt_s','T_liquid_tf_K','max_traj_error_K'}, conv);
sel = dts >= 2 & dts <= 30; pf = polyfit(log(dts(sel)), log(conv(sel,3)'), 1);
fprintf('Observed RK4 order (dt 2-30 s) = %.4f\n', pf(1));
Qliq = p.ML*p.CpL*(p.T0 - YB(end,1));
Qout = simpson_rule(0.20*p.CpC*(YB(:,2) - p.Tin), tB(2) - tB(1));
Qst = p.MJ*p.CpC*(YB(end,2) - p.Tin);
fprintf('Energy residual (Simpson) = %.6f J\n', Qliq - Qout - Qst);
UA = U_coef(300, 0.2, p, p.Uref)*p.A; C = 0.2*p.CpC; G = UA*C/(UA + C);
fprintf('Reduced-order UA_eff = %.4f W/K -> T_L(1800) = %.4f K\n', G, p.Tin + (p.T0 - p.Tin)*exp(-G*p.tf/(p.ML*p.CpL)));

%% 2. Deterministic optimum --------------------------------------------------
[rN, mN, EN] = nominal_optimum(p);
fprintf('Nominal optimum: %.4f rpm, %.4f kg/s, E = %.6f kWh\n', rN, mN, EN);

%% 3. Chance-constrained fixed setting: optimise on OPT, validate on VAL -----
OPT = sample_inputs(20000, 20261109, p);
VAL = sample_inputs(200000, 20261120, p);       % final-release validation seed (see CHANGELOG)
[rb, mb, Eb] = chance_constrained(alpha, OPT, p);
rec = ceil(rb/5 - 1e-9)*5;                    % pre-declared rounding rule (OPT data only)
x = mc_final_TL(rec, mb, VAL, p); k = sum(x <= p.Ttarget); [lo, hi, lb] = wilson(k, numel(x));
Ev = mean(mc_energy_kWh(rec, mb, VAL, p));
fprintf('Chance-constrained optimum (OPT): %.3f rpm @ %.3f kg/s, E = %.4f kWh\n', rb, mb, Eb);
fprintf('Recommended setting: %d rpm @ %.2f kg/s; VAL P(pass) = %.4f [%.4f, %.4f], one-sided LB = %.4f, E = %.4f kWh\n', ...
    rec, mb, k/numel(x), lo, hi, lb, Ev);
legacy = [125 0.49; rN mN; 300 0.50; 285 0.50; rec mb]; leg = zeros(size(legacy,1), 6);
for q = 1:size(legacy,1)
    xq = mc_final_TL(legacy(q,1), legacy(q,2), VAL, p); kq = sum(xq <= p.Ttarget); [l1, h1, b1] = wilson(kq, numel(xq));
    leg(q,:) = [legacy(q,:), kq/numel(xq), l1, h1, mean(mc_energy_kWh(legacy(q,1), legacy(q,2), VAL, p))];
end
write_csv(fullfile(outDir, 'designs_VAL.csv'), {'rpm','m_cmd','pass_frac','ci_lo','ci_hi','mean_energy_kWh'}, leg);
fprintf('V2 recommendation 285 rpm under V3 model: P(pass) = %.4f\n', leg(4,3));

%% 4. Adaptive policy --------------------------------------------------------
% Configuration (N1 = 200 rpm, t_obs = 600 s) is the one selected by the Python
% CAL scan (data/adaptive_policy/configuration_scan_CAL.csv). The margin is
% calibrated here independently on CAL, rounded up to the next 0.05 K by the
% pre-declared rule, and validated once on VAL. VAL2 is used only if the policy
% fails the acceptance rule and is modified by the pre-declared rule (+0.05 K).
rpm1 = 200; tobs = 600; sigma = 0.1;
CAL = sample_inputs(20000, 20261112, p);
preq = required_pass_fraction(alpha, 20000); dlo = 0; dhi = 3;
for it = 1:18
    mid = 0.5*(dlo + dhi);
    if adaptive_policy(CAL, rpm1, 0.5, tobs, mid, p, 31, sigma) >= preq, dhi = mid; else, dlo = mid; end
end
dset = ceil(dhi/0.05 - 1e-9)*0.05;            % pre-declared practical rounding (CAL data only)
[pa, Ea, cap, uaErr] = adaptive_policy(VAL, rpm1, 0.5, tobs, dset, p, 40, sigma);
[al, ah, ab] = wilson(pa*numel(Ea), numel(Ea)); stage = 1;
if ab < alpha                                 % pre-declared modification, fresh data
    dset = dset + 0.05; VAL2 = sample_inputs(200000, 20261121, p);
    [pa, Ea, cap, uaErr] = adaptive_policy(VAL2, rpm1, 0.5, tobs, dset, p, 41, sigma);
    [al, ah, ab] = wilson(pa*numel(Ea), numel(Ea)); stage = 2;
end
fprintf('Adaptive policy: calibrated delta %.4f K, rounded to %.2f K; VAL P(pass) = %.4f [%.4f, %.4f], LB %.4f, mean E %.4f kWh, capped %.4f, UA err sd %.4f (stage %d)\n', ...
    dhi, dset, pa, al, ah, ab, mean(Ea), mean(cap), std(uaErr), stage);
write_csv(fullfile(outDir, 'adaptive_validation.csv'), {'stage','delta_calibrated_K','delta_K','pass_frac','ci_lo','ci_hi','lb95','mean_energy_kWh','frac_capped'}, ...
    [stage, dhi, dset, pa, al, ah, ab, mean(Ea), mean(cap)]);
fprintf('Completed in %.1f s.\n', toc(t0));
end
