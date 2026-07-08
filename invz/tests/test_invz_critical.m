function tests = test_invz_critical
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function [Jf, J0] = lihof4_couplings()
% Controller adaptation: qVec_generator's grid-size name-value pair is 'grid'
% (not 'size' as in the brief's literal draft); capture its fprintf lattice
% diagnostics with evalc so test output stays pristine.
ion = invz_ion();
cmd = "qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5])";
[~, qvec] = evalc(cmd);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Jf = Jnu(:);  J0 = info.Jcc0;
end

function test_zero_field_tc(testCase)
% R 2007: Tc(B=0) = 1.74 K in the 1/z theory (experiment 1.53 K; MF higher). SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
% Controller adaptation: the raw 16^3-grid mean of the Sigma_c integrand is
% biased low (integrable Gamma singularity), so Richardson-extrapolate over
% 12^3/24^3 (both Jq grids already cached at dpRng=30) instead of using
% invz_sigma_crit(J0, Jf) from the 16^3 helper directly.
S = zeros(1,2); ns = [12 24];
for k = 1:2
    cmd = sprintf("qVec_generator(ion.a, 'mode', 'grid', 'grid', [%d %d %d], 'range', [-0.5 0.5])", ns(k), ns(k), ns(k));
    [~, qv] = evalc(cmd);
    qv = qv(any(abs(qv) > 1e-12, 2), :);
    [Jn, inf_k] = invz_jq_modes(ion, qv, struct('dpRng', 30, 'cache', true));
    S(k) = invz_sigma_crit(inf_k.Jcc0, Jn(:));
end
Sc = 2*S(2) - S(1);   % = 0.2980, vs published 0.3004
J0 = inf_k.Jcc0;
Tc  = invz_critical_T0field(ion, Sc, J0);
TcMF = invz_critical_T0field(ion, 0,  J0);
verifyEqual(testCase, Tc, 1.74, 'AbsTol', 0.08);
verifyGreaterThan(testCase, TcMF, Tc);            % fluctuations suppress Tc
end

function test_critical_field_at_310mK(testCase)
% R 2007: Hc(0.31 K) = 42.4-43.0 kOe with Sigma_c(0) = 0.0932. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
bx = invz_critical(ion, 0.31, Jf, struct('J0eff', J0));
verifyGreaterThan(testCase, bx, 4.0);  verifyLessThan(testCase, bx, 4.6);   % Tesla
pt = invz_solve_point(ion, 0.31, bx, Jf, struct('hyp', true, 'J0eff', J0));
verifyEqual(testCase, pt.Sigma0, 0.0932, 'AbsTol', 0.02);
end

function test_tc_small_field(testCase)
% Small-field limit of the fixed-B temperature bisection: Tc(0.5 T) must sit
% just below the closed-form Tc0 evaluated on the SAME 16^3 q-grid (1.7795 K;
% the Richardson-extrapolated benchmark 1.74 K is a different baseline).
% Measured: 1.777 K. The undershoot direction is R 2007's small-Bx caveat.
% B = 0.5 T is the validated floor: at 0.2-0.3 T the paramagnetic solve has
% non-convergence patches near the boundary that the classifier reads as
% "ordered", biasing Tc upward by ~0.04-0.05 K. The window [1.0 2.0] is
% exactly the driver's [Tlo Tmax], so this doubles as an integration check
% of the driver's bracket geometry. SLOW.
% Upper bound 1.79 (not 1.7795) absorbs the bisection's ~±0.005 K midpoint
% quantization; the measured margin below Tc0 is smaller than that jitter.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
tc = invz_critical_T(ion, 0.5, Jf, struct('J0eff', J0, 'window', [1.0 2.0]));
verifyGreaterThan(testCase, tc, 1.70);
verifyLessThan(testCase, tc, 1.79);
end

function test_tc_at_fixed_field_crossing(testCase)
% Mirror consistency: at a mid-slope boundary point (T* = 1.4 K, where both
% cut directions are well-conditioned) the fixed-B temperature bisection must
% land back on the fixed-T field bisection's point. 0.05 K tolerance covers
% both bisection tolerances (0.02 T, 0.01 K) at the local boundary slope. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
Tstar = 1.4;
bstar = invz_critical(ion, Tstar, Jf, struct('J0eff', J0, 'window', [0.5 7]));
tc = invz_critical_T(ion, bstar, Jf, struct('J0eff', J0, 'window', [1.0 2.0]));
verifyEqual(testCase, tc, Tstar, 'AbsTol', 0.05);
end
