function tests = test_invz_critical
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
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
% Small-field limit of the fixed-B temperature cut: Tc(0.5 T) must sit just below
% the closed-form Tc0 on the SAME 16^3 q-grid (1.7795 K), reflecting R 2007's
% small-Bx undershoot caveat. Uses an explicit window [1.0 2.0] to exercise the
% grid+interpolate path of invz_critical_T. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
tc = invz_critical_T(ion, 0.5, Jf, struct('J0eff', J0, 'window', [1.0 2.0]));
verifyGreaterThan(testCase, tc, 1.70);
verifyLessThan(testCase, tc, 1.79);
end

function test_tc_at_fixed_field_crossing(testCase)
% Mirror consistency: at a mid-slope boundary point (T* = 1.4 K, where both
% cut directions are well-conditioned) the fixed-B temperature cut must land
% back on the fixed-T field cut's point. 0.05 K tolerance covers the two
% solvers' crossing tolerances (0.02 T, 0.005 K) at the local boundary slope. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
Tstar = 1.4;
bstar = invz_critical(ion, Tstar, Jf, struct('J0eff', J0, 'window', [0.5 7]));
tc = invz_critical_T(ion, bstar, Jf, struct('J0eff', J0, 'window', [1.0 2.0]));
verifyEqual(testCase, tc, Tstar, 'AbsTol', 0.05);
end

function test_tc_boundary_is_smooth(testCase)
% Regression: the boundary classifier must count only CONVERGED points, or
% critical-slowing-down failures near the boundary get misread as "ordered"
% and Tc(B) becomes non-smooth. Asserted via RMS discrete second difference
% (SMOOTHNESS, not monotonicity, since a genuine re-entrant nose is allowed). SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
Bs = 0.4:0.2:1.8;                                   % low-field / high-T band
Tc = arrayfun(@(B) invz_critical_T(ion, B, Jf, struct('J0eff', J0)), Bs);
verifyTrue(testCase, all(isfinite(Tc)), 'every field must return a finite Tc');
d2 = Tc(1:end-2) - 2*Tc(2:end-1) + Tc(3:end);      % discrete curvature proxy
verifyLessThan(testCase, sqrt(mean(d2.^2)), 0.02);  % smooth (old code ~0.14)
end
