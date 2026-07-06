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
