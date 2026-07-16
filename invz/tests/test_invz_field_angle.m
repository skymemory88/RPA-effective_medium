function tests = test_invz_field_angle
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_transverse_baseline_anchor(testCase)
% Frozen x-field baseline (verified 2026-07-16, matches the Codex review probe):
% guards the seed/diagnostic edits against any accidental transverse change.
si = invz_single_ion(invz_ion(), 0.31, [4 0 0], struct('hyp', false));
verifyEqual(testCase, si.Jexp(1),        3.58424135285,  'AbsTol', 1e-9);
verifyEqual(testCase, si.Jexp(2),       -0.0689991117463,'AbsTol', 1e-10);
verifyEqual(testCase, si.E(2) - si.E(1), 0.369235620278, 'AbsTol', 1e-10);
verifyTrue(testCase, si.mf_converged);
verifyLessThan(testCase, si.mf_residual, 1e-12);
end

function test_sign_aware_branch_and_Fmf(testCase)
% Spec test 3 (second review finding 1, verified digit-for-digit): an explicit
% Bz selects the aligned branch; the VARIATIONAL F_mf ranks branches correctly
% (the naive shifted-spectrum comparison mis-ranks them and must not be used).
ion = invz_ion();  T = 0.31;
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
sp = invz_single_ion(ion, T, [2 0 +0.01], struct('hyp', false, 'order', true));
sm = invz_single_ion(ion, T, [2 0 -0.01], struct('hyp', false, 'order', true));
verifyGreaterThan(testCase, sp.Jexp(3), 0);
verifyLessThan(testCase, sm.Jexp(3), 0);
verifyLessThan(testCase, abs(sp.Jexp(3) + sm.Jexp(3)), 1e-10);      % exact Z2 mirror
% force the metastable anti-aligned branch with an explicit wrong-sign seed
sw = invz_single_ion(ion, T, [2 0 -0.01], struct('hyp', false, 'order', true, 'mz_seed', +1));
verifyGreaterThan(testCase, sw.Jexp(3), 0);                          % metastable branch reached
verifyLessThan(testCase, sm.F_mf, sw.F_mf);                          % aligned branch lower
verifyEqual(testCase, sm.F_mf, -21.4766457412, 'AbsTol', 1e-6);      % verified anchors
verifyEqual(testCase, sw.F_mf, -21.4696393612, 'AbsTol', 1e-6);
end

function test_mf_convergence_reporting(testCase)
ion = invz_ion();
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
si = invz_single_ion(ion, 0.31, [2 0 0.01], struct('hyp', false, 'order', true, 'mf_maxit', 1));
verifyFalse(testCase, si.mf_converged);
verifyEqual(testCase, si.mf_iters, 1);
verifyGreaterThan(testCase, si.mf_residual, 1e-12);
% hz_fixed mode: F_mf undefined by design (hz is imposed, not variational)
sf = invz_single_ion(ion, 0.31, [2 0 0], struct('hyp', false, 'hz_fixed', 5e-3));
verifyTrue(testCase, isnan(sf.F_mf));
verifyTrue(testCase, isfinite(sf.E0));
end
