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

function test_scalar_vs_vector_boundaries(testCase)
% Spec test 2: scalar B and [B 0 0] are literally the same solve at every
% scalar-accepting boundary (struct-exact equality: identical code path).
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';       % synthetic branch fixture (fast, no lattice sum)
o   = struct('J0eff', 6.4e-3);
verifyEqual(testCase, invz_twolevel(ion, 0.31, [5.5 0 0]), invz_twolevel(ion, 0.31, 5.5));
verifyEqual(testCase, invz_twolevel_ordered(ion, 0.31, [2 0 0], 5e-3), ...
                      invz_twolevel_ordered(ion, 0.31, 2, 5e-3));
pt1 = invz_solve_point(ion, 0.31, 5.5, Jnu, o);
pt2 = invz_solve_point(ion, 0.31, [5.5 0 0], Jnu, o);
verifyEqual(testCase, pt2, pt1);
[pa1, ph1] = invz_solve_auto(ion, 0.31, 5.5, Jnu, o);
[pa2, ph2] = invz_solve_auto(ion, 0.31, [5.5; 0; 0], Jnu, o);   % column form too
verifyEqual(testCase, ph2, ph1);
verifyEqual(testCase, pa2, pa1);
w = (0.1:0.1:0.4).';
o1 = invz_chi_realaxis(ion, 0.31, 5.5, pt1, w, struct('eta', 1e-3));
o2 = invz_chi_realaxis(ion, 0.31, [5.5 0 0], pt1, w, struct('eta', 1e-3));
verifyEqual(testCase, o2, o1);
end

function test_longitudinal_routing_threshold(testCase)
% Spec test 7: just below bz_tol -> transverse path (strict-PM two-level);
% just above -> forced moment-form with machine-readable branch metadata.
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o   = struct('J0eff', 6.4e-3);
[ptA, phA] = invz_solve_auto(ion, 0.31, [5.5 0 0.5e-9], Jnu, o);   % |Bz| <= 1e-9: dead band
verifyEqual(testCase, phA, 2);                                      % strict paramagnet
verifyEqual(testCase, ptA.tl.m, 0, 'AbsTol', 1e-3);                 % invz_twolevel gate honored
[ptB, phB] = invz_solve_auto(ion, 0.31, [5.5 0 2e-9], Jnu, o);      % above: moment route
verifyEqual(testCase, phB, 1);
verifyEqual(testCase, ptB.moment_branch, 'field_induced');
verifyTrue(testCase, ptB.is_ordered);                               % moment-form self-energy flag
end

function test_early_return_struct_completeness(testCase)
% Spec test 12: every early return of invz_solve_point_ordered carries the full
% declared field set, so invz_solve_auto / the map never probe a missing member.
ion  = invz_ion();
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
flds = {'m0','is_ordered','converged','Sigma0','crit','si','tl','moment_branch'};
% (a) spontaneous-mode paramagnetic early return (PM point: no spontaneous moment)
pta = invz_solve_point_ordered(ion, 1.0, 5.5, Jnu, struct('J0eff', 6.4e-3));
verifyFalse(testCase, pta.is_ordered);
verifyEqual(testCase, pta.moment_branch, 'none');
cellfun(@(f) verifyTrue(testCase, isfield(pta, f), ['missing ' f]), flds);
% (b) forced_moment with a crippled MF loop -> mf-gate early return
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
ptb = invz_solve_point_ordered(ion, 0.31, [2 0 0.01], Jnu, ...
      struct('J0eff', 6.4e-3, 'forced_moment', true, 'mf_maxit', 1));
verifyFalse(testCase, ptb.converged);
verifyEqual(testCase, ptb.moment_branch, 'field_induced');
cellfun(@(f) verifyTrue(testCase, isfield(ptb, f), ['missing ' f]), flds);
end

function test_solve_auto_returns_failed_pto(testCase)
% Second-review finding 2: a failed longitudinal solve returns the pto (with si)
% so invz_spectra_map can compute its RPA-only overlay -- never pt = [].
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
[ptc, phc] = invz_solve_auto(ion, 0.31, [2 0 0.01], Jnu, ...
             struct('J0eff', 6.4e-3, 'mf_maxit', 1));
verifyEqual(testCase, phc, 0);
verifyFalse(testCase, isempty(ptc));
verifyFalse(testCase, isempty(ptc.si));
end
