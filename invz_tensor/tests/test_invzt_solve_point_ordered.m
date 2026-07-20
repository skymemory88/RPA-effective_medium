function tests = test_invzt_solve_point_ordered
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
tc.TestData.ion = ion;
tc.TestData.lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
end

function test_deep_ordered_point(tc)
% Deep FM point: spontaneous moment, converged Sigma, stable (crit >= 0 within tol).
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], tc.TestData.lat, struct());
verifyTrue(tc, pt.is_ordered);
verifyTrue(tc, pt.converged);
verifyEqual(tc, pt.moment_branch, 'spontaneous');
verifyGreaterThan(tc, abs(pt.m0), 1.0);              % LiHoF4 FM moment is O(5) here
verifyTrue(tc, isfinite(pt.Sigma0) && all(isfinite(pt.Sigma)));
verifyTrue(tc, isfinite(pt.alpha_m));
verifyGreaterThan(tc, pt.crit, -1e-3);               % broken-symmetry state locally stable
verifyEqual(tc, pt.J0z_used, tc.TestData.lat.info.Jcc0);
verifyEqual(tc, pt.Jxx0_used, tc.TestData.lat.info.Jaa0);
verifyLessThan(tc, pt.sumrule_rel, 0.1);
fprintf('ordered a1 @ 3T: m0=%.4f Sigma0=%.5f crit=%.5f alpha_m=%.4g sumrule=%.3g\n', ...
    pt.m0, pt.Sigma0, pt.crit, pt.alpha_m, pt.sumrule_rel);
end

function test_converged_state_is_self_consistent(tc)
% Review P2-2: on a converged exit, the returned (Sigma, K, lambda, alpha_m) must
% describe the SAME state -- re-evaluating one medium+self-energy pass at pt.Sigma
% must reproduce pt.Sigma to the outer tolerance (check-before-mix loop ordering).
ion = tc.TestData.ion;  lat = tc.TestData.lat;
pt = invzt_solve_point_ordered(ion, 0.1, [3.0 0 0], lat, struct());
verifyTrue(tc, pt.converged);
sg2 = local_one_pass(ion, 0.1, [3.0 0 0], lat, pt);   % helper below: one g(Sigma) map step
verifyLessThan(tc, max(abs(sg2.Sigma - pt.Sigma)), 1e-6);   % >= tol_outer scale, not machine
end

function sg = local_one_pass(ion, T, B, lat, pt)
% One medium + moment-form self-energy evaluation AT pt.Sigma (no mixing): the
% fixed-point map g(Sigma) evaluated at the returned state. WHOLE-CC (2026-07-20
% amendment): no dominant/rest split -- mirrors the solver's own medium step.
[wn, wts, beta] = invz_matsubara(T, 40);
c0 = invz_chi0z(pt.si, T, 1i*wn, struct('elastic', true));
c0_cc = real(squeeze(c0(3,3,:)));
g = real(invz_g(pt.tl, 1i*wn));
ctil = c0 ./ reshape(1 + pt.Sigma, 1, 1, numel(wn));
Gcc = invzt_gcc_lattice(ctil, lat);
Gloc = -Gcc(:);
G0til = -(c0_cc ./ (1 + pt.Sigma));
K = 1 ./ Gloc - 1 ./ G0til;
lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg = invz_sigma_ordered(pt.tl, lam, K, g, beta);
end

function test_pm_early_return(tc)
% WELL above the bare-MF boundary (~5.0 T -- NOT the QPT at 4.65-4.70; measured
% m0 = 1.17 at 4.8 T, so 4.8 would NOT early-return) the order-mode MF relaxes
% to ~0: paramagnetic early return.
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [5.5 0 0], tc.TestData.lat, struct());
verifyFalse(tc, pt.is_ordered);
verifyFalse(tc, pt.converged);
verifyEqual(tc, pt.moment_branch, 'none');
verifyTrue(tc, isnan(pt.Sigma0) && isnan(pt.crit));
verifyTrue(tc, isempty(pt.tl));
verifyLessThan(tc, abs(pt.m0), 1e-2);
end

function test_mf_moment_persists_past_QPT(tc)
% Locks the P0-1 finding into the suite: the bare-MF moment is STILL nonzero at
% the 4.8 T PM anchor (measured m0 = 1.1717) -- which is exactly why phase
% selection must NOT be ordered-first (see invzt_solve_auto, Task 5).
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [4.8 0 0], tc.TestData.lat, struct());
verifyTrue(tc, abs(pt.m0) > 1.0);
fprintf('P0-1 lock: bare-MF m0(4.8 T) = %.4f (QPT is at 4.65-4.70 T)\n', pt.m0);
end

function test_longitudinal_rejected(tc)
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0.1], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
end

function test_full_dress_a3_rejected(tc)
% Full-dress 'a3' is PERMANENTLY rejected (136-state vertex budget-refused).
% This assertion stays valid after Task 7D adds 'a3d' to the mode gate.
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('mode', 'a3')), 'invzt:orderedMode');
end

function test_nlevels_std_only(tc)
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('nlevels', 'three')), 'invzt:orderedNlevels');
end

function test_split_knobs_rejected(tc)
% 2026-07-20 amendment: the ordered medium is WHOLE-CC (no dominant/rest split);
% the PM solver's 'Esplit'/'chi_rest' knobs have no meaning here and must fail loud.
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('Esplit', 0.4)), 'invzt:orderedSplitKnobs');
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('chi_rest', false)), 'invzt:orderedSplitKnobs');
end
