function tests = test_invzt_solver_contracts
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g = invzt_qgrid(4, 'halfopen');
tc.TestData.ion = ion;
tc.TestData.lat = invzt_jq_tensor(ion, g, struct('dpRng', 10, 'cache', false));
end

function test_pm_rejected_exit_is_transactional(tc)
% A capped solve returns the state that was actually evaluated: no final
% post-evaluation mix may make Sigma disagree with K/G/lambda.
pt = invzt_solve_point(tc.TestData.ion, 0.1, [3 0 0], tc.TestData.lat, ...
    struct('max_outer', 1));
verifyFalse(tc, pt.converged);
verifyEqual(tc, pt.outer_iters, 1);
verifyEqual(tc, pt.Sigma, zeros(size(pt.Sigma)), 'AbsTol', 0);
verifyTrue(tc, isfinite(pt.outer_residual) && pt.outer_residual > 0);
end

function test_ordered_rejected_exit_is_transactional(tc)
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3 0 0], ...
    tc.TestData.lat, struct('max_outer', 1));
verifyFalse(tc, pt.converged);
verifyEqual(tc, pt.outer_iters, 1);
verifyEqual(tc, pt.Sigma, zeros(size(pt.Sigma)), 'AbsTol', 0);
verifyTrue(tc, isfinite(pt.outer_residual) && pt.outer_residual > 0);
end

function test_pm_seed_validation(tc)
base = @() invzt_solve_point(tc.TestData.ion, 0.1, [3 0 0], ...
    tc.TestData.lat, struct('Sigma_seed', NaN));
verifyError(tc, base, 'invzt:SigmaSeed');
end

function test_auto_zero_field_fails_closed(tc)
[pt, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [0 0 0], ...
    tc.TestData.lat, struct());
verifyEmpty(tc, pt);
verifyEqual(tc, phase, 0);
verifyEqual(tc, di.para.err, 'invzt:a1ZeroField');
verifyEqual(tc, di.para.reason, 'invzt:a1ZeroField');
verifyFalse(tc, di.ordered.attempted);
verifyEqual(tc, di.ordered.reason, 'zero_field_unsupported');
end

function test_auto_rejects_phase_dependent_option_surface(tc)
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [3 0 0], ...
    tc.TestData.lat, struct('chi_rest', false)), 'invzt:autoSplitKnobs');
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [3 0 0], ...
    tc.TestData.lat, struct('nlevels', 'three')), 'invzt:autoNlevels');
end

function test_invalid_pm_sigma_provenance_does_not_seed_ordered(tc)
% Raising sigma_floor above every finite PM result forces the handoff-invalid
% path without fabricating a malformed point or accepting a rejected seed.
[~, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [5.5 0 0], ...
    tc.TestData.lat, struct('sigma_floor', 1e6));
verifyEqual(tc, phase, 0);
verifyTrue(tc, di.para.converged);
verifyFalse(tc, di.para.accepted);
verifyFalse(tc, di.ordered.attempted);
verifyEqual(tc, di.para.reason, 'sigma_below_floor');
verifyEqual(tc, di.ordered.reason, 'invalid_pm_handoff');
end

function test_auto_accepts_certified_modified_field_root(tc)
[pt, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [3 0 0], ...
    tc.TestData.lat, struct());
verifyEqual(tc, phase, 1);
verifyTrue(tc, di.ordered.accepted);
verifyEqual(tc, di.para.reason, 'unstable');
verifyEqual(tc, di.ordered.reason, 'accepted');
verifyTrue(tc, pt.converged);
verifyGreaterThan(tc, 1 + pt.Sigma0, 0);
verifyLessThan(tc, abs(pt.hmf_residual), 1e-10);
verifyLessThan(tc, pt.outer_residual, 1e-8);
end

function test_realaxis_rejects_bruteforce_lattice_mismatch(tc)
pt = invzt_solve_point(tc.TestData.ion, 0.1, [6 0 0], ...
    tc.TestData.lat, struct());
verifyTrue(tc, pt.converged);
verifyError(tc, @() invzt_chi_realaxis(tc.TestData.ion, 0.1, [6 0 0], ...
    pt, [0.005; 0.01], struct('qsel', [0.1 0 0], 'dpRng', 11, ...
    'cache', false, 'force_sigma0', true)), 'invzt:realaxisDipoleMismatch');
end
