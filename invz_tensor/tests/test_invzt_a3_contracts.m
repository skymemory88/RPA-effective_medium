function tests = test_invzt_a3_contracts
% Focused A3/A3d map, co-state, status, and rank-provenance contracts.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
tc.TestData.ion = invz_ion();
g = invzt_qgrid(4, 'halfopen');
tc.TestData.lat = invzt_jq_tensor(tc.TestData.ion, g, ...
    struct('dpRng', 10, 'cache', false));
tc.TestData.T = 0.1;
tc.TestData.Ecut = 0.05;
[wn, ~, ~] = invz_matsubara(tc.TestData.T, tc.TestData.Ecut);
tc.TestData.seed = complex(zeros(3, 3, numel(wn)));
end

function test_pm_one_step_is_explicit_and_co_state_consistent(tc)
for dress = {'dominant', 'full'}
    clear invzt_sigma_tensor
    opts = pm_opts(tc, dress{1});
    opts.evaluation_kind = 'one_step_reevaluation';
    pt = invzt_solve_point(tc.TestData.ion, tc.TestData.T, [0.5 0 0], ...
        tc.TestData.lat, opts);

    verifyEqual(tc, pt.evaluation_kind, 'one_step_reevaluation');
    verifyFalse(tc, pt.converged);
    verifyEqual(tc, pt.Vmat, tc.TestData.seed, 'AbsTol', 0);
    verifyEqual(tc, pt.a3_residual, pt.Vnext - pt.Vmat, 'AbsTol', 0);
    verifyEqual(tc, pt.outer_residual, max(abs(pt.a3_residual(:))), 'AbsTol', 0);
    verifyEqual(tc, pt.a3_residual_component_max, ...
        max(abs(pt.a3_residual), [], 3), 'AbsTol', 0);
    verifyEqual(tc, pt.a3_residual_frequency_max, ...
        reshape(max(max(abs(pt.a3_residual), [], 1), [], 2), [], 1), 'AbsTol', 0);

    clear invzt_sigma_tensor
    repeat = invzt_solve_point(tc.TestData.ion, tc.TestData.T, [0.5 0 0], ...
        tc.TestData.lat, opts);
    verifyEqual(tc, repeat.chi_til, pt.chi_til, 'AbsTol', 0);
    verifyEqual(tc, repeat.Kmat, pt.Kmat, 'AbsTol', 0);
    verifyEqual(tc, repeat.Vnext, pt.Vnext, 'AbsTol', 0);
end
end

function test_pm_capped_and_converged_statuses_are_distinct(tc)
clear invzt_sigma_tensor
opts = pm_opts(tc, 'dominant');
capped = invzt_solve_point(tc.TestData.ion, tc.TestData.T, [0.5 0 0], ...
    tc.TestData.lat, opts);
verifyEqual(tc, capped.evaluation_kind, 'rejected_capped_solve');
verifyFalse(tc, capped.converged);
verifyEqual(tc, capped.Vmat, tc.TestData.seed, 'AbsTol', 0);

opts.tol_outer = 1e9; % Contract fixture: the same one-step map passes its declared tolerance.
accepted = invzt_solve_point(tc.TestData.ion, tc.TestData.T, [0.5 0 0], ...
    tc.TestData.lat, opts);
verifyEqual(tc, accepted.evaluation_kind, 'converged_solve');
verifyTrue(tc, accepted.converged);
verifyTrue(tc, accepted.map_converged);
verifyEqual(tc, accepted.a3_residual, accepted.Vnext - accepted.Vmat, 'AbsTol', 0);
verifyLessThan(tc, accepted.outer_residual, opts.tol_outer);
end

function test_ordered_rank_and_reevaluation_provenance(tc)
clear invzt_sigma_tensor
opts = ordered_opts(tc, 4, 4);
pt = invzt_solve_point_ordered(tc.TestData.ion, tc.TestData.T, [3 0 0], ...
    tc.TestData.lat, opts);
verifyEqual(tc, pt.evaluation_kind, 'one_step_reevaluation');
verifyFalse(tc, pt.converged);
verifyEqual(tc, pt.Vmat, tc.TestData.seed, 'AbsTol', 0);
verifyEqual(tc, pt.vb.n_vertex, 4);
verifyEqual(tc, pt.a3d.n_vertex, 4);
verifyEqual(tc, pt.a3d.vertex.requested_rank, 4);
verifyEqual(tc, pt.a3d.vertex.realized_rank, 4);
verifyEqual(tc, pt.a3d.vertex.max_vertex_states, 4);
verifyEqual(tc, pt.a3d.residual, pt.a3d.Vnext - pt.a3d.Vmat, 'AbsTol', 0);
verifyEqual(tc, pt.outer_residual, pt.a3d.residual_inf, 'AbsTol', 0);

bad = ordered_opts(tc, 5, 4);
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, tc.TestData.T, ...
    [3 0 0], tc.TestData.lat, bad), 'invzt:orderedA3Budget');
end

function test_reevaluation_and_seed_inputs_fail_closed(tc)
bad = pm_opts(tc, 'dominant');
bad.Vmat_seed = NaN;
verifyError(tc, @() invzt_solve_point(tc.TestData.ion, tc.TestData.T, ...
    [0.5 0 0], tc.TestData.lat, bad), 'invzt:VmatSeed');

oo = struct('mode', 'a3d', 'Ecut', tc.TestData.Ecut, ...
    'reeval', struct('mode', 'a3d'), 'n_vertex', 4, 'max_vertex_states', 4);
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, tc.TestData.T, ...
    [3 0 0], tc.TestData.lat, oo), 'invzt:a3dReevalSeed');
end

function test_a3_remains_outside_zero_dispatch_and_realaxis_scope(tc)
opts = pm_opts(tc, 'dominant');
verifyError(tc, @() invzt_solve_point(tc.TestData.ion, tc.TestData.T, ...
    [0 0 0], tc.TestData.lat, opts), 'invzt:a1ZeroField');
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, tc.TestData.T, ...
    [0.5 0 0], tc.TestData.lat, struct('mode', 'a3')), 'invzt:autoMode');

opts.evaluation_kind = 'one_step_reevaluation';
pt = invzt_solve_point(tc.TestData.ion, tc.TestData.T, [0.5 0 0], ...
    tc.TestData.lat, opts);
verifyError(tc, @() invzt_chi_realaxis(tc.TestData.ion, tc.TestData.T, ...
    [0.5 0 0], pt, 0.01, struct()), 'invzt:realaxisMode');
end

function opts = pm_opts(tc, dress)
opts = struct('mode', 'a3', 'nlevels', 'three', 'Ecut', tc.TestData.Ecut, ...
    'dress', dress, 'max_outer', 1, 'Vmat_seed', tc.TestData.seed);
end

function opts = ordered_opts(tc, n_vertex, max_vertex_states)
opts = struct('mode', 'a3d', 'Ecut', tc.TestData.Ecut, 'max_outer', 1, ...
    'Vmat_seed', tc.TestData.seed, 'evaluation_kind', 'one_step_reevaluation', ...
    'n_vertex', n_vertex, 'max_vertex_states', max_vertex_states);
end
