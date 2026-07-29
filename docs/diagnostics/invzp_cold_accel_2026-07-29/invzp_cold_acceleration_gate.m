function R = invzp_cold_acceleration_gate(save_path)
%INVZP_COLD_ACCELERATION_GATE Reproduce the safeguarded cold-start discriminator.
%
% Frozen configuration:
%   T=0.10 K, legacy 16^3 brute-force/dpRng-30 couplings, resummed medium,
%   legacy_x transverse MF.  The 4.400 T comparisons are accelerated 0.70/200,
%   ordinary cold 0.70/1000, and QCP-down-equivalent warm 0.50/1000.  The
%   4.300 and 3.825 T classification traces use accelerated 0.70/1000.
%
% The returned/saved object is compact: it retains full stationary vectors for
% the three accepted 4.400 T states, node summaries, all proposal records, and
% only the last twelve iterations of each failed node. It deliberately does not
% save the full trace's repeated provenance payload.
if nargin < 1, save_path = ''; end

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

ion = invz_ion();
bz = struct('grid', [16 16 16], 'dpRng', 30, ...
            'cache', true, 'dipole', 'bruteforce');
[J, ci, Jaa0] = invz_bz_couplings(ion, bz);
J = J(:);
fresh_bz = bz;
fresh_bz.cache = false;
[Jfresh, cifresh, Jaa0fresh] = invz_bz_couplings(ion, fresh_bz);
Jfresh = Jfresh(:);
if ~isequaln(J, Jfresh) || ~isequaln(ci.Jcc0, cifresh.Jcc0) || ...
        ~isequaln(Jaa0, Jaa0fresh)
    error('invzp:coldAccelCoupling', ...
        'The cached and freshly generated frozen coupling fixtures differ.');
end
base = struct('hyp', true, 'J0eff', ci.Jcc0, 'Jxx0', Jaa0, ...
    'transverse_mf', 'legacy_x', 'static_medium', 'resummed', ...
    'ordered_mode', 'jensen');

fprintf('4.400 T accelerated 0.70/200 point\n');
oa = base;
oa.mix_outer = 0.70;
oa.max_outer = 200;
oa.cold_acceleration = 'signed_aitken1';
p_acc = invz_solve_point_ordered(ion, 0.10, 4.400, J, oa);

fprintf('4.400 T ordinary cold 0.70/1000 reference\n');
ol = base;
ol.mix_outer = 0.70;
ol.max_outer = 1000;
p_long = invz_solve_point_ordered(ion, 0.10, 4.400, J, ol);

fprintf('4.425 -> 4.400 T predictor-only warm 0.50/1000 reference\n');
ow = base;
ow.mix_outer = 0.50;
ow.max_outer = 1000;
p425 = invz_solve_point_ordered(ion, 0.10, 4.425, J, ow);
ow.hmf_seed = p425.hmf_prof.hmf_seed_out;
p_warm = invz_solve_point_ordered(ion, 0.10, 4.400, J, ow);

fprintf('4.400 T accelerated trace\n');
t440 = run_trace(ion, 4.400, J, ci.Jcc0, Jaa0, 200);
fprintf('4.300 T accelerated classification trace\n');
t430 = run_trace(ion, 4.300, J, ci.Jcc0, Jaa0, 1000);
fprintf('3.825 T accelerated classification trace\n');
t3825 = run_trace(ion, 3.825, J, ci.Jcc0, Jaa0, 1000);

R = struct();
R.schema = 'invzp_cold_acceleration_gate/v1';
R.created_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
R.config = struct('T', 0.10, 'grid', [16 16 16], 'dpRng', 30, ...
    'dipole', 'bruteforce', 'cache_used', true, ...
    'fresh_cache_bit_identity', true, 'static_medium', 'resummed', ...
    'transverse_mf', 'legacy_x', 'J0eff', ci.Jcc0, 'Jxx0', Jaa0, ...
    'coupling_digest', invz_exact_numeric_digest(J), ...
    'historical_expected_digest', ...
        'ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17', ...
    'historical_digest_match', strcmp(invz_exact_numeric_digest(J), ...
        'ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17'));
R.states_440 = struct('accelerated', compact_state(p_acc), ...
    'long_cold', compact_state(p_long), 'warm', compact_state(p_warm));
R.comparison_440 = struct( ...
    'accelerated_vs_long', state_distance(p_acc, p_long), ...
    'accelerated_vs_warm', state_distance(p_acc, p_warm));
R.trace_440 = compact_trace(t440);
R.trace_430 = compact_trace(t430);
R.trace_3825 = compact_trace(t3825);

if ~isempty(save_path)
    save(save_path, 'R', '-v7');
end
fprintf('done\n');
end

function tr = run_trace(ion, B, J, J0eff, Jxx0, max_outer)
so = struct('hyp', true, 'J0eff', J0eff, 'Jxx0', Jxx0, ...
    'transverse_mf', 'legacy_x', 'static_medium', 'resummed', ...
    'mix_outer', 0.70, 'max_outer', max_outer, ...
    'cold_acceleration', 'signed_aitken1');
tr = invz_ordered_trace(ion, 0.10, [B 0 0], J, struct('solve', so));
end

function s = compact_state(p)
s = struct('accepted', logical(p.converged), 'stable', logical(p.stable_1z), ...
    'hmf', p.hmf, 'm0', p.m0, 'Sigma0', p.Sigma0, 'K0', p.K(1), ...
    'final_resid', p.final_resid, 'D_uni', p.D_uni, 'Dq_min', p.Dq_min, ...
    'Sigma', p.Sigma, 'K', p.K, 'lambda', p.lambda, ...
    'predictor_acceleration', p.hmf_prof.predictor_acceleration);
end

function d = state_distance(a, b)
d = struct('Sigma_inf', max(abs(a.Sigma-b.Sigma)), ...
    'K_inf', max(abs(a.K-b.K)), ...
    'lambda_inf', max(abs(a.lambda-b.lambda)), ...
    'hmf_abs', abs(a.hmf-b.hmf), 'm0_abs', abs(a.m0-b.m0), ...
    'Sigma0_abs', abs(a.Sigma0-b.Sigma0), ...
    'K0_abs', abs(a.K(1)-b.K(1)));
end

function c = compact_trace(tr)
nodes = tr.nodes;
nb = numel(nodes);
blank = struct('id', NaN, 'phase', '', 'h', NaN, 'seed_kind', '', ...
    'accepted', false, 'outer_iters', NaN, 'term_reason', '', ...
    'resid_A', NaN, 'resid_B', NaN, 'resid_C', NaN, 'resid_D', NaN, ...
    'resid_static', NaN, 'Dq_abs_min', NaN, 'acceleration', struct());
ns = repmat(blank, 1, nb);
for k = 1:nb
    n = nodes(k);
    ns(k) = struct('id', n.id, 'phase', n.phase, 'h', n.h, ...
        'seed_kind', n.seed_kind, 'accepted', logical(n.accepted), ...
        'outer_iters', n.outer_iters, 'term_reason', n.term_reason, ...
        'resid_A', n.resid_A, 'resid_B', n.resid_B, ...
        'resid_C', n.resid_C, 'resid_D', n.resid_D, ...
        'resid_static', n.resid_static, 'Dq_abs_min', n.Dq_abs_min, ...
        'acceleration', n.acceleration);
end

iters = tr.iters;
keep = false(1, numel(iters));
if ~isempty(iters)
    keep = [iters.accel_attempted];
    failed_ids = [nodes(~[nodes.accepted]).id];
    for id = failed_ids
        ix = find([iters.node_id] == id);
        ix = ix(max(1, numel(ix)-11):end);
        keep(ix) = true;
    end
end
c = struct('result', tr.result, 'nodes', ns, ...
    'failed_count', nnz(~[nodes.accepted]), ...
    'proposal_count', nnz([iters.accel_attempted]), ...
    'selected_iters', iters(keep));
end
