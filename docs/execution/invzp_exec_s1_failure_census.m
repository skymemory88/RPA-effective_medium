function R = invzp_exec_s1_failure_census(Bx, save_path, solve_over)
%INVZP_EXEC_S1_FAILURE_CENSUS Per-node failure census for one ordered field column.
%
% Execution packet S1 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY: runs the production ordered node loop (invz_ordered_trace ->
% invz_hmf_ordered, ordered_mode = jensen, static_medium = resummed) at ONE
% transverse field and classifies every h-node by its proximate failure cause.
% It changes no solver behaviour and gates nothing.
%
%   Bx          transverse field in T (default 3.825 -- the diagnosed sliver,
%               invzp_convg_diagnosis.md SS5)
%   save_path   optional .mat path for the raw trace + census
%   solve_over  optional struct merged into the solve opts (e.g. mix_outer,
%               max_outer, static_medium, emt_static)
%
% Returns R with fields:
%   .meta     T, Bx, grid, dpRng, solve opts actually used, hmf_status, hstar
%   .tab      table, one row per traced node, with the acceptance verdict, the
%             term_reason / medium_status strings, the A-D block residuals and
%             the pole-proximity observables (Dq_min, Dq_abs_min, D_uni,
%             gstat_local_denom, Gstat, K0, r)
%   .census   struct: counts by (accepted, term_reason, medium_status)
%
% Fixed conditions mirror the production driver invz_run_spectra.m: T = 0.1 K,
% 16^3 coupling grid, dpRng 30, transverse_mf legacy_x, hyp true.
if nargin < 1 || isempty(Bx), Bx = 3.825; end
if nargin < 2, save_path = ''; end
if nargin < 3 || isempty(solve_over), solve_over = struct(); end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

T  = 0.1;
ion = invz_ion();
bz  = struct('grid', [16 16 16], 'dpRng', 30, 'cache', true, 'dipole', 'bruteforce');
[J, ci, Jaa0] = invz_bz_couplings(ion, bz);
J = J(:);

so = struct('hyp', true, 'J0eff', ci.Jcc0, 'Jxx0', Jaa0, ...
    'transverse_mf', 'legacy_x', 'static_medium', 'resummed', ...
    'mix_outer', 0.40, 'max_outer', 1000, 'cold_acceleration', 'none');
fn = fieldnames(solve_over);
for k = 1:numel(fn), so.(fn{k}) = solve_over.(fn{k}); end

t0 = tic;
trace = invz_ordered_trace(ion, T, [Bx 0 0], J, ...
    struct('solve', so, 'grid', bz.grid, 'dpRng', bz.dpRng, 'lattice_info', ci));
wall = toc(t0);

nodes = trace.nodes;
n = numel(nodes);
get = @(k, f, d) local_get(nodes(k), f, d);

id        = zeros(n,1);   h        = zeros(n,1);   accepted = false(n,1);
resA      = nan(n,1);     resB     = nan(n,1);     resC     = nan(n,1);
resD      = nan(n,1);     resStat  = nan(n,1);     iters    = nan(n,1);
Dq_min    = nan(n,1);     Dq_absmin= nan(n,1);     D_uni    = nan(n,1);
gdenom    = nan(n,1);     Gstat    = nan(n,1);     K0       = nan(n,1);
rr        = nan(n,1);     mm       = nan(n,1);     crit     = nan(n,1);
predictor = false(n,1);   hitmax   = false(n,1);
term      = cell(n,1);    medst    = cell(n,1);    seedk    = cell(n,1);

for k = 1:n
    id(k)        = get(k, 'id', k);
    h(k)         = get(k, 'h', NaN);
    accepted(k)  = logical(get(k, 'accepted', false));
    resA(k)      = get(k, 'resid_A', NaN);
    resB(k)      = get(k, 'resid_B', NaN);
    resC(k)      = get(k, 'resid_C', NaN);
    resD(k)      = get(k, 'resid_D', NaN);
    resStat(k)   = get(k, 'resid_static', NaN);
    iters(k)     = get(k, 'outer_iters', NaN);
    hitmax(k)    = logical(get(k, 'outer_hit_max', false));
    Dq_min(k)    = get(k, 'Dq_min', NaN);
    Dq_absmin(k) = get(k, 'Dq_abs_min', NaN);
    D_uni(k)     = get(k, 'D_uni', NaN);
    gdenom(k)    = get(k, 'gstat_local_denom', NaN);
    Gstat(k)     = get(k, 'Gstat', NaN);
    K0(k)        = get(k, 'K0', NaN);
    rr(k)        = get(k, 'r', NaN);
    mm(k)        = get(k, 'm', NaN);
    crit(k)      = get(k, 'crit', NaN);
    predictor(k) = logical(get(k, 'is_predictor', false));
    term{k}      = char(string(get(k, 'term_reason', '')));
    medst{k}     = char(string(get(k, 'medium_status', '')));
    seedk{k}     = char(string(get(k, 'seed_kind', '')));
end

tab = table(id, h, accepted, predictor, term, medst, seedk, iters, hitmax, ...
    resA, resB, resC, resD, resStat, Dq_min, Dq_absmin, D_uni, gdenom, ...
    Gstat, K0, rr, mm, crit);

keys = strcat(string(accepted), '|', string(term), '|', string(medst));
[uk, ~, ic] = unique(keys);
cnt = accumarray(ic, 1);
census = table(uk, cnt, 'VariableNames', {'accepted_term_medium', 'count'});

R = struct('meta', struct('T', T, 'Bx', Bx, 'grid', bz.grid, 'dpRng', bz.dpRng, ...
        'solve', so, 'hmf_status', trace.result.hmf_status, ...
        'hstar', trace.result.hstar, 'wall_s', wall, 'n_nodes', n, ...
        'n_accepted', nnz(accepted), 'n_failed', nnz(~accepted)), ...
    'tab', tab, 'census', census);

fprintf('=== S1 census: T = %.3f K, Bx = %.4f T ===\n', T, Bx);
fprintf('hmf_status = %s   hstar = %.10g   wall = %.1f s\n', ...
    trace.result.hmf_status, trace.result.hstar, wall);
fprintf('nodes = %d   accepted = %d   failed = %d\n', n, nnz(accepted), nnz(~accepted));
disp(census);
fprintf('--- failed nodes ---\n');
disp(tab(~accepted, :));
fprintf('--- accepted nodes ---\n');
disp(tab(accepted, :));

if ~isempty(save_path)
    save(save_path, 'trace', 'R', '-v7.3');
    fprintf('saved: %s\n', save_path);
end
end

function v = local_get(s, f, d)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
    if ~ischar(v) && ~isstring(v), v = v(1); end
else
    v = d;
end
end
