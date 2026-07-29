function result = invzp_exec_s3_twin_polish(save_path)
%INVZP_EXEC_S3_TWIN_POLISH Newton-polish diagnostic for the 4.05 T "twin roots" question.
%
% Execution packet S3 of docs/execution/invzp_plan_execution_diary.md
% (blind_convg_plan.md S4 "Resulting merged program", first immediate-diagnostics bullet;
% blind WP1c). MEASUREMENT ONLY. It changes no solver behaviour and promotes nothing.
%
% Background fact under test (invzp_convg_diagnosis.md S6.2): "At one common h near 4.05 T,
% two states pass the numerical A-D contract, with r ~ 0.768 and r ~ 0.822, separated by a
% state-space distance of about 1.037e-3." This script asks: are these two genuinely distinct
% high-accuracy numerical zeros of the coupled residual, or is one a pseudo-root that only
% passes because the acceptance tolerance (tol_outer = 1e-8) is loose?
%
% *** WHAT THIS SCRIPT DOES NOT ESTABLISH ***
% Newton-polishing a seed down to a tiny backward (residual) error is NOT interval
% certification. It shows that a floating-point iterate exists whose residual, in the node's
% own coordinates, is smaller than the reported separation; it does NOT bound the true root
% inside a verified enclosure and does NOT rule out that both iterates lie in the basin of one
% ill-conditioned root. A rigorous "this is THE unique zero in this box" claim requires an
% interval-Newton / Krawczyk enclosure, which this script does not run. Every "polished" root
% below remains an uncertified high-accuracy floating-point approximation, not a certified zero.
%
% *** WHAT THIS SCRIPT DOES NOT DO, BY CONSTRAINT ***
% No tolerance is loosened anywhere below: invz_ordered_node_newton's A-D accept gate
% (invz_ordered_residual, called internally) always uses its own tol_outer, and polishing only
% ever TIGHTENS it (1e-8 -> 1e-10 -> 1e-12 -> 1e-14). No pole floor is added. No masked/bare
% state is substituted for an accepted one. Seeds are drawn from a fixed, pre-registered,
% reproducible grid (documented below) — none is hand-tuned to land near r = 0.768 or r = 0.822.
%
% METHOD
%  Step 1 -- coexistence scan. At every h in a fixed scan grid (h0 from 1e-5 to 3e-2, the range
%    given in the task -- see hgrid below), invzp_enumerate_ordered_roots is run with the SAME
%    27 explicit seeds: 9 K0 values spanning linspace(min(Jnu_flat), max(Jnu_flat), 9) crossed
%    with 3 Sigma perturbation scales {0, 0.3, 1.0} x ctx.Jscale x randn (one fixed rng stream,
%    seeded once before seed construction so the whole grid is reproducible). Every seed is
%    solved by invz_ordered_node_newton and clustered by invzp_cluster_ordered_full_states'
%    independent [Sigma;K0] distance metric (sigma_scale = k0_scale = ctx.Jscale,
%    same_tol = 1e-6, distinct_tol = 1e-4 in that rescaled metric -- chosen because same-root
%    numerical replicates measured diam ~1e-8..1e-11 and cross-root gaps measured >~1e-2 in a
%    preliminary calibration run at this fixture, i.e. the tolerance band sits four-plus orders
%    of magnitude from both boundaries it separates, not tuned to manufacture a split).
%    The number of distinct accepted clusters and their r are recorded at every h.
%  Step 2 -- polish. At the h with the richest multi-cluster coexistence, the TWO
%    largest-support clusters (most seeds converging there -- the least flukey candidates) are
%    each re-solved from their own medoid state by invz_ordered_node_newton with maxit raised
%    to 80 and tol_outer swept over {1e-8,1e-10,1e-12,1e-14}; the tightest tol_outer at which
%    info.accepted remains true is reported as that root's "polish".
%  Step 3 -- separation. 2-norm/inf-norm distance between the two polished u=[Sigma;K0], and
%    the plain |Delta r|, compared against the two polished residual_inf values and against
%    the fixed A-block tolerance 1e-8.
%
% save_path (optional): .mat path; default docs/execution/out/s3_twin_polish.mat

if nargin < 1 || isempty(save_path)
    save_path = '';
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(root,'invz_common'));
addpath(fullfile(root,'invz_projected'));
addpath(fullfile(root,'docs','diagnostics','invzp_solver_stability_2026-07-27'));
addpath(fullfile(root,'docs','diagnostics','biased_smooth_r_2026-07-28'));

% ------------------------------------------------------------------------------------------
% Fixture (matches the production driver; matches docs/execution/invzp_exec_s1_failure_census.m
% and invzp_exec_s2_defactor_gate.m exactly).
% ------------------------------------------------------------------------------------------
T = 0.1;  Bx = 4.05;
ion = invz_ion();
bz  = struct('grid',[16 16 16],'dpRng',30,'cache',true,'dipole','bruteforce');
[J, ci, Jaa0] = invz_bz_couplings(ion, bz);
J = J(:);
ctx = invz_ordered_node_context(ion, T, [Bx 0 0], J, ...
    struct('J0eff', ci.Jcc0, 'Jxx0', Jaa0, 'hyp', true, 'transverse_mf', 'legacy_x'));
nw = numel(ctx.wn);

% ------------------------------------------------------------------------------------------
% Step 1 -- fixed, reproducible seed spread + h scan.
% ------------------------------------------------------------------------------------------
rng(20260729, 'twister');   % ONE fixed stream for the whole seed grid; not re-seeded per h.
K0grid = linspace(min(J), max(J), 9).';
sigScales = [0 0.3 1.0];
seeds = repmat(struct('id','','state',struct('Sigma',[],'K0s',[])), 0, 0);
sid = 0;
for kk = 1:numel(K0grid)
    for sc = sigScales
        sid = sid + 1;
        Sigma = sc * ctx.Jscale * randn(nw,1);
        seeds(sid).id = sprintf('s%02d_K0_%+.4g_sc%.2g', sid, K0grid(kk), sc); %#ok<AGROW>
        seeds(sid).state = struct('Sigma', Sigma, 'K0s', K0grid(kk));
    end
end
n_seed = numel(seeds);

spec = struct('cluster', struct('sigma_scale', ctx.Jscale, 'k0_scale', ctx.Jscale, ...
        'same_tol', 1e-6, 'distinct_tol', 1e-4), ...
    'newton', struct('maxit', 30));

hgrid = unique(sort([ ...
    logspace(-5,-3,5), ...
    0.0010:0.0005:0.0090, ...
    logspace(-2,log10(0.03),5) ]));
nh = numel(hgrid);

scan = repmat(struct('h',NaN,'status','','status_detail','', ...
    'nclusters',0,'r',[],'member_count',[],'wall_s',NaN), nh, 1);
reports = cell(nh,1);

fprintf('=== S3 step 1: coexistence scan, %d h values x %d seeds ===\n', nh, n_seed);
for k = 1:nh
    h = hgrid(k);
    t0 = tic;
    rep = invzp_enumerate_ordered_roots(ctx, h, seeds, spec);
    wall = toc(t0);
    reports{k} = rep;
    scan(k).h = h;
    scan(k).status = rep.status;
    scan(k).status_detail = rep.status_detail;
    scan(k).wall_s = wall;
    if ~isempty(rep.roots)
        scan(k).nclusters = numel(rep.roots);
        scan(k).r = arrayfun(@(x) x.info.r, rep.roots);
        scan(k).member_count = arrayfun(@(x) numel(x.member_seed_ids), rep.roots);
    end
    fprintf('  h=%.6g  (%.1fs)  status=%s  nclusters=%d  r=[%s]\n', ...
        h, wall, rep.status, scan(k).nclusters, sprintf('%.6g ', sort(scan(k).r)));
end

coexist_idx = find([scan.nclusters] >= 2);
if isempty(coexist_idx)
    result = struct('schema','invzp_s3_twin_polish/v1', 'reproduced', false, ...
        'scan', scan, 'conclusion', 'not reproduced: no coexistence point was found in the scanned range');
    fprintf('\n=== RESULT ===\nNo h in the scanned range produced >= 2 distinct accepted clusters.\n');
    fprintf('CONCLUSION: not reproduced: no coexistence point was found in the scanned range\n');
    if ~isempty(save_path), save(save_path, 'result', 'reports'); end
    return
end

% Best coexistence h: a two-stage, pre-specified (not data-reactive) rule.
%   Stage 1 (preferred): among h with >= 2 clusters whose member_count >= 2 (i.e. at least two
%     INDEPENDENT seeds converge there each -- not a one-off numerical fluke), pick the h that
%     MAXIMISES the r-gap between its two largest-support such clusters. A clean, well-supported,
%     widely-separated split is the strongest available test of "genuinely distinct roots";
%     picking on raw cluster COUNT instead would reward fragmentation (many singleton clusters)
%     over robustness, which is not the question being asked.
%   Stage 2 (fallback, only if no h qualifies for stage 1): revert to maximising nclusters,
%     tie-broken by proximity of any two clusters' r to the S6.2-reported {0.768, 0.822} pair.
%     This never happens silently -- result.selection_rule records which stage fired.
robust_gap = nan(size(scan));
for kk = coexist_idx
    mc = scan(kk).member_count;
    big = find(mc >= 2);
    if numel(big) >= 2
        rv = scan(kk).r(big);
        [~, ord] = sort(mc(big), 'descend');
        top2 = ord(1:2);
        robust_gap(kk) = abs(rv(top2(1)) - rv(top2(2)));
    end
end
if any(isfinite(robust_gap))
    selection_rule = 'stage1_robust_pair_max_gap';
    [~, best_k] = max(robust_gap);
else
    selection_rule = 'stage2_fallback_max_nclusters';
    best_k = coexist_idx(1);
    best_score = local_pair_score(scan(best_k).r);
    for kk = coexist_idx(2:end)
        if scan(kk).nclusters > scan(best_k).nclusters || ...
                (scan(kk).nclusters == scan(best_k).nclusters && ...
                 local_pair_score(scan(kk).r) < best_score)
            best_k = kk;
            best_score = local_pair_score(scan(kk).r);
        end
    end
end
h_best = hgrid(best_k);
rep_best = reports{best_k};
fprintf('\n=== Best coexistence (%s): h = %.8g, %d distinct clusters, r = [%s] ===\n', ...
    selection_rule, h_best, scan(best_k).nclusters, sprintf('%.6g ', sort(scan(best_k).r)));

% Pick the TWO largest-support clusters at h_best (the least flukey candidates -- support is
% measured before polishing, so this selection cannot be biased by the polish outcome).
memberCounts = arrayfun(@(x) numel(x.member_seed_ids), rep_best.roots);
[~, orderBySupport] = sort(memberCounts, 'descend');
pick = orderBySupport(1:min(2,numel(orderBySupport)));
rootA = rep_best.roots(pick(1));
if numel(pick) >= 2
    rootB = rep_best.roots(pick(2));
else
    rootB = [];
end
fprintf('  root A: support=%d (seed_ids %s)\n', numel(rootA.member_seed_ids), ...
    strjoin(rootA.member_seed_ids,','));
if ~isempty(rootB)
    fprintf('  root B: support=%d (seed_ids %s)\n', numel(rootB.member_seed_ids), ...
        strjoin(rootB.member_seed_ids,','));
end

% ------------------------------------------------------------------------------------------
% Step 2 -- polish the two selected roots as tightly as the kernel allows.
% ------------------------------------------------------------------------------------------
[node_best, node_meta] = invz_ordered_make_node(ctx, h_best);
if isempty(node_best)
    error('invzp:S3:NodeUnavailable', 'node construction failed at h_best=%.8g: %s', ...
        h_best, node_meta.status);
end

tolLadder = [1e-8 1e-10 1e-12 1e-14];
polishA = local_polish(node_best, rootA.state, tolLadder);
fprintf('\n--- Polish root A (seed_ids: %s, support=%d) ---\n', ...
    strjoin(rootA.member_seed_ids,','), numel(rootA.member_seed_ids));
local_print_polish(polishA);

if ~isempty(rootB)
    polishB = local_polish(node_best, rootB.state, tolLadder);
    fprintf('\n--- Polish root B (seed_ids: %s, support=%d) ---\n', ...
        strjoin(rootB.member_seed_ids,','), numel(rootB.member_seed_ids));
    local_print_polish(polishB);
else
    polishB = [];
    fprintf('\nOnly one distinct cluster had support >0 for a second root; no root B to polish.\n');
end

% ------------------------------------------------------------------------------------------
% Step 3 -- separation vs backward error.
% ------------------------------------------------------------------------------------------
result = struct('schema','invzp_s3_twin_polish/v1', 'reproduced', true, ...
    'fixture', struct('T',T,'Bx',Bx,'nw',nw,'Jscale',ctx.Jscale, ...
        'nq',numel(J),'Jcc0',ci.Jcc0,'Jaa0',Jaa0), ...
    'seed_spec', struct('rng','twister:20260729','K0grid',K0grid,'sigScales',sigScales, ...
        'cluster_spec',spec.cluster), ...
    'hgrid', hgrid, 'scan', scan, 'h_best', h_best, 'selection_rule', selection_rule, ...
    'rootA', rootA, 'rootB', rootB, 'polishA', polishA);
if ~isempty(polishB)
    result.polishB = polishB;
end

if isempty(polishB) || ~polishA.accepted || ~polishB.accepted
    result.conclusion = 'not reproduced: no coexistence point was found in the scanned range';
    if polishA.accepted && (isempty(polishB) || ~polishB.accepted)
        result.conclusion = ['not reproduced: only one of the two candidate roots at h_best ' ...
            'remained A-D-accepted under polishing'];
    end
    fprintf('\nCONCLUSION: %s\n', result.conclusion);
    if ~isempty(save_path), save(save_path, 'result', 'reports'); end
    return
end

uA = [polishA.state.Sigma(:); polishA.state.K0s];
uB = [polishB.state.Sigma(:); polishB.state.K0s];
d2   = norm(uA-uB, 2);
dinf = norm(uA-uB, Inf);
dr   = abs(polishA.info.r - polishB.info.r);
worst_resid = max(polishA.residual_inf_final, polishB.residual_inf_final);
factor_2   = d2 / worst_resid;
factor_inf = dinf / worst_resid;
factor_tol = dinf / 1e-8;

result.separation = struct('d2', d2, 'dinf', dinf, 'dr', dr, ...
    'worst_polished_residual_inf', worst_resid, ...
    'factor_over_worst_residual_2norm', factor_2, ...
    'factor_over_worst_residual_infnorm', factor_inf, ...
    'factor_over_1e8_tol', factor_tol);

fprintf('\n=== Step 3: separation vs backward error ===\n');
fprintf('  ||uA-uB||_2   = %.6g\n', d2);
fprintf('  ||uA-uB||_inf = %.6g\n', dinf);
fprintf('  |rA-rB|       = %.6g   (rA=%.8g, rB=%.8g)\n', dr, polishA.info.r, polishB.info.r);
fprintf('  worst polished residual_inf (backward error) = %.6g\n', worst_resid);
fprintf('  separation(inf) / worst backward error = %.6g x\n', factor_inf);
fprintf('  separation(inf) / 1e-8 A-block tolerance = %.6g x\n', factor_tol);

% Threshold is exactly the wording of the two allowed forms: "exceeds both polished backward
% errors" means the separation is larger than EITHER individual residual_inf, i.e. larger than
% their max; "within the polished backward errors" is the complementary (<=) case. No extra
% safety margin is inserted here -- factor_inf/factor_tol are printed regardless so the reader
% can see exactly how close to 1x the call was.
if dinf > worst_resid
    result.conclusion = sprintf( ...
        'distinct: the separation exceeds both polished backward errors by a factor of %.4g', ...
        factor_inf);
else
    result.conclusion = 'not distinguishable: the separation is within the polished backward errors';
end
fprintf('\nCONCLUSION: %s\n', result.conclusion);

if ~isempty(save_path)
    save(save_path, 'result', 'reports');
    fprintf('saved: %s\n', save_path);
end
end

% ================================================================================================
function score = local_pair_score(rvals)
if numel(rvals) < 2
    score = Inf;
    return
end
best = Inf;
for i = 1:numel(rvals)
    for j = i+1:numel(rvals)
        s = min(abs(rvals(i)-0.768)+abs(rvals(j)-0.822), ...
                abs(rvals(j)-0.768)+abs(rvals(i)-0.822));
        best = min(best, s);
    end
end
score = best;
end

% ================================================================================================
function polish = local_polish(node, seedState, tolLadder)
%LOCAL_POLISH Re-solve from seedState with a large iteration budget at progressively tighter
% tol_outer. Returns the result at the TIGHTEST tol_outer for which acceptance still holds
% (never loosens below 1e-8; if even 1e-8 fails to re-accept, that failure is reported honestly).
polish = struct('tol_ladder', tolLadder, 'attempts', repmat(struct( ...
    'tol_outer',NaN,'accepted',false,'reason','','residual_inf',NaN,'iterations',NaN, ...
    'rcond',NaN,'raw_rcond',NaN,'equilibrated_rcond',NaN,'pole_margin',NaN,'mean_margin',NaN, ...
    'r',NaN), 1, numel(tolLadder)), ...
    'accepted', false, 'tightest_tol_outer', NaN, 'state', seedState, 'info', struct());
for k = 1:numel(tolLadder)
    opts = struct('tol_outer', tolLadder(k), 'maxit', 80, 'rcond_min', 1e-12);
    [state, info] = invz_ordered_node_newton(node, seedState, opts);
    a = struct('tol_outer', tolLadder(k), 'accepted', info.accepted, ...
        'reason', info.reason, 'residual_inf', info.residual_inf, ...
        'iterations', info.iterations, 'rcond', info.rcond, ...
        'raw_rcond', info.raw_rcond, 'equilibrated_rcond', info.equilibrated_rcond, ...
        'pole_margin', info.pole_margin, 'mean_margin', info.mean_margin, 'r', info.r);
    polish.attempts(k) = a;
    if info.accepted
        polish.accepted = true;
        polish.tightest_tol_outer = tolLadder(k);
        polish.state = state;
        polish.info = info;
        polish.residual_inf_final = info.residual_inf;
        % Keep trying tighter tolerances; do not break, so the ladder is fully reported.
    end
end
if ~polish.accepted
    % Report the loosest (1e-8) attempt's raw numbers even though it did not (re-)accept, so the
    % failure itself is visible rather than silently empty.
    polish.state = seedState;
    polish.info = struct('accepted', false, 'reason', polish.attempts(1).reason, ...
        'r', polish.attempts(1).r);
    polish.residual_inf_final = polish.attempts(1).residual_inf;
end
end

% ================================================================================================
function local_print_polish(p)
fprintf('  tol_outer ladder:\n');
for k = 1:numel(p.attempts)
    a = p.attempts(k);
    fprintf(['    tol_outer=%.0e  accepted=%d  reason=%-24s  resid_inf=%.6g  iters=%2d  ' ...
        'rcond=%.3g  raw_rcond=%.3g  eq_rcond=%.3g  pole_margin=%.3g  mean_margin=%.3g  r=%.8g\n'], ...
        a.tol_outer, a.accepted, a.reason, a.residual_inf, a.iterations, ...
        a.rcond, a.raw_rcond, a.equilibrated_rcond, a.pole_margin, a.mean_margin, a.r);
end
if p.accepted
    fprintf('  TIGHTEST ACCEPTED tol_outer = %.0e  (residual_inf = %.6g, r = %.8g)\n', ...
        p.tightest_tol_outer, p.residual_inf_final, p.info.r);
else
    fprintf('  DID NOT RE-ACCEPT at any tol_outer in the ladder (loosest failure reason: %s)\n', ...
        p.attempts(1).reason);
end
end
