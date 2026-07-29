function A = invzp_exec_s1_iter_anatomy(mat_path, node_id)
%INVZP_EXEC_S1_ITER_ANATOMY Outer-iteration anatomy of one traced ordered node.
%
% Execution packet S1 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY, offline: reads a saved invz_ordered_trace .mat and
% characterises the outer Sigma<->K map history at one node so the failure mode
% can be classified as (a) a period-2 / period-p limit cycle, (b) a slow
% monotone crawl, or (c) an aperiodic wander between basins.
%
% Reports, on the outer-iteration residual resid_map (= dS, the raw map step):
%   - first/last/min values and the iteration of the minimum
%   - the mean log10 contraction rate over the last 200 iterations
%   - lag-1..lag-8 autocorrelation of the SIGNED increment of K0, and of
%     log10(resid_map), from which a dominant period is read off
%   - the count of sign alternations of diff(K0) over the tail
%   - the excursion range of gstat_local_denom, Gstat and Dq_abs_min
if nargin < 2 || isempty(node_id), node_id = []; end
S = load(mat_path, 'trace');
it = S.trace.iters;
nid = [it.node_id];
if isempty(node_id)
    % default: the node with the most outer iterations (the binding one)
    u = unique(nid);
    cnt = arrayfun(@(k) nnz(nid == k), u);
    [~, j] = max(cnt);
    node_id = u(j);
end
sel = it(nid == node_id);
n = numel(sel);
if n == 0, error('no iteration records for node %d', node_id); end

res  = [sel.resid_map].';
rst  = [sel.resid_static].';
K0   = [sel.K0].';
Gs   = [sel.Gstat].';
gd   = [sel.gstat_local_denom].';
dqa  = [sel.Dq_abs_min].';
dqn  = [sel.Dq_neg_count].';
yr   = [sel.y_rank].';

tail = max(1, n-199):n;
lr   = log10(max(res(tail), realmin));
p    = polyfit((1:numel(tail)).', lr, 1);

dK   = diff(K0);
tK   = dK(max(1, numel(dK)-199):end);
alt  = nnz(diff(sign(tK)) ~= 0);

lags = 1:8;
acK  = zeros(size(lags));
acR  = zeros(size(lags));
xK   = tK - mean(tK);
xR   = lr - mean(lr);
for k = lags
    acK(k) = sum(xK(1:end-k).*xK(1+k:end)) / max(sum(xK.^2), realmin);
    acR(k) = sum(xR(1:end-k).*xR(1+k:end)) / max(sum(xR.^2), realmin);
end
[~, pK] = max(acK);
[~, pR] = max(acR);

[rmin, imin] = min(res);

A = struct('node_id', node_id, 'n_iters', n, ...
    'resid_map_first', res(1), 'resid_map_last', res(end), ...
    'resid_map_min', rmin, 'resid_map_argmin', imin, ...
    'resid_static_last', rst(end), ...
    'log10_slope_per_iter_tail200', p(1), ...
    'sign_alternations_dK0_tail200', alt, ...
    'autocorr_dK0_lag1to8', acK, 'peak_lag_dK0', pK, ...
    'autocorr_log10resid_lag1to8', acR, 'peak_lag_log10resid', pR, ...
    'K0_tail_min', min(K0(tail)), 'K0_tail_max', max(K0(tail)), ...
    'Gstat_tail_min', min(Gs(tail)), 'Gstat_tail_max', max(Gs(tail)), ...
    'gdenom_tail_min', min(gd(tail)), 'gdenom_tail_max', max(gd(tail)), ...
    'gdenom_global_min', min(gd), ...
    'Dq_abs_min_tail_min', min(dqa(tail)), ...
    'Dq_neg_count_tail_range', [min(dqn(tail)) max(dqn(tail))], ...
    'y_rank_tail_range', [min(yr(tail)) max(yr(tail))]);

fprintf('=== iteration anatomy: node %d, %d outer iterations ===\n', node_id, n);
fprintf('resid_map  first %.6g  last %.6g  min %.6g @ it %d\n', ...
    res(1), res(end), rmin, imin);
fprintf('resid_static last %.6g\n', rst(end));
fprintf('log10(resid_map) slope over last 200 its: %+.3e per iteration\n', p(1));
fprintf('sign alternations of dK0 over last 200 steps: %d / %d\n', alt, numel(tK)-1);
fprintf('autocorr dK0        lags 1-8: %s  (peak lag %d)\n', mat2str(acK, 3), pK);
fprintf('autocorr log10resid lags 1-8: %s  (peak lag %d)\n', mat2str(acR, 3), pR);
fprintf('tail ranges: K0 [%.6g %.6g]  Gstat [%.6g %.6g]  gdenom [%.6g %.6g]\n', ...
    A.K0_tail_min, A.K0_tail_max, A.Gstat_tail_min, A.Gstat_tail_max, ...
    A.gdenom_tail_min, A.gdenom_tail_max);
fprintf('gdenom global min over the whole node history: %.6g\n', A.gdenom_global_min);
fprintf('Dq_abs_min tail min: %.6g   Dq_neg_count tail [%g %g]   y_rank tail [%g %g]\n', ...
    A.Dq_abs_min_tail_min, A.Dq_neg_count_tail_range(1), A.Dq_neg_count_tail_range(2), ...
    A.y_rank_tail_range(1), A.y_rank_tail_range(2));
fprintf('--- last 24 iterations (outer, resid_map, resid_static, K0, gdenom, Gstat) ---\n');
k0 = max(1, n-23);
for k = k0:n
    fprintf('%6d  %12.6g  %12.6g  %14.10g  %12.6g  %12.6g\n', ...
        sel(k).outer, res(k), rst(k), K0(k), gd(k), Gs(k));
end
end
