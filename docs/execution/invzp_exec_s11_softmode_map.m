function S = invzp_exec_s11_softmode_map(fields, save_path, opts)
%INVZP_EXEC_S11_SOFTMODE_MAP Map the finite-q soft mode across the field axis.
%
% Execution packet S11 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY. Supersedes the S9 damping-ladder census, which was stopped
% after 2 of 79 fields once E9b established that its question no longer mattered.
%
% WHY THIS REPLACES S9. S9 asked "which mix_outer rung closes each column", to
% justify an adaptive-damping production mode. E9/E9b withdrew that mode: a rung
% that closes a column only shows the iterate found a path AROUND the finite-q
% soft mode, never that the column is free of it, so the answer would not be acted
% on. S9 also cost up to SIX solver runs per field (1.0 T: 962.7 s; 1.5 T:
% 1554.7 s, both closing under no rung) to produce that unusable answer.
%
% The live question after E9b is where the soft mode BITES, and its observable is
% already recorded per node by the existing census:
%
%   finite Dq_absmin  <=>  node fails
%
% measured without exception in either direction over ~100 node evaluations across
% three iteration configurations and two fields (1.0 T: 11 of 11 failures carry it,
% 0 of 23 accepted; 3.825 T mix 0.40: 1 of 1, 0 of 33; 3.825 T mix 0.70: 22 of 22,
% 0 of 12). That needs ONE run per field, not six, and yields strictly more: which
% nodes, at which h, and how close D(q) = 1 + (J(q) - K0)*Gstat gets to zero.
%
% NOTE ON Dq_min VERSUS Dq_absmin. D(q) merely going negative somewhere is NOT the
% predictor: 27 of 34 nodes at 3.825 T have Dq_min < 0 INCLUDING accepted ones. It
% is D(q) approaching zero -- a finite Dq_absmin -- that corresponds to failure.
% Both are recorded here; do not substitute one for the other.
%
% WHAT THIS CANNOT SETTLE. It locates the soft mode in (Bx, h); it does not decide
% whether the instability is (a) a real competing ordering wavevector at low field,
% (b) an artefact of the 1/z truncation's q-dependence, or (c) the 1.5 T fold
% region of diagnosis Sec 6.1. That is the open physical question and needs the
% wavevector itself, not just the magnitude of the minimum.
%
%   fields      Bx in T (default: coarse low field + dense production window)
%   save_path   .mat written after EVERY field; resumable
%   opts        .mix_outer (0.30, the best single rung measured at 3.825 T)
%               .max_outer (1000)
%               .resume    (true)
if nargin < 1 || isempty(fields)
    fields = [1.0:0.25:2.75, 3.0:0.05:4.75, 5.0 5.5];
end
if nargin < 2, save_path = ''; end
if nargin < 3 || isempty(opts), opts = struct(); end
mix  = getf(opts, 'mix_outer', 0.30);
maxo = getf(opts, 'max_outer', 1000);
resume = getf(opts, 'resume', true);

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); addpath(fullfile(root, 'docs', 'execution'));

fields = fields(:).';
rows = struct('Bx', {}, 'status', {}, 'hstar', {}, 'n_nodes', {}, 'n_failed', {}, ...
    'n_soft', {}, 'soft_h_lo', {}, 'soft_h_hi', {}, 'dq_absmin_min', {}, ...
    'dq_min_min', {}, 'n_dqneg', {}, 'soft_eq_fail', {}, 'wall_s', {});
if resume && ~isempty(save_path) && isfile(save_path)
    L = load(save_path);
    if isfield(L, 'S') && isfield(L.S, 'rows') && ~isempty(L.S.rows)
        rows = L.S.rows;
        fprintf('resuming: %d field(s) already recorded\n', numel(rows));
    end
end

fprintf('=== S11 soft-mode map (mix_outer = %.2f, max_outer = %d) ===\n', mix, maxo);
fprintf('%d field(s)\n\n', numel(fields));
S = struct('mix_outer', mix, 'max_outer', maxo, 'rows', rows);

for Bx = fields
    if ~isempty(rows) && any(abs([rows.Bx] - Bx) < 1e-12), continue; end
    t0 = tic;
    so = struct('mix_outer', mix, 'max_outer', maxo);
    try
        txt = evalc('R = invzp_exec_s1_failure_census(Bx, '''', so);'); %#ok<NASGU>
        tb = R.tab;
        acc  = logical(tb.accepted);
        soft = isfinite(tb.Dq_absmin);          % the soft-mode marker
        row = struct('Bx', Bx, 'status', R.meta.hmf_status, 'hstar', R.meta.hstar, ...
            'n_nodes', height(tb), 'n_failed', nnz(~acc), 'n_soft', nnz(soft), ...
            'soft_h_lo', minor(tb.h(soft)), 'soft_h_hi', maxor(tb.h(soft)), ...
            'dq_absmin_min', minor(tb.Dq_absmin(soft)), ...
            'dq_min_min', minor(tb.Dq_min), 'n_dqneg', nnz(tb.Dq_min < 0), ...
            'soft_eq_fail', isequal(find(soft(:)).', find(~acc(:)).'), ...
            'wall_s', toc(t0));
    catch ME
        row = struct('Bx', Bx, 'status', ['error:' ME.identifier], 'hstar', NaN, ...
            'n_nodes', NaN, 'n_failed', NaN, 'n_soft', NaN, 'soft_h_lo', NaN, ...
            'soft_h_hi', NaN, 'dq_absmin_min', NaN, 'dq_min_min', NaN, ...
            'n_dqneg', NaN, 'soft_eq_fail', false, 'wall_s', toc(t0));
    end
    if isempty(S.rows), S.rows = row; else, S.rows(end+1) = row; end
    rows = S.rows;
    fprintf(['Bx=%.4f  %-13s hstar=%-12.7g failed=%2d/%2d soft=%2d  ' ...
             'soft_h=[%.3g,%.3g]  min|Dq|=%-10.4g Dqmin=%-9.4g  soft==fail:%d  (%.0f s)\n'], ...
        row.Bx, row.status, row.hstar, row.n_failed, row.n_nodes, row.n_soft, ...
        row.soft_h_lo, row.soft_h_hi, row.dq_absmin_min, row.dq_min_min, ...
        row.soft_eq_fail, row.wall_s);
    if ~isempty(save_path), save(save_path, 'S'); end
end

r = S.rows;
if ~isempty(r)
    softy = [r.n_soft] > 0;
    fprintf('\n--- summary over %d field(s) ---\n', numel(r));
    fprintf('  fields with a soft mode at any node : %d\n', nnz(softy));
    if any(softy)
        fprintf('  soft-mode field range              : %.4f to %.4f T\n', ...
            min([r(softy).Bx]), max([r(softy).Bx]));
        fprintf('  smallest min|D(q)| anywhere        : %.4g\n', ...
            min([r(softy).dq_absmin_min]));
    end
    % The identity is the load-bearing claim; report any field where it breaks
    % rather than assuming it holds everywhere it was not checked.
    bad = find(~[r.soft_eq_fail] & isfinite([r.n_failed]));
    if isempty(bad)
        fprintf('  "finite Dq_absmin <=> node fails" held at EVERY field measured\n');
    else
        fprintf('  *** identity BROKE at %d field(s): %s ***\n', ...
            numel(bad), num2str([r(bad).Bx]));
    end
    fprintf('  total wall %.0f s (%.0f s/field mean)\n', sum([r.wall_s]), mean([r.wall_s]));
end
if ~isempty(save_path), save(save_path, 'S'); fprintf('\nsaved: %s\n', save_path); end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function v = minor(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = min(x); end, end
function v = maxor(x), x = x(isfinite(x)); if isempty(x), v = NaN; else, v = max(x); end, end
