function O = invzp_exec_s10_overshoot(save_path)
%INVZP_EXEC_S10_OVERSHOOT Would a step-acceptance rule have prevented the excursion?
%
% Execution packet S10 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY, OFFLINE -- reads saved S10 traces, runs no solver.
%
% CONTEXT. Two candidate fixes for the 1.0 T mask were on the table.
%   (1) Pole-regular arithmetic (closure_coordinate = 'defactored'). REFUTED by
%       measurement: at mix_outer = 0.30 it returns the IDENTICAL failure set,
%       11/34, node-for-node [1 2 3 4 5 6 7 19 22 23 24]. The reason is visible in
%       the pole-locality data -- only ~2% of iterations have |d0| < 0.1, so the
%       iterate almost never sits NEAR the pole; it leaps across it. Making the
%       arithmetic finite at a point the iterate barely visits fixes nothing.
%   (2) Unbounded overshoot. The accepted solutions live at |d0| in [0.48, 1.61]
%       while the path reaches |d0| ~ 1.8e5 -- five orders of magnitude outside the
%       region containing the answer, with nothing bounding the step.
%
% THIS PROBE TESTS (2) BEFORE ANY CODE IS WRITTEN FOR IT. If overshoot is the
% binding mechanism, then the steps that carry the iterate outward must be steps on
% which the map residual GETS WORSE. A standard globalization (backtracking line
% search: accept a step only if |G(u)-u| decreases, else shorten it) would then have
% rejected exactly those steps. If instead the residual mostly decreases while |d0|
% explodes, the excursion is not a simple overshoot and a line search would not have
% caught it -- in which case this fix must be dropped too, not implemented hopefully.
%
% WHY A STEP RULE IS NOT A POLE FLOOR. Repo standards forbid pole floors, clipped
% denominators and masked-to-bare substitution. A backtracking line search is
% categorically different: it does not modify the equation, floor any denominator,
% or alter the residual whose zero is being sought. It changes only the PATH taken
% to a fixed point, so the fixed-point SET is provably unchanged -- which is exactly
% the property that makes it admissible here, and it must be stated in any proposal
% that follows from this measurement.
%
% Reports, per (alpha, node) over the S10 traces:
%   - fraction of outer steps on which |resid_map| INCREASED
%   - the same fraction restricted to steps that moved |d0| outward past 10
%   - the largest single-step growth factor in |resid_map|
%   - correlation between step-wise growth in |d0| and growth in |resid_map|
if nargin < 1, save_path = ''; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outd = fullfile(root, 'docs', 'execution', 'out');

f = dir(fullfile(outd, 's10_B1_a*.mat'));
f = f(~contains({f.name}, 'summary'));
if isempty(f), error('no S10 traces found in %s', outd); end

rows = struct('file', {}, 'node', {}, 'n', {}, 'frac_worse', {}, ...
    'frac_worse_far', {}, 'max_growth', {}, 'corr_d0_resid', {}, ...
    'd0_max', {}, 'resid_final', {});
fprintf('=== S10 overshoot probe: do outward steps make the residual worse? ===\n\n');
for k = 1:numel(f)
    L = load(fullfile(outd, f(k).name), 'trace');
    it = L.trace.iters;
    nid = [it.node_id];
    for c = unique(nid)
        sel = it(nid == c);
        if numel(sel) < 200, continue; end
        r  = abs([sel.resid_map].');
        d0 = abs([sel.gstat_local_denom].');
        good = isfinite(r) & isfinite(d0);
        r = r(good); d0 = d0(good);
        if numel(r) < 200, continue; end
        gr = r(2:end)./max(r(1:end-1), realmin);     % step-wise residual growth
        gd = d0(2:end)./max(d0(1:end-1), realmin);   % step-wise |d0| growth
        far = d0(2:end) > 10;                        % steps landing far outside
        lg = log10(max(gr, realmin)); ld = log10(max(gd, realmin));
        ok = isfinite(lg) & isfinite(ld);
        cc = NaN;
        if nnz(ok) > 10
            C = corrcoef(ld(ok), lg(ok));
            cc = C(1,2);
        end
        rows(end+1) = struct('file', f(k).name, 'node', c, 'n', numel(r), ...
            'frac_worse', mean(gr > 1), ...
            'frac_worse_far', ternary_num(any(far), mean(gr(far) > 1), NaN), ...
            'max_growth', max(gr), 'corr_d0_resid', cc, ...
            'd0_max', max(d0), 'resid_final', r(end)); %#ok<AGROW>
    end
end

fprintf('%-22s %-5s %-7s %-11s %-13s %-12s %-11s %s\n', 'file', 'node', 'n', ...
    'frac_worse', 'frac_worse_far', 'max_growth', 'corr(d0,r)', 'd0_max');
for i = 1:numel(rows)
    fprintf('%-22s %-5d %-7d %-11.4f %-13s %-12.3e %-11.3f %.3e\n', ...
        rows(i).file, rows(i).node, rows(i).n, rows(i).frac_worse, ...
        fmtnum(rows(i).frac_worse_far), rows(i).max_growth, ...
        rows(i).corr_d0_resid, rows(i).d0_max);
end

fw  = [rows.frac_worse];
fwf = [rows.frac_worse_far];  fwf = fwf(isfinite(fwf));
cc  = [rows.corr_d0_resid];   cc  = cc(isfinite(cc));
fprintf('\n--- aggregate over %d (file,node) pairs ---\n', numel(rows));
fprintf('  fraction of ALL steps that worsen the residual : median %.4f  [%.4f, %.4f]\n', ...
    median(fw), min(fw), max(fw));
if ~isempty(fwf)
    fprintf('  same, restricted to steps landing at |d0| > 10 : median %.4f  [%.4f, %.4f]\n', ...
        median(fwf), min(fwf), max(fwf));
end
if ~isempty(cc)
    fprintf('  corr( log10 d|d0|, log10 d|resid| )           : median %+.3f  [%+.3f, %+.3f]\n', ...
        median(cc), min(cc), max(cc));
end
fprintf(['\nREADING. A line search rejects a step when the residual grows. It can only\n' ...
         'help in proportion to how often that happens, and it targets the excursion\n' ...
         'only if outward steps are preferentially the worsening ones (frac_worse_far\n' ...
         'clearly above frac_worse, and/or a positive d0/resid growth correlation).\n' ...
         'If frac_worse is near 0.5 with no outward preference, the walk is a random\n' ...
         'wander rather than a systematic overshoot, and a line search is NOT the fix.\n']);

O = rows;
if ~isempty(save_path), save(save_path, 'O'); fprintf('\nsaved: %s\n', save_path); end
end

function v = ternary_num(c, a, b), if c, v = a; else, v = b; end, end
function s = fmtnum(x)
if isnan(x), s = '--'; else, s = sprintf('%.4f', x); end
end
