function P = invzp_exec_s10_pole_locality(save_path)
%INVZP_EXEC_S10_POLE_LOCALITY Is the Gstat pole ON the solution, or only on the path?
%
% Execution packet S10 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY, OFFLINE -- reads .mat files already produced by S1/S10 and runs
% no solver.
%
% THE QUESTION AND WHY IT DECIDES THE FIX. S10's eigenvalue probe showed the 1.0 T
% iterate crosses gstat_local_denom = 0 repeatedly. Two very different situations
% produce that, and they call for opposite responses:
%
%   (a) PATH-ONLY. The converged/accepted states sit well away from d0 = 0, and only
%       the ITERATION PATH wanders across it. Then the pole is an artefact of how the
%       solver travels, not of the solution, and the fix is algorithmic: pole-regular
%       arithmetic, a Newton step, or continuation -- any route that does not traverse
%       the singularity. Nothing physical is at stake.
%
%   (b) ON-SOLUTION. The accepted states themselves sit at small |d0|. Then
%       1 + Sigma0 + K0*G0inel0 -> 0 is a property of the state being sought -- a
%       divergent local static response of the effective medium -- and no change of
%       arithmetic makes it go away. That is a physics finding about the ordered
%       solution at low field, and it would have to be understood before any spectrum
%       computed there is trusted, however cleanly the solver converges.
%
% The discriminator: compare the distribution of |d0| at ACCEPTED nodes against the
% path excursion at FAILED nodes, at a field that converges (3.825 T) and at one that
% does not (1.0 T).
%
% Reads whichever of these exist:
%   docs/execution/out/s1_census_3p825.mat     (S1 census, field that converges)
%   docs/execution/out/s10_B1p0_a*.mat         (S10 traces, field that does not)
if nargin < 1, save_path = ''; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outd = fullfile(root, 'docs', 'execution', 'out');

P = struct('sets', []);
sets = {};

% ---- fields that converge: accepted-node |d0| from the S1 census table ------------
c1 = fullfile(outd, 's1_census_3p825.mat');
if isfile(c1)
    L = load(c1, 'R');
    tb = L.R.tab;
    d0 = tb.gdenom;
    acc = logical(tb.accepted);
    sets{end+1} = report('3.825 T (converges)', d0(acc), d0(~acc)); %#ok<AGROW>
end

% ---- field that does not converge: accepted vs failed, plus PATH excursion --------
% num2str(1.0) is '1', not '1.0', so the S10 traces are named s10_B1_a*.mat.
f = dir(fullfile(outd, 's10_B1_a*.mat'));
f = f(~contains({f.name}, 'summary'));
for k = 1:numel(f)
    fp = fullfile(outd, f(k).name);
    v = whos('-file', fp);
    names = {v.name};
    if ismember('R', names)
        L = load(fp, 'R');
        tb = L.R.tab;
        acc = logical(tb.accepted);
        sets{end+1} = report(sprintf('1.0 T node values [%s]', f(k).name), ...
            tb.gdenom(acc), tb.gdenom(~acc)); %#ok<AGROW>
    end
    if ismember('trace', names)
        L2 = load(fp, 'trace');
        gd = [L2.trace.iters.gstat_local_denom].';
        gd = gd(isfinite(gd));
        fprintf(['PATH excursion over ALL outer iterations [%s]: min %.4g  max %.4g  ' ...
                 'median|d0| %.4g  frac |d0| < 0.1: %.3f  sign changes: %d\n'], ...
            f(k).name, min(gd), max(gd), median(abs(gd)), mean(abs(gd) < 0.1), ...
            nnz(diff(sign(gd)) ~= 0));
    end
end

P.sets = sets;
if ~isempty(save_path), save(save_path, 'P'); fprintf('\nsaved: %s\n', save_path); end
end

function s = report(tag, d_acc, d_fail)
d_acc = d_acc(isfinite(d_acc));
d_fail = d_fail(isfinite(d_fail));
fprintf('\n--- %s ---\n', tag);
if ~isempty(d_acc)
    fprintf('  ACCEPTED nodes (n=%d): |d0| min %.4g  median %.4g  max %.4g;  frac |d0| < 0.1 = %.3f\n', ...
        numel(d_acc), min(abs(d_acc)), median(abs(d_acc)), max(abs(d_acc)), mean(abs(d_acc) < 0.1));
    fprintf('                        signed range [%.4g, %.4g]\n', min(d_acc), max(d_acc));
else
    fprintf('  ACCEPTED nodes: none\n');
end
if ~isempty(d_fail)
    fprintf('  FAILED   nodes (n=%d): |d0| min %.4g  median %.4g  max %.4g\n', ...
        numel(d_fail), min(abs(d_fail)), median(abs(d_fail)), max(abs(d_fail)));
else
    fprintf('  FAILED   nodes: none\n');
end
s = struct('tag', tag, 'd_accepted', d_acc, 'd_failed', d_fail);
end
