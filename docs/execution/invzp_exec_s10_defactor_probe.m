function S = invzp_exec_s10_defactor_probe(Bx, mixes, save_path)
%INVZP_EXEC_S10_DEFACTOR_PROBE Does the pole-regular coordinate clear the deep-ordered mask?
%
% Execution packet S10 of docs/execution/invzp_plan_execution_diary.md.
% MEASUREMENT ONLY. S10's eigenvalue probe showed the deep-ordered (1.0 T) failures
% are NOT a damping problem: the tail ratio never settles (IQR 0.9-2.5 on all ten
% commonly-failing nodes, so no local linear regime exists), while
% gstat_local_denom sweeps through ZERO with excursions to [-7550, +1450] and
% 550-900 sign flips per 1000 iterations. The iterate crosses its own Gstat pole.
%
% At that crossing the PRODUCTION DEFAULT path evaluates
%     Gq = Gs./(1 + (Jf - K0).*Gs)          (invz_emt_static_ordered.m:86,148)
% with Gs -> +-Inf, i.e. Inf/(+-Inf) = NaN. invz_gstat_ordered.m:35-38 documents
% exactly this and supplies the cure: the reassociation 1/Gtil0 = 1/Gstat - K0 is
% finite from both sides (Gtil0 -> -1/K0, r -> -G0bare*K0), and its own header warns
% the default arrangement "would turn a removable singularity into a node failure".
%
% closure_coordinate = 'defactored' selects that pole-regular arithmetic (it sets
% stable_form = true AND routes the q-average through
% invz_reciprocal_static_closure, whose weights scale./(z + Jf - K0) are regular at
% z = 0). S2 tested it only at 3.825 T, where the excursion is +-2.9, and found it
% changed nothing node-for-node -- a null result that is true THERE and does not
% generalise to a field where the excursion is 2500x larger.
%
% This probe runs matched factored/defactored pairs at one field and reports the
% failure sets. It changes no default and gates nothing.
%
% WHAT A POSITIVE RESULT WOULD AND WOULD NOT MEAN. Clearing these nodes would be a
% NUMERICAL fix to a REMOVABLE singularity -- it makes chi computable where it was
% masked. It would NOT certify the branch: E4 (1e-12 perturbations move nodes to
% different accepted states) and E8 (6-11 accepted clusters at one h) are untouched
% by it, so anything produced downstream stays branch-not-certified.
if nargin < 1 || isempty(Bx), Bx = 1.0; end
if nargin < 2 || isempty(mixes), mixes = [0.30 0.15]; end
if nargin < 3, save_path = ''; end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); addpath(fullfile(root, 'docs', 'execution'));

coords = {'factored', 'defactored'};
rows = struct('Bx', {}, 'mix_outer', {}, 'coord', {}, 'status', {}, ...
    'n_failed', {}, 'n_nodes', {}, 'hstar', {}, 'failed_ids', {}, 'wall_s', {});

fprintf('=== S10 defactored-coordinate probe, Bx = %.4f T ===\n\n', Bx);
for m = mixes(:).'
    for b = 1:numel(coords)
        c = struct('mix_outer', m, 'max_outer', 1000);
        c.emt_static = struct('closure_coordinate', coords{b});
        t0 = tic;
        txt = evalc('R = invzp_exec_s1_failure_census(Bx, '''', c);'); %#ok<NASGU>
        w = toc(t0);
        fid = find(~R.tab.accepted).';
        rows(end+1) = struct('Bx', Bx, 'mix_outer', m, 'coord', coords{b}, ...
            'status', R.meta.hmf_status, 'n_failed', R.meta.n_failed, ...
            'n_nodes', R.meta.n_nodes, 'hstar', R.meta.hstar, ...
            'failed_ids', fid, 'wall_s', w); %#ok<AGROW>
        fprintf('mix=%.2f coord=%-11s -> %-13s failed=%2d/%2d  hstar=%-12.8g (%.1f s)  ids=[%s]\n', ...
            m, coords{b}, R.meta.hmf_status, R.meta.n_failed, R.meta.n_nodes, ...
            R.meta.hstar, w, num2str(fid));
    end
end

fprintf('\n--- effect of the pole-regular coordinate, per matched pair ---\n');
for m = mixes(:).'
    f = rows([rows.mix_outer] == m & strcmp({rows.coord}, 'factored'));
    d = rows([rows.mix_outer] == m & strcmp({rows.coord}, 'defactored'));
    if isempty(f) || isempty(d), continue; end
    cleared = setdiff(f.failed_ids, d.failed_ids);
    added   = setdiff(d.failed_ids, f.failed_ids);
    fprintf('mix=%.2f : %d -> %d failures   cleared [%s]   newly failing [%s]\n', ...
        m, f.n_failed, d.n_failed, num2str(cleared), num2str(added));
    if isfinite(f.hstar) && isfinite(d.hstar)
        % Both coordinates are algebraically identical (S2 gate G1, rel err <= 1e-11),
        % so a LARGE h* difference here would not be roundoff -- it would mean the two
        % arithmetics landed on different roots, which is a branch-selection finding
        % and must not be reported as agreement.
        fprintf('         h* factored %.10g  defactored %.10g  |rel diff| %.3e\n', ...
            f.hstar, d.hstar, abs(d.hstar - f.hstar)/max(abs(f.hstar), realmin));
    end
end

S = rows;
if ~isempty(save_path), save(save_path, 'S'); fprintf('\nsaved: %s\n', save_path); end
end
