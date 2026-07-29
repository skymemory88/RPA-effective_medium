function S = invzp_exec_s9_field_ladder(fields, save_path, opts)
%INVZP_EXEC_S9_FIELD_LADDER Field-wide damping-ladder census (execution packet S9).
%
% Execution packet S9 of docs/execution/invzp_plan_execution_diary.md
% ("PRIORITY REVISION 2026-07-29"). MEASUREMENT ONLY. For each transverse field
% it runs the production ordered column at a ladder of mix_outer values and
% records the FIRST rung that closes the column, or that none does. It changes
% no solver behaviour, wires nothing into production, and gates nothing.
%
% WHAT THIS CAN AND CANNOT ESTABLISH. A column closed by a damping rung yields a
% spectrum at an A-D-accepted state. It does NOT yield a certified equilibrium
% phase label: diary E4 found a 1e-12 arithmetic perturbation moving nodes to
% different accepted states, and E8 found 6-11 accepted root clusters at a
% single h near 4.05 T. Anything produced downstream of this census must carry
% the label `converged-under-damping-ladder, branch-not-certified`.
%
%   fields      vector of Bx in T (default: the S9 pilot set spanning 3.0-6.0 T)
%   save_path   optional .mat path; written after EVERY field so a long run can
%               be inspected while running and resumed after an interrupt
%   opts        .ladder        mix_outer values, tried in order (see below)
%               .max_outer     iteration budget per rung (default 1000)
%               .full_ladder   false (default): stop at the first rung that
%                              closes. true: run EVERY rung at every field, so
%                              the h* values different rungs converge to can be
%                              compared. The comparison is the point -- if two
%                              rungs close to different h*, the ladder is
%                              selecting a branch, not just removing a mask.
%               .resume        true (default): if save_path exists, skip fields
%                              already recorded in it
%
% LADDER ORDER IS EVIDENCE-BASED, NOT THE ORDER WRITTEN IN THE DIARY SPEC. The
% spec listed {0.7, 0.5, 0.4, 0.3, 0.2, 0.15}, i.e. least damping first, on the
% usual assumption that weak damping is fastest when it works. Diary E1 measured
% the opposite at 3.825 T: mix_outer = 0.70 was both the WORST (17/34 nodes
% failed) and the SLOWEST (202.6 s, because every failing node burns the full
% max_outer budget), while mix_outer = 0.30 closed the whole column in 23.4 s.
% Descending order would therefore pay ~200 s per field to learn nothing before
% reaching the rung that works. The ladder below leads with the rung E1 measured
% as best and cheapest. This changes the cost and the reported "first rung"
% label; it does not change which fields are recoverable by damping at all,
% which is the measurement S9 exists to make.
%
% The fixed conditions (T = 0.1 K, 16^3 grid, dpRng 30, BRUTEFORCE dipole) are
% inherited unchanged from invzp_exec_s1_failure_census, so every number here is
% directly comparable to E1. Note this is deliberately NOT the user's current
% driver setting of dipoleBackend = 'ewald': the two backends differ by up to
% 1.14e-4 meV at finite q (diary E7/WP0), and mixing them inside one census
% would confound a backend change with a damping change.
if nargin < 1 || isempty(fields)
    fields = [3.00 3.30 3.60 3.825 4.05 4.35 4.70 5.00 5.40 6.00];
end
if nargin < 2, save_path = ''; end
if nargin < 3 || isempty(opts), opts = struct(); end

ladder      = getf(opts, 'ladder', [0.30 0.40 0.20 0.50 0.15 0.70]);
max_outer   = getf(opts, 'max_outer', 1000);
full_ladder = getf(opts, 'full_ladder', false);
resume      = getf(opts, 'resume', true);

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(root, 'docs', 'execution'));

fields = fields(:).';
done = struct('Bx', {}, 'rungs', {}, 'first_ok_mix', {}, 'hstar', {}, ...
    'closed', {}, 'wall_s', {});
if resume && ~isempty(save_path) && isfile(save_path)
    L = load(save_path);
    if isfield(L, 'S') && isfield(L.S, 'rows') && ~isempty(L.S.rows)
        done = L.S.rows;
        fprintf('resuming: %d field(s) already recorded in %s\n', numel(done), save_path);
    end
end

fprintf('=== S9 field-wide damping-ladder census ===\n');
fprintf('ladder mix_outer = [%s], max_outer = %d, mode = %s\n', ...
    num2str(ladder), max_outer, ternary(full_ladder, 'FULL (all rungs)', 'early-exit'));
fprintf('%d field(s): %s\n\n', numel(fields), num2str(fields));

S = struct('fields', fields, 'ladder', ladder, 'max_outer', max_outer, ...
    'full_ladder', full_ladder, 'rows', done);

for i = 1:numel(fields)
    Bx = fields(i);
    if ~isempty(done) && any(abs([done.Bx] - Bx) < 1e-12)
        fprintf('Bx = %.4f T  [already recorded, skipped]\n', Bx);
        continue
    end
    t0 = tic;
    rungs = struct('mix_outer', {}, 'status', {}, 'n_failed', {}, ...
        'n_nodes', {}, 'hstar', {}, 'wall_s', {}, 'failed_ids', {});
    first_ok = NaN;  hstar = NaN;
    fprintf('Bx = %.4f T\n', Bx);
    for k = 1:numel(ladder)
        so = struct('mix_outer', ladder(k), 'max_outer', max_outer);
        tk = tic;
        try
            txt = evalc('Rk = invzp_exec_s1_failure_census(Bx, '''', so);'); %#ok<NASGU>
            st  = Rk.meta.hmf_status;
            nf  = Rk.meta.n_failed;
            nn  = Rk.meta.n_nodes;
            hs  = Rk.meta.hstar;
            fid = find(~Rk.tab.accepted).';
        catch ME
            st = ['error:' ME.identifier]; nf = NaN; nn = NaN; hs = NaN; fid = [];
        end
        wk = toc(tk);
        rungs(end+1) = struct('mix_outer', ladder(k), 'status', st, ...
            'n_failed', nf, 'n_nodes', nn, 'hstar', hs, 'wall_s', wk, ...
            'failed_ids', fid); %#ok<AGROW>
        fprintf('   mix=%.2f -> %-14s failed=%s/%s  hstar=%.7g  (%.1f s)\n', ...
            ladder(k), st, numstr(nf), numstr(nn), hs, wk);
        if strcmp(st, 'ok')
            if isnan(first_ok), first_ok = ladder(k); hstar = hs; end
            if ~full_ladder, break; end
        elseif strcmp(st, 'no_bare_order')
            % NOT a convergence failure and NOT damping-dependent. invz_hmf_ordered.m:201
            % returns this when the BARE single-ion mean-field solve does not order
            % (|<Jz>| <= 1e-6), which is decided before any outer Sigma<->K iteration
            % runs -- so mix_outer cannot influence it, and the remaining rungs would
            % reproduce it identically. Measured in the 2026-07-29 pilot: at 5.0/5.4/6.0 T
            % all six rungs returned no_bare_order with 0/0 nodes and identical
            % sub-second wall times. Stop the ladder here.
            break
        end
    end
    row = struct('Bx', Bx, 'rungs', rungs, 'first_ok_mix', first_ok, ...
        'hstar', hstar, 'closed', ~isnan(first_ok), 'wall_s', toc(t0), ...
        'class', local_classify(rungs, first_ok, hstar));
    if isempty(S.rows), S.rows = row; else, S.rows(end+1) = row; end
    switch row.class
        case 'ordered'
            fprintf('   => CLOSED at mix = %.2f with an ORDERED root, hstar = %.7g  (%.1f s)\n\n', ...
                first_ok, hstar, row.wall_s);
        case 'pm_no_root'
            fprintf('   => converged at mix = %.2f, NO nonzero root (PM side)  (%.1f s)\n\n', ...
                first_ok, row.wall_s);
        case 'pm_no_bare_order'
            fprintf('   => PM: bare single ion does not order (mix-independent)  (%.1f s)\n\n', ...
                row.wall_s);
        otherwise
            fprintf('   => NO RUNG CONVERGED -- GENUINE CONVERGENCE FAILURE  (%.1f s)\n\n', ...
                row.wall_s);
    end
    if ~isempty(save_path), save(save_path, 'S'); end
end

% ---- summary ---------------------------------------------------------------
r = S.rows;
if ~isempty(r)
    cls = {r.class};
    ord = strcmp(cls, 'ordered');
    pmr = strcmp(cls, 'pm_no_root');
    pmb = strcmp(cls, 'pm_no_bare_order');
    fail = strcmp(cls, 'failed');
    % The four classes are reported SEPARATELY and never summed into one
    % "closed" rate. Only 'failed' is a convergence failure; both PM classes are
    % the physically correct answer that no ordered root exists at that field,
    % and folding them into a failure count would overstate the residual set.
    fprintf('--- summary over %d field(s) ---\n', numel(r));
    fprintf('   ordered root found        : %d\n', sum(ord));
    fprintf('   PM, converged, no root    : %d  (status ok, hstar NaN)\n', sum(pmr));
    fprintf('   PM, bare does not order   : %d  (mix-independent, not a failure)\n', sum(pmb));
    fprintf('   GENUINE convergence fail  : %d\n', sum(fail));
    if any(ord)
        fprintf('\n   rung that first closed each ORDERED column:\n');
        for k = 1:numel(ladder)
            n = sum(ord & ([r.first_ok_mix] == ladder(k)));
            if n > 0
                fprintf('      mix = %.2f : %d field(s)\n', ladder(k), n);
            end
        end
    end
    if any(fail)
        fprintf('\n   RESIDUAL FAILURE SET (no rung converges): %s\n', num2str([r(fail).Bx]));
    else
        fprintf('\n   RESIDUAL FAILURE SET (no rung converges): EMPTY\n');
    end
    if any(ord)
        fprintf('   ordered-root field range: %.4f to %.4f T\n', ...
            min([r(ord).Bx]), max([r(ord).Bx]));
    end
    fprintf('   total wall %.1f s (%.1f s/field mean)\n', sum([r.wall_s]), mean([r.wall_s]));
    if full_ladder
        fprintf('\n--- h* agreement across rungs that closed the same field ---\n');
        for i = 1:numel(r)
            ok = strcmp({r(i).rungs.status}, 'ok');
            if sum(ok) < 2, continue; end
            hs = [r(i).rungs(ok).hstar];
            fprintf('   Bx = %.4f T: %d rungs closed, h* spread = %.3e (values %s)\n', ...
                r(i).Bx, sum(ok), max(hs) - min(hs), num2str(hs, '%.7g '));
        end
    end
end
if ~isempty(save_path), save(save_path, 'S'); fprintf('\nsaved: %s\n', save_path); end
end

function c = local_classify(rungs, first_ok, hstar)
%LOCAL_CLASSIFY Four-way column outcome. Only 'failed' is a convergence failure.
%
%   ordered           some rung returned status 'ok' AND a finite nonzero root h*.
%                     This is the deliverable: chi is computable on an ordered branch.
%   pm_no_root        some rung returned 'ok' but h* is NaN. invz_hmf_ordered.m:380
%                     returns here when the bracket contains no nonzero root -- the
%                     column CONVERGED and the correct answer is "no ordered state".
%   pm_no_bare_order  every rung returned 'no_bare_order' (invz_hmf_ordered.m:201):
%                     the bare single ion does not order, decided before any outer
%                     iteration, hence independent of mix_outer.
%   failed            no rung converged for any other reason (node_failed etc.).
%                     THIS is the set the damping ladder exists to shrink.
if ~isnan(first_ok)
    if isfinite(hstar) && hstar ~= 0, c = 'ordered'; else, c = 'pm_no_root'; end
elseif ~isempty(rungs) && all(strcmp({rungs.status}, 'no_bare_order'))
    c = 'pm_no_bare_order';
else
    c = 'failed';
end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
function s = numstr(x)
if isnan(x), s = '--'; else, s = sprintf('%d', x); end
end
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
