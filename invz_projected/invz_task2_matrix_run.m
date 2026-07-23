function summary = invz_task2_matrix_run(cells, opts)
%INVZ_TASK2_MATRIX_RUN Executes + incrementally checkpoints a Task-2 discriminator-matrix cell
% list (stage-2c task 2b-driver; cells = invz_task2_matrix_enumerate.m's output). RESUMABLE:
% loads any existing opts.results_path, SKIPS cells whose cfg_id is already present, and
% appends+saves the results array to disk after EVERY newly-run cell -- so a crash/kill loses
% at most the in-flight cell, and a re-run continues from the last checkpoint. This function
% only RUNS + CHECKPOINTS: the prereg Sec. A/C/D/F classification/analysis is 2b-report's job
% (not performed here, beyond what invz_task2_run_config/invz_task2_classify already record
% per node).
%
% cells  struct array from invz_task2_matrix_enumerate.m (or any subset/superset sharing its
%        cfg_id/depends_on/resolve/cfg schema).
% opts.results_path  (default '.superpowers/sdd/task2_matrix_results.mat', a git-ignored
%        scratch path -- see invz_task2_matrix.m's header). NOTE for tests: ALWAYS override
%        this to a temporary path; never let a test touch the real controller checkpoint.
% opts.cell_filter   (default {} = run every cell in `cells`; else a cellstr of cfg_ids to
%        run). A filtered isolated cell's OWN prerequisite (c.depends_on) is pulled in and
%        run/skipped automatically even when NOT itself listed in cell_filter -- see
%        ensure_cell_done below -- so a narrow filter can never silently omit a needed
%        dependency and corrupt an isolated cfg's h_list/seed.
% opts.ion           (default invz_ion()) -- forwarded to invz_task2_resolve_cell_cfg.m for
%        coupling-variant materialization (G3/G4/G5/G6); should match the SAME ion every
%        cell's own cfg.ion already carries.
% opts.verbose       (default true) -- one fprintf line per newly-run cell.
%
% Error/failure handling (three independent, non-crashing paths, all checkpointed as a
% record and NEVER abort the loop):
%  (a) invz_task2_run_config's OWN error policy (stage-2c task 2a review-fix pass): a genuine
%      'invz:*' physics/numerics exception from the per-config solve is ALREADY absorbed
%      there into out.status='failed' -- this function does not re-implement that, it just
%      checkpoints whatever `out` invz_task2_run_config returns, success or failure alike.
%  (b) A cell's OWN prerequisite (c.depends_on) itself failed or produced no nodes (so an
%      isolated cell's h_list/seed cannot be derived): checkpointed here as a structured
%      failed-dependency record (invz_task2_resolve_cell_cfg's dep_err), WITHOUT ever calling
%      invz_task2_run_config for that cell.
%  (c) A coupling-variant materialization (G3/G4/G5/G6) itself raises a genuine 'invz:*'
%      exception (see invz_task2_resolve_cell_cfg.m): folded into the same failed-dependency
%      record; anything NON-'invz:*' (a driver bug) rethrows immediately, exactly as
%      invz_task2_run_config's own policy does for its own call chain.
%
% summary  struct('n_total_enumerated','n_requested','n_run','n_skipped','n_failed',
%          'n_checkpointed','results_path').
if nargin < 2, opts = struct(); end
results_path = getf(opts, 'results_path', fullfile('.superpowers', 'sdd', 'task2_matrix_results.mat'));
cell_filter  = getf(opts, 'cell_filter', {});
ion          = getf(opts, 'ion', invz_ion());
verbose      = getf(opts, 'verbose', true);

results = load_results(results_path);

if isempty(cell_filter)
    want_ids = {cells.cfg_id};
else
    want_ids = cell_filter;
end

n_run = 0;  n_skip = 0;  n_failed = 0;
for i = 1:numel(want_ids)
    [results, n_run, n_skip, n_failed] = ensure_cell_done(want_ids{i}, cells, results, ...
        results_path, ion, verbose, n_run, n_skip, n_failed);
end

summary = struct('n_total_enumerated', numel(cells), 'n_requested', numel(want_ids), ...
    'n_run', n_run, 'n_skipped', n_skip, 'n_failed', n_failed, ...
    'n_checkpointed', numel(results), 'results_path', results_path);
end

% =================================================================================================
function results = load_results(results_path)
results = empty_results();
if exist(results_path, 'file')
    S = load(results_path, 'results');
    if isfield(S, 'results') && ~isempty(S.results)
        results = S.results;
    end
end
end

% =================================================================================================
function results = empty_results()
results = struct('cfg_id', {}, 'group', {}, 'variant', {}, 'field_T', {}, 'cfg', {}, ...
    'out', {}, 'lattice_provenance', {});
end

% =================================================================================================
function [results, n_run, n_skip, n_failed] = ensure_cell_done(cfg_id, cells, results, ...
    results_path, ion, verbose, n_run, n_skip, n_failed)
%ENSURE_CELL_DONE Recursive, resumable single-cell executor. Skips (bumping n_skip) when
% cfg_id is already present in `results` (loaded from disk at the top of the run, and
% accumulating in-memory as this recurses/iterates) -- the resumability contract. Otherwise
% resolves its prerequisite FIRST (recursion on c.depends_on -- this is what lets a G1/G2
% isolated cell's h_list/seed2 be derived from its swept sibling regardless of processing
% order or which subset opts.cell_filter asked for: the dependency is pulled in and
% run/skipped even when not itself in the filter), builds its final cfg (isolated dependency
% + coupling-variant resolution, invz_task2_resolve_cell_cfg.m), runs it (or synthesizes a
% failed-dependency record if resolution itself failed), appends + immediately re-saves the
% ENTIRE results array to results_path (the incremental-checkpoint contract: a crash/kill
% loses at most this one in-flight cell), and returns the updated bookkeeping counters.
if any(strcmp({results.cfg_id}, cfg_id))
    n_skip = n_skip + 1;
    return;
end
idx = find(strcmp({cells.cfg_id}, cfg_id), 1);
if isempty(idx)
    error('invz:task2Matrix', 'invz_task2_matrix_run: unknown cfg_id ''%s'' -- not in the enumerated cell list.', cfg_id);
end
c = cells(idx);

if ~isempty(c.depends_on)
    [results, n_run, n_skip, n_failed] = ensure_cell_done(c.depends_on, cells, results, ...
        results_path, ion, verbose, n_run, n_skip, n_failed);
end

[cfg_resolved, provenance, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion);
if ~isempty(dep_err)
    out = failed_dependency_record(c, dep_err);
else
    out = invz_task2_run_config(cfg_resolved);
end

rec = struct('cfg_id', c.cfg_id, 'group', c.group, 'variant', c.variant, ...
    'field_T', c.field_T, 'cfg', cfg_resolved, 'out', out, 'lattice_provenance', provenance);
results(end+1) = rec; %#ok<AGROW>
n_run = n_run + 1;
if isfield(out, 'status') && strcmp(out.status, 'failed')
    n_failed = n_failed + 1;
end

save(results_path, 'results');
if verbose
    fprintf('[invz_task2_matrix] ran %s (group=%s, variant=%s, field=%.6g T)\n', ...
        c.cfg_id, c.group, c.variant, c.field_T);
end
end

% =================================================================================================
function out = failed_dependency_record(c, dep_err)
%FAILED_DEPENDENCY_RECORD A structured, invz_task2_run_config-shaped 'failed' record for a
% cell that could not even be ATTEMPTED because its own prerequisite/coupling-materialization
% failed (see invz_task2_resolve_cell_cfg.m) -- mirrors invz_task2_run_config.m's OWN
% out.status='failed' shape field-for-field (out.meta/out.nodes/out.summary/out.extra/
% out.status/out.err_id/out.err_msg) so downstream consumers (2b-report) can treat every
% checkpointed cell uniformly regardless of WHICH layer the failure occurred at.
Bx = c.cfg.Bx;
try
    Bx = invz_field_vec(Bx);
catch
    % leave Bx as given if it is not itself a valid field spec (defensive only)
end
out = struct();
out.meta = struct('T', c.cfg.T, 'Bx', Bx, 'J0eff', NaN, 'Jxx0', NaN, 'mode', c.cfg.mode, ...
    'label', c.cfg_id, 'grid', [], 'dpRng', NaN, 'is_synthetic', NaN, ...
    'solve_opts', getf(c.cfg, 'solve_opts', struct()), 'Jnu_hash', '');
out.nodes = struct([]);
out.summary = struct('n_nodes', 0, 'hstar', NaN, 'hmf_status', 'error', ...
    'replay_mismatch_any', false, 'mode', c.cfg.mode, 'n_stable', 0, 'n_marginal', 0, ...
    'n_unstable', 0, 'n_unconverged', 0, 'n_accepted', 0);
out.extra = struct('hstar', NaN, 'hmf_status', 'error', 'replay_mismatch_any', false, 'trc', []);
out.status  = 'failed';
out.err_id  = 'invz:task2MatrixDependency';
out.err_msg = dep_err;
end
