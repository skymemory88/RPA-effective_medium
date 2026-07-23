function summary = invz_task2_matrix(opts)
%INVZ_TASK2_MATRIX Top-level resumable Task-2 causal-discriminator matrix driver (stage-2c
% task 2b-driver; spec: .superpowers/sdd/task2b-driver-brief.md; frozen definitions:
% docs/invzp_task2_prereg.md Sec. A-H; global constraints: .superpowers/sdd/stage2c-context.md).
%
% Enumerates EVERY prereg Sec. G experiment-matrix cell (G1/G2 isolated-vs-swept/cold-vs-
% continued/multistart, G3 physical downsampling, G4 cardinality-controlled synthetic, G5
% histogram-matched synthetic, G6 grid-offset x dpRng ladder -- see
% invz_task2_matrix_enumerate.m) at the Sec. E frozen physical fields, runs each cell through
% the Task-2a per-config harness (invz_task2_run_config.m), and CHECKPOINTS every result
% incrementally to a resumable .mat file (invz_task2_matrix_run.m). THIS DRIVER ONLY RUNS +
% CHECKPOINTS -- it does not classify beyond what the harness already records per node
% (invz_task2_classify.m's stable/marginal/unstable/unconverged, D_uni/Dq_min/s, full state,
% seed provenance, closest +/-Dq q/branch): the Sec. A/C/D/F analysis and the committed report
% are the NEXT stage (2b-report)'s job.
%
% ---------------------------------------------------------------------------------------
% WHY three sibling helper layers exist below this one function (read before changing
% anything -- this is the part of the task with the most engineering judgement in it):
%
% invz_task2_run_config.m (Task-2a, NOT modified here) has exactly ONE hook for a non-
% default-real-lattice coupling set: cfg.couplings.Jnu_flat, which its resolve_couplings
% ALWAYS treats as fully synthetic (is_synthetic=true, qc/Jnu_unflat forced empty, NO override
% hook). That is fine for a genuinely synthetic fixture, but four of the six Sec. G groups
% here (G3 downsampling, G6 grid offsets) are REAL-lattice-derived and need q/branch
% provenance preserved (prereg Sec. D: "q/branch provenance is retained for the closest modes
% at every rung"; brief G3: "retain the flat->(q,branch) map so Task-0 provenance still
% resolves"). Modifying invz_task2_run_config.m to add such a hook is explicitly out of scope
% (task2b-driver-brief.md's scope fence: "Do NOT modify ... the Task-2a harness"). So this
% driver keeps the needed provenance ITSELF, alongside each checkpointed cell
% (results(k).lattice_provenance: qc/Jnu_unflat/nq, row-aligned with whatever Jnu_flat was fed
% through the synthetic-injection hook) rather than inside invz_task2_run_config's own
% out.nodes(:).qbranch_pos/neg (which will correctly, harmlessly read back as the NaN triple
% its own is_synthetic=true branch always produces for these cells -- a known, documented
% consequence of not touching that file, not a silent loss: 2b-report recovers q/branch
% identity via invz_ordered_trace_resolve(struct('is_synthetic',false,'Jnu_unflat',
% results(k).lattice_provenance.Jnu_unflat,'nq',results(k).lattice_provenance.nq),
% out.nodes(j).idx_pos_flat), then results(k).lattice_provenance.qc(q_idx,:) for the physical
% (h,k,l) -- see invz_task2_resolve_cell_cfg.m's header for the full construction, and
% invz_ordered_trace_resolve.m (Task 0, unmodified, called read-only here) for the resolver
% itself.
%
% Similarly, invz_task2_run_config's mode='isolated' takes cfg.h_list as a REQUIRED, already-
% concrete input -- there is no hook to say "solve whatever h a companion swept cell visited."
% Since prereg Sec. G item 1 / invz_task2_run_config's OWN run_isolated_mode docstring frame
% "isolated" as literally re-solving "the SAME h a swept run visited", the G1/G2 isolated
% cells' h_list (and, for the multistart variant, its seed) is deferred at enumeration time
% and filled in from the field's own swept sibling's CHECKPOINTED output immediately before
% the isolated cell is run (invz_task2_resolve_cell_cfg.m) -- a cross-cell dependency
% (cells(:).depends_on) that invz_task2_matrix_run.m resolves recursively and resumably (a
% filtered subset that asks for an isolated cell without its swept prerequisite still pulls
% the prerequisite in automatically; already-checkpointed prerequisites are simply skipped).
%
% ★ THE DELICATE PART (G6 grid offset / dipolar-Gamma handling): see
% invz_task2_couplings_shifted_grid.m's own header for the full derivation. Summary: a naive
% "shift the grid, then apply the OLD hardcoded q==0 Gamma filter" would silently fail to
% exclude a Gamma-equivalent point the shift reintroduces (proven, not assumed: the standard
% unshifted 16^3 grid never actually contains q=(0,0,0) -- 15 intervals is odd, so no linspace
% sample lands on exactly zero -- but a uniform half-of-the-axis-step shift provably DOES
% reintroduce exactly (0,0,0), verified numerically over the full 3D grid: 0 Gamma-equivalent
% points unshifted, exactly 1 after the half-step shift). This driver's helper avoids the bug
% by shifting FIRST and filtering with the GENERAL invz_is_gamma_equiv test SECOND (the same
% test invz_jq_modes.m itself uses for Lorentz-cavity placement), on the shifted coordinates,
% never the pre-shift ones -- see test_invz_task2_matrix.m's
% test_shifted_grid_half_step_drops_exactly_gamma and
% test_shifted_grid_naive_preshift_filter_would_miss_gamma for the concrete, numeric proof.
%
% opts.ion           (default invz_ion())
% opts.results_path  (default '.superpowers/sdd/task2_matrix_results.mat', git-ignored scratch
%                    -- see .gitignore's top-level '.superpowers/' entry). Callers that want to
%                    run a validation subset WITHOUT touching the real controller checkpoint
%                    MUST override this (test_invz_task2_matrix.m always does).
% opts.cell_filter   (default {} = the full 40-cell matrix; else a cellstr of cfg_ids -- see
%                    invz_task2_matrix_run.m).
% opts.verbose       (default true)
% opts.dry_run       (default false) build run_opts + enumerate, then RETURN WITHOUT running
%                    any cell. Exists so the DEFAULT (unfiltered, full-matrix) option path --
%                    the one the controller actually runs -- is reachable by a test without
%                    paying for 40 real solves; see the run_opts note below for the bug this
%                    guards. summary = struct('dry_run',true,'n_total_enumerated',..,'run_opts',..).
%
% summary  invz_task2_matrix_run.m's own return (n_total_enumerated/n_requested/n_run/
%          n_skipped/n_failed/n_checkpointed/results_path).
%
% NOTE (binding, stage-2c task 2b-driver scope fence): this function is capable of running
% the FULL matrix (opts.cell_filter = {}), which is expensive (many real 16^3-grid solves)
% and is explicitly NOT this task's job -- "the CONTROLLER runs the full (expensive) matrix as
% a background job" (task2b-driver-brief.md). This task's own validation
% (test_invz_task2_matrix.m) NEVER calls this top-level function unfiltered; it always passes
% a tiny, explicit opts.cell_filter (or calls invz_task2_matrix_enumerate/invz_task2_matrix_run
% directly) against a temporary opts.results_path, run in the foreground.
if nargin < 1, opts = struct(); end
ion = getf(opts, 'ion', invz_ion());
cells = invz_task2_matrix_enumerate(ion);
% run_opts is built FIELD-BY-FIELD, never via struct('cell_filter', <cell>, ...). MATLAB's
% struct() CONSTRUCTOR treats a cell value as an array generator, so struct('cell_filter', {})
% yields a 0x0 EMPTY STRUCT ARRAY -- not a scalar struct whose cell_filter is {} -- and every
% downstream getf(run_opts, ...) then dies with "Insufficient number of outputs". That is
% exactly the DEFAULT (unfiltered, full-matrix) path this driver exists to run, and no test
% reached it (the suite only ever called invz_task2_matrix_enumerate/_run directly, always with
% an explicit non-empty filter). Plain field assignment has no such gotcha: s.f = {} is fine.
run_opts = struct();
run_opts.results_path = getf(opts, 'results_path', fullfile('.superpowers', 'sdd', 'task2_matrix_results.mat'));
run_opts.cell_filter  = getf(opts, 'cell_filter', {});
run_opts.ion          = ion;
run_opts.verbose      = getf(opts, 'verbose', true);
if getf(opts, 'dry_run', false)
    summary = struct('dry_run', true, 'n_total_enumerated', numel(cells), 'run_opts', run_opts);
    return;
end
summary = invz_task2_matrix_run(cells, run_opts);
end
