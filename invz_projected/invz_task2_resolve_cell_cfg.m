function [cfg, provenance, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion)
%INVZ_TASK2_RESOLVE_CELL_CFG Fully resolves one enumerated cell's DEFERRED cfg fields into the
% concrete cfg invz_task2_run_config.m already understands natively, immediately before that
% cell is run (stage-2c task 2b-driver). Two independent kinds of deferral, both resolved here:
%
% (1) G1/G2 isolated-companion dependency (c.depends_on non-empty): cfg.h_list (and, for the
%     multistart variant, cfg.seed_list{1}) is filled in from the swept sibling's OWN
%     checkpointed `results` entry -- see invz_task2_derive_isolated_h_list.m. The multistart
%     seed (c.resolve.seed2_from_dep) is a deterministic, non-cold, non-random SECOND starting
%     guess: Sigma = zeros(Nw,1) (unperturbed -- isolates K0 as the one changed variable),
%     K0s = J0eff (the sibling's OWN J0eff, i.e. "start the K0 EMT self-consistency AT the bare
%     coupling value" -- the natural alternate anchor to K0=0, since D_uni=1+(J0eff-K0)*Gstat's
%     own (J0eff-K0) term is exactly zero there). This is an implementer's judgement call
%     (flagged for review, matching Task-2a's own precedent of flagging non-frozen choices):
%     invz_ordered_node_solve.m's OWN cold/warm dichotomy has no third "second cold" state to
%     borrow, so ANY genuinely different deterministic multistart probe must be warm-LABELLED
%     (info.seed_kind='warm', or 'cold_after_warm_fail' if it does not converge and the
%     solver's own cold_retry falls back -- both distinguishable from the plain 'isolated_cold'
%     variant's seed_kind='cold' in the checkpointed record either way).
%
% (2) Coupling-variant materialization (cfg.couplings.variant present): a declarative recipe
%     (G3 'downsample' / G4 'cardinality_synth' / G5 'histmatch_synth' / G6 'shifted_grid',
%     produced by invz_task2_matrix_enumerate.m) is turned into the concrete
%     struct('Jnu_flat',...,'J0eff',...,'Jxx0',...) invz_task2_run_config.m's resolve_couplings
%     already accepts as its ONE non-default-real-lattice input hook. provenance (non-empty
%     for the REAL-lattice-derived variants 'downsample'/'shifted_grid') carries qc/Jnu_unflat/
%     nq -- the q/branch identity invz_task2_run_config's own internal resolve_qbranch cannot
%     supply for a synthetic-injected Jnu_flat (its resolve_couplings hard-codes
%     is_synthetic=true whenever cfg.couplings.Jnu_flat is supplied, with no override hook --
%     see invz_task2_matrix.m's header for why this driver keeps provenance itself instead of
%     modifying that harness). 2b-report recovers q/branch identity for these cells' closest
%     +/-Dq modes via invz_ordered_trace_resolve(struct('is_synthetic',false,'Jnu_unflat',
%     provenance.Jnu_unflat,'nq',provenance.nq), out.nodes(k).idx_pos_flat), then
%     provenance.qc(q_idx,:) for the physical (h,k,l).
%
% Pure w.r.t. `results` (never mutated); dep_err is '' on success, else a human-readable
% reason resolution failed -- the caller (invz_task2_matrix_run.m) checkpoints a structured
% failed-dependency record instead of calling invz_task2_run_config, and does NOT abort the
% driver loop (mirrors the repo's invz:* error-absorption convention, applied at the driver
% level for this specific -- and, for the real matrix, expected to be rare/never-hit --
% failure mode: a prerequisite cell that itself errored or produced no nodes).
if nargin < 3 || isempty(ion), ion = invz_ion(); end
cfg = c.cfg;
provenance = struct();
dep_err = '';

if ~isempty(c.depends_on)
    di = find(strcmp({results.cfg_id}, c.depends_on), 1);
    if isempty(di)
        dep_err = sprintf(['prerequisite ''%s'' not found in results (should be impossible -- ' ...
            'invz_task2_matrix_run''s ensure_cell_done must run/skip it before its dependents).'], ...
            c.depends_on);
        return;
    end
    dep_out = results(di).out;
    if isfield(dep_out, 'status') && strcmp(dep_out.status, 'failed')
        dep_err = sprintf('prerequisite ''%s'' itself failed (%s: %s).', c.depends_on, ...
            dep_out.err_id, dep_out.err_msg);
        return;
    end
    if isempty(dep_out.nodes)
        dep_err = sprintf('prerequisite ''%s'' produced no nodes (cannot derive an isolated h_list).', ...
            c.depends_on);
        return;
    end
    try
        cfg.h_list = invz_task2_derive_isolated_h_list(dep_out);
    catch ME
        dep_err = sprintf('deriving h_list from prerequisite ''%s'' failed: %s', c.depends_on, ME.message);
        return;
    end
    if isfield(c, 'resolve') && isfield(c.resolve, 'seed2_from_dep') && c.resolve.seed2_from_dep
        J0eff_dep = dep_out.meta.J0eff;
        Ecut = getf(getf(cfg, 'solve_opts', struct()), 'Ecut', 40);
        [wn, ~, ~] = invz_matsubara(cfg.T, Ecut);
        cfg.seed_list = {struct('Sigma', zeros(numel(wn), 1), 'K0s', J0eff_dep)};
    end
end

cpl = getf(cfg, 'couplings', struct());
if isfield(cpl, 'variant') && ~isempty(cpl.variant)
    try
        [cfg, provenance] = local_materialize(cfg, cpl, ion);
    catch ME
        if startsWith(ME.identifier, 'invz:') && ~startsWith(ME.identifier, 'invz:task2')
            dep_err = sprintf('coupling materialization (variant ''%s'') raised %s: %s', ...
                cpl.variant, ME.identifier, ME.message);
        else
            rethrow(ME);
        end
    end
end
end

% =================================================================================================
function [cfg, provenance] = local_materialize(cfg, cpl, ion)
provenance = struct();
switch cpl.variant
    case 'downsample'
        [~, J0eff_full, Jxx0_full, qc_full, Jnu_unflat_full, ~, ~] = ...
            invz_task2_couplings_shifted_grid(ion, cpl.grid, cpl.dpRng, 'unshifted');
        [Jnu_flat, qc_sub, Jnu_unflat_sub, keep_idx] = invz_task2_couplings_downsample( ...
            Jnu_unflat_full, qc_full, cpl.stride);
        cfg.couplings = struct('Jnu_flat', Jnu_flat, 'J0eff', J0eff_full, 'Jxx0', Jxx0_full);
        provenance = struct('qc', qc_sub, 'Jnu_unflat', Jnu_unflat_sub, ...
            'nq', size(Jnu_unflat_sub, 1), 'keep_idx', keep_idx, 'construction', sprintf( ...
            'stride-%d subsample of the %s/dpRng=%d baseline q-grid (flat q-row order, all 4 branches kept per row).', ...
            cpl.stride, mat2str(cpl.grid), cpl.dpRng));

    case 'shifted_grid'
        [Jnu_flat, J0eff, Jxx0, qc, Jnu_unflat, ~, n_gamma_dropped] = ...
            invz_task2_couplings_shifted_grid(ion, cpl.grid, cpl.dpRng, cpl.shift_mode);
        cfg.couplings = struct('Jnu_flat', Jnu_flat, 'J0eff', J0eff, 'Jxx0', Jxx0);
        provenance = struct('qc', qc, 'Jnu_unflat', Jnu_unflat, 'nq', size(Jnu_unflat, 1), ...
            'n_gamma_dropped', n_gamma_dropped, 'construction', sprintf( ...
            '%s q-grid (%s, dpRng=%d); %d Gamma-equivalent point(s) excluded via invz_is_gamma_equiv.', ...
            cpl.shift_mode, mat2str(cpl.grid), cpl.dpRng, n_gamma_dropped));

    case 'cardinality_synth'
        base_vals = getf(cpl, 'base_vals', linspace(-2e-3, 6.0e-3, 24).');
        J0eff_ = getf(cpl, 'J0eff', 6.42e-3);
        [Jnu_flat, detail] = invz_task2_couplings_cardinality_synth(base_vals, cpl.n_total);
        cfg.couplings = struct('Jnu_flat', Jnu_flat, 'J0eff', J0eff_);
        provenance = struct('construction', sprintf(['cardinality-duplicated %d-pt fixture ' ...
            '(%d full tiles + %d padding entries = %d total, NO RNG).'], detail.n_base, ...
            detail.n_rep, detail.n_rem, cpl.n_total));

    case 'histmatch_synth'
        [Jnu_ref_flat, info_ref, Jaa0_ref] = invz_bz_couplings(ion, ...
            struct('grid', cpl.grid, 'dpRng', cpl.dpRng));
        Jnu_flat = invz_task2_couplings_histmatch_synth(Jnu_ref_flat, cpl.n_total);
        cfg.couplings = struct('Jnu_flat', Jnu_flat, 'J0eff', info_ref.Jcc0, 'Jxx0', Jaa0_ref);
        provenance = struct('construction', sprintf(['inverse-CDF histogram match to the real ' ...
            '%s/dpRng=%d Jnu distribution at matched cardinality %d (fixed quantile grid, NO RNG).'], ...
            mat2str(cpl.grid), cpl.dpRng, cpl.n_total));

    otherwise
        error('invz:task2Config', 'invz_task2_resolve_cell_cfg: unknown couplings.variant ''%s''.', cpl.variant);
end
end
