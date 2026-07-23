function tests = test_invz_task2_matrix
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% ===========================================================================================
% (1) ENUMERATION: full-matrix cell counts, group tags, uniqueness, determinism (stage-2c
% task 2b-driver gate: "Enumeration: assert the full cell list has the expected group counts
% and every cfg_id is unique and deterministic (calling the enumerator twice yields identical
% ids)"). PURE/no solves -- instantaneous.
% ===========================================================================================
function test_enumerate_cell_counts_and_group_tags(testCase)
cells = invz_task2_matrix_enumerate();
verifyEqual(testCase, numel(cells), 40, 'expected 40 total cells (12+12+2+2+12, brief tally)');

groups = {cells.group};
verifyEqual(testCase, nnz(strcmp(groups, 'G1G2')), 12);
verifyEqual(testCase, nnz(strcmp(groups, 'G3')),   12);
verifyEqual(testCase, nnz(strcmp(groups, 'G4')),    2);
verifyEqual(testCase, nnz(strcmp(groups, 'G5')),    2);
verifyEqual(testCase, nnz(strcmp(groups, 'G6')),   12);

% Every field_T is one of the 4 frozen prereg Sec. E values (Task-2a report table).
FIELDS = [1.173192, 2.581023, 3.754215, 2.85];
for k = 1:numel(cells)
    verifyTrue(testCase, any(abs(cells(k).field_T - FIELDS) < 1e-9), ...
        sprintf('cell %s has an unrecognized field_T=%.6g', cells(k).cfg_id, cells(k).field_T));
end
end

function test_enumerate_cfg_id_unique_and_stable_format(testCase)
cells = invz_task2_matrix_enumerate();
ids = {cells.cfg_id};
verifyEqual(testCase, numel(unique(ids)), numel(ids), 'every cfg_id must be unique');
for k = 1:numel(ids)
    verifyTrue(testCase, ischar(ids{k}) || isstring(ids{k}));
    % no whitespace/punctuation that would be unsafe as a MAT variable/filename component
    % (the field tag's leading 'F', e.g. 'F1173192', is the only uppercase letter used).
    verifyTrue(testCase, ~isempty(regexp(ids{k}, '^[a-zA-Z0-9_]+$', 'once')), ...
        sprintf('cfg_id ''%s'' is not a plain alphanumeric/underscore string', ids{k}));
end
end

function test_enumerate_deterministic(testCase)
% No Date/now/rand/tempname anywhere in the enumerator: calling it twice must reproduce the
% IDENTICAL cell list, in the SAME order, byte-for-byte on every field.
c1 = invz_task2_matrix_enumerate();
c2 = invz_task2_matrix_enumerate();
verifyEqual(testCase, {c1.cfg_id}, {c2.cfg_id});
verifyEqual(testCase, {c1.group}, {c2.group});
verifyEqual(testCase, {c1.variant}, {c2.variant});
verifyEqual(testCase, [c1.field_T], [c2.field_T]);
verifyEqual(testCase, {c1.depends_on}, {c2.depends_on});
for k = 1:numel(c1)
    verifyEqual(testCase, c1(k).cfg, c2(k).cfg, sprintf('cell %d (%s) cfg differs across calls', k, c1(k).cfg_id));
    verifyEqual(testCase, c1(k).resolve, c2(k).resolve);
end
end

function test_enumerate_isolated_cells_depend_on_field_matched_swept_cell(testCase)
cells = invz_task2_matrix_enumerate();
is_iso = strcmp({cells.variant}, 'isolated_cold') | strcmp({cells.variant}, 'isolated_seed2');
verifyEqual(testCase, nnz(is_iso), 8);   % 2 isolated variants x 4 fields
for k = find(is_iso)
    verifyFalse(testCase, isempty(cells(k).depends_on), ...
        sprintf('%s must declare a swept-cell dependency', cells(k).cfg_id));
    dep_idx = find(strcmp({cells.cfg_id}, cells(k).depends_on));
    verifyEqual(testCase, numel(dep_idx), 1);
    verifyEqual(testCase, cells(dep_idx).variant, 'swept');
    verifyEqual(testCase, cells(dep_idx).field_T, cells(k).field_T, ...
        'an isolated cell must depend on the SAME field''s own swept cell');
end
is_seed2 = strcmp({cells.variant}, 'isolated_seed2');
for k = find(is_seed2)
    verifyTrue(testCase, cells(k).resolve.seed2_from_dep);
end
is_cold = strcmp({cells.variant}, 'isolated_cold');
for k = find(is_cold)
    verifyFalse(testCase, cells(k).resolve.seed2_from_dep);
end
% self-contained groups (G3/G4/G5/G6) never depend on another cell.
is_self_contained = ismember({cells.group}, {'G3', 'G4', 'G5', 'G6'});
for k = find(is_self_contained)
    verifyTrue(testCase, isempty(cells(k).depends_on));
end
end

% ===========================================================================================
% (2) Bc_PM RE-DERIVATION (brief: "Re-derive Bc_PM via the same independent PM-mass path 2a
% used and assert it matches 4.692769, or use the frozen values citing 2a -- state which").
% invz_task2_matrix_enumerate.m cites the frozen values (fast/pure enumeration); THIS test
% independently re-derives Bc_PM via the SAME invz_critical PM-mass root-find Task-2a used and
% confirms both (a) it still matches 4.692769 and (b) the enumerator's own hardcoded field
% values agree with {0.25,0.55,0.80}*Bc_PM to the frozen table's own precision. Not
% INVZ_SLOW-gated: invz_critical's bracketing scan on the real lattice is ~6 s (Task-2a
% report), comfortably inside the "fast suite" norm (well under the brief's own ~1-2 min
% INVZ_SLOW threshold).
% ===========================================================================================
function test_bc_pm_rederivation_matches_frozen_fields(testCase)
ion = invz_ion();
[Jf, info, Jaa0] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30));
Bc_PM = invz_critical(ion, 0.1, Jf, struct('J0eff', info.Jcc0, 'Jxx0', Jaa0));
verifyEqual(testCase, Bc_PM, 4.692769, 'AbsTol', 1e-3);

frozen = [1.173192, 2.581023, 3.754215];
derived = [0.25, 0.55, 0.80] * Bc_PM;
verifyEqual(testCase, derived, frozen, 'AbsTol', 1e-3);
end

% ===========================================================================================
% (3) G6 SHIFTED-GRID HELPER -- the DELICATE dipolar-Gamma handling. Cheap: only coupling
% GENERATION (cached invz_jq_modes calls on the real 16^3 grid), no iterative solve at all.
% ===========================================================================================
function test_shifted_grid_unshifted_matches_bz_couplings(testCase)
% Regression proof that the reimplemented grid-generation recipe is bit-faithful to
% production: 'unshifted' must reproduce invz_bz_couplings/invz_jq_modes EXACTLY (same qc,
% hence same invz_jq_modes cache key, hence bit-identical outputs).
ion = invz_ion();
[Jnu_ref, info_ref, Jaa0_ref] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30));

[Jnu_flat, J0eff, Jxx0, qc, Jnu_unflat, ~, n_gamma_dropped] = ...
    invz_task2_couplings_shifted_grid(ion, [16 16 16], 30, 'unshifted');

verifyEqual(testCase, n_gamma_dropped, 0, ...
    'the standard 16^3 grid never contains an exact Gamma-equivalent point (odd interval count)');
verifyEqual(testCase, size(qc, 1), 4096);
verifyEqual(testCase, size(Jnu_unflat), [4096 4]);
verifyEqual(testCase, Jnu_flat, Jnu_ref, 'AbsTol', 0);
verifyEqual(testCase, J0eff, info_ref.Jcc0, 'AbsTol', 0);
verifyEqual(testCase, Jxx0, Jaa0_ref, 'AbsTol', 0);
end

function test_shifted_grid_half_step_drops_exactly_gamma(testCase)
% THE concrete, numeric proof the delicate Gamma handling actually fires (non-vacuously) on
% the half-step ladder rung: a uniform shift by half of the grid's own axis step reintroduces
% EXACTLY one Gamma-equivalent point -- (0,0,0) itself -- which the general
% invz_is_gamma_equiv-based filter must (and does) exclude.
ion = invz_ion();
[Jnu_flat, J0eff, Jxx0, qc, Jnu_unflat, latinfo, n_gamma_dropped] = ...
    invz_task2_couplings_shifted_grid(ion, [16 16 16], 30, 'half_step');

verifyEqual(testCase, n_gamma_dropped, 1);
verifyEqual(testCase, size(qc, 1), 4095);
verifyEqual(testCase, size(Jnu_unflat), [4095 4]);
verifyEqual(testCase, numel(Jnu_flat), 4095*4);
verifyEqual(testCase, J0eff, latinfo.Jcc0);
verifyEqual(testCase, Jxx0, latinfo.Jaa0);

% Minor #2 review-fix: a DIRECT self-check on the FUNCTION'S OWN returned qc (not merely
% inferred from n_gamma_dropped==1 and the row-count dropping by exactly one) -- the row
% actually excluded is self-evidently near-zero, on the qc this function actually hands back.
verifyFalse(testCase, any(all(abs(qc) < 1e-9, 2)), ...
    'the returned qc must contain no near-zero (Gamma-equivalent) row');

% independently reconstruct the PRE-filter shifted grid and confirm the dropped point really
% is (0,0,0) (not some other, unrelated point) -- proves the filter fired for the RIGHT reason.
qc0 = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
step = 1/15;   % linspace(-0.5,0.5,16) has 15 intervals
qc_shifted = qc0 + step/2;
is_zero = all(abs(qc_shifted) < 1e-9, 2);
verifyEqual(testCase, nnz(is_zero), 1);
end

function test_shifted_grid_naive_preshift_filter_would_miss_gamma(testCase)
% Concrete demonstration of the BUG this task's implementation avoids: filtering with the OLD
% hardcoded "q==0" check BEFORE shifting is a no-op (the unshifted grid never hits exactly
% zero), so the shift silently reintroduces an UNEXCLUDED Gamma point -- unlike this driver's
% own helper, which shifts FIRST and filters (generally, via invz_is_gamma_equiv) SECOND.
ion = invz_ion();
qc0 = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
step = 1/15;

naive_filtered_then_shifted = qc0(any(abs(qc0) > 1e-12, 2), :) + step/2;   % BUGGY order
verifyEqual(testCase, size(naive_filtered_then_shifted, 1), 4096, ...
    'the naive pre-shift filter must be a no-op here (nothing to drop before shifting)');

[~, ~, ~, qc_correct] = invz_task2_couplings_shifted_grid(ion, [16 16 16], 30, 'half_step');
verifyEqual(testCase, size(qc_correct, 1), 4095, ...
    'the correct (shift-then-filter, general Gamma test) construction must exclude the reintroduced point');
end

function test_shifted_grid_half_step_rejects_noncubic_grid(testCase)
% Minor #4 review-fix: 'half_step' builds ONE scalar delta = step/2 from grid(1) alone and adds
% it uniformly to h/k/l, which silently assumes a cubic grid. A non-cubic grid must be rejected
% explicitly (rather than silently mis-shifting two of the three axes by the wrong step) -- this
% task only ever uses [16 16 16], so a guard (not a per-axis delta) is the appropriate fix.
ion = invz_ion();
verifyError(testCase, ...
    @() invz_task2_couplings_shifted_grid(ion, [16 16 8], 30, 'half_step'), ...
    'invz:task2Config');
end

function test_resolve_cell_cfg_materializes_shifted_grid_variant(testCase)
ion = invz_ion();
c = struct('cfg', struct('ion', ion, 'T', 0.1, 'Bx', 2.85, 'mode', 'swept', ...
        'couplings', struct('variant', 'shifted_grid', 'shift_mode', 'half_step', ...
                             'grid', [16 16 16], 'dpRng', 30), 'label', 'x'), ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false));
results = struct('cfg_id', {}, 'group', {}, 'variant', {}, 'field_T', {}, 'cfg', {}, ...
    'out', {}, 'lattice_provenance', {});
[cfg_out, provenance, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion);
verifyTrue(testCase, isempty(dep_err));
verifyTrue(testCase, isfield(cfg_out.couplings, 'Jnu_flat'));
verifyEqual(testCase, numel(cfg_out.couplings.Jnu_flat), 4095*4);
verifyEqual(testCase, provenance.n_gamma_dropped, 1);
verifyEqual(testCase, provenance.nq, 4095);
verifyEqual(testCase, size(provenance.qc, 1), 4095);
end

% ===========================================================================================
% (4) G3 DOWNSAMPLING HELPER -- deterministic, provenance-preserving. Cheap (array slicing on
% an already-cached real lattice).
% ===========================================================================================
function test_downsample_preserves_rows_and_determinism(testCase)
% Review-fix: also asserts the anti-single-axis-collapse property directly -- G3's whole point
% (prereg Sec. G item 3: "isolates coupling density from distribution"). Every one of the three
% BZ axes must retain its FULL 16-value marginal range after decimation, at every stride: this
% is exactly what the earlier flat-stride construction violated (it collapsed the 2nd axis to
% 8/4/2 values at stride 2/4/8 while leaving the other two axes untouched at all 16).
ion = invz_ion();
grid0 = [16 16 16];
[~, ~, ~, qc_full, Jnu_unflat_full] = invz_task2_couplings_shifted_grid(ion, grid0, 30, 'unshifted');

for stride = [2 4 8]
    [Jnu_flat, qc_sub, Jnu_unflat_sub, keep_idx] = ...
        invz_task2_couplings_downsample(Jnu_unflat_full, qc_full, grid0, stride);
    n_expected = 4096/stride;
    verifyEqual(testCase, numel(keep_idx), n_expected, ...
        sprintf('stride %d must keep exactly nq_full/stride points', stride));
    verifyEqual(testCase, size(Jnu_unflat_sub, 1), n_expected);
    verifyEqual(testCase, Jnu_unflat_sub, Jnu_unflat_full(keep_idx, :), 'AbsTol', 0);
    verifyEqual(testCase, qc_sub, qc_full(keep_idx, :), 'AbsTol', 0);
    verifyEqual(testCase, Jnu_flat, Jnu_unflat_sub(:), 'AbsTol', 0);
    verifyEqual(testCase, keep_idx(1), 1);   % first q-row (h=k=l=0) always kept: sum==0 == 0 mod any stride

    % THE anti-single-axis-collapse assertion -- what the original flat-stride bug would have
    % failed: all three axes retain every one of their 16 unique reduced-coordinate values.
    for axis = 1:3
        verifyEqual(testCase, numel(unique(qc_sub(:, axis))), 16, ...
            sprintf('stride %d: axis %d lost marginal range (single-axis-collapse regression)', stride, axis));
    end

    % determinism: calling again reproduces the identical result.
    [Jnu_flat2, qc_sub2, Jnu_unflat_sub2, keep_idx2] = ...
        invz_task2_couplings_downsample(Jnu_unflat_full, qc_full, grid0, stride);
    verifyEqual(testCase, Jnu_flat, Jnu_flat2, 'AbsTol', 0);
    verifyEqual(testCase, qc_sub, qc_sub2, 'AbsTol', 0);
    verifyEqual(testCase, keep_idx, keep_idx2, 'AbsTol', 0);
end
end

function test_resolve_cell_cfg_materializes_downsample_variant(testCase)
ion = invz_ion();
c = struct('cfg', struct('ion', ion, 'T', 0.1, 'Bx', 2.85, 'mode', 'swept', ...
        'couplings', struct('variant', 'downsample', 'stride', 8, ...
                             'grid', [16 16 16], 'dpRng', 30), 'label', 'x'), ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false));
results = struct('cfg_id', {}, 'group', {}, 'variant', {}, 'field_T', {}, 'cfg', {}, ...
    'out', {}, 'lattice_provenance', {});
[cfg_out, provenance, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion);
verifyTrue(testCase, isempty(dep_err));
verifyEqual(testCase, numel(cfg_out.couplings.Jnu_flat), 4*512);
verifyEqual(testCase, provenance.nq, 512);
verifyEqual(testCase, size(provenance.qc, 1), 512);
% review-fix: end-to-end confirmation (not just the raw helper's own unit test) that all three
% axes keep their full 16-value marginal range through the FULL resolve_cell_cfg pathway.
for axis = 1:3
    verifyEqual(testCase, numel(unique(provenance.qc(:, axis))), 16);
end
end

function test_resolve_cell_cfg_absorbs_nontask2_materialization_error(testCase)
% Minor #1 review-fix: the coupling-materialization catch (invz_task2_resolve_cell_cfg.m, the
% try/catch around local_materialize) absorbs an 'invz:*' identifier that is NOT 'invz:task2*'
% into dep_err, and RETHROWS everything else -- exactly invz_task2_run_config.m's own
% task2_is_absorbable(id) predicate. NONE of the 4 concrete couplings.variant helpers this
% driver ships today (downsample/shifted_grid/cardinality_synth/histmatch_synth, nor their own
% callees qVec_generator/invz_is_gamma_equiv/invz_jq_modes' default non-ODD path) ever raises
% anything but 'invz:task2Config' or a plain non-'invz' MATLAB error -- i.e. the ABSORB branch
% is not reachable today via any genuine production error from those helpers (every error(...)
% call in that reachable call graph was inspected at implementation time: all are
% 'invz:task2*'). To exercise the ABSORB branch deterministically and cheaply, this test
% shadows qVec_generator -- a stable, non-task2 leaf utility reached via the 'downsample'
% variant's invz_task2_couplings_shifted_grid call, untouched by this review-fix pass -- on the
% MATLAB path with a throwaway stub that raises a genuine physics-class identifier
% (invz:transverseMF, one of the identifiers task2_is_absorbable's own docstring names as a
% "genuine production physics/numerics exception"). The shadow directory is added ahead of the
% real qVec_generator via addpath and removed in onCleanup (guaranteed even on test failure), so
% no other test observes it.
shadow_dir = tempname();
mkdir(shadow_dir);
cleanupObj = onCleanup(@() cleanup_shadow_path(shadow_dir)); %#ok<NASGU>
fid = fopen(fullfile(shadow_dir, 'qVec_generator.m'), 'w');
fprintf(fid, ['function varargout = qVec_generator(varargin) %%#ok<STOUT>\n' ...
    'error(''invz:transverseMF'', ''synthetic non-task2 materialization failure (task2b-driver review-fix absorb-branch test).'');\n' ...
    'end\n']);
fclose(fid);
addpath(shadow_dir, '-begin');
rehash path;

ion = invz_ion();
c = struct('cfg', struct('ion', ion, 'T', 0.1, 'Bx', 2.85, 'mode', 'swept', ...
        'couplings', struct('variant', 'downsample', 'stride', 8, ...
                             'grid', [16 16 16], 'dpRng', 30), 'label', 'x'), ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false));
results = struct('cfg_id', {}, 'group', {}, 'variant', {}, 'field_T', {}, 'cfg', {}, ...
    'out', {}, 'lattice_provenance', {});
[cfg_out, provenance, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion);

verifyFalse(testCase, isempty(dep_err), ...
    'a non-task2 invz:* materialization error must be ABSORBED into dep_err, not crash the driver');
verifyTrue(testCase, contains(dep_err, 'invz:transverseMF'));
verifyTrue(testCase, contains(dep_err, 'downsample'));
verifyEqual(testCase, cfg_out, c.cfg, 'on an absorbed failure, cfg must be returned unmodified');
verifyEqual(testCase, provenance, struct(), 'on an absorbed failure, provenance must stay the empty-struct default');
end

% ===========================================================================================
% (5) G4/G5 SYNTHETIC-COUPLING HELPERS -- deterministic, no RNG. Cheap.
% ===========================================================================================
function test_cardinality_synth_exact_count_and_values_from_base(testCase)
base24 = linspace(-2e-3, 6.0e-3, 24).';
[Jnu_flat, detail] = invz_task2_couplings_cardinality_synth(base24, 16384);
verifyEqual(testCase, numel(Jnu_flat), 16384);
verifyEqual(testCase, detail.n_base, 24);
verifyEqual(testCase, detail.n_rep, 682);
verifyEqual(testCase, detail.n_rem, 16);
verifyEqual(testCase, detail.n_rep*24 + detail.n_rem, 16384);
verifyTrue(testCase, all(ismember(Jnu_flat, base24)), 'every entry must come from base24 (no new values)');

Jnu_flat2 = invz_task2_couplings_cardinality_synth(base24, 16384);
verifyEqual(testCase, Jnu_flat, Jnu_flat2, 'AbsTol', 0);
end

function test_resolve_cell_cfg_materializes_cardinality_synth_variant(testCase)
ion = invz_ion();
base24 = linspace(-2e-3, 6.0e-3, 24).';
c = struct('cfg', struct('ion', ion, 'T', 0.1, 'Bx', 2.85, 'mode', 'swept', ...
        'couplings', struct('variant', 'cardinality_synth', 'n_total', 16384, ...
                             'base_vals', base24, 'J0eff', 6.42e-3), 'label', 'x'), ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false));
results = struct('cfg_id', {}, 'group', {}, 'variant', {}, 'field_T', {}, 'cfg', {}, ...
    'out', {}, 'lattice_provenance', {});
[cfg_out, ~, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion);
verifyTrue(testCase, isempty(dep_err));
verifyEqual(testCase, numel(cfg_out.couplings.Jnu_flat), 16384);
verifyEqual(testCase, cfg_out.couplings.J0eff, 6.42e-3);
end

function test_histmatch_synth_exact_resort_at_equal_cardinality(testCase)
% The degenerate, exactly-checkable property: at n_total == numel(real_ref), the inverse-CDF
% construction (fixed midpoint quantile grid) reproduces sort(real_ref) exactly.
ion = invz_ion();
[Jnu_ref, ~, ~] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30));
out = invz_task2_couplings_histmatch_synth(Jnu_ref, numel(Jnu_ref));
verifyEqual(testCase, sort(out), sort(Jnu_ref), 'AbsTol', 1e-12);

% smaller cardinality: no crash, stays within the reference's observed range, deterministic.
out_half = invz_task2_couplings_histmatch_synth(Jnu_ref, floor(numel(Jnu_ref)/2));
verifyEqual(testCase, numel(out_half), floor(numel(Jnu_ref)/2));
verifyGreaterThanOrEqual(testCase, min(out_half), min(Jnu_ref) - 1e-12);
verifyLessThanOrEqual(testCase, max(out_half), max(Jnu_ref) + 1e-12);
out_half2 = invz_task2_couplings_histmatch_synth(Jnu_ref, floor(numel(Jnu_ref)/2));
verifyEqual(testCase, out_half, out_half2, 'AbsTol', 0);
end

function test_resolve_cell_cfg_materializes_histmatch_synth_variant(testCase)
ion = invz_ion();
c = struct('cfg', struct('ion', ion, 'T', 0.1, 'Bx', 2.85, 'mode', 'swept', ...
        'couplings', struct('variant', 'histmatch_synth', 'n_total', 16384, ...
                             'grid', [16 16 16], 'dpRng', 30), 'label', 'x'), ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false));
results = struct('cfg_id', {}, 'group', {}, 'variant', {}, 'field_T', {}, 'cfg', {}, ...
    'out', {}, 'lattice_provenance', {});
[cfg_out, ~, dep_err] = invz_task2_resolve_cell_cfg(c, results, ion);
verifyTrue(testCase, isempty(dep_err));
verifyEqual(testCase, numel(cfg_out.couplings.Jnu_flat), 16384);
[~, info_ref] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30));
verifyEqual(testCase, cfg_out.couplings.J0eff, info_ref.Jcc0);
end

% ===========================================================================================
% (6) ISOLATED h_list DERIVATION HELPER -- pure, hand-built fixtures, no solves.
% ===========================================================================================
function test_derive_isolated_h_list_unique_stable_order(testCase)
swept_out.nodes = struct('h', {0.10, 0.20, 0.10, 0.30});
h_list = invz_task2_derive_isolated_h_list(swept_out);
verifyEqual(testCase, h_list, [0.10; 0.20; 0.30]);
verifyTrue(testCase, iscolumn(h_list));
end

function test_derive_isolated_h_list_errors_on_empty_nodes(testCase)
swept_out.nodes = struct([]);
verifyError(testCase, @() invz_task2_derive_isolated_h_list(swept_out), 'invz:task2Matrix');
end

% ===========================================================================================
% (7) MATRIX RUN: checkpoint/resume mechanics on a tiny, FAST, all-synthetic subset (stage-2c
% task 2b-driver gate: "Run + checkpoint + resume ... assert the results file is written
% incrementally and that a second run SKIPS the completed cells"). Uses the pinned 24-pt
% synthetic fixture (T=0.31, B=2.85, J0eff=6.42e-3 -- test_invz_task2_harness.m's own already-
% fast-proven config), NEVER the real 16^3 lattice, so this whole section stays fast/ungated.
% NEVER touches the real controller checkpoint path (.superpowers/sdd/task2_matrix_results.mat)
% -- always uses a fresh tempname()-based .mat path, cleaned up via onCleanup.
% ===========================================================================================
function test_matrix_run_checkpoint_incremental_and_resume(testCase)
[cfgA, cfgB] = tiny_synthetic_cfgs();
mk = @(id, cfg) struct('cfg_id', id, 'group', 'TEST', 'variant', 'x', 'field_T', 2.85, ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false), 'cfg', cfg);
cells = [mk('cellA', cfgA), mk('cellB', cfgB)];

results_path = [tempname() '.mat'];
cleanupObj = onCleanup(@() delete_if_exists(results_path)); %#ok<NASGU>

s1 = invz_task2_matrix_run(cells, struct('results_path', results_path, ...
    'cell_filter', {{'cellA'}}, 'verbose', false));
verifyEqual(testCase, s1.n_run, 1);
verifyEqual(testCase, s1.n_skipped, 0);
S = load(results_path, 'results');
verifyEqual(testCase, numel(S.results), 1);
verifyEqual(testCase, S.results(1).cfg_id, 'cellA');

% cellA must be SKIPPED (resume); cellB newly run+appended WITHOUT disturbing cellA's own
% already-checkpointed record -- the incremental-append contract.
s2 = invz_task2_matrix_run(cells, struct('results_path', results_path, ...
    'cell_filter', {{'cellA', 'cellB'}}, 'verbose', false));
verifyEqual(testCase, s2.n_run, 1);
verifyEqual(testCase, s2.n_skipped, 1);
S2 = load(results_path, 'results');
verifyEqual(testCase, numel(S2.results), 2);
verifyEqual(testCase, {S2.results.cfg_id}, {'cellA', 'cellB'});
verifyEqual(testCase, S2.results(1).out, S.results(1).out, ...
    'cellA''s checkpointed record must be untouched by the second run');

% full resume: re-running the SAME filter again must skip EVERYTHING.
s3 = invz_task2_matrix_run(cells, struct('results_path', results_path, ...
    'cell_filter', {{'cellA', 'cellB'}}, 'verbose', false));
verifyEqual(testCase, s3.n_run, 0);
verifyEqual(testCase, s3.n_skipped, 2);
end

function test_matrix_run_absorbs_failed_cell_and_continues(testCase)
[~, cfgB] = tiny_synthetic_cfgs();
ion = invz_ion();  Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
cfg_bad = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'isolated', 'h_list', 0.01, ...
    'seed_list', {{[]}}, 'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), ...
    'solve_opts', struct('transverse_mf', 'bogus_mode'), 'label', 'bad_cell');
mk = @(id, cfg) struct('cfg_id', id, 'group', 'TEST', 'variant', 'x', 'field_T', 2.85, ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false), 'cfg', cfg);
cells = [mk('bad_cell', cfg_bad), mk('good_cell', cfgB)];

results_path = [tempname() '.mat'];
cleanupObj = onCleanup(@() delete_if_exists(results_path)); %#ok<NASGU>

summary = invz_task2_matrix_run(cells, struct('results_path', results_path, 'verbose', false));
verifyEqual(testCase, summary.n_run, 2);
verifyEqual(testCase, summary.n_failed, 1);

S = load(results_path, 'results');
i_bad  = find(strcmp({S.results.cfg_id}, 'bad_cell'));
i_good = find(strcmp({S.results.cfg_id}, 'good_cell'));
verifyEqual(testCase, S.results(i_bad).out.status, 'failed');
verifyEqual(testCase, S.results(i_bad).out.err_id, 'invz:transverseMF');
verifyTrue(testCase, isempty(S.results(i_bad).out.nodes));
verifyFalse(testCase, isfield(S.results(i_good).out, 'status') && strcmp(S.results(i_good).out.status, 'failed'));
verifyTrue(testCase, S.results(i_good).out.nodes(1).accepted);
end

function test_matrix_run_absorbs_failed_dependency(testCase)
% The (expected-rare-in-the-real-matrix) case where an isolated cell's OWN prerequisite
% (swept sibling) itself failed: the dependent isolated cell must ALSO be checkpointed as a
% (driver-level) failed record, via invz:task2MatrixDependency, WITHOUT crashing the loop.
ion = invz_ion();  Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
cfg_swept_bad = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'swept', ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), ...
    'solve_opts', struct('transverse_mf', 'bogus_mode'), 'label', 'swept_bad');
cfg_iso = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'isolated', ...
    'h_list', [], 'seed_list', {{[]}}, ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'label', 'iso_dep_on_bad');
c1 = struct('cfg_id', 'swept_bad', 'group', 'TEST', 'variant', 'swept', 'field_T', 2.85, ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false), 'cfg', cfg_swept_bad);
c2 = struct('cfg_id', 'iso_dep_on_bad', 'group', 'TEST', 'variant', 'isolated_cold', ...
    'field_T', 2.85, 'depends_on', 'swept_bad', 'resolve', struct('seed2_from_dep', false), ...
    'cfg', cfg_iso);
cells = [c1, c2];

results_path = [tempname() '.mat'];
cleanupObj = onCleanup(@() delete_if_exists(results_path)); %#ok<NASGU>

summary = invz_task2_matrix_run(cells, struct('results_path', results_path, 'verbose', false));
verifyEqual(testCase, summary.n_run, 2);
verifyEqual(testCase, summary.n_failed, 2);

S = load(results_path, 'results');
i_swept = find(strcmp({S.results.cfg_id}, 'swept_bad'));
i_iso   = find(strcmp({S.results.cfg_id}, 'iso_dep_on_bad'));
verifyEqual(testCase, S.results(i_swept).out.status, 'failed');
verifyEqual(testCase, S.results(i_swept).out.err_id, 'invz:transverseMF');
verifyEqual(testCase, S.results(i_iso).out.status, 'failed');
verifyEqual(testCase, S.results(i_iso).out.err_id, 'invz:task2MatrixDependency');
end

function test_matrix_run_resolves_isolated_dependency_end_to_end(testCase)
% Full pipeline, on the cheap synthetic fixture: an isolated (multistart) cell's h_list AND
% seed are correctly derived from its swept sibling's OWN checkpointed output.
ion = invz_ion();  Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
swept_id = 'dep_swept';  iso_id = 'dep_isolated_seed2';
cfg_swept = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'swept', ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'label', swept_id);
cfg_iso = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'isolated', ...
    'h_list', [], 'seed_list', {{[]}}, ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'label', iso_id);
c1 = struct('cfg_id', swept_id, 'group', 'TEST', 'variant', 'swept', 'field_T', 2.85, ...
    'depends_on', '', 'resolve', struct('seed2_from_dep', false), 'cfg', cfg_swept);
c2 = struct('cfg_id', iso_id, 'group', 'TEST', 'variant', 'isolated_seed2', 'field_T', 2.85, ...
    'depends_on', swept_id, 'resolve', struct('seed2_from_dep', true), 'cfg', cfg_iso);
cells = [c1, c2];

results_path = [tempname() '.mat'];
cleanupObj = onCleanup(@() delete_if_exists(results_path)); %#ok<NASGU>

summary = invz_task2_matrix_run(cells, struct('results_path', results_path, 'verbose', false));
verifyEqual(testCase, summary.n_run, 2);
verifyEqual(testCase, summary.n_failed, 0);

S = load(results_path, 'results');
i_swept = find(strcmp({S.results.cfg_id}, swept_id));
i_iso   = find(strcmp({S.results.cfg_id}, iso_id));
swept_out = S.results(i_swept).out;
iso_rec   = S.results(i_iso);

expected_h = unique(reshape([swept_out.nodes.h], [], 1), 'stable');
verifyEqual(testCase, iso_rec.cfg.h_list, expected_h);

verifyEqual(testCase, iso_rec.cfg.seed_list{1}.K0s, swept_out.meta.J0eff);
verifyTrue(testCase, all(iso_rec.cfg.seed_list{1}.Sigma == 0));
[wn, ~, ~] = invz_matsubara(0.31, 40);
verifyEqual(testCase, numel(iso_rec.cfg.seed_list{1}.Sigma), numel(wn));

verifyEqual(testCase, numel(iso_rec.out.nodes), numel(expected_h));
end

% ===========================================================================================
% Local helpers (fixtures + cleanup)
% ===========================================================================================
function [cfgA, cfgB] = tiny_synthetic_cfgs()
%TINY_SYNTHETIC_CFGS Two fast, all-synthetic cfgs (the pinned 24-pt fixture, T=0.31/B=2.85/
% J0eff=6.42e-3 -- test_invz_task2_harness.m's own already-fast-proven config) for exercising
% the matrix RUNNER's own mechanics without ever touching the real 16^3 lattice. cfgB solves
% AT the independently-converged root pt.hmf (exactly test_invz_task2_harness.m's own
% test_isolated_mode_cold_matches_solve_point_ordered precedent) so it is GUARANTEED to
% accept -- not a guess at an arbitrary field value that might happen to land outside this
% fixture's convergent domain.
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;
pt = invz_solve_point_ordered(ion, 0.31, 2.85, Jnu, ...
    struct('J0eff', J0eff, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen'));
cfgA = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'swept', ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'label', 'cellA');
cfgB = struct('ion', ion, 'T', 0.31, 'Bx', 2.85, 'mode', 'isolated', 'h_list', pt.hmf, ...
    'seed_list', {{[]}}, 'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'label', 'cellB');
end

function delete_if_exists(p)
if exist(p, 'file'), delete(p); end
end

% ===========================================================================================
% REGRESSION (controller fix, stage-2c task 2b matrix-run pass): the TOP-LEVEL
% invz_task2_matrix was never called by any test -- the suite only exercised
% invz_task2_matrix_enumerate/_run directly, always with an explicit NON-EMPTY cell_filter.
% That left the DEFAULT (unfiltered, full-matrix) option path -- the exact path the controller
% runs -- unexercised, and it was broken: run_opts was built with
% struct(..., 'cell_filter', getf(opts,'cell_filter',{}), ...), and MATLAB's struct()
% CONSTRUCTOR treats a cell value as an array generator, so struct('cell_filter', {}) produced
% a 0x0 EMPTY STRUCT ARRAY rather than a scalar struct. Every downstream getf(run_opts,...)
% then failed with "Insufficient number of outputs from right hand side of equal sign"
% (observed: the real full-matrix launch died instantly at invz_task2_matrix_run.m:44).
% opts.dry_run reaches that construction without paying for 40 real solves.
% ===========================================================================================
function test_toplevel_default_opts_builds_scalar_run_opts(testCase)
s = invz_task2_matrix(struct('dry_run', true));      % NO cell_filter => the default path
verifyTrue(testCase, isstruct(s.run_opts) && isscalar(s.run_opts), ...
    'run_opts must be a SCALAR struct on the default path (struct(...,{}) empty-array gotcha)');
verifyTrue(testCase, iscell(s.run_opts.cell_filter) && isempty(s.run_opts.cell_filter), ...
    'default cell_filter must survive as an EMPTY CELL (= run every enumerated cell)');
verifyEqual(testCase, s.n_total_enumerated, 40, ...
    'the default path must still enumerate the full 40-cell matrix');
% the default checkpoint target is the controller''s real (git-ignored) results file
verifyTrue(testCase, ischar(s.run_opts.results_path) || isstring(s.run_opts.results_path));
verifyTrue(testCase, contains(char(s.run_opts.results_path), 'task2_matrix_results.mat'));
end

function cleanup_shadow_path(shadow_dir)
%CLEANUP_SHADOW_PATH onCleanup helper for test_resolve_cell_cfg_absorbs_nontask2_materialization_error:
% removes shadow_dir from the MATLAB path (if present) and deletes it, so the shadowed
% (throwaway) function definition can never leak into any other test.
if any(strcmp(strsplit(path, pathsep), shadow_dir))
    rmpath(shadow_dir);
    rehash path;
end
if exist(shadow_dir, 'dir'), rmdir(shadow_dir, 's'); end
end
