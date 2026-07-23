function cells = invz_task2_matrix_enumerate(ion)
%INVZ_TASK2_MATRIX_ENUMERATE Enumerates every prereg Sec. G experiment-matrix cell at the
% Sec. E frozen physical fields (stage-2c task 2b-driver). PURE, deterministic, NO I/O, NO
% solves, NO physics calls of any kind: returns a struct array ready for
% invz_task2_matrix_run.m. Calling this function twice yields an IDENTICAL list (same order,
% same cfg_id strings, same cfg contents) -- there is no Date/now/rand/tempname anywhere in
% this file, only fixed literals and integer/string arithmetic, matching the brief's
% requirement that cfg_id be a stable deterministic string with no timestamp/RNG dependence.
%
% Frozen physical fields (prereg Sec. E; the four Task-2a-derived/cited values, used VERBATIM
% here per that report's own table): F = [0.25, 0.55, 0.80]*Bc_PM (Bc_PM = 4.692769 T, from
% Task-2a's invz_critical PM-mass root-find on the real 16^3/dpRng=30 lattice) plus the
% existing 2.85 T defect anchor. This function CITES those already-derived numbers (frozen,
% not re-derived here -- re-deriving Bc_PM requires an actual invz_critical solve, which is
% emphatically not "pure/no I/O", and would make calling this enumerator twice needlessly
% repeat an ~6 s physics computation for a purely bookkeeping operation). The independent
% re-derivation + tolerance check the brief also asks for ("re-derive ... and assert it
% matches ... OR use the frozen values citing 2a -- state which") lives instead in
% test_invz_task2_matrix.m's own test_bc_pm_rederivation_matches_frozen_fields, run once,
% separately, as an explicit regression/confidence check -- NOT baked into every enumerate()
% call. So: THIS function uses the frozen literal values (citing Task-2a); the companion test
% file re-derives and asserts agreement. Both of the brief's offered options are exercised,
% at the appropriate cost/location for each.
%
% ion  (optional, default invz_ion()) -- threaded into every cell's cfg.ion so every solve
%      uses the identical ion parameters.
%
% cells  struct array, one element per matrix cell:
%   .cfg_id      (char) stable, unique, deterministic identifier (field/group/variant tag).
%   .group       (char) 'G1G2' | 'G3' | 'G4' | 'G5' | 'G6'.
%   .variant     (char) short human-readable sub-tag within the group.
%   .field_T     (double) the applied transverse field (T) this cell runs at.
%   .depends_on  (char) cfg_id of a prerequisite cell that MUST be run/checkpointed first
%                (currently: the G1/G2 isolated companions depend on their field's own G1
%                'swept' cell, for the h_list/multistart-seed derivation -- see
%                invz_task2_resolve_cell_cfg.m); '' when the cell is fully self-contained.
%   .resolve     (struct) what invz_task2_resolve_cell_cfg.m must still fill in at run time:
%                .seed2_from_dep (logical) -- construct the multistart seed from the
%                dependency's own J0eff (see invz_task2_resolve_cell_cfg.m's header).
%   .cfg         (struct) ready for invz_task2_run_config, MODULO the deferred fields above
%                (isolated cfg.h_list/cfg.seed_list, and any cfg.couplings.variant recipe --
%                both resolved by invz_task2_resolve_cell_cfg.m immediately before the cell
%                is actually run).
%
% Total cell count (matches the brief's own approximate tally): G1G2 12 + G3 12 + G4 2 +
% G5 2 + G6 12 = 40.
if nargin < 1 || isempty(ion), ion = invz_ion(); end

FIELDS = [1.173192, 2.581023, 3.754215, 2.85];   % T; prereg Sec. E (0.25/0.55/0.80*Bc_PM,
                                                  % Bc_PM=4.692769 T, + the 2.85 T defect
                                                  % anchor) -- Task-2a report table, verbatim.
T0    = 0.1;         % K; prereg Sec. E.
GRID0 = [16 16 16];  % baseline BZ grid (prereg "Notation": invz_bz_couplings default).
DPRNG0 = 30;         % baseline real-space dipole cutoff (prereg "Notation" baseline).

BASE24  = linspace(-2e-3, 6.0e-3, 24).';   % the pinned Task-2a/prereg synthetic fixture,
J0EFF24 = 6.42e-3;                         % VERBATIM (test_invz_task2_harness.m's own values).

build = {};   % cell array of scalar-struct cell-records; concatenated into a struct array once,
              % at the very end (mirrors invz_task2_run_config.m's own recs{}->[recs{:}] idiom).

% =============================================================================================
% G1/G2: isolated-vs-swept + cold-vs-continued + multistart (brief: ~3 x 4 = 12 cells).
% swept = the continuation sweep (invz_ordered_trace -> invz_hmf_ordered), the SAME baseline
% real-lattice construction used throughout this matrix's "full/unshifted/dp30" rung
% (G3's density ladder and G6's offset ladder both reuse THIS SAME cfg_id as their own
% baseline rung rather than duplicating it -- see the G3/G6 blocks below).
% isolated_cold / isolated_seed2 DEPEND on their field's own swept cell (depends_on) for the
% h-values to re-solve in isolation -- see invz_task2_derive_isolated_h_list.m's header for
% why "the swept run's own visited h's" is the only interpretation that answers "does the
% SAME node converge alone" regardless of whether the swept run itself masks or converges.
% =============================================================================================
for iF = 1:numel(FIELDS)
    F = FIELDS(iF);  ftag = field_tag(F);

    swept_id = sprintf('g1_swept_%s', ftag);
    cfg_swept = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'swept', ...
        'couplings', struct('grid', GRID0, 'dpRng', DPRNG0), 'label', swept_id);
    build{end+1} = make_cell(swept_id, 'G1G2', 'swept', F, '', cfg_swept); %#ok<AGROW>

    iso_cold_id = sprintf('g1_isolated_cold_%s', ftag);
    cfg_iso_cold = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'isolated', ...
        'h_list', [], 'seed_list', {{[]}}, ...
        'couplings', struct('grid', GRID0, 'dpRng', DPRNG0), 'label', iso_cold_id);
    build{end+1} = make_cell(iso_cold_id, 'G1G2', 'isolated_cold', F, swept_id, cfg_iso_cold); %#ok<AGROW>

    iso_seed2_id = sprintf('g1_isolated_seed2_%s', ftag);
    cfg_iso_seed2 = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'isolated', ...
        'h_list', [], 'seed_list', {{[]}}, ...   % placeholder; overwritten via resolve.seed2_from_dep
        'couplings', struct('grid', GRID0, 'dpRng', DPRNG0), 'label', iso_seed2_id);
    build{end+1} = make_cell(iso_seed2_id, 'G1G2', 'isolated_seed2', F, swept_id, cfg_iso_seed2, ...
        struct('seed2_from_dep', true)); %#ok<AGROW>
end

% =============================================================================================
% G3: physical downsampling, deterministic subsets of the 16^3 q-set at density {1/2,1/4,1/8},
% swept, at each field (brief: full = G1, so only the 3 NEW densities are enumerated here;
% ~3 x 4 = 12 cells). stride 2/4/8 -> density 1/2/1/4/1/8 (see
% invz_task2_couplings_downsample.m).
% =============================================================================================
STRIDES = [2 4 8];
for iF = 1:numel(FIELDS)
    F = FIELDS(iF);  ftag = field_tag(F);
    for is_ = 1:numel(STRIDES)
        stride = STRIDES(is_);
        cfg_id = sprintf('g3_ds%d_swept_%s', stride, ftag);
        cfg_ = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'swept', ...
            'couplings', struct('variant', 'downsample', 'stride', stride, ...
                                 'grid', GRID0, 'dpRng', DPRNG0), ...
            'label', cfg_id);
        build{end+1} = make_cell(cfg_id, 'G3', sprintf('downsample_stride%d', stride), F, '', cfg_); %#ok<AGROW>
    end
end

% =============================================================================================
% G4: cardinality-controlled synthetic (duplicate the pinned 24-pt fixture up to the real
% cardinality 16384), swept, at >= 1 representative condition (brief: ~2 cells;
% lattice-diagnostic only, NOT an existence endpoint per prereg Sec. E/G). Representative
% fields chosen here: the 2.85 T defect anchor (the SAME field the pinned 24-pt fixture is
% already validated at, test_invz_task2_harness.m) and 0.25*Bc_PM = 1.173192 T (the field
% closest to zero, for a low-field vs. mid/high-field representative pair) -- an
% implementer's choice (flagged for review), since the brief leaves "representative" open.
% =============================================================================================
G4_FIELD_IDX = [4, 1];
for k = 1:numel(G4_FIELD_IDX)
    iF = G4_FIELD_IDX(k);  F = FIELDS(iF);  ftag = field_tag(F);
    cfg_id = sprintf('g4_cardsynth_swept_%s', ftag);
    cfg_ = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'swept', ...
        'couplings', struct('variant', 'cardinality_synth', 'n_total', 16384, ...
                             'base_vals', BASE24, 'J0eff', J0EFF24), ...
        'label', cfg_id);
    build{end+1} = make_cell(cfg_id, 'G4', 'cardinality_synth', F, '', cfg_); %#ok<AGROW>
end

% =============================================================================================
% G5: histogram-matched synthetic (synthetic couplings whose histogram matches the REAL Jnu
% distribution shape at real cardinality), swept, at >= 1 representative condition (brief:
% ~2 cells; lattice-diagnostic only). SAME representative-field choice as G4, for a directly
% comparable G4-vs-G5 pair at each field (isolating "cardinality alone" vs. "cardinality +
% real shape" against the SAME real-lattice baseline).
% =============================================================================================
G5_FIELD_IDX = [4, 1];
for k = 1:numel(G5_FIELD_IDX)
    iF = G5_FIELD_IDX(k);  F = FIELDS(iF);  ftag = field_tag(F);
    cfg_id = sprintf('g5_histmatch_swept_%s', ftag);
    cfg_ = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'swept', ...
        'couplings', struct('variant', 'histmatch_synth', 'n_total', 16384, ...
                             'grid', GRID0, 'dpRng', DPRNG0), ...
        'label', cfg_id);
    build{end+1} = make_cell(cfg_id, 'G5', 'histmatch_synth', F, '', cfg_); %#ok<AGROW>
end

% =============================================================================================
% G6: grid offsets x dpRng ladder (brief: offset {unshifted,half-step} x dpRng {30,40},
% unshifted/dp30 = G1's own swept cell, so only the 3 NEW combos are enumerated here; swept,
% at each field; ~3 x 4 = 12 cells). See invz_task2_couplings_shifted_grid.m for the DELICATE
% half-step/dipolar-Gamma handling.
% =============================================================================================
G6_COMBOS = {{'unshifted', 40}, {'half_step', 30}, {'half_step', 40}};
for iF = 1:numel(FIELDS)
    F = FIELDS(iF);  ftag = field_tag(F);
    for ig = 1:numel(G6_COMBOS)
        sm = G6_COMBOS{ig}{1};  dp = G6_COMBOS{ig}{2};
        if strcmp(sm, 'half_step'), offtag = 'halfstep'; else, offtag = 'unshifted'; end
        cfg_id = sprintf('g6_off%s_dp%d_swept_%s', offtag, dp, ftag);
        cfg_ = struct('ion', ion, 'T', T0, 'Bx', F, 'mode', 'swept', ...
            'couplings', struct('variant', 'shifted_grid', 'shift_mode', sm, ...
                                 'grid', GRID0, 'dpRng', dp), ...
            'label', cfg_id);
        build{end+1} = make_cell(cfg_id, 'G6', sprintf('%s_dp%d', sm, dp), F, '', cfg_); %#ok<AGROW>
    end
end

cells = [build{:}];
end

% =================================================================================================
function c = make_cell(cfg_id, group, variant, field_T, depends_on, cfg, resolve_)
if nargin < 7 || isempty(resolve_), resolve_ = struct('seed2_from_dep', false); end
c = struct('cfg_id', cfg_id, 'group', group, 'variant', variant, 'field_T', field_T, ...
    'depends_on', depends_on, 'resolve', resolve_, 'cfg', cfg);
end

% =================================================================================================
function tag = field_tag(F)
%FIELD_TAG Deterministic, unique, filename/identifier-safe field tag: micro-Tesla integer
% encoding (all four prereg Sec. E fields are given to exactly 6 decimal places, so this
% round-trips exactly). E.g. 1.173192 -> 'F1173192'; 2.85 -> 'F2850000'.
tag = sprintf('F%07d', round(F * 1e6));
end
