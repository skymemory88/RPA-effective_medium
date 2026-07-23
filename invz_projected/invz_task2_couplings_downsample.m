function [Jnu_flat, qc_sub, Jnu_unflat_sub, keep_idx] = ...
    invz_task2_couplings_downsample(Jnu_unflat_full, qc_full, grid, stride)
%INVZ_TASK2_COUPLINGS_DOWNSAMPLE Deterministic, provenance-preserving density subsample of an
% already-computed real-lattice coupling set (stage-2c task 2b-driver; prereg
% docs/invzp_task2_prereg.md Sec. G item 3, brief G3: "deterministic subsets of the 16^3 q-set
% at density {1/2,1/4,1/8} ... isolates coupling density from distribution"; "subsample the
% q-grid / flattened Jnu by a fixed stride or fixed index set, documented; retain the
% flat->(q,branch) map so Task-0 provenance still resolves").
%
% Construction (REVIEW-FIX, stage-2c task 2b-driver review pass -- replaces an earlier flat-
% stride construction that confounded density with distribution, see "why not a flat stride"
% below): MODULAR decimation on the INTEGER per-axis BZ grid indices. Recover each q-row's
% integer axis indices (h,k,l), each in {0,...,grid(a)-1}, by inverting qVec_generator's own
% inclusive-linspace grid construction (generate_grid: qx=linspace(range(1),range(2),grid(1)),
% same for qy/qz, [QX,QY,QZ]=meshgrid(qx,qy,qz), qvec=[QX(:) QY(:) QZ(:)], range=[-0.5 0.5] --
% the SAME hardcoded range invz_bz_couplings.m/invz_task2_couplings_shifted_grid.m both use;
% confirmed against an independent index-space meshgrid at implementation time, not assumed --
% see task2b-driver-report.md's review-fix pass), then KEEP every q-row whose indices satisfy
% mod(h+k+l, stride) == 0. ALL 4 branches of every kept q-row are kept (rows are never split
% across branches), so q/branch identity is preserved exactly as before.
%
% Why not a flat stride (the bug this replaces): qVec_generator's grid-mode q-list is a
% COLUMN-MAJOR meshgrid(:) flatten, so a CONSTANT STRIDE on the flat row index (keep_idx =
% 1:stride:nq_full) keeps ALL 16 values of two axes and decimates ONLY the fastest-varying
% (2nd, "k") axis -- empirically: stride 2 -> a 16x8x16 sub-cuboid, stride 4 -> 16x4x16,
% stride 8 -> 16x2x16 (verified numerically). That is a single-axis collapse, NOT a uniform
% lower-density BZ sample: it substitutes a direction-biased distribution for the "lower
% density" leg and CONFOUNDS density with distribution, defeating this section's stated
% purpose (isolate density from distribution) for every G3 cell.
%
% Why the modular rule fixes it (EXACT properties, verified numerically, not merely argued):
% for a 16-point axis and stride in {2,4,8} (16 is exactly divisible by each), each axis's
% index modulo stride is EXACTLY uniform over its 16 values (4/2/... hits per residue class),
% so h+k+l mod stride is ALSO exactly uniform over the full 16^3 grid (the sum of independent
% uniform-on-Z_stride variables is itself uniform on Z_stride). Consequences:
%   (a) DENSITY: numel(keep_idx) == nq_full/stride identically (2048/1024/512 of 4096 for
%       stride 2/4/8) -- exact, no off-by-one, no rounding.
%   (b) MARGINALS: for EVERY fixed value of ANY one axis, the other two axes' indices modulo
%       stride are themselves uniform, so at least one (in fact exactly nq_full/(16*stride))
%       kept row exists at every one of that axis's 16 values -- i.e. each of the three BZ
%       axes retains its FULL 16-value marginal range after decimation, for every stride. This
%       is the anti-single-axis-collapse guarantee the flat-stride construction violated.
% Both properties are asserted in test_invz_task2_matrix.m's
% test_downsample_preserves_rows_and_determinism for stride in {2,4,8}.
%
% NO RNG anywhere -- purely a fixed integer index rule on data already in hand (plus a
% defensive round-trip residual check on the recovered indices, see below), deterministic and
% idempotent: the same inputs always produce byte-identical outputs.
%
% Jnu_unflat_full [nq_full x 4], qc_full [nq_full x 3]: the FULL baseline lattice's own
%   unflattened branch matrix / q-grid (row-aligned), e.g. from
%   invz_task2_couplings_shifted_grid(ion, grid, dpRng, 'unshifted') or invz_jq_modes directly.
%   MUST be the UNSHIFTED grid: the index recovery inverts qVec_generator's plain inclusive-
%   linspace(-0.5,0.5,grid(a)) construction, so a shifted qc_full would silently mis-recover
%   indices (caught by the residual check below, not silently accepted) -- G3 only ever calls
%   this helper with the 'unshifted' rung, by construction (see invz_task2_resolve_cell_cfg.m).
% grid      [1x3] positive-integer (each >= 2) BZ grid size that qc_full/Jnu_unflat_full were
%   built on (e.g. [16 16 16] -- the SAME grid argument passed to
%   invz_task2_couplings_shifted_grid/invz_bz_couplings). Threaded in explicitly by the caller
%   (not inferred from qc_full via e.g. numel(unique(.))) so a partially Gamma-filtered or
%   otherwise incomplete qc_full cannot silently under-count an axis; the nq_full == prod(grid)
%   check below enforces this link.
% stride (positive integer, e.g. 2/4/8): keep every q-row whose three integer axis indices sum
%   to a multiple of stride (see construction above).
%
% Jnu_flat        column-major flatten of Jnu_unflat_sub (matches invz_bz_couplings.m:17's own
%                 convention, so invz_ordered_trace_resolve's flat-index formula applies
%                 unchanged to a meta built from Jnu_unflat_sub/qc_sub).
% qc_sub          [nq_sub x 3], qc_full(keep_idx, :).
% Jnu_unflat_sub  [nq_sub x 4], Jnu_unflat_full(keep_idx, :).
% keep_idx        [nq_sub x 1] column vector, ASCENDING, of the ORIGINAL (full-lattice) row
%                 indices kept -- still the exact map back to the full-lattice row Task-0
%                 provenance needs (unchanged contract from the earlier flat-stride version).
if ~(isscalar(stride) && isfinite(stride) && stride == round(stride) && stride >= 1)
    error('invz:task2Config', 'invz_task2_couplings_downsample: stride must be a positive integer; got %s.', mat2str(stride));
end
if ~(isnumeric(grid) && isequal(size(grid), [1 3]) && all(isfinite(grid)) && all(grid == round(grid)) && all(grid >= 2))
    error('invz:task2Config', 'invz_task2_couplings_downsample: grid must be a 1x3 vector of integers >= 2; got %s.', mat2str(grid));
end
nq_full = size(Jnu_unflat_full, 1);
if size(qc_full, 1) ~= nq_full
    error('invz:task2Config', ['invz_task2_couplings_downsample: qc_full (%d rows) and ' ...
        'Jnu_unflat_full (%d rows) must be row-aligned.'], size(qc_full, 1), nq_full);
end
if nq_full ~= prod(grid)
    error('invz:task2Config', ['invz_task2_couplings_downsample: Jnu_unflat_full has %d rows but ' ...
        'grid %s implies %d q-points; qc_full/Jnu_unflat_full must be the FULL, un-Gamma-' ...
        'filtered unshifted grid (see header) for the axis-index recovery below to be valid.'], ...
        nq_full, mat2str(grid), prod(grid));
end

% --- recover integer per-axis grid indices by inverting qVec_generator's own inclusive
% linspace(range(1),range(2),grid(a)) construction (range = [-0.5 0.5], the SAME hardcoded
% convention invz_bz_couplings.m/invz_task2_couplings_shifted_grid.m both use throughout this
% task). The round-trip residual check immediately below RE-CONFIRMS the convention on every
% call (rather than trusting it silently) -- a wrong index recovery would silently re-break the
% density/marginal isolation this function exists to guarantee, so this fails LOUDLY instead. --
range = [-0.5 0.5];
step = (range(2) - range(1)) ./ (grid - 1);          % [1x3] per-axis step (isotropic for [16 16 16])
idx = round((qc_full - range(1)) ./ step);           % [nq_full x 3] integer, each column in 0..grid(a)-1
resid = qc_full - (range(1) + idx .* step);
max_resid = max(abs(resid(:)));
if max_resid > 1e-9
    error('invz:task2Config', ['invz_task2_couplings_downsample: qc_full does not lie on the ' ...
        'expected inclusive linspace(%g,%g,grid) grid (max index-recovery residual %.3g) -- ' ...
        'refusing to guess the axis-index convention (a wrong recovery would silently re-break ' ...
        'the density/marginal isolation this function exists to guarantee).'], range(1), range(2), max_resid);
end

keep_mask = mod(sum(idx, 2), stride) == 0;
keep_idx = find(keep_mask);
Jnu_unflat_sub = Jnu_unflat_full(keep_idx, :);
qc_sub = qc_full(keep_idx, :);
Jnu_flat = Jnu_unflat_sub(:);
end
