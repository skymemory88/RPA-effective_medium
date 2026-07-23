function [Jnu_flat, qc_sub, Jnu_unflat_sub, keep_idx] = ...
    invz_task2_couplings_downsample(Jnu_unflat_full, qc_full, stride)
%INVZ_TASK2_COUPLINGS_DOWNSAMPLE Deterministic, provenance-preserving density subsample of an
% already-computed real-lattice coupling set (stage-2c task 2b-driver; prereg
% docs/invzp_task2_prereg.md Sec. G item 3, brief G3: "deterministic subsets of the 16^3 q-set
% at density {1/2,1/4,1/8} ... isolates coupling density from distribution"; "subsample the
% q-grid / flattened Jnu by a fixed stride or fixed index set, documented; retain the
% flat->(q,branch) map so Task-0 provenance still resolves").
%
% Construction: a FIXED STRIDE on the flattened Q-ROW order, keep_idx = 1:stride:nq_full.
% Density 1/2 -> stride 2, 1/4 -> stride 4, 1/8 -> stride 8 (nq_full = 4096 for the standard
% 16^3 grid, divisible by 2/4/8 exactly, so numel(keep_idx) = nq_full/stride exactly). WHOLE
% q-rows are kept -- ALL 4 branches of every retained q, never split across branches -- so
% q/branch identity is preserved by construction: every kept row is an UNMUTATED copy of a
% specific row of the full-lattice Jnu_unflat/qc, and keep_idx doubles as the map back to that
% full-lattice row (invz_task2_run_config's own synthetic-Jnu_flat injection path cannot carry
% qc/Jnu_unflat through its harness -- see invz_task2_matrix.m's header for why the driver
% keeps this provenance itself rather than relying on invz_ordered_trace_resolve inside the
% harness). NO RNG anywhere -- purely a fixed integer index rule, deterministic and idempotent.
%
% Jnu_unflat_full [nq_full x 4], qc_full [nq_full x 3]: the FULL baseline lattice's own
%   unflattened branch matrix / q-grid (row-aligned), e.g. from
%   invz_task2_couplings_shifted_grid(ion, grid, dpRng, 'unshifted') or invz_jq_modes directly.
% stride (positive integer, e.g. 2/4/8): keep every stride-th q-row (1-based, first row always
%   kept).
%
% Jnu_flat        column-major flatten of Jnu_unflat_sub (matches invz_bz_couplings.m:17's own
%                 convention, so invz_ordered_trace_resolve's flat-index formula applies
%                 unchanged to a meta built from Jnu_unflat_sub/qc_sub).
% qc_sub          [nq_sub x 3], qc_full(keep_idx, :).
% Jnu_unflat_sub  [nq_sub x 4], Jnu_unflat_full(keep_idx, :).
% keep_idx        [nq_sub x 1] column vector of the ORIGINAL (full-lattice) row indices kept.
if ~(isscalar(stride) && isfinite(stride) && stride == round(stride) && stride >= 1)
    error('invz:task2Config', 'invz_task2_couplings_downsample: stride must be a positive integer; got %s.', mat2str(stride));
end
nq_full = size(Jnu_unflat_full, 1);
if size(qc_full, 1) ~= nq_full
    error('invz:task2Config', ['invz_task2_couplings_downsample: qc_full (%d rows) and ' ...
        'Jnu_unflat_full (%d rows) must be row-aligned.'], size(qc_full, 1), nq_full);
end
keep_idx = (1:stride:nq_full).';
Jnu_unflat_sub = Jnu_unflat_full(keep_idx, :);
qc_sub = qc_full(keep_idx, :);
Jnu_flat = Jnu_unflat_sub(:);
end
