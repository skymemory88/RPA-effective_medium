function Jnu_flat = invz_task2_couplings_histmatch_synth(real_ref, n_total)
%INVZ_TASK2_COUPLINGS_HISTMATCH_SYNTH Deterministic histogram-matched synthetic coupling set
% (stage-2c task 2b-driver; prereg docs/invzp_task2_prereg.md Sec. G item 5, brief G5:
% "synthetic couplings whose histogram matches the REAL Jnu distribution shape at real
% cardinality (deterministic construction, e.g. inverse-CDF sampling on a fixed quantile
% grid -- NO RNG). Lattice-diagnostic.").
%
% Implements exactly that recipe: the empirical CDF of real_ref (sorted values) is inverse-
% sampled at n_total FIXED, evenly-spaced quantiles u_k = (k-1/2)/n_total, k=1..n_total (the
% standard midpoint/plotting-position convention -- avoids the u=0/u=1 endpoints so no
% extrapolation past real_ref's observed min/max is ever needed), via linear interpolation of
% the empirical quantile function. NO RANDOM DRAW anywhere: a deterministic function
% evaluation only, so the same (real_ref, n_total) always reproduces bit-identical output.
%
% Degenerate, EXACTLY-CHECKABLE property when n_total == numel(real_ref) (the actual G5 use:
% matching the real lattice's OWN cardinality): u_k*n_real + 0.5 = k exactly for every k, so
% linear interpolation lands exactly on integer nodes and Jnu_flat reproduces sort(real_ref)
% EXACTLY (up to floating point) -- this is the correct, expected behaviour, not a bug: at
% equal cardinality, a histogram-matched, order-preserving-shape resample of a finite dataset
% IS that dataset (any relabelling of which physical q/branch got which value does not change
% the 1-D value histogram at all), and the DIFFERENCE this construction deliberately makes
% from the real lattice is that the (q,branch) <-> value PAIRING is discarded/decorrelated
% (assigned by sorted-value order here), isolating whether the raw value distribution alone
% (vs. its spatial/dispersion correlation across q) drives the physical masking behaviour.
% test_invz_task2_matrix.m's test_histmatch_synth_matches_histogram_shape_at_real_cardinality
% checks this property directly.
%
% real_ref  numeric vector, the REAL Jnu_flat distribution to match the shape of (e.g. the
%           real 16^3-grid, dpRng=30 lattice's own flattened Jnu).
% n_total   target cardinality (e.g. numel(real_ref) itself, for the real-cardinality-matched
%           G5 cells).
%
% Jnu_flat  [n_total x 1] column vector, histogram-shape-matched to real_ref.
real_sorted = sort(real_ref(:));
n_real = numel(real_sorted);
if n_real < 2
    error('invz:task2Config', 'invz_task2_couplings_histmatch_synth: real_ref needs >= 2 values.');
end
if ~(isscalar(n_total) && isfinite(n_total) && n_total == round(n_total) && n_total >= 1)
    error('invz:task2Config', 'invz_task2_couplings_histmatch_synth: n_total must be a positive integer; got %s.', mat2str(n_total));
end
n_total = round(n_total);
u = ((1:n_total).' - 0.5) / n_total;              % fixed quantile grid, deterministic, NO RNG
pos = u * n_real + 0.5;                           % map to the empirical CDF's index scale
pos = min(max(pos, 1), n_real);                   % clamp to the valid index range [1, n_real]
Jnu_flat = interp1((1:n_real).', real_sorted, pos, 'linear');
end
