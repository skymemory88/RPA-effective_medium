function res = invz_phase1_offset_spread(kind, values8, J_ref)
%INVZ_PHASE1_OFFSET_SPREAD Item 6, offset-sensitivity spread + pairwise agreement gate
% (docs/invzp_phase1_quadrature_prereg.md item 6): "Report the spread of the item-5 summaries
% across the eight {0,1/2}^3 offsets ... apply the pass/fail agreement gate (item-5 tolerances)."
% Spread is defined here as max-min of the (already normalized, for kind=='shape') statistic
% across the eight offsets; the agreement gate is applied PAIRWISE (all 28 pairs of 8 offsets)
% using the SAME invz_phase1_gate formula item 5's refinement uses -- "apply the pass/fail
% agreement gate (item-5 tolerances)" reuses those tolerances verbatim, not a separate formula.
%
% kind      'shape' (values8 already normalized: stat/J_ref or stat/J_ref^2) | 'energy' (raw meV).
% values8   [1x8] or [8x1], the statistic's value at each of the eight canonical offsets, in
%           invz_phase1_offsets() order (order does not affect max/min/spread/pairwise-gate, only
%           provenance/reporting).
% J_ref     required for kind=='energy' (passed through to invz_phase1_gate); pass [] for 'shape'.
%
% OUTPUT res (struct):
%   .max, .min, .spread     max(values8), min(values8), max-min.
%   .pairwise_pass          logical, true iff EVERY one of the 28 offset pairs agrees within the
%                           item-5 tolerance for this statistic.
%   .worst_diff, .worst_tol, .worst_margin   the single worst pair's |Dv|, tolerance, and
%                           |Dv|-tolerance (<=0 everywhere means PASS) -- report-only detail.
values8 = values8(:).';
if numel(values8) ~= 8 || ~all(isfinite(values8))
    error('invz:phase1Config', 'invz_phase1_offset_spread: values8 must be 8 finite scalars; got %s.', mat2str(values8));
end
if nargin < 3, J_ref = []; end

res.max    = max(values8);
res.min    = min(values8);
res.spread = res.max - res.min;

worst_margin = -Inf;  worst_diff = 0;  worst_tol = 0;  pairwise_pass = true;
for a = 1:8
    for b = (a+1):8
        [p, d, t] = invz_phase1_gate(kind, values8(a), values8(b), J_ref);
        if ~p, pairwise_pass = false; end
        margin = d - t;
        if margin > worst_margin
            worst_margin = margin;  worst_diff = d;  worst_tol = t;
        end
    end
end
res.pairwise_pass = pairwise_pass;
res.worst_diff   = worst_diff;
res.worst_tol    = worst_tol;
res.worst_margin = worst_margin;
end
