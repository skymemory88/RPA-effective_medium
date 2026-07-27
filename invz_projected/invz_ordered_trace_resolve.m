function [q_idx, branch_idx, Jq] = invz_ordered_trace_resolve(meta, flat_idx)
%INVZ_ORDERED_TRACE_RESOLVE Flat coupling index -> (q-index, branch-index, J(q)).
% Uses an invz_ordered_trace.m run's meta block (meta.nq, meta.Jnu_unflat) to invert the
% flattening invz_bz_couplings.m:17 performs (Jnu = Jnu(:), column-major over
% invz_jq_modes' [nq x 4] branch matrix): flat k -> q = mod(k-1,nq)+1,
% branch = floor((k-1)/nq)+1 (both 1-based). Kept here as ONE canonical implementation
% rather than re-derived ad hoc at each call site.
%
% Errors ('invz:orderedTraceResolve') when the run has no lattice provenance (synthetic
% Jnu_flat, meta.is_synthetic true) or when flat_idx is not a valid finite index (in
% particular: the per-iteration idx_pos_flat/idx_neg_flat trace fields are NaN when no
% mode of that sign exists -- callers must check isfinite(flat_idx) before resolving).
if meta.is_synthetic || isempty(meta.Jnu_unflat) || ~isfinite(meta.nq)
    error('invz:orderedTraceResolve', ['this run has no lattice provenance (synthetic ' ...
        'Jnu_flat) -- flat index %s cannot be resolved to (q, branch).'], mat2str(flat_idx));
end
if ~(isscalar(flat_idx) && isfinite(flat_idx) && flat_idx >= 1 ...
        && flat_idx <= numel(meta.Jnu_unflat) && flat_idx == round(flat_idx))
    error('invz:orderedTraceResolve', 'flat_idx must be an integer in [1, %d]; got %s.', ...
        numel(meta.Jnu_unflat), mat2str(flat_idx));
end
nq = meta.nq;
q_idx      = mod(flat_idx - 1, nq) + 1;
branch_idx = floor((flat_idx - 1) / nq) + 1;
Jq = meta.Jnu_unflat(q_idx, branch_idx);
end
