function Bc = invz_boundary_field(fields, isBelow, isAbove)
%INVZ_BOUNDARY_FIELD midpoint estimate of a phase boundary from per-field labels.
% Bc between the last field labelled "below" and the first field labelled "above";
% fields assumed ascending (the drivers use linspace). NaN when the sweep does not
% bracket the flip. The CALLER owns the predicates: pass only VALID labels (e.g. exclude
% suspect auto-PM points, which then WIDEN the bracket instead of anchoring it -- QCP
% Stage 1). Precision is half the gap between the two bracketing labelled fields:
% masked/suspect/other-labelled columns between them widen it beyond half the field
% step. Overlapping/nonmonotone label sequences (reentrant or noisy near-boundary
% labels) return NaN whenever the first above-label precedes the last below-label --
% refine the sweep or inspect manually; use invz_critical for a solver-grade Bc.
Bc = NaN;
kb = find(isBelow, 1, 'last');  ka = find(isAbove, 1, 'first');
if ~isempty(kb) && ~isempty(ka) && ka > kb
    Bc = 0.5*(fields(kb) + fields(ka));
end
end
