function [interval, status, Bc] = invz_boundary_interval(fields, orderedMask, pmMask, unknownMask)
%INVZ_BOUNDARY_INTERVAL Anchor-bracketed 1/z boundary reduction with an explicit status.
%   [interval, status, Bc] = invz_boundary_interval(fields, orderedMask, pmMask, unknownMask)
%   locates the LAST ordered anchor and the FIRST stable-PM anchor and reports the field
%   interval that brackets the 1/z boundary, together with how trustworthy that bracket is.
%
%   interval  [lo hi] = [fields(last ordered)  fields(first PM)];  [NaN NaN] when there is
%             no bracket.
%   status    'valid'        bracketed, and NO unknown/indeterminate column lies strictly
%                            between the two anchors -- the bracket is the sweep resolution.
%             'widened'      bracketed, but at least one unknown/indeterminate column lies
%                            BETWEEN the anchors, so the true boundary is only known to
%                            within the wider interval.
%             'unbracketed'  no anchor pair with the PM anchor above the ordered one
%                            (includes "nothing is labelled at all").
%             'invalid'      the three masks are not mutually exclusive: some column carries
%                            more than one label. That is a CALLER predicate defect, not a
%                            physics outcome, so no boundary is invented from it.
%   Bc        midpoint of interval on 'valid'/'widened'; NaN otherwise.
%
% WHY THIS IS NOT min/max OVER ALL UNKNOWN COLUMNS. An unknown column on the far side of the
% sweep says nothing about where the boundary is; folding it into the interval would widen a
% perfectly resolved bracket with unrelated information. Only INTERVENING unknowns widen.
%
% This is a SECOND helper, not a replacement for invz_boundary_field: that one remains the
% source of the historical scalar S.Bc_1z on the default 'resummed' path (G9), and this one
% feeds S.Bc_1z only under a strict static-medium scheme.
%
% INPUT CONTRACT (hard -- violations are wiring errors and throw 'invz:boundaryInterval'):
% fields must be finite, real and STRICTLY INCREASING, because "last ordered" and "first PM"
% are otherwise presentation-order artefacts; a caller wanting another display order must
% reorder after this reduction, not before it. All three masks must match numel(fields).
if ~isnumeric(fields) || ~isreal(fields) || isempty(fields) || ~all(isfinite(fields(:)))
    error('invz:boundaryInterval', 'fields must be a nonempty finite real vector.');
end
fields = fields(:).';
n = numel(fields);
if n > 1 && ~all(diff(fields) > 0)
    error('invz:boundaryInterval', ...
        'fields must be STRICTLY INCREASING; reorder before this reduction.');
end
masks = {orderedMask, pmMask, unknownMask};
names = {'orderedMask', 'pmMask', 'unknownMask'};
for k = 1:3
    m = masks{k};
    if ~(islogical(m) || isnumeric(m)) || numel(m) ~= n
        error('invz:boundaryInterval', '%s must be a %d-element logical mask.', names{k}, n);
    end
    masks{k} = logical(m(:)).';
end
ord = masks{1};  pm = masks{2};  unk = masks{3};

interval = [NaN NaN];  Bc = NaN;

% Mutually exclusive labels: each column carries exactly one reason string upstream, so an
% overlap can only mean the caller built its predicates wrong. Report it, never anchor on it.
if any((double(ord) + double(pm) + double(unk)) > 1)
    status = 'invalid';
    return;
end

kb = find(ord, 1, 'last');  ka = find(pm, 1, 'first');
if isempty(kb) || isempty(ka) || ka <= kb
    status = 'unbracketed';       % also the all-unknown / nothing-labelled case
    return;
end

interval = [fields(kb) fields(ka)];
Bc = 0.5*(fields(kb) + fields(ka));
if any(unk((kb+1):(ka-1)))
    status = 'widened';
else
    status = 'valid';
end
end
