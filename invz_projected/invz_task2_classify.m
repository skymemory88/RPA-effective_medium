function [class, s, D_uni, Dq_min] = invz_task2_classify(node)
%INVZ_TASK2_CLASSIFY Per-node stability classification (stage-2c task 2a; prereg
% docs/invzp_task2_prereg.md Sec. A, FROZEN 2026-07-23 -- values below are cited VERBATIM,
% not re-derived or tuned).
%
% Stability scalar: s = min(D_uni, min(Dq)).  Classification:
%   stable       s > +1e-3
%   marginal     |s| <= 1e-3
%   unstable     s < -1e-3
%   unconverged  separate classification: checker NOT accepted, OR non-finite
%                (prereg Sec. A table, literal fourth row)
%
% THE central rule (stage2c-context.md; prereg Sec. A prose): "Only checker-accepted
% exported states are classified stable/marginal/unstable. A node whose solve is not
% checker-accepted is unconverged, full stop -- it is NOT labelled 'unstable' no matter how
% negative its transient Dq/D_uni got during iteration. Negative transient iterates are
% diagnostic evidence, not proof of an unstable solution." This function enforces that:
% node.accepted is checked FIRST and unconditionally routes to 'unconverged' when false,
% before s is ever compared against the +-1e-3 thresholds.
%
% node (struct):
%   .accepted  (logical, REQUIRED) -- invz_ordered_residual's res.accepted verbatim (or the
%              equivalent invz_ordered_node_solve info.accepted).
%   .D_uni     (scalar, REQUIRED) -- the uniform static pole observable (checker's own
%              independent Block-B recomputation, res.D_uni / info.so.D_uni).
%   .Dq        (vector, OPTIONAL) -- the full per-mode static denominator array
%              (Dq = 1 + (J(q)-K0).*Gstat); min(Dq) is reduced HERE, NaN-propagating (see
%              robust_min below -- MATLAB's bare min() silently ignores individual NaN
%              entries, exactly the silent-masking failure mode
%              invz_ordered_residual.m's own robust_max_abs was written to avoid).
%   .Dq_min    (scalar, OPTIONAL, used only when .Dq is absent) -- an already-reduced
%              min(Dq) the caller computed itself; still finite-checked here.
%   Exactly one of .Dq / .Dq_min must be supplied.
%
% Returns the class plus the components (s, D_uni, min(Dq)) for reporting -- ALWAYS
% returned (even when class == 'unconverged'), since a non-accepted node's transient
% s/D_uni/Dq_min are diagnostic evidence Task 2's report must still show, just never used
% to justify a stable/marginal/unstable label.
THR = 1e-3;   % prereg Sec. A, FROZEN -- do not tune, do not reinterpret.

if ~isfield(node, 'accepted')
    error('invz:task2Classify', 'node.accepted is required.');
end
if ~isfield(node, 'D_uni')
    error('invz:task2Classify', 'node.D_uni is required.');
end
D_uni = node.D_uni;

has_Dq     = isfield(node, 'Dq') && ~isempty(node.Dq);
has_Dq_min = isfield(node, 'Dq_min');
if has_Dq
    Dq_min = robust_min(node.Dq(:));
elseif has_Dq_min
    Dq_min = node.Dq_min;
    if ~isfinite(Dq_min), Dq_min = NaN; end
else
    error('invz:task2Classify', 'node must carry either .Dq (vector) or .Dq_min (scalar).');
end

finite_ok = isfinite(D_uni) && isfinite(Dq_min);
if finite_ok
    s = min(D_uni, Dq_min);   % both operands finite: MATLAB's NaN-ignoring min() quirk cannot fire
else
    s = NaN;                 % never fabricate a numeric s from a missing/non-finite component
end

if ~node.accepted || ~finite_ok
    class = 'unconverged';   % prereg Sec. A: "checker NOT accepted, OR non-finite" -- either
    return;                  % condition alone is sufficient; s above is diagnostic only.
end

if s > THR
    class = 'stable';
elseif abs(s) <= THR
    class = 'marginal';
else
    class = 'unstable';
end
end

% =================================================================================================
function m = robust_min(v)
%ROBUST_MIN min(v), but NaN-propagating (mirrors invz_ordered_residual.m's robust_max_abs):
% MATLAB's bare min() silently ignores individual NaN entries in a vector rather than
% reporting the reduction as undefined -- exactly the silent-masking failure mode this
% diagnostic stage must not repeat. Any non-finite entry makes the whole reduction NaN.
if any(~isfinite(v))
    m = NaN;
else
    m = min(v);
end
end
