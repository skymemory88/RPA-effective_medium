function res = invz_pole_continuity(h, d, y, tol)
%INVZ_POLE_CONTINUITY Actual-path G17 crossing-continuity check (Task 18 Step 1; the
% preregistered algorithm, implemented verbatim -- constants quoted in
% docs/invzp_strict_medium_gate0_report.md SS1). G = -chi (meV^-1), ferromagnetic positive J.
%
% Detects every crossing of the LOCAL Gstat denominator d = prof.gstat_local_denom along the
% actual solved profile h = prof.hgrid, and checks that the observable y (prof.r or prof.crit)
% is CONTINUOUS through it -- i.e. that the singularity is removable in the integrand, not a
% genuine discontinuity of the actual solved path (prereg SS3(d) vs SS5).
%
% h, d, y (REQUIRED, real finite numeric vectors, same length >= 2, h strictly increasing):
%   ALIGNED, FINITE inputs -- the precondition, not something this function filters for. A
%   caller MUST drop any non-finite (h,d,y) triple before calling (e.g. a domain-failed node
%   where r/crit/gstat_local_denom are NaN); this function errors loudly on a violation rather
%   than silently producing a wrong slope from an operation involving NaN, or from a duplicate
%   abscissa. NEVER uses interp1: duplicate or exact-zero abscissae are excluded by
%   construction below, not delegated to a library routine that may special-case them silently.
% tol (REQUIRED, positive finite scalar): pole_cont_tol, the frozen relative-jump limit.
%
% A crossing is either:
%   sign_change  d(i) and d(i+1) both nonzero with opposite sign. h_cross is located by linear
%                interpolation of d between the two nodes. yL/yR are independent linear
%                extrapolants of y to h_cross from the two nodes strictly left (i-1,i) and
%                strictly right (i+1,i+2) of the bracket; jump_lr = |yL-yR|/max(1,|yL|,|yR|).
%   exact_zero   d(i) == 0 exactly. h_cross = h(i) (a DIRECT evaluation AT the crossing, the
%                strongest evidence available -- not a special case to be waived). That node is
%                EXCLUDED from both fits: yL/yR are independent linear extrapolants from the
%                nearest two nodes strictly left (i-2,i-1) and strictly right (i+1,i+2),
%                EXCLUDING node i itself. y_node = y(i) is then checked against BOTH
%                extrapolants: jump_lr = |yL-yR|/max(1,|yL|,|yR|), jump_node_L =
%                |y_node-yL|/max(1,|y_node|,|yL|), jump_node_R = |y_node-yR|/max(1,|y_node|,|yR|).
%                A finiteness-only rule would admit a finite but arbitrarily wrong y_node
%                sitting between two mutually consistent extrapolants -- all THREE differences
%                are required <= tol, not merely that y_node be finite.
% A crossing is UNRESOLVED when fewer than two nodes exist strictly on either side (an edge
% effect, e.g. the crossing sits within 1-2 nodes of either end of the profile); it is
% reported, never silently dropped.
%
% res.status (exactly one of):
%   'no_crossing'   d never changes sign and is never exactly 0 -- nothing to check.
%   'unresolved'    at least one crossing lacks two-sided node coverage. Gate 0 FAILS on this
%                   (prereg SS5): "not silently omitted."
%   'jump_exceeded' every crossing is resolved and every jump finite, but at least one jump
%                   (of ANY kind above) exceeds tol. This is the status the exact-zero
%                   "displaced but finite" case must produce -- see the header note above.
%   'ok'            every crossing is resolved and every jump is finite and <= tol.
% Precedence when several crossings disagree: unresolved > jump_exceeded > ok (any single bad
% crossing spoils the whole call, matching the prereg's "not silently omitted").
%
% res.crossings: struct array (fixed schema, one entry per detected crossing):
%   kind ('sign_change'|'exact_zero'), index (i for exact_zero, [i i+1] for sign_change),
%   h_cross, resolved (logical), yL, yR, jump_lr, and -- populated for 'exact_zero' only,
%   NaN otherwise so every entry shares one schema -- y_node, jump_node_L, jump_node_R.
% res.max_jump: max over every finite jump value reported above (NaN when there are no
%   crossings, or when every crossing is unresolved). Diagnostic only -- res.status is the gate.
if nargin < 4
    error('invz:poleContinuityInput', 'invz_pole_continuity requires (h, d, y, tol).');
end
name = {'h', 'd', 'y'};
vecs = {h, d, y};
n = numel(h);
for k = 1:3
    v = vecs{k};
    if ~(isnumeric(v) && isreal(v) && isvector(v) && numel(v) == n && all(isfinite(v(:))))
        error('invz:poleContinuityInput', ['%s must be a real finite numeric vector of length ' ...
            '%d (aligned with h); got class %s, numel=%d, all finite=%d.'], ...
            name{k}, n, class(v), numel(v), all(isfinite(v(:))));
    end
end
if ~(isnumeric(tol) && isscalar(tol) && isfinite(tol) && tol > 0)
    error('invz:poleContinuityInput', 'tol must be a positive finite scalar; got %s.', mat2str(tol));
end
if n < 2
    error('invz:poleContinuityInput', 'h must have at least 2 nodes; got %d.', n);
end
h = h(:).';  d = d(:).';  y = y(:).';    % canonical row orientation throughout
if any(diff(h) <= 0)
    error('invz:poleContinuityInput', ...
        'h must be strictly increasing (no duplicate/reordered abscissae); got diff(h) = %s.', ...
        mat2str(diff(h)));
end

blank = struct('kind', '', 'index', [], 'h_cross', NaN, 'resolved', false, ...
                'yL', NaN, 'yR', NaN, 'jump_lr', NaN, ...
                'y_node', NaN, 'jump_node_L', NaN, 'jump_node_R', NaN);
crossings = blank([]);   % 0x0 struct, fixed schema, grown below

for i = 1:n
    if d(i) == 0
        c = blank;  c.kind = 'exact_zero';  c.index = i;  c.h_cross = h(i);  c.y_node = y(i);
        if i-2 >= 1 && i+2 <= n
            c.resolved = true;
            c.yL = local_extrap(h(i-2), y(i-2), h(i-1), y(i-1), c.h_cross);
            c.yR = local_extrap(h(i+1), y(i+1), h(i+2), y(i+2), c.h_cross);
            c.jump_lr     = local_relgap(c.yL, c.yR);
            c.jump_node_L = local_relgap(c.y_node, c.yL);
            c.jump_node_R = local_relgap(c.y_node, c.yR);
        end
        crossings(end+1) = c; %#ok<AGROW>
    elseif i < n && d(i+1) ~= 0 && sign(d(i)) ~= sign(d(i+1))
        % both nonzero (d(i)==0 handled above; d(i+1)==0 is its OWN exact_zero crossing at the
        % next i, never double-counted here) with opposite sign -> a bracketed sign change.
        c = blank;  c.kind = 'sign_change';  c.index = [i i+1];
        c.h_cross = h(i) + (h(i+1) - h(i)) * (0 - d(i)) / (d(i+1) - d(i));
        if i-1 >= 1 && i+2 <= n
            c.resolved = true;
            c.yL = local_extrap(h(i-1), y(i-1), h(i), y(i), c.h_cross);
            c.yR = local_extrap(h(i+1), y(i+1), h(i+2), y(i+2), c.h_cross);
            c.jump_lr = local_relgap(c.yL, c.yR);
        end
        crossings(end+1) = c; %#ok<AGROW>
    end
end

if isempty(crossings)
    res = struct('status', 'no_crossing', 'crossings', crossings, 'max_jump', NaN);
    return;
end

all_jumps = [];
any_unresolved = false;
for k = 1:numel(crossings)
    c = crossings(k);
    if ~c.resolved
        any_unresolved = true;
        continue;
    end
    if strcmp(c.kind, 'exact_zero')
        all_jumps = [all_jumps, c.jump_lr, c.jump_node_L, c.jump_node_R]; %#ok<AGROW>
    else
        all_jumps = [all_jumps, c.jump_lr]; %#ok<AGROW>
    end
end

if any_unresolved
    status = 'unresolved';
elseif ~isempty(all_jumps) && (any(~isfinite(all_jumps)) || max(all_jumps) > tol)
    status = 'jump_exceeded';   % non-finite OR too-large jump on an otherwise "resolved" crossing
else
    status = 'ok';
end

if isempty(all_jumps)
    max_jump = NaN;
else
    max_jump = max(all_jumps);
end
res = struct('status', status, 'crossings', crossings, 'max_jump', max_jump);
end

% ---------------------------------------------------------------------------------------------
function ye = local_extrap(h1, y1, h2, y2, hq)
%LOCAL_EXTRAP Linear extrapolant through (h1,y1)-(h2,y2), evaluated at hq. h1 ~= h2 is
% guaranteed by the caller (strictly increasing h, distinct integer indices).
ye = y1 + (y2 - y1) * (hq - h1) / (h2 - h1);
end

% ---------------------------------------------------------------------------------------------
function r = local_relgap(a, b)
%LOCAL_RELGAP |a-b| normalised by max(1,|a|,|b|) (prereg SS5's shared normalisation).
r = abs(a - b) / max([1, abs(a), abs(b)]);
end
