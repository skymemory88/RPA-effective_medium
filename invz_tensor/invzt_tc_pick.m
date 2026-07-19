function [act, ka, kb, ncross] = invzt_tc_pick(cv)
%INVZT_TC_PICK Pure crossing/slide decision on ascending-T valid crit votes.
%   [act, ka, kb, ncross] = INVZT_TC_PICK(cv) inspects the VALID samples
%   (criticality votes in ascending-T order: cv > 0 paramagnet, < 0 ordered,
%   == 0 exactly critical) from one INVZT_CRITICAL_T window pass and decides:
%     act = 'zero'    cv(ka) == 0: the sample IS the boundary; kb = NaN.
%     act = 'bracket' voters ka, kb = ka+1 bracket the highest-T
%                     ordered->para crossing (cv(ka) < 0 < cv(kb)).
%     act = 'up'      the highest-T voter is ordered: the requested
%                     highest-T ordered->para crossing lies ABOVE the
%                     window (the physical high-T side is paramagnetic).
%                     This includes the lower para->ordered leg of a
%                     re-entrant region, which the projected classifier
%                     abandons (inherited gap, closed here).
%     act = 'down'    every voter is paramagnetic: boundary below.
%   ncross counts the boundary indicators in the window: adjacent STRICT
%   sign flips plus exactly-critical RUNS (a zero-run is one root, never
%   double-counted); ncross > 1 signals candidate re-entrance (caller warns).
%
%   Exact roots: a sampled cv == 0 is recognized and, when it is the
%   highest-T root, returned in preference to refining a lower interval
%   (the projected classifier mis-reads [-1, 0, 1] as two sign changes with
%   no returnable crossing -- inherited gap, closed here).
%
%   PURE: no solves, no state -- millisecond-testable on synthetic votes.
%   Precondition: cv is a nonempty numeric vector of finite values (the
%   caller passes only valid votes; invalid samples were already dropped).
%
%   See also INVZT_CRITICAL_T.
ka = NaN;  kb = NaN;
cv = cv(:).';                                   % row, defensively
z      = (cv == 0);
nzruns = sum(diff([false, z]) == 1);            % exactly-critical RUNS (roots)
strict = sum(cv(1:end-1).*cv(2:end) < 0);       % adjacent strict sign flips
ncross = strict + nzruns;
if cv(end) < 0
    act = 'up';                                 % top voter ordered: boundary above
    return;
end
iz  = find(z);                                  % exact roots at samples
upk = find(cv(1:end-1) < 0 & cv(2:end) > 0);    % strict ordered->para pairs
if ~isempty(iz) && (isempty(upk) || iz(end) > upk(end) + 1)
    act = 'zero';  ka = iz(end);                % highest-T root is an exact zero
elseif ~isempty(upk)
    act = 'bracket';  ka = upk(end);  kb = ka + 1;
else
    act = 'down';                               % top PM, no zeros, no up: all PM
end
end
