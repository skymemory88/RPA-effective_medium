function [c, ok, pt] = invzt_crit_at(ion, T, B, lat, opts)
%INVZT_CRIT_AT One tensor criticality sample: crit at (T, B) + sample VALIDITY.
%   [c, ok, pt] = INVZT_CRIT_AT(ion, T, B, lat, opts) solves one A1 tensor
%   point (INVZT_SOLVE_POINT, opts forwarded verbatim) and returns the
%   criticality scalar c = pt.crit plus ok, the three-part sample-validity
%   verdict (converged, finite crit, Sigma0 >= getf(opts,'sigma_floor',-0.5)
%   -- the floor single-sourced with INVZT_CRITICAL, rejecting the spurious
%   negative-Sigma fixed point). The third output pt is the solved point
%   struct (an EMPTY struct when the sample was absorbed by the catch below)
%   -- INVZT_CRITICAL_T's branch-tracked sampler reads pt.Sigma from it for
%   its rolling warm-start seed (execution finding E1).
%
%   ok is VALIDITY-only -- deliberately NO crit > 0 term. Each consumer
%   applies its own phase logic: INVZT_TC_PM_EXTRAP filters the PM points
%   itself (metastable ordered-side samples are FILTERED there, never
%   asserted on); INVZT_CRITICAL_T lets every valid sample vote by sign(c).
%
%   Physics signals (invz:degenerateDoublet, invz:orderedPhase,
%   invzt:a1ZeroField) are absorbed as c = NaN, ok = false; every other
%   identifier (invzt:mode, invzt:a1OddGamma, ...) is a MISCONFIGURATION and
%   rethrows -- absorbing it would silently bias a production sweep by
%   reading a config error as "ordered side" (mirrors the projected
%   invz_crit_at rule).
%
%   See also INVZT_SOLVE_POINT, INVZT_CRITICAL_T, INVZT_TC_PM_EXTRAP,
%   INVZ_CRIT_AT (projected reference).
try
    pt = invzt_solve_point(ion, T, B, lat, opts);
    c  = pt.crit;
    ok = pt.converged && isfinite(c) && pt.Sigma0 >= getf(opts, 'sigma_floor', -0.5);
catch err
    switch err.identifier
        case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
            c = NaN;  ok = false;      % phase/physics signal, not an error
            pt = struct();             % absorbed sample: no point to expose
        otherwise
            rethrow(err);              % misconfiguration: surface it
    end
end
end
