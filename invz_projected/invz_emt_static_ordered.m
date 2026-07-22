function [K0, Gstat, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu_flat, K0_seed, beta, J0eff, G0inel0, G0el0, opts)
%INVZ_EMT_STATIC_ORDERED Static-sector EMT closure for the ordered elastic propagator (SS3, P0-2/P0-A).
% The elastic G(0) of J 2.28-2.29 breaks the ordinary Dyson structure, so the closed-form
% direct solve of INVZ_EMT_SCALAR does not apply at w = 0 for m ~= 0. This function solves
% the scalar fixed point demanded by the EMT definitions (J 2.10-2.11 / HTML 16):
%   Gstat(K0)  = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0)
%   G(q,0)     = Gstat ./ (1 + (J(q) - K0).*Gstat)                       (HTML 17 insertion)
%   closure:   mean_q G(q,0) = Gstat   and   K0 = mean_q(J.*Gq)/mean_q(Gq)
% by damped iteration on K0. The static weights are the caller's: production passes the
% FULL-electronuclear split (round-2 P0-A), so at m = 0 (G0el0 -> 0) the fixed point
% coincides with invz_emt_scalar's direct solve ON THE FULL PROPAGATOR -- the PM solver's
% own K(0) (gate C1; the load-bearing loop-level identity is Task 3's Gate 6b). lam and
% Sigma0 are FIXED inputs here; the caller's outer loop refreshes them with the closed K0
% written into the K vector's w = 0 slot (see invz_hmf_ordered / the jensen solver mode).
% out: xi/h0/G0bare/Gtil0/r at the closed K0 (from invz_gstat_ordered), D_uni = the uniform
% static inverse response 1 + (J0eff - K0)*Gstat (pole observable), resid, iters, converged.
% opts (threaded from the caller as ONE NAMED STRUCT emt_static -- P1-F, never a bare
% struct() at call sites): resid_tol (default 1e-10, meV^-1 -- the PRIMARY convergence
% criterion AND the documented closure-residual acceptance threshold callers gate on),
% tol (default 0 -- an OPTIONAL absolute |dK0| stall floor; 0 means machine-resolution-only,
% see below), maxit (default 200), mix (default 0.5). out.converged is measured on the
% EXPORTED (K0, Gstat, resid) tuple (out.resid < rtol), so it can never disagree with
% out.resid.
% Amendment round 1 (controller resolution of the round-2 P0-2 BLOCKED escalation): the
% original |dK0|-only stopping test under-delivers closure by the map's local conditioning
% (resid/|dK0| ~ 1.8e5 at the gate-C2 operating point), so the residual itself -- the
% physical closure quantity -- became the primary convergence test.
% Amendment round 2 (round 1 was itself defective on two counts, both fixed here): (a) a
% fixed dK stall floor (1e-12) structurally preempted the residual criterion given the same
% ~1.8e5 conditioning, so the stall floor is now machine resolution, 4*eps(|K0|) (opts.tol
% is an optional additional absolute floor, default 0, i.e. inert unless the caller raises
% it); (b) the residual check now breaks BEFORE any K0 update, so the exported K0 is exactly
% the one that closed, and out.converged is recomputed post-loop from out.resid so it can
% never be a half-step out of sync with what was actually exported.
if nargin < 10, opts = struct(); end
rtol  = getf(opts, 'resid_tol', 1e-10);  % PRIMARY: closure-residual convergence (meV^-1)
tol   = getf(opts, 'tol', 0);            % optional ABSOLUTE |dK0| stall floor (0 = machine-only)
maxit = getf(opts, 'maxit', 200);
mix   = getf(opts, 'mix', 0.5);
Jf = Jnu_flat(:);
K0 = K0_seed;
for it = 1:maxit
    Gs = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
    Gq = Gs ./ (1 + (Jf - K0).*Gs);
    Gbar = mean(Gq);
    if abs(Gbar - Gs) < rtol, break; end % closed at the CURRENT K0 -- exported as-is
    K0_new = mean(Jf .* Gq) / Gbar;
    dK = abs(K0_new - K0);
    if dK < max(tol, 4*eps(abs(K0)))     % TRUE stall: no representable progress possible
        break;
    end
    K0 = K0 + mix*(K0_new - K0);
end
[Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
out = go;
out.D_uni = 1 + (J0eff - K0)*Gstat;
out.resid = abs(mean(Gq) - Gstat);
out.iters = it;
out.converged = out.resid < rtol;        % measured on the EXPORTED tuple
if ~out.converged
    warning('invz:emtStatic', 'static closure not converged after %d iterations: resid = %.3g', it, out.resid);
end
end
