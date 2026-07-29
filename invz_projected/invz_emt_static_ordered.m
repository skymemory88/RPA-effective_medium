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
% see below), maxit (default 200), mix (default 0.5), warn (default true -- emits the
% invz:emtStatic non-convergence warning; in-loop sweep callers pass false because they
% gate on out.converged and would otherwise flood the console). out.converged is measured
% on the EXPORTED (K0, Gstat, resid) tuple (out.resid < rtol), so it can never disagree
% with out.resid.
% G = -chi (meV^-1), ferromagnetic positive J.
% opts.static_medium ('resummed' default): under 'strict_1z_dyson_ref'/'strict_1z_bare_ref' this
% function does NOT iterate. K0 = Jbar - mu2*Gref once; the selected Dyson reference is
% Gref = (G0inel0+G0el0)/(1+Sigma0), while the explicit bare comparator omits that division.
% mix/maxit/tol and the invz:emtStatic warning are inapplicable and K0_seed is IGNORED.
% out.resid becomes the algebraic |K0 - Kstrict| check (zero for a correct call), out.iters = 0,
% and out.converged reports DOMAIN validity via out.medium_status
% ('ok' | 'ref_denom_nonpositive' | 'ref_denom_small' | 'nonfinite'). D_uni and Dq are still
% built from the physical Gstat and reported through Dq_min/Dq_max/Dq_neg_count.
% Strict mode requests invz_gstat_ordered's stable-form reassociation; resummed mode preserves
% the historical seven-argument arithmetic bitwise. opts.Jmom / opts.ref_margin as in
% invz_emt_scalar.
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
warn  = getf(opts, 'warn', true);            % emit the invz:emtStatic non-convergence warning
Jf = Jnu_flat(:);
smid = getf(opts, 'static_medium', 'resummed');
% opts.closure_coordinate ('factored' default | 'defactored'): which ARITHMETIC the
% resummed q-average is evaluated in. 'factored' is the historical seven-argument form,
% preserved bitwise. 'defactored' evaluates the SAME closure through
% INVZ_RECIPROCAL_STATIC_CLOSURE (Gq = 1/(1/Gstat + J(q) - K0)), which is an exact
% algebraic reassociation -- no floor, no broadening, no clipped denominator -- and which
% therefore evaluates the finite limit at the local Gstat pole where the factored form
% produces Inf/Inf = NaN. The CLOSURE CONDITION and the exported out.resid are unchanged
% (|mean_q Gq - Gstat| < resid_tol in both), so acceptance semantics do not move; only the
% iterates' behaviour on a pole crossing does. Default-off (blind_convg_plan.md SS2.6, SS4
% promotion boundary): it must not change a production phase label until its gates pass.
defactored = strcmp(getf(opts, 'closure_coordinate', 'factored'), 'defactored');
if ~ismember(getf(opts, 'closure_coordinate', 'factored'), {'factored', 'defactored'})
    error('invz:emtStatic', ...
        'opts.closure_coordinate must be ''factored'' or ''defactored''; got %s.', ...
        char(string(getf(opts, 'closure_coordinate', 'factored'))));
end
medium = struct('scheme', smid, 'status', 'not_applicable', 'ref', [], 'closure', []);
if strcmp(smid, 'resummed')
    K0 = K0_seed;
    for it = 1:maxit
        if defactored
            [Gs, gof] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, ...
                                           struct('stable_form', true));
            rc = invz_reciprocal_static_closure(gof, Jf, K0);
            Gbar = rc.Gbar;
            if abs(Gbar - Gs) < rtol, break; end % SAME closure test as the factored branch
            K0_new = rc.Jloc;
        else
            Gs = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
            Gq = Gs ./ (1 + (Jf - K0).*Gs);
            Gbar = mean(Gq);
            if abs(Gbar - Gs) < rtol, break; end % closed at the CURRENT K0 -- exported as-is
            K0_new = mean(Jf .* Gq) / Gbar;
        end
        dK = abs(K0_new - K0);
        if dK < max(tol, 4*eps(abs(K0)))     % TRUE stall: no representable progress possible
            break;
        end
        K0 = K0 + mix*(K0_new - K0);
    end
else
    % --- STRICT O(1/z): ONE-SHOT, NO ITERATION (spec SS1, SSB) -----------------------------
    % K0 = Jbar - mu2*Gref with Gref = G0bare/(1+Sigma0), a K0/lambda/xi-INDEPENDENT
    % reference. There is no fixed point here, so mix/maxit/tol and the invz:emtStatic
    % non-convergence warning are all inapplicable, and K0_seed is deliberately IGNORED: a
    % one-shot medium has no warm start, so a contaminated seed cannot propagate between nodes.
    % The resummed q-average this replaces carries the denominator that dies below Bc; its
    % feedback into K0 is what exceeds retained order.
    if ~isvector(Jnu_flat)
        error('invz:staticMedium', ['invz_emt_static_ordered: static_medium ''%s'' does not ' ...
            'support a [nJ,nw] coupling matrix in this phase (spec SS7.5): the pre-existing ' ...
            'Jnu_flat(:) flattening would average every frequency column into the static ' ...
            'q-average.'], smid);
    end
    if isfield(opts, 'Jmom') && ~isempty(opts.Jmom)
        mom = opts.Jmom;                          % threaded once per resolved point / node
    else
        mom = invz_coupling_moments(Jf);          % compatibility fallback
    end
    [Gref, refi] = invz_static_medium_reference(G0inel0 + G0el0, Sigma0, smid, ...
        struct('ref_margin', getf(opts, 'ref_margin', 1e-6)));
    [K0, clo] = invz_medium_moment_closure(Gref, mom, smid);
    medium.ref = refi;  medium.closure = clo;
    if strcmp(refi.status, 'ok') && strcmp(clo.status, 'ok')
        medium.status = 'ok';
    elseif ~strcmp(refi.status, 'ok')
        medium.status = refi.status;
    else
        medium.status = clo.status;
    end
    it = 0;                                       % nothing iterated
end
strict_mode = ~strcmp(medium.status, 'not_applicable');
if strict_mode && ~strcmp(medium.status, 'ok')
    % Do not feed an invalid reference/K0 into Jensen's local denominators. The caller will
    % stop before lambdas/Sigma consume it.
    Gstat = NaN;
    go = struct('xi', NaN, 'h0', NaN, 'G0bare', G0inel0 + G0el0, ...
                'Gtil0', NaN, 'r', NaN, 'gstat_local_denom', NaN, ...
                'G0inel0', G0inel0, 'G0el0', G0el0);   % same field set as invz_gstat_ordered
else
    if strict_mode || defactored
        [Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, ...
                                         struct('stable_form', true));
    else
        [Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
    end
end
if defactored && ~strict_mode
    Gbar_out = invz_reciprocal_static_closure(go, Jf, K0).Gbar;
else
    Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
    Gbar_out = mean(Gq);
end
out = go;
out.D_uni = 1 + (J0eff - K0)*Gstat;      % collective observable, built from the PHYSICAL Gstat
out.Dq_min = min(1 + (Jf - K0).*Gstat);
out.Dq_max = max(1 + (Jf - K0).*Gstat);
out.Dq_neg_count = nnz(1 + (Jf - K0).*Gstat <= 0);
out.iters = it;
out.medium = medium;  out.medium_status = medium.status;
out.omit_mu3 = NaN;  out.omit_cubic = NaN;  out.omit_max = NaN;
if strcmp(medium.status, 'not_applicable')
    out.resid = abs(Gbar_out - Gstat);   % resummed: the q-average closure residual
    out.converged = out.resid < rtol;    % measured on the EXPORTED tuple
    if warn && ~out.converged
        warning('invz:emtStatic', 'static closure not converged after %d iterations: resid = %.3g', it, out.resid);
    end
else
    % STRICT: the load-bearing residual is the ALGEBRAIC K0 check, not the closure of the
    % discarded resummation (spec SS4.4). Running that discarded inner solve here would restore
    % the very iteration and pole exposure this scheme removes, so it is not computed.
    % out.converged reports DOMAIN validity: there is no iteration to converge.
    out.resid = abs(K0 - medium.closure.Kstrict);
    out.converged = strcmp(medium.status, 'ok');
    out.omit_mu3 = medium.closure.omit_mu3;
    out.omit_cubic = medium.closure.omit_cubic;
    out.omit_max = medium.closure.omit_max;
end
end
