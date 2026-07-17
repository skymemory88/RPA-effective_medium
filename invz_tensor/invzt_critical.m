function [Bc, out] = invzt_critical(ion, T, lat, Brange, opts)
%INVZT_CRITICAL PM-side critical transverse field Bc at fixed T via sign(crit) bisection.
%   [Bc, out] = INVZT_CRITICAL(ion, T, lat, Brange, opts) locates the paramagnet /
%   ordered phase boundary in transverse field at fixed temperature T, using the A1
%   tensor point solver invzt_solve_point on the LOCKED lattice struct lat. Brange =
%   [Blo Bhi] (Tesla) brackets the search: Bhi the high-field paramagnet, Blo the
%   low-field ordered floor. Bisection returns Bc, the field where crit crosses 0.
%
%   PHASE KEY = sign(pt.crit), NOT pt.converged (Task-6 finding, LOCKED). The A1
%   Anderson solver CONVERGES metastable paramagnetic fixed points in the ORDERED
%   phase -- there pt.converged is true but pt.crit < 0. So the boundary is where
%   pt.crit changes sign (crit > 0 = PM, crit < 0 = ordered); convergence +
%   finiteness are only a per-sample VALIDITY precondition, never the phase verdict.
%
%   SAMPLE VALIDITY (all three required for a sample to cast a crit-signed vote):
%     (1) pt.converged, (2) isfinite(pt.crit), (3) pt.Sigma0 >= opts.sigma_floor.
%   Guard (3) rejects the SPURIOUS NEGATIVE-Sigma fixed point the Anderson map can
%   land on at some depths (Task-6 finding: Sigma should be >~ 0 in the paramagnet;
%   a strongly negative Sigma0 reports a misleading crit). An INVALID sample reads
%   as the ORDERED side (its crit is discarded, never defines the boundary) -- this
%   is the "non-convergence = phase signal" rule: away from the boundary the PM
%   solve converges cleanly, so a failure/near-singular RPA near the search floor is
%   the ordered signal (Task-4 carry: a near-singular RPA / NaN is caught, not a
%   crash). A genuine tensor-layer machinery error (invzt:*) is a MISCONFIGURATION,
%   not a phase signal, and rethrows (mirrors invz_crit_at's structural-error rule);
%   a physics signal from the shared engine (e.g. invz:degenerateDoublet) reads as
%   ordered.
%
%   WARM START: the converged pt.Sigma of each solve seeds the next (opts.Sigma_seed
%   threaded). T is fixed, so the Matsubara length is constant and the seed always
%   fits; the seed accelerates nearby solves and helps the ordered-side solves land
%   on the metastable PM branch (giving clean crit < 0 votes for the crossing).
%
%   opts (getf defaults; every field also FORWARDS to invzt_solve_point, which
%   ignores the control-only fields):
%     tol          0.02   field-bracket tolerance (T); Bc resolved to ~tol.
%     maxit        40      max bisection steps.
%     sigma_floor  -0.5    Sigma0 below this is spurious/unphysical -> sample INVALID
%                          (physical PM Sigma0 ~ +0.3..+0.5; well clear of the floor).
%   All other opts (odd, chi_rest, mode, Ecut, hyp, transverse_mf, ...) pass to the
%   solver unchanged, so a Bc scan runs in whatever ODD/chi_rest configuration the
%   caller selects.
%
%   Bc (crossing estimate): regula-falsi on crit between the bracketing CONVERGED
%   points when the ordered end is itself a valid converged crit < 0 sample (the
%   metastable-PM case); otherwise the tight-bracket midpoint (ordered end
%   non-converged -- the field edge of the converged paramagnet, matching the
%   projected invz_critical fallback).
%
%   out: .iters (per-sample struct array: .B, .side (+1 PM / -1 ordered), .valid,
%   .pt), .Bc, .T, .bracket = [Blo Bhi] (final), .method ('interp'|'midpoint').
%   out.iters(end).pt is guaranteed a CONVERGED PM point (the confirmed PM edge that
%   brackets Bc from above).
%
%   ERRORS invzt:bracket when Brange does not bracket a crossing: the top is not a
%   converged PM point (raise Bhi), or the floor is already PM (lower Blo).
%
%   See also INVZT_SOLVE_POINT, INVZT_TC_PM_EXTRAP, INVZ_CRITICAL (projected
%   field-root reference), INVZ_REFINE_CROSSING.
if nargin < 5, opts = struct(); end
if numel(Brange) ~= 2 || ~(Brange(2) > Brange(1))
    error('invzt:bracket', 'Brange must be [Blo Bhi] with Bhi > Blo (got %s).', mat2str(Brange));
end
Blo0 = Brange(1);  Bhi0 = Brange(2);
tol      = getf(opts, 'tol', 0.02);
maxit    = getf(opts, 'maxit', 40);
sigfloor = getf(opts, 'sigma_floor', -0.5);

iters = struct('B', {}, 'side', {}, 'valid', {}, 'pt', {});   % growable sample log
seed  = [];                                                   % warm-start Sigma

% --- top of window: must be a validated PARAMAGNET (valid AND crit > 0) ----------
[sHi, vHi, ptHi] = sample(Bhi0);
if ~(vHi && ptHi.crit > 0)
    error('invzt:bracket', ...
        ['Top of field window B = %.3f T at T = %.3f K is not a converged ' ...
         'paramagnet (valid = %d, crit = %.3g, Sigma0 = %.3g): raise the window top.'], ...
        Bhi0, T, vHi, ptHi.crit, ptHi.Sigma0);
end
Bhi = Bhi0;  ptPM = ptHi;                     % Bhi tracks the lowest-field PM point

% --- floor of window: must be ORDERED (side -1: crit < 0 or invalid) -------------
[sLo, vLo, ptLo] = sample(Blo0);              %#ok<ASGLU> vLo kept for readability
if sLo > 0
    error('invzt:bracket', ...
        ['Floor of field window B = %.3f T at T = %.3f K is already paramagnetic ' ...
         '(crit = %.3g > 0): no PM/ordered crossing in [%.3f, %.3f] T; lower the floor.'], ...
        Blo0, T, ptLo.crit, Blo0, Bhi0);
end
Blo = Blo0;  ptOrd = ptLo;                    % Blo tracks the highest-field ordered point

% --- bisection on sign(crit) -----------------------------------------------------
for it = 1:maxit
    if Bhi - Blo < tol, break; end
    Bm = 0.5*(Blo + Bhi);
    [sM, ~, ptM] = sample(Bm);
    if sM > 0
        Bhi = Bm;  ptPM = ptM;                % paramagnet: tighten from above
    else
        Blo = Bm;  ptOrd = ptM;               % ordered/invalid: tighten from below
    end
end

% --- crossing estimate -----------------------------------------------------------
critPM = ptPM.crit;                                       % > 0 at Bhi
ordValid = ptOrd.converged && isfinite(ptOrd.crit) && ...
           ptOrd.crit < 0 && ptOrd.Sigma0 >= sigfloor;
if ordValid
    critOrd = ptOrd.crit;                                 % < 0 at Blo
    Bc = Bhi - critPM*(Blo - Bhi)/(critOrd - critPM);     % regula-falsi on crit
    Bc = min(max(Bc, Blo), Bhi);                          % clamp into the bracket
    method = 'interp';
else
    Bc = 0.5*(Blo + Bhi);                                 % PM-edge midpoint fallback
    method = 'midpoint';
end

% --- guarantee out.iters ends on a CONVERGED PM point (the confirmed PM edge that
%     brackets Bc from above); the bisection may have ended on an ordered sample.
if isempty(iters) || ~(iters(end).valid && iters(end).pt.converged && iters(end).pt.crit > 0)
    iters(end+1) = struct('B', Bhi, 'side', 1, 'valid', true, 'pt', ptPM);
end

out = struct('iters', iters, 'Bc', Bc, 'T', T, 'bracket', [Blo Bhi], 'method', method);

% ================================ nested helper ===================================
    function [side, valid, pt] = sample(B)
    %SAMPLE Solve at field B (warm-started), classify, record into iters, update seed.
        o = opts;
        if ~isempty(seed), o.Sigma_seed = seed; end       % T fixed => length always fits
        try
            pt = invzt_solve_point(ion, T, B, lat, o);
        catch err
            if startsWith(err.identifier, 'invzt:')        % structural machinery error
                rethrow(err);                              % (config/lattice/mode) -> surface
            end
            % physics/numeric signal from the shared engine (e.g. invz:degenerateDoublet,
            % a near-singular linear solve): read as an INVALID (ordered) sample.
            pt = struct('converged', false, 'crit', NaN, 'Sigma0', NaN, ...
                'Sigma', [], 'err', err.identifier);
        end
        valid = isfield(pt, 'converged') && pt.converged && isfinite(pt.crit) && ...
                pt.Sigma0 >= sigfloor;
        if valid && pt.crit > 0
            side = 1;                                      % paramagnet
        else
            side = -1;                                     % ordered / invalid
        end
        iters(end+1) = struct('B', B, 'side', side, 'valid', valid, 'pt', pt);
        if isfield(pt, 'converged') && pt.converged && ~isempty(pt.Sigma) && all(isfinite(pt.Sigma))
            seed = pt.Sigma(:);                            % warm-start the next solve
        end
    end
end
