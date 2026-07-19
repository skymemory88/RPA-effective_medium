function [tc, out] = invzt_critical_T(ion, B, lat, opts)
%INVZT_CRITICAL_T Critical temperature at fixed field B (PM side, tensor A1).
%   [tc, out] = INVZT_CRITICAL_T(ion, B, lat, opts) locates the ordered/
%   paramagnet boundary in TEMPERATURE at fixed field -- the fixed-B cut that
%   crosses the near-vertical part of the boundary (T near Tc0) transversally,
%   where INVZT_CRITICAL's fixed-T field cut is ill-conditioned. B is a scalar
%   (transverse-along-a) or [Bx By Bz] (T); lat is the LOCKED lattice struct
%   from INVZT_JQ_TENSOR.
%
%   ALGORITHM (transplant of the projected INVZ_CRITICAL_T's 2026-07-09
%   rugged-boundary fix, with the inherited classifier gaps closed): sample
%   crit on a T-grid across the window, let ONLY VALID samples vote
%   (INVZT_CRIT_AT's three-part rule -- converged, finite, Sigma0 above the
%   sigma_floor; near the boundary the outer loop suffers critical slowing
%   down and an invalid sample must get NO vote, or the classifier latches
%   onto spurious sign flips), then act on INVZT_TC_PICK's pure decision: an
%   exactly-critical sample (crit == 0) is itself the boundary and is
%   returned directly; a bracketing ordered->para voter pair is refined by
%   regula-falsi (INVZ_REFINE_CROSSING, shared with the projected finders via
%   invz_common); otherwise the ADAPTIVE window moves toward the boundary.
%   NB unlike the projected solver, the tensor A1 map CONVERGES metastable PM
%   fixed points inside the ordered phase (crit < 0, valid) -- so votes exist
%   on both sides of the crossing even with odd on, where the projected T-cut
%   cannot bracket at all (ODD-LOG T2).
%
%   WINDOW: opts.window = [Tlo Thi] (K) is a HARD bound -- exactly one grid
%   pass, no sliding; if it contains no returnable root the function errors
%   invzt:bracket (widen the window). Without opts.window the ADAPTIVE mode
%   anchors the top at opts.Tc0 + 0.05 K, spans opts.width down, and moves
%   (up to 9 window attempts total) following INVZT_TC_PICK: top voter
%   ordered -> move up (this includes the lower para->ordered leg of a
%   re-entrant region -- the physical high-T side is paramagnetic, so the
%   requested highest-T crossing lies above); all voters PM -> move down; no
%   voters at all -> grow down keeping the top. The tensor branch has NO
%   zero-field closed form to fall back on (mode 'a1' forbids B = 0,
%   invzt:a1ZeroField), so adaptive mode REQUIRES a finite positive scalar
%   opts.Tc0 -- the driver computes the small-Bx proxy once
%   (INVZT_TC_PM_EXTRAP) and passes it (mirrors the projected invz:oddTc0
%   rule). Errors invzt:tcAnchor otherwise.
%
%   SEEDING (execution finding E1): any caller-supplied opts.Sigma_seed is
%   STRIPPED (the Matsubara vector length is T-dependent, so a caller's seed
%   fits at most ONE of the sampled temperatures) -- but the finder threads
%   its OWN rolling seed: the grid is sampled DESCENDING in T (paramagnet
%   first) and each solve is seeded with the previous valid Sigma
%   re-interpolated onto the new T's Matsubara grid (nested sampler below).
%   Without branch tracking, cold starts hop between A1 fixed-point branches
%   and the sampled crit(T) flickers sign (measured ncross = 5 and a
%   +0.20 K Tc error at the 8^3 round-trip gate) -- the same reason
%   INVZT_CRITICAL warm-starts its field bisection (its header: the seed
%   "helps the ordered-side solves land on the metastable PM branch").
%
%   RE-ENTRANCE: boundary indicators (INVZT_TC_PICK's ncross: strict sign
%   flips + exactly-critical runs) are ACCUMULATED across adaptive attempts
%   (within-window indicators never span windows: crossing pairs are
%   intra-window, and an exact zero at a window top returns immediately so
%   it is never re-counted) -- more than one across the whole scan warns
%   invzt:multipleCrossings and the highest-T root is returned (candidate
%   hyperfine re-entrance is physically reported at low field -- report,
%   never mask; second review R5: a down-leg seen in one window and the
%   upper crossing in a later window must still warn).
%
%   OPTIONS (getf defaults; every other field forwards to INVZT_SOLVE_POINT):
%     window    []     HARD [Tlo Thi] bound (K); validated (degenerate spans
%                      Thi - Tlo <= 1e-9 rejected); no sliding.
%     Tc0       --     zero-field Tc anchor (K) for the adaptive window;
%                      REQUIRED (finite positive scalar > the solve floor)
%                      when window is absent.
%     width     0.5    adaptive-window width (K); finite positive real scalar.
%     gridstep  1/30   coarse-grid step (K); finite positive real scalar.
%     tol       0.005  crossing refinement tolerance (K); finite positive
%                      real scalar.
%
%   out: .Tg/.c/.ok (the last SAMPLED window's samples, incl. invalid ones),
%   .window (the last SAMPLED [Tlo Thi] -- never a slide-mutated pair, second
%   review R2), .ncross (boundary indicators accumulated across the scan),
%   .B (validated field row).
%
%   ERRORS invzt:tcOpts (malformed width/gridstep/tol), invzt:tcWindow
%   (malformed or below-floor opts.window), invzt:tcAnchor (adaptive mode
%   without a usable Tc0), invzt:bracket (no returnable root: hard window
%   without a crossing, adaptive attempts exhausted, the window collapsed at
%   the solve floor, or a no-valid-voter window already at the floor --
%   nothing below left to expose, second review R3). Bracket messages always
%   report the last SAMPLED window: winS is captured only AFTER the collapse
%   guard confirms a grid will be evaluated, so a floor collapse reports the
%   previous attempt's really-sampled window (third review T2).
%
%   See also INVZT_CRITICAL (fixed-T field cut), INVZT_TC_PICK,
%   INVZT_CRIT_AT, INVZT_TC_PM_EXTRAP, INVZ_REFINE_CROSSING,
%   INVZ_CRITICAL_T (projected reference whose algorithm this transplants).
if nargin < 4, opts = struct(); end
B     = invz_field_vec(B);
width = getf(opts, 'width',    0.5);
gstep = getf(opts, 'gridstep', 1/30);
tol   = getf(opts, 'tol',      0.005);
Tmin  = 0.02;                                   % single-ion solve floor
posscal = @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x > 0;
if ~(posscal(width) && posscal(gstep) && posscal(tol))
    error('invzt:tcOpts', ['width, gridstep and tol must be finite positive ' ...
        'real scalars; got width = %s, gridstep = %s, tol = %s.'], ...
        invzt_str(width), invzt_str(gstep), invzt_str(tol));
end
if isfield(opts, 'Sigma_seed')                  % caller seeds fit ONE T only (see header)
    opts = rmfield(opts, 'Sigma_seed');
end
Ecut_used = getf(opts, 'Ecut', 40);             % must match the solver's Matsubara grid
wn_prev = [];  Sig_prev = [];                   % rolling branch-tracking seed state
f = @(T) sample(T);

hardwin = isfield(opts, 'window') && ~isempty(opts.window);
if hardwin
    win = opts.window;
    if ~(isnumeric(win) && isreal(win) && numel(win) == 2 && all(isfinite(win(:))) ...
            && win(2) - win(1) > 1e-9 && win(1) > 0)
        error('invzt:tcWindow', ...
            'opts.window must be finite real [Tlo Thi] with 0 < Tlo < Thi (span > 1e-9); got %s.', ...
            invzt_str(win));
    end
    if win(2) <= Tmin
        error('invzt:tcWindow', ...
            'opts.window = %s lies entirely at/below the %.3g K solve floor.', ...
            invzt_str(win), Tmin);
    end
    Tlo = win(1);  Thi = win(2);
else
    Tc0ok = isfield(opts, 'Tc0') && ~isempty(opts.Tc0) && posscal(opts.Tc0) ...
            && opts.Tc0 > Tmin;
    if ~Tc0ok
        error('invzt:tcAnchor', ['invzt_critical_T needs opts.window = [Tlo Thi] ' ...
            'or a finite positive scalar opts.Tc0 anchor > %.3g K (e.g. the ' ...
            'INVZT_TC_PM_EXTRAP small-Bx proxy): the tensor branch has no ' ...
            'zero-field closed form to fall back on.'], Tmin);
    end
    Thi = opts.Tc0 + 0.05;  Tlo = Thi - width;
end

nattempt = 1;  if ~hardwin, nattempt = 9; end   % hard window: ONE pass, no sliding
nseen = 0;                                      % boundary indicators ACROSS attempts (R5)
winS  = [Tlo Thi];                              % last SAMPLED window (reporting, R2)
for attempt = 1:nattempt
    Tlo = max(Tlo, Tmin);
    if Thi <= Tlo + 1e-9, break; end            % collapsed at the floor -> invzt:bracket
                                                % (winS keeps the PREVIOUS attempt's
                                                % really-sampled window -- T2)
    winS = [Tlo Thi];                           % captured only when a grid WILL run
    ng  = max(5, round((Thi - Tlo)/gstep) + 1);
    Tg  = linspace(Tlo, Thi, ng);
    c   = nan(1, ng);  ok = false(1, ng);
    for i = ng:-1:1                             % DESCENDING T: the first (cold)
        [c(i), ok(i)] = f(Tg(i));               % solve sits on the PM-most point;
    end                                         % lower samples ride the rolling seed
    Tv = Tg(ok);  cv = c(ok);                   % valid samples: the voters
    if isempty(cv)
        act = 'grow';  ncross = 0;              % nothing valid: keep top, grow down
    else
        [act, ka, kb, ncross] = invzt_tc_pick(cv);
    end
    % Re-entrance evidence accumulates ACROSS attempts (R5). No double count:
    % crossing pairs are intra-window, and an exact zero at a window top
    % returns immediately (so no zero is ever re-counted by a later window).
    nseen = nseen + ncross;
    switch act
        case {'zero', 'bracket'}
            if nseen > 1
                warning('invzt:multipleCrossings', ...
                    ['|B| = %.3f T: %d boundary indicators seen across the scan ' ...
                     '(possible re-entrance); returning the highest-T root of ' ...
                     'window [%.3f, %.3f] K.'], norm(B), nseen, winS(1), winS(2));
            end
            if strcmp(act, 'zero')
                tc = Tv(ka);                    % exactly-critical sample IS the boundary
            else
                tc = invz_refine_crossing(f, Tv(ka), cv(ka), Tv(kb), cv(kb), tol);
            end
            out = struct('Tg', Tg, 'c', c, 'ok', ok, 'window', winS, ...
                         'ncross', nseen, 'B', B);
            return;
    end
    if hardwin, break; end                      % ONE pass (R2): leave BEFORE any slide,
                                                % so the error reports the window the
                                                % user actually asked for and we sampled
    switch act
        case 'up'                               % top voter ordered: boundary above
            Tlo = Thi;  Thi = Thi + width;
        case 'down'                             % all voters PM: boundary below
            Thi = Tlo;  Tlo = Tlo - width;
        case 'grow'
            if winS(1) <= Tmin + 1e-9, break; end   % already at the floor: nothing
                                                    % below left to expose (R3 -- do
                                                    % not re-sample the identical grid)
            Tlo = Tlo - width;
    end
end
if hardwin
    error('invzt:bracket', ...
        ['|B| = %.3f T: no valid ordered/paramagnet root in the EXPLICIT window ' ...
         '[%.3f, %.3f] K (hard bound, no sliding): widen opts.window.'], ...
        norm(B), winS(1), winS(2));
else
    error('invzt:bracket', ...
        '|B| = %.3f T: no valid ordered/paramagnet root found (last sampled window [%.3f, %.3f] K).', ...
        norm(B), winS(1), winS(2));
end

% ------------------------------ nested sampler --------------------------------
    function [c1, ok1] = sample(T)
    %SAMPLE Branch-tracked criticality sample: INVZT_CRIT_AT + a rolling seed.
    %   Seeds each solve with the previous VALID Sigma re-interpolated onto
    %   this T's Matsubara grid (linear in the frequency values; the short
    %   tail beyond the previous grid's top is clamped to its last value --
    %   Sigma decays there). The state advances only on valid samples, so an
    %   absorbed/invalid solve never poisons the seed chain.
        oi = opts;
        if ~isempty(Sig_prev)
            wn_i = invz_matsubara(T, Ecut_used);
            sd = interp1(wn_prev, Sig_prev, wn_i, 'linear');
            sd(wn_i > wn_prev(end)) = Sig_prev(end);
            oi.Sigma_seed = sd(:);
        end
        [c1, ok1, pt] = invzt_crit_at(ion, T, B, lat, oi);
        if ok1 && isfield(pt, 'Sigma') && ~isempty(pt.Sigma) && all(isfinite(pt.Sigma))
            wn_prev  = invz_matsubara(T, Ecut_used);
            Sig_prev = pt.Sigma(:);
        end
    end
end
