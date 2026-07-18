# Tensor-branch temperature-cut finder: `invzt_critical_T` + two-regime driver

Date: 2026-07-18. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-18): user asked for fixed-field Tc(B) search "just as
in the projected branch" (their experience: T-cuts work better than field cuts
in many occasions), conditional on the tensor/projected core calculations being
fundamentally different — condition affirmatively established (full 12×12
tensor RPA + min-eig criticality vs scalar cc channel; measured +0.016 K
boundary difference at the identical grid). Approach A (full mirror of the
current projected algorithm, shared helpers) explicitly selected.

## Decisions

1. **Mirror the CURRENT projected algorithm, not naive bisection.** The
   projected `invz_critical_T` was rewritten on 2026-07-09 after a rugged-
   boundary failure (documented in the project memory and its header):
   near the boundary the EMT self-consistency suffers critical slowing down,
   so naive bisection latches onto spurious sign flips. The fix — sample
   `crit` on a T-grid, let only converged/finite samples vote, take the
   highest-T ordered→para sign change, refine by regula-falsi — is
   transplanted wholesale with tensor classifiers.
2. **Shared helpers** (Approach A): `invz_refine_crossing.m` moves to
   `invz_common/` (pure `git mv`, same precedent as `invz_peak_energy`);
   the driver's local `invzt_crit_at` is promoted to a module file
   `invz_tensor/invzt_crit_at.m` shared by the finder, the driver, and the
   proxy (exactly the role the projected `invz_crit_at.m` plays).
3. **Driver gains the projected two-regime interface**: `Ts` (fixed-T field
   cuts, low-T branch) + `Bs` (fixed-B temperature cuts, near-vertical
   branch), one flat two-kind `parfor`, merged T-sorted `phase_boundary`.
4. **A unified cross-model driver stays declined** (as in the drivers spec):
   the tensor side still has no ordered/FM solve. This work removes one of
   the two blockers noted there (the T-cut finder); revisit unification when
   FM lands.

## Verified facts this design rests on

- **The projected algorithm** (read in full from
  `invz_projected/invz_critical_T.m`, 108 lines): grid over `[Tlo Thi]`
  (`gridstep` 1/30 K), converged-only voters, highest-T ordered→para
  crossing, `invz_refine_crossing` refinement (`tol` 0.005 K), adaptive
  window (top `Tc0 + 0.05`, `width` 0.5 K, slides ≤ 8×), explicit
  `opts.window` override, `invz:multipleCrossings` warning on re-entrance
  (returns highest-T crossing), `invz:bracket` when nothing brackets.
  With `opts.odd` on and no `opts.Tc0` it errors `invz:oddTc0` (the
  fallback anchor would be wrong) — the pattern the tensor version
  generalizes, since the tensor branch has NO zero-field closed form at all.
- **`invz_refine_crossing` is generic** (closure `f(x) -> [value, ok]`,
  skips non-converged interior samples, falls back to linear interpolation)
  and has exactly two caller files, both in `invz_projected/`
  (`invz_critical.m:40,51`, `invz_critical_T.m:67`). Every test/driver
  reaching those callers already addpaths `invz_common` (established during
  the `invz_peak_energy` move — same call graph). No name collision in
  `invz_common/`.
- **The tensor A1 solver converges metastable PM fixed points inside the
  ordered phase** (LOCKED convention 7; Task-6 finding): valid samples with
  `crit < 0` exist below the boundary — so the converged-only voting grid
  gets votes on BOTH sides of the crossing. This is the opposite of the
  projected ODD-on situation (no metastable window → `invz_critical_T`
  cannot bracket there, ODD-LOG T2), so the tensor T-cut is expected to work
  with `odd` on — the smoke verifies this empirically.
- **No `Sigma_seed` warm-start across a T-scan**: the Matsubara vector
  length is T-dependent (`invzt_critical` seeds only because "T fixed =>
  length always fits", its header). A caller-supplied `Sigma_seed` must be
  stripped before the T-grid solves.
- **The tensor sample-validity rule** lives in the driver's local
  `invzt_crit_at` today: `ok = converged && isfinite(crit) &&
  Sigma0 >= getf(opts,'sigma_floor',-0.5)` — validity-only, deliberately no
  `crit > 0` term (each consumer applies its own phase logic: the proxy
  filters PM points itself; the T-cut finder votes by `sign(crit)`).
  Selective catch: `{invz:degenerateDoublet, invz:orderedPhase,
  invzt:a1ZeroField}` → `ok = false`; all else rethrows.
- **Byte-parity invariant**: the drivers spec keeps its driver code blocks
  byte-identical to the committed files. This spec's Component 4 becomes the
  governing copy of `invzt_run_phase_diagram.m`; the drivers spec's
  Component-1 block gets a one-line SUPERSEDED pointer.

## Components

### 1. `git mv invz_projected/invz_refine_crossing.m invz_common/invz_refine_crossing.m`

Zero content change. Gate: PROJECTED 143/0/19 (its two callers are exercised
by the critical-point tests), CORE unchanged, INTEROP unchanged.

### 2. Promote `invzt_crit_at` to `invz_tensor/invzt_crit_at.m`

The driver's local function moves verbatim into a module file with a proper
header; the driver's local-function section is deleted (call sites unchanged
— same name). Header must document: (a) the tensor analogue-of-`invz_crit_at`
role; (b) `ok` is VALIDITY-only (no `crit > 0`) and why — consumers apply
their own phase logic (`invzt_tc_pm_extrap` filters PM points itself;
`invzt_critical_T` votes by `sign(c)`); (c) the selective-catch contract
(physics signals absorbed, all else rethrows); (d) the single-sourced
`sigma_floor`.

```matlab
function [c, ok] = invzt_crit_at(ion, T, B, lat, opts)
%INVZT_CRIT_AT One tensor criticality sample: crit at (T, B) + sample VALIDITY.
%   [c, ok] = INVZT_CRIT_AT(ion, T, B, lat, opts) solves one A1 tensor point
%   (INVZT_SOLVE_POINT, opts forwarded verbatim) and returns the criticality
%   scalar c = pt.crit plus ok, the three-part sample-validity verdict
%   (converged, finite crit, Sigma0 >= getf(opts,'sigma_floor',-0.5) -- the
%   floor single-sourced with INVZT_CRITICAL, rejecting the spurious
%   negative-Sigma fixed point).
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
        otherwise
            rethrow(err);              % misconfiguration: surface it
    end
end
end
```

### 3. `invz_tensor/invzt_critical_T.m` (new)

```matlab
function [tc, out] = invzt_critical_T(ion, B, lat, opts)
%INVZT_CRITICAL_T Critical temperature at fixed field B (PM side, tensor A1).
%   [tc, out] = INVZT_CRITICAL_T(ion, B, lat, opts) locates the ordered/
%   paramagnet boundary in TEMPERATURE at fixed field -- the fixed-B cut that
%   crosses the near-vertical part of the boundary (T near Tc0) transversally,
%   where INVZT_CRITICAL's fixed-T field cut is ill-conditioned. B is a scalar
%   (transverse-along-a) or [Bx By Bz] (T); lat is the LOCKED lattice struct
%   from INVZT_JQ_TENSOR.
%
%   ALGORITHM (transplant of the projected INVZ_CRITICAL_T, 2026-07-09 rugged-
%   boundary fix): sample crit on a T-grid across the window, let ONLY VALID
%   samples vote (INVZT_CRIT_AT's three-part rule -- converged, finite,
%   Sigma0 above the sigma_floor; near the boundary the outer loop suffers
%   critical slowing down and an invalid sample must get NO vote, or the
%   classifier latches onto spurious sign flips), take the HIGHEST-T
%   ordered(-) -> para(+) sign change, and refine by regula-falsi
%   (INVZ_REFINE_CROSSING, shared with the projected finders via invz_common).
%   NB unlike the projected solver, the tensor A1 map CONVERGES metastable PM
%   fixed points inside the ordered phase (crit < 0, valid) -- so votes exist
%   on both sides of the crossing even with odd on, where the projected T-cut
%   cannot bracket at all (ODD-LOG T2).
%
%   WINDOW: opts.window = [Tlo Thi] (K) explicit; otherwise adaptive -- top
%   anchored at opts.Tc0 + 0.05 K, spanning opts.width down, slid up/down
%   (<= 8x) until it brackets. The tensor branch has NO zero-field closed
%   form to fall back on (mode 'a1' forbids B = 0, invzt:a1ZeroField), so the
%   adaptive path REQUIRES a finite opts.Tc0 -- the driver computes the
%   small-Bx proxy once (INVZT_TC_PM_EXTRAP) and passes it (mirrors the
%   projected invz:oddTc0 rule). Errors invzt:tcAnchor otherwise.
%
%   NO WARM START: any caller-supplied opts.Sigma_seed is STRIPPED before the
%   T-grid solves -- the Matsubara vector length is T-dependent, so a seed
%   from one temperature does not fit another (INVZT_CRITICAL may seed only
%   because its T is fixed).
%
%   RE-ENTRANCE: more than one valid sign change warns invzt:multipleCrossings
%   and returns the highest-T crossing (candidate hyperfine re-entrance is
%   physically reported at low field -- report, never mask).
%
%   OPTIONS (getf defaults; every other field forwards to INVZT_SOLVE_POINT):
%     window    []     explicit [Tlo Thi] (K); validated if given.
%     Tc0       --     zero-field Tc anchor (K) for the adaptive window;
%                      REQUIRED when window is absent.
%     width     0.5    adaptive-window width (K).
%     gridstep  1/30   coarse-grid step (K).
%     tol       0.005  crossing refinement tolerance (K).
%
%   out: .Tg/.c/.ok (last window's samples), .window (final [Tlo Thi]),
%   .ncross (valid sign changes found), .B (validated field row).
%
%   ERRORS invzt:tcWindow (malformed opts.window), invzt:tcAnchor (no window
%   and no finite Tc0), invzt:bracket (no valid crossing after sliding).
%
%   See also INVZT_CRITICAL (fixed-T field cut), INVZT_CRIT_AT,
%   INVZT_TC_PM_EXTRAP, INVZ_REFINE_CROSSING, INVZ_CRITICAL_T (projected
%   reference whose algorithm this transplants).
if nargin < 4, opts = struct(); end
B     = invz_field_vec(B);
width = getf(opts, 'width',    0.5);
gstep = getf(opts, 'gridstep', 1/30);
tol   = getf(opts, 'tol',      0.005);
Tmin  = 0.02;                                   % single-ion solve floor
if isfield(opts, 'Sigma_seed')                  % no warm start across T (see header)
    opts = rmfield(opts, 'Sigma_seed');
end
f = @(T) invzt_crit_at(ion, T, B, lat, opts);

if isfield(opts, 'window') && ~isempty(opts.window)
    win = opts.window;
    if ~(isnumeric(win) && isreal(win) && numel(win) == 2 && all(isfinite(win)) ...
            && win(2) > win(1) && win(1) > 0)
        error('invzt:tcWindow', ...
            'opts.window must be finite [Tlo Thi] with 0 < Tlo < Thi; got %s.', mat2str(win));
    end
    Tlo = win(1);  Thi = win(2);
else
    if ~isfield(opts, 'Tc0') || isempty(opts.Tc0) || ~(isnumeric(opts.Tc0) && isfinite(opts.Tc0))
        error('invzt:tcAnchor', ['invzt_critical_T needs opts.window = [Tlo Thi] or a ' ...
            'finite opts.Tc0 anchor (e.g. the INVZT_TC_PM_EXTRAP small-Bx proxy): the ' ...
            'tensor branch has no zero-field closed form to fall back on.']);
    end
    Thi = opts.Tc0 + 0.05;  Tlo = Thi - width;
end

for slide = 0:8
    Tlo = max(Tlo, Tmin);
    ng  = max(5, round((Thi - Tlo)/gstep) + 1);
    Tg  = linspace(Tlo, Thi, ng);
    c   = nan(1, ng);  ok = false(1, ng);
    for i = 1:ng
        [c(i), ok(i)] = f(Tg(i));
    end
    Tv = Tg(ok);  cv = c(ok);                   % valid samples: the voters
    if numel(cv) >= 2
        sc = find(diff(sign(cv)) ~= 0);         % ordered(-) <-> para(+) transitions
        if numel(sc) > 1
            warning('invzt:multipleCrossings', ...
                ['|B| = %.3f T: %d valid sign changes in [%.3f, %.3f] K ' ...
                 '(possible re-entrance); returning the highest-T crossing.'], ...
                norm(B), numel(sc), Tlo, Thi);
        end
        up = sc(sign(cv(sc)) < 0 & sign(cv(sc+1)) > 0);   % low-T ordered -> high-T para
        if ~isempty(up)
            k  = up(end);                       % highest-T ordered->para crossing
            tc = invz_refine_crossing(f, Tv(k), cv(k), Tv(k+1), cv(k+1), tol);
            out = struct('Tg', Tg, 'c', c, 'ok', ok, 'window', [Tlo Thi], ...
                         'ncross', numel(sc), 'B', B);
            return;
        end
    end
    % No valid crossing in this window: slide toward where it must be.
    % (Check isempty first: all([] > 0) is true in MATLAB.)
    if isempty(cv)                              % nothing valid: keep top, grow down
        Tlo = Tlo - width;
    elseif all(cv > 0)                          % window all paramagnet: Tc below
        Thi = Tlo;  Tlo = Tlo - width;
    elseif all(cv < 0)                          % window all ordered: Tc above
        Tlo = Thi;  Thi = Thi + width;
    else
        break;                                  % mixed signs but no ord->para: give up
    end
end
error('invzt:bracket', ...
    '|B| = %.3f T: no valid paramagnet/ordered crossing found (last window [%.3f, %.3f] K).', ...
    norm(B), Tlo, Thi);
end
```

### 4. `invz_tensor/invzt_run_phase_diagram.m` — two-regime rewiring

Full replacement (this block becomes the governing byte-parity copy;
the drivers spec's Component-1 block gets a SUPERSEDED pointer):

```matlab
%INVZT_RUN_PHASE_DIAGRAM  Full-tensor 1/z PM-side phase boundary: Bc(T) + Tc(B).
%
% TWO-REGIME SEARCH, mirroring the projected invz_run_phase_diagram: fixed-T
% FIELD CUTS (invzt_critical, low-T branch, one Bc per Ts entry) plus fixed-B
% TEMPERATURE CUTS (invzt_critical_T, the near-vertical branch approaching
% Tc(0) where a field cut crosses the boundary at a glancing angle and is
% ill-conditioned -- one Tc per Bs entry). Both kinds run in one flat parfor.
%
% SCOPE: PM side only -- there is still no ordered-phase tensor solve, so
% nothing below the boundary is computed. Unlike the projected T-cut (which
% cannot bracket with ODD on, ODD-LOG T2), the tensor A1 solver converges
% metastable PM points inside the ordered phase, so the T-cut brackets with
% odd on/off alike.
%
% ERROR POLICY: the sweep absorbs ONLY per-point invzt:bracket (a genuine
% no-crossing outcome once Brange/Bs/Twin are preflighted below); every other
% identifier that ESCAPES the finders rethrows. The finders' own internal
% sampler additionally classifies shared-engine physics signals as
% invalid/ordered votes (their committed policy, documented in their headers).

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs ------------------------------------------------------------------
Ts     = linspace(0.4, 1.4, 11);       % fixed-T FIELD-CUT grid (K); low-T branch.
Bs     = linspace(0.25, 1.5, 6);       % fixed-B TEMPERATURE-CUT fields (T); the
                                       % near-vertical branch. [] disables T-cuts.
Twin   = [];                           % [] -> adaptive T-cut window anchored at the
                                       % small-Bx proxy Tc0 (computed below); or an
                                       % explicit [Tlo Thi] (K) forwarded to
                                       % invzt_critical_T (skips the anchor).
Brange = [0.05 6.0];                   % field-cut [Blo Bhi] bracket (T). Blo > 0:
                                       % mode 'a1' forbids exact zero transverse field.
gridN  = 16;  gridConv = 'halfopen';   % invzt_qgrid(gridN, gridConv)
dpRng  = 30;                           % invzt_jq_tensor coupling-sum range
useParallel = true;                    % false -> force serial
solve_opts  = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', 'dress', 'full');
                                       % sigma_floor may be added here too; defaults to
                                       % invzt_critical's own -0.5 (single-sourced)

show_proxy_anchor     = true;          % plot the small-Bx Tc(0) proxy marker. The proxy
                                       % is COMPUTED regardless whenever Bs is nonempty
                                       % and Twin is [] (it anchors the adaptive T-cut
                                       % window); this knob only controls the marker.
Bproxy                = 0.05;          % transverse field (T) for that proxy
Ts_proxy              = 1.40:1/30:2.00;% proxy extrapolation grid (K); must contain >= 2
                                       % converged PM points at Bproxy (i.e. reach above
                                       % Tc0 ~ 1.56 K at production grids).
show_projected_anchor = false;         % OPT-IN cross-model comparator: the PROJECTED
                                       % closed-form Tc(0) (invz_odd_zero_field).
                                       % Off by default because it (a) puts
                                       % invz_projected on the path, (b) requires
                                       % ion.demag == 0, (c) solves 7 variants per
                                       % grid, and (d) is a DIFFERENT ODD treatment
                                       % (projected Tier-1 chi_perp-mediated dJ) from
                                       % this driver's structural tensor `odd` flag --
                                       % a comparator, NOT the same quantity.
% -----------------------------------------------------------------------------

% Preflights (invzt:bracket doubles as the finders' arg-validation id; absorbing
% it per point without these checks would turn a typo into a silent all-NaN sweep).
assert(isnumeric(Brange) && isreal(Brange) && numel(Brange) == 2 && ...
    all(isfinite(Brange)) && Brange(2) > Brange(1) && Brange(1) > 0, ...
    'invzt_run_phase_diagram:Brange', ...
    'Brange must be finite [Blo Bhi] with 0 < Blo < Bhi (mode ''a1'' forbids B = 0); got %s.', ...
    mat2str(Brange));
assert(isempty(Bs) || (isnumeric(Bs) && isreal(Bs) && isvector(Bs) && ...
    all(isfinite(Bs)) && all(Bs > 0)), 'invzt_run_phase_diagram:Bs', ...
    'Bs must be empty or a finite positive vector (mode ''a1'' forbids B = 0); got %s.', ...
    mat2str(Bs));
assert(isempty(Twin) || (isnumeric(Twin) && isreal(Twin) && numel(Twin) == 2 && ...
    all(isfinite(Twin)) && Twin(2) > Twin(1) && Twin(1) > 0), ...
    'invzt_run_phase_diagram:Twin', ...
    'Twin must be [] (adaptive) or finite [Tlo Thi] with 0 < Tlo < Thi; got %s.', mat2str(Twin));

if show_projected_anchor
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_projected'));
    assert(ion.demag == 0, 'invzt:anchorDemag', ...
        ['show_projected_anchor requires ion.demag == 0 (invz_odd_blocks is ' ...
         'intrinsic-only). Set ion.demag = 0 or disable the anchor.']);
end

g   = invzt_qgrid(gridN, gridConv);
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

% ---- small-Bx proxy Tc(0): plot marker AND the adaptive T-cut window anchor --
Tc0_proxy = NaN;
need_anchor = ~isempty(Bs) && isempty(Twin);
if show_proxy_anchor || need_anchor
    critfun = @(T) invzt_crit_at(ion, T, [Bproxy 0 0], lat, solve_opts);
    try
        Tc0_proxy = invzt_tc_pm_extrap(critfun, Ts_proxy);
    catch err
        if ~strcmp(err.identifier, 'invzt:tcNoWindow'), rethrow(err); end
        fprintf('  Tc(0) proxy: no PM extrapolation window on Ts_proxy (skipped)\n');
    end
end
if need_anchor
    assert(isfinite(Tc0_proxy), 'invzt_run_phase_diagram:tcAnchor', ...
        ['The T-cut jobs need an adaptive-window anchor but the small-Bx proxy did ' ...
         'not resolve on Ts_proxy. Extend Ts_proxy above Tc0, or set an explicit Twin.']);
end
topts = solve_opts;                       % T-cut finder opts: anchor or window
if ~isempty(Twin), topts.window = Twin; else, topts.Tc0 = Tc0_proxy; end

% ---- one flat parfor over both cut kinds -------------------------------------
nWorkers = 0;
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

nT = numel(Ts);  nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];              % one independent 1-D root find per job
kind = [ones(1, nT) 2*ones(1, nB)];    % 1 = Bc(T) at fixed T, 2 = Tc(B) at fixed B
res  = nan(1, nT + nB);
parfor (k = 1:nT+nB, nWorkers)
    t0 = tic;  v = jobs(k);  val = NaN;
    % val-then-assign: keeps the sliced output unconditionally assigned.
    try
        if kind(k) == 1
            val = invzt_critical(ion, v, lat, Brange, solve_opts);
        else
            val = invzt_critical_T(ion, v, lat, topts);
        end
    catch err
        % invzt:bracket is a genuine per-point outcome (inputs preflighted
        % above, so this window/bracket simply lacks a crossing). Anything
        % else is a MISCONFIGURATION and rethrows.
        if ~strcmp(err.identifier, 'invzt:bracket'), rethrow(err); end
        if kind(k) == 1
            fprintf('  T = %.2f K: no bracket in [%.2f %.2f] T (left NaN)\n', ...
                    v, Brange(1), Brange(2));
        else
            fprintf('  B = %.2f T: no valid T-crossing (left NaN)\n', v);
        end
    end
    if kind(k) == 1
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', ...
                k, nT+nB, v, val, toc(t0));
    else
        fprintf('  [%2d/%2d] B = %.2f T  ->  Tc = %.3f K   (%.0f s)\n', ...
                k, nT+nB, v, val, toc(t0));
    end
    res(k) = val;
end
Bc  = res(1:nT);                       % low-T branch:  Bc at each Ts
TcB = res(nT+1:end);                   % high-T branch: Tc at each Bs

% ---- optional projected closed-form Tc(0) comparator --------------------------
Tc0_closed = NaN;
if show_projected_anchor
    % Track this driver's own odd setting -- comparing a tensor ODD-on curve
    % against a projected ODD-off anchor (or vice versa) would be misleading.
    anchor_mode = 'off';
    if ~isfield(solve_opts, 'odd') || ~isequal(solve_opts.odd, false), anchor_mode = 'full'; end
    % Same NOMINAL N and dipole cutoff as this driver's lat -- but NB the
    % projected engine always builds its own legacy endpoint-inclusive
    % [-0.5, 0.5] mesh (a different quadrature convention from 'halfopen'),
    % on a different model with a different ODD treatment: a cross-model
    % COMPARATOR, never the same quantity.
    Tc0_closed = invz_odd_zero_field(ion, struct('mode', anchor_mode, ...
        'grids', {{gridN}}, 'dpRng', dpRng, 'cache', true));
end

figure; hold on;
plot(Bc, Ts, 'o-', 'DisplayName', 'tensor A1: field-cut B_c(T)');
if nB > 0
    plot(Bs, TcB, 's-', 'DisplayName', 'tensor A1: T-cut T_c(B)');
end
if show_proxy_anchor && isfinite(Tc0_proxy)
    plot(0, Tc0_proxy, 'k^', 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('tensor small-B_x proxy T_c(0), %.2f T', Bproxy));
end
if isfinite(Tc0_closed)
    plot(0, Tc0_closed, 'ks', 'MarkerFaceColor', 'c', ...
         'DisplayName', 'projected closed-form T_c(0) (comparator; legacy-inclusive mesh)');
end
xlabel('B_c (T)'); ylabel('T (K)');
title('LiHoF_4 full-tensor 1/z phase boundary (PM side: field cuts + T cuts)');
legend('Location', 'southwest');

% Merged boundary, T-sorted, finite points only. Columns [T(K) B(T)]. Near the
% regime join both branches can contribute a point, so the curve is not
% strictly single-valued in T there (same note as the projected driver).
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);
```

Structural notes: the local `invzt_crit_at` function at the bottom of the
current driver is DELETED (Component 2 provides it as a module file; the
`critfun` closure line is unchanged). The `Ts` default drops its old
`> 1.5 K` tail (those points were guaranteed no-bracket NaNs; the T-cut
branch now covers that region properly).

### 5. Tests — `invz_tensor/tests/test_invzt_critical_T.m` (new)

Standard CORE `setupOnce` boilerplate (invz_tensor, repo root, invz_common).
Three test functions:

- `test_anchor_required` (fast, no solves — the guard fires before any
  compute; `lat` may be a dummy struct):
  `verifyError(@() invzt_critical_T(ion, 1.0, struct(), struct()), 'invzt:tcAnchor')`.
- `test_window_validated` (fast): reversed window
  `struct('window', [1.8 1.0])` → `verifyError(..., 'invzt:tcWindow')`.
- `test_tcut_matches_field_cut_slow` (`INVZ_SLOW`-gated, the projected
  crossing-consistency pattern): at 8³/`dpRng` 10–15, `odd` false for speed,
  T* = 1.4 K: `Bstar = invzt_critical(ion, 1.4, lat, [0.05 6], o)`;
  `tc = invzt_critical_T(ion, Bstar, lat, setfield(o, 'window', [1.0 1.8]))`
  must return 1.4 within 0.05 K. Validates the mirror against the
  already-validated field finder with no new reference data.

CORE becomes **52 passed / 0 failed / 2 incomplete** (50 + the two fast
tests; the slow one joins the A4 ladder test as filtered-incomplete).

### 6. Driver smoke (verification, not a committed test)

Same-directory `sed` copy: `Ts` → `[1.2 1.4]`, `Bs` → `[0.5 1.5]`,
`gridN` → 8, `dpRng` → 15, `useParallel` → false,
`show_projected_anchor` → true (exercises the comparator once).
`Ts_proxy` stays (8³ Tc0 ≈ 1.53 per the drivers-plan smoke — inside the grid).
Assertions: ≥ 1 finite `Bc`; **both** `TcB` finite (this is the empirical
proof the tensor T-cut brackets with `odd = true`, which the projected branch
cannot do); `TcB(1) > TcB(2)` (Tc decreases with field — physical);
`phase_boundary` has ≥ 3 rows; a figure exists. Delete the smoke copy after.
Expected T-cut values: both lie between 1.40 K and the proxy Tc0 ≈ 1.53 K at
this coarse grid (the 8³ boundary has Bc(1.4 K) = 1.916 T > 1.5 T, so both
fields are ordered at 1.4 K and PM at Tc0) — Tc(0.5 T) near ≈ 1.50–1.53 K,
Tc(1.5 T) near ≈ 1.42–1.50 K. Magnitudes reported, only the ordering
`TcB(0.5 T) > TcB(1.5 T)` gated.

## Error handling

- Preflights before compute: `Brange`/`Bs`/`Twin` at the driver
  (`invzt_run_phase_diagram:Brange/:Bs/:Twin`), window/anchor at the finder
  (`invzt:tcWindow`/`invzt:tcAnchor`), proxy-anchor availability
  (`invzt_run_phase_diagram:tcAnchor`).
- Per-point: only `invzt:bracket` absorbed (both cut kinds), all else
  rethrows; the sampler's selective-catch policy is unchanged (Component 2).
- `invzt:multipleCrossings` is a warning (re-entrance is physics to report,
  never mask); `invzt:tcNoWindow` from the proxy is caught, the anchor
  assert then decides whether it matters.

## Testing / verification

1. Component 1 first: `git mv`, then PROJECTED (143/0/19) + CORE (50/0/1)
   + INTEROP (8/0/2) — all unchanged.
2. Components 2–3 + 5: fast tests RED→GREEN where meaningful (the two error
   guards are new behavior — write test first, watch it fail for the right
   reason: `invzt_critical_T` undefined), then CORE 52/0/2.
3. Slow gate: `INVZ_SLOW=1` run of the new test file only (minutes at 8³).
4. Component 4 + smoke (Component 6). Suites re-run at the end:
   CORE 52/0/2, INTEROP 8/0/2, PROJECTED 143/0/19.
5. Docs: README §2 callout (T-cut no longer missing — reword the scope
   sentence: still no ordered/FM solve; the near-vertical region is now
   covered by fixed-B cuts), §2.1 recipe note, module map rows for
   `invzt_critical_T` + `invzt_crit_at`, architecture table A1 row mention;
   drivers spec Component-1 SUPERSEDED pointer; drivers spec/plan
   quickstart counts stay (they cite their own point-in-time gates —
   only the README's living quickstart count updates to 52/0/2).

## Success criteria

- `invzt_critical_T` returns Tc(B) on the tensor branch with `odd` on and
  off; the slow crossing-consistency test round-trips Bc↔Tc within 0.05 K.
- The driver produces a merged two-branch boundary in one run; smoke shows
  both T-cut points finite and ordered correctly.
- `invz_refine_crossing` lives in `invz_common/`; all three suites at their
  expected counts (CORE 52/0/2, INTEROP 8/0/2, PROJECTED 143/0/19).
- README §2 and the module map reflect the new capability.

## Out of scope

- Ordered-phase tensor solve (still the FM blocker to driver unification).
- The projected `invz_critical_T`'s smoothness regression
  (`test_tc_boundary_is_smooth`) — a production-sweep-scale test; the tensor
  branch defers it until someone runs a production T-cut sweep.
- Adaptive per-field window *learning* (each T-cut is independent; the
  projected driver's per-field self-adaptation via slide is inherited as-is).
- Warm-starting across T (documented impossibility with T-dependent
  Matsubara grids; revisit only with an interpolating seed scheme).
