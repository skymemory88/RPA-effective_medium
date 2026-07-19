# Tensor-branch temperature-cut finder: `invzt_critical_T` + two-regime driver

Date: 2026-07-18. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-18): user asked for fixed-field Tc(B) search "just as
in the projected branch" (their experience: T-cuts work better than field cuts
in many occasions), conditional on the tensor/projected core calculations being
fundamentally different — condition affirmatively established (full 12×12
tensor RPA + min-eig criticality vs scalar cc channel; measured +0.016 K
boundary difference at the identical grid). Approach A (full mirror of the
current projected algorithm, shared helpers) explicitly selected.
**REV 2 (2026-07-18)**, after the pre-execution external review
`critical-t-review_by_Codex.md` (findings F1–F8, all verified against source
and MATLAB R2025a behavior before acceptance) — dispositions in the
"Review dispositions" section; the Component code blocks below are the
post-review, governing versions.
**REV 3 (2026-07-19)**, after the second-pass review (same file, R1–R6, all
verified before acceptance): the seed-strip test was a FALSE PROOF (the
solver silently ignores wrong-length seeds — now a poisonous length-matched
NaN seed asserting the endpoint stays valid); hard windows no longer mutate
their reported window (break before the slide switch); the no-valid `grow`
path terminates at the solve floor instead of re-sampling the identical
grid; the driver preflights `Btol`/`Ttol`/`Twidth`/`Tgridstep` and Twin's
floor; re-entrance evidence accumulates ACROSS adaptive attempts (`nseen`);
`sigma_floor` gets direct threading coverage.
**REV 4 (2026-07-19)**, after the third-pass review (T1–T3): the hard-window
message assertion is rewritten as an explicit try/catch (`verifyError` does
NOT return the caught `MException` — reviewer-reproduced `<missing>`); `winS`
is captured only AFTER the floor-collapse guard so a collapse reports the
previous really-sampled window; degenerate hard-window spans are rejected;
the hard-window diagnostic and the adaptive floor termination gain CHEAP
fast tests via the B = 0 all-invalid trick (pre-lattice `invzt:a1ZeroField`
makes every sample invalid — no physics solves); wording/id-list drift
cleaned (the seed is length-MATCHED and poisonous-VALUED, not
"incompatible").
**REV 5 (2026-07-19, execution finding E1)**: the Task-2 integration gates
caught a real algorithm gap none of the three reviews saw — with all
transcriptions byte-verified and every synthetic test green, the sampled
`crit(T)` at fixed field flipped sign FIVE times (round-trip returned
1.596 K vs the field-cut's 1.4 K; the odd-on gate exceeded the zero-field
Tc). Root cause: cold-started solves hop between A1 fixed-point branches —
the field-cut finder only works because it WARM-STARTS every sample
(`invzt_critical`'s header: the seed "helps the ordered-side solves land on
the metastable PM branch"). Fix, mirroring that precedent: the T-grid is
sampled DESCENDING in T (paramagnet first) with a rolling seed
re-interpolated onto each new T's Matsubara grid (`invzt_crit_at` gains a
third `pt` output to expose the converged `Sigma`); caller seeds are still
stripped.

## Decisions

1. **Mirror the CURRENT projected algorithm's core idea, with its inherited
   classifier gaps closed.** The projected `invz_critical_T` was rewritten on
   2026-07-09 after a rugged-boundary failure: near the boundary the EMT
   self-consistency suffers critical slowing down, so naive bisection latches
   onto spurious sign flips. Transplanted: valid-samples-only T-grid voting,
   highest-T ordered→para crossing, regula-falsi refinement, adaptive
   Tc0-anchored window. Closed here (review F1; the projected code shares
   these defects verbatim — noted as a possible projected-side follow-up, not
   touched in this work): an exactly-sampled root (`crit == 0`) is recognized
   and returned; a window whose top voter is ordered always slides UP (this
   includes the lower para→ordered leg of a re-entrant region, which the
   projected classifier abandons); zeros are never double-counted as two
   crossings.
2. **Pure decision core.** The vote-classification/slide policy lives in a
   pure function `invzt_tc_pick(cv)` (no solves, no state) so every branch is
   millisecond-testable on synthetic votes (review F1/F3).
3. **Explicit window is a HARD bound** (review F4): one grid pass, no
   sliding, `invzt:bracket` if it contains no returnable root. Only the
   adaptive (`Tc0`-anchored) mode slides — up to 9 window attempts.
4. **Shared helpers** (Approach A): `invz_refine_crossing.m` moves to
   `invz_common/` (pure `git mv`, same precedent as `invz_peak_energy`);
   the driver's local `invzt_crit_at` is promoted to a module file shared by
   the finder, the driver, and the proxy.
5. **Driver gains the projected two-regime interface** with **namespaced
   tolerances** (review F5): `Ts` (field cuts) + `Bs` (temperature cuts),
   knobs `Btol` (tesla, → `invzt_critical` `tol`) and `Ttol`/`Twidth`/
   `Tgridstep` (kelvin, → `invzt_critical_T`), merged into per-finder opts
   at the call boundary — never one shared `tol` in two different units.
6. **Strict input validation with safe formatting** (review F2): all
   validators use `invzt_str`, never raw `mat2str` (which throws on structs
   while constructing the intended error). NB `invzt_str`'s CURRENT
   fallthrough is itself raw `mat2str`, so it too throws on structs — it is
   HARDENED in this work (Component 2b: a class/size placeholder branch for
   anything `mat2str` cannot format), with a fast regression assert.
7. **A unified cross-model driver stays declined** (as in the drivers spec):
   the tensor side still has no ordered/FM solve. This work removes one of
   the two blockers noted there (the T-cut finder); revisit unification when
   FM lands.

## Verified facts this design rests on

- **The projected algorithm** (read in full from
  `invz_projected/invz_critical_T.m`, 108 lines): grid over `[Tlo Thi]`
  (`gridstep` 1/30 K), converged-only voters, highest-T ordered→para
  crossing, `invz_refine_crossing` refinement (`tol` 0.005 K), adaptive
  window (top `Tc0 + 0.05`, `width` 0.5 K), explicit `opts.window`
  override (which there still slides — the hard-bound contract here is a
  deliberate improvement), `invz:multipleCrossings` warning,
  `invz:bracket` when nothing brackets. With `opts.odd` on and no
  `opts.Tc0` it errors `invz:oddTc0` — the pattern the tensor version
  generalizes, since the tensor branch has NO zero-field closed form at all.
- **The projected classifier's gaps (review F1, confirmed by
  hand-simulation):** for votes `[-1, 0, 1]`, `diff(sign(cv))` reports two
  changes but the strict `-→+` predicate matches neither pair — an
  exactly-sampled root produces `invz:bracket` plus a spurious
  multiple-crossings warning; for `[+, +, -, -]` (the lower leg of a
  re-entrant region) it breaks instead of sliding up toward the physical
  high-T paramagnet. Both fixed in `invzt_tc_pick`.
- **MATLAB R2025a validation behavior (review F2, confirmed):**
  `mat2str(struct())` throws `MATLAB:mat2str:NumericInput` (so an error
  message built with it masks the intended identifier); a vector on the
  right of `&&` throws `MATLAB:nonLogicalConditional`; `isfinite(1+1i)` is
  `true`. Hence: scalar/real/positive validators and `invzt_str` formatting.
  Also the `struct('f', [a b])` constructor trap: array values create struct
  ARRAYS — tests build malformed-opts fixtures by field assignment. And
  `invzt_str.m` (read in full): its non-char branch is raw `mat2str(x)`, so
  it throws on structs/cells/objects — hardened in Component 2b before any
  validator relies on it.
- **The projected `invz_critical_T` was NOT modified** — it shares F1's
  exact-zero and re-entrant-lower-leg gaps verbatim; flagged as a candidate
  follow-up, out of scope here (validated code, own test surface).
- **`invz_refine_crossing` is generic** (closure `f(x) -> [value, ok]`,
  skips non-converged interior samples via midpoint retry, falls back to
  linear interpolation) and has exactly two caller files, both in
  `invz_projected/` (`invz_critical.m:40,51`, `invz_critical_T.m:67`).
  Every test/driver reaching those callers already addpaths `invz_common`
  (established during the `invz_peak_energy` move — same call graph). No
  name collision in `invz_common/`. Fast direct coverage of the moved
  helper post-move: the fast interop critical-parity test exercises the
  projected field-cut path, plus a NEW direct unit test in CORE (the
  substantive projected T-cut tests are `INVZ_SLOW`-gated — review F8).
- **The tensor A1 solver converges metastable PM fixed points inside the
  ordered phase** (LOCKED convention 7): valid `crit < 0` votes exist below
  the boundary — the voting grid sees both sides of the crossing. This is
  the opposite of the projected ODD-on situation (no metastable window →
  `invz_critical_T` cannot bracket there, ODD-LOG T2), so the tensor T-cut
  works with `odd` on — gated by a committed `INVZ_SLOW` test, not just the
  throwaway smoke (review F3).
- **No `Sigma_seed` warm-start across a T-scan**: the Matsubara vector
  length is T-dependent (`invzt_critical` seeds only because "T fixed =>
  length always fits", its header). A caller-supplied `Sigma_seed` is
  stripped; the slow round-trip test proves the strip with a seed that is
  length-MATCHED to the window's top endpoint but poisonous in VALUE (NaN)
  — a wrong-LENGTH seed is silently ignored by the solver's
  `numel == nwn` guard and proves nothing (third-review wording fix).
- **`invzt_crit_at`'s guard paths are pre-lattice**: `invzt_solve_point`
  validates `opts.mode` and the zero-transverse-field guard BEFORE touching
  `lat`, so sampler-contract tests run in milliseconds with dummy structs.
- **The tensor sample-validity rule**: `ok = converged && isfinite(crit) &&
  Sigma0 >= getf(opts,'sigma_floor',-0.5)` — validity-only, deliberately no
  `crit > 0` term (each consumer applies its own phase logic: the proxy
  filters PM points itself; the T-cut finder votes by `sign(crit)`).
  Selective catch: `{invz:degenerateDoublet, invz:orderedPhase,
  invzt:a1ZeroField}` → `ok = false`; all else rethrows.
- **Byte-parity invariant**: this spec's Component 5 becomes the governing
  copy of `invzt_run_phase_diagram.m`; the drivers spec's Component-1 block
  gets a one-line SUPERSEDED pointer.
- **The worktree carries the user's own WIP** (at spec time:
  `invz_tensor/invzt_run_spectra.m` with in-progress production knobs, plus
  other pre-existing modifications). Execution takes a fresh
  `git status --short` checkpoint at every task start and stages only
  named paths — no hard-coded dirt list (review F7).

## Components

### 1. `git mv invz_projected/invz_refine_crossing.m invz_common/invz_refine_crossing.m`

Zero content change. Gate: PROJECTED 143/0/19 (unchanged — the frozen
baseline count is a repo invariant, so the new direct helper test lives in
CORE, not in the projected suite), CORE unchanged at this task, INTEROP
unchanged.

### 2. `invz_tensor/invzt_crit_at.m` (promotion)

The driver's local function moves verbatim into a module file with a proper
header; the driver's local-function section is deleted in Component 5 (call
sites unchanged — same name; until then the identical local legally shadows
the module file inside the script).

```matlab
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
```

### 2b. Harden `invz_tensor/invzt_str.m` (safe formatter — never throws)

```matlab
function s = invzt_str(x)
%INVZT_STR Compact string form of x for error messages -- never throws.
%   s = INVZT_STR(x) is char(x) when x is a char row or a scalar string,
%   mat2str(x) for numeric/logical/string arrays, and a class/size
%   placeholder (e.g. '<1x1 struct>') for anything mat2str cannot format --
%   so an error MESSAGE can always be built, whatever malformed value
%   triggered it (mat2str itself throws on structs/cells/objects, which
%   would mask the intended error identifier; 2026-07-18 T-cut review F2).
%   Shared error-message helper for the invz_tensor drivers -- replaces the
%   per-file local_str / local_conv_str copies.
if ischar(x) || (isstring(x) && isscalar(x))
    s = char(x);
elseif isnumeric(x) || islogical(x) || isstring(x)
    s = mat2str(x);
else
    s = sprintf('<%s %s>', strjoin(cellstr(string(size(x))), 'x'), class(x));
end
end
```

Behavior-additive: char/string/numeric/logical inputs format exactly as
before; only inputs that previously THREW now return a placeholder. Fast
regression assert (`invzt_str(struct())` → `'<1x1 struct>'`) rides in the
validation test.

### 3. `invz_tensor/invzt_tc_pick.m` (new — the pure decision core)

```matlab
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
```

Policy notes (each is a committed fast test case): `[-1 1]` → bracket(1,2);
`[1 -1 1]` → bracket(2,3) with ncross 2 (re-entrance warn); `[1 2 3]` →
down; `[-1 -2]` → up; `[1 1 -1 -1]` → up (re-entrant lower leg — the
projected code breaks here); `[-1 0 1]` → zero(2) with ncross 1 (the
projected code errors here); `[-1 0 0 1]` → zero(3), ncross 1; `[1 0 -1]`
→ up (ordered above the zero ⇒ the true highest boundary is above);
`[-1 1 -1 0 1]` → zero(4), ncross 3; singletons `1`/`-1`/`0` →
down/up/zero(1).

### 4. `invz_tensor/invzt_critical_T.m` (new)

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
%   +0.15 K Tc error at the 8^3 round-trip gate) -- the same reason
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
```

### 5. `invz_tensor/invzt_run_phase_diagram.m` — two-regime rewiring

Full replacement (this block is the governing byte-parity copy; the drivers
spec's Component-1 block gets a SUPERSEDED pointer):

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
% OPTION NAMESPACES: solve_opts carries solver physics knobs shared by every
% point solve; the two finders get their OWN control knobs (Btol in tesla vs
% Ttol/Twidth/Tgridstep in kelvin -- both finders name their option 'tol',
% in different units, so they are merged deliberately at the call boundary,
% never shared).
%
% ERROR POLICY: the sweep absorbs ONLY per-point invzt:bracket (a genuine
% no-crossing outcome once Ts/Bs/Twin/Brange are preflighted below); every
% other identifier that ESCAPES the finders rethrows. The finders' own
% internal sampler additionally classifies shared-engine physics signals as
% invalid/ordered votes (their committed policy, documented in their headers).

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs ------------------------------------------------------------------
Ts     = linspace(0.4, 1.4, 11);       % fixed-T FIELD-CUT grid (K); low-T branch.
                                       % [] disables field cuts.
Bs     = linspace(0.25, 1.5, 6);       % fixed-B TEMPERATURE-CUT fields (T); the
                                       % near-vertical branch. [] disables T-cuts.
Twin   = [];                           % [] -> adaptive T-cut window anchored at the
                                       % small-Bx proxy Tc0 (computed below); or an
                                       % explicit HARD [Tlo Thi] bound (K) forwarded
                                       % to invzt_critical_T (no sliding; skips the
                                       % proxy anchor).
Brange = [0.05 6.0];                   % field-cut [Blo Bhi] bracket (T). Blo > 0:
                                       % mode 'a1' forbids exact zero transverse field.
Btol   = 0.02;                         % field-cut bracket tolerance (TESLA) ->
                                       % invzt_critical opts.tol
Ttol   = 0.005;                        % T-cut refinement tolerance (KELVIN) ->
                                       % invzt_critical_T opts.tol
Twidth = 0.5;                          % T-cut adaptive-window width (K)
Tgridstep = 1/30;                      % T-cut coarse-grid step (K)
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
% it per point without these checks would turn a typo into a silent all-NaN
% sweep). invzt_str, not mat2str: mat2str throws on structs while BUILDING the
% intended error message.
assert(isnumeric(Brange) && isreal(Brange) && numel(Brange) == 2 && ...
    all(isfinite(Brange)) && Brange(2) > Brange(1) && Brange(1) > 0, ...
    'invzt_run_phase_diagram:Brange', ...
    'Brange must be finite [Blo Bhi] with 0 < Blo < Bhi (mode ''a1'' forbids B = 0); got %s.', ...
    invzt_str(Brange));
assert(isempty(Ts) || (isnumeric(Ts) && isreal(Ts) && isvector(Ts) && ...
    all(isfinite(Ts)) && all(Ts > 0)), 'invzt_run_phase_diagram:Ts', ...
    'Ts must be empty or a finite positive real vector; got %s.', invzt_str(Ts));
assert(isempty(Bs) || (isnumeric(Bs) && isreal(Bs) && isvector(Bs) && ...
    all(isfinite(Bs)) && all(Bs > 0)), 'invzt_run_phase_diagram:Bs', ...
    'Bs must be empty or a finite positive real vector (mode ''a1'' forbids B = 0); got %s.', ...
    invzt_str(Bs));
assert(isempty(Twin) || (isnumeric(Twin) && isreal(Twin) && numel(Twin) == 2 && ...
    all(isfinite(Twin)) && Twin(2) > Twin(1) && Twin(1) > 0 && Twin(2) > 0.02), ...
    'invzt_run_phase_diagram:Twin', ...
    ['Twin must be [] (adaptive) or a finite HARD [Tlo Thi] bound with 0 < Tlo < Thi ' ...
     'reaching above the 0.02 K solve floor; got %s.'], invzt_str(Twin));
% Finder control knobs: validate BEFORE the lattice build and proxy (the
% finders would reject or misuse them only inside a sweep job, after the
% expensive setup -- second review R4). invzt_critical does NOT validate its
% opts.tol at all, so a bad Btol would silently corrupt the bisection.
posk = @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x > 0;
assert(posk(Btol), 'invzt_run_phase_diagram:Btol', ...
    'Btol must be a finite positive real scalar (TESLA); got %s.', invzt_str(Btol));
assert(posk(Ttol), 'invzt_run_phase_diagram:Ttol', ...
    'Ttol must be a finite positive real scalar (KELVIN); got %s.', invzt_str(Ttol));
assert(posk(Twidth), 'invzt_run_phase_diagram:Twidth', ...
    'Twidth must be a finite positive real scalar (KELVIN); got %s.', invzt_str(Twidth));
assert(posk(Tgridstep), 'invzt_run_phase_diagram:Tgridstep', ...
    'Tgridstep must be a finite positive real scalar (KELVIN); got %s.', invzt_str(Tgridstep));

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

% Per-finder opts, merged deliberately at the call boundary (both finders name
% their tolerance 'tol' -- tesla for the field cut, kelvin for the T-cut).
bopts = solve_opts;  bopts.tol = Btol;
topts = solve_opts;  topts.tol = Ttol;  topts.width = Twidth;  topts.gridstep = Tgridstep;
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
            val = invzt_critical(ion, v, lat, Brange, bopts);
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

### 6. Tests — `invz_tensor/tests/test_invzt_critical_T.m` (new)

Standard CORE `setupOnce` boilerplate. **Five fast** test functions (pure
core policy ×2, validation, sampler contract, refine-crossing helper) and
**two `INVZ_SLOW`** integration gates (odd-off Bc↔Tc round-trip with a
length-MATCHED poisonous `Sigma_seed` and `out`-field assertions; a
committed odd-on T-cut). Full test code lives in the implementation plan
(Task 2 Step 1); the acceptance matrix:

| Area | Fast cases (synthetic, no solves) | Slow cases (INVZ_SLOW, 8³) |
|---|---|---|
| Crossing policy | `[-1 1]` bracket; `[1 -1 1]` highest + ncross 2; all-PM down; all-ordered up; re-entrant lower leg `[1 1 -1 -1]` up; singletons | Bc↔Tc round-trip, `AbsTol` 0.05 K |
| Exact zeros | `[-1 0 1]` zero, ncross 1; zero-run `[-1 0 0 1]`; `[1 0 -1]` up; zero above bracket `[-1 1 -1 0 1]`; lone `0` | — |
| Validation | `invzt_str(struct())` placeholder (the hardened formatter); `invzt:tcAnchor` (missing/vector/complex/negative/below-floor `Tc0`); `invzt:tcWindow` (reversed/struct/below-floor); `invzt:tcOpts` (zero gridstep, negative tol, Inf width) — malformed-opts fixtures built by FIELD ASSIGNMENT (the `struct('f',[a b])` constructor makes struct arrays); CHEAP bracket diagnostics via the B = 0 all-invalid trick (pre-lattice `invzt:a1ZeroField`, dummy `lat`, no solves): hard-window `invzt:bracket` message reports the EXACT user window, and the adaptive no-valid grow path terminates at the floor (T1/T3) | hard window `[1.0 1.8]` used by the round-trip; deep-ordered hard window errors `invzt:bracket` whose MESSAGE reports the unmutated user window `[0.250, 0.450]` (one pass, no slide-mutated diagnostics — F4 + R2 end-to-end, in the odd-on test) |
| Sampler | zero-field absorbed (`ok=false`, `c=NaN`); bad mode rethrows `invzt:mode` — both pre-lattice, dummy `lat` | valid votes on both signs (implicit in round-trip); `sigma_floor = Inf` invalidates a finite-crit converged point (threading + rejection, one solve — R6) |
| Helper | `invz_refine_crossing`: linear bracket to `1e-3`; interior dead-zone with midpoint recovery/interp fallback | — |
| Seed | — | POISONOUS NaN seed, length-matched to the top endpoint's Matsubara count (a wrong-length seed is silently ignored by the solver's `numel == nwn` guard — a 7-element seed proves NOTHING, second review R1); with the strip, `out.ok(end)` stays a valid voter |
| odd-on | — | `Tc(1.5 T)` finite in `[1.35, 1.60]`, window `[1.2 1.7]` |
| out struct | — | fields `Tg/c/ok/window/ncross/B` present, sizes consistent |

CORE becomes **55 passed / 0 failed / 3 incomplete** (50 + the five fast
tests; the two slow tests join the A4 ladder as filtered-incomplete).

### 7. Driver smoke (verification, not a committed test)

Same-directory `sed` copy: `Ts` → `[1.2 1.4]`, `Bs` → `[0.5 1.5]`,
`gridN` → 8, `dpRng` → 15, `useParallel` → false,
`show_projected_anchor` → true. `Ts_proxy` stays (8³ Tc0 ≈ 1.53 sits inside
it). Assertions: both `Bc` finite (warm-cache knowns ≈ 2.148/1.916 T);
**both `TcB` finite with `odd = true`** (the driver default — the empirical
two-branch demonstration); `TcB(1) > TcB(2)`;
`TcB(1) <= Tc0_proxy + 0.02`; finite comparator; `phase_boundary` ≥ 3 rows;
figure exists. Expected T-cut values: both between 1.40 K and the proxy
Tc0 ≈ 1.53 K at this coarse grid (the 8³ boundary has Bc(1.4 K) = 1.916 T >
1.5 T) — magnitudes reported, only the ordering gated. Durations are
estimates, never acceptance criteria (review F8). Delete the smoke copy
after.

## Error handling

- Preflights before compute: `Brange`/`Ts`/`Bs`/`Twin` at the driver
  (`invzt_run_phase_diagram:*`, `invzt_str`-formatted), controls/window/
  anchor at the finder (`invzt:tcOpts`/`invzt:tcWindow`/`invzt:tcAnchor`),
  proxy-anchor availability (`invzt_run_phase_diagram:tcAnchor`).
- Per-point: only `invzt:bracket` absorbed (both cut kinds), all else
  rethrows; the sampler's selective-catch policy is unchanged (Component 2).
- Hard window: one pass, `invzt:bracket` with a widen-the-window message.
  Adaptive: up to 9 window attempts; collapse at the `Tmin` floor
  terminates into `invzt:bracket` (no infinite resampling).
- `invzt:multipleCrossings` is a warning (re-entrance is physics to report,
  never mask); `invzt:tcNoWindow` from the proxy is caught, the anchor
  assert then decides whether it matters.

## Testing / verification

1. Component 1 first: `git mv`, then PROJECTED (143/0/19) + CORE (50/0/1)
   + INTEROP (8/0/2) — all unchanged.
2. Components 2–4 + 6, TDD: fast tests written first and RED
   (`MATLAB:UndefinedFunction` ≠ the expected identifiers), then the three
   module files land, then fast GREEN, then the `INVZ_SLOW` gate on the new
   file (7/0/0), then CORE **55/0/3** + INTEROP 8/0/2.
3. Component 5 + smoke (Component 7). Suites re-run at the end:
   CORE 55/0/3, PROJECTED 143/0/19.
4. Docs (review F6 — full sweep, not spot edits): README §1 quickstart
   count → 55/0/3 with the incompletes named; §2 "Drivers" callout — BOTH
   the parenthetical ("PM-side field-cut boundary…") and the scope sentence
   rewritten for two-regime search; §2.1 heading and intro mention both
   finders, with a one-line T-cut example added to the recipe; the
   zero-field note corrected (`invzt_tc_pm_extrap` is the small-Bx PROXY —
   the only true-B=0 route in scope is the projected closed form); §7's
   slow-gate sentence names all three `INVZ_SLOW` tests; module-map rows
   for `invzt_critical_T`, `invzt_tc_pick`, `invzt_crit_at`; architecture
   A1 row. Structure check asserts the stale phrases are gone (not just
   one). Drivers spec Component-1 gets the SUPERSEDED pointer.
5. Every task starts with a fresh `git status --short` checkpoint; only
   named paths are staged (the worktree carries user WIP — at spec time
   `invz_tensor/invzt_run_spectra.m` — which must never be staged or
   reverted).

## Success criteria

- `invzt_critical_T` returns Tc(B) on the tensor branch with `odd` on and
  off; the slow round-trip and odd-on gates pass; every `invzt_tc_pick`
  policy case from the acceptance matrix passes in milliseconds.
- Malformed `Tc0`/window/controls produce their documented identifiers
  (never `MATLAB:nonLogicalConditional` or a `mat2str` throw).
- The driver produces a merged two-branch boundary in one run; smoke shows
  both T-cut points finite (odd on) and ordered correctly.
- `invz_refine_crossing` lives in `invz_common/` with direct CORE coverage;
  suites at CORE 55/0/3, INTEROP 8/0/2, PROJECTED 143/0/19.
- README reflects the new capability with no stale field-cut-only or
  proxy-as-true-B=0 phrasing.

## Review dispositions (`critical-t-review_by_Codex.md`, F1–F8)

- **F1 (High) — accepted, verified by hand-simulation.** Exact-zero and
  re-entrant-lower-leg failures confirmed; both inherited verbatim from the
  projected `invz_critical_T` (flagged as a candidate projected-side
  follow-up, deliberately NOT touched here). Fixed via the pure
  `invzt_tc_pick` core with committed synthetic tests for every policy case.
- **F2 (High) — accepted, MATLAB behavior confirmed** (`mat2str(struct())`
  throws; vector through `&&` throws; `isfinite(1+1i)` true). Strict
  scalar/real/positive validators, new `invzt:tcOpts` identifier,
  `invzt_str` everywhere (the module already owned the safe formatter),
  below-floor window rejection, floor-collapse termination, driver `Ts`
  preflight.
- **F3 (High) — accepted.** The pure core enables the full fast acceptance
  matrix; sampler-contract tests ride the verified pre-lattice guard paths;
  direct `invz_refine_crossing` test added in CORE (PROJECTED's frozen
  count is a repo invariant); seed-strip proven by a length-MATCHED
  poisonous seed in the slow round-trip; odd-on T-cut is now a COMMITTED `INVZ_SLOW` test,
  not just the throwaway smoke.
- **F4 (Medium) — accepted, reviewer's preferred option**: explicit window
  = hard bound, no sliding; adaptive mode slides. Documented in the finder
  header and the driver's `Twin` knob.
- **F5 (Medium) — accepted, collision confirmed** (`tol` = 0.02 T vs
  0.005 K under one name). Driver knobs `Btol`/`Ttol`/`Twidth`/`Tgridstep`
  merged into `bopts`/`topts` at the call boundary; `Ts` preflight added.
- **F6 (Medium) — accepted**, including the reviewer's catch that the
  README's zero-field note mislabeled the small-Bx proxy as a true-B=0
  route. Full doc sweep in Testing item 4 with stale-phrase assertions.
- **F7 (Low) — accepted, confirmed by `git status`** (the user's own
  in-progress edits to `invzt_run_spectra.m` were present; the earlier
  review file had been replaced). Hard-coded dirt lists dropped in favor of
  per-task status checkpoints + exact-path staging.
- **F8 (Low) — accepted**: "up to 9 window attempts" wording; Task-1
  coverage claim corrected (fast helper coverage = interop field-cut parity
  + the new direct CORE test); durations documented as estimates.

## Second-pass review dispositions (R1–R6, 2026-07-19)

- **R1 (High) — accepted, confirmed at source** (`invzt_solve_point.m:274`:
  `numel(opts.Sigma_seed) == nwn` guard; measured nwn 75/54/43 at
  1.0/1.4/1.8 K, so a 7-element seed is silently ignored and the old test
  passed with the strip deleted). Fixed: poisonous NaN seed length-matched
  to the top endpoint via `invz_matsubara(window(2), 40)`; `out.ok(end)`
  asserted.
- **R2 (Medium) — accepted, confirmed by re-reading the rev-2 loop**: the
  slide switch mutated `[Tlo Thi]` before a hard-window exit, so the error
  reported a window never sampled. Fixed: `if hardwin, break; end` BEFORE
  the slide switch + all reporting through the last-sampled `winS`; the
  odd-on test asserts the message contains `[0.250, 0.450]`.
- **R3 (Medium) — accepted**: the no-valid `grow` path re-sampled the
  identical floor-clamped grid for the remaining attempt budget. Fixed:
  terminate when the sampled `Tlo` is already at the floor. (The reviewer's
  optional pure window-transition helper was declined as over-engineering
  for a four-line switch; the floor guard is FAST-covered via the B = 0
  all-invalid trick — third review T3.)
- **R4 (Medium) — accepted**: `Btol`/`Ttol`/`Twidth`/`Tgridstep` preflighted
  as finite positive real scalars with stable driver identifiers, and Twin
  must reach above the 0.02 K solve floor — all before the lattice build
  (NB `invzt_critical` does not validate its `opts.tol` at all).
- **R5 (Medium) — accepted**: boundary indicators now accumulate across
  adaptive attempts (`nseen`) and feed both the warning and `out.ncross`;
  no-double-count argument documented in-code (crossing pairs are
  intra-window; a zero at a window top returns immediately).
- **R6 (Low) — accepted**: direct `sigma_floor` coverage added to the slow
  round-trip (one extra solve: `sigma_floor = Inf` must yield finite `c`
  with `ok = false`).

## Third-pass review dispositions (T1–T3, 2026-07-19)

- **T1 (High) — accepted, reviewer-reproduced**: MATLAB's `verifyError`
  outputs are the evaluated function's outputs (`<missing>` when it throws),
  never the caught `MException` — the rev-3 message assertion failed on
  CORRECT behavior. Fixed: explicit try/catch asserting identifier AND
  message; plus the reviewer's cheap variant adopted (B = 0 + dummy `lat`
  makes every sample invalid pre-lattice, exercising the hard-window
  diagnostic with zero physics solves).
- **T2 (Medium) — accepted, confirmed in the rev-3 block**: `winS` was
  assigned before the collapse guard, so a floor collapse reported a
  never-sampled degenerate pair as the "last sampled window". Fixed:
  capture after the guard; degenerate hard-window spans (≤ 1e-9) rejected
  at validation.
- **T3 (Low) — accepted**: adaptive floor termination now fast-covered
  (same B = 0 trick, `Tc0` just above the floor); the plan's error-policy
  list gains the four driver control ids; every "incompatible seed" phrase
  corrected to length-MATCHED / poisonous-VALUED.

## Execution findings

- **E1 (2026-07-19, caught by the Task-2 integration gates — the reason
  those gates exist):** byte-verified transcriptions + all-green synthetic
  tests, yet `crit(T)` at fixed field flipped sign 5× (round-trip 1.596 K
  vs 1.4 K; odd-on above the zero-field Tc). Cold-started samples hop
  between A1 fixed-point branches; the field-cut finder avoids this only
  via its warm-start seed. Fixed with the descending branch-tracked rolling
  seed (Component 4's nested sampler; `invzt_crit_at` third output).
  Caller-seed stripping and its poisonous-seed test are unchanged (the
  chain's FIRST/top sample is cold, so an unstripped caller seed still
  poisons exactly that endpoint).

## Out of scope

- Ordered-phase tensor solve (still the FM blocker to driver unification).
- Fixing the projected `invz_critical_T`'s own exact-zero /
  re-entrant-lower-leg gaps (F1's twins) — validated code with its own test
  surface; flagged to the user as a candidate follow-up.
- The projected smoothness regression (`test_tc_boundary_is_smooth`) — a
  production-sweep-scale test; defer until someone runs a production tensor
  T-cut sweep.
- Warm-starting across T (documented impossibility with T-dependent
  Matsubara grids; revisit only with an interpolating seed scheme).
