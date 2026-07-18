# Tensor-branch drivers: `invzt_run_phase_diagram` / `invzt_run_spectra`

Date: 2026-07-18. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-18); scope and the peak-picker sharing decision were
each confirmed via an explicit multiple-choice pick — see Decisions below.
REVISED twice on 2026-07-18: first after a source-verification pass
(subagent review) that found a runtime-fatal defect in the first draft's
`Tc0`-proxy handle; second after an external review
(`invzt_driver_review_by_Codex.md`, findings F1–F8) whose blocker F1 —
the explicit-q branch of `invzt_chi_realaxis` real-projects the response,
making every q-path χ'' map identically zero — adds a prerequisite
Component 0 below. Full finding dispositions: "Review amendments" section.

## Decisions (user-approved, govern over any conflicting prose below)

1. **Scope: tensor-only drivers, PM-side, as-is.** Ship two new scripts in
   `invz_tensor/`, honestly scoped to what `invzt_solve_point` /
   `invzt_critical` / `invzt_chi_realaxis` support today. No new tensor
   solver *physics* (no ordered-phase branch, no `invzt_critical_T`) as a
   prerequisite. One narrow, non-physics module fix IS in scope
   (Component 0: the `invzt_chi_realaxis` explicit-q complex-response
   contract — a one-line projection bug plus a scope guard, gated by new
   CORE regression tests). `invz_projected`'s own drivers and module code
   stay unmodified except for the one `git mv` in item 2.
2. **`invz_peak_energy.m` moves to `invz_common/`** (`git mv`, zero logic
   change) so both `invz_projected` and the new `invzt_run_spectra.m` call
   the one shared copy, matching the precedent already set by the 16
   single-ion functions moved there for the same reason.

## Problem

`invz_projected` ships `invz_run_phase_diagram.m` and `invz_run_spectra.m` —
scripts a user just runs to get a phase diagram or a susceptibility
spectrum. `invz_tensor` has no equivalents: only the point-solver primitives
a driver would call (`invzt_solve_point`, `invzt_critical`,
`invzt_chi_realaxis`). A user on the tensor branch currently has to hand-roll
the sweep loop and plotting every time (see `invz_tensor/README.html` §2,
added as a stopgap in the prior session).

Porting the projected drivers verbatim is not possible: the tensor branch
has no ordered/FM solve at all (`invzt_solve_point` only ever returns a
paramagnetic fixed point — LOCKED convention 7) and no temperature-cut root
finder (`invzt_critical_T` does not exist, only the field-cut
`invzt_critical`). Both gaps are real physics/numerics work, not driver
plumbing, and are out of scope here (Decision 1).

## Solution overview

One prerequisite module fix (Component 0), then two new files,
`invz_tensor/invzt_run_phase_diagram.m` and `invz_tensor/invzt_run_spectra.m`,
structurally mirroring their projected counterparts but reduced to the
PM-side field-cut-only physics the tensor branch actually has, plus
tensor-only knobs (`mode`/`nlevels`; `dress` only where A3 is reachable)
the projected driver has no equivalent of. One shared numeric helper
(`invz_peak_energy.m`) moves to `invz_common/` so `invzt_run_spectra.m` can
use it without a projected dependency.

Both drivers stay **self-contained by default**: neither requires
`invz_projected` on the path unless the user opts into the cross-model
comparison anchor (Component 2, `show_projected_anchor`).

## Verified facts this design rests on

Checked against source (two verification passes + the Codex review); no
longer assumptions:

- **The `invz_peak_energy.m` move is safe.** All ten tests that invoke
  `invz_spectra_map` / `invz_spectra_qpath` / `invz_chi_tensor_ref` share an
  identical `setupOnce` that does
  `addpath(fullfile(here, '..', '..', 'invz_common'))`, as do the three
  driver callers (`invz_run_spectra.m`, `invz_run_ip_scan.m`,
  `invz_run_tensor_ref.m`) and the tensor interop test
  `tests/interop/test_invzt_rpa_parity.m` (which reaches it transitively via
  `invz_chi_tensor_ref`). `invz_common/` has no name collision.
- **`imag(out.chi_uniform(3,3,:))` is the correct, already-positive χ''**
  — no sign flip. The projected chain un-flips its internal `G`-language at
  `chit = (-G0)./(1+Sw)`, and `invz_spectra_map.m:143` builds the plotted
  `S.chiz` as `imag(o.chi_cc_q(1,:)).'` with no extra minus; the tensor
  `chi_uniform` chain (`invzt_chi0_split` → `invzt_chi_rpa` → `u'*X*u`)
  introduces no additional sign. `test_invzt_rpa_parity.m` compares both
  un-negated.
- **`out.chi_cc_q` is (pre-fix) real by construction — Codex F1,
  confirmed.** The explicit-q branch accumulates
  `acc = acc + real(Xq(3*(s-1)+3, 3*(s-1)+3, iq))` — a transplant of
  `invzt_gcc_lattice`'s Matsubara-axis pattern (where `real()` is a
  legitimate noise-clean) onto the real axis, where it deletes the entire
  dissipative part. No production code consumes tensor `chi_cc_q` today;
  the only CORE test on it asserts size/finiteness (zeros pass — Codex F2),
  and the interop test's `imag(op.chi_cc_q(1,:))` is the **projected**
  function's complex output, which is exactly the contract the tensor side
  must match. Component 0 fixes this.
- **`invzt_tc_pm_extrap` does NOT guard its `critfun`.** It calls
  `[c(i), okv(i)] = critfun(Tg(i))` bare (line 46); any exception propagates
  and aborts the whole call. The projected analogue's effective critfun,
  `invz_crit_at.m`, supplies its own *selective* try/catch. The driver's
  local `invzt_crit_at` does the same — see Component 1.
- **`invzt_solve_point` throws on a realistic sweep.** It contains no
  try/catch. Sweep-plausible identifiers: `invzt:a1ZeroField` (transverse
  field below 1e-6 T), and — propagated from `invz_common/invz_twolevel.m`
  — `invz:degenerateDoublet` (splitting < 1e-4 meV) and `invz:orderedPhase`
  (|m| > 1e-3). A *failed PM fixed point* deep in the ordered phase does
  NOT throw: it returns `pt.converged = false` cleanly. It records the
  physics layer in `pt.mode` (`invzt_solve_point.m:370`) — the hook for
  Component 0's scope guard.
- **`invzt:bracket` doubles as `invzt_critical`'s argument-validation
  identifier** (Codex F5, confirmed at `invzt_critical.m:61-63`): a
  malformed or reversed `Brange` raises the same id as a genuine
  no-crossing window. A driver that absorbs `invzt:bracket` per point must
  therefore preflight `Brange` itself, or a reversed bracket becomes a
  silent all-NaN sweep.
- **`invzt_critical`'s nested `sample()` absorbs every non-`invzt:*`
  exception as an invalid/ordered vote** (Codex F7, confirmed at
  `invzt_critical.m:130-142`) — its committed T7 "non-convergence = phase
  signal" policy. The drivers' fail-loud guarantee is therefore scoped:
  *the driver rethrows every non-bracket error that escapes
  `invzt_critical`; the finder's own broader internal classification stands
  as designed.* `sigma_floor` is an existing `invzt_critical` option
  (default −0.5) — the drivers read the same name from `solve_opts` so the
  floor is single-sourced.
- **`invz_odd_zero_field` requires `ion.demag == 0`** for every `mode`
  (guard in its callee `invz_odd_blocks.m`, fires before any mode branch),
  solves **seven** variants per grid at `mode='full'`, and always builds
  its own **legacy endpoint-inclusive `[-0.5, 0.5]` mesh** via
  `qVec_generator` (Codex F6) — a different quadrature convention from the
  tensor driver's `'halfopen'` grid even at equal nominal N.
- **Strict Γ in an explicit q-list is accepted, not rejected** —
  `invzt_jq_tensor` assembles Γ-equivalent points with the Lorentz cavity
  (`invz_is_gamma_equiv` gate), i.e. a strict-uniform-mode page, which is
  the strict-q=0 *observable*, not the q→0⁺ intrinsic limit (the projected
  README §1.6 distinction). Convention decision (Codex F8): dispersion
  q-paths in the drivers, recipes, and smoke tests are **Γ-excluded**;
  every example starts at finite q.

## Components

### 0. Prerequisite fix: `invzt_chi_realaxis` explicit-q complex contract

Two changes to `invz_tensor/invzt_chi_realaxis.m`, both TDD-gated by new
CORE tests in `invz_tensor/tests/test_invzt_chi_realaxis.m`:

**(a) Keep the explicit-q response complex** (fixes Codex F1). In the
explicit-q accumulation loop, drop the `real()` projection and preallocate
complex:

```matlab
% before (bug -- Matsubara-pattern transplant, deletes chi'' on the real axis):
out.chi_cc_q = zeros(nq, nw);
...
acc = acc + real(Xq(3*(s-1)+3, 3*(s-1)+3, iq));

% after:
out.chi_cc_q = complex(zeros(nq, nw));
...
acc = acc + Xq(3*(s-1)+3, 3*(s-1)+3, iq);
```

Update the function header's `chi_cc_q` description: COMPLEX per-q
sublattice-averaged site-diagonal cc response; `imag()` is the positive
χ'' — contract parity with the projected `invz_chi_realaxis`'s `chi_cc_q`
(whose `imag` the interop suite already consumes).

**(b) Enforce the documented A1-only scope** (fixes Codex F3 at the callee,
so every caller is protected):

```matlab
ptmode = getf(pt, 'mode', 'a1');
if ~strcmp(ptmode, 'a1')
    error('invzt:realaxisMode', ['invzt_chi_realaxis is the A1 scalar-Sigma ' ...
        'continuation ONLY (LOCKED scope; see the SCOPE box in this header): ' ...
        'got pt.mode = ''%s''. Re-solve the point with mode ''a1''. Continuing ' ...
        'the A2/A3 matrix objects is an open item (README section 10).'], ptmode);
end
```

New CORE tests (regression + guard; exact code in the implementation plan):
- `test_qsel_explicit_q_complex_response` — explicit q-list at a known PM
  point: `~isreal(out.chi_cc_q)` and positive dissipative weight
  (`max(imag) > 1e-6`) on the physical (full-Σ) call; the exact
  χ''(ω>0) non-negativity gate runs on a second `force_sigma0 = true`
  call (the bare-RPA limit, manifestly passive) — AS-EXECUTED AMENDMENT:
  the gate cannot run on the full-Σ call because the frozen-Kw A1
  continuation has a pre-existing near-resonance negative-χ'' artifact
  (−312.9 beside a +652.5 peak at the test point; present identically in
  the untouched `chi_uniform` path AND in the projected reference, so
  inherited, not introduced). The function header carries the caveat;
  "imag() = χ''" is the repository's no-extra-sign-flip convention, with
  passivity established only in the bare-RPA limit.
- `test_realaxis_rejects_non_a1_point` — a synthetic `pt.mode = 'a2'`
  override of a solved a1 point fed to `invzt_chi_realaxis` must raise
  `invzt:realaxisMode` (AS-EXECUTED AMENDMENT per second-review R6: the
  override isolates the guard from A2 solver behavior).
- `test_qsel_explicit_q_odd_mask` (second-review R1) — with
  `pt.odd = false` and `force_sigma0`, the explicit-q response must equal
  (1e-12) a manual reconstruction through `invzt_odd_mask(latq.Jt)`, and
  differ >1% from the unmasked route.

**(c) Honor `pt.odd` at explicit q** (second-review R1, added post-execution):
the odd=false Cartesian-off-diagonal zeroing rule is extracted from
`invzt_solve_point` into the shared helper `invz_tensor/invzt_odd_mask.m`
(byte-identical semantics, behavior-neutral refactor) and applied by the
continuation when `~pt.odd`: the explicit-q `latq.Jt` (the load-bearing
site — unmasked it gave a 17.2% response error) and the `'gamma'` `Jpage`
(a C2 no-op at Γ, ~1e-19, noted in-code); the Cartesian-diagonal
`'gamma_uniform'` `Jd` needs no mask by construction.

CORE suite expectation: **49 passed / 0 failed / 1 incomplete** after (a)+(b)
(47 + those 2 tests), **50 / 0 / 1** after (c) — the final gate.

### 1. `invz_tensor/invzt_run_phase_diagram.m`

```matlab
%INVZT_RUN_PHASE_DIAGRAM  Full-tensor 1/z PM-side phase boundary Bc(T).
%
% SCOPE (read this before trusting the plot): FIELD-CUTS ONLY. This driver
% bisects Bc at each fixed T via invzt_critical. There is no temperature-cut
% finder on the tensor branch (no invzt_critical_T), so the near-vertical
% part of the boundary approaching Tc(0) is NOT resolved -- those Ts points
% simply fail their bracket and come back NaN. That region is left BLANK by
% design, the same way invz_run_phase_diagram's ODD-overlay branch leaves its
% own near-Tc0 region blank. There is also no ordered-phase tensor solve, so
% nothing below the boundary is computed at all.
%
% ERROR POLICY: the sweep absorbs ONLY per-point invzt:bracket (a genuine
% no-crossing outcome once Brange itself has been preflighted below); every
% other identifier that ESCAPES invzt_critical rethrows. invzt_critical's
% own internal sampler additionally classifies shared-engine physics signals
% (and other non-invzt errors) as ordered votes -- its committed T7 policy,
% documented in its header, not overridden here.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs ------------------------------------------------------------------
Ts     = linspace(1.0, 2.2, 25);       % temperature grid (K)
Brange = [0.05 6.0];                   % [Blo Bhi] bisection bracket (T). Blo > 0:
                                       % mode 'a1' forbids exact zero transverse field.
gridN  = 16;  gridConv = 'halfopen';   % invzt_qgrid(gridN, gridConv)
dpRng  = 30;                           % invzt_jq_tensor coupling-sum range
useParallel = true;                    % false -> force serial
solve_opts  = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', 'dress', 'full');
                                       % sigma_floor may be added here too; defaults to
                                       % invzt_critical's own -0.5 (single-sourced below)

show_proxy_anchor     = true;          % tensor-native small-Bx Tc(0) proxy
                                       % (invzt_tc_pm_extrap). Costs ~numel(Ts)
                                       % EXTRA serial point solves -- set false to skip.
Bproxy                = 0.05;          % transverse field (T) for that proxy
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

% Brange preflight (invzt:bracket doubles as invzt_critical's arg-validation
% id; absorbing it per point without this check would turn a reversed/typo'd
% bracket into a silent all-NaN sweep).
assert(isnumeric(Brange) && isreal(Brange) && numel(Brange) == 2 && ...
    all(isfinite(Brange)) && Brange(2) > Brange(1) && Brange(1) > 0, ...
    'invzt_run_phase_diagram:Brange', ...
    'Brange must be finite [Blo Bhi] with 0 < Blo < Bhi (mode ''a1'' forbids B = 0); got %s.', ...
    mat2str(Brange));

if show_projected_anchor
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_projected'));
    assert(ion.demag == 0, 'invzt:anchorDemag', ...
        ['show_projected_anchor requires ion.demag == 0 (invz_odd_blocks is ' ...
         'intrinsic-only). Set ion.demag = 0 or disable the anchor.']);
end

g   = invzt_qgrid(gridN, gridConv);
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

nWorkers = 0;
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

nT = numel(Ts);
Bc = nan(1, nT);
parfor (k = 1:nT, nWorkers)
    t0 = tic;  val = NaN;
    % val-then-assign (NOT `Bc(k) = ...` inside the try): keeps the sliced
    % output variable unconditionally assigned on every iteration, matching
    % invz_run_phase_diagram's own parfor pattern.
    try
        val = invzt_critical(ion, Ts(k), lat, Brange, solve_opts);
    catch err
        % invzt:bracket here is a genuine per-point outcome (Brange was
        % preflighted above, so this window simply lacks a crossing at this
        % T -- expected near Tc0, where the boundary is near-vertical).
        % Anything else is a MISCONFIGURATION that would otherwise silently
        % produce a whole sweep of NaNs, so it rethrows.
        if ~strcmp(err.identifier, 'invzt:bracket'), rethrow(err); end
        fprintf('  T = %.2f K: no bracket in [%.2f %.2f] T (left NaN)\n', ...
                Ts(k), Brange(1), Brange(2));
    end
    fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT, Ts(k), val, toc(t0));
    Bc(k) = val;
end

% ---- zero-field anchors ------------------------------------------------------
Tc0_proxy = NaN;  Tc0_closed = NaN;
if show_proxy_anchor
    % Pass the FULL Ts: invzt_tc_pm_extrap filters to converged-PM points
    % itself and extrapolates the two lowest, so no hand-picked cutoff.
    critfun = @(T) invzt_crit_at(ion, T, [Bproxy 0 0], lat, solve_opts);
    try
        Tc0_proxy = invzt_tc_pm_extrap(critfun, Ts);
    catch err
        if ~strcmp(err.identifier, 'invzt:tcNoWindow'), rethrow(err); end
        fprintf('  Tc(0) proxy: no PM extrapolation window on this Ts grid (skipped)\n');
    end
end
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
if isfinite(Tc0_proxy)
    plot(0, Tc0_proxy, 'k^', 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('tensor small-B_x proxy T_c(0), %.2f T', Bproxy));
end
if isfinite(Tc0_closed)
    plot(0, Tc0_closed, 'ks', 'MarkerFaceColor', 'c', ...
         'DisplayName', 'projected closed-form T_c(0) (comparator; legacy-inclusive mesh)');
end
xlabel('B_c (T)'); ylabel('T (K)');
title('LiHoF_4 full-tensor 1/z phase boundary (PM side, field-cuts only)');
legend('Location', 'southwest');

phase_boundary = sortrows([Ts(:) Bc(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);

% =============================================================================
function [c, ok] = invzt_crit_at(ion, T, B, lat, opts)
% Tensor analogue of invz_projected's invz_crit_at: one criticality sample,
% with the SELECTIVE try/catch invzt_tc_pm_extrap does not provide itself.
% Physics signals (degenerate doublet, ordered toy doublet, zero-field guard)
% become ok = false; structural/config errors rethrow.
try
    pt = invzt_solve_point(ion, T, B, lat, opts);
    c  = pt.crit;
    % Same three-part sample-validity rule as invzt_critical, with the SAME
    % single-sourced floor (invzt_critical's opts.sigma_floor, default -0.5).
    ok = pt.converged && isfinite(c) && pt.Sigma0 >= getf(opts, 'sigma_floor', -0.5);
catch err
    switch err.identifier
        case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
            c = NaN;  ok = false;      % phase/physics signal, not an error
        otherwise
            rethrow(err);              % invzt:mode, invzt:a1OddGamma, ... = misconfiguration
    end
end
end
```

### 2. `invz_tensor/invzt_run_spectra.m`

```matlab
%INVZT_RUN_SPECTRA  Full-tensor 1/z chi''_cc spectra vs transverse field, or along a q-path.
%
% SCOPE: PARAMAGNETIC SIDE ONLY. invzt_solve_point has no ordered-phase branch,
% so unlike invz_run_spectra this driver cannot sweep ACROSS Bc to show the soft
% mode hardening on both sides. Fields that land on the ordered side (or fail the
% sample-validity rule) are masked to NaN with a console note.
%
% solve_opts.mode must be 'a1' (enforced below AND by invzt_chi_realaxis's own
% invzt:realaxisMode guard): the real-axis continuation exists only for the A1
% scalar Sigma. There is consequently no 'dress' knob here -- it is A3-only.
%
% ERROR POLICY: field sweep -- per-field selective masking (physics signals and
% invalid samples become NaN columns with a console note; everything else
% rethrows). q-path view -- FAIL LOUD by design: the whole product hinges on the
% single Bq point, so an invalid point raises invzt:qpathNotPM (with converged/
% crit/Sigma0) instead of drawing an empty map, and any physics throw from the
% one solve propagates for the same reason.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs (mirroring invz_run_spectra's names where the concept carries over) --
T      = 1.6;                          % K
fields = linspace(0.3, 4.5, 40);       % T -- keep on the PM side (see SCOPE above)
w      = linspace(0, 0.6, 400).';      % meV, uniform spacing (invz_peak_energy assumes it)
eta    = 5e-3;                         % meV, real-axis Lorentzian HWHM
sliceMax  = 6;                         % field count at/below which line slices are drawn
peak_wmin = 0.02;                      % meV -- lower bound for the peak pick. NON-ZERO by
                                       % default to exclude the low-frequency hyperfine
                                       % line, mirroring invz_spectra_map's opts.peak_wmin;
                                       % set 0 to pick over the whole w grid.
theta_c = 0.0;  phi_ab = 0.0;          % field-direction tilt knobs (deg). Ported as-is:
                                       % invzt_solve_point already takes a full [Bx By Bz]
                                       % and forwards transverse_mf.
transverse_mf = 'legacy_x';            % 'legacy_x' | 'none' | 'vector_ab'
gridN = 16; gridConv = 'halfopen'; dpRng = 30;
useParallel = true;
solve_opts = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', ...
                    'transverse_mf', transverse_mf);

% ---- q-path view: set qpath non-empty to switch views ------------------------
qpath = [];                            % [] = field sweep; [nq x 3] r.l.u. = q-path view.
                                       % Keep q-paths GAMMA-EXCLUDED: a strict q = [0 0 0]
                                       % is assembled with the Lorentz cavity (the strict-
                                       % uniform observable), NOT the q->0+ intrinsic limit
                                       % a dispersion plot wants -- start at finite q.
% qpath = [linspace(0.1, 0.5, 41).' zeros(41, 2)];   % example: toward (0.5,0,0)
Bq    = 2.0;                           % T -- fixed field for the q-path view
wq    = linspace(0, 0.6, 400).';       % meV -- own frequency grid for the q-path view
% -----------------------------------------------------------------------------

if ~strcmp(getf(solve_opts, 'mode', 'a1'), 'a1')
    error('invzt_run_spectra:mode', ['invzt_run_spectra requires solve_opts.mode = ' ...
        '''a1'' (got ''%s''): invzt_chi_realaxis is the A1 scalar-Sigma continuation ' ...
        'ONLY. A2/A3 real-axis continuation is an open item (README section 10).'], ...
        char(getf(solve_opts, 'mode', 'a1')));
end
if phi_ab ~= 0 && strcmp(transverse_mf, 'legacy_x')
    error('invzt_run_spectra:transverseMF', ...
        ['phi_ab = %.3g deg needs the vector transverse mean field: set transverse_mf ' ...
         'to ''vector_ab'' (or ''none'' for a bare CF+Zeeman diagnostic). legacy_x is ' ...
         'x-only and C4-inconsistent for rotated fields.'], phi_ab);
end
if ~isempty(qpath)
    % F8 convention preflight (R4, 2026-07-18 second Codex review): dispersion
    % q-paths must exclude strict Gamma (and any Gamma-equivalent row) -- see
    % the qpath knob comment above. Runs BEFORE any lattice/solve work.
    if ~(isnumeric(qpath) && isreal(qpath) && ismatrix(qpath) && size(qpath, 2) == 3 ...
            && all(isfinite(qpath(:))))
        error('invzt_run_spectra:qpath', ['qpath must be an [nq,3] real finite numeric ' ...
            'array of r.l.u. coordinates (see the qpath knob comment above); got a %s ' ...
            '%s.'], class(qpath), mat2str(size(qpath)));
    end
    if any(invz_is_gamma_equiv(qpath, ion.tau))
        error('invzt_run_spectra:qpathGamma', ['qpath must exclude Gamma-equivalent rows: ' ...
            'a strict q = [0 0 0] (or any Gamma-equivalent point) is assembled with the ' ...
            'Lorentz cavity -- the strict-uniform observable, NOT the q->0+ intrinsic limit ' ...
            'a dispersion plot wants (see the qpath knob comment above). Start the path at ' ...
            'finite q.']);
    end
end
dhat   = [cosd(theta_c)*cosd(phi_ab), cosd(theta_c)*sind(phi_ab), sind(theta_c)];
sfloor = getf(solve_opts, 'sigma_floor', -0.5);   % single-sourced with invzt_critical

g   = invzt_qgrid(gridN, gridConv);
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

if ~isempty(qpath)
    % ---------------- q-path dispersion at one fixed field --------------------
    pt = invzt_solve_point(ion, T, Bq*dhat, lat, solve_opts);
    if ~(pt.converged && isfinite(pt.crit) && pt.crit > 0 && pt.Sigma0 >= sfloor)
        error('invzt:qpathNotPM', ['q-path point B = %.2f T, T = %.2f K is not a ' ...
            'valid PM sample (converged = %d, crit = %.4g, Sigma0 = %.4g): the whole ' ...
            'q-path product hinges on this one point, so failing loudly beats ' ...
            'drawing an empty map. Raise Bq, raise T, or check the knobs.'], ...
            Bq, T, pt.converged, pt.crit, pt.Sigma0);
    end
    out = invzt_chi_realaxis(ion, T, Bq*dhat, pt, wq, ...
            struct('qsel', qpath, 'dpRng', dpRng, 'eta', eta));
    chipp_q = imag(out.chi_cc_q);                 % [nq, nw] positive chi'' (Component 0)
    figure; imagesc(wq, 1:size(chipp_q, 1), chipp_q); set(gca, 'YDir', 'normal');
    xlabel('\omega (meV)'); ylabel('q index along path'); colorbar;
    title(sprintf('tensor 1/z \\chi''''_{cc}(q,\\omega), T = %.2f K, B = %.2f T', T, Bq));
    Epeak_q = invz_peak_energy(chipp_q.', wq, peak_wmin);   % columns must be per-q
    figure; plot(1:numel(Epeak_q), Epeak_q, 'o-');
    xlabel('q index along path'); ylabel('E_{peak} (meV)');
else
    % ---------------- field sweep at the uniform mode -------------------------
    nWorkers = 0;
    if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end
    nf = numel(fields);
    chipp = nan(numel(w), nf);
    parfor (k = 1:nf, nWorkers)
        col = nan(numel(w), 1);
        try
            pt = invzt_solve_point(ion, T, fields(k)*dhat, lat, solve_opts);
            % Same three-part sample-validity rule as invzt_critical (converged,
            % finite positive crit, single-sourced Sigma0 floor) -- so the
            % spurious negative-Sigma fixed point invzt_critical warns about
            % never reaches the plot.
            if pt.converged && isfinite(pt.crit) && pt.crit > 0 && pt.Sigma0 >= sfloor
                o = invzt_chi_realaxis(ion, T, fields(k)*dhat, pt, w, ...
                        struct('qsel', 'gamma_uniform', 'eta', eta));
                col = squeeze(imag(o.chi_uniform(3, 3, :)));   % already positive chi''
            else
                fprintf('  B = %.2f T: ordered / non-converged / invalid (masked)\n', fields(k));
            end
        catch err
            switch err.identifier
                case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
                    fprintf('  B = %.2f T: %s (masked)\n', fields(k), err.identifier);
                otherwise
                    rethrow(err);
            end
        end
        chipp(:, k) = col;
    end

    if nf <= sliceMax
        figure; hold on;  co = lines(nf);
        for k = 1:nf
            plot(w, chipp(:, k), '-', 'Color', co(k, :), ...
                 'DisplayName', sprintf('%.2f T', fields(k)));
        end
        xlabel('\omega (meV)'); ylabel('\chi''''_{cc}');
        title(sprintf('tensor 1/z, T = %.2f K', T)); legend show;
    else
        figure; imagesc(fields, w, chipp); set(gca, 'YDir', 'normal');
        xlabel('B (T)'); ylabel('\omega (meV)'); colorbar;
        title(sprintf('tensor 1/z \\chi''''_{cc}(B,\\omega), T = %.2f K', T));
    end

    Epeak = invz_peak_energy(chipp, w, peak_wmin);
    figure; plot(fields, Epeak, 'o-');
    xlabel('B (T)'); ylabel('E_{peak} (meV)');
    title(sprintf('\\chi''''_{cc} peak energy vs field, T = %.2f K', T));
    % Gaps are CENSORED peaks (boundary max / non-positive column) or masked
    % ordered points -- do not interpolate over them.
end
```

### 3. `invz_peak_energy.m` → `invz_common/`

`git mv invz_projected/invz_peak_energy.m invz_common/invz_peak_energy.m`,
no content changes. Verified safe (see "Verified facts"): every test and
driver that reaches it already has `invz_common` on the path, and there is no
name collision. No caller edits are expected; the PROJECTED regression gate
confirms it.

## Error handling

Summary of the conventions used above, all inherited from existing module
practice rather than newly invented:

- **Configuration preflights fail loud before any compute**: malformed
  `Brange` (`invzt_run_phase_diagram:Brange`), non-a1 `solve_opts.mode`
  (`invzt_run_spectra:mode` at the driver, `invzt:realaxisMode` at the
  callee), `phi_ab ~= 0` under `legacy_x`
  (`invzt_run_spectra:transverseMF`), anchor with `ion.demag ~= 0`
  (`invzt:anchorDemag`).
- **Per-point boundary failure** (`invzt:bracket`, post-preflight):
  reported, left `NaN`, sweep continues. Every other identifier *that
  escapes `invzt_critical`* rethrows. (Inside `invzt_critical`, the nested
  sampler's committed policy additionally reads shared-engine physics
  signals — and other non-`invzt:*` exceptions — as ordered votes; that
  internal classification is out of the drivers' hands and stands as
  designed.)
- **Field-sweep physics signals** (`invz:degenerateDoublet`,
  `invz:orderedPhase`, `invzt:a1ZeroField`) and **invalid-but-non-throwing
  samples** (`pt.converged == false`, non-finite/non-positive `crit`,
  `Sigma0` below the single-sourced `sigma_floor`): masked column with a
  console note, never an error.
- **q-path invalid point**: deliberate FAIL LOUD (`invzt:qpathNotPM` with
  `converged`/`crit`/`Sigma0` in the message) — one invalid `Bq`
  invalidates the whole product, so no empty map is ever drawn; physics
  throws from that single solve propagate for the same reason.
- **`invzt:tcNoWindow`** from the `Tc0` proxy: caught, anchor skipped, rest
  of the plot still renders.

## Testing / verification

The drivers themselves get no MATLAB unit tests — matching precedent (the
projected `invz_run_*.m` have none; drivers are verified by running them).
Component 0 DOES get unit tests (it is module code). Verification:

1. Component 0 first, TDD: the complex-response regression test must FAIL
   against the current `real()` code, then pass after the fix; the
   `invzt:realaxisMode` guard test likewise. CORE becomes 49/0/1.
2. Run both drivers at reduced size (`gridN = 8`, `numel(Ts) ~ 4`,
   `numel(fields) ~ 5`, `useParallel = false`) via `matlab -batch` from the
   repo root; confirm they complete without error and produce non-empty
   figures. Exercise the spectra q-path branch too.
3. Smoke assertions are SEMANTIC, not just structural (Codex F2): the
   q-path smoke requires positive spectral weight
   (`max(chipp_q(:)) > 1e-6`) and at least one finite `Epeak_q` (if all
   censored, widen `wq` — do not weaken the gate); the field-sweep smoke
   requires at least one finite column AND at least one masked column at a
   temperature below the boundary top.
4. Exercise the opt-in anchor once (`show_projected_anchor = true`) to
   confirm the conditional `addpath` + `ion.demag` assert behave.
5. Re-run CORE (49/0/1 after Component 0(a)+(b), 50/0/1 after 0(c) — the final gate) and PROJECTED (143/0/19)
   after the `invz_peak_energy.m` move and at the end.
6. Update `invz_tensor/README.html` §2 to point at the real scripts,
   keeping the recipes as the "under the hood" explanation, and make its
   q-path recipe Γ-excluded (start at 0.1) so spec, plan, README, and smoke
   are identical (Codex F8).

Full-resolution production sweeps (the real `gridN = 16`+ / dense-`Ts` runs)
are left to the user to kick off, same as the projected driver's own
"production sweep... LEFT TO THE USER" note.

## Success criteria

- `invzt_chi_realaxis` returns a complex `chi_cc_q` with positive
  dissipative weight at a PM point, rejects non-a1 points with
  `invzt:realaxisMode`, honors `pt.odd` at explicit q via `invzt_odd_mask`,
  and the three new CORE tests pass (final gate 50/0/1).
- `invzt_run_phase_diagram.m` and `invzt_run_spectra.m` exist in
  `invz_tensor/`, run end-to-end at reduced size with no runtime errors
  (both branches of the spectra driver), and produce the plots described —
  with the q-path smoke showing genuinely nonzero χ''.
- Both drivers run with `invz_projected` absent from the path under their
  default knobs.
- `invz_peak_energy.m` lives in `invz_common/`; PROJECTED (143/0/19) and
  INTEROP (8/0/2) unchanged.
- `invz_tensor/README.html` §2 reflects the new scripts and the Γ-excluded
  q-path convention.

## Review amendments (Codex review, 2026-07-18)

Dispositions of `invzt_driver_review_by_Codex.md` F1–F8, each verified
against source before acceptance:

- **F1 (Blocker) — accepted.** `chi_cc_q` real-projection confirmed at the
  explicit-q accumulation. → Component 0(a) + regression test. The spec's
  earlier "verified fact" about the sign convention covered only
  `chi_uniform`; corrected above.
- **F2 (High) — accepted.** Structural-only smoke/unit gates would have
  certified the F1 zero-map. → semantic assertions (Testing item 3) + the
  Component 0 regression test. The reviewer's optional manual
  `invzt_chi_rpa` cross-projection was replaced by a χ''(ω>0)
  non-negativity physics gate — decisive without reimplementing the
  function's internals in the test.
- **F3 (High) — accepted, both layers.** Driver preflight
  (`invzt_run_spectra:mode`) + callee guard (`invzt:realaxisMode`,
  Component 0(b)); `dress` removed from the spectra driver's knobs (A3-only,
  dead under forced a1; still exposed in the phase driver where A3 runs).
- **F4 (Medium) — accepted.** q-path validity now includes the `Sigma0`
  floor; policy decision: FAIL LOUD (`invzt:qpathNotPM`) for the q-path,
  per-column masking only for the multi-field sweep; floor single-sourced
  as `solve_opts.sigma_floor` (existing `invzt_critical` option name,
  default −0.5) across both drivers and the proxy helper.
- **F5 (Medium) — accepted.** `Brange` preflight added before the `parfor`.
  Splitting `invzt_critical`'s arg-validation identifier from its
  no-crossing identifier is noted as a possible module follow-up, not done
  here (out of driver scope).
- **F6 (Medium) — accepted as a wording/labeling fix.** The anchor comment
  and plot legend now state: same nominal N and dipole cutoff, projected
  legacy-inclusive mesh, different model and ODD treatment — comparator
  only. Running the tensor curve on `'legacy_inclusive'` when the anchor is
  enabled was considered and declined: the anchor is already declared a
  cross-model comparator, and quadrature is the smallest of its three
  differences.
- **F7 (Medium) — spec qualified; callee change declined.** The drivers'
  fail-loud guarantee is now explicitly scoped to errors that escape
  `invzt_critical`. Narrowing `invzt_critical`'s internal absorbed set
  would alter its committed T7 "non-convergence = phase signal" policy —
  a module-behavior change with its own test surface, out of scope for a
  driver task.
- **F8 (Low) — accepted.** Convention decision: dispersion q-paths are
  Γ-excluded (strict Γ is *accepted* by `invzt_jq_tensor` but assembles the
  Lorentz-cavity strict-uniform page — the strict-q=0 observable, not the
  q→0⁺ limit); spec, plan, README recipe, and smoke now all start q-paths
  at finite q, and the driver comment documents the distinction. Renaming
  the spec to drop `-design` declined — the suffix is this repo's
  established spec-filename convention.

## Second-review amendments (Codex review round 2, R1–R6, 2026-07-18)

The second review was written against `266e799` (pre-execution); dispositions
verified against the executed state:

- **R1 (High) — accepted, fixed** (commit `1ae0c7f`): explicit-q continuation
  ignored `pt.odd = false` (17.2% response error). Shared `invzt_odd_mask`
  helper + regression test; Component 0(c) above.
- **R2 (High) — already resolved during execution**: the very failure the
  reviewer measured (−312.9/+652.5) was hit by the Task-1 implementer and
  resolved by moving the causality gate to the `force_sigma0` bare-RPA call
  (min = 0.319 ≥ 0 there — passivity established empirically in that limit,
  so the scoped claim stays rather than dropping the gate entirely as the
  reviewer suggested). Header carries the caveat; the reviewer's observation
  that the projected reference also shows negative lobes (−496/+671) is
  recorded as a shared known limitation.
- **R3 (High) — already resolved during execution**: the invalid
  `runtests('dir/file/testname')` selector form was hit by the Task-1
  implementer; working forms are the `.m`-suffixed whole-file run or
  `runtests('...file.m', 'ProcedureName', 'name')`. The plan's three
  commands are corrected in its execution addendum.
- **R4 (Medium) — accepted, fixed** (commit `0157477`): `qpath` preflight
  (`invzt_run_spectra:qpath` / `invzt_run_spectra:qpathGamma`) using the
  module's own `invz_is_gamma_equiv` gate — stricter than the suggested
  `abs(q) <= 1e-12` (catches Γ-equivalent points). Component-2 block above
  is regenerated byte-true from the committed driver.
- **R5 (Medium) — two bullets accepted** (README quickstart count now 50/0/1;
  SESSION-note items for the mode guard and the odd-off mismatch retired as
  resolved); **third bullet refuted with evidence**: the A1 proxy-Tc
  comparison IS grid-matched — both sides used the identical
  `legacy_inclusive` 16³ mesh (`tests/interop/test_invzt_critical_parity.m:64`
  explicitly passes `'legacy_inclusive'`; `tests/invzt_anchors.m:128` pins
  `grid = 'legacy_inclusive'`; `invz_odd_zero_field` builds the byte-identical
  `qVec_generator([-0.5 0.5])` mesh). The reviewer's premise assumed the
  tensor default (`'halfopen'`) was used; it was explicitly overridden. The
  README headline wording stands. (The driver's anchor caveat from F6 remains
  correct — the DRIVER defaults to `'halfopen'`.)
- **R6 (Low) — accepted, fixed** (commit `1ae0c7f`): mode-guard test uses a
  synthetic `pt.mode` override, no A2 solve.

## Out of scope

- `invzt_critical_T` (temperature-cut root finder) and any ordered-phase
  tensor solver — the actual blockers to full parity; explicitly deferred
  per Decision 1, not attempted as part of this task.
- Any change to `invz_projected`'s own drivers beyond the one `git mv`.
- Unifying the two branches' drivers behind a shared `invz_common` entry
  point / model toggle — considered and explicitly declined for now (the
  adapter work is only worth it once the tensor branch has FM + T-cut
  support; revisit then).
- Narrowing `invzt_critical`'s internal non-`invzt:*` absorption, or
  splitting its arg-validation identifier from `invzt:bracket` (Codex
  F5/F7 module-side suggestions) — noted as possible module follow-ups.
- A2/A3 real-axis continuation (the reason for the a1-only guards).
- Warm-starting `Sigma_seed` across sweep points (the projected drivers do
  not do this either; `invzt_critical` warm-starts only *within* one
  bisection).
