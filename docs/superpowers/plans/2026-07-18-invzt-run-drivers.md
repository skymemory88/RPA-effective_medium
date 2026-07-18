# Tensor-Branch Drivers (`invzt_run_phase_diagram` / `invzt_run_spectra`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the full-tensor branch the same "just run it" driver scripts the projected branch has — a PM-side phase-diagram sweep and a χ''_cc spectra driver — scoped honestly to what the tensor solvers support today, after first fixing the `invzt_chi_realaxis` explicit-q complex-response bug that would make every q-path χ'' map identically zero.

**Architecture:** One prerequisite module fix (drop a `real()` projection in `invzt_chi_realaxis`'s explicit-q branch + add an A1-only scope guard, TDD-gated), then two new standalone MATLAB *scripts* in `invz_tensor/` mirroring the projected `invz_run_*.m` pattern (knob block at top, `parfor` sweep, figures + workspace exports, local helpers at the bottom). One 30-line numeric helper (`invz_peak_energy.m`) moves from `invz_projected/` to `invz_common/`. `invz_tensor/README.html` §2 is updated last.

**Tech Stack:** MATLAB R2025a; existing module primitives `invzt_qgrid`, `invzt_jq_tensor`, `invzt_solve_point`, `invzt_critical`, `invzt_tc_pm_extrap`, `invzt_chi_realaxis`; shared engine in `invz_common/`.

**Spec (source of record):** `docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md` (rev 2, incorporating the Codex review dispositions F1–F8). The driver code blocks below are copied from that spec verbatim; if they diverge, the spec governs. The external review itself is `invzt_driver_review_by_Codex.md` (repo root, user-provided — do not commit or modify it).

## Global Constraints

- MATLAB is at `/Applications/MATLAB_R2025a.app/bin/matlab`; every `matlab -batch` command runs **from the repository root**, and the repo path contains spaces — always quote the binary path exactly as shown.
- Repository root: `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion` (all `cd` commands below mean this directory).
- Both drivers must run with `invz_projected` **absent from the path** under their default knobs (`show_projected_anchor = false`). Only the opt-in anchor may `addpath` it.
- New drivers are **scripts with local functions at the end** (same as the projected drivers). Do not turn them into functions.
- Error-handling conventions (from the spec, non-negotiable):
  - Configuration preflights fail loud before any compute (`Brange`, `solve_opts.mode`, `phi_ab`/`transverse_mf`, anchor demag).
  - Phase-diagram sweep absorbs **only** `invzt:bracket` per point; every other identifier that escapes `invzt_critical` rethrows.
  - Physics signals `invz:degenerateDoublet`, `invz:orderedPhase`, `invzt:a1ZeroField` become `ok = false` / a masked column in the field sweep; the q-path branch instead FAILS LOUD by design (`invzt:qpathNotPM`).
  - `parfor` sliced outputs use the val-then-assign pattern (never assign the sliced variable inside `try`).
  - The `Sigma0` validity floor is single-sourced as `getf(solve_opts, 'sigma_floor', -0.5)` (the existing `invzt_critical` option name) everywhere it appears.
- Suite expectations: CORE `runtests('invz_tensor/tests')` → **47 / 0 / 1 before Task 1, 49 / 0 / 1 from Task 1 on**; PROJECTED `runtests('invz_projected/tests')` → **143 / 0 / 19**; INTEROP `runtests('invz_tensor/tests/interop')` → **8 / 0 / 2**.
- Dispersion q-paths are **Γ-excluded** everywhere (driver example, README recipe, smoke tests): strict `[0 0 0]` is accepted by `invzt_jq_tensor` but assembles the Lorentz-cavity strict-uniform page, not the q→0⁺ limit.
- The worktree contains **unrelated pre-existing modifications** (`docs/SESSION-2026-07-16-invz-odd-mainbody.md`, `framework_revision_suggestion.txt` deletion, `jensen_1z_framework.html`, part of `invz_projected/README.html`, and `invzt_driver_review_by_Codex.md`). Never `git add -A` / `git commit -a`; stage only the exact paths named in each commit step.
- `imag(...)` of the returned susceptibilities is the positive χ'' — **no sign flip anywhere** (`chi_uniform` verified against `invz_spectra_map.m:143` and `test_invzt_rpa_parity.m`; `chi_cc_q` gains the same contract in Task 1).

---

### Task 1: Fix `invzt_chi_realaxis` — complex explicit-q response + A1-only guard

**Files:**
- Modify: `invz_tensor/invzt_chi_realaxis.m` (explicit-q accumulation ~lines 170-185; header doc; new guard after the opts-parsing block ~line 112)
- Test: `invz_tensor/tests/test_invzt_chi_realaxis.m` (two new local test functions)

**Interfaces:**
- Consumes: `pt.mode` (recorded by `invzt_solve_point.m:370` for every mode); `getf(struct, name, default)` from `invz_common/`.
- Produces: `out.chi_cc_q` `[nq, nw]` **complex** (was real-projected) — `imag()` is the positive χ''; error identifier `invzt:realaxisMode` raised for any `pt.mode ~= 'a1'`. Task 4's spectra driver relies on both.

- [ ] **Step 1: Add the failing regression test**

Append to `invz_tensor/tests/test_invzt_chi_realaxis.m` (function-based test file; new local functions are auto-collected):

```matlab
function test_qsel_explicit_q_complex_response(testCase)
% F1 regression (Codex review 2026-07-18): the explicit-q response must keep
% its imaginary (dissipative) part on the real axis. The pre-fix code
% real()-projected each site-diagonal element (a Matsubara-pattern transplant
% from invzt_gcc_lattice, where real() is a legitimate noise-clean), making
% chi_cc_q identically real and every q-path chi'' map zero.
ion = invz_ion();  T = 1.6;  B = [2 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
w = (0.05:0.05:0.55).';
qlist = [0.25 0 0; 0.1 0.2 0.3];
out = invzt_chi_realaxis(ion, T, B, pt, w, ...
    struct('odd', false, 'qsel', qlist, 'dpRng', 10, 'cache', false));
verifyFalse(testCase, isreal(out.chi_cc_q), ...
    'chi_cc_q lost its imaginary (dissipative) part');
mx = max(imag(out.chi_cc_q(:)));
verifyGreaterThan(testCase, mx, 1e-6, ...
    'no positive dissipative weight anywhere on the q list');
% chi''(w>0) >= 0 for a passive response's site-diagonal average -- a physics
% gate that needs no reimplementation of the function's internals.
verifyGreaterThan(testCase, min(imag(out.chi_cc_q(:))), -1e-6*mx, ...
    'negative chi'''' beyond eta-broadening tolerance');
end
```

- [ ] **Step 2: Run it to verify it fails against the current code**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/test_invzt_chi_realaxis/test_qsel_explicit_q_complex_response'); disp(r); assert(r.Failed == 1)"
```
Expected: FAIL on the `isreal` verification (chi_cc_q is real by construction today).

- [ ] **Step 3: Fix the accumulation (two one-line edits in `invzt_chi_realaxis.m`)**

In the explicit-q setup:
```matlab
% old:
    out.chi_cc_q = zeros(nq, nw);
% new:
    out.chi_cc_q = complex(zeros(nq, nw));
```
In the per-q accumulation loop:
```matlab
% old:
                acc = acc + real(Xq(3*(s-1)+3, 3*(s-1)+3, iq));
% new:
                acc = acc + Xq(3*(s-1)+3, 3*(s-1)+3, iq);
```

- [ ] **Step 4: Update the function header's `chi_cc_q` contract**

In the header's explicit-qvec description, change:
```matlab
%       per-q, sublattice-averaged site-diagonal cc response out.chi_cc_q
%       [nq,nw] (mirrors INVZT_GCC_LATTICE's per-q diag4 mean, evaluated one
%       q at a time rather than BZ-averaged together).
```
to:
```matlab
%       per-q, sublattice-averaged site-diagonal cc response out.chi_cc_q
%       [nq,nw], COMPLEX -- imag() is the positive chi'' (contract parity
%       with the projected INVZ_CHI_REALAXIS's chi_cc_q). The sublattice
%       diag4 mean mirrors INVZT_GCC_LATTICE's pattern but WITHOUT its
%       real(): that projection is Matsubara-axis-only and would delete the
%       whole dissipative part here (fixed 2026-07-18, Codex review F1).
```
And in the Returns block, change:
```matlab
%     out.chi_cc_q    [nq,nw]  : per-q sublattice-averaged cc response;
```
to:
```matlab
%     out.chi_cc_q    [nq,nw]  : COMPLEX per-q sublattice-averaged cc
%                                response (imag = chi'' >= 0 for w > 0);
```

- [ ] **Step 5: Run the regression test to verify it passes**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/test_invzt_chi_realaxis/test_qsel_explicit_q_complex_response'); disp(r); assertSuccess(r)"
```
Expected: PASS.

- [ ] **Step 6: Add the failing mode-guard test**

Append to the same test file:

```matlab
function test_realaxis_rejects_non_a1_point(testCase)
% invzt_chi_realaxis is the A1 scalar-Sigma continuation ONLY (LOCKED scope):
% an A2/A3 point must be rejected loudly -- its alpha/lambda/K fields are
% matrix or diagnostic objects that would otherwise produce silent garbage
% or all-NaN spectra (Codex review F3).
ion = invz_ion();  T = 1.6;  B = [2 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt2 = invzt_solve_point(ion, T, B, lat, struct('odd', false, 'mode', 'a2'));
w = (0.1:0.1:0.3).';
verifyError(testCase, ...
    @() invzt_chi_realaxis(ion, T, B, pt2, w, struct('odd', false)), ...
    'invzt:realaxisMode');
end
```

Run to verify it fails (no guard exists yet):
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/test_invzt_chi_realaxis/test_realaxis_rejects_non_a1_point'); disp(r); assert(r.Failed == 1)"
```

- [ ] **Step 7: Add the guard to `invzt_chi_realaxis.m`**

Insert after the `force_sigma0 = ...` assignment (end of the opts-parsing block) and before the `if ~isequal(odd_req, pt.odd)` check:

```matlab
ptmode = getf(pt, 'mode', 'a1');
if ~strcmp(ptmode, 'a1')
    error('invzt:realaxisMode', ['invzt_chi_realaxis is the A1 scalar-Sigma ' ...
        'continuation ONLY (LOCKED scope; see the SCOPE box in this header): ' ...
        'got pt.mode = ''%s''. Re-solve the point with mode ''a1''. Continuing ' ...
        'the A2/A3 matrix objects is an open item (README section 10).'], ptmode);
end
```

Run the guard test again — expected: PASS.

- [ ] **Step 8: Run the full CORE suite (count grows by the two new tests) + INTEROP**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```
Expected: **49 passed / 0 failed / 1 incomplete** and **8 / 0 / 2**. (The pre-existing `test_qsel_explicit_qlist_resolves_per_q` still passes: `isfinite` on complex checks both parts.)

- [ ] **Step 9: Commit**

```bash
git add invz_tensor/invzt_chi_realaxis.m invz_tensor/tests/test_invzt_chi_realaxis.m
git commit -m "fix(invzt): keep chi_cc_q complex on the real axis + A1-only guard

The explicit-q branch of invzt_chi_realaxis real()-projected the
site-diagonal response -- a Matsubara-pattern transplant from
invzt_gcc_lattice that deleted the entire dissipative part, so any q-path
chi'' map built from imag(chi_cc_q) was identically zero (Codex review
F1, confirmed at source). Drop the projection, preallocate complex, and
document imag() = positive chi'' (contract parity with the projected
invz_chi_realaxis). Also enforce the documented A1-only scope at run time
(invzt:realaxisMode on pt.mode ~= 'a1', Codex F3) instead of letting A2/A3
points produce silent garbage. Two new CORE regression tests; CORE 49/0/1.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Move `invz_peak_energy.m` to `invz_common/`

**Files:**
- Move: `invz_projected/invz_peak_energy.m` → `invz_common/invz_peak_energy.m` (pure `git mv`, zero content change)
- Test: existing suites only

**Interfaces:**
- Consumes: nothing new.
- Produces: `E = invz_peak_energy(chi, w, wmin)` reachable from `invz_common/` — `chi` is `[nw, nB]` (one spectrum per **column**), `w` the uniform frequency grid, `wmin` the lower peak-pick bound; returns `E` `[1, nB]` with censored (NaN) parabolic-refined peak positions. Task 4's spectra driver calls it exactly like this. Existing projected callers keep working unmodified (every caller already has `invz_common` on the path — spec "Verified facts").

- [ ] **Step 1: Confirm clean starting state for the paths this task touches**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
git status --short invz_projected/invz_peak_energy.m invz_common/
```
Expected: no output.

- [ ] **Step 2: Move the file**

```bash
git mv invz_projected/invz_peak_energy.m invz_common/invz_peak_energy.m
git status --short | grep peak_energy
```
Expected: `R  invz_projected/invz_peak_energy.m -> invz_common/invz_peak_energy.m`

- [ ] **Step 3: Confirm no stale path references anywhere**

```bash
grep -rn "invz_projected/invz_peak_energy" --include="*.m" --include="*.html" --include="*.md" . | grep -v "docs/superpowers" | grep -v "invzt_driver_review"
```
Expected: no output.

- [ ] **Step 4: Run the PROJECTED regression gate (the suite that exercises the moved file)**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
```
Expected: **143 / 0 / 19**.

- [ ] **Step 5: Run the CORE and INTEROP tensor suites**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```
Expected: **49 / 0 / 1** and **8 / 0 / 2**.

- [ ] **Step 6: Commit (the rename is already staged by `git mv`; stage nothing else)**

```bash
git commit -m "refactor: move invz_peak_energy to invz_common (shared by projected + tensor drivers)

Pure git mv, zero content change. All callers (invz_spectra_map,
invz_spectra_qpath, invz_chi_tensor_ref) already addpath invz_common.
PROJECTED 143/0/19, CORE 49/0/1, INTEROP 8/0/2 unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `invzt_run_phase_diagram.m`

**Files:**
- Create: `invz_tensor/invzt_run_phase_diagram.m`
- Temporary (deleted before commit): `invz_tensor/smoke_invzt_run_phase_diagram.m`

**Interfaces:**
- Consumes: `invzt_qgrid(n, conv)`, `invzt_jq_tensor(ion, g, opts)` → `lat`; `[Bc, out] = invzt_critical(ion, T, lat, Brange, opts)` (raises `invzt:bracket` BOTH for malformed `Brange` and for a valid window with no crossing — hence the preflight); `pt = invzt_solve_point(ion, T, B, lat, opts)` (fields used: `pt.crit`, `pt.converged`, `pt.Sigma0`); `[tc, used] = invzt_tc_pm_extrap(critfun, Tg)` where `critfun: T -> [crit, ok]` (throws `invzt:tcNoWindow`; provides NO try/catch around critfun — the local `invzt_crit_at` supplies it); optional opt-in: `Tc = invz_odd_zero_field(ion, opts)` from `invz_projected/` (requires `ion.demag == 0`; always builds its own legacy endpoint-inclusive mesh).
- Produces: workspace variables `Bc` `[1, numel(Ts)]`, `Tc0_proxy`, `Tc0_closed` (scalars, NaN when skipped), `phase_boundary` `[nfinite, 2]` columns `[T(K) Bc(T)]`; one figure. Local function `[c, ok] = invzt_crit_at(ion, T, B, lat, opts)` (script-local, not on the path).

- [ ] **Step 1: Write `invz_tensor/invzt_run_phase_diagram.m` with production knobs**

Exact content — the spec's Component-1 code block, copied verbatim (knob lines must stay byte-identical to the `sed` patterns in Step 2). See `docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md` §"1. `invz_tensor/invzt_run_phase_diagram.m`" — the file content is that code block, starting at the `%INVZT_RUN_PHASE_DIAGRAM` header line and ending with the closing `end` of `invzt_crit_at`. Key structural points the implementer must not lose:
- The `Brange` preflight `assert` sits between the knob block and the `show_projected_anchor` conditional (before any compute).
- The `parfor` catch absorbs only `err.identifier == 'invzt:bracket'`.
- `invzt_crit_at`'s validity line reads the floor via `getf(opts, 'sigma_floor', -0.5)`.
- The anchor comment + legend label state the projected engine's legacy-inclusive mesh (comparator, never the same quantity).

- [ ] **Step 2: Generate the reduced-size smoke copy (same directory, so the `mfilename`-relative `addpath` lines stay valid)**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
sed -e "s/^Ts     = linspace(1.0, 2.2, 25);/Ts     = [1.2 1.4 2.0 2.2];/" \
    -e "s/^gridN  = 16;/gridN  = 8;/" \
    -e "s/^dpRng  = 30;/dpRng  = 15;/" \
    -e "s/^useParallel = true;/useParallel = false;/" \
    -e "s/^show_projected_anchor = false;/show_projected_anchor = true;/" \
    invz_tensor/invzt_run_phase_diagram.m > invz_tensor/smoke_invzt_run_phase_diagram.m
grep -n "^Ts \|^gridN\|^dpRng\|^useParallel\|^show_projected_anchor" invz_tensor/smoke_invzt_run_phase_diagram.m
```
Expected grep output shows the five overridden lines.

Smoke-point rationale: T = 1.2, 1.4 K bracket and return finite Bc (~2–4.5 T); T = 2.0, 2.2 K sit above Tc0 → `invzt:bracket` → the catch path prints "no bracket" and leaves NaN. The proxy anchor finds its two lowest converged-PM points at 2.0/2.2 K. `show_projected_anchor = true` exercises the conditional `addpath` + demag assert + `invz_odd_zero_field` at nominal N = 8 / dpRng 15.

- [ ] **Step 3: Run the smoke script with assertions**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "run('invz_tensor/smoke_invzt_run_phase_diagram.m'); fprintf('finite Bc: %d of %d | Tc0_proxy = %.4f | Tc0_closed = %.4f\n', sum(isfinite(Bc)), numel(Bc), Tc0_proxy, Tc0_closed); assert(sum(isfinite(Bc)) >= 1, 'smoke: no finite Bc at all'); assert(any(isnan(Bc)), 'smoke: expected at least one no-bracket NaN'); assert(isfinite(Tc0_closed), 'smoke: projected anchor did not resolve'); assert(~isempty(phase_boundary), 'smoke: empty phase_boundary'); assert(~isempty(findall(0, 'Type', 'figure')), 'smoke: no figure created'); disp('SMOKE OK')"
```
Expected: per-point progress lines; "no bracket" notes for T = 2.0/2.2; finite Bc for T = 1.2/1.4; finite `Tc0_proxy` (≈1.5–1.7 K, coarse — value not gated; a distance-sanity flag print from `invzt_tc_pm_extrap` is reported-not-gated) or a printed skip; finite `Tc0_closed` (≈1.5–1.8 K, single coarse 8³ legacy-inclusive grid); `SMOKE OK`. Duration ~2–10 min cold cache, ~1–4 min warm.

- [ ] **Step 4: Delete the smoke copy**

```bash
rm invz_tensor/smoke_invzt_run_phase_diagram.m
git status --short invz_tensor/
```
Expected: `?? invz_tensor/invzt_run_phase_diagram.m` (plus the pre-existing `?? invz_tensor/README.html` / ` D invz_tensor/README.md` — leave those; Task 5 handles them).

- [ ] **Step 5: Commit (stage only the driver)**

```bash
git add invz_tensor/invzt_run_phase_diagram.m
git commit -m "feat(invzt): phase-diagram driver -- PM-side field-cut Bc(T) sweep

Mirrors invz_run_phase_diagram scoped to what the tensor branch supports:
invzt_critical field-cuts only (no T-cut finder exists), Brange
preflighted before the sweep (invzt:bracket doubles as the finder's
arg-validation id), invzt:bracket absorbed per point (all else that
escapes invzt_critical rethrows), val-then-assign parfor slicing,
tensor-native small-Bx Tc(0) proxy via a selective-catch invzt_crit_at
local (mirrors invz_crit_at; Sigma0 floor single-sourced as sigma_floor),
opt-in projected closed-form Tc(0) comparator (off by default; conditional
addpath + demag assert; legacy-inclusive-mesh caveat in comment + legend).
Smoke-verified at 8^3/dpRng-15 on [1.2 1.4 2.0 2.2] K: 2 bracketed,
2 no-bracket NaN, both anchors resolved.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `invzt_run_spectra.m`

**Files:**
- Create: `invz_tensor/invzt_run_spectra.m`
- Temporary (deleted before commit): `invz_tensor/smoke_invzt_run_spectra_a.m`, `invz_tensor/smoke_invzt_run_spectra_b.m`

**Interfaces:**
- Consumes: `invzt_qgrid` / `invzt_jq_tensor` → `lat`; `pt = invzt_solve_point(...)` (fields: `pt.converged`, `pt.crit`, `pt.Sigma0`); `out = invzt_chi_realaxis(ion, T, B, pt, w, opts)` with `opts.qsel = 'gamma_uniform'` → `out.chi_uniform` `[3,3,nw]`, or `opts.qsel = [nq,3] qvec` → `out.chi_cc_q` `[nq,nw]` **complex** (Task 1) — both in the positive-χ sign convention, A1-only (guarded); `E = invz_peak_energy(chi, w, wmin)` from `invz_common/` (**Task 2**).
- Produces: field-sweep branch — workspace `chipp` `[nw, nfields]`, `Epeak` `[1, nfields]`, figures; q-path branch — `chipp_q` `[nq, nw]`, `Epeak_q` `[1, nq]`, figures.

- [ ] **Step 1: Write `invz_tensor/invzt_run_spectra.m` with production knobs**

Exact content — the spec's Component-2 code block, copied verbatim (knob lines byte-identical to the Step-2/3 `sed` patterns). See `docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md` §"2. `invz_tensor/invzt_run_spectra.m`". Key structural points:
- The `solve_opts.mode` preflight (`invzt_run_spectra:mode`) and the `phi_ab`/`transverse_mf` guard sit before the lattice build.
- `solve_opts` has NO `dress` field (A3-only knob; dead under forced a1).
- `sfloor = getf(solve_opts, 'sigma_floor', -0.5)` feeds both the field-sweep validity check and the q-path check.
- The q-path branch FAILS LOUD via `invzt:qpathNotPM` (message carries `converged`, `crit`, `Sigma0`) instead of masking.
- The example `qpath` comment starts at 0.1 with the Γ-exclusion note (strict Γ = Lorentz-cavity strict-uniform page, not q→0⁺).

- [ ] **Step 2: Generate and run smoke A (field-sweep branch, with masked ordered points)**

At T = 1.4 K the boundary sits near ~2.2–2.5 T, so fields 0.5/1.5 T land on the ordered side (masked-note path) while 3.0/4.0/4.5 T are PM:

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
sed -e "s/^T      = 1.6;/T      = 1.4;/" \
    -e "s/^fields = linspace(0.3, 4.5, 40);/fields = [0.5 1.5 3.0 4.0 4.5];/" \
    -e "s/^gridN = 16;/gridN = 8;/" \
    -e "s/dpRng = 30;/dpRng = 15;/" \
    -e "s/^useParallel = true;/useParallel = false;/" \
    invz_tensor/invzt_run_spectra.m > invz_tensor/smoke_invzt_run_spectra_a.m

"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "run('invz_tensor/smoke_invzt_run_spectra_a.m'); nfin = sum(any(isfinite(chipp), 1)); fprintf('finite columns: %d of %d | finite Epeak: %d\n', nfin, size(chipp, 2), sum(isfinite(Epeak))); assert(nfin >= 1, 'smoke A: no finite spectra at all'); assert(nfin < size(chipp, 2), 'smoke A: expected at least one masked (ordered-side) field'); assert(max(chipp(:), [], 'omitnan') > 1e-6, 'smoke A: no positive spectral weight'); assert(~isempty(findall(0, 'Type', 'figure')), 'smoke A: no figure'); disp('SMOKE A OK')"
```
Expected: masked notes for the low fields, finite spectra with positive weight for the high fields, `SMOKE A OK`. With 5 fields (≤ `sliceMax`) the line-slice view is exercised; the production 40-field default exercises the colormap branch. Duration ~1–4 min warm.

- [ ] **Step 3: Generate and run smoke B (q-path branch) — SEMANTIC gates (Codex F2)**

Same reductions plus a 5-point Γ-excluded q-path and a safely-PM `Bq` (3.5 T at T = 1.4 K):

```bash
sed -e "s/^T      = 1.6;/T      = 1.4;/" \
    -e "s/^gridN = 16;/gridN = 8;/" \
    -e "s/dpRng = 30;/dpRng = 15;/" \
    -e "s/^useParallel = true;/useParallel = false;/" \
    -e "s/^qpath = \[\];/qpath = [linspace(0.1, 0.5, 5).' zeros(5, 2)];/" \
    -e "s/^Bq    = 2.0;/Bq    = 3.5;/" \
    invz_tensor/invzt_run_spectra.m > invz_tensor/smoke_invzt_run_spectra_b.m

"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "run('invz_tensor/smoke_invzt_run_spectra_b.m'); fprintf('chipp_q: [%d x %d] | max chi'''' = %.4g | finite Epeak_q: %d of %d\n', size(chipp_q, 1), size(chipp_q, 2), max(chipp_q(:)), sum(isfinite(Epeak_q)), numel(Epeak_q)); assert(max(chipp_q(:)) > 1e-6, 'smoke B: no positive spectral weight -- the F1 zero-map regression'); assert(sum(isfinite(Epeak_q)) >= 1, 'smoke B: every q column censored -- widen wq, do not weaken this gate'); assert(~isempty(findall(0, 'Type', 'figure')), 'smoke B: no figure'); disp('SMOKE B OK')"
```
Expected: the 3.5 T point passes the PM validity check, `chipp_q` is `[5 x 400]` with genuinely positive χ'' (this is the assertion that would have caught F1 — zero-weight output fails it), at least one uncensored peak, `SMOKE B OK`. If every `Epeak_q` is censored, widen `wq` in the smoke copy and re-run — do not weaken the gate.

- [ ] **Step 4: Delete both smoke copies**

```bash
rm invz_tensor/smoke_invzt_run_spectra_a.m invz_tensor/smoke_invzt_run_spectra_b.m
git status --short invz_tensor/
```
Expected: `?? invz_tensor/invzt_run_spectra.m` (plus the pre-existing README items handled in Task 5).

- [ ] **Step 5: Commit (stage only the driver)**

```bash
git add invz_tensor/invzt_run_spectra.m
git commit -m "feat(invzt): spectra driver -- PM-side chi''_cc field sweep + q-path view

Mirrors invz_run_spectra scoped to the tensor branch: A1-only real-axis
continuation (mode preflight at the driver + invzt:realaxisMode at the
callee; no dress knob -- A3-only), PM-side only (ordered/invalid fields
masked via the three-part sample-validity rule with the single-sourced
sigma_floor; physics signals caught selectively, all else rethrows),
q-path view fails loud (invzt:qpathNotPM) since one invalid Bq
invalidates the whole product, theta_c/phi_ab tilt knobs ported, peak
pick via the shared invz_common/invz_peak_energy with nonzero peak_wmin
default. Dispersion q-paths Gamma-excluded by convention. Smoke-verified
both branches at 8^3/dpRng-15, T = 1.4 K, with semantic gates (positive
spectral weight + uncensored peak on the q-path).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: README §2 update + final gates

**Files:**
- Modify: `invz_tensor/README.html` (§2 callout + recipe intros + q-path recipe Γ-exclusion + two stale §-refs)
- Commit also (pre-existing, related, currently uncommitted): `invz_tensor/README.html` (whole file is untracked — the earlier session's MD→HTML conversion), `invz_tensor/README.md` (its deletion), `invz_tensor/SESSION-2026-07-18-invz-tensor-full.md` (one-line `.md`→`.html` reference fix)
- Do **NOT** commit: `invz_projected/README.html` (pre-existing user modifications + one link fix — flag to the user instead), `docs/SESSION-2026-07-16-invz-odd-mainbody.md`, `jensen_1z_framework.html`, `framework_revision_suggestion.txt` deletion, `Data/*.fig`, `LiHoF4_soft_mode_thermodynamics_v2.html`, `.claude/`, `invzt_driver_review_by_Codex.md`

**Interfaces:**
- Consumes: the two driver filenames from Tasks 3–4 — referenced verbatim in the HTML.
- Produces: documentation only.

- [ ] **Step 1: Replace §2's "No turnkey driver" callout**

Find the callout directly under `<h2 id="recipes">` (single `<div>` on one line):

```html
<div class="callout caution"><span class="tag">No turnkey driver</span> <code>invz_projected</code> ships <code>invz_run_phase_diagram.m</code> and <code>invz_run_spectra.m</code> &mdash; scripts you just run. <code>invz_tensor</code> does <strong>not</strong> have equivalents yet: it ships the point-solver primitives (<code>invzt_solve_point</code>, <code>invzt_critical</code>, <code>invzt_chi_realaxis</code>) that a driver would call, but no driver script that sweeps them over T/B and plots the result. The two recipes below are exactly what such a driver would do &mdash; write them to a script under <code>invz_tensor/</code> (or the repo root) if you want to reuse them; <code>invzt_run_ladder.m</code> is <strong>not</strong> that driver &mdash; it walks the A0&ndash;A4 basis/mode ladder at one fixed validation point (&sect;4/&sect;9 A4), it does not sweep T or B.</div>
```

Replace with (one line):

```html
<div class="callout note"><span class="tag">Drivers</span> Two runnable driver scripts ship with the module &mdash; <code>invzt_run_phase_diagram.m</code> (PM-side field-cut boundary B<sub>c</sub>(T) with zero-field anchors) and <code>invzt_run_spectra.m</code> (field-sweep &chi;&Prime;<sub>cc</sub> spectra at q&rarr;0, or a q-path dispersion view via the <code>qpath</code> knob) &mdash; the tensor analogues of the projected <code>invz_run_phase_diagram.m</code> / <code>invz_run_spectra.m</code>, with knobs at the top of each script (including the tensor-only <code>mode</code>/<code>nlevels</code> axes; the spectra driver requires <code>mode = 'a1'</code>, the only rung with a real-axis continuation, and dispersion q-paths are &Gamma;-excluded &mdash; strict q = 0 computes the Lorentz-cavity strict-uniform observable, not the q&rarr;0<sup>+</sup> limit). Both are honestly scoped to PM-side physics: there is no ordered/FM tensor solve and no temperature-cut finder (<code>invzt_critical_T</code>), so the near-vertical region of the boundary approaching T<sub>c</sub>(0) is left blank by design (&sect;10 "Open items"). <code>invzt_run_ladder.m</code> is <strong>not</strong> a sweep driver &mdash; it walks the A0&ndash;A4 basis/mode ladder at one fixed validation point (&sect;5/&sect;9 A4). The recipes below are what these drivers do under the hood &mdash; the starting point for custom sweeps.</div>
```

- [ ] **Step 2: Point the two recipe subsections at the drivers + fix stale refs**

In §2.1, replace:

```html
<p><code>invzt_critical</code> (&sect;4) finds one Bc at one fixed T by bisection. Build the boundary curve by calling it on a temperature grid, exactly as <code>invz_run_phase_diagram</code> does one level up for the projected module:</p>
```

with:

```html
<p><code>invzt_critical</code> (&sect;5) finds one Bc at one fixed T by bisection. <code>invzt_run_phase_diagram.m</code> packages the loop below (plus parfor, the small-B<sub>x</sub> proxy anchor, and an opt-in projected closed-form comparator); the inline version, for custom sweeps:</p>
```

In §2.2, replace:

```html
<p>Converge a point with <code>invzt_solve_point</code>, then continue it to the real axis with <code>invzt_chi_realaxis</code> over whatever frequency grid you want (mirrors the projected <code>invz_chi_realaxis</code> call inside <code>invz_run_spectra</code>/<code>invz_spectra_map</code>):</p>
```

with:

```html
<p>Converge a point with <code>invzt_solve_point</code>, then continue it to the real axis with <code>invzt_chi_realaxis</code> over whatever frequency grid you want &mdash; <code>invzt_run_spectra.m</code> packages exactly this over a field sweep (with ordered-side masking and the shared <code>invz_peak_energy</code> peak pick):</p>
```

Make the §2.2 q-path recipe Γ-excluded: in its `<pre>` block, replace

```
qpath = [linspace(0, 0.5, 41).' zeros(41,2)];   % Gamma -> (0.5,0,0)
```

with

```
qpath = [linspace(0.1, 0.5, 41).' zeros(41,2)]; % finite-q -> (0.5,0,0); keep Gamma-excluded
```

(if the exact spacing/comment differs, locate the recipe's `qpath =` line by grep and apply the same 0→0.1 start + Γ-exclusion comment).

Fix the two stale `&sect;9 "Open items"` references inside §2 (Open items is §10 in the final numbering) — they are the only two `&sect;9 "Open items"` occurrences in the file; change both to `&sect;10 "Open items"`.

- [ ] **Step 3: Check §2 renders and its structure is consistent**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
python3 -c "
import re
s = open('invz_tensor/README.html').read()
assert 'No turnkey driver' not in s, 'old callout still present'
assert 'invzt_run_phase_diagram.m' in s and 'invzt_run_spectra.m' in s
assert '&sect;9 \"Open items\"' not in s, 'stale section-9 Open-items ref remains'
for tag in ['div','table','ul','ol','pre']:
    o = len(re.findall(r'<'+tag+r'[ >]', s)); c = len(re.findall(r'</'+tag+'>', s))
    assert o == c, f'{tag}: {o} open vs {c} close'
print('README checks OK')
"
open invz_tensor/README.html
```
Expected: `README checks OK`; visual check of §2 in the browser.

- [ ] **Step 4: Run the final suite gates**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
```
Expected: **49 / 0 / 1** and **143 / 0 / 19**.

- [ ] **Step 5: Commit the README work (deliberately includes the earlier session's MD→HTML conversion — stage exactly these three paths)**

```bash
git add invz_tensor/README.html invz_tensor/README.md invz_tensor/SESSION-2026-07-18-invz-tensor-full.md
git diff --cached --name-status
git commit -m "docs(invzt): README.md -> README.html (projected-style) + section-2 driver quickstart

Converts the module reference to the same HTML style as
invz_projected/README.html (earlier-session conversion, committed here),
adds section 2 'Getting a phase diagram or a susceptibility spectrum'
pointing at the new invzt_run_phase_diagram.m / invzt_run_spectra.m
drivers with the inline recipes kept as the under-the-hood explanation
(q-path recipe made Gamma-excluded per the F8 convention decision), and
fixes the SESSION note's README.md reference.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
The `git diff --cached --name-status` check must list exactly: `A	invz_tensor/README.html`, `D	invz_tensor/README.md`, `M	invz_tensor/SESSION-2026-07-18-invz-tensor-full.md`. If anything else appears staged, unstage it before committing.

- [ ] **Step 6: Flag the deliberately-uncommitted files to the user**

Report at handoff: (a) `invz_projected/README.html` carries one line from this work (the `invz_tensor/README.md` → `README.html` link fix in §7) on top of pre-existing modifications that predate this task — left uncommitted so those user edits aren't swept into an unrelated commit; (b) `invzt_driver_review_by_Codex.md` is the user's review file, untouched and uncommitted.

---

## Execution notes

- Task order: 1 → 2 → 3 → 4 → 5, strictly. Task 4 requires Task 1 (complex `chi_cc_q` + guard) and Task 2 (`invz_peak_energy` on the `invz_common` path). Task 3 is independent of both but keeping the order lets its smoke run warm the 8³/dpRng-15 caches Task 4 reuses.
- All smoke runs use the **same-directory temp-copy** pattern (`sed` a reduced-knob copy next to the original, run, delete): the drivers' `mfilename('fullpath')`-relative `addpath` lines only work from inside `invz_tensor/`, and the production knobs (16³ grid, 25-point Ts, 40-field sweep) are hours of compute that stay "left to the user" per repo precedent.
- Production-knob defaults are deliberately committed **unexecuted at full scale** — identical precedent to `invz_run_phase_diagram.m`, whose committed default sweep is also an hours-long run.
- Smoke gates are SEMANTIC where it matters (Codex F2): positive spectral weight and an uncensored peak on the q-path, not merely "finite" (zeros are finite).
