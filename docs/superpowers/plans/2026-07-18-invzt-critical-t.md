# Tensor T-Cut Finder (`invzt_critical_T`) + Two-Regime Driver Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the tensor branch the projected branch's fixed-field critical-temperature search — `invzt_critical_T` (Tc at fixed B) — and rewire `invzt_run_phase_diagram.m` into the projected-style two-regime sweep (`Ts` field cuts + `Bs` temperature cuts).

**Architecture:** A pure decision core (`invzt_tc_pick`) classifies valid-sample votes (fixing the projected classifier's inherited exact-zero and re-entrant-lower-leg gaps), consumed by `invzt_critical_T` (valid-samples-only T-grid, highest-T ordered→para root, regula-falsi refinement, hard explicit window OR adaptive Tc0-anchored sliding). Shared helpers: `invz_refine_crossing.m` moves to `invz_common/`; the driver's local `invzt_crit_at` is promoted to a module file. The driver gains a flat two-kind `parfor` with namespaced per-finder tolerances.

**Tech Stack:** MATLAB R2025a; existing primitives `invzt_solve_point`, `invzt_critical`, `invzt_tc_pm_extrap`, `invzt_qgrid`, `invzt_jq_tensor`, `invzt_str`; shared engine in `invz_common/`.

**Spec (source of record):** `docs/superpowers/specs/2026-07-18-invzt-critical-t-design.md` **rev 4** (post `critical-t-review_by_Codex.md` rounds 1–3 — dispositions F1–F8, R1–R6, and T1–T3 baked in). All code blocks below are copied from it verbatim; if they diverge, the spec governs.

## Global Constraints

- MATLAB is at `/Applications/MATLAB_R2025a.app/bin/matlab`; every `matlab -batch` command runs **from the repository root**; the repo path contains spaces — always quote the binary path exactly as shown. Generous timeouts (600000 ms). All durations quoted below are estimates, never acceptance criteria.
- Repository root: `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion` (all `cd` commands below mean this directory).
- Single-test selector form that works on this install: `runtests('<file>.m', 'ProcedureName', '<name>')` — the `dir/file/testname` slash form is INVALID here.
- Suite expectations: CORE `runtests('invz_tensor/tests')` → **50 / 0 / 1 before Task 2, 55 / 0 / 3 from Task 2 on**; INTEROP `runtests('invz_tensor/tests/interop')` → **8 / 0 / 2**; PROJECTED `runtests('invz_projected/tests')` → **143 / 0 / 19** (frozen — no new tests go in the projected suite).
- **Byte-parity invariant:** the committed `invz_tensor/invzt_run_phase_diagram.m` must be byte-identical to the spec rev 4's Component-5 fenced ```matlab block (Task 4 marks the drivers spec's Component-1 block SUPERSEDED).
- Error policy (spec, non-negotiable): driver preflights fail loud before compute (`invzt_run_phase_diagram:Brange/:Ts/:Bs/:Twin/:Btol/:Ttol/:Twidth/:Tgridstep/:tcAnchor`, all `invzt_str`-formatted — never `mat2str`, which throws on structs while building the message); the sweep absorbs **only** per-point `invzt:bracket`; `invzt_crit_at`'s selective catch stays exactly `{invz:degenerateDoublet, invz:orderedPhase, invzt:a1ZeroField}` → `ok=false`, all else rethrows; `parfor` sliced outputs use val-then-assign.
- `invzt_critical_T`: explicit `opts.window` is a HARD bound (one pass, no sliding, and the exit happens BEFORE any slide so errors report the window actually sampled — R2); adaptive mode requires a finite positive scalar `opts.Tc0` (`invzt:tcAnchor`), makes up to 9 window attempts, terminates the no-valid `grow` path at the 0.02 K floor (R3), and accumulates re-entrance indicators across attempts (`nseen` — R5); `width`/`gridstep`/`tol` are validated finite positive real scalars (`invzt:tcOpts`); caller-supplied `Sigma_seed` is stripped (proven by a length-MATCHED poisonous seed — a wrong-length seed is silently ignored by the solver, R1).
- **Worktree hygiene (per-task checkpoint, no hard-coded dirt list):** start EVERY task with `git status --short`, note the unrelated dirty/untracked paths present at that moment, and stage only the exact paths named in the task's commit step. The worktree carries the user's own WIP — at plan time this includes `invz_tensor/invzt_run_spectra.m` (in-progress production knobs) — which must NEVER be staged, reverted, or edited. Never `git add -A` / `git commit -a`.
- MATLAB test-fixture trap: `struct('f', [a b])` creates a struct ARRAY — malformed-opts fixtures in tests are built by field assignment (`o = struct(); o.window = [1.8 1.0];`).

---

### Task 1: Move `invz_refine_crossing.m` to `invz_common/`

**Files:**
- Move: `invz_projected/invz_refine_crossing.m` → `invz_common/invz_refine_crossing.m` (pure `git mv`, zero content change)
- Test: existing suites only (direct unit coverage of the helper arrives in Task 2's CORE test file — the projected suite count stays frozen)

**Interfaces:**
- Consumes: nothing new.
- Produces: `bx = invz_refine_crossing(f, Ba, ca, Bb, cb, tol)` reachable from `invz_common/` — regula-falsi between CONVERGED bracket ends; `f` is a closure `x -> [value, ok]`; a non-converged interior sample triggers a midpoint retry, then a linear-interpolation fallback. Task 2's `invzt_critical_T` calls it exactly like the projected finders do. Existing projected callers (`invz_critical.m:40,51`, `invz_critical_T.m:67`) keep working unmodified — every test/driver reaching them already addpaths `invz_common`; the fast direct exercise after the move is the interop critical-parity test (the substantive projected T-cut tests are `INVZ_SLOW`-gated).

- [ ] **Step 1: Worktree checkpoint + clean state for the touched paths**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
git status --short
git status --short invz_projected/invz_refine_crossing.m invz_common/
```
Expected: the second command prints nothing (neither touched path has pending changes); note the unrelated dirty paths from the first for later staging discipline.

- [ ] **Step 2: Move the file**

```bash
git mv invz_projected/invz_refine_crossing.m invz_common/invz_refine_crossing.m
git status --short | grep refine_crossing
```
Expected: `R  invz_projected/invz_refine_crossing.m -> invz_common/invz_refine_crossing.m`

- [ ] **Step 3: Confirm no stale path references**

```bash
grep -rn "invz_projected/invz_refine_crossing" --include="*.m" --include="*.html" --include="*.md" . | grep -v "docs/superpowers" | grep -v "review_by_Codex"
```
Expected: no output.

- [ ] **Step 4: Run the PROJECTED gate**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
```
Expected: **143 passed / 0 failed / 19 incomplete**.

- [ ] **Step 5: Run CORE and INTEROP**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```
Expected: **50 / 0 / 1** and **8 / 0 / 2** (the interop critical-parity test is the fast path through the moved helper).

- [ ] **Step 6: Commit (the rename is already staged by `git mv`; stage nothing else)**

```bash
git commit -m "refactor: move invz_refine_crossing to invz_common (shared by projected + tensor finders)

Pure git mv, zero content change. Its two callers (invz_critical,
invz_critical_T) and every test/driver reaching them already addpath
invz_common (same call-graph argument as the invz_peak_energy move);
fast post-move exercise via the interop critical-parity test, direct
unit coverage lands with the tensor T-cut test file (CORE, keeping the
frozen PROJECTED count untouched). PROJECTED 143/0/19, CORE 50/0/1,
INTEROP 8/0/2 unchanged.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `invzt_tc_pick` + `invzt_crit_at` + `invzt_critical_T` + tests

**Files:**
- Create: `invz_tensor/invzt_tc_pick.m` (spec Component 3 — the pure decision core)
- Create: `invz_tensor/invzt_crit_at.m` (spec Component 2 — the promoted sampler)
- Create: `invz_tensor/invzt_critical_T.m` (spec Component 4)
- Modify: `invz_tensor/invzt_str.m` (spec Component 2b — harden the formatter fallthrough: raw `mat2str` throws on structs, masking the intended error id)
- Create: `invz_tensor/tests/test_invzt_critical_T.m` (5 fast + 2 slow test functions)
- NOT touched: `invz_tensor/invzt_run_phase_diagram.m` — its identical local `invzt_crit_at` legally shadows the new module file inside the script until Task 3 removes it (MATLAB scoping: local functions win); behavior is identical either way because the code is identical.

**Interfaces:**
- Consumes: `invz_refine_crossing` from `invz_common/` (Task 1); `invzt_solve_point` (fields `pt.crit`, `pt.converged`, `pt.Sigma0`; validates `opts.mode` and the zero-transverse-field guard BEFORE touching `lat` — the sampler tests exploit this with dummy structs; uses a caller seed ONLY when `numel(opts.Sigma_seed) == nwn`); `invz_matsubara(T, Ecut)` (for the length-matched seed fixture); `invz_field_vec`, `getf` from `invz_common/`; `invzt_str`.
- Produces: `[act, ka, kb, ncross] = invzt_tc_pick(cv)` — `act ∈ {'zero','bracket','up','down'}`, indices into the voter vector, re-entrance count. `[c, ok] = invzt_crit_at(ion, T, B, lat, opts)` — validity-only `ok`. `[tc, out] = invzt_critical_T(ion, B, lat, opts)` — opts `window` (hard) | `Tc0` (adaptive, required) | `width`/`gridstep`/`tol`; errors `invzt:tcOpts`/`invzt:tcWindow`/`invzt:tcAnchor`/`invzt:bracket`; warns `invzt:multipleCrossings`; `out` fields `.Tg/.c/.ok/.window/.ncross/.B`. Task 3's driver calls all three.

- [ ] **Step 1: Worktree checkpoint, then write the test file**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && git status --short
```

Create `invz_tensor/tests/test_invzt_critical_T.m`:

```matlab
function tests = test_invzt_critical_T
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

% ---------- fast: pure crossing/slide policy (synthetic votes, no solves) -----

function test_pick_crossing_policy(testCase)
% Simple bracket
[act, ka, kb, nc] = invzt_tc_pick([-1 1]);
verifyEqual(testCase, act, 'bracket');
verifyEqual(testCase, [ka kb], [1 2]);
verifyEqual(testCase, nc, 1);
% Multiple crossings: HIGHEST ordered->para pair wins, re-entrance counted
[act, ka, kb, nc] = invzt_tc_pick([1 -1 1]);
verifyEqual(testCase, act, 'bracket');
verifyEqual(testCase, [ka kb], [2 3]);
verifyEqual(testCase, nc, 2);
% All-PM window: boundary below
act = invzt_tc_pick([1 2 3]);
verifyEqual(testCase, act, 'down');
% All-ordered window: boundary above
act = invzt_tc_pick([-1 -2]);
verifyEqual(testCase, act, 'up');
% Re-entrant LOWER leg (para below, ordered above): the physical high-T side
% is paramagnetic, so the highest ordered->para crossing is ABOVE -- must
% move up, not give up (the projected classifier's inherited gap).
act = invzt_tc_pick([1 1 -1 -1]);
verifyEqual(testCase, act, 'up');
% Singletons
verifyEqual(testCase, invzt_tc_pick(1),  'down');
verifyEqual(testCase, invzt_tc_pick(-1), 'up');
end

function test_pick_exact_zero(testCase)
% A sampled crit == 0 IS the boundary (the projected classifier mis-reads
% [-1, 0, 1] as two sign changes with no returnable crossing).
[act, ka, ~, nc] = invzt_tc_pick([-1 0 1]);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 2);
verifyEqual(testCase, nc, 1);                    % one boundary, not two
% Zero RUN counted once
[act, ka, ~, nc] = invzt_tc_pick([-1 0 0 1]);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 3);
verifyEqual(testCase, nc, 1);
% Zero with ordered ABOVE it: the true highest boundary is above -> up
act = invzt_tc_pick([1 0 -1]);
verifyEqual(testCase, act, 'up');
% Zero ABOVE a strict crossing: the zero is the higher root
[act, ka, ~, nc] = invzt_tc_pick([-1 1 -1 0 1]);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 4);
verifyEqual(testCase, nc, 3);
% Lone exactly-critical voter
[act, ka] = invzt_tc_pick(0);
verifyEqual(testCase, act, 'zero');
verifyEqual(testCase, ka, 1);
end

% ---------- fast: input-contract validation (guards fire pre-compute) --------

function test_anchor_and_window_validation(testCase)
% lat is a dummy struct throughout -- every guard fires before any compute.
% Malformed-opts fixtures are built by FIELD ASSIGNMENT: the struct('f',[a b])
% constructor would create a struct ARRAY, not a field holding a vector.
% First, the hardened safe formatter itself: raw mat2str throws on structs
% (masking the intended error id); invzt_str must return a placeholder.
verifyEqual(testCase, invzt_str(struct()), '<1x1 struct>');
ion = invz_ion();
% invzt:tcAnchor -- adaptive mode needs a USABLE Tc0
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), struct()), 'invzt:tcAnchor');
o = struct(); o.Tc0 = [1 2];      % vector: must NOT surface MATLAB:nonLogicalConditional
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
o = struct(); o.Tc0 = 1 + 1i;     % complex: isfinite(1+1i) is true -- needs isreal
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
o = struct(); o.Tc0 = -1;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
o = struct(); o.Tc0 = 0.01;       % at/below the 0.02 K solve floor
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcAnchor');
% invzt:tcWindow -- explicit window malformed / nonnumeric / below-floor
o = struct(); o.window = [1.8 1.0];
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcWindow');
o = struct(); o.window = struct();  % nonnumeric: message building must not mat2str-throw
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcWindow');
o = struct(); o.window = [0.005 0.015];
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcWindow');
% invzt:tcOpts -- numerical controls must be finite positive real scalars
o = struct(); o.Tc0 = 1.5; o.gridstep = 0;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcOpts');
o = struct(); o.Tc0 = 1.5; o.tol = -1;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcOpts');
o = struct(); o.Tc0 = 1.5; o.width = Inf;
verifyError(testCase, @() invzt_critical_T(ion, 1.0, struct(), o), 'invzt:tcOpts');
% CHEAP bracket diagnostics (third review T1/T3): with B = 0 every sample is
% invalidated by the PRE-LATTICE invzt:a1ZeroField guard (absorbed by
% invzt_crit_at), so whole windows run in milliseconds on a dummy lat.
% (a) Hard window, one pass: the invzt:bracket MESSAGE must report EXACTLY
% the user's window -- never a slide-mutated pair. NB explicit try/catch:
% verifyError does NOT return the caught MException (its outputs are the
% evaluated function's outputs -- <missing> on a throw).
o = struct(); o.window = [1.1 1.4];
didThrow = false;
try
    invzt_critical_T(ion, 0, struct(), o);
catch ME
    didThrow = true;
    verifyEqual(testCase, ME.identifier, 'invzt:bracket');
    verifySubstring(testCase, ME.message, '[1.100, 1.400]');
end
verifyTrue(testCase, didThrow, 'expected invzt:bracket was not thrown');
% (b) Adaptive no-valid grow path terminates at the 0.02 K floor (R3) --
% Tc0 just above the floor makes the first window floor-clamped already, so
% the scan must stop instead of re-sampling the identical grid.
o = struct(); o.Tc0 = 0.05;
didThrow = false;
try
    invzt_critical_T(ion, 0, struct(), o);
catch ME
    didThrow = true;
    verifyEqual(testCase, ME.identifier, 'invzt:bracket');
end
verifyTrue(testCase, didThrow, 'expected invzt:bracket was not thrown');
end

% ---------- fast: sampler contract (pre-lattice guard paths, dummy lat) ------

function test_crit_at_contract(testCase)
% Physics signals absorbed as invalid samples; config errors rethrow. Both
% paths fire inside invzt_solve_point BEFORE any lattice access, so a dummy
% lat keeps these at millisecond cost.
ion = invz_ion();
[c, ok] = invzt_crit_at(ion, 1.6, [0 0 0], struct(), struct());   % zero transverse field
verifyFalse(testCase, ok);
verifyTrue(testCase, isnan(c));                                   % invzt:a1ZeroField absorbed
o = struct(); o.mode = 'bogus';
verifyError(testCase, @() invzt_crit_at(ion, 1.6, [2 0 0], struct(), o), ...
    'invzt:mode');                                                % misconfiguration rethrows
end

% ---------- fast: the moved invz_common helper, directly ----------------------

function test_refine_crossing_helper(testCase)
% Regula-falsi between converged bracket ends; non-converged interior ->
% midpoint retry; total interior failure -> linear-interpolation fallback.
froot = 1.3;
f1 = @(x) deal(2*(x - froot), true);
bx = invz_refine_crossing(f1, 1.0, 2*(1.0 - froot), 1.8, 2*(1.8 - froot), 1e-4);
verifyEqual(testCase, bx, froot, 'AbsTol', 1e-3);
% Interior dead zone straddling the root: midpoint retries + the final
% linear interpolation must still land the root (f is linear, so the
% fallback is exact up to the surviving bracket).
f2 = @(x) deal(2*(x - froot), ~(x > 1.25 && x < 1.35));
bx = invz_refine_crossing(f2, 1.0, 2*(1.0 - froot), 1.8, 2*(1.8 - froot), 1e-4);
verifyEqual(testCase, bx, froot, 'AbsTol', 0.06);
end

% ---------- slow: integration gates (8^3, warm caches) ------------------------

function test_tcut_matches_field_cut_slow(testCase)
% Crossing consistency (the projected branch's own validation pattern): at
% T* = 1.4 K find B* with the validated field-cut finder, then the T-cut at
% B* must return T*. Also proves the Sigma_seed strip (a length-MATCHED,
% poisonous-VALUED seed is passed -- unstripped, the matching endpoint
% consumes the NaNs and turns invalid) and the out-struct contract.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
g   = invzt_qgrid(8, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
o   = struct('odd', false);                             % lighter; round-trip is self-consistent
Bstar = invzt_critical(ion, 1.4, lat, [0.05 6], o);
oT = o;
oT.window = [1.0 1.8];                                  % HARD window containing T* = 1.4
% Deliberately POISONOUS seed, length-MATCHED to the top endpoint's Matsubara
% count (Ecut default 40): a wrong-length seed would be silently ignored by
% invzt_solve_point's numel(Sigma_seed) == nwn guard and prove nothing about
% the strip (second review R1). Unstripped, the 1.8 K endpoint consumes the
% NaN seed and turns invalid; stripped, it starts from the zero seed and
% stays a converged PM voter.
wnHi = invz_matsubara(oT.window(2), 40);
oT.Sigma_seed = nan(numel(wnHi), 1);
[tc, out] = invzt_critical_T(ion, Bstar, lat, oT);
verifyEqual(testCase, tc, 1.4, 'AbsTol', 0.05);
verifyTrue(testCase, all(isfield(out, {'Tg', 'c', 'ok', 'window', 'ncross', 'B'})));
verifyEqual(testCase, numel(out.c), numel(out.Tg));
verifyEqual(testCase, numel(out.ok), numel(out.Tg));
verifyGreaterThanOrEqual(testCase, out.ncross, 1);
verifyTrue(testCase, out.ok(end), ...
    'top-endpoint sample invalid: the poisonous Sigma_seed was not stripped');
% sigma_floor threading + rejection (second review R6): an impossible floor
% must invalidate a converged point whose crit is still finite (one solve).
oF = o;  oF.sigma_floor = Inf;
[cF, okF] = invzt_crit_at(ion, 1.8, [Bstar 0 0], lat, oF);
verifyTrue(testCase, isfinite(cF));
verifyFalse(testCase, okF);
fprintf('T-cut round-trip: Bc(1.4 K) = %.4f T -> Tc = %.4f K (ncross = %d)\n', ...
        Bstar, tc, out.ncross);
end

function test_tcut_odd_on_slow(testCase)
% The capability the projected T-cut lacks: bracketing with ODD on (the
% tensor A1 solver converges metastable PM points inside the ordered phase,
% so valid votes exist on both sides of the crossing -- ODD-LOG T2 records
% that the projected finder cannot bracket there at all).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
g   = invzt_qgrid(8, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
oT = struct('odd', true);
oT.window = [1.2 1.7];
tc = invzt_critical_T(ion, 1.5, lat, oT);
verifyGreaterThan(testCase, tc, 1.35);   % loose physical band at 8^3: the odd-on
verifyLessThan(testCase, tc, 1.60);      % boundary has Bc(1.4 K) = 1.916 T > 1.5 T
fprintf('odd-on T-cut: Tc(1.5 T) = %.4f K\n', tc);
% Hard window entirely inside the ordered phase: ONE pass, no sliding ->
% invzt:bracket with the widen-the-window message (F4's contract end-to-end).
oT2 = struct('odd', true);
oT2.window = [0.25 0.45];
% Physical deep-ordered hard window: one pass, no sliding, and the error
% reports the window the user asked for (R2). Explicit try/catch: verifyError
% does NOT return the caught MException (third review T1).
didThrow = false;
try
    invzt_critical_T(ion, 1.5, lat, oT2);
catch ME
    didThrow = true;
    verifyEqual(testCase, ME.identifier, 'invzt:bracket');
    verifySubstring(testCase, ME.message, '[0.250, 0.450]');
end
verifyTrue(testCase, didThrow, 'expected invzt:bracket was not thrown');
end
```

- [ ] **Step 2: Run the fast tests to verify they fail where they should**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/test_invzt_critical_T.m'); disp(r); assert(sum([r.Failed]) == 4 && sum([r.Passed]) == 1)"
```
Expected: **4 failed, 1 passed, 2 incomplete**. The two `invzt_tc_pick` tests, the validation test (its first assert hits the UN-hardened `invzt_str`, which still throws on structs), and the sampler test all fail (`MATLAB:UndefinedFunction` / thrown formatter ≠ expected identifiers); `test_refine_crossing_helper` PASSES already — `invz_refine_crossing` exists since Task 1, so that test is new direct coverage, not a TDD RED.

- [ ] **Step 3: Create `invz_tensor/invzt_tc_pick.m`**

Exact content — the spec rev 4's Component-3 code block, verbatim (from `function [act, ka, kb, ncross] = invzt_tc_pick(cv)` to its final `end`). Extract the fenced block with sed/awk and `diff` against your written file — expect byte-identical.

- [ ] **Step 4: Create `invz_tensor/invzt_crit_at.m` and harden `invz_tensor/invzt_str.m`**

`invzt_crit_at.m`: exact content — the spec rev 4's Component-2 code block, verbatim. Same extraction + `diff` verification.

`invzt_str.m`: replace its body with the spec rev 4's Component-2b code block, verbatim (the change is behavior-additive: char/string/numeric/logical inputs format exactly as before; only inputs that previously THREW — structs/cells/objects — now return a `'<1x1 struct>'`-style placeholder).

- [ ] **Step 5: Create `invz_tensor/invzt_critical_T.m`**

Exact content — the spec rev 4's Component-4 code block, verbatim (~190 lines). Same extraction + `diff` verification. Load-bearing structure to sanity-check after transcription: control-knob validation (`invzt:tcOpts`) BEFORE the `Sigma_seed` strip and `f` closure; hard-window branch (`invzt:tcWindow`, incl. the below-floor check) vs adaptive branch (`invzt:tcAnchor` requires finite positive scalar `Tc0 > Tmin`); `nattempt = 1` for hard windows vs 9 adaptive attempts; the floor-collapse `break`; `winS` (last SAMPLED window) captured each attempt AFTER the floor-collapse guard (T2 -- a collapse reports the previous really-sampled window) and used in ALL reporting; `nseen` accumulating `ncross` across attempts and feeding the warning + `out.ncross` (R5); the `invzt_tc_pick` return switch (`zero` returns `Tv(ka)` directly, `bracket` refines); `if hardwin, break; end` BEFORE the slide switch (R2 — the one-pass contract keeps the reported window unmutated); the `grow` branch's floor termination (R3); the two distinct `invzt:bracket` messages built from `winS`.

- [ ] **Step 6: Run the fast tests to verify they pass**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/test_invzt_critical_T.m'); disp(r); assertSuccess(r)"
```
Expected: **5 passed / 0 failed / 2 incomplete**.

- [ ] **Step 7: Run the slow integration gates once**

```bash
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/test_invzt_critical_T.m'); disp(r); assertSuccess(r)"
```
Expected: **7 passed / 0 failed / 0 incomplete**, with both printed lines. `Bc(1.4 K)` with `odd=false` should land noticeably ABOVE the odd-on 1.916 T (ODD suppresses ordering) — value reported, not gated; the gated quantities are the round-trip (`Tc` within 0.05 K of 1.4) and the odd-on band. Estimated ~10–20 min warm cache (≈ 40 + 20 point solves at 8³). If a gate fails, STOP and report with the printed values — do not widen tolerances or windows.

- [ ] **Step 8: Run full CORE + INTEROP**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```
Expected: **55 / 0 / 3** (50 + the five fast tests; the two slow ones join the A4 ladder as filtered-incomplete) and **8 / 0 / 2**.

- [ ] **Step 9: Commit**

```bash
git add invz_tensor/invzt_tc_pick.m invz_tensor/invzt_crit_at.m invz_tensor/invzt_critical_T.m invz_tensor/invzt_str.m invz_tensor/tests/test_invzt_critical_T.m
git commit -m "feat(invzt): T-cut finder invzt_critical_T + pure invzt_tc_pick core + shared invzt_crit_at

Transplant of the projected invz_critical_T's rugged-boundary algorithm
(valid-samples-only T-grid voting -- near the boundary critical slowing
down makes invalid samples poisonous voters -- highest-T ordered->para
root, regula-falsi via the shared invz_refine_crossing) with the
inherited classifier gaps CLOSED in a pure, synthetically-tested
decision core (invzt_tc_pick): exact crit==0 samples are recognized as
roots (projected mis-reads [-1,0,1] as unbracketable), zero runs count
once, and a window whose top voter is ordered always moves UP (incl.
the re-entrant lower leg the projected code abandons). Explicit window
is a HARD bound (one pass); adaptive mode (up to 9 attempts) REQUIRES a
finite positive scalar opts.Tc0 (invzt:tcAnchor -- no tensor zero-field
closed form). Strict invzt:tcOpts/tcWindow validation with invzt_str
formatting (mat2str throws on structs); Sigma_seed stripped (Matsubara
length is T-dependent, strip proven by a length-MATCHED poisonous seed in the slow
round-trip). invzt_crit_at promoted from the driver's local (validity-
only ok; consumers apply their own phase logic); invzt_str hardened
(its raw-mat2str fallthrough threw on structs, masking the very error
ids it was formatting -- now a class/size placeholder, behavior-additive).
Tests: 5 fast (pure policy incl. re-entrance/zeros, validation
identifiers + the formatter regression, sampler contract on pre-lattice
guard paths, direct invz_refine_crossing coverage) RED->GREEN + 2
INVZ_SLOW gates (odd-off Bc<->Tc round-trip 0.05 K + length-MATCHED
poisonous-seed strip proof (R1) + sigma_floor threading (R6); odd-on
T-cut incl. the hard-window no-crossing invzt:bracket contract with its
unmutated-window diagnostic (R2) -- the capability the projected finder
lacks, ODD-LOG T2). Hard windows exit before any slide (R2) with winS
captured after the collapse guard (T2); message diagnostics asserted
via explicit try/catch (verifyError returns <missing> on a throw, T1);
the no-valid grow path terminates at the solve floor (R3), fast-covered
along with the hard-window diagnostic by the B=0 all-invalid trick
(T1/T3, zero physics solves); re-entrance indicators accumulate across
attempts into nseen/out.ncross (R5). CORE 55/0/3.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Driver rewiring — two-regime sweep + smoke

**Files:**
- Modify: `invz_tensor/invzt_run_phase_diagram.m` (FULL replacement with the spec rev 4's Component-5 block; the local `invzt_crit_at` at the bottom is deleted — Task 2's module file takes over, same name, same code)
- Temporary (deleted before commit): `invz_tensor/smoke_invzt_run_phase_diagram.m`

**Interfaces:**
- Consumes: `invzt_critical(ion, T, lat, Brange, bopts)` (its `tol` fed from the `Btol` knob, tesla); `[tc, out] = invzt_critical_T(ion, B, lat, topts)` (its `tol`/`width`/`gridstep` fed from `Ttol`/`Twidth`/`Tgridstep`, kelvin; `window` or `Tc0` per the `Twin` knob); `invzt_crit_at`, `invzt_tc_pm_extrap` (throws `invzt:tcNoWindow`); `invzt_str` for assert messages.
- Produces: workspace `Bc` `[1, numel(Ts)]`, `TcB` `[1, numel(Bs)]`, `Tc0_proxy`, `Tc0_closed`, `phase_boundary` `[nfinite, 2]` columns `[T(K) B(T)]`; one figure with both branches.

- [ ] **Step 1: Worktree checkpoint, then replace the driver**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && git status --short
```
(Reminder: `invz_tensor/invzt_run_spectra.m` is the user's WIP — do not touch.)

The new file content is the spec rev 4's Component-5 fenced ```matlab block, verbatim — extract and `diff`, expect byte-identical (it is the governing byte-parity copy from now on). Key structural points to sanity-check: `Btol`/`Ttol`/`Twidth`/`Tgridstep` knobs with the units-in-caps comments; `Ts` preflight alongside `Brange`/`Bs`/`Twin` (Twin additionally required to reach above the 0.02 K solve floor) plus finite-positive-real-scalar preflights for `Btol`/`Ttol`/`Twidth`/`Tgridstep` (R4 — `invzt_critical` does not validate its `opts.tol` at all), all `invzt_str`-formatted; proxy computed before the sweep when `show_proxy_anchor || (~isempty(Bs) && isempty(Twin))` with the `invzt_run_phase_diagram:tcAnchor` assert; `bopts`/`topts` assembled at the call boundary; flat two-kind `parfor` with val-then-assign and `invzt:bracket`-only absorption; T-cut plot branch; merged `phase_boundary`; NO local function at the bottom.

- [ ] **Step 2: Generate the reduced-size smoke copy (same directory)**

```bash
sed -e "s/^Ts     = linspace(0.4, 1.4, 11);/Ts     = [1.2 1.4];/" \
    -e "s/^Bs     = linspace(0.25, 1.5, 6);/Bs     = [0.5 1.5];/" \
    -e "s/^gridN  = 16;/gridN  = 8;/" \
    -e "s/^dpRng  = 30;/dpRng  = 15;/" \
    -e "s/^useParallel = true;/useParallel = false;/" \
    -e "s/^show_projected_anchor = false;/show_projected_anchor = true;/" \
    invz_tensor/invzt_run_phase_diagram.m > invz_tensor/smoke_invzt_run_phase_diagram.m
grep -n "^Ts \|^Bs \|^gridN\|^dpRng\|^useParallel\|^show_projected_anchor" invz_tensor/smoke_invzt_run_phase_diagram.m
```
Expected: all six overridden lines shown. (`Ts_proxy` stays 1.40:1/30:2.00 — the 8³ proxy Tc0 ≈ 1.53 K sits inside it, established by the drivers-plan smoke.)

- [ ] **Step 3: Run the smoke with assertions**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "run('invz_tensor/smoke_invzt_run_phase_diagram.m'); fprintf('Bc = [%s] | TcB = [%s] | proxy %.4f | closed %.4f\n', num2str(Bc, '%.3f '), num2str(TcB, '%.4f '), Tc0_proxy, Tc0_closed); assert(sum(isfinite(Bc)) == 2, 'smoke: both field cuts must bracket at 8^3 (known ~2.15/1.92 T)'); assert(all(isfinite(TcB)), 'smoke: BOTH T-cuts must resolve -- the odd-on two-branch demonstration'); assert(TcB(1) > TcB(2), 'smoke: Tc must decrease with field'); assert(TcB(1) <= Tc0_proxy + 0.02, 'smoke: Tc(0.5 T) must not exceed the zero-field proxy (+refinement tol)'); assert(isfinite(Tc0_closed), 'smoke: projected comparator did not resolve'); assert(size(phase_boundary, 1) >= 3, 'smoke: merged boundary too small'); assert(~isempty(findall(0, 'Type', 'figure')), 'smoke: no figure'); disp('SMOKE OK')"
```
Expected: `Bc ≈ [2.148 1.916]` (warm-cache knowns), both `TcB` finite with `Tc(0.5 T)` ≈ 1.50–1.53 K > `Tc(1.5 T)` ≈ 1.42–1.50 K (both between 1.40 K and the proxy Tc0 ≈ 1.53 — the 8³ odd-on boundary has Bc(1.4 K) = 1.916 T > 1.5 T), finite proxy + comparator, `SMOKE OK`. Estimated ~10–25 min serial, warm caches (estimate, not a gate). The T-cut success with `odd = true` (the driver default) is the two-branch demonstration; if either `TcB` is NaN, STOP and report with the console output — do not widen windows or toggle `odd`.

- [ ] **Step 4: Delete the smoke copy and verify byte-parity**

```bash
rm invz_tensor/smoke_invzt_run_phase_diagram.m
git status --short invz_tensor/
python3 - <<'EOF'
spec = open("docs/superpowers/specs/2026-07-18-invzt-critical-t-design.md", encoding="utf-8").read()
marker = "### 5. `invz_tensor/invzt_run_phase_diagram.m` — two-regime rewiring"
i = spec.index(marker); i = spec.index("```matlab\n", i) + len("```matlab\n")
j = spec.index("\n```", i)
block = spec[i:j]
file = open("invz_tensor/invzt_run_phase_diagram.m", encoding="utf-8").read().rstrip("\n")
print("byte-parity:", "OK" if block == file else "MISMATCH")
EOF
```
Expected: ` M invz_tensor/invzt_run_phase_diagram.m` (plus the unrelated user-WIP paths noted at the checkpoint); `byte-parity: OK`.

- [ ] **Step 5: Commit (stage only the driver)**

```bash
git add invz_tensor/invzt_run_phase_diagram.m
git commit -m "feat(invzt): two-regime phase-diagram driver -- Ts field cuts + Bs temperature cuts

Mirrors the projected invz_run_phase_diagram interface: fixed-T field
cuts (invzt_critical) plus fixed-B temperature cuts (invzt_critical_T)
in one flat two-kind parfor, merged T-sorted phase_boundary. Per-finder
option namespaces (review F5): Btol in TESLA and Ttol/Twidth/Tgridstep
in KELVIN, merged into bopts/topts at the call boundary -- both finders
name their option 'tol' in different units, so it is never shared. The
small-Bx proxy Tc0 moves BEFORE the sweep (it anchors the adaptive
T-cut window; invzt_run_phase_diagram:tcAnchor if the T-cuts need it
and it fails; explicit HARD Twin bypasses it). Ts preflighted alongside
Brange/Bs/Twin (Twin must reach above the 0.02 K solve floor), and the
four finder-control knobs Btol/Ttol/Twidth/Tgridstep preflighted as
finite positive real scalars before any expensive work (R4 --
invzt_critical does not validate its opts.tol), all invzt_str-formatted. Ts default drops its old
>1.5 K tail (guaranteed no-bracket NaNs -- the T-cut branch covers that
region properly now). Local invzt_crit_at deleted (module file from the
previous commit takes over). Smoke-verified at 8^3/dpRng-15 with odd
on: both T-cuts bracket (the projected T-cut cannot with ODD on),
Tc(0.5 T) > Tc(1.5 T), both below the proxy Tc0.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Docs + final gates

**Files:**
- Modify: `invz_tensor/README.html` (quickstart count, §2 callout — parenthetical AND scope sentence, §2.1 heading/intro/example/zero-field note, §7 slow-gate sentence, architecture table A1 row, module map ×3 rows)
- Modify: `docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md` (SUPERSEDED pointer on its Component-1 block)

**Interfaces:**
- Consumes: names/behavior from Tasks 2–3 verbatim.
- Produces: documentation only.

- [ ] **Step 1: Worktree checkpoint, then README quickstart count**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && git status --short
```

In `invz_tensor/README.html`, replace:

```html
<p>&rarr; expected <strong>50 passed / 0 failed / 1 incomplete</strong> (the incomplete is the <code>INVZ_SLOW</code>-gated A4 production-ladder test; the count includes the three 2026-07-18 regression tests added with the drivers &mdash; complex explicit-q response, A1-scope guard, odd-mask parity).</p>
```

with:

```html
<p>&rarr; expected <strong>55 passed / 0 failed / 3 incomplete</strong> (the incompletes are the <code>INVZ_SLOW</code>-gated A4 production-ladder test, the T-cut crossing-consistency round-trip, and the odd-on T-cut gate; the count includes the 2026-07-18 regression tests added with the drivers and the T-cut finder &mdash; complex explicit-q response, A1-scope guard, odd-mask parity, the <code>invzt_tc_pick</code> policy/zero suites, validation identifiers, sampler contract, and the <code>invz_refine_crossing</code> helper test).</p>
```

- [ ] **Step 2: §2 "Drivers" callout — BOTH stale phrases**

In the callout under `<h2 id="recipes">`, replace the parenthetical:

```html
<code>invzt_run_phase_diagram.m</code> (PM-side field-cut boundary B<sub>c</sub>(T) with zero-field anchors)
```

with:

```html
<code>invzt_run_phase_diagram.m</code> (PM-side boundary via fixed-T field cuts <em>and</em> fixed-B temperature cuts, with zero-field anchors)
```

and replace the scope sentence:

```html
Both are honestly scoped to PM-side physics: there is no ordered/FM tensor solve and no temperature-cut finder (<code>invzt_critical_T</code>), so the near-vertical region of the boundary approaching T<sub>c</sub>(0) is left blank by design (&sect;10 "Open items").
```

with:

```html
Both are honestly scoped to PM-side physics: there is no ordered/FM tensor solve. The phase-diagram driver searches <strong>both ways</strong> &mdash; fixed-T field cuts (<code>invzt_critical</code>, the low-T branch) and fixed-B temperature cuts (<code>invzt_critical_T</code>, added 2026-07-18: the near-vertical branch approaching T<sub>c</sub>(0), where a field cut crosses the boundary at a glancing angle; the tensor T-cut also brackets with <code>odd</code> on, which the projected T-cut cannot).
```

- [ ] **Step 3: §2.1 heading, intro, T-cut example, zero-field note**

Replace the heading:

```html
<h3>2.1 &nbsp;Phase diagram &mdash; Bc(T) via a loop over <code>invzt_critical</code></h3>
```

with:

```html
<h3>2.1 &nbsp;Phase diagram &mdash; Bc(T) field cuts and Tc(B) temperature cuts</h3>
```

Replace the intro sentence:

```html
<p><code>invzt_critical</code> (&sect;5) finds one Bc at one fixed T by bisection. <code>invzt_run_phase_diagram.m</code> packages the loop below (plus parfor, the small-B<sub>x</sub> proxy anchor, and an opt-in projected closed-form comparator); the inline version, for custom sweeps:</p>
```

with:

```html
<p><code>invzt_critical</code> (&sect;5) finds one Bc at one fixed T; <code>invzt_critical_T</code> finds one Tc at one fixed B (valid-samples-only T-grid + interpolation; explicit <code>opts.window</code> is a HARD bound, adaptive mode needs <code>opts.Tc0</code>). <code>invzt_run_phase_diagram.m</code> packages both as a two-regime sweep (<code>Ts</code> + <code>Bs</code> knobs, one flat parfor, per-finder tolerances <code>Btol</code>/<code>Ttol</code>, the small-B<sub>x</sub> proxy anchor, an opt-in projected comparator); the inline field-cut version, for custom sweeps &mdash; the one-line T-cut equivalent is e.g. <code>tc = invzt_critical_T(ion, 1.0, lat, struct('window', [1.0 1.8]))</code>:</p>
```

Replace the zero-field note (the reviewer's catch — the proxy is NOT a true-B=0 route):

```html
For the true zero-field endpoint <code>Tc(B=0)</code> specifically (not a small-<code>Bx</code> proxy), use <code>invzt_tc_pm_extrap</code> (linear crit&rarr;0 extrapolation from the two lowest converged-PM points on a fixed T grid) or the projected closed form <code>invz_odd_zero_field</code> &mdash; A3's true-zero-field solve is deferred (&sect;10 "Open items").
```

with:

```html
For the zero-field endpoint <code>Tc(B=0)</code>: the tensor-native route is the small-<code>Bx</code> PROXY <code>invzt_tc_pm_extrap</code> (linear crit&rarr;0 extrapolation from the two lowest converged-PM points at <code>Bproxy</code> &mdash; an estimate near B=0, never a true B=0 solve); the only true-zero-field result in the tensor scope is the projected closed form <code>invz_odd_zero_field</code> &mdash; the tensor true-B=0 solve stays deferred (&sect;10 "Open items").
```

(If any old-string differs slightly from the live file, locate it by grep and keep the NEW strings exactly as given.)

- [ ] **Step 4: §7 slow-gate sentence, architecture row, module map**

In §7 (suite topology), replace:

```html
this covers the A4 production ladder run (see "Open items").
```

with:

```html
this covers the A4 production ladder run (see "Open items"), the T-cut crossing-consistency round-trip, and the odd-on T-cut gate (test_invzt_critical_T).
```

In the §3 architecture table's A1 row, replace `<code>invzt_critical</code>, <code>invzt_tc_pm_extrap</code>` with `<code>invzt_critical</code>, <code>invzt_critical_T</code>, <code>invzt_tc_pm_extrap</code>`.

In the §4 module map, directly after the row

```html
<tr><td><code>invzt_critical.m</code></td><td>T7</td><td>PM-side critical transverse field Bc via <code>sign(crit)</code> bisection</td></tr>
```

insert:

```html
<tr><td><code>invzt_critical_T.m</code></td><td>drivers (2026-07-18)</td><td>PM-side critical temperature Tc(B) at fixed field &mdash; valid-samples-only T-grid, highest-T ordered&rarr;para root via <code>invzt_tc_pick</code>, regula-falsi refinement; explicit window = HARD bound, adaptive window REQUIRES <code>opts.Tc0</code> (no tensor zero-field closed form); strips <code>Sigma_seed</code> (Matsubara length is T-dependent)</td></tr>
<tr><td><code>invzt_tc_pick.m</code></td><td>drivers (2026-07-18)</td><td>PURE crossing/slide decision core on valid votes (zero-recognition, re-entrance counting, top-voter-ordered &rarr; move up) &mdash; closes the projected classifier's inherited exact-zero and re-entrant-lower-leg gaps; synthetically tested</td></tr>
<tr><td><code>invzt_crit_at.m</code></td><td>drivers (2026-07-18)</td><td>Shared one-sample criticality probe: crit + three-part validity (validity-only <code>ok</code> &mdash; consumers apply their own phase logic); selective physics-signal catch, all else rethrows</td></tr>
```

- [ ] **Step 5: SUPERSEDED pointer in the drivers spec**

In `docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md`, directly after the heading line ``### 1. `invz_tensor/invzt_run_phase_diagram.m` `` insert:

```markdown
**SUPERSEDED (2026-07-18, T-cut extension):** the governing byte-parity copy
of this driver now lives in
`docs/superpowers/specs/2026-07-18-invzt-critical-t-design.md` Component 5
(two-regime `Ts` + `Bs` sweep with per-finder tolerances). The block below is
the as-first-shipped field-cut-only version, kept for the execution record.
```

- [ ] **Step 6: README structure check (stale-phrase assertions) + final gates**

```bash
python3 -c "
import re
s = open('invz_tensor/README.html').read()
for stale in ['no temperature-cut finder',
              'field-cut boundary B<sub>c</sub>(T) with zero-field anchors',
              'For the true zero-field endpoint',
              '50 passed / 0 failed / 1 incomplete']:
    assert stale not in s, 'stale phrase remains: ' + stale
for fresh in ['invzt_critical_T', 'invzt_tc_pick', 'invzt_crit_at',
              '55 passed / 0 failed / 3 incomplete']:
    assert fresh in s, 'expected phrase missing: ' + fresh
for tag in ['div','table','ul','ol','pre']:
    o = len(re.findall(r'<'+tag+r'[ >]', s)); c = len(re.findall(r'</'+tag+'>', s))
    assert o == c, f'{tag}: {o} open vs {c} close'
print('README checks OK')
"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
```
Expected: `README checks OK`; **55 / 0 / 3**; **143 / 0 / 19**.

- [ ] **Step 7: Commit (stage exactly these two files)**

```bash
git add invz_tensor/README.html "docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md"
git diff --cached --name-status
git commit -m "docs(invzt): README reflects the T-cut finder + two-regime driver

Quickstart count 50->55/0/3 with the incompletes named; section-2
callout (parenthetical AND scope sentence) and 2.1 heading/intro/
zero-field note rewritten for two-regime search -- incl. the review
catch that invzt_tc_pm_extrap was mislabeled a true-B=0 route (it is
the small-Bx proxy; only the projected closed form is true zero-field
in tensor scope); slow-gate sentence names all three INVZ_SLOW tests;
architecture/module-map rows for invzt_critical_T + invzt_tc_pick +
invzt_crit_at; drivers spec Component-1 block marked SUPERSEDED by the
T-cut spec's Component 5 (the governing byte-parity copy).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
The `git diff --cached --name-status` must list exactly: `M	invz_tensor/README.html`, `M	docs/superpowers/specs/2026-07-18-invzt-run-drivers-design.md`.

---

## Execution notes

- Task order: 1 → 2 → 3 → 4, strictly (2 needs 1's helper on the path; 3 needs 2's module files; 4 documents 2–3).
- The 8³/dpRng-15 caches are warm from the drivers-plan smokes; all durations are estimates, never gates — on an apparent hang, capture and report the console tail instead of killing early.
- Production sweeps at the committed knobs (16³, 11 field cuts + 6 T-cuts) remain "left to the user", per repo precedent.
- During Task 2 the driver's local `invzt_crit_at` deliberately coexists with the new module file (locals shadow path functions inside the script; the code is identical) — flag it in neither task as a defect; Task 3 removes the local.
- The projected `invz_critical_T` carries the same F1 classifier gaps this work fixes on the tensor side — deliberately untouched (validated code, own test surface); noted in the spec as a candidate follow-up for the user.
