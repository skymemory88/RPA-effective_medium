# c-axis Field-Misalignment (Tilt) Knob Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a c-axis field-tilt knob (`theta_c`) to `invz_run_spectra` by threading a full field vector `[Bx By Bz]` through the 1/z solve chain, routing longitudinal fields through the moment-form self-energy with sign-aware branch selection.

**Architecture:** A tiny normalizer (`invz_field_vec`) maps scalar fields to `[B 0 0]` so every existing caller is bit-for-bit unchanged; `invz_solve_auto` gains a `bz_tol` dead band that routes `|Bz| > bz_tol` exclusively through `invz_solve_point_ordered` with a `forced_moment` flag (the induced moment is physical; `invz_sigma_ordered` reduces to `invz_sigma` as `m→0`). The spectra drivers gain a `field_dir`/`solve_opts` API, and a Σ=0 scalar-vs-3×3-tensor reference quantifies the omitted cross-channel error to set the supported tilt range.

**Tech Stack:** MATLAB R2025a, `matlab.unittest` function-based tests (`functiontests(localfunctions)`).

**Spec:** `docs/superpowers/specs/2026-07-16-invz-field-angle-design.md` (read it once before starting; it has the physics background and the review-resolution history).

## Global Constraints

- Repo path contains spaces — ALWAYS quote: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "..."`. Run all MATLAB commands from the repo root (`.../Simulation/invZ expansion`).
- Fast suite must pass after EVERY task: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"` (baseline: all pass, a few skip as "incomplete" without `INVZ_SLOW`/removed-src — that is normal; `assertSuccess` tolerates skips, not failures).
- Slow suite after Tasks 8, 9, 10: same command prefixed `INVZ_SLOW=1 `.
- **Non-negotiable:** at `theta_c = 0` / scalar field arguments, every solver output is bit-for-bit unchanged and every published benchmark test passes untouched. NEVER edit published benchmark values in existing tests.
- Do NOT modify existing test files; all new tests go in the three new test files named below.
- New error identifiers: `invz:fieldVec`, `invz:solveOpts`, `invz:fieldDir`, `invz:fields`. All thrown errors keep `invz:*` prefixes (the `invz_solve_auto` absorb-policy depends on it).
- `bz_tol` default is `1e-9` (Tesla) everywhere it appears.
- Commit style: `type(invz): summary` (see recent history: `feat(invz):`, `test(invz):`, `docs(invz):`).

## File Structure

| File | Status | Responsibility |
|---|---|---|
| `invz/invz_field_vec.m` | create | scalar/3-vector field normalization (single validation point) |
| `invz/invz_single_ion.m` | modify | sign-aware `mz_seed`, `mf_converged/iters/residual`, `E0`, `F_mf` |
| `invz/invz_twolevel.m`, `invz/invz_twolevel_ordered.m` | modify | accept scalar-or-vector B (one line each) |
| `invz/invz_solve_point.m` | modify | normalize B once, pass vector down |
| `invz/invz_solve_point_ordered.m` | modify | `forced_moment`, sign-retry, complete early returns, `moment_branch` |
| `invz/invz_solve_auto.m` | modify | `bz_tol` dead band + longitudinal routing |
| `invz/invz_chi_realaxis.m` | modify | vector-capable paramagnet fallback (one line) |
| `invz/invz_spectra_map.m` | modify | `field_dir`/`bz_tol`/`solve_opts` API, dead band pre-parfor, longitudinal failure contract, metadata |
| `invz/invz_spectra_qpath.m` | modify | vector B, `S.Bvec`/`S.Bmag`, `mat2str` error, same opts contract |
| `invz/invz_run_spectra.m` | modify | `theta_c` knob, `dhat`, labels |
| `invz/invz_plot_spectra_map.m` | modify | axis label `|B| (T)` |
| `invz/invz_chi_tensor_ref.m` | create | Σ=0 scalar-vs-tensor cross-channel reference |
| `invz/invz_run_tensor_ref.m` | create | measurement driver: eps_spec/dE_peak vs angle table |
| `invz/tests/test_invz_field_vec.m` | create | Task 1 tests |
| `invz/tests/test_invz_field_angle.m` | create | Tasks 2-4 solver-level tests |
| `invz/tests/test_invz_field_angle_spectra.m` | create | Tasks 5-7 spectra/label tests + fast mirror test |
| `invz/tests/test_invz_field_angle_slow.m` | create | Task 8 slow tests (INVZ_SLOW-gated) |
| `invz/tests/test_invz_tensor_ref.m` | create | Task 9 reference tests |
| `docs/SESSION-2026-07-16-field-angle.md` | create | measured tensor-reference table + supported angle range |
| `invz/README.html` | modify | function reference + knob documentation |

---

### Task 1: `invz_field_vec` normalizer

**Files:**
- Create: `invz/invz_field_vec.m`
- Test: `invz/tests/test_invz_field_vec.m`

**Interfaces:**
- Produces: `B = invz_field_vec(B)` — scalar `b → [b 0 0]`; any real finite numeric 3-element vector (row or column) → `1x3` row; everything else errors with identifier `invz:fieldVec`. Idempotent on valid input. Every later task calls this at its public boundary.

- [ ] **Step 1: Write the failing test**

Create `invz/tests/test_invz_field_vec.m`:

```matlab
function tests = test_invz_field_vec
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_scalar_maps_to_transverse_row(testCase)
verifyEqual(testCase, invz_field_vec(4), [4 0 0]);
verifyEqual(testCase, invz_field_vec(0), [0 0 0]);
verifyEqual(testCase, invz_field_vec(-2.5), [-2.5 0 0]);   % signed amplitude passes through
end

function test_vectors_normalize_to_row(testCase)
verifyEqual(testCase, invz_field_vec([1 2 3]), [1 2 3]);
verifyEqual(testCase, invz_field_vec([1; 2; 3]), [1 2 3]); % column -> row
verifyEqual(testCase, invz_field_vec(invz_field_vec(5)), [5 0 0]);  % idempotent
end

function test_invalid_inputs_error(testCase)
bad = {NaN, [1 Inf 0], 1+2i, [], [1 2], [1 2 3 4], 'abc', true};
for k = 1:numel(bad)
    verifyError(testCase, @() invz_field_vec(bad{k}), 'invz:fieldVec');
end
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_vec.m'); assertSuccess(results)"
```
Expected: FAIL — `invz_field_vec` undefined.

- [ ] **Step 3: Write the implementation**

Create `invz/invz_field_vec.m`:

```matlab
function B = invz_field_vec(B)
%INVZ_FIELD_VEC Normalize a field argument to a 1x3 row [Bx By Bz] in Tesla.
% Scalar b -> [b 0 0] (the historical transverse-along-a convention); any real,
% finite, numeric 3-element vector (row or column) -> 1x3 row, values unchanged.
% Anything else (NaN/Inf, complex, empty, wrong length, non-numeric) errors with
% identifier 'invz:fieldVec'. Idempotent, so boundaries may normalize freely.
if ~isnumeric(B) || ~isreal(B) || ~all(isfinite(B(:)))
    error('invz:fieldVec', 'Field must be real, finite, and numeric.');
end
if isscalar(B)
    B = [B 0 0];
elseif numel(B) == 3
    B = reshape(B, 1, 3);
else
    error('invz:fieldVec', ...
        'Field must be a scalar or a 3-element vector; got %d element(s).', numel(B));
end
end
```

(Note: `all(isfinite([])) == true`, so the empty case falls through to the `numel` error — covered by the test.)

- [ ] **Step 4: Run test to verify it passes**

Same command as Step 2. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_field_vec.m invz/tests/test_invz_field_vec.m
git commit -m "feat(invz): invz_field_vec scalar/3-vector field normalizer"
```

---

### Task 2: `invz_single_ion` — sign-aware seed, MF convergence state, `E0`/`F_mf`

**Files:**
- Modify: `invz/invz_single_ion.m` (line 19 seed default; new output fields after line 83)
- Test: `invz/tests/test_invz_field_angle.m` (new file)

**Interfaces:**
- Consumes: nothing new.
- Produces (additive `si` fields, no existing caller changes):
  `si.mf_converged` (logical), `si.mf_iters` (int), `si.mf_residual` (double, final `dmf`),
  `si.E0` (double, unshifted ground energy — `si.E` stays shifted),
  `si.F_mf` (double, variational MF free energy; `NaN` in `hz_fixed` mode).
  New default: in `order` mode with `B(3) < 0`, `mz_seed` defaults to `-1` (explicit `opts.mz_seed` still wins).

- [ ] **Step 1: Write the failing tests**

Create `invz/tests/test_invz_field_angle.m`:

```matlab
function tests = test_invz_field_angle
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_transverse_baseline_anchor(testCase)
% Frozen x-field baseline (verified 2026-07-16, matches the Codex review probe):
% guards the seed/diagnostic edits against any accidental transverse change.
si = invz_single_ion(invz_ion(), 0.31, [4 0 0], struct('hyp', false));
verifyEqual(testCase, si.Jexp(1),        3.58424135285,  'AbsTol', 1e-9);
verifyEqual(testCase, si.Jexp(2),       -0.0689991117463,'AbsTol', 1e-10);
verifyEqual(testCase, si.E(2) - si.E(1), 0.369235620278, 'AbsTol', 1e-10);
verifyTrue(testCase, si.mf_converged);
verifyLessThan(testCase, si.mf_residual, 1e-12);
end

function test_sign_aware_branch_and_Fmf(testCase)
% Spec test 3 (second review finding 1, verified digit-for-digit): an explicit
% Bz selects the aligned branch; the VARIATIONAL F_mf ranks branches correctly
% (the naive shifted-spectrum comparison mis-ranks them and must not be used).
ion = invz_ion();  T = 0.31;
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
sp = invz_single_ion(ion, T, [2 0 +0.01], struct('hyp', false, 'order', true));
sm = invz_single_ion(ion, T, [2 0 -0.01], struct('hyp', false, 'order', true));
verifyGreaterThan(testCase, sp.Jexp(3), 0);
verifyLessThan(testCase, sm.Jexp(3), 0);
verifyLessThan(testCase, abs(sp.Jexp(3) + sm.Jexp(3)), 1e-10);      % exact Z2 mirror
% force the metastable anti-aligned branch with an explicit wrong-sign seed
sw = invz_single_ion(ion, T, [2 0 -0.01], struct('hyp', false, 'order', true, 'mz_seed', +1));
verifyGreaterThan(testCase, sw.Jexp(3), 0);                          % metastable branch reached
verifyLessThan(testCase, sm.F_mf, sw.F_mf);                          % aligned branch lower
verifyEqual(testCase, sm.F_mf, -21.4766457412, 'AbsTol', 1e-6);      % verified anchors
verifyEqual(testCase, sw.F_mf, -21.4696393612, 'AbsTol', 1e-6);
end

function test_mf_convergence_reporting(testCase)
ion = invz_ion();
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
si = invz_single_ion(ion, 0.31, [2 0 0.01], struct('hyp', false, 'order', true, 'mf_maxit', 1));
verifyFalse(testCase, si.mf_converged);
verifyEqual(testCase, si.mf_iters, 1);
verifyGreaterThan(testCase, si.mf_residual, 1e-12);
% hz_fixed mode: F_mf undefined by design (hz is imposed, not variational)
sf = invz_single_ion(ion, 0.31, [2 0 0], struct('hyp', false, 'hz_fixed', 5e-3));
verifyTrue(testCase, isnan(sf.F_mf));
verifyTrue(testCase, isfinite(sf.E0));
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle.m'); assertSuccess(results)"
```
Expected: FAIL — `mf_converged` field does not exist (and sign test fails on the `+` branch for `Bz < 0`).

- [ ] **Step 3: Implement**

In `invz/invz_single_ion.m`, replace line 19:

```matlab
mzsd = 1.0;         if isfield(opts,'mz_seed'), mzsd = opts.mz_seed; end
```

with:

```matlab
mzsd = 1.0;         if order && B(3) < 0, mzsd = -1.0; end
                    % sign-aware default: an explicit longitudinal field selects the
                    % ALIGNED branch (the +1 seed can trap the metastable mirror state
                    % below the transition; spec 2026-07-16, review finding 2)
if isfield(opts,'mz_seed'), mzsd = opts.mz_seed; end
```

After line 83 (`si.JzJz_fluct = jz2 - si.Jexp(3)^2;`), insert before the final `end`:

```matlab
si.mf_converged = converged;
si.mf_iters     = it;
si.mf_residual  = dmf;
si.E0 = E(1);                                  % unshifted ground energy (si.E stays shifted)
% Variational MF free energy (branch diagnostic, spec SR1): the 0.5*h*<J> terms
% restore the -1/2 J <J>^2 double counting; equals Fsite + hx^2/(2Jxx0) + hz^2/(2J0z)
% at a self-consistent point. Undefined (NaN) under hz_fixed: hz is imposed, not a MF.
Fsite = E(1) - log(sum(exp(-beta*(E - E(1)))))/beta;
if hzfix
    si.F_mf = NaN;
else
    si.F_mf = Fsite + 0.5*(hx*si.Jexp(1) + hz*si.Jexp(3));
end
```

Also add to the header comment (after the mode descriptions):

```matlab
% Returned diagnostics: si.mf_converged/mf_iters/mf_residual (mean-field loop state),
% si.E0 (unshifted ground energy), si.F_mf (variational MF free energy; NaN with hz_fixed).
% In order mode the mz_seed default is sign-aware: -1 when B(3) < 0 (aligned branch).
```

- [ ] **Step 4: Run to verify pass**

Same command as Step 2. Expected: PASS (3 tests).

- [ ] **Step 5: Run the full fast suite (regression gate)**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```
Expected: PASS. The seed edit is inert for every existing caller (all pass `B(3) = 0`).

- [ ] **Step 6: Commit**

```bash
git add invz/invz_single_ion.m invz/tests/test_invz_field_angle.m
git commit -m "feat(invz): sign-aware mz_seed + MF convergence/E0/F_mf diagnostics in invz_single_ion"
```

---

### Task 3: Field-vector plumbing at the five `[Bx 0 0]` leaf sites

**Files:**
- Modify: `invz/invz_twolevel.m:7`, `invz/invz_twolevel_ordered.m:13`, `invz/invz_solve_point.m:16`, `invz/invz_solve_point_ordered.m:34`, `invz/invz_chi_realaxis.m:37`
- Test: `invz/tests/test_invz_field_angle.m` (append)

**Interfaces:**
- Consumes: `invz_field_vec` (Task 1).
- Produces: every one of these functions now accepts scalar OR 3-vector `Bx`, normalizing once at its own boundary (idempotent, so nesting is safe). No signature changes, no output changes for scalar input.

- [ ] **Step 1: Write the failing test** (append to `test_invz_field_angle.m`)

```matlab
function test_scalar_vs_vector_boundaries(testCase)
% Spec test 2: scalar B and [B 0 0] are literally the same solve at every
% scalar-accepting boundary (struct-exact equality: identical code path).
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';       % synthetic branch fixture (fast, no lattice sum)
o   = struct('J0eff', 6.4e-3);
verifyEqual(testCase, invz_twolevel(ion, 0.31, [5.5 0 0]), invz_twolevel(ion, 0.31, 5.5));
verifyEqual(testCase, invz_twolevel_ordered(ion, 0.31, [2 0 0], 5e-3), ...
                      invz_twolevel_ordered(ion, 0.31, 2, 5e-3));
pt1 = invz_solve_point(ion, 0.31, 5.5, Jnu, o);
pt2 = invz_solve_point(ion, 0.31, [5.5 0 0], Jnu, o);
verifyEqual(testCase, pt2, pt1);
[pa1, ph1] = invz_solve_auto(ion, 0.31, 5.5, Jnu, o);
[pa2, ph2] = invz_solve_auto(ion, 0.31, [5.5; 0; 0], Jnu, o);   % column form too
verifyEqual(testCase, ph2, ph1);
verifyEqual(testCase, pa2, pa1);
w = (0.1:0.1:0.4).';
o1 = invz_chi_realaxis(ion, 0.31, 5.5, pt1, w, struct('eta', 1e-3));
o2 = invz_chi_realaxis(ion, 0.31, [5.5 0 0], pt1, w, struct('eta', 1e-3));
verifyEqual(testCase, o2.chi_cc_q, o1.chi_cc_q);
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle.m'); assertSuccess(results)"
```
Expected: FAIL — `invz_twolevel` builds `[Bx 0 0]` from a vector `Bx` (dimension error or wrong Hamiltonian).

- [ ] **Step 3: Implement (five one-line normalizations)**

`invz/invz_twolevel.m` — replace line 7:
```matlab
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false, 'Jxx0', Jxx0));
```
with:
```matlab
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', false, 'Jxx0', Jxx0));
```

`invz/invz_twolevel_ordered.m` — replace line 13 the same way (keep the `hz_fixed` option):
```matlab
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', false, 'hz_fixed', hz, 'Jxx0', Jxx0));
```

`invz/invz_solve_point.m` — after the opts block (after line 13), insert:
```matlab
Bx = invz_field_vec(Bx);                       % scalar -> [Bx 0 0]; 3-vector passes through
```
and replace line 16's `[Bx 0 0]` with `Bx`:
```matlab
si  = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0));
```
(line 19 `invz_twolevel(ion, T, Bx, ...)` is already vector-safe after the edit above).

`invz/invz_solve_point_ordered.m` — after the opts block (after line 26), insert:
```matlab
Bx = invz_field_vec(Bx);                       % scalar -> [Bx 0 0]; 3-vector passes through
```
and replace line 34's `[Bx 0 0]` with `Bx`:
```matlab
si = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'order', true, 'J0z', J0eff, 'Jxx0', Jxx0));
```

`invz/invz_chi_realaxis.m` — replace line 37:
```matlab
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', hyp, 'Jxx0', Jxx0));   % paramagnet
```
with:
```matlab
si = invz_single_ion(ion, T, invz_field_vec(Bx), struct('hyp', hyp, 'Jxx0', Jxx0));   % paramagnet
```

Update each touched function's header comment: change "at one (T, Bx) point" style text to note "`Bx`: scalar (transverse, historical) or `[Bx By Bz]` vector (T)".

- [ ] **Step 4: Run to verify pass, then full fast suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle.m'); assertSuccess(results)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```
Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_twolevel.m invz/invz_twolevel_ordered.m invz/invz_solve_point.m invz/invz_solve_point_ordered.m invz/invz_chi_realaxis.m invz/tests/test_invz_field_angle.m
git commit -m "feat(invz): vector-capable field at the five [Bx 0 0] leaf sites"
```

---

### Task 4: `forced_moment` + longitudinal routing (`invz_solve_point_ordered`, `invz_solve_auto`)

**Files:**
- Modify: `invz/invz_solve_point_ordered.m`, `invz/invz_solve_auto.m`
- Test: `invz/tests/test_invz_field_angle.m` (append)

**Interfaces:**
- Consumes: Task 2 (`si.mf_converged`, sign-aware seed), Task 3 (vector plumbing).
- Produces:
  - `invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts)` gains `opts.forced_moment` (logical, default false) and forwards `opts.mz_seed`/`opts.mf_maxit`/`opts.mf_mix` to `invz_single_ion`. EVERY return path now carries the full field set: `m0, is_ordered, converged, Sigma0, crit, si, tl, moment_branch` (`tl = []` on early returns; `moment_branch ∈ {'spontaneous','field_induced','none'}` — `'none'` only on the spontaneous-mode paramagnetic early return).
  - `invz_solve_auto(ion, T, Bx, Jnu_flat, opts)` gains `opts.bz_tol` (default `1e-9` T): `|Bz| <= bz_tol` is zeroed and routed through today's transverse logic verbatim; `|Bz| > bz_tol` routes ONLY to the moment solver with `forced_moment = true`, and on failure returns the failed `pto` (not `[]`) whenever it carries a nonempty `si`, so the map can build an RPA-only overlay.

- [ ] **Step 1: Write the failing tests** (append to `test_invz_field_angle.m`)

```matlab
function test_longitudinal_routing_threshold(testCase)
% Spec test 7: just below bz_tol -> transverse path (strict-PM two-level);
% just above -> forced moment-form with machine-readable branch metadata.
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o   = struct('J0eff', 6.4e-3);
[ptA, phA] = invz_solve_auto(ion, 0.31, [5.5 0 0.5e-9], Jnu, o);   % |Bz| <= 1e-9: dead band
verifyEqual(testCase, phA, 2);                                      % strict paramagnet
verifyEqual(testCase, ptA.tl.m, 0, 'AbsTol', 1e-3);                 % invz_twolevel gate honored
[ptB, phB] = invz_solve_auto(ion, 0.31, [5.5 0 2e-9], Jnu, o);      % above: moment route
verifyEqual(testCase, phB, 1);
verifyEqual(testCase, ptB.moment_branch, 'field_induced');
verifyTrue(testCase, ptB.is_ordered);                               % moment-form self-energy flag
end

function test_early_return_struct_completeness(testCase)
% Spec test 12: every early return of invz_solve_point_ordered carries the full
% declared field set, so invz_solve_auto / the map never probe a missing member.
ion  = invz_ion();
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
flds = {'m0','is_ordered','converged','Sigma0','crit','si','tl','moment_branch'};
% (a) spontaneous-mode paramagnetic early return (PM point: no spontaneous moment)
pta = invz_solve_point_ordered(ion, 1.0, 5.5, Jnu, struct('J0eff', 6.4e-3));
verifyFalse(testCase, pta.is_ordered);
verifyEqual(testCase, pta.moment_branch, 'none');
cellfun(@(f) verifyTrue(testCase, isfield(pta, f), ['missing ' f]), flds);
% (b) forced_moment with a crippled MF loop -> mf-gate early return
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
ptb = invz_solve_point_ordered(ion, 0.31, [2 0 0.01], Jnu, ...
      struct('J0eff', 6.4e-3, 'forced_moment', true, 'mf_maxit', 1));
verifyFalse(testCase, ptb.converged);
verifyEqual(testCase, ptb.moment_branch, 'field_induced');
cellfun(@(f) verifyTrue(testCase, isfield(ptb, f), ['missing ' f]), flds);
end

function test_solve_auto_returns_failed_pto(testCase)
% Second-review finding 2: a failed longitudinal solve returns the pto (with si)
% so invz_spectra_map can compute its RPA-only overlay -- never pt = [].
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
ws = warning('off', 'invz:mfNotConverged');
restore = onCleanup(@() warning(ws));
[ptc, phc] = invz_solve_auto(ion, 0.31, [2 0 0.01], Jnu, ...
             struct('J0eff', 6.4e-3, 'mf_maxit', 1));
verifyEqual(testCase, phc, 0);
verifyFalse(testCase, isempty(ptc));
verifyFalse(testCase, isempty(ptc.si));
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle.m'); assertSuccess(results)"
```
Expected: FAIL — no `moment_branch` field / longitudinal `[5.5 0 2e-9]` reaches the strict-PM `invz_twolevel` gate.

- [ ] **Step 3: Implement `invz_solve_point_ordered.m`**

After the opts block, extend option handling (keep the Task-3 `Bx = invz_field_vec(Bx);` line above this):

```matlab
fmom = isfield(opts,'forced_moment') && opts.forced_moment;
```

Replace the single-ion call block (lines 30-41 of the current file: the `si = ...` through the early return) with:

```matlab
% Ordered mean-field solve (full electronuclear space): spontaneous moment m0 and field hz.
% J0z is the SAME cc coupling J(0) used by the criticality and the RPA/1z denominator.
% forced_moment (spec 2026-07-16): with an explicit longitudinal Bx(3) the moment is
% field-induced -- the spontaneous |m0| > mtol gate is bypassed and branch alignment
% with the applied Bz is enforced (sign-aware seed + one mirrored retry).
siopts = struct('hyp', hyp, 'order', true, 'J0z', J0eff, 'Jxx0', Jxx0);
for f = {'mz_seed', 'mf_maxit', 'mf_mix'}                  % diagnostic pass-throughs (tests)
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
end
si = invz_single_ion(ion, T, Bx, siopts);
branch = 'spontaneous';  if fmom, branch = 'field_induced'; end
if fmom && si.mf_converged && abs(si.Jexp(3)) > 1e-10 && sign(si.Jexp(3)) ~= sign(Bx(3))
    % converged onto the metastable anti-aligned branch: one mirrored retry
    siopts.mz_seed = -sign(si.Jexp(3));
    si2 = invz_single_ion(ion, T, Bx, siopts);
    if si2.mf_converged && sign(si2.Jexp(3)) == sign(Bx(3))
        si = si2;
    else
        warning('invz:branchMismatch', ...
            'Anti-aligned moment persists at Bz = %.3g T after mirrored retry.', Bx(3));
        pt = early_return(si.Jexp(3), si, branch);
        return;
    end
end
if fmom && ~si.mf_converged
    pt = early_return(si.Jexp(3), si, branch);             % MF gate (second review finding 6)
    return;
end
m0 = si.Jexp(3);
pt.m0 = m0;
pt.is_ordered = fmom || abs(m0) > mtol;
if ~pt.is_ordered
    pt = early_return(m0, si, 'none');                     % paramagnetic point: use invz_solve_point
    return;
end
```

At the end of the function, after `pt.outer_iters = outer;`, add:

```matlab
pt.moment_branch = branch;
```

And add the local function at the bottom of the file (after the main function's `end`):

```matlab
% -------------------------------------------------------------------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (spec: callers never
% probe a missing member; tl = [] flags "no two-level params were built").
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch);
end
```

Update the header: document `opts.forced_moment`, `pt.moment_branch`, and that `pt.is_ordered` means strictly "uses the moment-form self-energy" (spontaneous FM **or** field-induced).

- [ ] **Step 4: Implement `invz_solve_auto.m`** (full replacement — the file is short)

```matlab
function [pt, phase, di] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Ordered-first 1/z solve at one (T, B) point, paramagnetic fallback.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
%
% Transverse route (|Bz| <= opts.bz_tol, default 1e-9 T; the component is ZEROED in
% the dead band): tries invz_solve_point_ordered (spontaneous moment, phase = 1),
% then invz_solve_point (paramagnet, phase = 2) -- identical to the historical logic.
%
% Longitudinal route (|Bz| > bz_tol): the moment is field-induced at every (T, B), so
% the solve goes EXCLUSIVELY through invz_solve_point_ordered with forced_moment = true
% (the strict-paramagnet solver would reject m ~= 0; invz_sigma_ordered reduces to the
% paramagnet form as m -> 0). phase = 1 on acceptance; on failure phase = 0 and pt is
% the failed pto WHENEVER it carries a nonempty si (RPA-overlay fallback for the map),
% pt = [] only when no usable single-ion state exists.
%
% Error policy: only invz:* identifiers are absorbed into di; anything else rethrows.
if nargin < 5, opts = struct(); end
bz_tol = getf(opts, 'bz_tol', 1e-9);
B = invz_field_vec(Bx);
if abs(B(3)) <= bz_tol, B(3) = 0; end          % dead band: exactly transverse below tolerance
pt = [];  phase = 0;  di = struct('ordered_err', '', 'para_err', '');

if B(3) ~= 0
    oo = opts;  oo.forced_moment = true;
    try
        pto = invz_solve_point_ordered(ion, T, B, Jnu_flat, oo);
        if pto.is_ordered && pto.converged && isfinite(pto.Sigma0) && pto.si.mf_converged
            pt = pto;  phase = 1;
        elseif ~isempty(pto.si)
            pt = pto;                          % failed, but si (and maybe tl) support an overlay
        end
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        di.ordered_err = err.identifier;
    end
    return;
end

try
    pto = invz_solve_point_ordered(ion, T, B, Jnu_flat, opts);
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
        pt = pto;  phase = 1;  return;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    di.ordered_err = err.identifier;
end
try
    pt = invz_solve_point(ion, T, B, Jnu_flat, opts);
    if pt.converged && isfinite(pt.Sigma0)
        phase = 2;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    di.para_err = err.identifier;
    pt = [];
end
end
```

- [ ] **Step 5: Run to verify pass, then full fast suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle.m'); assertSuccess(results)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```
Expected: both PASS. (The transverse branch is byte-identical logic, so `test_invz_solve_auto`'s existing phase-0 and rethrow tests still pass; `getf` is already on the path for all callers since `invz/` is one flat directory.)

- [ ] **Step 6: Commit**

```bash
git add invz/invz_solve_point_ordered.m invz/invz_solve_auto.m invz/tests/test_invz_field_angle.m
git commit -m "feat(invz): forced_moment longitudinal route with bz_tol dead band and complete return contracts"
```

---

### Task 5: `invz_spectra_map` — `field_dir` API, option propagation, failure contract, metadata

**Files:**
- Modify: `invz/invz_spectra_map.m`
- Test: `invz/tests/test_invz_field_angle_spectra.m` (new file)

**Interfaces:**
- Consumes: Tasks 1-4 (`invz_field_vec`, `bz_tol` routing, `forced_moment`, failed-`pto` return).
- Produces: `invz_spectra_map(ion, T, fields, w, opts)` gains
  `opts.field_dir` (default `[1 0 0]`; nonzero finite real 3-vector, normalized internally; error `invz:fieldDir`),
  `opts.bz_tol` (default `1e-9`; resolved ONCE — used for the pre-parfor dead band, forwarded in `sopts`, and used by `one_field`),
  `opts.solve_opts` (struct, default empty; merged into `sopts`; fields `J0eff`/`Jxx0`/`hyp` reserved → error `invz:solveOpts`).
  `fields` must be nonnegative (error `invz:fields`). Returns `S.field_dir` (normalized 1x3) and `S.Bvec` (`nB x 3`, the vectors actually used, dead band applied). `S.phase` semantics: `1` = moment-form (spontaneous FM or field-induced), `2` = strict paramagnet, `0` = masked.

- [ ] **Step 1: Write the failing tests**

Create `invz/tests/test_invz_field_angle_spectra.m`:

```matlab
function tests = test_invz_field_angle_spectra
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
% Synthetic-coupling fixture (matches test_invz_spectra_map: fast, no lattice sum).
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', ...
            'info', struct('Jcc0', 6.4e-3), 'verbose', false);
end

function test_default_dir_equals_explicit_x(testCase)
% Spec test 2 (map form): default field_dir is exactly [1 0 0].
ion = invz_ion();  T = 0.31;  fields = [2.0 5.5];  w = (0.02:0.02:0.6).';
S1 = invz_spectra_map(ion, T, fields, w, fixture());
fx = fixture();  fx.field_dir = [1 0 0];
S2 = invz_spectra_map(ion, T, fields, w, fx);
verifyEqual(testCase, S2.chiz, S1.chiz);
verifyEqual(testCase, S2.chirpa, S1.chirpa);
verifyEqual(testCase, S1.field_dir, [1 0 0]);
verifyEqual(testCase, S1.Bvec, [fields(:) zeros(2, 2)]);
end

function test_api_validation(testCase)
ion = invz_ion();  w = (0.1:0.1:0.5).';
fx = fixture();  fx.solve_opts = struct('Jxx0', 1e-3);            % reserved field
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, fx), 'invz:solveOpts');
fx = fixture();  fx.field_dir = [0 0 0];
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, fx), 'invz:fieldDir');
verifyError(testCase, @() invz_spectra_map(ion, 0.31, -1.0, w, fixture()), 'invz:fields');
end

function test_deadband_zeroes_tiny_tilt(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [1 0 1e-12];       % Bz = 5.5e-12 T << bz_tol
S = invz_spectra_map(ion, 0.31, 5.5, w, fx);
verifyEqual(testCase, S.Bvec(3), 0);               % dead band applied BEFORE the solve
verifyEqual(testCase, S.phase, 2);                 % genuinely transverse: strict PM
end

function test_tilted_sweep_is_moment_form(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(1) 0 sind(1)];
S = invz_spectra_map(ion, 0.31, [2.0 5.5], w, fx);
verifyEqual(testCase, S.phase, [1 1]);             % field-induced moment at BOTH fields
verifyTrue(testCase, all(isfinite(S.chiz(:))));
verifyEqual(testCase, S.Bvec, [2.0; 5.5] * [cosd(1) 0 sind(1)], 'RelTol', 1e-12);
end

function test_failed_longitudinal_masks_not_crashes(testCase)
% Spec test 8 (review finding 5): a non-converged longitudinal point must yield a
% masked 1/z column plus an RPA-only overlay -- and must NOT abort the parfor by
% falling through to the strict-paramagnet invz_twolevel gate.
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(1) 0 sind(1)];
fx.solve_opts = struct('max_outer', 1);            % cripple the outer EMT loop
S = invz_spectra_map(ion, 0.31, 2.0, w, fx);
verifyEqual(testCase, S.phase, 0);
verifyTrue(testCase, all(~isfinite(S.chiz(:))));   % 1/z column masked
verifyTrue(testCase, any(isfinite(S.chirpa(:))));  % RPA overlay from the failed pto
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle_spectra.m'); assertSuccess(results)"
```
Expected: FAIL — `field_dir`/`S.Bvec` unknown.

- [ ] **Step 3: Implement**

In `invz/invz_spectra_map.m`, after the existing `getf` option block (after line 40), add:

```matlab
fdir  = getf(opts, 'field_dir', [1 0 0]);
bztol = getf(opts, 'bz_tol', 1e-9);
sxtra = getf(opts, 'solve_opts', struct());
if any(isfield(sxtra, {'J0eff', 'Jxx0', 'hyp'}))
    error('invz:solveOpts', 'solve_opts fields J0eff/Jxx0/hyp are reserved (driver-owned).');
end
if ~isnumeric(fdir) || ~isreal(fdir) || numel(fdir) ~= 3 || ~all(isfinite(fdir)) || norm(fdir(:)) == 0
    error('invz:fieldDir', 'field_dir must be a nonzero finite real 3-vector.');
end
fdir = reshape(fdir, 1, 3) / norm(fdir);
if any(fields < 0)
    error('invz:fields', 'fields are sweep magnitudes |B| and must be nonnegative.');
end
```

After `fields = fields(:).';   w = w(:);` add the dead-banded vector table (same rule as `invz_solve_auto`, applied ONCE so `S.Bvec` is exactly what the solves use — second review finding 2):

```matlab
BvecM = fields(:) * fdir;                        % [nB x 3] actual solve fields
BvecM(abs(BvecM(:, 3)) <= bztol, 3) = 0;         % dead band: identical rule to invz_solve_auto
```

Replace the `sopts` construction inside `one_field` by building it once in the parent (replace nothing yet — see the full `one_field` below) and change the parfor body call to:

```matlab
sopts = sxtra;
sopts.hyp = hyp;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;

parfor (k = 1:nB, nWorkers)
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k)] = ...
        one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol);
    if verbose
        ph = {'masked (no converged solve)', 'moment-form (FM or field-induced)', 'paramagnet'};
        fprintf('  |B| = %5.2f T : %-34s Sigma0 = %s\n', fields(k), ph{phaseC(k)+1}, num2str(Sig0(k)));
    end
end
```

In the `S` packing add:

```matlab
S.field_dir = fdir;  S.Bvec = BvecM;
```

Replace the entire `one_field` local function with:

```matlab
% -------------------------------------------------------------------------------------------
function [chiz, chirpa, Sigma0, phase] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol)
%ONE_FIELD chi''_cc(omega) at one field via the shared ordered-first solve
% (invz_solve_auto); phase = 1 (moment-form: spontaneous FM or field-induced),
% 2 (strict PM), 0 (no accepted solution -> masked 1/z column).
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp);

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- moment-form branch (FM or induced) ---
    o  = invz_chi_realaxis(ion, T, B, pt, w, copts);   % reuses pt.si (moment-form eigenstates)
    chiz = imag(o.chi_cc_q(1, :)).';
    pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                 'K', [], 'is_ordered', true, 'si', pt.si);
    c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;   % share bare cc
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    Sigma0 = pt.Sigma0;
    return;
end

if abs(B(3)) > bztol
    % --- longitudinal failure: NEVER the strict-paramagnet overlay (its m = 0 gate
    % would raise invz:orderedPhase and abort the parfor -- review finding 5). If the
    % failed moment-branch pt still carries valid si/tl, compute the RPA-only overlay
    % from the ordered-style pt0; otherwise leave the whole column masked.
    if ~isempty(pt) && ~isempty(pt.si) && isfield(pt, 'tl') && ~isempty(pt.tl)
        pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                     'K', [], 'is_ordered', true, 'si', pt.si);
        c0opts = copts;  c0opts.npass = 1;
        try
            o0 = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        catch err
            if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        end
    end
    if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
    return;
end

% --- transverse paramagnetic side: unchanged historical logic --------------------------
if phase == 2 && ~isempty(pt), tl0 = pt.tl;  si0 = pt.si;
else, tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));  si0 = []; end
chi0cc = [];
try
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    c0opts = copts;  c0opts.npass = 1;  c0opts.si = si0;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    chi0cc = o0.chi0cc_w;                          % share the bare cc with the 1/z call
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
end
if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    copts1 = copts;
    if ~isempty(chi0cc), copts1.chi0cc_w = chi0cc; end
    o = invz_chi_realaxis(ion, T, B, pt, w, copts1);
    chiz = imag(o.chi_cc_q(1, :)).';
end
end
```

Update the file header: document the three new opts, `S.field_dir`/`S.Bvec`, and reword the phase legend to "1 = moment-form (spontaneous FM below Bc, or field-induced under a longitudinal tilt — a rounded crossover, no sharp Bc), 2 = strict paramagnet, 0 = masked".

- [ ] **Step 4: Run to verify pass, then full fast suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle_spectra.m'); assertSuccess(results)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```
Expected: both PASS (`test_invz_spectra_map` still passes: default `field_dir` reproduces the scalar behavior exactly).

- [ ] **Step 5: Commit**

```bash
git add invz/invz_spectra_map.m invz/tests/test_invz_field_angle_spectra.m
git commit -m "feat(invz): spectra_map field_dir API, solve_opts/bz_tol propagation, longitudinal failure contract"
```

---

### Task 6: `invz_spectra_qpath` — vector field, metadata, safe error formatting

**Files:**
- Modify: `invz/invz_spectra_qpath.m`
- Test: `invz/tests/test_invz_field_angle_spectra.m` (append)

**Interfaces:**
- Consumes: Tasks 1-4.
- Produces: `B` accepts scalar-or-3-vector; new `opts.bz_tol` / `opts.solve_opts` (same contract and reserved-field rule as the map); returns `S.Bvec` (1x3, dead band applied — the vector used) and `S.Bmag = norm(S.Bvec)`; `S.B` keeps returning the normalized vector for continuity. `invz:noSolution` message formats the field with `mat2str`.

- [ ] **Step 1: Write the failing tests** (append to `test_invz_field_angle_spectra.m`)

```matlab
function test_qpath_scalar_vs_vector_and_metadata(testCase)
ion = invz_ion();  w = (0.1:0.1:0.5).';
qp  = [1 0 0; 1.5 0 0; 2 0 0];
fx  = fixture();  fx.dpRng = 8;                    % small real-space cutoff: fast path sums
S1 = invz_spectra_qpath(ion, 0.31, 5.5, qp, w, fx);
S2 = invz_spectra_qpath(ion, 0.31, [5.5 0 0], qp, w, fx);
verifyEqual(testCase, S2.chiz, S1.chiz);
verifyEqual(testCase, S2.Epeak, S1.Epeak);
verifyEqual(testCase, S1.Bvec, [5.5 0 0]);
verifyEqual(testCase, S1.Bmag, 5.5);
end

function test_qpath_error_message_wellformed_for_vector(testCase)
% Review finding 7: a 3-vector B through a scalar %.3f recycles the format string;
% the message must instead carry one mat2str-rendered vector.
ion = invz_ion();  w = (0.1:0.1:0.5).';
qp  = [1 0 0; 2 0 0];
% (1.9 K, 0.04 T): PM two-level raises invz:degenerateDoublet -> phase 0 -> noSolution,
% thrown BEFORE any path lattice sum (cheap).
try
    invz_spectra_qpath(ion, 1.9, [0.04 0 0], qp, w, fixture());
    verifyTrue(testCase, false, 'expected invz:noSolution');
catch err
    verifyEqual(testCase, err.identifier, 'invz:noSolution');
    verifyTrue(testCase, contains(err.message, mat2str([0.04 0 0], 4)));
end
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle_spectra.m'); assertSuccess(results)"
```
Expected: FAIL — no `S.Bvec` field; vector `B` produces a malformed message.

- [ ] **Step 3: Implement**

In `invz/invz_spectra_qpath.m`, after the `getf` option block (after line 49), add:

```matlab
bztol = getf(opts, 'bz_tol', 1e-9);
sxtra = getf(opts, 'solve_opts', struct());
if any(isfield(sxtra, {'J0eff', 'Jxx0', 'hyp'}))
    error('invz:solveOpts', 'solve_opts fields J0eff/Jxx0/hyp are reserved (driver-owned).');
end
B = invz_field_vec(B);
if abs(B(3)) <= bztol, B(3) = 0; end             % same dead band as invz_solve_auto
```

Replace the solve call (line 65) with:

```matlab
sopts = sxtra;
sopts.hyp = true;  sopts.J0eff = Jcc0;  sopts.Jxx0 = Jaa0;  sopts.bz_tol = bztol;
[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);
```

Replace the error (lines 66-70) with:

```matlab
if phase == 0
    error('invz:noSolution', ...
        ['No converged 1/z solution at T = %.3f K, B = %s T ' ...
         '(near-degenerate doublet, critical band, or non-converged moment branch).'], ...
        T, mat2str(B, 4));
end
```

In the `S` packing (line 117), extend:

```matlab
S.w = w;  S.T = T;  S.B = B;  S.Bvec = B;  S.Bmag = norm(B);  S.phase = phase;  S.info = info;  S.Jq = Jq;
```

Header: document scalar-or-vector `B`, the two new opts, `S.Bvec`/`S.Bmag`, and reword the phase note to "(1 = moment-form solve, 2 = strict-PM solve)". The `phase == 1` / `else` overlay branch needs no code change: under a longitudinal `B`, `phase` is 1 or the function has already errored, so the strict-PM `invz_twolevel` line is unreachable with `Bz ~= 0`. State that in a comment above the `if phase == 1` block:

```matlab
% Under a longitudinal B (|Bz| > bz_tol) invz_solve_auto returns phase 1 or the error
% above already fired -- the strict-PM else-branch below is only reachable transversely.
```

- [ ] **Step 4: Run to verify pass, then full fast suite**

Same two commands as Task 5 Step 4. Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_spectra_qpath.m invz/tests/test_invz_field_angle_spectra.m
git commit -m "feat(invz): spectra_qpath vector field, Bvec/Bmag metadata, mat2str noSolution message"
```

---

### Task 7: Driver knob `theta_c` + axis labels

**Files:**
- Modify: `invz/invz_run_spectra.m`, `invz/invz_plot_spectra_map.m:32`
- Test: `invz/tests/test_invz_field_angle_spectra.m` (append — label test only; the driver script itself is exercised manually)

**Interfaces:**
- Consumes: `opts.field_dir` (Task 5), vector `B` in qpath (Task 6).
- Produces: driver knob `theta_c` (deg); `dhat = [cosd(theta_c) 0 sind(theta_c)]`; both plots label the sweep axis `|B| (T)`.

- [ ] **Step 1: Write the failing label test** (append to `test_invz_field_angle_spectra.m`)

```matlab
function test_plot_map_axis_label(testCase)
% Review finding 7: the sweep axis is a magnitude; 'B_x (T)' is wrong under a tilt.
S = struct('fields', [1 2], 'w', (0:0.1:1).');
f = figure('Visible', 'off');
restore = onCleanup(@() close(f));
ax = axes(f);
invz_plot_spectra_map(ax, S, rand(11, 2), 'ttl');
verifyEqual(testCase, ax.XLabel.String, '|B| (T)');
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle_spectra.m'); assertSuccess(results)"
```
Expected: FAIL — label is `B_x (T)`.

- [ ] **Step 3: Implement**

`invz/invz_plot_spectra_map.m` line 32: change `'B_x (T)'` to `'|B| (T)'`.

`invz/invz_run_spectra.m`:

(a) After the `sliceMax` knob (line 45), insert:

```matlab
theta_c = 0;                         % deg -- tilt of the field OUT of the transverse ab-plane
                                     % toward c (the Ising axis): models experimental field
                                     % misalignment. The direction is FIXED across the sweep;
                                     % x-axes stay the total magnitude |B| (the longitudinal
                                     % component is |B|*sind(theta_c)). theta_c = 0 reproduces
                                     % the pure transverse benchmark exactly. Convention matches
                                     % LiReF4_MF_Yikai theta at phi = 0 (spec 2026-07-16).
                                     % SCALAR-STAGE VALIDITY: Sigma dresses the cc channel only;
                                     % exact at theta_c = 0, uncontrolled O(theta^2) cross-
                                     % channel error beyond the tensor-referenced small-tilt
                                     % range (invz_run_tensor_ref); a longitudinal component
                                     % turns the sharp transition into a rounded crossover.
                                     % Azimuth (phi_ab) + full-tensor propagation: deferred.
```

(b) After the `eScale` switch block (after line 72), insert:

```matlab
dhat = [cosd(theta_c) 0 sind(theta_c)];          % unit field direction (ac-plane)
tiltStr = '';
if theta_c ~= 0, tiltStr = sprintf(', \\theta_c = %.2g\\circ', theta_c); end
```

(c) q-path view — line 83 and line 96: replace `Bq` / `Bq(k)` in the `invz_spectra_qpath(...)` calls with `Bq*dhat` / `Bq(k)*dhat` (titles keep the magnitude `Bq`); append `tiltStr` inside the two `sprintf` titles of the single-`Bq` branch and the multi-`Bq` title, e.g.:

```matlab
sprintf('1/z FM-mode \\chi''''_{cc}, T = %.2f K, B = %.2f T%s', T, Bq, tiltStr)
```

(d) Field-sweep view — line 113: add the direction to the map opts:

```matlab
S = invz_spectra_map(ion, T, fields, wMeV, ...
        struct('parallel', useParallel, 'eta', eta, 'field_dir', dhat));
```

(e) Peak plot — line 134: `xlabel('B_x (T)');` → `xlabel('|B| (T)');` and append `tiltStr` to the peak-plot and slice-view titles (`sprintf('T = %.2f K%s', T, tiltStr)`).

(f) Header comment: one added line under the view list — "theta_c (deg) tilts the swept field toward c to model misalignment; see the knob comment for validity."

- [ ] **Step 4: Run to verify pass, then full fast suite**

Same two commands as Task 5 Step 4. Expected: both PASS.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_run_spectra.m invz/invz_plot_spectra_map.m invz/tests/test_invz_field_angle_spectra.m
git commit -m "feat(invz): theta_c misalignment knob in invz_run_spectra; |B| axis labels"
```

---

### Task 8: Symmetry and continuity tests (fast mirror + slow crossover/continuity)

**Files:**
- Test: `invz/tests/test_invz_field_angle_spectra.m` (append: fast mirror test)
- Test: `invz/tests/test_invz_field_angle_slow.m` (new, INVZ_SLOW-gated)

**Interfaces:** consumes everything above; produces no new code — pure validation.

- [ ] **Step 1: Add the fast ±Bz mirror test** (append to `test_invz_field_angle_spectra.m`)

```matlab
function test_pm_bz_mirror_symmetry(testCase)
% Spec test 4: chi''_cc is even in Bz (Z2). Metric per second-review refinement 2.
ion = invz_ion();  w = (0:0.02:0.5).';
fp = fixture();  fp.field_dir = [cosd(1) 0 +sind(1)];
fm = fixture();  fm.field_dir = [cosd(1) 0 -sind(1)];
Sp = invz_spectra_map(ion, 0.31, 5.0, w, fp);
Sm = invz_spectra_map(ion, 0.31, 5.0, w, fm);
verifyEqual(testCase, Sp.phase, 1);
verifyEqual(testCase, Sm.phase, 1);
a = Sp.chiz(:);  b = Sm.chiz(:);
verifyTrue(testCase, all(isfinite(a)) && all(isfinite(b)));
verifyLessThan(testCase, max(abs(a - b)) / max(max(abs(a)), 1e-12), 1e-8);
end
```

- [ ] **Step 2: Create the slow test file**

Create `invz/tests/test_invz_field_angle_slow.m`:

```matlab
function tests = test_invz_field_angle_slow
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', ...
            'info', struct('Jcc0', 6.4e-3), 'verbose', false);
end

function test_crossover_continuity(testCase)
% Spec test 5: under a 0.5-degree tilt the former sharp transition is a rounded
% crossover -- no censored/masked points across the old Bc, sum rule intact.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'slow test: set INVZ_SLOW=1');
ion = invz_ion();  T = 0.31;  w = (0:0.02:0.5).';
fields = 4.6:0.05:5.3;
fx = fixture();  fx.field_dir = [cosd(0.5) 0 sind(0.5)];
S = invz_spectra_map(ion, T, fields, w, fx);
verifyTrue(testCase, all(S.phase == 1));                     % moment-form everywhere
verifyTrue(testCase, all(isfinite(S.Epeak)));                % no censored crossover gap
% sum rule at every field (pt-level; same fixture and route as the map)
sopts = struct('hyp', true, 'J0eff', fx.info.Jcc0, 'Jxx0', ion.Jxx0);
for k = 1:numel(fields)
    [pt, phase] = invz_solve_auto(ion, T, fields(k)*fx.field_dir, fx.Jnu, sopts);
    verifyEqual(testCase, phase, 1, sprintf('field %.2f', fields(k)));
    verifyLessThan(testCase, pt.sumrule_rel, 5e-2, sprintf('field %.2f', fields(k)));
end
end

function test_theta_to_zero_continuity(testCase)
% Spec test 6 (amended 2026-07-16, recorded justification -- see spec):
%   B = 6 T (PM at zero tilt): chi_cc even in Bz -> flat 1e-6 bound (measured
%     4.2e-7 at 1e-3 deg). Exercises forced moment-form -> strict-PM reduction.
%   B = 3 T (spontaneous FM at zero tilt; 2 T sits in a fixture-specific
%     non-convergence island B in [1,2] T): the single-domain response is
%     LINEAR in Bz (aligned branch breaks Z2; soft-mode coefficient ~ delta/eta,
%     measured 7.7e-3 at 1e-3 deg) -- a flat 1e-6 bound is the wrong physics.
%     Assert continuity as linear scaling: no O(1) jump at the spontaneous ->
%     forced routing boundary.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'slow test: set INVZ_SLOW=1');
ion = invz_ion();  T = 0.31;  w = (0:0.02:0.5).';
% --- PM side: B = 6 T, flat bound ---
S0 = invz_spectra_map(ion, T, 6.0, w, fixture());
ft = fixture();  ft.field_dir = [cosd(1e-3) 0 sind(1e-3)];
St = invz_spectra_map(ion, T, 6.0, w, ft);
a = S0.chiz(:);  b = St.chiz(:);
verifyTrue(testCase, all(isfinite(a)) && all(isfinite(b)), 'B = 6');
verifyLessThan(testCase, max(abs(b - a)) / max(abs(a)), 1e-6, 'B = 6');
% --- FM side: B = 3 T, linear-scaling assertion ---
S0 = invz_spectra_map(ion, T, 3.0, w, fixture());
angs = [1e-3 1e-4];  d = nan(1, 2);
for k = 1:2
    ft = fixture();  ft.field_dir = [cosd(angs(k)) 0 sind(angs(k))];
    St = invz_spectra_map(ion, T, 3.0, w, ft);
    a = S0.chiz(:);  b = St.chiz(:);
    verifyTrue(testCase, all(isfinite(a)) && all(isfinite(b)), sprintf('B = 3, %g deg', angs(k)));
    d(k) = max(abs(b - a)) / max(abs(a));
end
verifyLessThan(testCase, d(1), 3e-2);
r = d(1) / d(2);                                   % pure linear response -> 10
verifyGreaterThan(testCase, r, 6);
verifyLessThan(testCase, r, 15);
end
```

- [ ] **Step 3: Run the fast file, then the slow file explicitly**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle_spectra.m'); assertSuccess(results)"
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_field_angle_slow.m'); assertSuccess(results)"
```
Expected: both PASS. **If the mirror bound (1e-8) or continuity bound (1e-6) fails by a small factor:** do NOT silently loosen — record the measured value and the loosening justification in the commit message (spec: "refined only with recorded justification"); investigate first whether the miss is a real asymmetry (e.g. a missed seed mirror) rather than float noise.

- [ ] **Step 4: Full fast suite, then commit**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
git add invz/tests/test_invz_field_angle_spectra.m invz/tests/test_invz_field_angle_slow.m
git commit -m "test(invz): Bz mirror, crossover continuity, theta->0 reduction tests"
```

---

### Task 9: Σ=0 scalar-vs-tensor reference (`invz_chi_tensor_ref`) + supported-angle measurement

**Files:**
- Create: `invz/invz_chi_tensor_ref.m`, `invz/invz_tilt_err.m`, `invz/invz_run_tensor_ref.m`
- Create: `invz/tests/test_invz_tensor_ref.m`, `docs/SESSION-2026-07-16-field-angle.md`

**Interfaces:**
- Consumes: `invz_field_vec`, `invz_single_ion` (order mode, sign-aware), `invz_chi0z`, `invz_peak_energy`, `getf`.
- Produces: `R = invz_chi_tensor_ref(ion, T, Bvec, w, opts)` with fields
  `R.chi_sc`, `R.chi_ten` (`nw x 1` χ''_cc), `R.eps_spec` (floored spectral L2 relative error), `R.Epeak_sc`, `R.Epeak_ten`, `R.dE_peak`, `R.w`, `R.B`.
  Opts: `eta` (5e-3), `Jsel` (`ion.J0eff`), `Jaa0` (`ion.Jxx0`), `hyp` (true), `peak_wmin` (0.05).

- [ ] **Step 1: Write the failing structural tests**

Create `invz/tests/test_invz_tensor_ref.m`:

```matlab
function tests = test_invz_tensor_ref
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_cross_channel_hierarchy_at_zero_tilt(testCase)
% Measured physics (Task-9 blocked round, amended 2026-07-16): at theta = 0,
% B || a, the yz cross channel of chi0 is SYMMETRY-ALLOWED (B64s) and large,
% while xz is strongly suppressed (measured yz/zz = 0.183, xz/zz = 2.8e-3 at
% 6 T, hyp = false). The scalar-vs-tensor baseline discrepancy at zero tilt is
% therefore real and finite -- it is REPORTED as a baseline, never gated as a
% tilt error.
ion = invz_ion();
si = invz_single_ion(ion, 0.31, [6 0 0], struct('hyp', false));
z = (0:0.01:0.5).' + 1i*5e-3;
c0 = invz_chi0z(si, 0.31, z, struct('elastic', false));
xz = max(abs(squeeze(c0(1,3,:))));  yz = max(abs(squeeze(c0(2,3,:))));
zz = max(abs(squeeze(c0(3,3,:))));
verifyGreaterThan(testCase, yz, 10*xz);          % yz-dominated cross-channel hierarchy
verifyLessThan(testCase, xz, 0.02*zz);           % xz strongly suppressed at theta = 0
w = (0:0.01:0.5).';
R = invz_chi_tensor_ref(ion, 0.31, [6 0 0], w, struct('hyp', false, 'eta', 0.02));
verifyGreaterThan(testCase, R.eps_spec, 0);      % finite yz-driven baseline ...
verifyLessThan(testCase, R.eps_spec, 1);         % ... bounded (sanity)
end

function test_tilt_metric_wellformed(testCase)
% eps_tilt (invz_tilt_err) is the GATED metric: the error in the tilt-induced
% change, baseline-differenced so the theta-independent yz discrepancy drops out.
ion = invz_ion();
w = (0:0.01:0.5).';
o = struct('hyp', false, 'eta', 0.02);
R0 = invz_chi_tensor_ref(ion, 0.31, [6 0 0], w, o);
R1 = invz_chi_tensor_ref(ion, 0.31, 6*[cosd(1) 0 sind(1)], w, o);
e1 = invz_tilt_err(R1, R0);
verifyGreaterThan(testCase, e1, 0);
verifyLessThan(testCase, e1, 2);
end

function test_demag_guard(testCase)
ion = invz_ion();  ion.demag = 1;
verifyError(testCase, ...
    @() invz_chi_tensor_ref(ion, 0.31, [6 0 0], (0:0.1:0.4).'), 'invz:tensorRef');
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_tensor_ref.m'); assertSuccess(results)"
```
Expected: FAIL — `invz_chi_tensor_ref` undefined.

- [ ] **Step 3: Implement the reference**

Create `invz/invz_chi_tensor_ref.m`:

```matlab
function R = invz_chi_tensor_ref(ion, T, Bvec, w, opts)
%INVZ_CHI_TENSOR_REF Sigma = 0 scalar-vs-tensor RPA cross-check at one (T, B).
% Quantifies the cross-channel (xz/yz/xy) error the scalar cc pipeline omits under
% a c-axis tilt (spec 2026-07-16 SR3). Both sides share ONE single-ion state
% (order mode, sign-aware seed) and ONE bare chi0 tensor from invz_chi0z; they
% differ ONLY in propagation:
%   scalar:  chi_cc = chi0_cc / (1 - Jsel*chi0_cc)              (cc channel only)
%   tensor:  chi    = chi0 * inv(I3 - J*chi0),  J = diag(Jaa0, Jaa0, Jsel)
% NOTE: in the SPONTANEOUSLY ordered phase (m0 ~= 0 at theta = 0) the cross
% channels are nonzero even without tilt -- a pre-existing property of the scalar
% pipeline that this reference makes visible; expect eps_spec(theta=0) > 0 there.
% Intrinsic response only: asserts ion.demag == 0.
%
% R fields: chi_sc, chi_ten [nw x 1] (chi''_cc), eps_spec (floored spectral L2
% relative error), Epeak_sc/Epeak_ten/dE_peak (invz_peak_energy, opts.peak_wmin
% default 0.05 meV to skip the hyperfine line), w, B.
if nargin < 5, opts = struct(); end
eta   = getf(opts, 'eta', 5e-3);
Jsel  = getf(opts, 'Jsel', ion.J0eff);
Jaa0  = getf(opts, 'Jaa0', ion.Jxx0);
hyp   = getf(opts, 'hyp', true);
wmin  = getf(opts, 'peak_wmin', 0.05);
if ion.demag ~= 0
    error('invz:tensorRef', 'reference defined for the intrinsic response (ion.demag = 0).');
end
B  = invz_field_vec(Bvec);
si = invz_single_ion(ion, T, B, struct('hyp', hyp, 'order', true, 'J0z', Jsel, 'Jxx0', Jaa0));
w  = w(:);
z  = w + 1i*eta;
c0 = invz_chi0z(si, T, z, struct('elastic', false));
J  = diag([Jaa0, Jaa0, Jsel]);
nw = numel(z);
chi_sc = zeros(nw, 1);  chi_ten = zeros(nw, 1);
for k = 1:nw
    X = c0(:, :, k);
    chi_sc(k)  = imag(X(3, 3) / (1 - Jsel*X(3, 3)));
    Xt = X / (eye(3) - J*X);                     % chi0 * inv(I - J*chi0), full 3x3
    chi_ten(k) = imag(Xt(3, 3));
end
R.chi_sc = chi_sc;  R.chi_ten = chi_ten;  R.w = w;  R.B = B;
floorv = 1e-12 * max(abs(chi_ten)) * sqrt(nw);   % guards the metric at spectral zeros
R.eps_spec  = norm(chi_sc - chi_ten, 2) / max(norm(chi_ten, 2), floorv);
R.Epeak_sc  = invz_peak_energy(chi_sc,  w, wmin);
R.Epeak_ten = invz_peak_energy(chi_ten, w, wmin);
R.dE_peak   = abs(R.Epeak_sc - R.Epeak_ten);
% Peak-observable amplitude error (the GATED intensity metric; L2 lineshape
% metrics are positional-artifact-dominated for sharp lines -- see spec sec. 7).
% Same wmin mask as the peak search: the gate must measure the ELECTRONIC mode,
% not a sub-wmin hyperfine feature.
msk = w >= wmin;
R.amp_sc  = max(chi_sc(msk));
R.amp_ten = max(chi_ten(msk));
R.eps_amp = abs(R.amp_sc - R.amp_ten) / max(R.amp_ten, floorv);
end
```

- [ ] **Step 3b: Create the tilt-metric helper**

Create `invz/invz_tilt_err.m`:

```matlab
function eps = invz_tilt_err(R, R0)
%INVZ_TILT_ERR Tilt-referenced spectral error of the scalar cc chain.
% Measures how well the scalar pipeline captures the TILT-INDUCED change (the
% spec's accuracy statement), by differencing against the theta = 0 reference
% R0 at the same field -- this removes the theta-independent baseline
% discrepancy from the symmetry-allowed yz cross channel (B64s), which exists
% at zero tilt and is not a tilt error.
Dsc  = R.chi_sc  - R0.chi_sc;
Dten = R.chi_ten - R0.chi_ten;
floorv = 1e-12 * max(abs(R.chi_ten)) * sqrt(numel(R.chi_ten));
eps = norm(Dsc - Dten, 2) / max(norm(Dten, 2), floorv);
end
```

- [ ] **Step 4: Run structural tests to verify pass**

Same command as Step 2. Expected: PASS (3 tests). The hierarchy test's loose
bounds encode measured physics (yz/zz = 0.183, xz/zz = 2.8e-3); if they fail by
orders of magnitude, suspect the chi0 evaluation, not the bounds.

- [ ] **Step 5: Write the measurement driver**

Create `invz/invz_run_tensor_ref.m`:

```matlab
%INVZ_RUN_TENSOR_REF Sigma=0 scalar-vs-tensor cross-channel error vs tilt angle.
% Two-layer metric (spec section 7, amended after the Task-9 measurement):
%   eps_spec (invz_chi_tensor_ref): RAW discrepancy incl. the theta-independent
%     baseline from the symmetry-allowed yz cross channel (B64s; measured
%     yz/zz = 0.183 at 6 T even at theta = 0). REPORTED, never gated.
%   eps_tilt (invz_tilt_err): error in the TILT-INDUCED change, differenced
%     against the theta = 0 reference at the same field. Diagnostic only.
% Comparison eta = 0.02 meV (4 pts/HWHM on the 0.005 grid): at the production
% eta = 5e-3 the L2 norm is dominated by sub-linewidth peak misalignment (a
% metric instability, not physics).
% GATE (final, user-approved): peak observables only -- every L2 lineshape
% metric is dominated by the zero-tilt yz peak-offset artifact (delta0/eta;
% verified 0.11 vs 2.2/20 at 6T, ~0.28 vs 7.7/20 at 2T). eps_spec/eps_tilt are
% reported diagnostics; lineshape fidelity is explicitly not claimed.
%   supported(theta > 0) <=> eps_amp <= 0.10 AND dE_peak <= max(0.02*Epeak_ten, eta)
%   at EVERY tested field. Copy the printed table into
%   docs/SESSION-2026-07-16-field-angle.md and the constants of
%   invz/tests/test_invz_tensor_ref.m (reproducibility assertion).
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
ion = invz_ion();
T       = 0.1;                       % K -- the spectra-driver default temperature
angles  = [0 0.5 1 2 5];             % deg (spec grid; 0 = baseline row, not gated)
fieldsB = [2 4.95 6];                % T: ordered / near-crossover / paramagnetic
w       = (0:0.005:0.6).';           % meV
eta     = 0.02;                      % meV comparison broadening (metric stability)

% live couplings, as in the production drivers (cached lattice sum)
[qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
qc = qc(any(abs(qc) > 1e-12, 2), :);
[~, info] = invz_jq_modes(ion, qc, struct('dpRng', 30, 'cache', true));
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end
ropts = struct('Jsel', info.Jcc0, 'Jaa0', Jaa0, 'eta', eta);

fprintf('%8s %8s %12s %12s %12s %12s %10s %10s\n', 'theta', '|B| (T)', 'eps_spec', 'eps_tilt', 'eps_amp', 'dE_peak', 'Ep_sc', 'Ep_ten');
supported = true(size(angles));
for ib = 1:numel(fieldsB)
    B = fieldsB(ib);
    R0 = invz_chi_tensor_ref(ion, T, [B 0 0], w, ropts);   % theta = 0 reference
    for ia = 1:numel(angles)
        a = angles(ia);
        if a == 0
            R = R0;  et = 0;
        else
            R = invz_chi_tensor_ref(ion, T, B*[cosd(a) 0 sind(a)], w, ropts);
            et = invz_tilt_err(R, R0);
        end
        ok = a == 0 || (R.eps_amp <= 0.10 && ...
             ( (isnan(R.dE_peak) && isnan(R.Epeak_sc) == isnan(R.Epeak_ten)) || ...
               R.dE_peak <= max(0.02*R.Epeak_ten, eta) ));
        if a > 0, supported(ia) = supported(ia) && ok; end
        verdict = {'FAIL', 'ok'};
        fprintf('%8.2f %8.2f %12.4g %12.4g %12.4g %12.4g %10.4g %10.4g   %s\n', ...
                a, B, R.eps_spec, et, R.eps_amp, R.dE_peak, R.Epeak_sc, R.Epeak_ten, verdict{ok + 1});
    end
end
fprintf('Supported tilt range (peak-observable criterion): theta_c <= %.2g deg\n', ...
        max([0, angles(supported & angles > 0)]));
```

- [ ] **Step 6: Run the measurement and record the results**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "run('invz/invz_run_tensor_ref.m')"
```

This step CAPTURES DATA (the numbers cannot be known before the run):
1. Create `docs/SESSION-2026-07-16-field-angle.md` containing: the printed table verbatim, the support criterion, the resulting supported `theta_c` range, and one paragraph noting the expected `eps_spec(θ=0) > 0` at the ordered field (pre-existing scalar-pipeline property, per the `invz_chi_tensor_ref` header).
2. Append the reproducibility test to `test_invz_tensor_ref.m`, substituting the measured `eps_spec` values from the table into `expected`:

```matlab
function test_reproducibility_of_logged_table(testCase)
% Slow: re-measures three logged (angle, field) points and asserts the eps_tilt
% values match docs/SESSION-2026-07-16-field-angle.md to 1% (reproducibility,
% NOT a size target).
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'slow test: set INVZ_SLOW=1');
ion = invz_ion();
w = (0:0.005:0.6).';
o = struct('eta', 0.02);      % default couplings -- must match the logged convention
pts          = {0.5, 6;  2, 6;  5, 4.95};        % {theta_deg, B} spot checks
expected_amp  = [NaN; NaN; NaN];  % <- REPLACE with measured eps_amp (GATED metric)
expected_tilt = [NaN; NaN; NaN];  % <- REPLACE with measured eps_tilt (diagnostic)
                                  %    from the first run BEFORE committing
                                  %    (committing NaNs = task incomplete)
for k = 1:size(pts, 1)
    a = pts{k, 1};  B = pts{k, 2};
    R0 = invz_chi_tensor_ref(ion, 0.1, [B 0 0], w, o);
    R  = invz_chi_tensor_ref(ion, 0.1, B*[cosd(a) 0 sind(a)], w, o);
    verifyEqual(testCase, R.eps_amp, expected_amp(k), 'RelTol', 0.01);
    verifyEqual(testCase, invz_tilt_err(R, R0), expected_tilt(k), 'RelTol', 0.01);
end
end
```

(Note: the spot-check calls use default `Jsel`/`Jaa0` — re-run those three points once with defaults to fill `expected`, or pass the live couplings in the test exactly as the driver does; either is fine as long as the test and its recorded constants use the SAME couplings. State which was used in the session doc.)

- [ ] **Step 7: Run all tensor-ref tests, full fast suite, commit**

```bash
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_tensor_ref.m'); assertSuccess(results)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
git add invz/invz_chi_tensor_ref.m invz/invz_run_tensor_ref.m invz/tests/test_invz_tensor_ref.m docs/SESSION-2026-07-16-field-angle.md
git commit -m "feat(invz): Sigma=0 scalar-vs-tensor reference; measured supported tilt range"
```

---

### Task 10: Documentation + final gates

**Files:**
- Modify: `invz/README.html`
- Modify: `docs/SESSION-2026-07-16-field-angle.md` (extend)

- [ ] **Step 1: README.html updates**

Match the existing HTML structure (the file has a function-reference section and numbered feature sections; mirror the markup of neighboring entries):

1. Function reference: add entries for `invz_field_vec` ("scalar/3-vector field normalization; every solver boundary accepts `[Bx By Bz]`"), `invz_chi_tensor_ref` + `invz_run_tensor_ref` ("Σ=0 scalar-vs-tensor cross-channel reference; sets the supported tilt range").
2. Spectra-driver section: a short paragraph — "`theta_c` tilts the swept field toward c (misalignment). Longitudinal fields route through the moment-form self-energy (`forced_moment`; sign-aware branch seeding); the sharp transition becomes a rounded crossover, `S.phase = 1` meaning moment-form (spontaneous OR field-induced). Scalar-stage validity: cc-channel-only dressing, exact at zero tilt, tensor-referenced supported range `theta_c <= <measured value> deg` (see `docs/SESSION-2026-07-16-field-angle.md`). Azimuth and full-tensor propagation are deferred (spec §8)."
3. Scope section: move "longitudinal/tilted field" from out-of-scope to in-scope-(scalar-stage, c-tilt only), citing the spec.

Insert the measured supported range from Task 9's table — no placeholder text.

- [ ] **Step 2: Extend the session doc**

Append to `docs/SESSION-2026-07-16-field-angle.md`: the file/function map of this feature (table from this plan's File Structure), the new option/metadata names (`field_dir`, `bz_tol`, `solve_opts`, `moment_branch`, `S.Bvec`), and the deferred items (azimuth + two-channel MF; A0/A1) with a pointer to the spec.

- [ ] **Step 3: Final verification gates**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```
Expected: both PASS (slow run includes the existing INVZ_SLOW benchmarks — this is the spec's bit-for-bit non-negotiable gate).

Manual smoke (optional but recommended; needs a display-less figure export): set `theta_c = 1;` and `fields = linspace(4.4, 5.6, 25);` in `invz_run_spectra.m`, run it once, and confirm the map shows a rounded crossover (no masked band at the former Bc) — then revert the knob edits to `theta_c = 0` and the original `fields`.

- [ ] **Step 4: Commit**

```bash
git add invz/README.html docs/SESSION-2026-07-16-field-angle.md
git commit -m "docs(invz): README + session notes for the theta_c misalignment knob"
```
