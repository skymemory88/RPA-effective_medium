# Demag Knob Rework + q-Path Spectra (R2007 Fig. 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the sample-shape (demagnetization) correction physically correct and user-toggleable in both run scripts, and add a q-path spectrum mode reproducing Rønnow et al. PRB 75, 054426 (2007) Fig. 3.

**Architecture:** The ordering-channel coupling (`info.Jcc0`, `Jnu`) becomes demag-invariant (per R2007 the demag field cancels from the critical condition); the shape enters only through (a) a new strict-uniform observable correction `info.Jshape_cc` applied inside `invz_chi_realaxis`, and (b) a new live transverse coupling `info.Jaa0` (demag-aware) threaded as `opts.Jxx0` through every single-ion solve. A new `invz_spectra_qpath` reuses the once-solved 1/z medium plus the vectorized `Jsel` machinery of `invz_chi_realaxis` to evaluate chi''(q,w) along any r.l.u. path.

**Tech Stack:** MATLAB (R2025a), matlab.unittest `functiontests` suite in `invz/tests/`, no toolboxes required (Parallel Computing Toolbox optional).

---

## Context and physics rationale (read before executing)

- **Working directory / repo root:** `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion`. All paths below are relative to it. Current branch: `invz-1z-lihof4`. MATLAB invocation (from repo root):
  `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"`
  Single file: replace `'invz/tests'` with e.g. `'invz/tests/test_invz_jq_modes.m'`. Slow tests additionally need env `INVZ_SLOW=1`.
- **The bug being fixed** (Code_review_byCodex.md Finding #1, High): the current working tree folds the demag term into `info.Jcc0`, which feeds the mean field, the RPA denominator, AND the 1/z critical condition. R2007 (p. 054426-2/3) states verbatim: *"the demagnetization field cancels out in the determination of the critical condition, no matter the shape of the individual domains in the ordered phase … the ordering is going to occur at wave vectors infinitesimally different from zero, not at q ≡ 0."* Note R2007 says this **against** a contrary suggestion by Chakraborty et al. (their ref. 17) — cite it as "per R2007" not "per settled theory".
- **Where demag legitimately DOES act** (both per R2007's own usage):
  1. The strict-uniform **measured observable**: `chi_meas = chi_int / (1 + Jshape*chi_int)` with `Jshape = 4*dm_cc` (the factor 4 is the uniform-mode projection of a scalar-broadcast term onto the 4-site sublattice vector — same reason `info.Jcc0 = Jcc0_dipole + 4*J12`). Algebraically this reproduces the old shifted-pole result `chit/(1-(Jcc0-Jshape)*chit)` exactly, so the observable is unchanged relative to the old convention — only criticality/Sigma/mean-field change.
  2. The **transverse (aa) channel**: `hx = Jxx0*<Jx>` is not a critical/order-parameter channel (no q→0+ cancellation applies); a shape term there is the internal-vs-applied transverse field correction. Decision (user-approved): `info.Jaa0` is demag-aware directly. Consequence: with demag on, `Bc(T)` shifts **only** via this transverse channel (physical: boundary vs applied field for that shape); `Tc(B=0)` never shifts (at B=0, `<Jx>=0`).
- **Lorentz term:** stays unconditionally on, no knob (user-approved). It is the mandatory cavity term of the dipole-sum split, not a physical toggle.
- **Expected small numeric drift:** drivers switch the transverse coupling from the hardcoded `ion.Jxx0 = 3.512e-3` meV to the live `info.Jaa0 ≈ 3.5104e-3` meV (dpRng=30). ~0.05% shift in the transverse channel → sub-mK/mT boundary drift even at demag=0. Intended (self-consistency); slow-test tolerances absorb it. `ion.Jxx0` remains the fallback default inside `invz_single_ion` so direct callers are unaffected.
- **Cache invalidation:** `invz_jq_modes` results change schema (new `info.Jaa0`, `info.Jshape_cc`; demag removed from `Jcc0`). The cache key gets a schema version (filename prefix `jq2_`), so all old `invz/cache/jq_*.mat` files are simply never matched again (safe to leave or delete).
- **q-path spectra are intrinsic:** a finite-q probe (neutrons) never sees the shape correction; at Γ-equivalent path points the relevant limit is q→0+. So `invz_spectra_qpath` must NOT apply `Jshape`. Sanity anchors from the 4-site structure factor: `(2,0,0)` IS Γ-equivalent (structure factor 4 → Lorentz applies, max branch = `info.Jcc0`), `(1,0,0)` is NOT (structure factor 0). This matches R2007 Fig. 3: the mode softens toward `(2,0,0)` and is stiff at `(1,0,0)`.
- **Peak picking near criticality:** the hyperfine interaction adds a low-frequency pole (~0.01 meV, R2007 Fig. 2 discussion). Fig. 3 tracks the crystal-field doublet excitation, so peak extraction excludes `w < opts.peak_wmin` (default 0.05 meV).
- **Out of scope (do not attempt):** off-diagonal dipolar (ODD/full-tensor) terms; an aa-channel observable; the `hyp=false` real-axis mismatch (Codex #10); BZ-grid endpoint bias (Codex #2); causality/positivity of the real-axis continuation (Codex #3); the FM/PM branch-selection boundary (Codex #4).
- **Dirty working tree warning:** the tree already contains uncommitted work (eUnit driver rewrite, ordered-phase files, the demag draft this plan reworks). Task 0 checkpoints it so per-task diffs stay reviewable. Never `git add -A` after Task 0 — always add explicit paths. Never commit `invz/invz_run_spectra.asv` (editor autosave).

### File map

| File | Action | Role |
|---|---|---|
| `invz/invz_jq_modes.m` | modify | single source of J(q): demag-free `Jcc0`/`Jnu`, new `info.Jaa0`, `info.Jshape_cc`, cache v2 |
| `invz/invz_twolevel.m`, `invz/invz_twolevel_ordered.m` | modify | optional `opts.Jxx0` forwarding |
| `invz/invz_solve_point.m`, `invz/invz_solve_point_ordered.m` | modify | forward `opts.Jxx0` |
| `invz/invz_chi_realaxis.m` | modify | forward `opts.Jxx0`; apply `opts.Jshape` observable correction |
| `invz/invz_solve_auto.m` | create | shared ordered-first/para-fallback solve (extracted from `invz_spectra_map`) |
| `invz/invz_spectra_map.m` | modify | thread `Jaa0`/`Jshape`; `one_field` uses `invz_solve_auto` |
| `invz/invz_spectra_qpath.m` | create | chi''(q,w) along a path at fixed (T,B) |
| `invz/invz_plot_spectra_qpath.m` | create | colormap panel for the q-path view |
| `invz/invz_run_spectra.m` | modify | corrected demag knob comment; qpath/Bq/dispScale knobs + third view |
| `invz/invz_run_phase_diagram.m` | modify | demag knob comment; hoist `Jxx0 = info.Jaa0` into solver opts |
| `invz/invz_ion.m` | modify | demag comment corrected |
| `invz/tests/test_invz_jq_modes.m` | modify | rewrite demag test (invariance), `Jaa0` assertions |
| `invz/tests/test_invz_single_ion.m` | modify | add `Jxx0` override test |
| `invz/tests/test_invz_chi_observable.m` | modify | add observable-rescale identity test |
| `invz/tests/test_invz_demag_invariance.m` | create | Tc0/Sigma_c invariance; pinned-Jxx0 crit invariance (SLOW) |
| `invz/tests/test_invz_spectra_qpath.m` | create | structural + Γ-limit anchors for the q-path |
| `invz/README.html` | modify | demag semantics, new functions, qpath quick-start |
| `.gitignore` | create/modify | ignore `*.asv` |

---

### Task 0: Checkpoint the dirty working tree

The tree has pre-existing uncommitted work this plan builds on. Commit it as a single checkpoint so every later task produces a clean, reviewable diff.

**Files:** everything currently modified/untracked EXCEPT `invz/invz_run_spectra.asv`.

- [ ] **Step 1: Ignore MATLAB autosaves**

Create (or append to) `.gitignore` at repo root:

```gitignore
*.asv
```

- [ ] **Step 2: Review what will be committed**

```bash
git status
git diff --stat
```

Confirm nothing looks like a secret/credential. Expected: modified `invz/*.m`, `invz/tests/*.m`, `invz/README.html`, `docs/SESSION-2026-07-07-invz-handoff.md`, `jensen_1z_framework.html`, deleted `invz/README.md` and the HoF3 PDF, plus untracked `Code_review_byCodex.md`, `Data/`, `References/`, new `invz/*.m` ordered-phase files, HTML docs, `verify_*.py`.

Note: `References/` and `Data/` may contain large binary PDFs/data. If any single file exceeds ~50 MB, leave that file uncommitted and mention it in the task report instead.

- [ ] **Step 3: Commit the checkpoint**

```bash
git add -A
git reset -- invz/invz_run_spectra.asv
git commit -m "wip: checkpoint working tree before demag rework (eUnit driver, ordered-phase files, demag draft, refs)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

- [ ] **Step 4: Verify the fast test suite is green at the baseline**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: all pass (~32 passed, ~8 filtered as SLOW). If the baseline is already red, STOP and report — do not proceed on a broken baseline.

---

### Task 1: `invz_jq_modes` — demag-free ordering channel, `info.Jaa0`, `info.Jshape_cc`, cache v2

**Files:**
- Modify: `invz/invz_jq_modes.m`
- Test: `invz/tests/test_invz_jq_modes.m`

- [ ] **Step 1: Rewrite the demag test to assert the NEW physics (failing first)**

In `invz/tests/test_invz_jq_modes.m`, replace the entire `test_demag_shape_term` function (currently lines 30–54) with:

```matlab
function test_demag_shape_term(testCase)
% Ordering-channel couplings are demag-INVARIANT: per R2007 the demagnetization
% field cancels from the critical condition (ordering occurs at q -> 0+, not the
% strict uniform mode), so Jcc0 and the branch spectrum Jnu never see the sample
% shape. The shape enters only through (a) info.Jshape_cc, the strict-uniform
% OBSERVABLE correction, and (b) info.Jaa0, the transverse (non-critical) channel.
ion0 = invz_ion();                                  % demag = 0 default
qs = [0 0 0; 1 0 0; 0.25 0 0.1];                    % Gamma, non-Gamma zone point, generic q
[J0nu, i0] = invz_jq_modes(ion0, qs, struct('dpRng',30,'cache',false));
verifyEqual(testCase, i0.Jshape_cc, 0);             % off: no shape correction
verifyEqual(testCase, i0.Jaa0, i0.Jaa0_dipole + 4*ion0.J12);

C = invz_const();  lorz4 = 4 * (4*pi/(3*ion0.Vc)*C.gfac);  % uniform-mode Lorentz share

% sphere (Nz = Nx = 1/3): ordering channel must NOT move; observable and
% transverse channels must. For a sphere dm_cc = dm_aa = lorz, so
% Jshape_cc = 4*lorz and the transverse dipole share drops by the same amount.
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
[JSnu, iS] = invz_jq_modes(ionS, qs, struct('dpRng',30,'cache',false));
verifyEqual(testCase, JSnu, J0nu);                  % branch spectrum bit-identical
verifyEqual(testCase, iS.Jcc0, i0.Jcc0);            % criticality coupling bit-identical
verifyEqual(testCase, iS.Jcc0_dipole, i0.Jcc0_dipole);
verifyEqual(testCase, iS.Jshape_cc, lorz4, 'RelTol', 1e-12);
verifyEqual(testCase, iS.Jaa0_dipole, i0.Jaa0_dipole - lorz4, 'RelTol', 1e-9);
verifyEqual(testCase, iS.Jaa0, iS.Jaa0_dipole + 4*ionS.J12);

% c-axis needle (Nz = 0, Nx = 1/2): no cc shape term; transverse drops by
% 4*dm_aa = 4*gfac*(4pi/Vc)*(1/2) = 1.5*lorz4.
ionN = invz_ion();  ionN.demag = 1;  ionN.alpha = 0;
[~, iN] = invz_jq_modes(ionN, [0 0 0], struct('dpRng',30,'cache',false));
verifyEqual(testCase, iN.Jshape_cc, 0, 'AbsTol', 1e-15);
verifyEqual(testCase, iN.Jcc0, i0.Jcc0);
verifyEqual(testCase, iN.Jaa0_dipole, i0.Jaa0_dipole - 1.5*lorz4, 'RelTol', 1e-9);

% cache must distinguish demag settings via the SHAPE fields (Jcc0 no longer differs)
opts = struct('dpRng',10,'cache',true);
[~, c0] = invz_jq_modes(invz_ion(), [0 0 0], opts);
ionS2 = invz_ion();  ionS2.demag = 1;  ionS2.alpha = 1;
[~, cS] = invz_jq_modes(ionS2, [0 0 0], opts);
verifyGreaterThan(testCase, abs(cS.Jshape_cc - c0.Jshape_cc), 1e-6);
verifyEqual(testCase, cS.Jcc0, c0.Jcc0);
end
```

Also, in `test_gamma_point_constants` (after the existing `info.Jcc0` assertion at line 18), add one line:

```matlab
verifyEqual(testCase, info.Jaa0*1e3, 3.512e-3*1e3, 'RelTol', 0.03);  % live transverse J(0) matches R2007
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_jq_modes.m'); assertSuccess(results)"
```

Expected: FAIL — `info` has no field `Jshape_cc`/`Jaa0`, and `iS.Jcc0 ~= i0.Jcc0` under the current draft.

- [ ] **Step 3: Implement in `invz/invz_jq_modes.m`**

Four edits (line numbers refer to the current working tree):

(a) Replace the knob comment block at lines 69–72 with:

```matlab
% Sample-shape terms: the Lorentz cavity +4pi/(3Vc) is ALWAYS added at the uniform mode (it is
% the mandatory cavity term of the dipole-sum split, not a physical toggle). The demagnetization
% correction (ion.demag/ion.alpha, default off) NEVER touches the ordering channel: per R2007 it
% cancels from the critical condition (ordering occurs at q -> 0+, not strict q = 0), so Jnu and
% info.Jcc0 are demag-invariant. The shape is exported instead as (a) info.Jshape_cc, the
% strict-uniform OBSERVABLE correction applied downstream in invz_chi_realaxis
% (chi_meas = chi/(1 + Jshape_cc*chi)), and (b) demag-aware info.Jaa0, the transverse
% (non-critical) mean-field coupling (internal-vs-applied transverse field).
```

(b) Line 104 (per-q Γ branch): change

```matlab
        Jcc = Jcc + lorz - dm_cc;                    % uniform-mode shape term: Lorentz cavity - demag
```

to

```matlab
        Jcc = Jcc + lorz;                            % uniform-mode Lorentz cavity (demag-invariant)
```

(c) Line 112 (Γ diagnostic): change

```matlab
Jcc0d = -squeeze(C.gfac*dip0(3,3,:,:)) + lorz - dm_cc;
```

to

```matlab
Jcc0d = -squeeze(C.gfac*dip0(3,3,:,:)) + lorz;
```

Leave line 113 (`Jaa0d = ... + lorz - dm_aa;`) exactly as it is.

(d) After line 118 (`info.Jcc0 = info.Jcc0_dipole + 4*ion.J12;`), add:

```matlab
info.Jaa0      = info.Jaa0_dipole + 4*ion.J12;   % transverse J(0), demag-aware (meV)
info.Jshape_cc = 4*dm_cc;                        % strict-uniform observable correction (meV); 0 when demag = 0
```

(e) Cache schema version — line 85–86: change

```matlab
pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; demag; alpha];
key = sprintf('jq_%d_%s_%s.mat', dpRng, hash_vec(qvec(:)), hash_vec(pkey));
```

to

```matlab
pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; demag; alpha; 2];   % trailing 2 = cache schema v2
key = sprintf('jq2_%d_%s_%s.mat', dpRng, hash_vec(qvec(:)), hash_vec(pkey));
```

(f) Header docstring: delete the now-false sentence at lines 49–51 (*"No demagnetizing/sample-shape term is included here … out of scope for this Gamma-point diagnostic."*) and replace with:

```matlab
%   The demagnetizing/sample-shape term is deliberately EXCLUDED from Jnu and info.Jcc0
%   (it cancels in the critical condition per R2007); it is exported separately as
%   info.Jshape_cc (strict-uniform observable correction) and inside info.Jaa0
%   (transverse channel). See the knob block below.
```

- [ ] **Step 4: Run the test to verify it passes**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_jq_modes.m'); assertSuccess(results)"
```

Expected: PASS (all 4 local tests).

- [ ] **Step 5: Commit**

```bash
git add invz/invz_jq_modes.m invz/tests/test_invz_jq_modes.m
git commit -m "fix(invz): demag-invariant ordering channel; export Jaa0 + Jshape_cc (R2007 critical-condition cancellation)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Optional `opts.Jxx0` in `invz_twolevel` / `invz_twolevel_ordered`

**Files:**
- Modify: `invz/invz_twolevel.m`, `invz/invz_twolevel_ordered.m`
- Test: `invz/tests/test_invz_single_ion.m`

- [ ] **Step 1: Write the failing test**

Append to `invz/tests/test_invz_single_ion.m`:

```matlab
function test_jxx0_override(testCase)
% opts.Jxx0 must control the transverse mean field end-to-end: with Jxx0 = 0 the
% converged hx is exactly 0, and the two-level params must inherit the override.
ion = invz_ion();
s1 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false));
verifyGreaterThan(testCase, abs(s1.hx), 1e-6);                 % default: finite MF
s0 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'Jxx0', 0));
verifyEqual(testCase, s0.hx, 0);
t1 = invz_twolevel(ion, 0.31, 4);
t0 = invz_twolevel(ion, 0.31, 4, struct('Jxx0', 0));
verifyGreaterThan(testCase, abs(t1.Delta - t0.Delta), 1e-9);   % override reached the doublet
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_single_ion.m'); assertSuccess(results)"
```

Expected: FAIL with "Too many input arguments" on the 4-arg `invz_twolevel` call (the `invz_single_ion` half already passes — it reads `opts.Jxx0` today).

- [ ] **Step 3: Implement**

`invz/invz_twolevel.m` — change the signature and single-ion call:

```matlab
function tl = invz_twolevel(ion, T, Bx, opts)
%INVZ_TWOLEVEL Electronic two-level (split doublet) parameters for the Jensen self-energy.
% opts.Jxx0 (optional): transverse MF coupling forwarded to invz_single_ion (default ion.Jxx0).
if nargin < 4, opts = struct(); end
Jxx0 = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
C  = invz_const();
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false, 'Jxx0', Jxx0));
```

(rest of the function body unchanged).

`invz/invz_twolevel_ordered.m` — same pattern with the 5th argument:

```matlab
function tl = invz_twolevel_ordered(ion, T, Bx, hz, opts)
```

after the existing header comment add:

```matlab
% opts.Jxx0 (optional): transverse MF coupling forwarded to invz_single_ion (default ion.Jxx0).
if nargin < 5, opts = struct(); end
Jxx0 = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
```

and change line 12 to:

```matlab
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false, 'hz_fixed', hz, 'Jxx0', Jxx0));
```

- [ ] **Step 4: Run test to verify it passes**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_single_ion.m'); assertSuccess(results)"
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_twolevel.m invz/invz_twolevel_ordered.m invz/tests/test_invz_single_ion.m
git commit -m "feat(invz): optional opts.Jxx0 forwarding in the two-level helpers

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Thread `Jxx0` through the solvers; observable `Jshape` correction in `invz_chi_realaxis`

**Files:**
- Modify: `invz/invz_solve_point.m`, `invz/invz_solve_point_ordered.m`, `invz/invz_chi_realaxis.m`
- Test: `invz/tests/test_invz_chi_observable.m`

- [ ] **Step 1: Write the failing test**

Append to `invz/tests/test_invz_chi_observable.m`:

```matlab
function test_demag_observable_rescale(testCase)
% The measured strict-uniform response chi_meas = chi_int/(1 + Jshape*chi_int) must
% equal the OLD-convention shifted pole chit/(1 - (Jcc0 - Jshape)*chit) exactly
% (algebraic identity), and Jshape = 0 must be a byte-identical no-op.
ion = invz_ion();
T = 0.31;  Bx = 5.0;  w = (0:0.01:0.4).';
tl0 = invz_twolevel(ion, T, Bx);
pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
Jcc0 = ion.J0eff;  Jsh = 1e-3;
oi = invz_chi_realaxis(ion, T, Bx, pt0, w, struct('Jsel', Jcc0, 'npass', 1));
om = invz_chi_realaxis(ion, T, Bx, pt0, w, struct('Jsel', Jcc0, 'npass', 1, 'Jshape', Jsh));
chit  = oi.chi0cc_w(:).' ./ (1 + oi.Sigma_w(:).');
expct = chit ./ (1 - (Jcc0 - Jsh)*chit);
verifyEqual(testCase, om.chi_cc_q(1,:), expct, 'RelTol', 1e-10);
verifyEqual(testCase, oi.Jshape, 0);
verifyEqual(testCase, om.Jshape, Jsh);
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_chi_observable.m'); assertSuccess(results)"
```

Expected: FAIL — output struct has no field `Jshape` and no correction is applied.

- [ ] **Step 3: Implement**

`invz/invz_chi_realaxis.m` — three edits:

(a) After the `Jsel` opt (line 16), add:

```matlab
Jxx0   = ion.Jxx0; if isfield(opts,'Jxx0'),   Jxx0   = opts.Jxx0;   end
Jshape = 0;        if isfield(opts,'Jshape'), Jshape = opts.Jshape; end
```

(b) Line 24 (paramagnet fallback single-ion): change to

```matlab
    si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true, 'Jxx0', Jxx0));   % paramagnet
```

(c) After the `chi_cc_q` fill loop (after line 55), add:

```matlab
if Jshape ~= 0
    % Sample-shape correction for the STRICT-UNIFORM measured observable only:
    % chi_meas = chi_int/(1 + Jshape*chi_int)  (demag-limited: the soft mode
    % saturates at 1/Jshape instead of diverging). Callers evaluating a finite-q
    % path (intrinsic probe) must NOT pass Jshape.
    out.chi_cc_q = out.chi_cc_q ./ (1 + Jshape*out.chi_cc_q);
end
out.Jshape = Jshape;
```

Also extend the header docstring (after line 12) with:

```matlab
% opts.Jxx0   (ion.Jxx0)  transverse MF coupling for the internally built single-ion state.
% opts.Jshape (0)         strict-uniform demag observable correction (use info.Jshape_cc);
%                         applied in place to chi_cc_q. Leave 0 for finite-q (intrinsic) paths.
```

`invz/invz_solve_point.m` — after line 8 (`J0eff = ...`), add:

```matlab
Jxx0  = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
```

change line 15 to:

```matlab
si  = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', hyp, 'Jxx0', Jxx0));
```

change line 18 to:

```matlab
tl  = invz_twolevel(ion, T, Bx, struct('Jxx0', Jxx0));
```

`invz/invz_solve_point_ordered.m` — after line 21 (`J0eff = ...`), add:

```matlab
Jxx0  = ion.Jxx0;  if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
```

change line 34 to:

```matlab
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', hyp, 'order', true, 'J0z', J0eff, 'Jxx0', Jxx0));
```

change line 45 to:

```matlab
tl  = invz_twolevel_ordered(ion, T, Bx, si.hz, struct('Jxx0', Jxx0));
```

- [ ] **Step 4: Run the new test AND the full fast suite** (this task touches the hot path)

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: PASS. Defaults are unchanged (`Jxx0` falls back to `ion.Jxx0`, `Jshape` to 0), so no existing test may move.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_solve_point.m invz/invz_solve_point_ordered.m invz/invz_chi_realaxis.m invz/tests/test_invz_chi_observable.m
git commit -m "feat(invz): thread opts.Jxx0 through solvers; strict-uniform demag observable correction (opts.Jshape)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Extract `invz_solve_auto`; rewire `invz_spectra_map`

**Files:**
- Create: `invz/invz_solve_auto.m`
- Modify: `invz/invz_spectra_map.m`
- Test: existing `invz/tests/test_invz_spectra_map.m` (must stay green; no new test here — Task 5's tests exercise `invz_solve_auto` again)

- [ ] **Step 1: Create `invz/invz_solve_auto.m`**

```matlab
function [pt, phase] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Ordered-first 1/z solve at one (T, Bx) point, paramagnetic fallback.
% Tries invz_solve_point_ordered; a converged spontaneous-moment solution returns with
% phase = 1 (ferromagnet). Otherwise invz_solve_point is attempted: phase = 2 on a
% converged paramagnetic solve. phase = 0 means no usable 1/z solution: pt is then the
% non-converged paramagnetic pt when one was produced (its Sigma0 may still be of
% diagnostic value), or [] when both solves errored (e.g. invz:degenerateDoublet near
% Bx -> 0 in both phases). opts passes through to both solvers (hyp, J0eff, Jxx0, ...).
%
% Option-A caveat (see invz_solve_point_ordered): the FM/PM handoff happens at the bare
% MEAN-FIELD boundary, which sits slightly above the 1/z critical field.
if nargin < 5, opts = struct(); end
pt = [];  phase = 0;
try
    pto = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts);
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
        pt = pto;  phase = 1;  return;
    end
catch
    % fall through to the paramagnetic solve
end
try
    pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
    if pt.converged && isfinite(pt.Sigma0)
        phase = 2;
    end
catch
end
end
```

- [ ] **Step 2: Rewire `invz/invz_spectra_map.m`**

(a) After line 52 (`Jcc0 = info.Jcc0;`), add (the `isfield` guards keep the synthetic-info test path working):

```matlab
Jaa0   = ion.Jxx0;  if isfield(info, 'Jaa0'),      Jaa0   = info.Jaa0;      end
Jshape = 0;         if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end
```

(b) Change the `one_field` call (lines 64–65) to:

```matlab
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k)] = ...
        one_field(ion, T, fields(k), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp);
```

(c) Replace the whole `one_field` local function (lines 79–121) with:

```matlab
function [chiz, chirpa, Sigma0, phase] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp)
%ONE_FIELD chi''_cc(omega) at one field via the shared ordered-first solve
% (invz_solve_auto); returns phase = 1 (FM), 2 (PM), or 0 (no solution -> NaN columns).
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
sopts = struct('hyp', hyp, 'J0eff', Jcc0, 'Jxx0', Jaa0);
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape);

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- ferromagnetic (ordered) branch ---
    o  = invz_chi_realaxis(ion, T, B, pt, w, copts);
    chiz = imag(o.chi_cc_q(1, :)).';
    pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                 'K', [], 'is_ordered', true, 'si', pt.si);
    c0opts = copts;  c0opts.npass = 1;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    Sigma0 = pt.Sigma0;
    return;
end

% --- paramagnetic side: bare-RPA overlay first (needs only the two-level params),
% so a non-converged 1/z point still gets its RPA column ---
try
    tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    c0opts = copts;  c0opts.npass = 1;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
catch
end
if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    o = invz_chi_realaxis(ion, T, B, pt, w, copts);
    chiz = imag(o.chi_cc_q(1, :)).';
end
end
```

(d) In the file header (around lines 10–12), extend the returned-maps description:

```matlab
%     When ion.demag ~= 0, both maps are the strict-uniform MEASURED observable
%     (demag-corrected via info.Jshape_cc; the soft mode saturates at 1/Jshape_cc
%     instead of diverging); with demag = 0 they are the intrinsic response.
```

- [ ] **Step 3: Run the full fast suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: PASS. `test_invz_spectra_map` injects a synthetic `info` without `Jaa0`/`Jshape_cc` — the fallbacks must make behavior identical to before this task.

- [ ] **Step 4: Commit**

```bash
git add invz/invz_solve_auto.m invz/invz_spectra_map.m
git commit -m "refactor(invz): extract invz_solve_auto; thread Jaa0/Jshape through the spectra map

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: `invz_spectra_qpath` + plot helper (R2007 Fig. 3 engine)

**Files:**
- Create: `invz/invz_spectra_qpath.m`, `invz/invz_plot_spectra_qpath.m`
- Test: `invz/tests/test_invz_spectra_qpath.m`

- [ ] **Step 1: Write the failing test**

Create `invz/tests/test_invz_spectra_qpath.m`:

```matlab
function tests = test_invz_spectra_qpath
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_qpath_structure_and_gamma_limit(testCase)
% Structural contract + physics anchors of the q-path spectrum (R 2007 Fig 3 engine):
%  - shapes [nw x nq]; arc-length coordinate s starts at 0;
%  - (2,0,0) IS Gamma-equivalent for the 4-site basis (structure factor 4): the max
%    coupling branch there equals the uniform-mode info.Jcc0 of the same dpRng;
%  - (1,0,0) is NOT Gamma-equivalent (structure factor 0): weaker coupling, so the
%    mode sits HIGHER in energy and the dispersion falls toward the zone centre.
ion = invz_ion();
T = 0.31;  B = 5.5;                        % paramagnetic side: fast, well-converged
w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);             % synthetic medium (as in test_invz_spectra_map)
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 1.5 0 0; 2 0 0];
S = invz_spectra_qpath(ion, T, B, qpath, w, ...
                       struct('Jnu', Jnu, 'info', info, 'dpRng', 10));

verifySize(testCase, S.chiz,   [numel(w) 3]);
verifySize(testCase, S.chirpa, [numel(w) 3]);
verifyEqual(testCase, S.s, [0 0.5 1], 'AbsTol', 1e-12);
verifyEqual(testCase, S.phase, 2);
verifyTrue(testCase, all(isfinite(S.Epeak)));

% Gamma-equivalence anchor at matching dpRng
[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
verifyEqual(testCase, S.Jq(3), iref.Jcc0, 'RelTol', 1e-9);
verifyGreaterThan(testCase, iref.Jcc0, S.Jq(1));         % (1,0,0) couples more weakly
verifyGreaterThan(testCase, S.Epeak(1), S.Epeak(3));     % mode softens toward Gamma
end

function test_qpath_demag_invariant(testCase)
% A finite-q probe is intrinsic: the path spectra must be demag-invariant when the
% transverse channel is held at the same value (compare shapes with identical Jaa0
% by keeping demag's only other entry point -- info.Jaa0 -- out via synthetic info).
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
T = 0.31;  B = 5.5;  w = (0.05:0.05:0.5).';
info = struct('Jcc0', 6.4e-3);             % same synthetic medium for both shapes
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 2 0 0];
o = struct('Jnu', Jnu, 'info', info, 'dpRng', 10);
S0 = invz_spectra_qpath(ion0, T, B, qpath, w, o);
SS = invz_spectra_qpath(ionS, T, B, qpath, w, o);
verifyEqual(testCase, SS.chiz, S0.chiz);   % bit-identical: no shape leak on the path
verifyEqual(testCase, SS.Jq,   S0.Jq);
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_spectra_qpath.m'); assertSuccess(results)"
```

Expected: FAIL with "Undefined function 'invz_spectra_qpath'".

- [ ] **Step 3: Create `invz/invz_spectra_qpath.m`**

```matlab
function S = invz_spectra_qpath(ion, T, B, qpath, w, opts)
%INVZ_SPECTRA_QPATH chi''_cc(q, omega) along a reciprocal-space path at fixed (T, B).
%   S = invz_spectra_qpath(ion, T, B, qpath, w) computes the dissipative cc susceptibility
%   along the q-path (nq x 3, r.l.u.) at fixed temperature T (K) and transverse field B (T),
%   for the 1/z-renormalized and bare-RPA theories -- cf. R 2007 Fig 3 (dispersion along
%   (h,0,0) at 0.31 K). The 1/z medium (Sigma, K, lambda) is solved ONCE at (T, B) on the
%   BZ-integration grid (ordered-first via invz_solve_auto, so both phases work); the path
%   susceptibility then follows from the same single-site response chit(w) via the pole
%   formula chi(q, w) = chit/(1 - J_nu(q) chit), evaluated at the selected coupling branch
%   along the path (opts.branch, default 4 = largest/softest; at Gamma-equivalent q the
%   max branch IS the uniform mode = info.Jcc0).
%
%   NO sample-shape (demag) correction is applied along the path: a finite-q probe
%   (neutrons) measures the INTRINSIC response, and at Gamma-equivalent path points the
%   relevant limit is q -> 0+, where the demagnetizing field cancels (R 2007). The demag
%   knob still reaches these spectra indirectly through the transverse mean field
%   (info.Jaa0 -> opts Jxx0 of the solvers).
%
%   Returns:
%     S.chiz, S.chirpa     [nw x nq]  1/z and bare-RPA chi''_cc(q, w)
%     S.Epeak, S.Epeak_rpa [1 x nq]   peak (mode) energy per q, searched at w >= peak_wmin
%     S.Jq                 [1 x nq]   selected coupling branch along the path (meV)
%     S.s                  [1 x nq]   arc-length path coordinate (r.l.u. from qpath(1,:))
%     S.qpath, S.w, S.T, S.B, S.info, S.phase (1 = FM, 2 = PM solve used)
%
%   opts fields (all optional):
%     .hyp (true), .grid ([16 16 16]), .dpRng (30), .eta (5e-3)   as in invz_spectra_map
%     .branch (4)        which of the 4 ascending-sorted coupling branches to follow
%     .peak_wmin (0.05)  meV; excludes the low-frequency hyperfine pole (R 2007 Fig 2)
%                        from the peak search so Epeak tracks the doublet mode
%     .Jnu, .info        precomputed BZ-grid branches / info (skips the lattice sum; tests)

if nargin < 6, opts = struct(); end
hyp    = getf(opts, 'hyp', true);
grid   = getf(opts, 'grid', [16 16 16]);
dpRng  = getf(opts, 'dpRng', 30);
eta    = getf(opts, 'eta', 5e-3);
branch = getf(opts, 'branch', 4);
wmin   = getf(opts, 'peak_wmin', 0.05);

w = w(:);

if isfield(opts, 'Jnu') && isfield(opts, 'info')
    Jnu = opts.Jnu(:);   info = opts.info;
else
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5]);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    [Jnu, info] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', true));
    Jnu = Jnu(:);
end
Jcc0 = info.Jcc0;
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end

% one medium solve at (T, B) -- FM below the (bare-MF) boundary, PM above
[pt, phase] = invz_solve_auto(ion, T, B, Jnu, struct('hyp', hyp, 'J0eff', Jcc0, 'Jxx0', Jaa0));
if phase == 0
    error('invz:noSolution', ...
        ['No converged 1/z solution at T = %.3f K, B = %.3f T ' ...
         '(near-degenerate doublet or critical band).'], T, B);
end

% coupling branches along the path (same convention and cache as the BZ grid)
Jpath = invz_jq_modes(ion, qpath, struct('dpRng', dpRng, 'cache', true));
Jq = Jpath(:, branch).';                          % ascending sort; 4 = largest (softest pole)

% path spectra from ONE real-axis evaluation each: Jsel is vectorized over q, and the
% single-site chit(w) is q-independent. Intrinsic: no Jshape here (see header).
copts = struct('Jsel', Jq, 'eta', eta, 'Jxx0', Jaa0);
o = invz_chi_realaxis(ion, T, B, pt, w, copts);
chiz = imag(o.chi_cc_q).';                        % [nw x nq]

if phase == 1
    pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                 'K', [], 'is_ordered', true, 'si', pt.si);
else
    tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
end
c0opts = copts;  c0opts.npass = 1;
o0 = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
chirpa = imag(o0.chi_cc_q).';

S = struct();
S.qpath = qpath;
S.s = [0 cumsum(vecnorm(diff(qpath, 1, 1), 2, 2)).'];
S.w = w;  S.T = T;  S.B = B;  S.phase = phase;  S.info = info;  S.Jq = Jq;
S.chiz = chiz;  S.chirpa = chirpa;
S.Epeak     = peak_energy(chiz,   w, wmin);
S.Epeak_rpa = peak_energy(chirpa, w, wmin);
end

% -------------------------------------------------------------------------------------------
function E = peak_energy(chi, w, wmin)
%PEAK_ENERGY per-column peak position of chi''(w) restricted to w >= wmin.
E = nan(1, size(chi, 2));
mask = w >= wmin;
wm = w(mask);
for k = 1:size(chi, 2)
    c = chi(mask, k);
    if any(isfinite(c))
        [~, i] = max(c);
        E(k) = wm(i);
    end
end
end

% -------------------------------------------------------------------------------------------
function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
```

- [ ] **Step 4: Create `invz/invz_plot_spectra_qpath.m`**

```matlab
function invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit)
%INVZ_PLOT_SPECTRA_QPATH Render one chi''(q, omega) colormap panel from invz_spectra_qpath.
%   invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit) draws `chi` ([nw x nq], e.g.
%   S.chiz or S.chirpa) against the arc-length path coordinate S.s (x) and frequency S.w
%   (y), and overlays the peak dispersion `Epeak` (pass S.Epeak or S.Epeak_rpa) in white.
%   Colour conventions match invz_plot_spectra_map: log10 scale spanning three decades
%   below the robust (99.5th-percentile) peak; NaN columns transparent on grey; present-
%   but-negative chi'' floored into the darkest colour. As there, the caller pre-scales
%   S.w and Epeak when plotting in GHz (see invz_run_spectra's eUnit knob).
if nargin < 5, ttl = ''; end
if nargin < 6, eUnit = 'meV'; end

finiteMask = isfinite(chi);
Z = log10(max(chi, realmin));
im = imagesc(ax, S.s, S.w, Z);
set(im, 'AlphaData', double(finiteMask));
set(ax, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'Layer', 'top');

pos = chi(finiteMask & chi > 0);
if ~isempty(pos)
    x = sort(pos(:));
    hi = x(max(1, min(numel(x), ceil(0.995*numel(x)))));
    clim(ax, [log10(hi/1e3) log10(hi)]);
end
colormap(ax, turbo);
hold(ax, 'on');  plot(ax, S.s, Epeak, 'w.-', 'MarkerSize', 8);  hold(ax, 'off');
xlabel(ax, sprintf('s along path from Q = [%g %g %g] (r.l.u.)', S.qpath(1,:)));
ylabel(ax, sprintf('\\omega (%s)', eUnit));
title(ax, ttl);
cb = colorbar(ax);  cb.Label.String = 'log_{10} \chi''''_{cc}';
end
```

- [ ] **Step 5: Run test to verify it passes**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_spectra_qpath.m'); assertSuccess(results)"
```

Expected: PASS (2 tests). If `test_qpath_structure_and_gamma_limit` fails on the `S.Jq(3) == iref.Jcc0` anchor, check `branch`/sort order before touching tolerances: `Jnu` columns are ascending, so column 4 must be the uniform mode at Γ-equivalent q.

- [ ] **Step 6: Commit**

```bash
git add invz/invz_spectra_qpath.m invz/invz_plot_spectra_qpath.m invz/tests/test_invz_spectra_qpath.m
git commit -m "feat(invz): q-path spectra engine + plot helper (R2007 Fig 3)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Driver knobs — `invz_run_spectra.m` (qpath view) and `invz_run_phase_diagram.m`

**Files:**
- Modify: `invz/invz_run_spectra.m` (full rewrite below), `invz/invz_run_phase_diagram.m`
- Test: manual smoke run (drivers are scripts; the underlying functions are already covered)

- [ ] **Step 1: Replace `invz/invz_run_spectra.m` with**

```matlab
%INVZ_RUN_SPECTRA chi''_cc spectra vs transverse field (uniform mode) or along a q-path.
%
% One driver, three views:
%   qpath = []  (default)  -- field sweep at the uniform mode (q = 0):
%     numel(fields) <= sliceMax  ->  1D line-slice overlay (1/z solid, RPA dashed, one
%                                    colour per field)                      cf. R 2007 Fig 2
%     numel(fields) >  sliceMax  ->  2D field-vs-frequency colormap, 1/z and RPA panels
%                                                            cf. R 2007 Fig 2 / Kovacevic Fig 3d
%   qpath = [nq x 3] r.l.u.  -- q-path view at fixed field(s) Bq:          cf. R 2007 Fig 3
%     numel(Bq) == 1  ->  2D path-vs-frequency colormaps (1/z + RPA), peak overlay
%     numel(Bq) >  1  ->  E_peak(q) dispersion overlay, one colour per field
%
% The field sweep is done by invz_spectra_map (both phases: FM below Bc, PM above; S.phase
% labels each field 1=FM/2=PM/0=masked); the q-path by invz_spectra_qpath (one 1/z medium
% solve per field, then the pole formula along the path). Both result structs replot for
% free, e.g. invz_plot_spectra_map(gca, S, S.chiz, '1/z').
% Save with:  print(gcf, 'spectra.png', '-dpng', '-r150');

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

ion = invz_ion();
% ion.demag = 1; ion.alpha = 0.5;    % OPTIONAL sample-shape (demagnetization) knob; default off.
%   demag = 0 (default): intrinsic response -- the R 2007 benchmark. demag ~= 0 with spheroid
%   aspect ratio alpha (1 sphere, 0 c-needle, Inf disk) adds the sample shape CONSISTENTLY:
%     - the strict-uniform (q = 0) observable chi''_cc is demag-corrected via info.Jshape_cc
%       (chi_meas = chi/(1 + Jshape*chi)): the soft mode saturates instead of diverging;
%     - the transverse mean field uses the shape-corrected info.Jaa0 (internal vs applied field);
%     - the ordering-channel coupling info.Jcc0/Jnu -- and with it Bc and Tc -- is demag-
%       INVARIANT (R 2007: the demagnetizing field cancels from the critical condition,
%       because ordering occurs at q -> 0+, not strict q = 0).
%   q-path spectra are intrinsic (finite-q probe) and see demag only via info.Jaa0.
T = 0.31;                             % K
% fields = [3.6 4.2 4.8 5.4 6.0];       % few -> slices;  many -> colormap
fields = linspace(3,6.5,151);
w = (0:0.002:0.45).';               % meV
eta = 1e-3;                          % real-axis line width: Lorentzian HWHM in meV (5e-3 meV ~ 1.2 GHz).
                                     % Lower -> sharper peaks (resolves the sub-6-GHz hyperfine lines),
                                     % but keep eta >~ the w-spacing above or the peaks alias.
sliceMax = 6;                         % field count at/below which the line-slice view is used
useParallel = true;                  % true -> parfor over fields (Parallel Computing Toolbox)
eUnit = 'meV';                       % 'meV' or 'GHz' -- plotting only; computation always runs in meV

% ---- q-path view (R 2007 Fig 3): set qpath non-empty to switch views --------------------
qpath = [];                          % [] = field-sweep views; [nq x 3] r.l.u. = q-path view
% qh = linspace(1, 2, 51).';  qpath = [qh zeros(numel(qh), 2)];  % (1,0,0)->(2,0,0), Fig 3 path
Bq = 4.24;                           % field(s), T, for the q-path view. One value -> colormaps;
                                     % several -> E_peak(q) overlay (R 2007 Fig 3: [3.6 4.24 6.0])
dispScale = 1;                       % dispersion display scale factor; R 2007 scales the
                                     % calculated energies by 1.15 to match experiment (Fig 3)

C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end

if ~isempty(qpath)
    % ---------------- q-path view at fixed field(s) (R 2007 Fig 3) ----------------
    if numel(Bq) == 1
        S = invz_spectra_qpath(ion, T, Bq, qpath, w, struct('eta', eta));
        Splot = S;   % display-only copy; the solve above always ran in meV
        Splot.w = S.w*eScale;  Splot.Epeak = S.Epeak*eScale;  Splot.Epeak_rpa = S.Epeak_rpa*eScale;
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);
        invz_plot_spectra_qpath(ax1, Splot, Splot.chiz, Splot.Epeak, ...
            sprintf('1/z, T = %.2f K, B = %.2f T', T, Bq), eUnit);
        ax2 = subplot(1, 2, 2);
        invz_plot_spectra_qpath(ax2, Splot, Splot.chirpa, Splot.Epeak_rpa, ...
            sprintf('RPA, T = %.2f K, B = %.2f T', T, Bq), eUnit);
    else
        figure; hold on;  co = lines(numel(Bq));
        for k = 1:numel(Bq)
            Sk = invz_spectra_qpath(ion, T, Bq(k), qpath, w, struct('eta', eta));
            plot(Sk.s, Sk.Epeak*eScale*dispScale,     '-',  'Color', co(k, :), ...
                 'DisplayName', sprintf('1/z, %.2f T', Bq(k)));
            plot(Sk.s, Sk.Epeak_rpa*eScale*dispScale, '--', 'Color', co(k, :), ...
                 'DisplayName', sprintf('RPA, %.2f T', Bq(k)));
        end
        xlabel(sprintf('s along path from Q = [%g %g %g] (r.l.u.)', qpath(1, :)));
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('T = %.2f K  (dispScale = %.2f)', T, dispScale));  legend show;
        % cf. R 2007 Fig 3: for the (1,0,0)->(2,0,0) path, Q = (1+s, 0, 0); their plotted
        % theory lines are the calculated energies scaled by 1.15 (set dispScale = 1.15).
    end
else
    % ---------------- field-sweep views at the uniform mode ----------------
    S = invz_spectra_map(ion, T, fields, w, struct('parallel', useParallel, 'eta', eta));

    if numel(fields) <= sliceMax
        figure; hold on;  co = lines(numel(fields));
        for k = 1:numel(fields)
            plot(w*eScale, S.chiz(:, k),   '-',  'Color', co(k, :), 'DisplayName', sprintf('1/z, %.2f T', fields(k)));
            plot(w*eScale, S.chirpa(:, k), '--', 'Color', co(k, :), 'DisplayName', sprintf('RPA, %.2f T', fields(k)));
        end
        xlabel(eLabel);  ylabel('\chi''''_{cc}');  title(sprintf('T = %.2f K', T));  legend show;
    else
        Splot = S;  Splot.w = S.w * eScale;    % display-only copy; solve above always ran in meV
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);  invz_plot_spectra_map(ax1, Splot, Splot.chiz,   sprintf('1/z, T = %.2f K', T), eUnit);
        ax2 = subplot(1, 2, 2);  invz_plot_spectra_map(ax2, Splot, Splot.chirpa, sprintf('RPA, T = %.2f K', T), eUnit);
    end
end
```

- [ ] **Step 2: Edit `invz/invz_run_phase_diagram.m`**

(a) Immediately after line 43 (`ion = invz_ion();`), insert the knob comment:

```matlab
% ion.demag = 1;  ion.alpha = 0.25;  % OPTIONAL sample-shape (demagnetization) knob; default off.
%   Effect here: the boundary shifts ONLY through the transverse channel -- the internal
%   transverse field uses the shape-corrected coupling info.Jaa0 (below), i.e. Bc(T) is then
%   quoted vs the APPLIED field for that sample shape. The ordering-channel coupling
%   J0 = info.Jcc0 and the branch spectrum Jnu are demag-INVARIANT by construction (R 2007:
%   the demagnetizing field cancels from the critical condition, since ordering occurs at
%   q -> 0+, not strict q = 0) -- so at B = 0 the knob does not move Tc0 at all. demag = 0
%   (default) is the intrinsic / internal-field boundary that matches the R 2007 benchmark.
```

(b) After line 48 (`J0 = info.Jcc0; ...`), add:

```matlab
Jxx0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jxx0 = info.Jaa0; end   % live transverse J(0)
% (demag-aware; at demag = 0 it differs from the hardcoded ion.Jxx0 by <0.1% -- the live
% dipole sum supersedes the pasted constant). Tc0 below needs no Jxx0: at B = 0, <Jx> = 0.
```

(c) Change the two solver calls to pass it through (`opts` flows into `invz_solve_point` in both root finders):

Line 86:

```matlab
            val = invz_critical(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'window', [0.1 6]));
```

Line 93:

```matlab
            val = invz_critical_T(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'Tc0', Tc0));
```

- [ ] **Step 3: Smoke-run the q-path view (cheap settings)**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath('invz'); addpath('.'); ion = invz_ion(); w = (0.02:0.02:0.6).'; qh = linspace(1,2,9).'; S = invz_spectra_qpath(ion, 0.31, 5.5, [qh zeros(9,2)], w, struct('grid',[8 8 8],'dpRng',10)); disp(S.Epeak); assert(all(isfinite(S.Epeak)) && S.Epeak(1) > S.Epeak(end))"
```

Expected: prints 9 finite peak energies, decreasing from (1,0,0) toward (2,0,0); assertion passes. (First run pays one 8^3 lattice sum + one 1/z solve: about a minute.)

- [ ] **Step 4: Run the full fast suite** (guards against driver-adjacent regressions)

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_run_spectra.m invz/invz_run_phase_diagram.m
git commit -m "feat(invz): demag knob in both drivers; q-path (R2007 Fig 3) view in invz_run_spectra

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Demag-invariance regression tests (the coverage gap Codex flagged)

**Files:**
- Create: `invz/tests/test_invz_demag_invariance.m`

- [ ] **Step 1: Write the tests** (they should pass immediately — they are regression locks for Tasks 1–3; if any fails, a previous task is wrong: STOP and fix there)

```matlab
function tests = test_invz_demag_invariance
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function q = small_grid()
% 6^3 grid, Gamma excluded: invariance below is exact, so grid quality is irrelevant.
cmd = "qVec_generator(invz_ion().a, 'mode', 'grid', 'grid', [6 6 6], 'range', [-0.5 0.5])";
[~, q] = evalc(cmd);
q = q(any(abs(q) > 1e-12, 2), :);
end

function test_tc0_and_sigma_crit_demag_invariant(testCase)
% R 2007: the demagnetizing field cancels from the critical condition. Sigma_c and the
% zero-field Tc read ONLY Jnu/Jcc0, so with the shape knob on (sphere) they must be
% BIT-IDENTICAL to the intrinsic values. (At B = 0 the transverse channel is inert:
% <Jx> = 0, so the demag-aware info.Jaa0 cannot enter either.)
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
q = small_grid();
[J0nu, i0] = invz_jq_modes(ion0, q, struct('dpRng', 5, 'cache', false));
[JSnu, iS] = invz_jq_modes(ionS, q, struct('dpRng', 5, 'cache', false));
verifyEqual(testCase, JSnu, J0nu);
verifyEqual(testCase, iS.Jcc0, i0.Jcc0);
Sc0 = invz_sigma_crit(i0.Jcc0, J0nu(:));
ScS = invz_sigma_crit(iS.Jcc0, JSnu(:));
verifyEqual(testCase, ScS, Sc0);
verifyEqual(testCase, invz_critical_T0field(ionS, ScS, iS.Jcc0), ...
                      invz_critical_T0field(ion0, Sc0, i0.Jcc0));
end

function test_crit_demag_enters_only_via_transverse_channel(testCase)
% With the transverse coupling PINNED to the intrinsic value (opts.Jxx0), the full 1/z
% criticality at finite field must be demag-invariant -- the ordering channel carries no
% shape dependence. Unpinned (Jxx0 = demag-aware info.Jaa0), crit legitimately moves:
% that is the internal-vs-applied transverse field correction. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
q = small_grid();
[J0nu, i0] = invz_jq_modes(ion0, q, struct('dpRng', 5, 'cache', false));
[JSnu, iS] = invz_jq_modes(ionS, q, struct('dpRng', 5, 'cache', false));
T = 0.31;  B = 5.0;
p0 = invz_solve_point(ion0, T, B, J0nu(:), struct('J0eff', i0.Jcc0, 'Jxx0', ion0.Jxx0));
pS = invz_solve_point(ionS, T, B, JSnu(:), struct('J0eff', iS.Jcc0, 'Jxx0', ion0.Jxx0));
verifyEqual(testCase, pS.crit, p0.crit);          % pinned: bit-identical
pU = invz_solve_point(ionS, T, B, JSnu(:), struct('J0eff', iS.Jcc0, 'Jxx0', iS.Jaa0));
verifyGreaterThan(testCase, abs(pU.crit - p0.crit), 1e-9);   % unpinned: physical shift
end
```

- [ ] **Step 2: Run the fast test**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_demag_invariance.m'); assertSuccess(results)"
```

Expected: PASS (1 passed, 1 filtered as SLOW).

- [ ] **Step 3: Run the SLOW test once**

```bash
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_demag_invariance.m'); assertSuccess(results)"
```

Expected: PASS (2 passed). Runtime: a few minutes (two full 1/z solves at 0.31 K on a 6^3 grid plus one repeat).

- [ ] **Step 4: Commit**

```bash
git add invz/tests/test_invz_demag_invariance.m
git commit -m "test(invz): lock phase-boundary demag invariance (ordering channel) + transverse-only entry

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Documentation sync + final verification

**Files:**
- Modify: `invz/invz_ion.m`, `invz/README.html`

- [ ] **Step 1: Fix the `invz_ion.m` knob comment**

Replace lines 14–17 of `invz/invz_ion.m` (the block starting `% Sample-shape (demagnetization) correction ...`) with:

```matlab
% Sample-shape (demagnetization) knob. demag = 0 (default): intrinsic couplings, the R2007
% benchmark. demag ~= 0: the ellipsoid shape (ellipsoid_demagn(alpha)) enters ONLY as
%   (a) info.Jshape_cc -- strict-uniform observable correction applied in invz_chi_realaxis
%       (chi_meas = chi/(1 + Jshape_cc*chi)), and
%   (b) demag-aware info.Jaa0 -- the transverse mean-field channel (internal vs applied field).
% The ordering channel (info.Jcc0, Jnu) and hence Bc/Tc(ordering) are demag-INVARIANT: per
% R2007 the demagnetizing field cancels from the critical condition (ordering at q -> 0+).
```

(keep lines 18–19, the `ion.demag = 0;` / `ion.alpha = 1;` defaults, unchanged).

- [ ] **Step 2: Update `invz/README.html`**

Locate the demag callout near line 130–136 (the `Shape scope` note and the displayed `J(0) -> J(0) + Lorentz - demag` equation). Replace that whole demag block so it states:

1. The corrected equation set: `Jnu`/`info.Jcc0` carry the Lorentz cavity only; `info.Jshape_cc = 4*(4π/Vc)*gfac*demag*N_zz(α)` corrects the strict-uniform observable in `invz_chi_realaxis` (`chi_meas = chi/(1+Jshape_cc*chi)` — algebraically identical to the old shifted pole `chit/(1-(Jcc0-Jshape_cc)*chit)`); `info.Jaa0 = Jaa0_dipole(demag-aware) + 4*J12` feeds the transverse mean field via `opts.Jxx0`.
2. Phase-boundary behavior: `Tc(B=0)` demag-invariant exactly; `Bc(T)` shifts only via the transverse channel (applied-vs-internal field), never via the ordering coupling — "per Rønnow et al. 2007 (who note this contradicts a suggestion of Chakraborty et al.)".
3. Remove the stale "wire info.Jaa0 into the transverse MF" TODO (it is now wired) and the stale reference to an `info.Jaa0` field that "doesn't exist" — it exists now.
4. Add `invz_spectra_qpath` / `invz_plot_spectra_qpath` / `invz_solve_auto` rows to the function table, and a short q-path quick-start (Fig. 3 path snippet from the driver, note `dispScale = 1.15` comparison and that path spectra are intrinsic/no demag).
5. Note the cache-format bump: old `invz/cache/jq_*.mat` files are orphaned (prefix now `jq2_`) and can be deleted.

This is prose editing of an HTML file — match the existing `<div class="callout note">`/table markup style already used in the file.

- [ ] **Step 3: Full verification — fast suite, then slow suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: fast suite all green (now ~37+ passed). Slow suite all green (~6–10 min; `test_critical_field_at_310mK` re-derives Bc with the live `Jxx0` — the published-range bounds [4.0, 4.6] T absorb the <0.1% transverse-coupling drift).

- [ ] **Step 4: Commit**

```bash
git add invz/invz_ion.m invz/README.html
git commit -m "docs(invz): demag semantics (observable + transverse channels), q-path quick-start, cache-bump note

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

- [ ] **Step 5: Report**

Summarize for the user: what changed physically (demag no longer moves the ordering-channel criticality; observable + transverse channels carry the shape), the two driver knobs, the Fig. 3 workflow (`qpath` + `Bq = [3.6 4.24 6.0]` + `dispScale = 1.15`), the expected small `Jxx0` drift, and the cache invalidation.
