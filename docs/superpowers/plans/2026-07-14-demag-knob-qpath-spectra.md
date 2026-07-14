# Demag Knob Rework + Exploratory q-Path Dispersion (R2007 Fig. 3 Trends) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Revision 2 (2026-07-14):** incorporates the Codex plan review (`q-spectrum_plan_review_by_Codex.md`). Key changes vs rev 1: Task 0 no longer blanket-commits the working tree (baseline test first, user-approval gate for a scoped checkpoint); a direction-aware Γ-limit guard (`invz_jq_path`) replaces raw truncated-sum path couplings, with the reviewer's h→2 cutoff regression; the q-path view gets its own frequency grid (0.85 meV) and a censoring/interpolating peak picker; `invz_solve_auto` catches only `invz:*` errors and rethrows the rest; `hyp` is threaded through `invz_chi_realaxis`; the Task 2 test gets its missing `end` and ordered-path coverage; all demag wording follows one four-bullet formulation; the feature is consistently labeled *exploratory branch-resolved dispersion*, not a quantitative Fig. 3 reproduction.

**Goal:** Make the sample-shape (demagnetization) correction physically correct and user-toggleable in both run scripts, and add an **exploratory branch-resolved q-path energy/susceptibility workflow** for comparison with the **trends** in Rønnow et al. PRB 75, 054426 (2007) Fig. 3.

**Architecture:** The ordering-channel coupling (`info.Jcc0`, `Jnu`) becomes demag-invariant (per R2007 the demag field cancels from the critical condition); the shape enters only through (a) a new strict-uniform observable correction `info.Jshape_cc` applied inside `invz_chi_realaxis`, and (b) a new live transverse coupling `info.Jaa0` (demag-aware) threaded as `opts.Jxx0` through every single-ion solve. A new `invz_jq_path` evaluates coupling branches along a path with a direction-aware Γ-limit guard against the truncated dipole sum; `invz_spectra_qpath` reuses the once-solved 1/z medium plus the vectorized `Jsel` machinery of `invz_chi_realaxis` to evaluate chi''(q,w) along the path.

**Tech Stack:** MATLAB (R2025a), matlab.unittest `functiontests` suite in `invz/tests/`, no toolboxes required (Parallel Computing Toolbox optional).

---

## Context and physics rationale (read before executing)

- **Working directory / repo root:** `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion`. All paths below are relative to it. Current branch: `invz-1z-lihof4`. MATLAB invocation (from repo root):
  `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"`
  Single file: replace `'invz/tests'` with e.g. `'invz/tests/test_invz_jq_modes.m'`. Slow tests additionally need env `INVZ_SLOW=1`.
- **The bug being fixed** (Code_review_byCodex.md Finding #1, High): the current working tree folds the demag term into `info.Jcc0`, which feeds the mean field, the RPA denominator, AND the 1/z critical condition. R2007 (p. 054426-2/3) states verbatim: *"the demagnetization field cancels out in the determination of the critical condition, no matter the shape of the individual domains in the ordered phase … the ordering is going to occur at wave vectors infinitesimally different from zero, not at q ≡ 0."* Note R2007 says this **against** a contrary suggestion by Chakraborty et al. (their ref. 17) — cite it as "per R2007" not "per settled theory".

### The canonical demag wording (use this everywhere — comments, README, docstrings)

1. `info.Jcc0`, `Jnu`, and the ordering-channel contribution to criticality are **demag-invariant**.
2. `Tc(B=0)` is **exactly demag-invariant** because the transverse moment vanishes there.
3. `Bc(T)` expressed versus **applied** transverse field **can shift** through the demag-aware `info.Jaa0` (internal-vs-applied transverse field relation).
4. q-path calculations omit the longitudinal strict-uniform `Jshape_cc` transform (a finite-q probe measures the intrinsic longitudinal response), but **can still change** through the transverse applied/internal-field relation (`info.Jaa0`).

Never write "Bc and Tc are demag-invariant" or "q-path spectra have no demag" — both are wrong under this architecture.

- **Where demag legitimately acts** (both per R2007's own usage):
  1. The strict-uniform **measured observable**: `chi_meas = chi_int / (1 + Jshape*chi_int)` with `Jshape = 4*dm_cc` (the factor 4 is the uniform-mode projection of a scalar-broadcast term onto the 4-site sublattice vector — same reason `info.Jcc0 = Jcc0_dipole + 4*J12`). Algebraically this reproduces the old shifted-pole result `chit/(1-(Jcc0-Jshape)*chit)` exactly, so the observable is unchanged relative to the old convention — only criticality/Sigma/mean-field change.
  2. The **transverse (aa) channel**: `hx = Jxx0*<Jx>` is not a critical/order-parameter channel (no q→0+ cancellation applies); a shape term there is the internal-vs-applied transverse field correction. Decision (user-approved): `info.Jaa0` is demag-aware directly.
- **Lorentz term:** stays unconditionally on, no knob (user-approved). It is the mandatory cavity term of the dipole-sum split, not a physical toggle.
- **Truncated dipole sum near Γ-equivalent points (review Finding 1):** `MF_dipole`'s sharp real-space cutoff cannot resolve `|q − G| ≲ 1/(dpRng·a)`. Approaching `(2,0,0)` the computed max branch collapses (measured: 0.00156 meV at h=1.999, dpRng 30, vs the correct ≈0.00642) and then jumps at the exact endpoint where `is_gamma_equiv` adds the Lorentz term. **Guard convention adopted:** within a trust radius `ksnap = 2.5·2π/(dpRng·min‖a_i‖)` of a Γ-equivalent point G, replace the truncated branches by the exact **directional long-wavelength limit** — eigenvalues of `J_reg(Γ) + gfac·(4π/Vc)·(1/3 − k̂_z²)` (scalar sublattice broadcast), where `k̂` is the Cartesian direction of `q − G` (Cohen–Keffer nonanalytic term). For any in-plane approach (`k̂_z = 0`, the Fig. 3 path) this limit coincides with the uniform-mode Lorentz value, i.e. with `info.Jcc0`'s convention — so the endpoint is continuous by construction. Points at exactly G use the local path direction. This guard lives ONLY in the new `invz_jq_path` (path evaluation); the BZ-grid path through `invz_jq_modes` is untouched (grid-quadrature quality is Codex #2, out of scope). Near **non**-Γ-equivalent integer points (e.g. `(1,0,0)`) other (staggered) branches keep truncated-sum quality — disclosed as exploratory, not guarded.
- **Exploratory status (review Finding 3) — say this in docstrings, plot titles, README:** the q-path output is a *branch-resolved susceptibility*: one sorted coupling eigenvalue per q fed to the scalar pole formula. It is NOT neutron scattering intensity (no eigenvector/sublattice-interference weights, no polarization factor, no magnetic form factor). Branch index = sorted-eigenvalue position per q; branch identity is not tracked through crossings. The engine also inherits three known unresolved issues: biased closed-grid BZ quadrature (Codex #2), possible negative spectral weight in the real-axis continuation (Codex #3), and the FM/PM handoff at the bare mean-field boundary rather than the 1/z boundary (Codex #4, retained by `invz_solve_auto`). Energy-trend comparison with R2007 Fig. 3 is reasonable (Fig. 3 is itself an energy comparison); quantitative reproduction is NOT claimed.
- **Expected small numeric drift:** drivers switch the transverse coupling from the hardcoded `ion.Jxx0 = 3.512e-3` meV to the live `info.Jaa0 ≈ 3.5104e-3` meV (dpRng=30). ~0.05% shift in the transverse channel → sub-mK/mT boundary drift even at demag=0. Intended (self-consistency); slow-test tolerances absorb it. `ion.Jxx0` remains the fallback default inside `invz_single_ion` so direct callers are unaffected.
- **Cache invalidation:** `invz_jq_modes` results change schema (new `info.Jaa0`, `info.Jshape_cc`; demag removed from `Jcc0`). The cache key gets a schema version (filename prefix `jq2_`), so all old `invz/cache/jq_*.mat` files are simply never matched again (safe to leave or delete).
- **Error handling in new code (review Finding 8):** catch ONLY `invz:*`-prefixed identifiers (expected numerical conditions like `invz:degenerateDoublet`); record them in a diagnostic; rethrow everything else. Pre-existing broad catches elsewhere are Codex #7, out of scope.
- **Out of scope (do not attempt):** Ewald summation (the Γ-limit guard above is the sanctioned deferral); off-diagonal dipolar (ODD/full-tensor) terms; an aa-channel observable; BZ-grid endpoint bias (Codex #2); causality/positivity of the real-axis continuation (Codex #3); the FM/PM branch-selection boundary (Codex #4); eigenvector-overlap branch tracking.
- **Dirty working tree (review Finding 4):** the tree contains ~32 changed/untracked entries including work unrelated to this plan (PDFs, Data/, verification scripts, doc edits). NEVER `git add -A`. Every task commit lists explicit paths of files that task touched. Task 0 runs the baseline test BEFORE any commit and gates any checkpoint on explicit user approval of an exact path list. Commit-message trailers (`Co-Authored-By`) follow the session tooling convention; drop them if the user says so.

### File map

| File | Action | Role |
|---|---|---|
| `invz/invz_jq_modes.m` | modify | single source of J(q): demag-free `Jcc0`/`Jnu`, new `info.Jaa0`, `info.Jshape_cc`, cache v2 |
| `invz/invz_twolevel.m`, `invz/invz_twolevel_ordered.m` | modify | optional `opts.Jxx0` forwarding |
| `invz/invz_solve_point.m`, `invz/invz_solve_point_ordered.m` | modify | forward `opts.Jxx0` |
| `invz/invz_chi_realaxis.m` | modify | forward `opts.Jxx0` and `opts.hyp`; apply `opts.Jshape` observable correction |
| `invz/invz_solve_auto.m` | create | shared ordered-first/para-fallback solve; `invz:*`-only catches + diagnostics |
| `invz/invz_spectra_map.m` | modify | thread `Jaa0`/`Jshape`/`hyp`; `one_field` uses `invz_solve_auto` |
| `invz/invz_jq_path.m` | create | path couplings with direction-aware Γ-limit guard; dual path coordinates |
| `invz/invz_spectra_qpath.m` | create | chi''(q,w) along a path at fixed (T,B); censoring peak picker |
| `invz/invz_plot_spectra_qpath.m` | create | colormap panel for the q-path view |
| `invz/invz_run_spectra.m` | modify | corrected demag knob comment; qpath/Bq/wq/dispScale knobs + third view |
| `invz/invz_run_phase_diagram.m` | modify | demag knob comment; hoist `Jxx0 = info.Jaa0` into solver opts |
| `invz/invz_ion.m` | modify | demag comment corrected (canonical wording) |
| `invz/tests/test_invz_jq_modes.m` | modify | rewrite demag test (invariance), `Jaa0` assertions |
| `invz/tests/test_invz_single_ion.m` | modify | add `Jxx0` override test (para + ordered paths) |
| `invz/tests/test_invz_chi_observable.m` | modify | add observable-rescale identity test |
| `invz/tests/test_invz_solve_auto.m` | create | expected-error swallowing vs unexpected-error rethrow |
| `invz/tests/test_invz_demag_invariance.m` | create | Tc0/Sigma_c invariance; pinned-Jxx0 crit invariance (SLOW) |
| `invz/tests/test_invz_spectra_qpath.m` | create | structural + Γ-limit anchors; Γ-approach cutoff regression; peak censoring |
| `invz/README.html` | modify | demag semantics (canonical wording), new functions, qpath quick-start, exploratory caveats |
| `.gitignore` | create/modify | ignore `*.asv` |

---

### Task 0: Baseline verification and scoped checkpoint (user-gated)

- [ ] **Step 1: Run the fast test suite BEFORE touching anything**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: all pass (~32 passed, ~8 filtered as SLOW). If the baseline is red, STOP and report — do not proceed and do not commit anything.

- [ ] **Step 2: Ignore MATLAB autosaves and commit only that**

Create (or append to) `.gitignore` at repo root:

```gitignore
*.asv
```

```bash
git add .gitignore
git commit -m "chore: ignore MATLAB autosave (*.asv) files

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

- [ ] **Step 3: USER APPROVAL GATE — scoped checkpoint of pre-existing edits to plan-touched files**

The following files that this plan modifies already carry uncommitted changes (dirty) or are untracked, so a task commit would otherwise silently bundle pre-existing work:

- dirty: `invz/invz_jq_modes.m`, `invz/invz_ion.m`, `invz/invz_run_spectra.m`, `invz/invz_run_phase_diagram.m`, `invz/invz_chi_realaxis.m`, `invz/tests/test_invz_jq_modes.m`, `invz/tests/test_invz_single_ion.m`, `invz/README.html`
- untracked (plan dependencies): `invz/invz_solve_point_ordered.m`, `invz/invz_twolevel_ordered.m`, `invz/invz_sigma_ordered.m`, `invz/invz_spectra_map.m`, `invz/invz_plot_spectra_map.m`, `invz/tests/test_invz_spectra_map.m`, `invz/tests/test_invz_ordered_phase.m`, `invz/tests/test_invz_sigma_ordered.m`

ASK THE USER (AskUserQuestion or equivalent pause): "May I commit exactly these N files as a pre-plan checkpoint so per-task diffs stay clean? Files outside this list (PDFs, Data/, verify scripts, unrelated docs) stay uncommitted." Show `git diff --stat -- <list>` first.

- If approved: `git add <exact list above>` (no `-A`), review `git status`, commit as
  `wip(invz): checkpoint pre-plan state of files this plan modifies (demag draft, ordered-phase, eUnit driver)`.
- If declined: proceed on the dirty tree; each task commit still uses explicit paths, and its message must note when it subsumes pre-existing uncommitted hunks in those files.

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

Six edits (line numbers refer to the current working tree):

(a) Replace the knob comment block at lines 69–72 with:

```matlab
% Sample-shape terms: the Lorentz cavity +4pi/(3Vc) is ALWAYS added at the uniform mode (it is
% the mandatory cavity term of the dipole-sum split, not a physical toggle). The demagnetization
% correction (ion.demag/ion.alpha, default off) NEVER touches the ordering channel: per R2007 it
% cancels from the critical condition (ordering occurs at q -> 0+, not strict q = 0), so Jnu,
% info.Jcc0, and Tc(B=0) are demag-invariant. The shape is exported instead as (a) info.Jshape_cc,
% the strict-uniform OBSERVABLE correction applied downstream in invz_chi_realaxis
% (chi_meas = chi/(1 + Jshape_cc*chi)), and (b) demag-aware info.Jaa0, the transverse
% (non-critical) mean-field coupling -- through which Bc(T) vs APPLIED field can still shift
% (internal-vs-applied transverse field relation).
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

(e) Cache schema version — lines 85–86: change

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

- [ ] **Step 5: Commit (explicit paths)**

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

Append to `invz/tests/test_invz_single_ion.m` (note: the file terminates every local function with `end` — this one included):

```matlab
function test_jxx0_override(testCase)
% opts.Jxx0 must control the transverse mean field end-to-end, on BOTH the
% paramagnetic and ordered paths: with Jxx0 = 0 the converged hx is exactly 0,
% and the two-level helpers must inherit the override.
ion = invz_ion();
% paramagnetic path
s1 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false));
verifyGreaterThan(testCase, abs(s1.hx), 1e-6);                 % default: finite MF
s0 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'Jxx0', 0));
verifyEqual(testCase, s0.hx, 0);
t1 = invz_twolevel(ion, 0.31, 4);
t0 = invz_twolevel(ion, 0.31, 4, struct('Jxx0', 0));
verifyGreaterThan(testCase, abs(t1.Delta - t0.Delta), 1e-9);   % override reached the doublet
% ordered path (params match the existing ordered-branch test at (0.5 K, 2 T))
so  = invz_single_ion(ion, 0.5, [2 0 0], struct('hyp', false, 'order', true));
so0 = invz_single_ion(ion, 0.5, [2 0 0], struct('hyp', false, 'order', true, 'Jxx0', 0));
verifyGreaterThan(testCase, abs(so.hx), 1e-6);
verifyEqual(testCase, so0.hx, 0);
to1 = invz_twolevel_ordered(ion, 0.5, 2, so.hz);
to0 = invz_twolevel_ordered(ion, 0.5, 2, so.hz, struct('Jxx0', 0));
verifyGreaterThan(testCase, abs(to1.Delta - to0.Delta), 1e-9);
end
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

- [ ] **Step 5: Commit (explicit paths)**

```bash
git add invz/invz_twolevel.m invz/invz_twolevel_ordered.m invz/tests/test_invz_single_ion.m
git commit -m "feat(invz): optional opts.Jxx0 forwarding in the two-level helpers

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Thread `Jxx0`/`hyp` through the solvers; observable `Jshape` correction in `invz_chi_realaxis`

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
hyp    = true;     if isfield(opts,'hyp'),    hyp    = opts.hyp;    end
```

(b) Line 24 (paramagnet fallback single-ion): change to

```matlab
    si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', hyp, 'Jxx0', Jxx0));   % paramagnet
```

(c) After the `chi_cc_q` fill loop (after line 55), add:

```matlab
if Jshape ~= 0
    % Sample-shape correction for the STRICT-UNIFORM measured observable only:
    % chi_meas = chi_int/(1 + Jshape*chi_int)  (demag-limited: the soft mode
    % saturates at 1/Jshape instead of diverging). Callers evaluating a finite-q
    % path (intrinsic longitudinal probe) must NOT pass Jshape.
    out.chi_cc_q = out.chi_cc_q ./ (1 + Jshape*out.chi_cc_q);
end
out.Jshape = Jshape;
```

Also extend the header docstring (after line 12) with:

```matlab
% opts.Jxx0   (ion.Jxx0)  transverse MF coupling for the internally built single-ion state.
% opts.hyp    (true)      hyperfine manifold for that internal state; pass the caller's hyp so
%                         the real-axis chi0 matches the Matsubara medium's Hilbert space.
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

Expected: PASS. Defaults are unchanged (`Jxx0` falls back to `ion.Jxx0`, `Jshape` to 0, `hyp` to true), so no existing test may move.

- [ ] **Step 5: Commit (explicit paths)**

```bash
git add invz/invz_solve_point.m invz/invz_solve_point_ordered.m invz/invz_chi_realaxis.m invz/tests/test_invz_chi_observable.m
git commit -m "feat(invz): thread opts.Jxx0/hyp through solvers; strict-uniform demag observable correction (opts.Jshape)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: `invz_solve_auto` (typed error handling) + rewire `invz_spectra_map`

**Files:**
- Create: `invz/invz_solve_auto.m`
- Modify: `invz/invz_spectra_map.m`
- Test: `invz/tests/test_invz_solve_auto.m` (new); existing `invz/tests/test_invz_spectra_map.m` must stay green

- [ ] **Step 1: Write the failing test**

Create `invz/tests/test_invz_solve_auto.m`:

```matlab
function tests = test_invz_solve_auto
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_expected_invz_errors_become_phase0(testCase)
% Known numerical conditions (invz:* identifiers) are absorbed as "no solution"
% with the identifier preserved in the diagnostic. T = 1.9 K > Tc0, Bx = 0.05 T:
% the ordered solve relaxes to no moment (no error), then the paramagnetic
% two-level helper raises invz:degenerateDoublet.
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[pt, phase, di] = invz_solve_auto(ion, 1.9, 0.05, Jnu, struct('J0eff', 6.4e-3));
verifyEqual(testCase, phase, 0);
verifyEmpty(testCase, pt);
verifyEqual(testCase, di.para_err, 'invz:degenerateDoublet');
end

function test_unexpected_errors_rethrow(testCase)
% Programming/API defects must NOT be converted into phase = 0: a malformed ion
% struct (missing field) raises a non-invz error that must propagate.
ion = rmfield(invz_ion(), 'J');
Jnu = linspace(-2e-3, 6.0e-3, 24).';
verifyError(testCase, ...
    @() invz_solve_auto(ion, 0.31, 5.5, Jnu, struct('J0eff', 6.4e-3)), ...
    'MATLAB:nonExistentField');
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_solve_auto.m'); assertSuccess(results)"
```

Expected: FAIL with "Undefined function 'invz_solve_auto'".

- [ ] **Step 3: Create `invz/invz_solve_auto.m`**

```matlab
function [pt, phase, diag] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Ordered-first 1/z solve at one (T, Bx) point, paramagnetic fallback.
% Tries invz_solve_point_ordered; a converged spontaneous-moment solution returns with
% phase = 1 (ferromagnet). Otherwise invz_solve_point is attempted: phase = 2 on a
% converged paramagnetic solve. phase = 0 means no usable 1/z solution: pt is then the
% non-converged paramagnetic pt when one was produced (its Sigma0 may still be of
% diagnostic value), or [] when a solver raised an expected numerical condition.
% opts passes through to both solvers (hyp, J0eff, Jxx0, ...).
%
% Error policy: ONLY invz:* identifiers (expected numerical conditions, e.g.
% invz:degenerateDoublet near Bx -> 0) are absorbed; their identifiers are returned in
% diag.ordered_err / diag.para_err. Any other exception is a programming/API defect and
% is RETHROWN rather than converted into "no solution".
%
% Option-A caveat (see invz_solve_point_ordered): the FM/PM handoff happens at the bare
% MEAN-FIELD boundary, which sits slightly above the 1/z critical field.
if nargin < 5, opts = struct(); end
pt = [];  phase = 0;  diag = struct('ordered_err', '', 'para_err', '');
try
    pto = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts);
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
        pt = pto;  phase = 1;  return;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.ordered_err = err.identifier;
end
try
    pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
    if pt.converged && isfinite(pt.Sigma0)
        phase = 2;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.para_err = err.identifier;
    pt = [];
end
end
```

- [ ] **Step 4: Rewire `invz/invz_spectra_map.m`**

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
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp);

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

% --- paramagnetic side: bare-RPA overlay first (needs only the two-level params), so a
% non-converged 1/z point still gets its RPA column. invz:degenerateDoublet is the one
% expected condition here (Bx -> 0); anything else is a defect and propagates.
try
    tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    c0opts = copts;  c0opts.npass = 1;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
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

- [ ] **Step 5: Run the new test and the full fast suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: PASS. `test_invz_spectra_map` injects a synthetic `info` without `Jaa0`/`Jshape_cc` — the fallbacks must make behavior identical to before this task.

- [ ] **Step 6: Commit (explicit paths)**

```bash
git add invz/invz_solve_auto.m invz/invz_spectra_map.m invz/tests/test_invz_solve_auto.m
git commit -m "refactor(invz): extract invz_solve_auto (invz:*-only catches, diagnostics); thread Jaa0/Jshape/hyp through the spectra map

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: `invz_jq_path` (Γ-limit guard) + `invz_spectra_qpath` + plot helper

**Files:**
- Create: `invz/invz_jq_path.m`, `invz/invz_spectra_qpath.m`, `invz/invz_plot_spectra_qpath.m`
- Test: `invz/tests/test_invz_spectra_qpath.m`

- [ ] **Step 1: Write the failing tests**

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

function test_gamma_approach_regression(testCase)
% Review finding (blocker): the sharply truncated MF_dipole sum collapses on the
% approach to the Gamma-equivalent (2,0,0) (max branch fell to ~0.0016 meV at
% h = 1.999, dpRng 30, vs the correct ~0.0064) and then jumped at the endpoint.
% The direction-aware snap must give a smooth, cutoff-stable approach: for this
% in-plane path (khat_z = 0) the directional limit is the uniform-mode Lorentz
% value, so the endpoint equals info.Jcc0 by construction.
ion = invz_ion();
hs = [1.90 1.96 1.98 1.99 1.999 2.0].';
qpath = [hs zeros(numel(hs), 2)];
P20 = invz_jq_path(ion, qpath, struct('dpRng', 20, 'cache', false));
P40 = invz_jq_path(ion, qpath, struct('dpRng', 40, 'cache', false));
J20 = P20.Jnu(:, 4);  J40 = P40.Jnu(:, 4);          % max branch (ascending sort)
[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 40, 'cache', false));
verifyEqual(testCase, J40(end), iref.Jcc0, 'RelTol', 1e-9);   % endpoint = uniform mode
% no dip/jump: every point within 3% of the endpoint on this fine approach
verifyGreaterThan(testCase, min(J40)/J40(end), 0.97);
verifyGreaterThan(testCase, min(J20)/J20(end), 0.97);
% cutoff stability: dpRng 20 vs 40 agree pointwise within 2%
verifyEqual(testCase, J20, J40, 'RelTol', 0.02);
% the guard actually engaged near the endpoint, and reports it
verifyTrue(testCase, any(P40.snapped) && P40.snapped(end));
% dual path coordinates: index (r.l.u.) and Cartesian (Ang^-1)
verifyEqual(testCase, P40.s(end), 0.10, 'AbsTol', 1e-12);
verifyEqual(testCase, P40.s_cart(end), 2*pi*0.10/ion.a(1,1), 'RelTol', 1e-9);
end

function test_qpath_structure_and_gamma_limit(testCase)
% Structural contract + physics anchors of the q-path spectrum:
%  - shapes [nw x nq]; index path coordinate starts at 0;
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

[~, iref] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
verifyEqual(testCase, S.Jq(3), iref.Jcc0, 'RelTol', 1e-9);
verifyGreaterThan(testCase, iref.Jcc0, S.Jq(1));         % (1,0,0) couples more weakly
verifyGreaterThan(testCase, S.Epeak(1), S.Epeak(3));     % mode softens toward Gamma
end

function test_qpath_peak_censoring(testCase)
% Review finding (blocker): a maximum in the first/last usable frequency bin means
% the true peak lies outside the window and must be censored (NaN), not reported.
ion = invz_ion();
T = 0.31;  B = 5.5;
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 2 0 0];
o = struct('Jnu', Jnu, 'info', info, 'dpRng', 10);
% (a) clipped window: the mode at 5.5 T sits above 0.15 meV at every path point
wclip = (0.02:0.02:0.15).';
Sc = invz_spectra_qpath(ion, T, B, qpath, wclip, o);
verifyTrue(testCase, all(isnan(Sc.Epeak)));
% (b) peak_wmin at/above the window top: nothing usable, all censored, no error
w = (0.02:0.02:0.6).';
o2 = o;  o2.peak_wmin = 1.0;
Sm = invz_spectra_qpath(ion, T, B, qpath, w, o2);
verifyTrue(testCase, all(isnan(Sm.Epeak)));
end

function test_qpath_no_ordering_channel_shape_leak(testCase)
% With the transverse channel held identical (synthetic info carries no Jaa0, so
% both runs fall back to ion.Jxx0), the path spectra must be bit-identical for
% demag on vs off: the guard and the couplings carry no ordering-channel shape
% dependence. (The REAL demag pathway into q-path spectra is info.Jaa0 only.)
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
T = 0.31;  B = 5.5;  w = (0.05:0.05:0.5).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
qpath = [1 0 0; 2 0 0];
o = struct('Jnu', Jnu, 'info', info, 'dpRng', 10);
S0 = invz_spectra_qpath(ion0, T, B, qpath, w, o);
SS = invz_spectra_qpath(ionS, T, B, qpath, w, o);
verifyEqual(testCase, SS.chiz, S0.chiz);
verifyEqual(testCase, SS.Jq,   S0.Jq);
end
```

- [ ] **Step 2: Run test to verify it fails**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_spectra_qpath.m'); assertSuccess(results)"
```

Expected: FAIL with "Undefined function 'invz_jq_path'" / "'invz_spectra_qpath'".

- [ ] **Step 3: Create `invz/invz_jq_path.m`**

```matlab
function P = invz_jq_path(ion, qpath, opts)
%INVZ_JQ_PATH Coupling branches along a q-path with a direction-aware Gamma-limit guard.
%   P = invz_jq_path(ion, qpath) evaluates the 4 cc coupling branches (invz_jq_modes
%   convention, ascending sort) at each row of qpath (r.l.u.), then REPLACES values inside
%   a trust radius of any Gamma-equivalent point G by the exact directional long-wavelength
%   limit. Rationale: MF_dipole's sharp real-space cutoff cannot resolve |q - G| below
%   ~1/(dpRng*a) -- the raw branch collapses on approach to G and then jumps at the exact
%   endpoint where the Lorentz term is added (verified: max branch 0.0016 meV at h = 1.999
%   vs the correct 0.0064 at dpRng 30). The directional limit is
%       eig( J_reg(Gamma) + gfac*(4*pi/Vc)*(1/3 - khat_z^2) )   [scalar sublattice broadcast]
%   with khat the Cartesian direction of q - G (Cohen-Keffer nonanalytic term). For any
%   in-plane approach (khat_z = 0, e.g. the (h,0,0) R2007 Fig-3 path) this equals the
%   uniform-mode Lorentz value, i.e. the info.Jcc0 convention -- continuous endpoint by
%   construction. Points at exactly G use the LOCAL PATH DIRECTION; a single-point path
%   at G defaults to the in-plane (uniform-mode) convention.
%
%   Scope: the guard covers Gamma-EQUIVALENT points only. Near non-Gamma-equivalent
%   integer points (e.g. (1,0,0), structure factor 0) the staggered branches retain
%   truncated-sum quality -- exploratory. Branch index = sorted-eigenvalue position per q;
%   branch identity is NOT tracked through crossings.
%
%   Returns:
%     P.Jnu     [nq x 4]  branch couplings (meV), guard applied
%     P.snapped [nq x 1]  logical: true where the directional limit replaced the raw sum
%     P.s       [1 x nq]  cumulative path distance in INDEX (r.l.u.) coordinates
%     P.s_cart  [1 x nq]  cumulative path distance in Cartesian reciprocal Ang^-1
%     P.ksnap   scalar    trust radius used (Ang^-1)
%
%   opts: .dpRng (30), .cache (true)  -- forwarded to invz_jq_modes for the raw sums;
%         .snapfac (2.5)              -- trust radius = snapfac*2*pi/(dpRng*min ||a_i||).
if nargin < 3, opts = struct(); end
dpRng    = 30;  if isfield(opts,'dpRng'),   dpRng    = opts.dpRng;   end
useCache = ~isfield(opts,'cache') || opts.cache;
snapfac  = 2.5; if isfield(opts,'snapfac'), snapfac  = opts.snapfac; end
C  = invz_const();
nq = size(qpath, 1);

Jnu = invz_jq_modes(ion, qpath, struct('dpRng', dpRng, 'cache', useCache));

Brec  = 2*pi*inv(ion.a).';                 % reciprocal basis rows: k_cart = q_rlu * Brec
ksnap = snapfac * 2*pi / (dpRng * min(vecnorm(ion.a, 2, 2)));
snapped = false(nq, 1);
Greg = [];                                  % regular Gamma-point matrix, computed lazily once
for iq = 1:nq
    G = round(qpath(iq, :));
    if ~is_gamma_equiv(G, ion.tau), continue; end
    k = (qpath(iq, :) - G) * Brec;
    if norm(k) >= ksnap, continue; end
    if norm(k) < 1e-12                      % exactly at G: use the local path direction
        if iq < nq, dq = qpath(iq+1, :) - qpath(iq, :);
        elseif iq > 1, dq = qpath(iq, :) - qpath(iq-1, :);
        else, dq = [0 0 0];
        end
        k = dq * Brec;
        if norm(k) < 1e-12, k = [1 0 0] * Brec; end   % single point at G: in-plane default
    end
    kz2 = (k(3) / norm(k))^2;
    if isempty(Greg)
        dip0 = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
        ex0  = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
        Greg = -squeeze(C.gfac*dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));
    end
    Jm = Greg + C.gfac*(4*pi/ion.Vc)*(1/3 - kz2);     % directional nonanalytic broadcast
    Jm = (Jm + Jm')/2;
    Jnu(iq, :) = sort(real(eig(Jm))).';
    snapped(iq) = true;
end

P.Jnu = Jnu;  P.snapped = snapped;  P.ksnap = ksnap;
P.s      = [0 cumsum(vecnorm(diff(qpath, 1, 1),        2, 2)).'];
P.s_cart = [0 cumsum(vecnorm(diff(qpath, 1, 1) * Brec, 2, 2)).'];
end

function tf = is_gamma_equiv(q, tau)
tf = abs(real(sum(exp(2i*pi*(tau*q.'))))/size(tau,1) - 1) < 1e-9;
end
```

- [ ] **Step 4: Create `invz/invz_spectra_qpath.m`**

```matlab
function S = invz_spectra_qpath(ion, T, B, qpath, w, opts)
%INVZ_SPECTRA_QPATH Branch-resolved chi''_cc(q, omega) along a q-path at fixed (T, B).
%   EXPLORATORY: this is a branch susceptibility -- ONE sorted coupling eigenvalue per q
%   fed to the scalar pole formula. It is NOT neutron scattering intensity (no eigenvector/
%   sublattice-interference weights, no polarization factor, no magnetic form factor), and
%   it inherits the known open issues of the medium solve (closed-grid BZ quadrature, the
%   real-axis continuation's possible negative weight, and the bare-MF FM/PM handoff of
%   invz_solve_auto). Energy TRENDS are comparable with R 2007 Fig 3; quantitative
%   reproduction is not claimed.
%
%   S = invz_spectra_qpath(ion, T, B, qpath, w) computes chi''_cc along qpath (nq x 3,
%   r.l.u.) at fixed T (K) and transverse field B (T), for the 1/z and bare-RPA theories.
%   The 1/z medium (Sigma, K, lambda) is solved ONCE at (T, B) on the BZ-integration grid
%   (ordered-first via invz_solve_auto); the path susceptibility then follows from the same
%   single-site response chit(w) via chi(q, w) = chit/(1 - J_nu(q) chit), with J_nu(q) from
%   invz_jq_path (direction-aware Gamma-limit guard; see its header).
%
%   Demag semantics (canonical): the strict-uniform Jshape_cc transform is NOT applied
%   here -- a finite-q probe measures the intrinsic longitudinal response, and at Gamma-
%   equivalent path points the relevant limit is q -> 0+ where the demagnetizing field
%   cancels (R 2007). The spectra CAN still change with ion.demag through the transverse
%   applied/internal-field relation (demag-aware info.Jaa0 -> solver Jxx0).
%
%   Returns:
%     S.chiz, S.chirpa     [nw x nq]  1/z and bare-RPA chi''_cc(q, w)
%     S.Epeak, S.Epeak_rpa [1 x nq]   censored, parabolic-refined peak energy per q (NaN
%                                     when the maximum is non-positive/non-finite or sits
%                                     in the first/last usable bin -- i.e. the true peak
%                                     lies outside the sampled window)
%     S.Jq      [1 x nq]   selected coupling branch along the path (meV)
%     S.snapped [1 x nq]   true where invz_jq_path replaced the raw truncated sum
%     S.s       [1 x nq]   path distance in INDEX (r.l.u.) coordinates
%     S.s_cart  [1 x nq]   path distance in Cartesian reciprocal Ang^-1
%     S.qpath, S.w, S.T, S.B, S.info, S.phase (1 = FM, 2 = PM solve used)
%
%   opts fields (all optional):
%     .grid ([16 16 16]), .dpRng (30), .eta (5e-3)   as in invz_spectra_map
%     .branch (4)        which ascending-sorted coupling branch to follow (sorted index,
%                        NOT a tracked mode identity through crossings)
%     .snapfac (2.5)     Gamma-limit trust-radius factor (see invz_jq_path)
%     .peak_wmin (0.05)  meV; excludes the low-frequency hyperfine pole (R 2007 Fig 2)
%                        from the peak search so Epeak tracks the doublet mode
%     .Jnu, .info        precomputed BZ-grid branches / info (skips the lattice sum; tests)
%   The hyperfine manifold is always included (electronuclear solve + matching real-axis
%   state); there is deliberately no .hyp option on this API.

if nargin < 6, opts = struct(); end
grid    = getf(opts, 'grid', [16 16 16]);
dpRng   = getf(opts, 'dpRng', 30);
eta     = getf(opts, 'eta', 5e-3);
branch  = getf(opts, 'branch', 4);
snapfac = getf(opts, 'snapfac', 2.5);
wmin    = getf(opts, 'peak_wmin', 0.05);

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
[pt, phase] = invz_solve_auto(ion, T, B, Jnu, struct('hyp', true, 'J0eff', Jcc0, 'Jxx0', Jaa0));
if phase == 0
    error('invz:noSolution', ...
        ['No converged 1/z solution at T = %.3f K, B = %.3f T ' ...
         '(near-degenerate doublet or critical band).'], T, B);
end

% guarded coupling branches along the path
P  = invz_jq_path(ion, qpath, struct('dpRng', dpRng, 'cache', true, 'snapfac', snapfac));
Jq = P.Jnu(:, branch).';

% path spectra from ONE real-axis evaluation each: Jsel is vectorized over q, and the
% single-site chit(w) is q-independent. Intrinsic: no Jshape here (see header).
copts = struct('Jsel', Jq, 'eta', eta, 'Jxx0', Jaa0, 'hyp', true);
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
S.qpath = qpath;  S.s = P.s;  S.s_cart = P.s_cart;  S.snapped = P.snapped(:).';
S.w = w;  S.T = T;  S.B = B;  S.phase = phase;  S.info = info;  S.Jq = Jq;
S.chiz = chiz;  S.chirpa = chirpa;
S.Epeak     = peak_energy(chiz,   w, wmin);
S.Epeak_rpa = peak_energy(chirpa, w, wmin);
end

% -------------------------------------------------------------------------------------------
function E = peak_energy(chi, w, wmin)
%PEAK_ENERGY per-column peak of chi''(w) at w >= wmin, parabolic-refined; CENSORED (NaN)
% when the maximum is non-positive/non-finite or sits in the first/last usable bin (a
% boundary maximum means the true peak lies outside the sampled window -- reporting the
% grid edge would fabricate a flat dispersion). Assumes uniform w spacing.
E = nan(1, size(chi, 2));
mask = w >= wmin;
wm = w(mask);
if numel(wm) < 3, return; end
dw = wm(2) - wm(1);
for k = 1:size(chi, 2)
    c = chi(mask, k);
    [cmax, i] = max(c);
    if ~isfinite(cmax) || cmax <= 0 || i == 1 || i == numel(c), continue; end
    d = (c(i-1) - c(i+1)) / (2*(c(i-1) - 2*c(i) + c(i+1)));   % parabolic vertex offset
    if ~isfinite(d) || abs(d) > 1, d = 0; end
    E(k) = wm(i) + d*dw;
end
end

% -------------------------------------------------------------------------------------------
function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
```

- [ ] **Step 5: Create `invz/invz_plot_spectra_qpath.m`**

```matlab
function invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit)
%INVZ_PLOT_SPECTRA_QPATH Render one branch chi''(q, omega) colormap panel (exploratory).
%   invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit) draws `chi` ([nw x nq], e.g.
%   S.chiz or S.chirpa) against the INDEX-coordinate path distance S.s (x; r.l.u. Miller
%   coordinates, not Cartesian reciprocal distance -- S.s_cart carries that) and frequency
%   S.w (y), and overlays the censored peak dispersion `Epeak` (pass S.Epeak or
%   S.Epeak_rpa) in white; censored (NaN) peak points simply leave gaps. This is a BRANCH
%   susceptibility, not neutron intensity (no structure-factor/form-factor weights).
%   Colour conventions match invz_plot_spectra_map: log10 scale spanning three decades
%   below the robust (99.5th-percentile) peak; NaN transparent on grey; present-but-
%   negative chi'' floored into the darkest colour. As there, the caller pre-scales S.w
%   and Epeak when plotting in GHz (see invz_run_spectra's eUnit knob).
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
xlabel(ax, sprintf('s along path from Q = [%g %g %g] (index r.l.u.)', S.qpath(1,:)));
ylabel(ax, sprintf('\\omega (%s)', eUnit));
title(ax, ttl);
cb = colorbar(ax);  cb.Label.String = 'log_{10} \chi''''_{cc}';
end
```

- [ ] **Step 6: Run test to verify it passes**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests/test_invz_spectra_qpath.m'); assertSuccess(results)"
```

Expected: PASS (4 tests). Debug hints if not:
- `S.Jq(3) ~= iref.Jcc0`: check the ascending sort (column 4 = max branch) and that the endpoint snap used `kz2 = 0` (path direction along x).
- Γ-regression cutoff mismatch at `h = 1.90`: at dpRng 20 that point IS snapped (trust radius 0.125 r.l.u.) while at dpRng 40 it is raw (radius 0.0625) — the 2% cross-cutoff tolerance covers the measured 1.6% gap. Do not silently widen tolerances; investigate first.

- [ ] **Step 7: Commit (explicit paths)**

```bash
git add invz/invz_jq_path.m invz/invz_spectra_qpath.m invz/invz_plot_spectra_qpath.m invz/tests/test_invz_spectra_qpath.m
git commit -m "feat(invz): exploratory q-path dispersion engine with direction-aware Gamma-limit guard + censoring peak picker

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
%   qpath = [nq x 3] r.l.u.  -- EXPLORATORY branch-resolved q-path view at fixed field(s)
%     Bq, for comparison with the TRENDS in R 2007 Fig 3 (branch susceptibility, not
%     neutron intensity; see invz_spectra_qpath header for the inherited caveats):
%     numel(Bq) == 1  ->  2D path-vs-frequency colormaps (1/z + RPA), censored peak overlay
%     numel(Bq) >  1  ->  E_peak(q) dispersion overlay, one colour per field
%
% The field sweep is done by invz_spectra_map (both phases: FM below Bc, PM above; S.phase
% labels each field 1=FM/2=PM/0=masked); the q-path by invz_spectra_qpath (one 1/z medium
% solve per field, then the pole formula along the guarded path couplings). Both result
% structs replot for free, e.g. invz_plot_spectra_map(gca, S, S.chiz, '1/z').
% Save with:  print(gcf, 'spectra.png', '-dpng', '-r150');

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

ion = invz_ion();
% ion.demag = 1; ion.alpha = 0.5;    % OPTIONAL sample-shape (demagnetization) knob; default off.
%   demag = 0 (default): intrinsic response -- the R 2007 benchmark. demag ~= 0 with spheroid
%   aspect ratio alpha (1 sphere, 0 c-needle, Inf disk) adds the sample shape as follows:
%     - info.Jcc0, Jnu, and the ordering-channel contribution to criticality are
%       demag-INVARIANT (R 2007: the demagnetizing field cancels from the critical
%       condition; ordering occurs at q -> 0+, not strict q = 0);
%     - Tc(B = 0) is EXACTLY demag-invariant (the transverse moment vanishes there);
%     - Bc(T) vs APPLIED transverse field CAN shift through the demag-aware transverse
%       coupling info.Jaa0 (internal-vs-applied field relation);
%     - the strict-uniform (q = 0) observable chi''_cc is demag-corrected via info.Jshape_cc
%       (chi_meas = chi/(1 + Jshape*chi)): the soft mode saturates instead of diverging;
%     - q-path spectra omit the Jshape_cc transform (finite-q probe = intrinsic longitudinal
%       response) but still see demag through info.Jaa0.
T = 0.31;                             % K
% fields = [3.6 4.2 4.8 5.4 6.0];       % few -> slices;  many -> colormap
fields = linspace(3,6.5,151);
w = (0:0.002:0.45).';               % meV -- field-sweep views
eta = 1e-3;                          % real-axis line width: Lorentzian HWHM in meV (5e-3 meV ~ 1.2 GHz).
                                     % Lower -> sharper peaks (resolves the sub-6-GHz hyperfine lines),
                                     % but keep eta >~ the w-spacing above or the peaks alias.
sliceMax = 6;                         % field count at/below which the line-slice view is used
useParallel = true;                  % true -> parfor over fields (Parallel Computing Toolbox)
eUnit = 'meV';                       % 'meV' or 'GHz' -- plotting only; computation always runs in meV

% ---- q-path view (R 2007 Fig 3 trends): set qpath non-empty to switch views -------------
qpath = [];                          % [] = field-sweep views; [nq x 3] r.l.u. = q-path view
% qh = linspace(1, 2, 51).';  qpath = [qh zeros(numel(qh), 2)];  % (1,0,0)->(2,0,0), Fig 3 path
Bq = 4.24;                           % field(s), T, for the q-path view. One value -> colormaps;
                                     % several -> E_peak(q) overlay (R 2007 Fig 3: [3.6 4.24 6.0])
wq = (0:0.004:0.85).';               % meV -- q-path grid. Fig 3 reaches ~0.75 meV near h = 1 at
                                     % 60 kOe (after their 1.15 scaling); 0.85 avoids clipping,
                                     % which the censoring peak picker would flag as NaN.
dispScale = 1;                       % dispersion display scale factor; R 2007 scales the
                                     % calculated energies by 1.15 to match experiment (Fig 3)

C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end

if ~isempty(qpath)
    % ---------------- exploratory q-path view at fixed field(s) ----------------
    if numel(Bq) == 1
        S = invz_spectra_qpath(ion, T, Bq, qpath, wq, struct('eta', eta));
        Splot = S;   % display-only copy; the solve above always ran in meV
        Splot.w = S.w*eScale;  Splot.Epeak = S.Epeak*eScale;  Splot.Epeak_rpa = S.Epeak_rpa*eScale;
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);
        invz_plot_spectra_qpath(ax1, Splot, Splot.chiz, Splot.Epeak, ...
            sprintf('1/z branch \\chi''''_{cc} (exploratory), T = %.2f K, B = %.2f T', T, Bq), eUnit);
        ax2 = subplot(1, 2, 2);
        invz_plot_spectra_qpath(ax2, Splot, Splot.chirpa, Splot.Epeak_rpa, ...
            sprintf('RPA branch \\chi''''_{cc} (exploratory), T = %.2f K, B = %.2f T', T, Bq), eUnit);
    else
        figure; hold on;  co = lines(numel(Bq));
        for k = 1:numel(Bq)
            Sk = invz_spectra_qpath(ion, T, Bq(k), qpath, wq, struct('eta', eta));
            plot(Sk.s, Sk.Epeak*eScale*dispScale,     '-',  'Color', co(k, :), ...
                 'DisplayName', sprintf('1/z, %.2f T', Bq(k)));
            plot(Sk.s, Sk.Epeak_rpa*eScale*dispScale, '--', 'Color', co(k, :), ...
                 'DisplayName', sprintf('RPA, %.2f T', Bq(k)));
        end
        xlabel(sprintf('s along path from Q = [%g %g %g] (index r.l.u.)', qpath(1, :)));
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('branch dispersion (exploratory), T = %.2f K, dispScale = %.2f', T, dispScale));
        legend show;
        % cf. R 2007 Fig 3 TRENDS: for the (1,0,0)->(2,0,0) path, Q = (1+s, 0, 0); their
        % plotted theory lines are the calculated energies scaled by 1.15 (dispScale = 1.15).
        % Gaps in the lines are CENSORED peaks (mode outside the wq window) -- widen wq,
        % do not interpolate over them.
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

(a) Immediately after line 43 (`ion = invz_ion();`), insert the knob comment (canonical wording):

```matlab
% ion.demag = 1;  ion.alpha = 0.25;  % OPTIONAL sample-shape (demagnetization) knob; default off.
%   - info.Jcc0, Jnu, and the ordering-channel contribution to criticality are demag-
%     INVARIANT (R 2007: the demagnetizing field cancels from the critical condition,
%     since ordering occurs at q -> 0+, not strict q = 0).
%   - Tc(B = 0) is EXACTLY demag-invariant: the transverse moment vanishes there.
%   - Bc(T) vs APPLIED transverse field CAN shift through the demag-aware transverse
%     coupling info.Jaa0 (hoisted into Jxx0 below): internal-vs-applied field relation.
%   demag = 0 (default) is the intrinsic / internal-field boundary matching the R 2007
%   benchmark.
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

- [ ] **Step 3: Smoke-run the q-path engine (cheap settings, endpoint-dense)**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath('invz'); addpath('.'); ion = invz_ion(); w = (0.02:0.02:0.85).'; qh = [linspace(1,1.9,7) 1.95 1.99 2.0].'; S = invz_spectra_qpath(ion, 0.31, 5.5, [qh zeros(numel(qh),2)], w, struct('grid',[8 8 8],'dpRng',10)); disp([qh.'; S.Jq; S.Epeak]); assert(all(isfinite(S.Epeak)) && S.Epeak(1) > S.Epeak(end)); assert(max(abs(diff(S.Jq(end-2:end)))) < 0.1*S.Jq(end))"
```

Expected: prints path points, guarded couplings, and finite peak energies decreasing from (1,0,0) toward (2,0,0); both assertions pass (the second confirms no endpoint dip/jump survives the guard). First run pays one 8^3 lattice sum + one 1/z solve: about a minute.

- [ ] **Step 4: Run the full fast suite** (guards against driver-adjacent regressions)

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: PASS.

- [ ] **Step 5: Commit (explicit paths; note pre-existing hunks if Task 0's checkpoint was declined)**

```bash
git add invz/invz_run_spectra.m invz/invz_run_phase_diagram.m
git commit -m "feat(invz): demag knob in both drivers; exploratory q-path view (R2007 Fig 3 trends) in invz_run_spectra

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
% 6^3 grid, Gamma excluded: the invariances below are exact, so grid quality is irrelevant.
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
% that is the internal-vs-applied transverse field correction (canonical bullet 3). SLOW.
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

Expected: PASS (2 passed). Runtime: a few minutes (three full 1/z solves at 0.31 K on a 6^3 grid).

- [ ] **Step 4: Commit (explicit paths)**

```bash
git add invz/tests/test_invz_demag_invariance.m
git commit -m "test(invz): lock phase-boundary demag invariance (ordering channel) + transverse-only entry

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Documentation sync + final verification

**Files:**
- Modify: `invz/invz_ion.m`, `invz/README.html`

- [ ] **Step 1: Fix the `invz_ion.m` knob comment (canonical wording)**

Replace lines 14–17 of `invz/invz_ion.m` (the block starting `% Sample-shape (demagnetization) correction ...`) with:

```matlab
% Sample-shape (demagnetization) knob. demag = 0 (default): intrinsic couplings, the R2007
% benchmark. demag ~= 0: the ellipsoid shape (ellipsoid_demagn(alpha)) enters ONLY as
%   (a) info.Jshape_cc -- strict-uniform observable correction applied in invz_chi_realaxis
%       (chi_meas = chi/(1 + Jshape_cc*chi)), and
%   (b) demag-aware info.Jaa0 -- the transverse mean-field channel.
% Consequences: info.Jcc0/Jnu and the ordering-channel criticality are demag-INVARIANT
% (R2007: the demagnetizing field cancels from the critical condition; ordering at q -> 0+);
% Tc(B=0) is exactly demag-invariant (<Jx> = 0 there); Bc(T) vs APPLIED field can still
% shift through (b) -- the internal-vs-applied transverse field relation.
```

(keep lines 18–19, the `ion.demag = 0;` / `ion.alpha = 1;` defaults, unchanged).

- [ ] **Step 2: Update `invz/README.html`**

Locate the demag callout near line 130–136 (the `Shape scope` note and the displayed `J(0) -> J(0) + Lorentz - demag` equation). Replace that whole demag block so it states:

1. The corrected equation set: `Jnu`/`info.Jcc0` carry the Lorentz cavity only; `info.Jshape_cc = 4*(4π/Vc)*gfac*demag*N_zz(α)` corrects the strict-uniform observable in `invz_chi_realaxis` (`chi_meas = chi/(1+Jshape_cc*chi)` — algebraically identical to the old shifted pole `chit/(1-(Jcc0-Jshape_cc)*chit)`); `info.Jaa0 = Jaa0_dipole(demag-aware) + 4*J12` feeds the transverse mean field via `opts.Jxx0`.
2. Phase-boundary behavior in the CANONICAL four-bullet wording from this plan's context section, attributed as "per Rønnow et al. 2007 (who note this contradicts a suggestion of Chakraborty et al.)".
3. Remove the stale "wire info.Jaa0 into the transverse MF" TODO (it is now wired) and the stale reference to an `info.Jaa0` field that "doesn't exist" — it exists now.
4. Add `invz_jq_path` / `invz_spectra_qpath` / `invz_plot_spectra_qpath` / `invz_solve_auto` rows to the function table, and a short q-path quick-start (Fig. 3 path snippet from the driver, `wq` up to 0.85 meV, `dispScale = 1.15` comparison). The quick-start MUST carry the exploratory framing: branch-resolved susceptibility (not neutron intensity), sorted-index branch selection (no crossing tracking), the Γ-limit guard and its trust radius, censored peaks, and the inherited open issues (BZ quadrature, real-axis positivity, bare-MF phase handoff).
5. Note the cache-format bump: old `invz/cache/jq_*.mat` files are orphaned (prefix now `jq2_`) and can be deleted.

This is prose editing of an HTML file — match the existing `<div class="callout note">`/table markup style already used in the file.

- [ ] **Step 3: Full verification — fast suite, then slow suite**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: fast suite all green (now ~40+ passed). Slow suite all green (~6–10 min; `test_critical_field_at_310mK` re-derives Bc with the live `Jxx0` — the published-range bounds [4.0, 4.6] T absorb the <0.1% transverse-coupling drift).

- [ ] **Step 4: Commit (explicit paths)**

```bash
git add invz/invz_ion.m invz/README.html
git commit -m "docs(invz): canonical demag semantics; exploratory q-path quick-start; cache-bump note

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

- [ ] **Step 5: Report**

Summarize for the user: what changed physically (demag no longer moves the ordering-channel criticality; observable + transverse channels carry the shape, per the canonical four bullets), the two driver knobs, the exploratory q-path workflow (`qpath` + `Bq = [3.6 4.24 6.0]` + `wq` to 0.85 meV + `dispScale = 1.15`), the Γ-limit guard and what it does/doesn't cover, the expected small `Jxx0` drift, the cache invalidation, and which review findings were addressed vs deliberately deferred (Ewald, ODD terms, Codex #2/#3/#4).
