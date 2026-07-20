# invzt Ordered-Side (FM) a1 + dominant-vertex a3d Implementation Plan — REVISED after the two 2026-07-19 Codex reviews

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the full-tensor branch an ordered-phase (FM) solve for mode `a1` end-to-end (point solver → real-axis spectra → **stability-based** phase dispatcher → rewired `invzt_run_spectra`), so the spectrum driver spans the quantum phase transition and shows the electro-nuclear excitation modes on BOTH sides — plus **`a3d`: the full-response, fixed-rank field-adapted dominant-vertex ordered a3** (Matsubara-only, rank-16 field-adapted vertex basis, compact `cc;cc` storage, hard performance budget; Task 7 stages 7A–7E). Full-dress 136-state ordered a3 remains out of scope (budget-refused).

**Architecture:** Mirror the projected branch's ordered stack (`invz_solve_point_ordered` / `invz_sigma_ordered` / `invz_twolevel_ordered` / the ordered branch of `invz_chi_realaxis`) onto the tensor lattice engine (`invzt_gcc_lattice` 12×12 RPA). The moment-form self-energy pieces move to `invz_common/` (`git mv`, behavior-neutral). The ordered single-ion state comes from the EXISTING `invz_single_ion(..., 'order', true)` machinery with tensor-native couplings (`J0z = lat.info.Jcc0`, `Jxx0 = lat.info.Jaa0`). **Phase selection is stability-based, NOT ordered-first** (review P0-1): the PM leg's `crit > 0` is the tensor QPT criterion and decides FIRST; the ordered leg is consulted only when the PM sample is invalid.

**Tech Stack:** MATLAB R2025a (`MB="/Applications/MATLAB_R2025a.app/bin/matlab -batch"`), repo test convention `functiontests(localfunctions)` + `runtests`; python3 + mpmath for the vertex oracle (`verify_tensor_vertex.py`).

**Scope decisions:**
1. (user, 2026-07-19) **Transverse only:** spontaneous-moment route at `|Bz| <= bz_tol` only; longitudinal raises `invzt:orderedLongitudinal`; no `forced_moment` port.
2. (user, 2026-07-19, round 2 — supersedes the round-1 deferral; basis definition amended 2026-07-20 after the 7B falsification, see 7B) **Ordered a3 = `a3d`, the full-response dominant-vertex hybrid:** full 136-state ordered single-ion response and lattice medium; the four-point vertex computed ONLY in a FIXED-COUNT 16-state FIELD-ADAPTED dominant basis (the 16 lowest ordered mean-field eigenstates — NOT a moving `E < Esplit` cut, which drifts 13→8 states, and NOT the zero-field `e2xI8` content, which captures ~0% of the ordered cc response, both measured); only the `cc;cc` vertex stored (`[nwn,nl]`, no zero-padded `[9,3,3,nwn,nl]` 1.4 GB array); frequency work tiled and budget-guarded (`invzt:orderedA3Budget`); **Matsubara-only** (`invzt:realaxisMode` keeps rejecting it); production 0.1 K runs under `INVZ_SLOW`, never CORE. `mode='a3d'` in `invzt_solve_point_ordered` (Task 7D); full dress and `mode='a3'` are rejected (`invzt:orderedMode`). The dispatcher classifies the phase with `a1` ONLY — `a3d` is evaluated after classification, never to decide the phase.
3. Documentation names this "full-response, dominant-vertex ordered `a3d`" — never "full tensor ordered a3" — so it cannot be mistaken for the refused 136-state full-dress calculation.
4. The ordered Jensen moment-form (eqs 37–38) is an **approximation diagnostic** for `a3d`, not an exact gate: the derivation replaces `g(ωm±ωn)` by its static value (exact at n=0 only) and treats the elastic sector in a resummed static approximation (`G0 = -M²g - m²h`, J 2.7/2.8). Measured residuals of the naive bridge: 3.13e-5 (n=0) but 4.89e-2 (n=1), 8.43e-2 (n=2) at βΔ=14. Gates are static/asymptotic only; finite-frequency comparisons are REPORT.

## Global Constraints

- Tensor CORE suite (`invz_tensor/tests`) runs with `invz_projected` absent from the path; parity tests live in `invz_tensor/tests/interop/` only.
- Physics numbers are *reported, never tuned*; hard gates are exact identities, parity reductions, symmetry, convergence, reproducibility.
- Coupling threading (v3): ordered solve consumes `lat.info.Jcc0` / `lat.info.Jaa0`; record `pt.J0z_used` / `pt.Jxx0_used`.
- Elastic conventions: Matsubara medium = elastic ON (`invzt_chi0_split` default); real-axis rebuild = elastic OFF.
- Commit style: `feat(invzt):` / `test(invzt):` / `docs(invzt):` / `refactor(invzt):`. **Every commit stages explicit paths — never `git add -A`** (review P1-6: the worktree carries unrelated user work).
- **Suite invariant (review P1-5):** the repo baseline is CORE 47 passed / 0 failed / **1 Incomplete** (`INVZ_SLOW`-gated), INTEROP 8/0/2, PROJECTED 143/0/**19** — so suite gates assert **zero failures**, never all-passed: `assert(~any([r.Failed]))`, printing `{r([r.Incomplete]).Name}` for eyeball comparison against the known allowlist. Final acceptance (Task 8) checks the Incomplete set matches the allowlist exactly and runs the slow suite with `INVZ_SLOW=1`.
- Measured anchors (T = 0.1 K, gridN 16 halfopen, dpRng 30, odd on, legacy_x; verified 2026-07-19): tensor QPT (crit sign change) between **4.65 and 4.70 T**; bare-MF ordered moment persists to ~5.0 T (`m0` = 1.6553 / 1.5109 / 1.1717 at 4.65 / 4.70 / 4.80 T — MF boundary ≠ QPT, the P0-1 driver); PM solver **converges to spurious negative-crit points below Bc** (4.0 T: converged, crit = −0.186); PM soft mode 0.2595 meV @ 4.8 T → 0.3435 meV @ 6.0 T.
- All commands assume repo root = `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion`.

---

### Task 1: Move the moment-form engine to `invz_common/` (behavior-neutral)

**Files:**
- Move: `invz_projected/invz_sigma_ordered.m` → `invz_common/invz_sigma_ordered.m`
- Move: `invz_projected/invz_twolevel_ordered.m` → `invz_common/invz_twolevel_ordered.m`
- Modify: `invz_projected/README.html` (both functions are currently documented as projected-module residents — 3 references; review P1-6)

**Interfaces:**
- Produces: `sig = invz_sigma_ordered(tl, lam, K, g, beta)` (fields `alpha, alpha_m, gamma, Sigma`; requires `lam = invz_lambdas(K, g, wts, beta, [1 2 3])`); `tl = invz_twolevel_ordered(ion, T, Bx, hz, opts)` (PM `tl` fields plus nonzero-allowed `m`, `h0`, `hz`). Tasks 3–4 and 7 call both from `invz_tensor` (CORE path).

- [ ] **Step 1: Verify both functions are projected-free**

Run: `grep -nE "invz_emt_scalar|invz_jq|invz_odd|invz_chiperp|invz_solve" invz_projected/invz_sigma_ordered.m invz_projected/invz_twolevel_ordered.m`
Expected: no output.

- [ ] **Step 2: git mv both files; update the projected README**

```bash
git mv invz_projected/invz_sigma_ordered.m invz_common/invz_sigma_ordered.m
git mv invz_projected/invz_twolevel_ordered.m invz_common/invz_twolevel_ordered.m
```

Then edit `invz_projected/README.html`: move the two function entries to wherever it lists `invz_common` residents (or annotate them "moved to invz_common 2026-07-19, one source of truth"), matching how the 2026-07-17 single-ion-engine move is documented.

- [ ] **Step 3: Projected suite — zero failures, baseline incompletes unchanged**

Run: `$MB "cd('invz_projected/tests'); r = runtests(pwd); disp(table(r)); fprintf('failed=%d incomplete=%d\n', nnz([r.Failed]), nnz([r.Incomplete])); assert(~any([r.Failed]))"`
Expected: `failed=0 incomplete=19` (the known `INVZ_SLOW` allowlist — compare names if the count moves).

- [ ] **Step 4: Tensor CORE — zero failures**

Run: `$MB "cd('invz_tensor/tests'); r = runtests(pwd); disp(table(r)); assert(~any([r.Failed]))"`
Expected: `failed=0 incomplete=1`.

- [ ] **Step 5: Commit (path-scoped — the worktree carries unrelated user work)**

```bash
git add invz_common/invz_sigma_ordered.m invz_common/invz_twolevel_ordered.m invz_projected/invz_sigma_ordered.m invz_projected/invz_twolevel_ordered.m invz_projected/README.html
git commit -m "refactor(invz): move moment-form engine (sigma_ordered, twolevel_ordered) to invz_common -- branch-agnostic, behavior-neutral"
```

---

### Task 2: Extract the static criticality helper `invzt_crit_static.m`

**Files:**
- Create: `invz_tensor/invzt_crit_static.m`
- Modify: `invz_tensor/invzt_solve_point.m` (local function `local_crit` → thin call)

**Interfaces:**
- Produces: `[crit, cmass, arank] = invzt_crit_static(ctil0, JtGamma, rank_tol)` — EXACTLY the current `local_crit` body. Consumed by `invzt_solve_point` (existing) and `invzt_solve_point_ordered` (Task 3).

- [ ] **Step 1: Create the new file** — copy `local_crit` from `invzt_solve_point.m:396-410` verbatim:

```matlab
function [crit, cmass, arank] = invzt_crit_static(ctil0, JtGamma, rank_tol)
%INVZT_CRIT_STATIC Static Gamma-point criticality of a renormalized local 3x3.
%   [crit, cmass, arank] = INVZT_CRIT_STATIC(ctil0, JtGamma, rank_tol) returns
%   the min real eigenvalue of (M+M')/2, M = I - S*JtGamma*S, with S the
%   rank-clipped PSD square root of C12 = kron(eye(4), herm(ctil0)) (Hermitian
%   eigendecomposition, NOT sqrtm). crit shares the zeros of I - C12*JtGamma on
%   the active subspace; crit > 0 in a locally stable phase (PM above Bc, and
%   the converged FM state away from Bc). cmass = clipped eigenvalue mass,
%   arank = active rank. Extracted verbatim from INVZT_SOLVE_POINT local_crit
%   (2026-07-19) so the ordered solver shares one implementation.
if nargin < 3, rank_tol = 1e-12; end
C12 = kron(eye(4), (ctil0 + ctil0')/2);
[U, D] = eig((C12 + C12')/2);
d = real(diag(D));
clip = d < rank_tol;
cmass = sum(abs(d(clip)));
arank = sum(~clip);
d(clip) = 0;
S = U * diag(sqrt(max(d, 0))) * U';
M = eye(size(S,1)) - S*JtGamma*S;
crit = min(real(eig((M + M')/2)));
end
```

- [ ] **Step 2: Replace `local_crit` in `invzt_solve_point.m`** — delete the local function (lines 396–410) and repoint both call sites:

```matlab
[crit, crit_clipped_mass, crit_active_rank] = invzt_crit_static(ctil0, lat.JtGamma, rank_tol);
if strcmp(mode, 'a3')
    crit_add = invzt_crit_static(ctil0_add, lat.JtGamma, rank_tol);
    resum_spread_crit = crit - crit_add;
end
```

- [ ] **Step 3: CORE + INTEROP — zero failures (bit-identical extraction gate)**

Run: `$MB "cd('invz_tensor/tests'); r = runtests(pwd); assert(~any([r.Failed]))"` then `$MB "cd('invz_tensor/tests/interop'); r = runtests(pwd); assert(~any([r.Failed]))"`
Expected: zero failures both (pure code motion); incompletes at baseline (1 / 2).

- [ ] **Step 4: Commit**

```bash
git add invz_tensor/invzt_crit_static.m invz_tensor/invzt_solve_point.m
git commit -m "refactor(invzt): extract invzt_crit_static from invzt_solve_point -- shared with the ordered solver"
```

---

### Task 3: `invzt_solve_point_ordered.m` — the tensor a1 FM point solver

**EXECUTION AMENDMENT (2026-07-20, user-approved — supersedes the split lines in the code below):** the ordered medium is **WHOLE-CC** — no dominant/rest split. Measured during execution: the `E < Esplit` cut selects a drifting 13→8-state set carrying only **47.6% of the cc weight** at 0.1 K ordered points, which made the 5e-3 parity gate structurally unsatisfiable (dSigma0 = 8.6e-3, dalpha_m = 5.8e-2, with `chi0cc0` bit-identical and `dm0 = 0` — pure medium-convention divergence). The projected ordered reference itself renormalizes the FULL electronuclear `G0` (no split); the dominant/rest rule is the PM/E1 prescription. Therefore:
- Replace the `invzt_chi0_split` block with `c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true))`; `c0_cc = real(squeeze(c0(3,3,:)))`; `chi0cc0 = c0_cc(1)`.
- Loop: `ctil = c0 ./ denom(Sigma)`; `G0til = -(c0_cc ./ (1 + Sigma))`; K bookkeeping unchanged. `ctil0 = c0(:,:,1)/(1 + Sigma(1))`.
- REMOVE the `Esplit`/`chi_rest` options; passing either errors `invzt:orderedSplitKnobs` (fail-loud: the ordered solver is whole-cc by design). `pt.chi_rest = true` kept for interface parity (whole response renormalized), `pt.mspec = []`.
- `odd` stays, default true (production); the parity test's tensor leg passes `struct('odd', false)` only.
- The test file's `local_one_pass` helper uses the same whole-cc map; add a small test asserting `invzt:orderedSplitKnobs`. All gates unchanged (crit > −1e-3, sumrule < 0.1, parity 5e-3).
- The fixed dominant basis (e2xI8) is a3d's vertex concern (Task 7B) — NOT a1's.

**Files:**
- Create: `invz_tensor/invzt_solve_point_ordered.m`
- Test: `invz_tensor/tests/test_invzt_solve_point_ordered.m`
- Test: `invz_tensor/tests/interop/test_invzt_solve_ordered_parity.m`

**Interfaces:**
- Consumes: `invz_single_ion(ion,T,B,struct('hyp',hyp,'order',true,'J0z',J0z,'Jxx0',Jxx0,'transverse_mf',tmf))`; `invz_twolevel_ordered(ion,T,B,si.hz,struct('Jxx0',Jxx0,'transverse_mf',tmf))` (Task 1); `invzt_chi0_split`; `invzt_gcc_lattice`; `invz_lambdas(K,g,wts,beta,[1 2 3])`; `invz_sigma_ordered` (Task 1); `invzt_odd_mask`; `invzt_crit_static` (Task 2).
- Produces: `pt = invzt_solve_point_ordered(ion, T, B, lat, opts)`. Accepted path: the `invzt_solve_point` field surface (`Sigma0, Sigma, alpha, lambda(3x1), K, G, tl, si, chi0cc0, crit, crit_clipped_mass, crit_active_rank, sumrule_rel, converged, outer_iters, diag4_spread, mode='a1', odd, chi_rest, mspec, Jxx0_used, lat` provenance) PLUS `m0`, `alpha_m`, `is_ordered=true`, `moment_branch='spontaneous'`, `J0z_used`. Early returns carry the projected-parity fixed set: `m0, is_ordered=false, converged=false, Sigma0=NaN, crit=NaN, si, tl=[], moment_branch`. `opts.mode` other than `'a1'` errors `invzt:orderedMode` (P0-2 scope). Consumed by Tasks 4–6.

- [ ] **Step 1: Write the failing tests** — `invz_tensor/tests/test_invzt_solve_point_ordered.m`:

```matlab
function tests = test_invzt_solve_point_ordered
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
tc.TestData.ion = ion;
tc.TestData.lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
end

function test_deep_ordered_point(tc)
% Deep FM point: spontaneous moment, converged Sigma, stable (crit >= 0 within tol).
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], tc.TestData.lat, struct());
verifyTrue(tc, pt.is_ordered);
verifyTrue(tc, pt.converged);
verifyEqual(tc, pt.moment_branch, 'spontaneous');
verifyGreaterThan(tc, abs(pt.m0), 1.0);              % LiHoF4 FM moment is O(5) here
verifyTrue(tc, isfinite(pt.Sigma0) && all(isfinite(pt.Sigma)));
verifyTrue(tc, isfinite(pt.alpha_m));
verifyGreaterThan(tc, pt.crit, -1e-3);               % broken-symmetry state locally stable
verifyEqual(tc, pt.J0z_used, tc.TestData.lat.info.Jcc0);
verifyEqual(tc, pt.Jxx0_used, tc.TestData.lat.info.Jaa0);
verifyLessThan(tc, pt.sumrule_rel, 0.1);
fprintf('ordered a1 @ 3T: m0=%.4f Sigma0=%.5f crit=%.5f alpha_m=%.4g sumrule=%.3g\n', ...
    pt.m0, pt.Sigma0, pt.crit, pt.alpha_m, pt.sumrule_rel);
end

function test_converged_state_is_self_consistent(tc)
% Review P2-2: on a converged exit, the returned (Sigma, K, lambda, alpha_m) must
% describe the SAME state -- re-evaluating one medium+self-energy pass at pt.Sigma
% must reproduce pt.Sigma to the outer tolerance (check-before-mix loop ordering).
ion = tc.TestData.ion;  lat = tc.TestData.lat;
pt = invzt_solve_point_ordered(ion, 0.1, [3.0 0 0], lat, struct());
verifyTrue(tc, pt.converged);
sg2 = local_one_pass(ion, 0.1, [3.0 0 0], lat, pt);   % helper below: one g(Sigma) map step
verifyLessThan(tc, max(abs(sg2.Sigma - pt.Sigma)), 1e-6);   % >= tol_outer scale, not machine
end

function sg = local_one_pass(ion, T, B, lat, pt)
% One medium + moment-form self-energy evaluation AT pt.Sigma (no mixing): the
% fixed-point map g(Sigma) evaluated at the returned state.
[wn, wts, beta] = invz_matsubara(T, 40);
nwn = numel(wn);
[cdom, crest] = invzt_chi0_split(pt.si, T, 1i*wn, struct('Esplit', 0.4653));
cdom_cc  = real(squeeze(cdom(3,3,:)));
crest_cc = real(squeeze(crest(3,3,:)));
g = real(invz_g(pt.tl, 1i*wn));
ctil = cdom ./ reshape(1 + pt.Sigma, 1, 1, nwn) + crest;
Gcc = invzt_gcc_lattice(ctil, lat);
Gloc = -Gcc(:);
G0til = -(cdom_cc ./ (1 + pt.Sigma) + crest_cc);
K = 1 ./ Gloc - 1 ./ G0til;
lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg = invz_sigma_ordered(pt.tl, lam, K, g, beta);
end

function test_pm_early_return(tc)
% WELL above the bare-MF boundary (~5.0 T -- NOT the QPT at 4.65-4.70; measured
% m0 = 1.17 at 4.8 T, so 4.8 would NOT early-return) the order-mode MF relaxes
% to ~0: paramagnetic early return.
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [5.5 0 0], tc.TestData.lat, struct());
verifyFalse(tc, pt.is_ordered);
verifyFalse(tc, pt.converged);
verifyEqual(tc, pt.moment_branch, 'none');
verifyTrue(tc, isnan(pt.Sigma0) && isnan(pt.crit));
verifyTrue(tc, isempty(pt.tl));
verifyLessThan(tc, abs(pt.m0), 1e-2);
end

function test_mf_moment_persists_past_QPT(tc)
% Locks the P0-1 finding into the suite: the bare-MF moment is STILL nonzero at
% the 4.8 T PM anchor (measured m0 = 1.1717) -- which is exactly why phase
% selection must NOT be ordered-first (see invzt_solve_auto, Task 5).
pt = invzt_solve_point_ordered(tc.TestData.ion, 0.1, [4.8 0 0], tc.TestData.lat, struct());
verifyTrue(tc, abs(pt.m0) > 1.0);
fprintf('P0-1 lock: bare-MF m0(4.8 T) = %.4f (QPT is at 4.65-4.70 T)\n', pt.m0);
end

function test_longitudinal_rejected(tc)
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0.1], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
end

function test_full_dress_a3_rejected(tc)
% Full-dress 'a3' is PERMANENTLY rejected (136-state vertex budget-refused).
% This assertion stays valid after Task 7D adds 'a3d' to the mode gate.
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('mode', 'a3')), 'invzt:orderedMode');
end

function test_nlevels_std_only(tc)
verifyError(tc, @() invzt_solve_point_ordered(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('nlevels', 'three')), 'invzt:orderedNlevels');
end
```

- [ ] **Step 2: Run to verify failure**

Run: `$MB "cd('invz_tensor/tests'); r = runtests('test_invzt_solve_point_ordered'); disp(table(r))"`
Expected: FAIL — `invzt_solve_point_ordered` undefined.

- [ ] **Step 3: Implement `invz_tensor/invzt_solve_point_ordered.m`**

```matlab
function pt = invzt_solve_point_ordered(ion, T, B, lat, opts)
%INVZT_SOLVE_POINT_ORDERED Tensor a1 solve at one FERROMAGNETIC (T, B) point.
%   pt = INVZT_SOLVE_POINT_ORDERED(ion, T, B, lat, opts) is the ordered-phase
%   counterpart of INVZT_SOLVE_POINT mode 'a1': the single-ion problem is solved
%   with the longitudinal ORDERING mean field (spontaneous m0 = <Jz>, hz =
%   J0z*<Jz> with J0z = lat.info.Jcc0 -- tensor-native threading, NOT ion.J0eff),
%   the self-energy uses the moment form (INVZ_SIGMA_ORDERED, HTML eqs 37-38,
%   lambda_{1,2,3}), and the medium step is the SAME tensor lattice engine as
%   the PM solver (INVZT_GCC_LATTICE 12x12 RPA + LOCKED K bookkeeping).
%
%   SCOPE: TRANSVERSE fields only (spontaneous route; invzt:orderedLongitudinal
%   otherwise -- no forced_moment port, 2026-07-19 decision). Modes: 'a1' (this
%   task) and, once Task 7D lands, 'a3d' (full-response dominant-vertex,
%   Matsubara-only); full-dress 'a3' is PERMANENTLY rejected (invzt:orderedMode,
%   136-state vertex budget-refused). nlevels 'std' only
%   (invzt:orderedNlevels). Option-A parity with the projected ordered solver:
%   m0 is the BARE mean-field order parameter -- it onsets at the MF boundary
%   (~5.0 T at 0.1 K, MEASURED 2026-07-19), well above the tensor crit = 0 QPT
%   (4.65-4.70 T). A converged ordered point is therefore NOT evidence of being
%   in the FM phase; phase assignment belongs to INVZT_SOLVE_AUTO's
%   stability-based rule, never to this solver.
%
%   Mixing: plain damped mixing, CHECK-BEFORE-MIX (same loop ordering as the PM
%   solver, review P2-2): on every exit -- converged or not -- the returned
%   (Sigma, K, G, lambda, alpha, alpha_m) describe the SAME evaluated state.
%
%   OPTIONS (getf defaults): Ecut 40, hyp true, transverse_mf 'legacy_x',
%   mix_outer 0.7, tol_outer 1e-8, max_outer 80, m_tol 1e-2, bz_tol 1e-9,
%   odd true, chi_rest true, Esplit 0.4653, rank_tol 1e-12,
%   mz_seed/mf_maxit/mf_mix forwarded to invz_single_ion.
%
%   Returns the INVZT_SOLVE_POINT pt surface plus m0, alpha_m, is_ordered,
%   moment_branch ('spontaneous' | 'none'), J0z_used. Early returns (PM
%   relaxation |m0| <= m_tol; MF non-convergence) carry the projected-parity
%   fixed set: m0, is_ordered=false, converged=false, Sigma0=NaN, crit=NaN,
%   si, tl=[], moment_branch. Always check pt.is_ordered AND pt.converged.
%
%   See also INVZT_SOLVE_POINT, INVZ_SOLVE_POINT_ORDERED (projected reference),
%   INVZ_SIGMA_ORDERED, INVZ_TWOLEVEL_ORDERED, INVZT_GCC_LATTICE,
%   INVZT_CRIT_STATIC, INVZT_SOLVE_AUTO.
if nargin < 5, opts = struct(); end
mode   = getf(opts, 'mode', 'a1');
% At Task-3 time only 'a1' exists; Task 7D extends this gate to accept 'a3d'
% (full-response dominant-vertex, Matsubara-only). 'a3' (full dress) stays
% rejected PERMANENTLY: the 136-state vertex is budget-refused (A4 ladder gate).
if ~ismember(char(mode), {'a1'})            % Task 7D: {'a1','a3d'}
    error('invzt:orderedMode', ['invzt_solve_point_ordered implements mode ''a1'' ' ...
        '(and ''a3d'' once Task 7D lands) -- got ''%s''. Full-dress ordered ''a3'' ' ...
        'is permanently out of scope: the 136-state vertex is budget-refused.'], char(mode));
end
Ecut   = getf(opts, 'Ecut', 40);
hyp    = getf(opts, 'hyp', true);
tmf    = getf(opts, 'transverse_mf', 'legacy_x');
mixo   = getf(opts, 'mix_outer', 0.7);
tolo   = getf(opts, 'tol_outer', 1e-8);
maxo   = getf(opts, 'max_outer', 80);
mtol   = getf(opts, 'm_tol', 1e-2);
bztol  = getf(opts, 'bz_tol', 1e-9);
Esplit = getf(opts, 'Esplit', 0.4653);
rank_tol = getf(opts, 'rank_tol', 1e-12);
chi_rest = ~isfield(opts, 'chi_rest') || isempty(opts.chi_rest) || ~isequal(opts.chi_rest, false);
odd = ~isfield(opts, 'odd') || isempty(opts.odd) || ~isequal(opts.odd, false);
nlevels = getf(opts, 'nlevels', 'std');
if ~strcmp(char(nlevels), 'std')
    error('invzt:orderedNlevels', ['invzt_solve_point_ordered supports nlevels = ' ...
        '''std'' only (full electronuclear single ion); got ''%s''.'], char(nlevels));
end

B = invz_field_vec(B);
if abs(B(3)) > bztol
    error('invzt:orderedLongitudinal', ['invzt_solve_point_ordered is the ' ...
        'TRANSVERSE (spontaneous-moment) route only: got Bz = %.3g T > bz_tol = ' ...
        '%.3g T. The forced-moment longitudinal route is not ported to the ' ...
        'tensor branch (2026-07-19 scope decision).'], B(3), bztol);
end
B(3) = 0;                                            % dead band: exactly transverse

J0z  = lat.info.Jcc0;                                % tensor-native uniform cc coupling
Jxx0 = lat.info.Jaa0;                                % v3 threading (parity with PM solver)

[wn, wts, beta] = invz_matsubara(T, Ecut);
nwn = numel(wn);

% --- ordered single-ion mean-field fixed point (full electronuclear space) -------
siopts = struct('hyp', hyp, 'order', true, 'J0z', J0z, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mz_seed', 'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), siopts.(f{1}) = opts.(f{1}); end
end
si = invz_single_ion(ion, T, B, siopts);
m0 = si.Jexp(3);
if ~si.mf_converged
    pt = early_return(m0, si, 'spontaneous');
    return;
end
if abs(m0) <= mtol
    pt = early_return(m0, si, 'none');               % paramagnetic point: use invzt_solve_point
    return;
end

% --- two-level driver at the converged ordering field hz; chi0 split -------------
tl = invz_twolevel_ordered(ion, T, B, si.hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
g  = real(invz_g(tl, 1i*wn));
[cdom, crest, mspec] = invzt_chi0_split(si, T, 1i*wn, struct('Esplit', Esplit));
cdom_cc  = real(squeeze(cdom(3,3,:)));
crest_cc = real(squeeze(crest(3,3,:)));
chi0cc0  = cdom_cc(1) + crest_cc(1);
if ~chi_rest
    crest    = zeros(size(crest));
    crest_cc = zeros(nwn, 1);
end

% --- odd = false mask on a COPY of lat.Jt (same rule as the PM solver) -----------
lat_eff = lat;
if ~odd
    lat_eff.Jt = invzt_odd_mask(lat.Jt);
end

% --- C2 assertion at Gamma (same as PM solver) -----------------------------------
JG = lat.JtGamma;
odd_ca = JG(3:3:12, 1:3:12);  odd_cb = JG(3:3:12, 2:3:12);
oddmag = max(abs([odd_ca(:); odd_cb(:)]));
if ~(oddmag < 1e-10*abs(lat.info.Jcc0))
    error('invzt:a1OddGamma', ['lat.JtGamma ODD (c<->a,b) blocks do not vanish ' ...
        '(max %.3e >= 1e-10*|Jcc0| = %.3e): C2 symmetry violated at Gamma.'], ...
        oddmag, 1e-10*abs(lat.info.Jcc0));
end

% --- outer moment-form loop: damped mixing, CHECK-BEFORE-MIX (P2-2) --------------
Sigma = zeros(nwn, 1);
denom = @(s) reshape(1 + s, 1, 1, nwn);
converged = false;
Gloc = nan(nwn, 1);  K = nan(nwn, 1);  diag4 = nan(4, nwn);
lam = nan(3, 1);  sg = struct('alpha', NaN, 'alpha_m', NaN, 'gamma', nan(nwn,1), 'Sigma', nan(nwn,1));
for outer = 1:maxo
    ctil = cdom ./ denom(Sigma) + crest;             % [3,3,nwn] renormalized chi0
    [Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);
    Gloc  = -Gcc(:);
    G0til = -(cdom_cc ./ (1 + Sigma) + crest_cc);
    K = 1 ./ Gloc - 1 ./ G0til;                      % LOCKED K bookkeeping (same as PM)
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);    % moment form needs lambda_3
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);   % HTML eqs 37-38
    dS  = max(abs(sg.Sigma - Sigma));
    if dS < tolo, converged = true; break; end       % exit BEFORE mixing: state consistent
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
end

ctil0 = cdom(:,:,1) / (1 + Sigma(1)) + crest(:,:,1);
[crit, crit_clipped_mass, crit_active_rank] = invzt_crit_static(ctil0, lat.JtGamma, rank_tol);

% --- assemble pt (invzt_solve_point surface + ordered extras) --------------------
pt.m0     = m0;
pt.is_ordered = true;
pt.moment_branch = 'spontaneous';
pt.Sigma0 = Sigma(1);
pt.Sigma  = Sigma;
pt.alpha  = sg.alpha;
pt.alpha_m = sg.alpha_m;
pt.lambda = lam;
pt.K = K;
pt.G = Gloc;
pt.tl = tl;
pt.si = si;
pt.chi0cc0 = chi0cc0;
pt.crit = crit;
pt.crit_clipped_mass = crit_clipped_mass;
pt.crit_active_rank  = crit_active_rank;
pt.sumrule_rel = abs(sum(wts.*Gloc)/beta + si.JzJz_fluct) / max(abs(si.JzJz_fluct), 1e-12);
phys_finite = all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(Gloc));
pt.converged = converged && phys_finite;
pt.outer_iters = outer;
pt.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
pt.mode = 'a1';
pt.odd  = odd;
pt.chi_rest = chi_rest;
pt.mspec = mspec;
pt.Jxx0_used = Jxx0;
pt.J0z_used  = J0z;
pt.nlevels = 'std';
pt.lat = struct('qvec_hash', hash_vec(lat.qvec(:)), 'conv', lat.conv, ...
    'JtGamma', lat.JtGamma, 'info', lat.info, 'Jxx0_used', Jxx0);
end

% ------------------------------- local helpers ----------------------------------
function pt = early_return(m0, si, branch)
%EARLY_RETURN Complete field set for every non-accepted exit (projected parity).
pt = struct('m0', m0, 'is_ordered', false, 'converged', false, 'Sigma0', NaN, ...
            'crit', NaN, 'si', si, 'tl', [], 'moment_branch', branch);
end

function h = hash_vec(v)
% Weak grid-identity hash, same formula as invz_cache_key / invzt_solve_point.
v = v(:);
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
```

- [ ] **Step 4: Run the CORE tests, verify pass**

Run: `$MB "cd('invz_tensor/tests'); r = runtests('test_invzt_solve_point_ordered'); disp(table(r)); assert(~any([r.Failed]))"`
Expected: all 8 pass. Record the printed `m0/Sigma0/crit` anchors in the commit body.

- [ ] **Step 5: Write the interop parity test** — `invz_tensor/tests/interop/test_invzt_solve_ordered_parity.m`:

```matlab
function tests = test_invzt_solve_ordered_parity
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', '..'));
addpath(fullfile(here, '..', '..', '..', 'invz_common'));
addpath(fullfile(here, '..', '..', '..', 'invz_projected'));
end

function test_ordered_no_odd_parity_live(testCase)
% Ordered-phase analog of the PM no-ODD live parity: both legs on the SAME grid
% and couplings (J0eff = info.Jcc0, Jxx0 = info.Jaa0 -> identical single-ion MF
% fixed point, hence identical m0 by construction). The tensor medium differs
% from the scalar one only by named residuals (cross-Cartesian chi0 elements,
% dominant-mask vs whole-cc division), so gate Sigma0/alpha_m at the PM
% frozen-baseline scale and m0 tight.
ion = invz_ion();  T = 0.1;  Bx = 3.0;
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
[Jnu, info] = invz_jq_modes(ion, g.qvec, struct('dpRng', 15, 'cache', true));
pj = invz_solve_point_ordered(ion, T, [Bx 0 0], Jnu(:), ...
    struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0));
pt = invzt_solve_point_ordered(ion, T, [Bx 0 0], lat, ...
    struct('odd', false, 'chi_rest', false));
verifyTrue(testCase, pj.is_ordered && pj.converged);
verifyTrue(testCase, pt.is_ordered && pt.converged);
verifyLessThan(testCase, abs(pt.m0 - pj.m0), 1e-6);
verifyLessThan(testCase, abs(pt.Sigma0 - pj.Sigma0), 5e-3);
verifyLessThan(testCase, abs(pt.alpha_m - pj.alpha_m), 5e-3);
fprintf(['INTEROP ordered no-ODD live: dm0 = %.3e, dSigma0 = %.3e, dalpha_m = %.3e, ' ...
    'dcrit = %.3e (m0 %.4f; Sigma0 tensor %.5f / proj %.5f)\n'], ...
    pt.m0 - pj.m0, pt.Sigma0 - pj.Sigma0, pt.alpha_m - pj.alpha_m, ...
    pt.crit - pj.crit, pt.m0, pt.Sigma0, pj.Sigma0);
end
```

Note (review P2-2 interaction): the projected ordered solver mixes-then-checks, this solver checks-then-mixes — a tolerance-sized (≤ tol_outer/mixo) systematic on converged exits, far below the 5e-3 gate.

- [ ] **Step 6: Run interop, verify pass**

Run: `$MB "cd('invz_tensor/tests/interop'); r = runtests('test_invzt_solve_ordered_parity'); disp(table(r)); assert(~any([r.Failed]))"`
Expected: pass. If `dSigma0` exceeds 5e-3: STOP — investigate which named residual grew (print `chi0cc0` both legs first); do not widen the gate.

- [ ] **Step 7: Commit**

```bash
git add invz_tensor/invzt_solve_point_ordered.m invz_tensor/tests/test_invzt_solve_point_ordered.m invz_tensor/tests/interop/test_invzt_solve_ordered_parity.m
git commit -m "feat(invzt): ordered-phase a1 point solver -- moment-form Sigma on the tensor lattice engine, projected-parity gated"
```

---

### Task 4: Ordered branch in `invzt_chi_realaxis.m`

**Files:**
- Modify: `invz_tensor/invzt_chi_realaxis.m` (continuation block, lines 149–159; header)
- Test: extend `invz_tensor/tests/test_invzt_chi_realaxis.m`

**Interfaces:**
- Consumes: `pt` from `invzt_solve_point_ordered` (`is_ordered`, `alpha_m`, `tl.m`, ordered `si`) or `invzt_solve_point` (PM, unchanged).
- Produces: same `out` contract, now with the moment-form continuation for ordered points:
  `Sw = (pt.alpha - pt.alpha_m) + (gamma_w - (2*tl.m^2/tl.M2)*gamma0) .* g`, `gamma_w = pref*(pt.lambda(1) - (1-n01^2)*Kw)`, `gamma0` at `K0 = pt.K(1)`. (With the frozen static `Kw = pt.K(1)` seeding, `gamma_w == gamma0` identically; both symbols are kept for parity with the projected `realaxis_sigma`, whose `npass`/`Jfull` re-solve makes them differ — a possible future port.)

**Review P1-1 is the driver of this task's test design:** the CURRENT code silently applies the PM formula to an ordered `pt` (it never touches `alpha_m`), so broad spectral-shape tests can pass against the WRONG formula. The gate below tests `out.Sigma_w` against the exact ordered expression AND asserts that expression differs materially from the PM one at the anchor — so the pre-change code MUST fail.

- [ ] **Step 1: Write the failing tests** — append to `test_invzt_chi_realaxis.m` (reuse the file's existing setup/ion/lat plumbing):

```matlab
function test_ordered_sigma_w_exact_formula(tc)
% THE ordered-continuation gate (review P1-1). out.Sigma_w must equal the
% moment-form expression assembled from the SAME pt fields; and that expression
% must differ MATERIALLY from the PM expression here, so the pre-change code
% (which applies the PM formula to an ordered pt without error) FAILS this test.
ion = tc.TestData.ion;  lat = tc.TestData.lat;  T = 0.1;
pt = invzt_solve_point_ordered(ion, T, [3.5 0 0], lat, struct());
verifyTrue(tc, pt.is_ordered && pt.converged);        % anchor MUST hold (not assume)
w = [0; 0.10; 0.31; 0.45];  eta = 2e-3;  z = w + 1i*eta;
o = invzt_chi_realaxis(ion, T, [3.5 0 0], pt, w, struct('qsel', 'gamma_uniform', 'eta', eta));
tl = pt.tl;  g = invz_g(tl, z);  pref = tl.M2/tl.n01^2;
gamma0 = pref*(pt.lambda(1) - (1 - tl.n01^2)*pt.K(1));     % frozen Kw: gamma_w == gamma0
Sw_ord = (pt.alpha - pt.alpha_m) + gamma0*(1 - 2*tl.m^2/tl.M2) .* g;
Sw_pm  = pt.alpha + gamma0 .* g;
verifyGreaterThan(tc, max(abs(Sw_ord - Sw_pm)), 1e-4);     % formulas differ materially here
verifyEqual(tc, o.Sigma_w, Sw_ord, 'AbsTol', 1e-12, 'RelTol', 1e-10);   % exact algebra
fprintf('ordered Sigma_w gate: max|ord-pm| = %.4g, m=%.4f, alpha_m=%.4g\n', ...
    max(abs(Sw_ord - Sw_pm)), tl.m, pt.alpha_m);
end

function test_ordered_point_spectra(tc)
% Broad ordered-spectrum sanity: finite, non-negative up to the frozen-Kw
% caveat, soft mode inside (0.05, 0.6) meV.
ion = tc.TestData.ion;  lat = tc.TestData.lat;  T = 0.1;
pt = invzt_solve_point_ordered(ion, T, [3.5 0 0], lat, struct());
verifyTrue(tc, pt.is_ordered && pt.converged);
w = linspace(0, 0.6, 601).';
o = invzt_chi_realaxis(ion, T, [3.5 0 0], pt, w, struct('qsel', 'gamma_uniform', 'eta', 2e-3));
c = squeeze(imag(o.chi_uniform(3,3,:)));
verifyTrue(tc, all(isfinite(c)));
Ep = invz_peak_energy(c, w, 0.05);
verifyTrue(tc, isfinite(Ep) && Ep > 0.05 && Ep < 0.6);
verifyGreaterThan(tc, min(c), -0.05*max(c));
fprintf('ordered realaxis @ 3.5T: Epeak=%.4f meV, max chi''''=%.4g, min=%.3g\n', Ep, max(c), min(c));
end

function test_ordered_mode_softens_toward_Bc(tc)
% FM-side soft-mode direction: the mode SOFTENS approaching Bc (4.65-4.70 T)
% from below -- E(3.5 T) > E(4.4 T) > 0.
ion = tc.TestData.ion;  lat = tc.TestData.lat;  T = 0.1;
w = linspace(0, 0.6, 601).';
E = nan(1, 2);  Bs = [3.5 4.4];
for k = 1:2
    pt = invzt_solve_point_ordered(ion, T, [Bs(k) 0 0], lat, struct());
    verifyTrue(tc, pt.is_ordered && pt.converged);
    o = invzt_chi_realaxis(ion, T, [Bs(k) 0 0], pt, w, struct('qsel', 'gamma_uniform', 'eta', 2e-3));
    E(k) = invz_peak_energy(squeeze(imag(o.chi_uniform(3,3,:))), w, 0.05);
end
verifyGreaterThan(tc, E(1), E(2));
fprintf('FM soft mode: E(3.5T)=%.4f > E(4.4T)=%.4f meV\n', E(1), E(2));
end

function test_ordered_force_sigma0_bare_rpa(tc)
% BRANCH-INDEPENDENT regression only (review P1-1: force_sigma0 bypasses BOTH
% formulas, so this cannot gate the ordered continuation): bare RPA of the
% ORDERED chi0 stays non-negative.
ion = tc.TestData.ion;  lat = tc.TestData.lat;  T = 0.1;
pt = invzt_solve_point_ordered(ion, T, [3.5 0 0], lat, struct());
verifyTrue(tc, pt.is_ordered && pt.converged);
w = linspace(0, 0.6, 301).';
o = invzt_chi_realaxis(ion, T, [3.5 0 0], pt, w, ...
    struct('qsel', 'gamma_uniform', 'eta', 2e-3, 'force_sigma0', true));
c = squeeze(imag(o.chi_uniform(3,3,:)));
verifyGreaterThan(tc, min(c), -1e-10*max(abs(c)));
end
```

- [ ] **Step 2: Run to verify failure — the exact-formula gate must fail pre-change**

Run: `$MB "cd('invz_tensor/tests'); r = runtests('test_invzt_chi_realaxis'); disp(table(r))"`
Expected: `test_ordered_sigma_w_exact_formula` FAILS on the `verifyEqual(o.Sigma_w, Sw_ord, ...)` line (the current code returns the PM expression, which the test just proved differs by > 1e-4). The other new tests may pass or fail — only the formula gate is the load-bearing pre-change failure.

- [ ] **Step 3: Implement the ordered branch** — replace the continuation block (lines 149–159):

```matlab
% --- A1 scalar-Sigma continuation (Kw-seeding pattern; PM = projected non-ordered
%     branch, ordered = projected HTML-eq-37 moment form, both frozen static K) ---
ordered = isfield(pt, 'is_ordered') && pt.is_ordered;
tl   = pt.tl;
g    = invz_g(tl, z);
pref = tl.M2 / tl.n01^2;
Kw   = pt.K(1) * ones(nw, 1);
if force_sigma0
    Sw = zeros(nw, 1);
elseif ordered
    gamma_w = pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw);
    gamma0  = pref*(pt.lambda(1) - (1 - tl.n01^2)*pt.K(1));
    Sw = (pt.alpha - pt.alpha_m) + (gamma_w - (2*tl.m^2/tl.M2)*gamma0) .* g;
else
    Sw = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw) .* g;
end
```

ONE more code change (2026-07-20 amendment, matching Task 3's whole-cc medium): the ordered branch's ctil rebuild must be WHOLE-CC too — where the PM branch keeps its `invzt_chi0_split`-based `cdom/(1+Sw) + crest` rebuild unchanged, the ordered branch instead rebuilds

```matlab
c0w  = invz_chi0z(pt.si, T, z, struct('elastic', false));   % full local chi0, elastic OFF
ctil = c0w ./ reshape(1 + Sw, 1, 1, nw);
```

(the exact tensor analog of the projected ordered `chit = -G0/(1+Sw)`; elastic OFF on the real axis per the global convention). qsel/odd-mask/RPA/projection stay branch-independent. Header: SCOPE box → "A1 scalar-Sigma continuation, PM AND ordered branches; A2/A3 continuation remains the open item"; delete the line-6 claim that no ordered points exist; document the ordered whole-cc rebuild.

- [ ] **Step 4: Full CORE — zero failures**

Run: `$MB "cd('invz_tensor/tests'); r = runtests(pwd); disp(table(r)); assert(~any([r.Failed]))"`
Expected: zero failures; the PM `Sw` expression is verbatim unchanged, so all pre-existing realaxis tests hold.

- [ ] **Step 5: Commit**

```bash
git add invz_tensor/invzt_chi_realaxis.m invz_tensor/tests/test_invzt_chi_realaxis.m
git commit -m "feat(invzt): ordered-phase moment-form real-axis continuation, exact-formula gated"
```

---

### Task 5: `invzt_solve_auto.m` — STABILITY-BASED phase dispatcher (review P0-1 redesign)

**Files:**
- Create: `invz_tensor/invzt_solve_auto.m`
- Test: `invz_tensor/tests/test_invzt_solve_auto.m`

**Interfaces:**
- Consumes: `invzt_solve_point` (existing), `invzt_solve_point_ordered` (Task 3).
- Produces: `[pt, phase, di] = invzt_solve_auto(ion, T, B, lat, opts)` — `phase` 1 = ordered accepted, 2 = PM accepted, 0 = neither. `di` carries STRUCTURED per-leg diagnostics (review P2-1): `di.para` and `di.ordered`, each `struct('attempted',L,'accepted',L,'converged',L,'m0',D,'crit',D,'Sigma0',D,'err',S)` (`m0` = NaN on the PM leg). Consumed by Task 6.

**Selection rule (P0-1):** the bare-MF ordered leg keeps `m0 > m_tol` up to ~5.0 T while the tensor QPT sits at 4.65–4.70 T (measured, see Global Constraints) — so ordered-first selection would misassign the whole [Bc, ~5.0 T] band to the FM phase. The PM leg's three-part validity rule (`converged && crit > 0 && Sigma0 >= sigma_floor`) IS the tensor QPT criterion, so it decides FIRST; the ordered leg is consulted only when the PM sample is invalid, and must itself pass `converged && si.mf_converged && isfinite(Sigma0) && crit > crit_tol_ordered`. The ordered stability criterion `crit > crit_tol_ordered` (default −1e-3, the Task-3 deep-ordered gate scale) is the explicit justification the review demanded: `invzt_crit_static` at the ordered renormalized `ctil0` measures local stability of the broken-symmetry state; a converged ordered point failing it is masked (phase 0) with both diagnostics preserved, never labeled a physical phase.

**Input policy (round-2 P1-5):** validation happens AT DISPATCHER ENTRY, before either leg: the field is normalized and `|Bz| > bz_tol` raises `invzt:orderedLongitudinal` immediately (round-1 draft applied the guard only inside the ordered leg, so a PM-valid point silently accepted a longitudinal field), and `opts.mode` other than `'a1'` raises `invzt:autoMode` (`a3d` is never used to classify the phase — it is evaluated AFTER classification, on accepted ordered points only). The catch policy is an ENUMERATED ALLOWLIST of recoverable physics outcomes — `invz:degenerateDoublet`, `invzt:a1ZeroField`, `invz:orderedPhase` — everything else (config, shape, unsupported-mode, invariant violations) rethrows; the round-1 blanket `invz:*`/`invzt:*` absorption would have converted malformed options and lattice errors into phase 0.

- [ ] **Step 1: Write the failing tests** — `invz_tensor/tests/test_invzt_solve_auto.m`:

```matlab
function tests = test_invzt_solve_auto
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
ion = invz_ion();
g   = invzt_qgrid(16, 'halfopen');
tc.TestData.ion = ion;
tc.TestData.lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
end

function test_below_Bc_is_phase1(tc)
[pt, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0], tc.TestData.lat, struct());
verifyEqual(tc, phase, 1);
verifyTrue(tc, pt.is_ordered && pt.converged);
verifyTrue(tc, di.para.attempted && ~di.para.accepted);   % PM leg was consulted and rejected
end

function test_at_4p8_is_phase2_HARD(tc)
% THE P0-1 gate: 4.8 T is on the tensor-PM side (crit = +0.023 measured) even
% though the bare-MF ordered leg still has m0 = 1.17 there. Stability-based
% selection MUST return PM. An ordered-first implementation fails here.
[pt, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [4.8 0 0], tc.TestData.lat, struct());
verifyEqual(tc, phase, 2);
verifyTrue(tc, pt.converged && pt.crit > 0);
verifyFalse(tc, di.ordered.attempted);    % PM valid -> ordered leg never consulted
end

function test_above_Bc_is_phase2(tc)
[pt, phase] = invzt_solve_auto(tc.TestData.ion, 0.1, [5.5 0 0], tc.TestData.lat, struct());
verifyEqual(tc, phase, 2);
verifyTrue(tc, pt.converged && pt.crit > 0);
end

function test_spurious_pm_rejected_deterministic(tc)
% Review P1-2: force the ordered leg to early-return DETERMINISTICALLY
% (m_tol = Inf -> |m0| <= m_tol always -> is_ordered = false), leaving the PM
% leg untouched. At 4.0 T the PM leg CONVERGES with crit = -0.186 (measured):
% an implementation missing the crit > 0 rule returns phase 2 here and fails.
[~, phase, di] = invzt_solve_auto(tc.TestData.ion, 0.1, [4.0 0 0], tc.TestData.lat, ...
    struct('m_tol', Inf));
verifyEqual(tc, phase, 0);
verifyTrue(tc, di.para.attempted && ~di.para.accepted);
verifyTrue(tc, di.para.converged);                   % the spurious point DID converge...
verifyLessThan(tc, di.para.crit, 0);                 % ...and was rejected on stability
verifyTrue(tc, di.ordered.attempted && ~di.ordered.accepted);
end

function test_boundary_bracket(tc)
% Phase boundary consistency (P0-1): scanning across the measured crit bracket,
% no PM label at or below 4.65 T, PM at and above 4.75 T, and no re-entrant
% ordered label above the first PM field.
ion = tc.TestData.ion;  lat = tc.TestData.lat;
Bs = 4.55:0.05:4.85;  ph = zeros(size(Bs));
for k = 1:numel(Bs)
    [~, ph(k)] = invzt_solve_auto(ion, 0.1, [Bs(k) 0 0], lat, struct());
end
fprintf('boundary scan: B = [%s], phase = [%s]\n', num2str(Bs, '%.2f '), num2str(ph, '%d '));
verifyTrue(tc, ~any(ph(Bs <= 4.651) == 2));          % no PM at/below the lower bracket edge
verifyTrue(tc, all(ph(Bs >= 4.749) == 2));           % PM at/above the upper bracket edge
i2 = find(ph == 2, 1);
verifyTrue(tc, isempty(i2) || ~any(ph(i2:end) == 1));  % monotone: no FM re-entry
end

function test_longitudinal_rejected_at_entry(tc)
% Round-2 P1-5: the Bz guard is at DISPATCHER ENTRY -- it must fire even where
% the PM leg would be valid (5.5 T), not only on the ordered side. The round-1
% draft guarded only inside the ordered leg, so a PM-valid point silently
% accepted a longitudinal field.
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [5.5 0 1e-6], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0.1], ...
    tc.TestData.lat, struct()), 'invzt:orderedLongitudinal');
end

function test_invalid_mode_rejected(tc)
% a3d/a3 never classify the phase; a non-a1 mode is a caller error, not phase 0.
verifyError(tc, @() invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0], ...
    tc.TestData.lat, struct('mode', 'a3d')), 'invzt:autoMode');
end

function test_config_errors_rethrown_not_masked(tc)
% Round-2 P1-5: the catch allowlist must NOT convert configuration errors into
% phase 0. A malformed lattice struct raises a non-allowlisted identifier,
% which must propagate out of the dispatcher.
badlat = struct('nothing', 1);           % no Jt/JtGamma/info fields
try
    invzt_solve_auto(tc.TestData.ion, 0.1, [3.0 0 0], badlat, struct());
    verifyTrue(tc, false, 'malformed lattice must throw, not return phase 0');
catch err
    verifyFalse(tc, ismember(err.identifier, ...
        {'invz:degenerateDoublet', 'invzt:a1ZeroField', 'invz:orderedPhase'}));
end
end
```

- [ ] **Step 2: Run to verify failure**

Run: `$MB "cd('invz_tensor/tests'); r = runtests('test_invzt_solve_auto'); disp(table(r))"`
Expected: FAIL — `invzt_solve_auto` undefined.

- [ ] **Step 3: Implement `invz_tensor/invzt_solve_auto.m`**

```matlab
function [pt, phase, di] = invzt_solve_auto(ion, T, B, lat, opts)
%INVZT_SOLVE_AUTO Stability-based tensor a1 phase dispatcher at one (T, B).
%   [pt, phase, di] = INVZT_SOLVE_AUTO(ion, T, B, lat, opts) assigns the phase
%   by STABILITY, not by moment onset (review P0-1, 2026-07-19): the bare-MF
%   ordered leg keeps m0 > 0 well onto the PM side (measured 0.1 K: m0 = 1.17
%   at 4.8 T vs the tensor crit = 0 QPT at 4.65-4.70 T), so the projected
%   ordered-first pattern would misassign [Bc, ~5.0 T]. Here the PM leg decides
%   FIRST via the three-part validity rule -- its crit > 0 IS the tensor QPT
%   criterion -- and the ordered leg is consulted only when the PM sample is
%   invalid:
%     phase 2 : PM   accepted -- converged, finite crit > 0, Sigma0 >= sigma_floor.
%     phase 1 : PM invalid AND ordered accepted -- is_ordered, converged,
%               si.mf_converged, finite Sigma0, crit > crit_tol_ordered
%               (default -1e-3): INVZT_CRIT_STATIC at the ordered renormalized
%               ctil0 measures local stability of the broken-symmetry state; a
%               converged ordered point failing it is NOT labeled a phase.
%     phase 0 : neither accepted -- near-Bc Option-A window or solver failure;
%               pt carries the last usable attempt ([] if none), di says why.
%
%   TRANSVERSE ONLY: invzt:orderedLongitudinal is RETHROWN, not absorbed (no
%   tensor forced-moment route -- caller error, not a phase outcome).
%
%   di (review P2-1): STRUCTURED per-leg diagnostics -- di.para / di.ordered,
%   each with attempted, accepted, converged, m0 (NaN for PM), crit, Sigma0,
%   err (exception identifier, '' if returned normally). A leg that returned
%   but was rejected is fully described -- rejection reasons are read off the
%   numbers (e.g. converged = true, crit < 0 = the spurious below-Bc PM point).
%
%   INPUT POLICY (validated AT ENTRY, before either leg -- round-2 P1-5):
%     - B is normalized (invz_field_vec); |Bz| > bz_tol raises
%       invzt:orderedLongitudinal HERE, so a PM-valid point can never silently
%       accept a longitudinal field.
%     - opts.mode must be 'a1' (invzt:autoMode otherwise): the phase is NEVER
%       classified by a3d -- evaluate a3d separately on accepted ordered points.
%
%   Error policy: ENUMERATED ALLOWLIST -- only the recoverable physics outcomes
%   invz:degenerateDoublet, invzt:a1ZeroField, invz:orderedPhase are absorbed
%   into .err; every other identifier (config, shape, unsupported mode,
%   invariant violations) RETHROWS.
%
%   See also INVZT_SOLVE_POINT, INVZT_SOLVE_POINT_ORDERED, INVZ_SOLVE_AUTO
%   (projected reference -- NOTE its ordered-first rule is NOT ported, P0-1).
if nargin < 5, opts = struct(); end
sfloor = getf(opts, 'sigma_floor', -0.5);            % single-sourced validity floor
ctolo  = getf(opts, 'crit_tol_ordered', -1e-3);      % ordered stability tolerance
bztol  = getf(opts, 'bz_tol', 1e-9);

% --- entry validation (round-2 P1-5): field and mode, BEFORE any solve ----------
B = invz_field_vec(B);
if abs(B(3)) > bztol
    error('invzt:orderedLongitudinal', ['invzt_solve_auto is the TRANSVERSE ' ...
        '(spontaneous-moment) route only: got Bz = %.3g T > bz_tol = %.3g T. ' ...
        'No tensor forced-moment route exists (2026-07-19 scope decision).'], B(3), bztol);
end
B(3) = 0;                                            % dead band: exactly transverse
if ~strcmp(char(getf(opts, 'mode', 'a1')), 'a1')
    error('invzt:autoMode', ['invzt_solve_auto classifies the phase with mode ' ...
        '''a1'' ONLY (got ''%s''): a3d is an after-classification evaluation on ' ...
        'accepted ordered points, never the phase criterion.'], ...
        char(getf(opts, 'mode', 'a1')));
end
RECOVERABLE = {'invz:degenerateDoublet', 'invzt:a1ZeroField', 'invz:orderedPhase'};
leg0 = struct('attempted', false, 'accepted', false, 'converged', false, ...
              'm0', NaN, 'crit', NaN, 'Sigma0', NaN, 'err', '');
di = struct('para', leg0, 'ordered', leg0);
pt = [];  phase = 0;

% --- PM leg first: its crit > 0 is the QPT criterion ----------------------------
di.para.attempted = true;
try
    ptp = invzt_solve_point(ion, T, B, lat, opts);
    di.para.converged = ptp.converged;
    di.para.crit = ptp.crit;  di.para.Sigma0 = ptp.Sigma0;
    if ptp.converged && isfinite(ptp.crit) && ptp.crit > 0 && ptp.Sigma0 >= sfloor
        di.para.accepted = true;
        pt = ptp;  phase = 2;  return;
    end
    pt = ptp;                                        % keep for diagnostics
catch err
    if ~ismember(err.identifier, RECOVERABLE), rethrow(err); end
    di.para.err = err.identifier;
end

% --- ordered leg: only when the PM sample is invalid ----------------------------
di.ordered.attempted = true;
try
    pto = invzt_solve_point_ordered(ion, T, B, lat, opts);
    di.ordered.converged = pto.converged;
    di.ordered.m0 = pto.m0;  di.ordered.crit = pto.crit;  di.ordered.Sigma0 = pto.Sigma0;
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0) ...
            && pto.si.mf_converged && pto.crit > ctolo
        di.ordered.accepted = true;
        pt = pto;  phase = 1;  return;
    end
    if ~isempty(pto.si), pt = pto; end               % last usable attempt wins the pt slot
catch err
    if ~ismember(err.identifier, RECOVERABLE), rethrow(err); end
    di.ordered.err = err.identifier;
end
end
```

Cost note (conscious, not hidden): deep-ordered fields pay a failed PM attempt first (up to `max_outer` iterations). At the driver's 101 fields under parfor this is acceptable; if deep-ordered PM attempts ever dominate wall-clock, add a coarse crit-based pre-screen in a LATER change — do not reorder the legs, correctness owns the order.

- [ ] **Step 4: Run tests, verify pass**

Run: `$MB "cd('invz_tensor/tests'); r = runtests('test_invzt_solve_auto'); disp(table(r)); assert(~any([r.Failed]))"`
Expected: all 8 test functions pass (review minor: count the functions, don't trust prose — currently 8 after the P1-5 additions; the boundary-bracket scan is ~9 two-leg solves — minutes, not hours; if it exceeds ~10 min, gate it behind `INVZ_SLOW` and add its name to the CORE Incomplete allowlist in Task 8).

- [ ] **Step 5: Commit**

```bash
git add invz_tensor/invzt_solve_auto.m invz_tensor/tests/test_invzt_solve_auto.m
git commit -m "feat(invzt): stability-based phase dispatcher -- PM crit decides first; bare-MF moment never assigns the phase (review P0-1)"
```

---

### Task 6: Rewire `invzt_run_spectra.m` across the QPT (+ knob fixes from the 2026-07-19 diagnosis)

**Files:**
- Modify: `invz_tensor/invzt_run_spectra.m`

**Interfaces:**
- Consumes: `invzt_solve_auto` (Task 5), `invzt_chi_realaxis` (Task 4), `invz_peak_energy`.
- Produces: driver output only (figures + workspace `chipp`, `phasev`, `critv`, `Epeak`).

Diagnosis being fixed (measured 2026-07-19): (i) `w` topped out at 0.025 meV while the soft mode sits at 0.26–0.34 meV — 97% of the spectral weight above the window; (ii) `eta = 2e-5` < step 1.25e-4 → aliasing (10× peak-height error); (iii) fields 4.0–4.65 T masked by the PM-only solver. (i)+(ii) came from porting the projected GHz-unit window as meV.

- [ ] **Step 1: Replace the knob block** (current lines 26–42):

```matlab
% ---- knobs (mirroring invz_run_spectra's names where the concept carries over) --
T = 0.1;                            % K
fields = linspace(3.0, 7.0, 101);   % T -- SPANS the QPT (crit = 0 bracket 4.65-4.70 T at
                                    % 0.1 K): stability-based auto solve assigns FM below,
                                    % PM above (invzt_solve_auto; NOT the bare-MF moment,
                                    % which persists to ~5.0 T -- review P0-1)
w = linspace(0, 0.6, 601).';        % meV, uniform (invz_peak_energy assumes it). Window must
                                    % contain the soft mode: 0.26-0.34 meV on the PM side at
                                    % 4.8-6 T (measured 2026-07-19) -- the old 0-0.025 meV
                                    % port of the projected GHz window cut 97% of the weight.
eta = 2e-3;                         % meV, real-axis Lorentzian HWHM. MUST stay >= the w step
                                    % (1e-3 here) or peaks alias between grid points (measured
                                    % 10x peak-height error at eta/step = 0.16); guarded below.
sliceMax = 6;                       % field count at/below which line slices are drawn
peak_wmin = 0.05;                   % meV -- excludes the low-frequency hyperfine line so
                                    % Epeak tracks the SOFT MODE (the 0.003-0.009 meV
                                    % hyperfine feature dominates the raw max below ~5 T).
theta_c = 0.0;  phi_ab = 0.0;       % tilt knobs (deg). theta_c ~= 0 gives Bz ~= 0 ->
                                    % invzt:orderedLongitudinal (no tensor forced-moment
                                    % route; 2026-07-19 scope).
transverse_mf = 'legacy_x';         % 'legacy_x' | 'none' | 'vector_ab'
gridN = 16; gridConv = 'halfopen'; dpRng = 30;
useParallel = true;
solve_opts = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', ...
                    'transverse_mf', transverse_mf);
```

- [ ] **Step 2: Add the eta/step guard** after the knob block (and its `wq` twin in the q-path branch):

```matlab
if eta < (w(2) - w(1))
    error('invzt_run_spectra:etaStep', ['eta = %.3g meV < w step = %.3g meV: a ' ...
        'Lorentzian narrower than the sampling aliases between grid points ' ...
        '(measured 10x peak-height error at eta/step = 0.16, 2026-07-19). Raise ' ...
        'eta or refine w.'], eta, w(2) - w(1));
end
```

- [ ] **Step 3: Replace the field-sweep parfor body** (current lines 111–152) — two-phase solve with structured-diagnostic logging (review P2-1):

```matlab
    % ---------------- field sweep at the uniform mode, ACROSS Bc ---------------
    nWorkers = 0;
    if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end
    nf = numel(fields);
    chipp  = nan(numel(w), nf);
    phasev = zeros(1, nf);          % 1 = ordered, 2 = PM, 0 = masked
    critv  = nan(1, nf);            % branch stability of the ACCEPTED leg (review re-plan #4)
    parfor (k = 1:nf, nWorkers)
        col = nan(numel(w), 1);  ph = 0;  ck = NaN;
        try
            [pt, ph, di] = invzt_solve_auto(ion, T, fields(k)*dhat, lat, solve_opts);
            if ph > 0
                ck = pt.crit;
                o = invzt_chi_realaxis(ion, T, fields(k)*dhat, pt, w, ...
                        struct('qsel', 'gamma_uniform', 'eta', eta));
                col = squeeze(imag(o.chi_uniform(3, 3, :)));
            else
                % Structured per-leg outcomes (P2-1): a leg that RETURNED but was
                % rejected is described by its numbers, not a blank error string.
                fprintf(['  B = %.2f T: masked -- PM(att=%d conv=%d crit=%.4g err=%s) ' ...
                         'ORD(att=%d conv=%d m0=%.4g crit=%.4g err=%s)\n'], fields(k), ...
                    di.para.attempted, di.para.converged, di.para.crit, di.para.err, ...
                    di.ordered.attempted, di.ordered.converged, di.ordered.m0, ...
                    di.ordered.crit, di.ordered.err);
            end
        catch err
            switch err.identifier
                case {'invz:degenerateDoublet', 'invzt:a1ZeroField'}
                    fprintf('  B = %.2f T: %s (masked)\n', fields(k), err.identifier);
                otherwise
                    rethrow(err);
            end
        end
        chipp(:, k) = col;  phasev(k) = ph;  critv(k) = ck;
    end
    fprintf('phase summary: %d ordered / %d PM / %d masked of %d fields\n', ...
        sum(phasev == 1), sum(phasev == 2), sum(phasev == 0), nf);
```

Keep the plotting blocks; add to the colormap title `across B_c`; after the `Epeak` plot add `if any(phasev == 1) && any(phasev == 2), xline(fields(find(phasev == 2, 1)), '--', 'B_c'); end`.

- [ ] **Step 4: Update the q-path branch** — replace the single PM solve + `invzt:qpathNotPM` check (lines 93–100):

```matlab
    [pt, phaseq, diq] = invzt_solve_auto(ion, T, Bq*dhat, lat, solve_opts);
    if phaseq == 0
        error('invzt:qpathInvalid', ['q-path point B = %.2f T, T = %.2f K failed both ' ...
            'legs [PM: conv=%d crit=%.4g err=%s | ORD: conv=%d m0=%.4g crit=%.4g err=%s]: ' ...
            'the whole q-path product hinges on this one point, so failing loudly beats ' ...
            'drawing an empty map.'], Bq, T, diq.para.converged, diq.para.crit, ...
            diq.para.err, diq.ordered.converged, diq.ordered.m0, diq.ordered.crit, ...
            diq.ordered.err);
    end
```

and append `(phase %d)` provenance to the q-path figure titles. Replace the header SCOPE box (currently "PARAMAGNETIC SIDE ONLY") with the two-phase description naming `invzt_solve_auto`, the stability-based rule, and the transverse-only constraint.

- [ ] **Step 5: Acceptance run — ONLY meaningful after Task 5's boundary tests are green (review re-plan #4)**

Run: `$MB "run('invz_tensor/invzt_run_spectra.m'); saveas(figure(1), 'Data/invzt_spectra_acrossBc_T0p1.fig'); exit"`
Expected: phase summary with ordered fields below the 4.65–4.70 bracket, PM above, masked band (phase 0) confined to the near-Bc Option-A window; soft mode softening into Bc from BOTH sides and re-hardening; `Epeak` V-shaped around the `xline`. REPORT the masked-band width in the commit body — it is the measured Option-A boundary gap, a finding, not a defect.

- [ ] **Step 6: Full CORE — zero failures**

Run: `$MB "cd('invz_tensor/tests'); r = runtests(pwd); disp(table(r)); assert(~any([r.Failed]))"`

- [ ] **Step 7: Commit**

```bash
git add invz_tensor/invzt_run_spectra.m Data/invzt_spectra_acrossBc_T0p1.fig
git commit -m "feat(invzt): spectra driver spans the QPT -- stability-based auto solve + window/eta/peak_wmin fixes from the 2026-07-19 diagnosis"
```

---

### Task 7: Ordered `a3d` — full-response, dominant-vertex, Matsubara-only (staged 7A→7E)

**(Round-2 scope: the dominant-basis route is selected. Full-dress 136-state ordered a3 stays permanently out of scope. Every stage below is its own test/commit cycle; a stage that fails its gate STOPs the task — record in ODD-LOG, do not proceed on a broken foundation.)**

**Honesty rule for this task:** the plan author has read `invzt_sigma_tensor.m` only in excerpts (the dominant branch, lines 122–145, and the output assembly) and has NOT read `invzt_gamma4.m`/`invzt_vertex4.m` internals or `invzt_rung_basis.m` in full. Every step that touches those files starts with "read the file"; the specs below fix the CONTRACTS (shapes, identifiers, gates), not the internals.

---

#### 7A — Independent prerequisites: `eps_el` fix + ordered oracle rows + honest Jensen diagnostic

**Files:**
- Modify: `invz_tensor/invzt_solve_point.m` (a3 branch, `eps_el` at line ~259)
- Create: `invz_tensor/tests/test_invzt_eps_el.m` (self-contained — no dependence on unread helpers)
- Modify: `verify_tensor_vertex.py` (TWO scalar-beta ordered systems; round-2 P1-6 kills the `"beta": [14,3]` sketch)
- Regenerate: `invz_tensor/tests/fixtures/vertex_oracle.json`
- Modify: `invz_tensor/tests/test_invzt_vertex.m` (ordered V-row test; ordered Jensen STATIC gate + finite-frequency REPORT; `verifyNotEmpty` row-count protection added to the existing PM V-row test too)

- [ ] **7A Step 1: `eps_el` failing test** — `invz_tensor/tests/test_invzt_eps_el.m`, self-contained:

```matlab
function tests = test_invzt_eps_el
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_eps_el_is_elastic_only(tc)
% Constraint 7: eps_el must use the ELASTIC-ONLY static cc weight c_d, not the
% full equal-time variance JzJz_fluct (review P1-4: the old formula was an
% upper bound). The three-state toy has zero diagonal Mz in the doublet ->
% c_d ~ 0 while JzJz_fluct = O(M^2) > 0: the two definitions differ by orders
% of magnitude, so the old formula CANNOT pass this test.
ion = invz_ion();  T = 2.0;  Bx = 0.5;
g   = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
pt  = invzt_solve_point(ion, T, [Bx 0 0], lat, struct('mode', 'a3', 'nlevels', 'three'));
verifyTrue(tc, pt.converged);
[wn, ~, beta] = invz_matsubara(T, 40);
c0_el   = invz_chi0z(pt.si, T, 1i*wn(1), struct('elastic', true));
c0_inel = invz_chi0z(pt.si, T, 1i*wn(1), struct('elastic', false));
c_d_ref = real(c0_el(3,3,1) - c0_inel(3,3,1)) / beta;
verifyEqual(tc, pt.c_d, c_d_ref, 'AbsTol', 1e-12, 'RelTol', 1e-10);
verifyEqual(tc, pt.eps_el, beta*abs(pt.K(1))*c_d_ref, 'AbsTol', 1e-12, 'RelTol', 1e-10);
verifyLessThan(tc, pt.c_d, 0.5*pt.si.JzJz_fluct);   % materially below the old upper bound
fprintf('eps_el fix: c_d=%.4g vs JzJz_fluct=%.4g (old formula overcounted %.3gx)\n', ...
    pt.c_d, pt.si.JzJz_fluct, pt.si.JzJz_fluct/max(abs(pt.c_d), 1e-300));
end
```

Run to verify failure (`pt.c_d` undefined), then implement in `invzt_solve_point.m` — replace the line-259 formula:

```matlab
    % constraint 7: c_d = ELASTIC-ONLY static cc weight (chi0's z~0 elastic term
    % enters as beta*c_d, read off the elastic-on/off difference of the SAME
    % chi0 convention the solve consumes). Review P1-4 (2026-07-19): the
    % previous formula used the FULL equal-time variance si.JzJz_fluct -- an
    % upper bound; prior logged eps_el values are upper bounds (ODD-LOG, Task 8).
    c0_el   = invz_chi0z(si, T, 1i*wn(1), struct('elastic', true));
    c0_inel = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
    c_d     = real(c0_el(3,3,1) - c0_inel(3,3,1)) / beta;
    eps_el  = beta * abs(K(1)) * c_d;
```

Add `pt.c_d = c_d;` next to `pt.eps_el` (NaN on a1/a2 branches, alongside the other a3-only fields). Re-run: pass. CORE zero-failures. Commit (path-scoped): `fix(invzt): eps_el uses elastic-only c_d (constraint 7) -- prior values were upper bounds`.

- [ ] **7A Step 2: ordered oracle systems (python).** READ the `jensen_2lvl` emission block (`verify_tensor_vertex.py` ~lines 700–724) first, then add TWO systems with SCALAR beta (round-2 P1-6 — `sysdata` expects scalar `S.beta`), reusing the block's own operator/row machinery verbatim with the operator swapped:

```python
# Ordered two-level (split doublet, diagonal moment +-m, off-diagonal M): the
# minimal system whose vertex carries the moment/elastic paths of J 2.26-2.29.
# TWO scalar-beta systems (sysdata expects scalar beta):
#   _b14: beta*Delta = 14 -> 1 - n01^2 = 3.3261e-6 (elastic sector ~dead)
#   _b3 : beta*Delta = 3  -> elastic sector alive (REPORT rows for the
#         elastic-resummation ambiguity, constraint 7)
# Operator centered by the THERMAL MEAN exactly as jensen_2lvl centers its
# operator; the system dict records m as the UNcentered diagonal (+m on the
# ground state), matching tl.m = ground-state diagonal element.
for tag_b, beta_o in (("b14", mpf(14)), ("b3", mpf(3))):
    # ... transcribe the jensen_2lvl block with:
    #     E = [0, Delta_o], operator [[+m_o, M_o], [M_o, -m_o]] (then centered),
    #     Delta_o = 1.0, M_o = 0.9, m_o = 0.6, same kernel K(l) and Lmax,
    #     rows tagged f"V;cc;jensen2lvlord_{tag_b};subtract={conv}",
    #     plus F/Gamma rows with the standard 4-slot tags (auto-consumed by
    #     test_dense_F_Gamma_vs_oracle)
    fx["systems"][f"jensen_2lvl_ordered_{tag_b}"] = {
        "beta": float(beta_o), "E": [0.0, 1.0], "Delta": 1.0, "M2": 0.81, "m": 0.6}
```

Run `python3 verify_tensor_vertex.py`; fixture diff must be ADDITIVE only (existing rows bit-identical modulo provenance).

- [ ] **7A Step 3: MATLAB ordered V-row test.** In `test_invzt_vertex.m`: clone `test_dense_V_jensen_vs_oracle` (lines 77–102) → `test_dense_V_jensen_ordered_vs_oracle`, once per system (`_b14`, `_b3`): `S = fx.systems.jensen_2lvl_ordered_b14`, operator `X = centered([S.m sqrt(S.M2); sqrt(S.M2) -S.m], pcc)`, filter `startsWith(row.tags, 'V;cc;jensen2lvlord_b14')`, gate 1e-9. Add `verifyGreaterThan(tc, nrow, 0)` row-count protection to BOTH the new tests AND the existing PM `test_dense_V_jensen_vs_oracle` (round-2 P1-6: an empty filter must fail, not vacuously pass).

- [ ] **7A Step 4: the ordered Jensen DIAGNOSTIC (replaces the round-1 "exact bridge" — round-2 P0-2).** The moment-form Σ (eqs 37–38 / J 2.26–2.29) is derived with `g(ωm±ωn) → g(ωm)` (exact at n=0 ONLY) and a resummed-static elastic sector; the ordered bare propagator is `G0 = -M²g - m²h`, `h = β(1−n01²)δ_{n,0}` (J 2.7/2.8). Measured residuals of the naive `V == (−M²g)·Σ` bridge: 3.13e-5 (n=0), 4.89e-2 (n=1), 8.43e-2 (n=2) at βΔ=14. So: **gate n=0 only, with the elastic term included and a tolerance derived from the omitted-term scale; REPORT n>0 and all βΔ=3 rows.** Append to `test_invzt_vertex.m`:

```matlab
function test_jensen_ordered_static_diagnostic(testCase)
% Ordered moment-form comparison (J 2.26-2.29) -- STATIC GATE + FINITE-FREQUENCY
% REPORT (round-2 P0-2: the moment form is n=0-exact/finite-n-approximate BY
% DERIVATION -- g(wm+-wn) -> g(wm) -- so equality at n>0 must NOT be gated).
% This is an approximation diagnostic for a3d, NOT a correctness gate on the
% dense vertex (that gate is the oracle-row test above).
Delta = 1.0;  M2 = 0.81;  m = 0.6;  W = 1.3;
for ib = 1:2
    betas = [14, 3];  beta = betas(ib);
    n01 = tanh(beta*Delta/2);
    tl = struct('Delta', Delta, 'M2', M2, 'm', m, 'n01', n01, 'g0', 2*n01/Delta);
    Lmax = 400; kk = (0:Lmax).'; wts = [1; 2*ones(Lmax, 1)];
    wln = 2*pi*kk/beta;  Kint = 0.37 ./ (1 + (wln/W).^2);
    gint = invz_g(tl, 1i*wln);
    lam = invz_lambdas(Kint, gint, wts, beta, [1 2 3]);
    next = [0; 1; 2];  wext = 2*pi*next/beta;
    Kext = 0.37 ./ (1 + (wext/W).^2);
    gext = invz_g(tl, 1i*wext);
    sig = invz_sigma_ordered(tl, lam, Kext, gext, beta);
    h0  = beta*(1 - n01^2);                      % elastic term (J 2.8), n = 0 slot only
    G0  = -M2*gext;  G0(1) = G0(1) - m^2*h0;     % ordered bare propagator (J 2.7)
    Vjen = G0 .* sig.Sigma;
    % independent: exact four-point contraction, m != 0 (centered operator)
    Ecc = [0; Delta];  pcc = boltz(Ecc, beta);
    X   = centered([m sqrt(M2); sqrt(M2) -m], pcc);
    es  = struct('E', Ecc, 'p', pcc);  ops = struct('c', X);
    Kf  = @(ri, si, l) kjen(ri, si, l, beta, W);
    opts = struct('stage', 'V', 'Lmax', 400);
    opts.comps = {'c'};  opts.ext = {{'c', 'c'}};
    out = invzt_vertex4(es, ops, Kf, next.', beta, opts);
    Vme = out.val(1, :).';
    res = abs(Vme - Vjen) ./ max(abs(Vjen), 1e-300);
    if ib == 1
        % n = 0 STATIC GATE, low-elastic-weight regime. Tolerance derived from
        % the omitted elastic/static-resummation scale (30x margin over the
        % omitted-term ratio; floor 1e-4). Measured naive residual: 3.13e-5.
        tol0 = max(1e-4, 30 * m^2*h0 / (M2*abs(gext(1))));
        verifyLessThan(testCase, res(1), tol0);
        % derivation-shape check: the n=0 identity is MUCH better than the
        % n>0 static-value approximation (their measured ratio is ~1e3).
        verifyLessThan(testCase, res(1), 0.1*res(2));
        fprintf('[ordered static gate bD=14] res(n=0)=%.3e (tol %.3e); REPORT n=1: %.3e, n=2: %.3e\n', ...
            res(1), tol0, res(2), res(3));
    else
        fprintf(['[ordered diagnostic bD=3, REPORT ONLY] res = [%.3e %.3e %.3e] -- ' ...
            'elastic-sector resummation ambiguity (constraint 7), a3d input\n'], res);
    end
end
end
```

If the bD=14 n=0 gate fails at O(1): first suspect the centering/moment convention (`tl.m` = UNcentered ground-state diagonal vs the thermally-centered vertex operator) — resolve against J 2.7/2.26–2.29 in `jensen_1z_framework.html`, record the resolution in the test header. A residual persisting after convention resolution = STOP + record.

- [ ] **7A Step 5: run + commit**

Run: `$MB "cd('invz_tensor/tests'); r = runtests({'test_invzt_vertex','test_invzt_eps_el'}); disp(table(r)); assert(~any([r.Failed]))"`

```bash
git add invz_tensor/invzt_solve_point.m invz_tensor/tests/test_invzt_eps_el.m verify_tensor_vertex.py invz_tensor/tests/fixtures/vertex_oracle.json invz_tensor/tests/test_invzt_vertex.m
git commit -m "feat(invzt): a3d prerequisites -- eps_el elastic-only c_d, ordered oracle systems, static-gated Jensen diagnostic"
```

---

#### 7B — Fixed dominant vertex basis: constructor + ordered basis diagnostics

**EXECUTION AMENDMENT (2026-07-20, user-approved — supersedes the e2xI8 content definition below):** the fixed ZERO-FIELD `e2xI8` basis is **experimentally falsified**: measured `chi_share` = 0.00000–0.00013 across 3.0–4.65 T (overlap deficit 1.56–1.87) — the zero-field Ising doublet has ⟨+|Jz|−⟩ = 0 and the true low-energy eigenstates at 3–5 T transverse field have rotated far from the zero-field content, so the fixed subspace sees ~none of the cc response; the field-adapted lowest-16 captures **97.7%**. The vertex basis is therefore **FIXED-COUNT (16), FIELD-ADAPTED**: the 16 lowest eigenstates of the ordered mean-field Hamiltonian at each point — a direct truncation of `si_full` (no projector machinery). Constant dimension resolves the round-2 P0-3 count-drift objection; content evolution across the sweep is gated by continuity diagnostics. Revised contract: `vb = invzt_ordered_vertex_basis(ion, T, si_full, opts)` (T added — populations need beta) with `vb.E` [16×1] (= `si_full.E(1:16)`, already ground-shifted), `vb.Mz` [16×16] (= `si_full.Mz(1:16,1:16)`), `vb.p` [16×1] (Boltzmann at T, renormalized over the 16), `vb.V16` (= `si_full.V(:,1:16)`, for cross-field subspace diagnostics), `vb.n_full`, `vb.n_vertex = 16`, `vb.chi_share` (static cc ratio, elastic ON both legs), `vb.gap_16_17 = si_full.E(17) - si_full.E(16)` (level-crossing monitor). Test gates: at 3.0 T — `n_vertex == 16`, `chi_share > 0.9`, `gap_16_17` finite positive (REPORT its value); continuity over B = [3.0 3.5 4.0 4.4 4.65] — `chi_share` jumps < 0.05, adjacent-field subspace overlap `min(svd(V16(Bk)'*V16(Bk+1)))` > 0.9 (REPORT the values), E finite. The falsification measurement is recorded in the 7B report + ODD-LOG (Task 8).

**REFINEMENT (2026-07-20, Codex T7B feedback, verified + user-directed):** the object is the **fixed-RANK, field-adapted spectral subspace** (rank invariant under rotations within the 16; only a closing 16/17 gap threatens the truncation). Contract refinements: (i) populations from the converged state — `vb.p_mass = sum(si_full.P(1:16))`, `vb.p = si_full.P(1:16)/vb.p_mass` (NOT recomputed Boltzmann; `p_mass` is an approximation diagnostic, measured 1.000000, soft-gated > 0.99); (ii) carry **`vb.Mx`/`vb.My`** too — the a3d hybrid needs the full 3×3 truncated response tensor `chi_dom`; (iii) isolation diagnostics recorded, not the vacuous `gap > 0`: `vb.cluster_width = E(16) − E(1)`, `vb.gap_ratio = gap_16_17/max(cluster_width, eps)` (measured 2.71–3.66), `vb.gap_kBT = gap_16_17/(kB·T)`; (iv) **`vb.var_share`** = truncated/full equal-time connected Jz variance (measured 0.6650 at 3.0 T, 0.8381 at 4.65 T — the static `chi_share` ≈ 0.98 does NOT prove vertex convergence: high-energy states are denominator-suppressed statically but enter the four-point vertex as virtual intermediates); (v) `vb.min_subspace_overlap`/`vb.projector_distance = sqrt(1 − overlap²)` when a previous `vb` is supplied for continuity; (vi) the vertex operator centering (in 7C/7D consumption) uses the NORMALIZED dominant populations: `Mz_centered = vb.Mz − real(sum(vb.p .* diag(vb.Mz)))·eye(16)`; (vii) optional Procrustes alignment of `V16` to a previous subspace is gauge only — vertex results must not depend on it (YAGNI until needed). CONSEQUENCE for 7C/7E: a **direct `Vcc` basis-convergence study is REQUIRED before production** — see the 7C benchmark step and the 7E gate table.

**Files:**
- Modify: `invz_tensor/invzt_rung_basis.m` (document + test `'e2xI8'`; the label regexp already admits it — verify the constructor actually builds the 2-dim electronic ground doublet ⊗ I8 = 16 correctly)
- Create: `invz_tensor/invzt_ordered_vertex_basis.m`
- Test: `invz_tensor/tests/test_invzt_ordered_vertex_basis.m`

**Why fixed (round-2 P0-3, measured):** the instantaneous-energy cut `E < 0.4653` selects 13/13/11/10/9/8 states at B = 3.0/3.5/4.0/4.4/4.65/4.8 T — a moving basis that introduces cutoff jumps in a field sweep. The dominant basis is defined by CONTENT: the zero-field electronic ground doublet tensored with the complete nuclear space, dimension 16, field-independent.

**Interfaces:**
- Produces: `vb = invzt_ordered_vertex_basis(ion, si_full, opts)` with fields:
  `vb.P` [136×16] projector (from `invzt_rung_basis(ion,'e2xI8')`), `vb.E` [16×1] and `vb.Mz` [16×16] — the ordered mean-field Hamiltonian and Jz REBUILT IN the fixed basis at the CONVERGED (hx, hz) of `si_full` (one mean field, two representations: NO separate MF iteration in the reduced space), `vb.p` [16×1] normalized populations, `vb.n_full = 136`, `vb.n_vertex = 16`, `vb.overlap_deficit` = `norm((eye(136) - P*P')*V16, 'fro')` with `V16` the 16 lowest `si_full` eigenvectors, `vb.chi_share` = `chi0_cc^vertex(0) / chi0_cc^full(0)` (static susceptibility captured by the basis).

- [ ] **7B Step 1: failing tests** — `test_invzt_ordered_vertex_basis.m`: (i) `invzt_rung_basis(ion,'e2xI8')` returns `dim_actual == 16` with an orthonormal projector; (ii) at B = 3.0 T ordered `si_full`, `vb.n_vertex == 16`, `vb.overlap_deficit < 0.05`, `vb.chi_share > 0.9` (REPORT the numbers; the two thresholds are direction gates — a fixed basis that captures < 90% of the static cc response cannot honestly drive a3d, STOP and reconsider); (iii) **field-continuity**: over B = [3.0 3.5 4.0 4.4 4.65], `vb.n_vertex` constant (16 by construction), `vb.E` finite, and `max |chi_share(B_{k+1}) - chi_share(B_k)| < 0.05` (no cutoff-induced jumps — the failure mode of the `Esplit` cut).
- [ ] **7B Step 2: implement** — read `invzt_rung_basis.m` and `invzt_solve_point.m`'s `local_rung_si` (lines 451–522, the existing reduced-basis Hamiltonian pattern) first; `invzt_ordered_vertex_basis` mirrors `local_rung_si`'s projected-operator construction but takes (hx, hz) FIXED from `si_full` instead of iterating.
- [ ] **7B Step 3: run, pass, commit** — `feat(invzt): fixed e2xI8 ordered vertex basis with continuity diagnostics (a3d 7B)`.

---

#### 7C — Compact `cc;cc` vertex path: memory fix, tiling, budget guards, benchmark

**Files:**
- Modify: `invz_tensor/invzt_sigma_tensor.m` (dominant path, lines 122–145 + contraction)
- Test: extend `invz_tensor/tests/test_invzt_a3_threestate.m` (compact-vs-old equality) + new budget-guard tests

**Measured problem (round-2 P1-1/P1-2):** line 141 zero-pads the compact `G4dom` into `complex(zeros(npair,nc,nc,nwn,nl))` ≈ **1.42 GB** at 0.1 K (nwn=740, nl=1479) though only the `cc;cc` slot is used; the dense walk is ~`6·N⁴·nwn·nl` ≈ 1.9e11 kernel evaluations at N=13. GPU is NOT an escape on this machine (`gpuDeviceCount('all') = 0` on Apple silicon — CPU path is the only verified path; optional CUDA backend deferred to a machine that has one).

- [ ] **7C Step 1: compact storage.** In dominant mode: keep `G4cc = G4dom` as `[nwn, nl]`, never allocate the 5-D array, and contract directly against `Kcc(l)` (read `contract_vertex`, lines 222–253, first; add a compact branch `Vcc(n) = (1/2beta) * sum_l wl(l) * Kcc(l) * G4cc(n,l)` with the same l-weighting the general branch uses). Regression test BEFORE removing the old path: on the T = 2 K three-state fixture, compact vs old dominant contraction equal to 1e-12 (both `Vmat(3,3,:)` and the downstream `chi_til`).
- [ ] **7C Step 2: tiling + budget guards.** Evaluate the `(n,l)` grid in `l`-tiles of `opts.tile_nl` (default 128) columns so peak temporary memory is bounded; EVERY tiling parameter enters the `gamma4_cached` cache key. Add guards, all failing with `invzt:orderedA3Budget` BEFORE any allocation: `max_vertex_states` (default 16), `max_vertex_work` (default 5e11 kernel evals, computed as `6*N^4*nwn*nl`), `max_vertex_bytes` (default 4e9, computed from the compact shapes + tile temporaries). Tests: a synthetic over-budget call (e.g. `max_vertex_work = 1e3`) errors with the right identifier and allocates nothing.
- [ ] **7C Step 3: benchmark gate.** Script step (recorded in ODD-LOG, not CORE): run N = 8, 10, 13, 16 on a REDUCED grid (T = 1 K or Ecut = 10) and fit the measured `N⁴·nwn·nl` scaling; extrapolate to the 0.1 K production grid. **Hard gate: extrapolated wall-clock ≤ 12 h** (the repo's existing budget anchor). If it fails: STOP — the follow-on options (on-the-fly contraction, Matsubara compression/tail treatment, transition/time-domain factorization) become the next plan; do NOT promise the production run.
- [ ] **7C Step 3b: VERTEX-BASIS CONVERGENCE study (required before production — 2026-07-20 refinement).** The rank-16 basis captures ≈ 98% of the static cc response but only 0.665–0.838 of the equal-time Jz variance (measured) — high-energy states enter the vertex as virtual intermediates, so static capture is NOT vertex convergence. On a reduced frequency grid at ≥ 1 ordered anchor (e.g. 0.1 K / 3 T with reduced Ecut, and the affordable 1 K anchor), compute the compact `Vcc` at `Nv = 16` vs `Nv = 24` (and `Nv = 32` if it fits the budget guard) and REPORT the static and frequency-dependent relative differences. Acceptance judgment (7E row): `Vcc(16)` vs `Vcc(24)` agreeing at the few-percent scale supports production at rank 16; a large spread is a STOP — the sanctioned fallbacks are (1) enlarging the field-adapted cluster through the next robust energy gap, or (2) a response-augmented/block-Krylov basis seeded by the low-energy cluster and the omitted `Q·Jz·P` directions. NEVER fall back to the zero-field projector, a moving `Esplit` cutoff, or an untracked ranking of individual eigenstates.
- [ ] **7C Step 4: commit** — `perf(invzt): compact cc-only dominant vertex + tiling + invzt:orderedA3Budget guards (a3d 7C)`.

---

#### 7D — `invzt_solve_point_ordered(..., mode='a3d')`

**EXECUTION AMENDMENT (2026-07-20, user-required — supersedes the step-2/3 wiring below):** the a3d fixed point iterates V on the **COMPLETE HYBRID MAP** — the defining invariant is that the returned `Vcc` is generated by the returned `Kmat`:

```text
chi_dom_til(V) = dress[chi_16, V]                       (dominant Dyson dressing, cc channel)
chi_H(V)      = chi_full + (chi_dom_til(V) - chi_16)
K_H(V)        = EMT[chi_H(V)]                            (invzt_emt_matrix on lat_eff)
V_new(n)      = (1/2beta) sum_l K_H,cc(l) * Gamma16_cc;cc(n,l)
```

Rationale (user, 2026-07-20): a fixed point closed on the isolated 16-state system (`V = Γ₁₆K₁₆[χ̃₁₆(V)]`) with the hybrid assembled afterward returns a `Vcc` NOT generated by the returned `Kmat` — violating the advertised full-response hybrid. The omitted 119-state response modifies K (including transverse/ODD channels) at the SAME vertex order — not automatically O(1/z²) — and `var_share` = 0.665 at 3 T says spectator feedback is not negligible. The expensive Γ¹⁶ stays rank-16/compact and K-independent (one-time build); only the cheap dress→EMT→contract cycle iterates.

Implementation: extend `invzt_sigma_tensor` with `opts.chi_base` (fixed `[3,3,nwn]`, default zeros = bit-identical legacy) added to its internally-dressed response BEFORE the EMT step; a3d passes arg-1 `si = si_vb` (16-state) + `chi_base = chi_full − chi_16`. Returned `chi_til` = the hybrid `chi_H`; `pt.Kmat` = the hybrid-map medium; crit from `chi_H(:,:,1)`.

HARD GATES (all in CORE at the affordable anchor): (1) **self-generation** — `Vcc == contract_cc(Gamma16, pt.Kmat)` within the fixed-point tolerance; (2) **complete-map reeval** reproduces `Vcc`/`Kmat`/`chi_til`; (3) **reduction identity** — with `chi_base = 0` (equivalently `chi_full == chi_16`) the hybrid solver reduces EXACTLY to the truncated inner solver; (4) REPORT the as-built-truncated vs hybrid difference at the affordable anchor(s). The as-built truncated solve is retained ONLY as a warm-start/diagnostic comparator — it is NOT a3d (it would have to be renamed "self-consistent 16-state vertex model embedded post hoc" to ship as a theory, which is weaker than the stated full-response hybrid).

**Files:**
- Modify: `invz_tensor/invzt_solve_point_ordered.m` (mode gate + a3d branch)
- Modify: `invz_tensor/invzt_sigma_tensor.m` (accept an EXPLICIT dominant basis: `opts.dom_basis = struct('E','p','Mz')` from 7B's `vb` — when present, the vertex uses it and NEVER infers membership from `Esplit`; fixes round-2 P1-4 for the ordered route)
- Test: extend `invz_tensor/tests/test_invzt_solve_point_ordered.m`

**The a3d approximation, stated once (round-2 P1-3 — ONE consistent response everywhere):**

```text
chi_til = chi_full + (chi_dom_til - chi_dom)      % "dominant dressed, rest bare"
```

`chi_full` = bare ordered local response (all 136 states); `chi_dom` = bare response of the 16-state vertex basis; `chi_dom_til` = the vertex-dressed dominant response from the 7C compact solve. THIS `chi_til` — not two different maps — feeds (i) the EMT medium loop, (ii) the returned `pt.chi_til`, and (iii) `invzt_crit_static(herm_real(chi_til(:,:,1)), lat.JtGamma)`. `Sigma_cc_equiv`, if reported, is defined against the MATCHING dominant propagator (`V_cc/G0_dom_cc`) and labeled DIAGNOSTIC.

**Output surface (round-2 P0-1, honest):** `Vcc [nwn,1]`, `Kmat`, `chi_til [3,3,nwn]`, `Sigma_cc_equiv` (diagnostic), `vb` basis diagnostics (n_full/n_vertex/overlap_deficit/chi_share), budget/convergence data, `eps_el`/`c_d` (7A formula on the ordered si), `crit` from the hybrid `chi_til`. **`alpha`, `alpha_m`, `lambda` = NaN — they are Jensen moment-form objects the dense vertex does not produce; fabricating them is forbidden.** `pt.mode = 'a3d'`; the mode gate becomes `ismember(char(mode), {'a1','a3d'})`.

- [ ] **7D Step 1: failing tests** — append to `test_invzt_solve_point_ordered.m` (small/affordable grid for CORE — e.g. T = 1 K, `Ecut` reduced, gridN 8, all knobs printed; the production 0.1 K run is `INVZ_SLOW` only):

```matlab
function test_ordered_a3d_point(tc)
% a3d Matsubara solve on an AFFORDABLE anchor: converged physical objects,
% honest surface (no fabricated Jensen fields), consistent criticality.
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point_ordered(ion, 1.0, [2.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
if ~pt.is_ordered, fprintf('anchor relaxed PM -- pick a lower B\n'); end
verifyTrue(tc, pt.is_ordered && pt.converged);
verifyTrue(tc, all(isfinite(pt.Vcc)) && all(isfinite(pt.chi_til(:))));
verifyTrue(tc, isnan(pt.alpha) && isnan(pt.alpha_m));        % NOT fabricated
verifyEqual(tc, pt.vb.n_vertex, 16);
verifyTrue(tc, isfinite(pt.eps_el) && isfinite(pt.c_d));
% criticality consistency (7E gate): crit is recomputed from the RETURNED chi_til
cr = invzt_crit_static(real((pt.chi_til(:,:,1) + pt.chi_til(:,:,1)')/2), lat8.JtGamma);
verifyEqual(tc, pt.crit, cr, 'AbsTol', 1e-10);
% real-axis stays rejected
verifyError(tc, @() invzt_chi_realaxis(ion, 1.0, [2.0 0 0], pt, ...
    linspace(0, 0.6, 11).', struct()), 'invzt:realaxisMode');
fprintf('a3d anchor: m0=%.4f crit=%.5f chi_share=%.4f eps_el=%.3g\n', ...
    pt.m0, pt.crit, pt.vb.chi_share, pt.eps_el);
end

function test_a3d_vs_a1_approximation_control(tc)
% 7E approximation-control gate: a3d vs ordered a1 at the same affordable
% anchor -- REPORT the spread (the beyond-Jensen content + basis truncation),
% and require same-sign, same-order crit (a wildly different crit means a bug,
% not physics).
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
p1 = invzt_solve_point_ordered(ion, 1.0, [2.0 0 0], lat8, struct('Ecut', 10));
p3 = invzt_solve_point_ordered(ion, 1.0, [2.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
verifyTrue(tc, p1.converged && p3.converged);
verifyEqual(tc, sign(p3.crit), sign(p1.crit));
verifyLessThan(tc, abs(p3.crit - p1.crit), max(0.5*abs(p1.crit), 0.05));
fprintf('a3d vs a1: crit %.5f / %.5f, dSigma_cc(0) = %.4g\n', ...
    p3.crit, p1.crit, real(p3.Sigma_cc_equiv(1)) - p1.Sigma0);
end

function test_a3d_fixed_point_consistency(tc)
% 7E gate: one re-evaluation of the a3d map at the returned state reproduces
% Vcc/Kmat/chi_til within tolerance (the a3d analog of the a1 check-before-mix
% test). Implement via a solver diagnostic hook: opts.reeval = pt returns the
% one-pass images; spec fixed here, wiring lives in the a3d branch.
ion = tc.TestData.ion;
g8  = invzt_qgrid(8, 'halfopen');
lat8 = invzt_jq_tensor(ion, g8, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point_ordered(ion, 1.0, [2.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10));
re = invzt_solve_point_ordered(ion, 1.0, [2.0 0 0], lat8, ...
    struct('mode', 'a3d', 'Ecut', 10, 'reeval', pt));
verifyLessThan(tc, max(abs(re.Vcc - pt.Vcc)), 1e-6);
verifyLessThan(tc, max(abs(re.chi_til(:) - pt.chi_til(:))), 1e-6);
end
```

- [ ] **7D Step 2: implement the a3d branch** — mode gate extended to `{'a1','a3d'}`; the a1 front half (guards, ordered MF `si_full`, tl, splits, odd mask, C2 assert) reused verbatim; then: `vb = invzt_ordered_vertex_basis(ion, T, si, struct())` (fixed-count field-adapted, 7B amendment) → compact 7C vertex on `vb` (passed as `opts.dom_basis` — never `Esplit`) → hybrid `chi_til` → EMT loop (`invzt_emt_matrix` on `chi_til`) iterated with check-before-mix → crit from the returned `chi_til(:,:,1)` via `invzt_crit_static` → assemble the honest surface above (`alpha/alpha_m/lambda = NaN`). `opts.reeval` hook: skip the loop, evaluate the map ONCE at `reeval`'s state, return the images.
- [ ] **7D Step 3: production-grid run under `INVZ_SLOW` only** — a single 0.1 K/3 T a3d point behind the `INVZ_SLOW` env gate (goes on the CORE Incomplete allowlist), preceded by the 7C budget guard check; REPORT its numbers to ODD-LOG.
- [ ] **7D Step 4: run, pass, commit** — `feat(invzt): a3d ordered solve -- full-response dominant-vertex hybrid, one consistent chi_til (7D)`.

---

#### 7E — Acceptance gates (all must hold before Task 8 documents a3d)

| Gate | Requirement |
|---|---|
| Vertex correctness | compact `cc;cc` equals the general/old dominant contraction on small fixtures (1e-12) |
| Ordered correctness | centered m≠0 vertex matches the independent oracle rows (1e-9) |
| Jensen comparison | STATIC n=0 gate with derived tolerance only; n>0 and βΔ=3 residuals REPORTED |
| Basis stability | fixed 16-state vertex basis; no cutoff jumps over the 3.0–4.65 T scan (7B continuity test) |
| Fixed-point consistency | `reeval` reproduces `Vcc/Kmat/chi_til` within 1e-6 |
| Criticality consistency | returned `crit` recomputed from returned `chi_til(:,:,1)` (1e-10) |
| Approximation control | a3d vs ordered a1 spread reported at the affordable anchor; same-sign crit |
| Performance | 7C reduced-grid extrapolation ≤ 12 h before any production 0.1 K promise |
| Vertex-basis convergence | `Vcc(Nv=16)` vs `Vcc(Nv=24)` (32 if affordable) compared on a reduced grid at ordered anchor(s) BEFORE production; static + frequency-dependent diffs reported; divergence → sanctioned fallbacks only (next-gap enlargement / response-augmented basis) |
| Suite placement | structural tests in CORE (affordable anchors); production run `INVZ_SLOW` only |

---

### Task 8: Documentation + logs + final acceptance

**Files:**
- Modify: `invz_tensor/README.html`
- Modify: `docs/ODD-LOG.md`
- Create: `invz_tensor/SESSION-2026-07-19-invzt-ordered-side.md`

- [ ] **Step 1: README** — mode table: `a1` gains "ordered: solver + real-axis + stability-based dispatcher (2026-07-19)"; new row `a3d` — "**full-response, fixed-rank field-adapted dominant-vertex ordered a3d** (NEVER call it 'full tensor ordered a3' — round-2 minor): Matsubara-only, rank-16 field-adapted vertex basis (zero-field e2xI8 falsified, chi_share ~ 0; Esplit cut drifts; scorecard 48%/0%/97.7%), compact cc-only storage, `invzt:orderedA3Budget` guards; real-axis OPEN"; `a3` row notes full-dress ordered permanently budget-refused. `invzt_run_spectra` scope text: spans Bc, stability-based phase rule, transverse-only, knob rationale with the three measured failure numbers. Correct the `invzt_chi0_split` header claim that `E < Esplit` selects exactly ndom = 16 with hyp — measured FALSE for ordered 0.1 K states (13→8 across the sweep); the doc must scope that claim to the PM states it was written for.

- [ ] **Step 2: ODD-LOG §A-ordered** — implemented items, suite status, headline numbers (m0/Sigma0/crit at 3 T; interop dSigma0/dalpha_m; the 4.8 T stability-selection lock; boundary bracket; ordered static-gate n=0 residual + the n>0/βΔ=3 REPORT rows; **eps_el correction note: prior logged eps_el values were upper bounds** (full variance, P1-4); a3d affordable-anchor numbers (crit vs a1, chi_share, overlap_deficit) + the 7C benchmark scaling fit and extrapolated production wall-clock; driver phase summary + masked-band width), grid conv + dpRng provenance.

- [ ] **Step 3: Session doc** — the driver-knob diagnosis numbers; the P0-1 m0-persistence finding and the stability-based rule; the P0-2 infeasibility numbers (740 wn, ~1.4 GB Gamma4, budget-refused 136-basis) and the follow-on split; scope decisions with dates.

- [ ] **Step 4: Final acceptance — full suites with the Incomplete allowlist (review P1-5)**

Run: `$MB "cd('invz_tensor/tests'); r = runtests(pwd); inc = {r([r.Incomplete]).Name}; disp(inc); assert(~any([r.Failed]))"` — the printed Incomplete set must be exactly the known allowlist (baseline 1 `INVZ_SLOW` name, plus the 7D `INVZ_SLOW` production a3d point, plus `test_boundary_bracket` IF it was slow-gated in Task 5).
Run the projected suite the same way (allowlist: the 19 known names). Where budget permits, run the slow suite once: `INVZ_SLOW=1` prefixed to the same commands, expecting zero failures and an empty (or reduced) Incomplete set.

- [ ] **Step 5: Commit**

```bash
git add invz_tensor/README.html docs/ODD-LOG.md invz_tensor/SESSION-2026-07-19-invzt-ordered-side.md
git commit -m "docs(invzt): ordered-side a1 + a3 prerequisites -- README scope, ODD-LOG A-ordered entry, session log"
```

---

## Review Disposition (2026-07-19 Codex review — all findings verified against the codebase before adoption)

| Finding | Verdict after verification | Disposition |
|---|---|---|
| P0-1 ordered-first misassigns [Bc, ~5.0 T] | CONFIRMED (m0 probe reproduced bit-identically: 1.6553/1.5109/1.1717 at 4.65/4.70/4.80 T) | Task 5 redesigned: stability-based PM-first rule; hard 4.8 T gate; boundary-bracket scan; ordered stability criterion `crit > -1e-3` explicit |
| P0-2 ordered a3 infeasible at `std`/0.1 K | CONFIRMED (README ladder gate budget-refuses `e17xI8`; Gamma4 ≈ 1.4 GB at 740 wn) | Ordered-a3 SOLVE removed from this plan (mode `'a1'` only, `invzt:orderedMode`); Task 7 reduced to the two prerequisites; follow-on plan owns the basis-vs-optimization decision |
| P1-1 real-axis tests don't gate the formula | CONFIRMED (current code applies the PM formula to ordered pt silently) | Task 4: exact `Sigma_w` formula gate that MUST fail pre-change; `verifyTrue` anchors; force_sigma0 demoted to branch-independent regression |
| P1-2 nondeterministic rejection test | CONFIRMED | Task 5: `m_tol = Inf` deterministic variant; asserts `di.para.converged && di.para.crit < 0` |
| P1-3 oracle step not executable; wrong target | CONFIRMED (accessors are `TestData.fx`/`row.tags`/`value_re,value_im`; PM bridge is `test_jensen_bridge`) | Task 7 rebuilt on the real patterns; the load-bearing test is now the CONTRACTED ordered bridge `V == G0 .* Sigma_ordered` incl `alpha_m`, exact-gated in the dead-elastic regime, REPORT at moderate `h0` |
| P1-4 `eps_el` uses full variance | CONFIRMED (`invzt_solve_point.m:259`) | Task 7 Step 1: elastic-only `c_d` from chi0 elastic-on/off; `pt.c_d`; ODD-LOG upper-bound note |
| P1-5 suite commands fail at baseline | CONFIRMED (CORE 47/0/1, PROJECTED 143/0/19) | All gates now `assert(~any([r.Failed]))` + Incomplete-allowlist check at final acceptance; `INVZ_SLOW=1` slow run where budget permits |
| P1-6 `git add -A` captures user work | CONFIRMED (dirty worktree) | All commits path-scoped; `invz_projected/README.html` update added to Task 1 |
| P2-1 `di` loses returned-but-rejected outcomes | CONFIRMED | Task 5: structured per-leg `di` (attempted/accepted/converged/m0/crit/Sigma0/err); Task 6 logs it |
| P2-2 mix-after-evaluate inconsistency | CONFIRMED (my draft copied the projected ordered loop, which shares the flaw) | Task 3: check-before-mix (PM-solver ordering) + self-consistency test; parity-gate impact noted (≤ tol/mixo, far below 5e-3) |

## Review Disposition — round 2 (2026-07-19 Codex re-review; all findings verified before adoption)

| Finding | Verdict after verification | Disposition |
|---|---|---|
| P0-1 plan excluded the selected a3d deliverable | CONFIRMED (re-review header records the user's dominant-basis decision) | a3d in scope: Task 7 restaged 7A–7E; mode gate `{'a1','a3d'}` at 7D, full-dress `'a3'` permanently rejected; dispatcher stays a1-only (`invzt:autoMode`); Matsubara-only; honest surface (no fabricated `alpha/alpha_m/lambda`) |
| P0-2 ordered Jensen bridge not an exact identity | CONFIRMED (HTML: `G0 = -M²g - m²h`, J 2.7/2.8; "g(ωm±ωn) → static value — exact at n=0, approximate at finite frequency"; my βΔ=14 elastic weight was wrong: 3.3261e-6, not 6.6e-7) | 7A Step 4: static n=0 gate with the elastic term in `G0` and a derived tolerance (30× omitted-term scale, floor 1e-4) + shape check res(0) < 0.1·res(1); n>0 and βΔ=3 rows REPORT-only; never blocks a3d |
| P0-3 `E < Esplit` unstable ordered basis | CONFIRMED bit-identically (ndom 13/13/11/10/9/8 at 3.0–4.8 T) | 7B: fixed `e2xI8` 16-state content-defined basis; `si_full` + `vb` split; continuity test; `invzt_chi0_split` doc claim corrected (Task 8) |
| P1-1 dominant path allocates the 1.42 GB zero-padded G4 | CONFIRMED (`invzt_sigma_tensor.m:141`) | 7C Step 1: compact `G4cc [nwn,nl]` + direct Kcc contraction; compact-vs-old 1e-12 regression before removing the old path |
| P1-2 whole-grid work ~1.9e11 kernel evals; GPU unavailable | CONFIRMED (arithmetic checks; `gpuDeviceCount('all') = 0` on this machine) | 7C: l-tiling (params in cache key), `invzt:orderedA3Budget` guards (states/work/bytes), N = 8/10/13/16 reduced-grid benchmark with a hard ≤ 12 h extrapolation gate; production run `INVZ_SLOW` only; CPU is the verified path |
| P1-3 medium and criticality use different maps | CONFIRMED (dominant crit uses `cdom/(1+Σ)+crest` while the medium resums full `G0`) | 7D: one hybrid `chi_til = chi_full + (chi_dom_til - chi_dom)` for EMT + returned state + crit; crit-recompute consistency test (1e-10); `Sigma_cc_equiv` vs the MATCHING dominant propagator, labeled diagnostic |
| P1-4 `Esplit` not forwarded to the vertex | CONFIRMED (`st_opts` carries rank_tol/mix/tol/max/seed/dress only) | 7D: explicit `opts.dom_basis` (from 7B's `vb`) — the ordered vertex never infers membership from `Esplit`; single-sourced provenance |
| P1-5 dispatcher catch too broad; Bz guard bypassable | CONFIRMED (round-1 draft guarded Bz only in the ordered leg; blanket invz:*/invzt:* catch) | Task 5: entry validation (Bz + mode) before either leg; enumerated RECOVERABLE allowlist; three new tests (longitudinal-at-PM-valid-field, invalid mode, malformed lattice rethrown) |
| P1-6 multi-beta fixture incompatible; placeholder helper | CONFIRMED (`sysdata` expects scalar `S.beta`; `local_a3_anchor` did not exist) | 7A: two scalar-beta systems `_b14`/`_b3`; self-contained `test_invzt_eps_el.m` (no invented helpers); `verifyNotEmpty` row-count protection added to PM + ordered V-row tests |
| Minor: "all 7 pass" count; naming; field-surface distinction | CONFIRMED | Count fixed ("count the functions, don't trust prose"); docs mandate "full-response, dominant-vertex ordered a3d"; a3d output surface separates meaningful fields from NaN'd Jensen fields |

## Self-Review (revised plan)

1. **Spec coverage:** across-QPT spectra deliverable = Tasks 3–6; a3d deliverable = Task 7 stages 7A–7E with the 7E gate table; user scope decisions (transverse-only; a3d dominant-basis; Matsubara-only; naming) encoded at the solver (`invzt:orderedLongitudinal`, `invzt:orderedMode`, `invzt:autoMode`, `invzt:orderedA3Budget`) and in the docs task. Every round-1 and round-2 "Required revision/change" maps to a task change (two disposition tables above).
2. **Placeholder scan:** Task 7 retains explicitly-flagged "read the file first" steps for `invzt_gamma4.m`/`invzt_vertex4.m`/`invzt_rung_basis.m` internals and the `jensen_2lvl` emission block — contracts (shapes, identifiers, gates) are fixed here; internals are read at implementation time. The 7D `reeval` hook is specced by its contract (skip loop, one-pass images) with wiring left to the a3d branch implementation. No TBD/TODO otherwise.
3. **Type consistency:** `pt` surface consistent across Tasks 3–6; a3d surface separated (7D: `Vcc/Kmat/chi_til/Sigma_cc_equiv/vb` meaningful, `alpha/alpha_m/lambda` NaN); `invzt_crit_static` signature consistent (Tasks 2/3/7D); `di` schema identical in Tasks 5/6; `lam` 3×1 + `invz_sigma_ordered(tl, lam, K, g, beta)` matches the verified signature; frozen-Kw `gamma_w == gamma0` identity stated in Task 4's interface and test comment; `vb` field names consistent between 7B and 7D.

**Known risks:** (a) the 7A static gate may still expose an m≠0 centering-convention subtlety — STOP and resolve against J 2.7/2.26–2.29, record, never tune; (b) the near-Bc Option-A window (bare-MF boundary ≈ 5.0 T vs QPT 4.65–4.70 T) may leave a visible masked band in the Task 6 sweep — the measured boundary-scheme gap, REPORTED; (c) the 7C benchmark may fail the ≤ 12 h extrapolation gate — then the production 0.1 K a3d run is NOT promised and the compression/factorization options become the follow-on (the affordable-anchor CORE tests still land); (d) `e2xI8` may capture < 90% of the ordered static cc response at some field (7B chi_share gate) — STOP and reconsider the basis before 7C/7D spend; (e) `invz_twolevel_ordered` at very low transverse field relies on hz for the doublet splitting — fine at B ≥ 3 T tested here; a future B → 0 sweep should add a dedicated test.
