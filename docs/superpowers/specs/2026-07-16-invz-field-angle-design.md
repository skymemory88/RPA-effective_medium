# c-axis field-misalignment (tilt) knob for `invz_run_spectra`

Date: 2026-07-16. Branch: `invz-1z-lihof4`. Status: design approved by user;
**revised 2026-07-16** after external review, **second revision 2026-07-16**
after the follow-up review (`field-angle-plan-review_by_Codex.md`, both rounds;
numerical findings verified against the code — see "Review resolutions" at the
end).
Scope: **scalar stage, c-axis tilt only** (`phi_ab = 0`); azimuthal support and
tensor propagation are deferred follow-ups (§8).

## Problem

`invz_run_spectra.m` sweeps a scalar field magnitude and treats it as a purely
**transverse** field along the crystal **a**-axis: every solver in the chain
turns the scalar `Bx` into the literal vector `[Bx 0 0]` (ordering axis = **c** =
z). Real LiHoF4 experiments have a small **misalignment** of the nominally
transverse field toward **c**; the resulting longitudinal component `Bz` rounds
the sharp quantum phase transition into a crossover and is a known source of
experiment/theory discrepancy. The user wants a field-angle knob in the 1/z
spectra driver, following the convention of `LiReF4_MF_Yikai.m` (located at
`/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/
Matlab/Simulation/Mean Field/LiReF4/LiReF4_MF_Yikai.m`, lines 60-66 — the
**formula** (`Hx = B·cosθ·cosφ, Hz = B·sinθ`) is authoritative; that file's
line-60 comment "angle from c-axis" contradicts its own formula and should be
ignored. Companion files there are `MF_Chi.m` (capital C) and `RPA.m`; the
reviewer confirmed `RPA.m` inverts the full 3×3 response, with a
diagonal-projected interaction).

## Physics background — why longitudinal field is *not* trivial in 1/z

At the **MF** and **RPA** levels a longitudinal field is trivial: both are
**full-tensor** methods. MF diagonalizes `H(B)` for any `B = [Bx By Bz]` and
self-consistifies all of `⟨Jx⟩,⟨Jy⟩,⟨Jz⟩`; RPA builds the full 3×3 `χ0(ω)` from
all transitions and inverts `χ = [χ0⁻¹ − 𝒥(q)]⁻¹`. A tilt just changes the
numbers in the same matrices.

The **1/z self-energy is a different kind of object.** Jensen's `Σ(ω)` is a
*scalar* built from one transition — the field-split Ising doublet — and it
renormalizes **only the cc (c-axis) component**: the code computes
`χ̃0_cc = χ0_cc/(1+Σ)` and never forms `χ0_xx`, `χ0_xz`, or any off-diagonal.
Every ingredient (`α`, `γ`, the `λp` Matsubara sums, the ordered `α_m`, the
elastic sector) is a projection onto the doublet parameters
`(Δ, M²=|⟨0|Jz|1⟩|², n01, m=⟨0|Jz|0⟩)`.

Consequences for a tilted field:

- **The moment-carrying machinery is already the general case.** A longitudinal
  `Bz` adds `−gL·µB·Bz·Jz` to `H0`, inducing `⟨Jz⟩ ≠ 0` at every field —
  structurally identical to the ordered branch's longitudinal molecular field
  `hz`. `invz_sigma_ordered` (the m≠0 self-energy) reduces exactly to
  `invz_sigma` as `m→0` (`alpha_m → 0`); `invz_chi0z` already folds the induced
  moment into the elastic term via `Jexp`; the sum rule already uses the
  variance `JzJz_fluct = ⟨Jz²⟩−⟨Jz⟩²`. The induced-moment case flows through
  the existing ordered path with no new self-energy.

- **Accuracy statement (corrected per review).** `χ_cc` is even in the tilt
  angle, so the entire tilt effect starts at O(θ²) — and the scalar port keeps
  only *part* of the O(θ²) coefficient (the nonperturbative doublet
  re-splitting and induced-moment effects) while omitting another same-order
  part (the `χ_zx·χ_xz` cross-channel contribution that a full tensor inversion
  would mix in; `χ_xz` is odd in θ and enters `χ_cc` quadratically). The
  defensible claim is therefore: **exact at zero tilt; uncontrolled relative
  error in the tilt-induced change, beginning at O(θ²)**. Near the former
  critical point the response to the conjugate field is non-analytic, so Taylor
  counting does not apply there at all — the supported angle range must be
  established **numerically** (see the Σ=0 tensor-reference test, §7), not by
  power counting. Mitigating expectation (to be measured, not assumed): the
  omitted cross terms are Van Vleck-scale because `⟨0|Jx|1⟩ ≈ 0` within the
  Ising doublet.

## Design decisions (resolved with the user; revised after review)

1. **Scope this stage to the c-axis tilt** (`theta_c`), the physically requested
   experimental misalignment. `phi_ab` support is **descoped** (review finding
   1): with `B64s ≠ 0` an x-field already induces a small perpendicular
   `⟨Jy⟩ ≈ −0.069` at 4 T (verified numerically) with no feedback channel, so a
   "rigorous azimuth" would require a two-channel (hx, hy) transverse mean
   field that is exactly C4-consistent — and that is **incompatible with a
   bit-for-bit x-field baseline**. The transverse MF in this stage remains the
   **legacy x-only approximation**, now documented as such. Two-channel MF +
   azimuth is a deferred follow-up (§8).
2. **Longitudinal routing through the moment-form self-energy** with
   sign-aware branch selection (review finding 2, verified: the default
   `mz_seed = +1` converges to the metastable anti-aligned branch for
   `Bz < 0`).
3. **The ODD extension does not help here.** `odd_implementation_plan.html`
   Tiers 1-2 are deliberately scalar-cc (orthogonal to an external longitudinal
   field). Only its Appendix A route is relevant: **A0 is the rigorous
   tensor-RPA layer; A1 remains a projected, dominant-transition 1/z
   approximation** (wording per review).
4. **Staged delivery:** scalar c-tilt now; azimuth/two-channel MF and A0(+A1)
   tensor propagation later (§8).

## Solution overview

Thread a full field **vector** through the solve chain in place of the hardcoded
`[Bx 0 0]`, backward-compatibly: a helper maps a scalar to `[B 0 0]` and passes
a 3-vector through. The driver builds `Bvec = |B|·[cosθc 0 sinθc]`. When
`|Bz|` exceeds a routing tolerance the solve goes exclusively through the
moment-carrying path with the spontaneous-moment gate bypassed
(`forced_moment`); below the tolerance the longitudinal component is **zeroed**
and today's transverse logic runs verbatim. All changes are inert at
`theta_c = 0`: the existing suite and published benchmarks must reproduce
bit-for-bit.

## Components

### 1. New helper `invz/invz_field_vec.m`

```matlab
function B = invz_field_vec(B)
%INVZ_FIELD_VEC Normalize a field argument to a 1x3 row [Bx By Bz] in Tesla.
```

Contract (review finding 8): input must be numeric, real, finite, and either a
scalar (→ `[B 0 0]`, the historical convention) or exactly three elements (row
or column, normalized to a `1x3` row). Anything else — NaN/Inf, complex, empty,
wrong length — errors with stable identifier `invz:fieldVec`. Used at every
site that currently writes `[Bx 0 0]`.

### 2. `invz/invz_single_ion.m` — branch seeding + MF convergence reporting

**No `hy` channel** (descoped, decision 1). Changes:

- **Sign-aware seed:** in `order` mode, default
  `mz_seed = sign(B(3))` when `B(3) ≠ 0` (else `+1` as today), so an explicit
  longitudinal field selects the aligned branch. `opts.mz_seed` still
  overrides. Verified failure this fixes: `B = [2 0 −0.01]` T, seed `+1` →
  `⟨Jz⟩ = +4.815` (metastable); seed `−1` → `−4.86686`, the exact Z2 mirror of
  the `+Bz` result.
- **MF convergence surfaced** (review finding 6): return `si.mf_converged`
  (logical), `si.mf_iters`, `si.mf_residual` (final `dmf`); keep the existing
  warning. No behavior change for existing callers (additive fields).
- **Energy diagnostics** (second review finding 1): return `si.E0` (the
  unshifted ground-state energy `E(1)` — `si.E` stays shifted as today) and
  `si.F_mf`, the **variational** mean-field free energy
  `F_MF = −kT·ln Tr e^{−βH} + hx²/(2·Jxx0) + hz²/(2·J0z)` (the double-counting
  correction restores `−½J⟨J⟩²` per interaction channel; use the unshifted
  spectrum). The naive shifted-spectrum comparison is **wrong**: at
  `T = 0.31 K`, `B = [2 0 −0.01]` it mis-ranks the branches
  (−5.15e-8 vs −2.74e-8 meV), while `F_mf` correctly puts the aligned branch
  lower: −21.47664574 vs −21.46963936 meV (verified; matches review).

### 3. Field-vector plumbing (the 5 `[Bx 0 0]` leaf sites)

Each becomes `invz_field_vec(B)` and accepts scalar-or-vector `B`:

- `invz_twolevel.m:7` — keep the `m=0` gate (only reached when the longitudinal
  component has been zeroed by the routing dead band, so `m = 0` holds exactly).
- `invz_twolevel_ordered.m:13` — rebuild the doublet with the full vector + the
  fixed MF `hz`. No double-counting (external `Bz` in `H0`; `hz` is the MF
  piece). The `hz` handed over comes from the sign-selected branch of §2, so
  branch consistency is automatic.
- `invz_solve_point.m:16,19`, `invz_solve_point_ordered.m:34,45`,
  `invz_chi_realaxis.m:37` — pass the vector through.

`invz_critical.m`, `invz_critical_T.m`, `invz_critical_T0field.m`,
`invz_run_phase_diagram.m` keep passing scalars → unchanged. The angle knob is
**not** added to the phase-diagram driver (no sharp boundary under `Bz`).

### 4. Longitudinal routing — `invz_solve_auto.m` + `invz_solve_point_ordered.m`

**Single routing rule** (review finding 8): named option `opts.bz_tol`
(default `1e-9` T). In `invz_solve_auto`, after `Bvec = invz_field_vec(B)`:

- `|Bvec(3)| <= bz_tol`: **zero the component** (`Bvec(3) = 0`) and run today's
  logic verbatim (ordered-first, paramagnet fallback). The dead band can never
  reach the strict-paramagnet `m` gate with a material moment because the
  component is exactly zero.
- `|Bvec(3)| > bz_tol`: route **exclusively** to `invz_solve_point_ordered`
  with `opts.forced_moment = true`. Never attempt `invz_solve_point` (its
  two-level gate rejects `m ≠ 0`).

`forced_moment` semantics in `invz_solve_point_ordered` (review finding 6 —
explicitly non-circular, in order):

1. Bypass the early `abs(m0) > mtol` paramagnetic return.
2. Require `si.mf_converged` (else return non-converged immediately).
3. Assert `sign(si.Jexp(3)) == sign(Bvec(3))`; on mismatch, re-solve once with
   the mirrored seed and keep the field-aligned solution; if still anti-aligned,
   return non-converged with a warning.
4. Run the outer EMT⇆Σ loop as today; final acceptance uses the same
   finite-`Sigma0` and medium-convergence checks as the existing route.
5. Only then set the moment-form flags: `pt.is_ordered = true` (documented
   strictly as "uses the moment-form self-energy") **plus** new machine-readable
   `pt.moment_branch = 'spontaneous' | 'field_induced'`.

**Return contracts** (second review finding 2 + refinement 4):

- **Every** return path of `invz_solve_point_ordered` — including the early
  paramagnetic return and the new step-2/3 failures — populates the same field
  set: `is_ordered`, `converged`, `Sigma0`, `crit`, `si`, `tl` (may be `[]`),
  `m0`, `moment_branch`. Callers never probe a missing struct member.
- On a longitudinal failure (`phase = 0`), `invz_solve_auto` returns the failed
  ordered-style `pto` (not `[]`) whenever it carries valid `si`/`tl`, so the
  map can implement its RPA-only fallback. `pt = []` only when no usable
  single-ion state exists.

### 5. `invz/invz_spectra_map.m` — direction API, failure contract, metadata

**API** (review finding 4 — `fields` already means a list of scalar sweep
values, so a 3-element field input would be ambiguous):

```matlab
opts.field_dir = [1 0 0];   % unit-normalized internally; nonzero finite real 3-vector
```

`fields` stays a vector of **nonnegative magnitudes** (validated); internally an
`nB x 3` array `Bvec(k,:) = fields(k)*dhat` is formed once, before the parfor,
**with the dead-band normalization already applied** (`|Bvec(k,3)| <= bz_tol`
zeroed) so that `S.Bvec` is exactly the vectors the solves use — the same rule
`invz_solve_auto` applies, resolved from the same option (second review finding
2; no requested-vs-used ambiguity). Returned metadata: `S.fields` (magnitudes),
`S.field_dir` (normalized), `S.Bvec` (vectors used).

**Solver-option propagation** (second review finding 2 — today both spectra
functions build a fresh `sopts` with only `hyp/J0eff/Jxx0`, so the failure-
injection tests would be unimplementable):

- New `opts.bz_tol` on both spectra functions, resolved **once** and used for
  (a) the pre-parfor dead-band normalization above, (b) the `sopts` passed to
  `invz_solve_auto`, and (c) `one_field`'s longitudinal-failure branch. One
  value, three consumers, no duplicated defaults.
- New `opts.solve_opts` (struct, default empty) on both spectra functions:
  merged into `sopts` after the live couplings are set. The fields `J0eff`,
  `Jxx0`, `hyp` are **reserved** (error `invz:solveOpts` if present — the
  driver owns them); everything else (`max_outer`, `mf_maxit`, `mix_outer`,
  `Ecut`, `m_tol`, `mz_seed`, ...) passes through to the point solvers. This is
  what lets tests force `max_outer = 1` / `mf_maxit = 1` through the map.

**Longitudinal failure contract** (review finding 5 — verified crash path:
`one_field` line 113 calls `invz_twolevel` *outside* the try block, so a
`phase = 0` longitudinal point would raise `invz:orderedPhase` and abort the
parfor): when `|Bz| > bz_tol`, `one_field` never falls through to the
strict-paramagnet overlay. `phase == 1` → unchanged ordered handling. Otherwise:
if the failed moment-branch `pt` still carries valid `si`/`tl`, compute the
RPA-only overlay from the ordered-style `pt0` and leave the 1/z column masked;
else mask the whole column. `S.phase` stays 0 there.

**Labels/docs** (review finding 7): header and phase strings distinguish
"spontaneous FM" from "field-induced moment (crossover)"; no code path may
label a `Bz ≠ 0` point as FM/PM phases of a sharp transition.

### 6. `invz/invz_spectra_qpath.m` — vector field + safe formatting

`B` accepts scalar-or-3-vector via `invz_field_vec`; returns `S.Bvec` (the
vector **used**, i.e. after the same `bz_tol` dead-band normalization) and
`S.Bmag`. Gains the same `opts.bz_tol` / `opts.solve_opts` contract as the map
(§5). The `invz:noSolution` error at lines 67-69 currently formats `B` with a
single `%.3f` — MATLAB recycles the format over a 3-vector, producing a
malformed message; rewrite with `mat2str`-style formatting (review finding 7).

### 7. Driver `invz/invz_run_spectra.m` + plot labels

One new knob near the existing ones:

```matlab
theta_c = 0;   % deg: tilt of the field OUT of the transverse ab-plane toward c (Ising axis)
```

`dhat = [cosd(theta_c) 0 sind(theta_c)]`, fixed across the sweep; fed to the
field-sweep view via `opts.field_dir` and to the q-path view as `Bq*dhat`.
`theta_c = 0` reproduces today's `[|B| 0 0]` exactly. Header documents: the
convention (matches `LiReF4_MF_Yikai.m` theta at phi = 0); that the sweep
x-axis is the total magnitude (c-component `|B|·sinθc`); the legacy x-only
transverse MF; the corrected accuracy statement; and the deferred follow-ups.

Label sweep (review finding 7): `invz_plot_spectra_map.m:32` and
`invz_run_spectra.m:134` change `'B_x (T)'` → `'|B| (T)'` (with the direction in
the title when `theta_c ≠ 0`).

**Σ=0 tensor-reference validation** (review finding 3; metric per second
review finding 3): a test/driver utility builds the full 3×3 Cartesian RPA from
`invz_chi0z` (already available) with the same diagonal couplings
`diag(Jaa0, Jaa0, Jsel)` and compares its `χ''_cc` against the scalar-chain
`Σ=0` result over representative `(B, ω)` at `theta_c ∈ {0, 0.5, 1, 2, 5}` deg.
Conditions and metrics are fixed **before** the measurement:

- Reference conditions: `ion.demag = 0` (intrinsic response on both sides);
  the **full 3×3** inversion — all cross channels retained, including any
  `yz`/`xy` blocks allowed by `B64s`, not an xz-only sub-block.
- Spectral metric, per (angle, field):
  `eps_spec = ||χ''_sc − χ''_ten||₂ / max(||χ''_ten||₂, 1e-12·max|χ''_ten|·√nw)`
  (L2 over the ω grid — robust at spectral zeros, unlike a pointwise relative
  max).
- Peak metric: `dE_peak = |Epeak_sc − Epeak_ten|`, censored peaks compared as
  in `invz_peak_energy` (both-NaN passes, one-sided NaN fails).
- **Support criterion** (spec defaults, adjustable at review): an angle is
  supported when `eps_spec <= 5%` AND `dE_peak <= max(0.02·Epeak_ten, eta)` at
  every tested field. The resulting supported range is stated in the
  README/session log; the logged numbers also get a 1% reproducibility
  assertion (slow test).

This does not require the 12×12 A0 build.

## Backward compatibility

- Scalar field anywhere → `[B 0 0]` → identical Hamiltonian, identical path.
- `theta_c = 0` → routing dead band zeroes nothing (component already 0) →
  today's exact branch logic; sign-aware seed inert at `B(3) = 0`.
- Phase-diagram / critical drivers and all their tests untouched.

**Non-negotiable:** with `theta_c = 0` the full existing suite and every
published benchmark (Σc, Tc(0), Hc, soft-mode minimum) reproduce bit-for-bit
(no model change remains in this scope — the hy channel was descoped).

## Testing

New tests under `invz/tests/` (fast unless noted). Explicit grids/tolerances
(review test-plan items 2-8); values below are the spec's commitments, refined
only with recorded justification:

1. **Regression:** full suite green at defaults; `invz_single_ion` output at
   `(T, [Bx 0 0])` bitwise-unchanged (fields added, values identical).
2. **Field-vector contract:** scalar `B` vs `[B 0 0]` identical at every
   scalar-accepting boundary (`invz_twolevel`, `invz_solve_point`,
   `invz_solve_point_ordered`, `invz_solve_auto`, `invz_chi_realaxis`,
   `invz_spectra_qpath`); for `invz_spectra_map` (whose third argument stays a
   magnitude list) the equivalent check is default `field_dir` vs explicit
   `[1 0 0]` (second-review refinement 3); row vs column 3-vectors identical;
   NaN/Inf/complex/empty/wrong-length inputs error with `invz:fieldVec`.
3. **Branch selection (finding 2, corrected per second review finding 1):** at
   `T = 0.31 K`, `B = [2 0 ∓0.01]`: `sign(⟨Jz⟩) == sign(Bz)` on both sides;
   `±Bz` single-ion states exact Z2 mirrors (`|⟨Jz⟩(+) + ⟨Jz⟩(−)| < 1e-10`);
   **variational** `si.F_mf` of the aligned branch strictly lower (verified
   margin at this point: 7.0e-3 meV — the naive shifted-spectrum comparison
   mis-ranks and must not be used).
4. **±Bz mirror of spectra:** `χ''_cc` at `θc = ±1°`, `T = 0.31 K`,
   `B = 5.0 T`, `w = (0:0.02:0.5)'` meV; explicit metric (second-review
   refinement 2): `max|a−b| / max(max|a|, 1e-12) < 1e-8`.
5. **Crossover continuity:** same `(T, w)`, `θc = 0.5°`, fields
   `4.6:0.05:5.3` T: no NaN `Epeak` in the crossover window, and
   `pt.sumrule_rel < 5e-2` at every field (same order as the existing
   ordered-phase tolerance).
6. **θc → 0 continuity:** at `T = 0.31 K`,
   `max_w |χ''(θc=10⁻³ deg) − χ''(0)| / max_w χ''(0) < 1e-6` at **two** fields:
   `B = 2 T` (spontaneous-ordered at zero tilt: moment-form → moment-form) and
   `B = 6 T` (paramagnetic at zero tilt: exercises the forced moment-form →
   strict-PM formula reduction; second-review refinement 1).
7. **Routing tolerance:** `|Bz|` just below `bz_tol` → transverse path
   (`pt.moment_branch = 'spontaneous'` or PM); just above → moment path
   (`'field_induced'`).
8. **Longitudinal failure masking (finding 5):** a deliberately non-converged
   longitudinal map point (e.g. `max_outer = 1`) yields a masked/RPA-only
   column, `S.phase = 0`, and the parfor completes.
9. **MF convergence gate:** a forced `mf_maxit = 1` longitudinal point returns
   non-converged (`si.mf_converged` respected in acceptance).
10. **Metadata/labels:** `S.field_dir`, `S.Bvec`, `S.Bmag` present and
    consistent; map axis label `|B| (T)`; qpath error message well-formed for
    vector `B`.
11. **Σ=0 tensor reference (finding 3, slow):** `eps_spec` and `dE_peak` (§7
    metrics) at the §7 angle grid; test asserts reproducibility of the logged
    numbers (1%), not their size — the §7 support criterion applied to the
    numbers sets the documented supported angle range.
12. **Struct completeness (second-review refinement 4):** force each early
    return of `invz_solve_point_ordered` (paramagnetic `m0 < mtol`, MF
    non-convergence, sign-mismatch after retry) and assert the full declared
    field set is present each time.

## 8. Deferred follow-ups

1. **Azimuthal field + two-channel transverse MF.** Requires iterating
   `hx` AND `hy` for every direction (including `[Bx 0 0]`, where `B64s`
   induces `⟨Jy⟩ ≠ 0`), which changes the default x-field model and demands
   benchmark revalidation — exactly-C4 symmetry restored in exchange. Design
   options recorded in the review; do not resurrect the rejected `By ≠ 0`
   guard (angle-discontinuous, C4-inconsistent).
2. **Full-tensor A0(+A1)** (`odd_implementation_plan.html` Appendix A, ODD
   blocks zeroed): **A0** — rigorous tensor-RPA layer, `[12,12,nq]`
   Cartesian⊗sublattice `𝒥(q)` against the full 3×3 `χ0` from `invz_chi0z`;
   captures the `χ0_xz` cross-channel exactly. **A1** — projected 1/z bridge
   (`χ̃0 = χ_dom/(1+Σ_c) + χ_rest`), still a dominant-transition
   approximation, not a rigorous tensor 1/z. A0+A1 would supersede this
   stage's scalar routing for arbitrary tilt and subsume retardation.

## Review resolutions (field-angle-plan-review_by_Codex.md, 2026-07-16)

| Finding | Status | Resolution |
|---|---|---|
| 1 hy-guard vs C4 | Verified (`⟨Jy⟩ = −0.0690` at `[4 0 0]`, matches review) | Azimuth descoped; legacy x-only MF documented; two-channel MF deferred (§8.1) |
| 2 negative-Bz branch | Verified digit-for-digit | Sign-aware seed + sign assertion + mirrored retry (§2, §4); free-energy check in tests |
| 3 O(tilt²) overclaim | Accepted | Accuracy statement corrected; Σ=0 3×3 reference test added (§7); A0/A1 wording fixed |
| 4 map API ambiguity | Accepted | `opts.field_dir` + `S.field_dir/S.Bvec` (§5) |
| 5 overlay crash path | Verified (`one_field` line 113 outside try) | Explicit longitudinal `phase = 0` contract, no `invz_twolevel` fallthrough (§5) |
| 6 forced_moment semantics | Accepted | Non-circular 5-step acceptance; `si.mf_converged`; `pt.moment_branch` (§2, §4) |
| 7 labels/formatting | Verified (3 sites) | `\|B\| (T)` labels, vector-safe qpath error, crossover wording (§5-§7) |
| 8 helper/threshold rigor | Accepted | Precise `invz_field_vec` contract; `opts.Jyy0` dropped (moot after descope); single named `bz_tol` rule with dead-band zeroing (§1, §4) |
| 9 absent sources | Partial pushback | `LiReF4_MF_Yikai.m` exists outside the repo — absolute path now cited; `MF_chi.m`/`RPA.m` marked user-reported |

Second review (same file, top section), all accepted:

| Finding | Status | Resolution |
|---|---|---|
| SR1 free-energy formula | Verified digit-for-digit (F_mf −21.4696/−21.4766; naive comparison mis-ranks) | Variational `si.F_mf` + `si.E0` added (§2); test 3 corrected |
| SR2 option/dead-band propagation | Verified (`sopts` carries only `hyp/J0eff/Jxx0` today) | `opts.bz_tol` single-resolution rule + `opts.solve_opts` pass-through with reserved fields (§5-§6); pre-parfor dead-band normalization makes `S.Bvec` the used vectors; `invz_solve_auto` returns failed `pto` with valid `si`/`tl`; safe early-return field set (§4) |
| SR3 supported-angle metric | Accepted | `eps_spec` (spectral L2 with floor) + censoring-aware `dE_peak`; support = `eps_spec <= 5%` AND `dE_peak <= max(2%, eta)`; `demag = 0`, full 3×3 incl. `yz` (§7) |
| Refinements 1-5 | Accepted | `B = 6 T` continuity point (test 6); explicit mirror metric (test 4); map API test via `field_dir` (test 2); struct-completeness test 12; `MF_Chi.m` case + formula-authoritative citation |
