# invz_tensor — Full-Tensor 1/z Expansion (A0→A4)

Full-tensor (3 Cartesian × 4 sublattice, multi-level) implementation of Jensen's 1/z
effective-medium expansion for LiHoF4, in contrast to `../invz_projected/`, which projects
onto the scalar cc (Ising) channel with a two-level self-energy. The tensor route exists to
give the projected module's opt-in ODD (off-diagonal dipolar) extension (`invz_projected/README.html`
§1.9) a **rigorous completion**: a genuine multi-level, multi-Cartesian 1/z self-energy that
must reduce to the projected ODD equations E1/E4/E5 as a designed cross-validation
(framework `jensen_1z_framework.html` §11.8).

Implementation plan (the source of record for every design decision below):
`../docs/superpowers/plans/2026-07-17-invz-tensor-full.md` (v3). It supersedes
`../docs/superpowers/plans/2026-07-16-invz-tensor-odd.md` (SUPERSEDED banner atop that file);
see §"Status" below. Physics base: `../jensen_1z_framework.html` §§0–11 (scalar theory + ODD),
`../odd_implementation_plan.html` Appendix A (A0–A2 staging), `../three_level_1z_extension.html`
(A3 exact multi-level vertex), Yang–Wang PRB 10, 4714 (1974) & PRB 12, 1057 (1975) (SBO
diagrammatics). Full measurement record: `../docs/ODD-LOG.md` §A0–§A4.

**Status (2026-07-18): A0→A4 implemented and cross-validated.** All five stages below are
in place; the emergent-boundary test (framework §11.8) passes in two sectors (exact scalar
identity + matched-truncation ODD collapse). See "Headline results" and "Open items".

## Running the suites (quickstart)

**Prerequisites:** MATLAB **R2025a** (the module uses `pagemtimes`/`pagemldivide`); run every command **from the repository root** (the path contains spaces — keep the whole command quoted exactly as shown); **Python is not needed at run time** (the vertex oracle is the committed JSON fixture — MATLAB never calls Python).

**Fast check — the default gate (CORE, ~1–2 min warm-cache).** The tensor unit tests, run with `invz_projected` *off* the path:
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
```
→ expected **47 passed / 0 failed / 1 incomplete** (the incomplete is the `INVZ_SLOW`-gated A4 production-ladder test).

**Interop parity (optional) — live cross-checks against `invz_projected`:**
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```
→ expected **8 / 0 / 2**.

**Projected regression (must stay identical to the pre-tensor baseline):**
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
```
→ expected **143 / 0 / 19**.

**Full run including the slow gates** (`INVZ_SLOW=1`; adds the 16³/dpRng-30 reproducibility gate and the slow spectra/critical benchmarks — minutes; NB the A4 `e6` production rung is ~10 h if reached, see "Open items"). Prefix any of the above, e.g.:
```
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
```

The first cold run of a grid-dependent test recomputes the dipole/coupling sums into `invz_tensor/cache/` (gitignored; a few minutes at 16³/dpRng-30); subsequent runs are warm. Detailed path-isolation / boilerplate / grid-contract rationale is under **"Suite topology"** below.

**Run one computation (not a test)** — a mode-`'a1'` point solve, from the repo root:
```matlab
addpath(pwd); addpath('invz_common'); addpath('invz_tensor');   % root: MF_dipole/exchange/qVec_generator; shared engine; module
ion = invz_ion();
g   = invzt_qgrid(8, 'halfopen');                                % Γ-excluded 8³ grid
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
pt  = invzt_solve_point(ion, 2.0, [0.5 0 0], lat, struct('mode', 'a1', 'odd', true));
fprintf('crit = %+.4f (>0 = PM),  Sigma0 = %.4f,  converged = %d\n', pt.crit, pt.Sigma0, pt.converged);
```
Swap `'mode','a3'` for the genuine tensor self-energy, add `'nlevels','e3'` for an A4 rung, `'dress','dominant'` for the E1-matched truncation — see **"The mode ladder"**.

## Architecture — staged A0 → A4

| Stage | What it is | Key functions |
|---|---|---|
| **A0** | Tensor-RPA parity layer: cached `[12,12,nq]` coupling tensor, page-wise 12×12 RPA, tensor-owned anchors | `invzt_qgrid`, `invzt_jq_tensor`, `invzt_chi_rpa`, `invzt_gcc_lattice` |
| **A1** | Projected-1/z bridge: dominant-sector Σ inside the full tensor RPA (no `invz_chiperp`/`invz_odd_deltaJ` needed — the tensor RPA carries transverse/ODD mediation automatically) | `invzt_chi0_split`, `invzt_solve_point` (mode `'a1'`), `invzt_critical`, `invzt_tc_pm_extrap`, `invzt_chi_realaxis` |
| **A2** | Matrix effective medium: DIRECT closure `K = ctil^-1 − chibar^-1` (no iteration) | `invzt_emt_matrix`, `invzt_solve_point` (mode `'a2'`) |
| **A3** | Genuine tensor 1/z self-energy from the exact component-labelled four-point vertex | `invzt_kernels`, `invzt_vertex4`, `invzt_gamma4`, `invzt_vertex3`, `invzt_threestate`, `invzt_sigma_tensor`, `invzt_solve_point` (mode `'a3'`); oracle: `verify_tensor_vertex.py` |
| **A4** | Basis-defined state-space ladder: how the A3 physics evolves as CF/nuclear basis content grows | `invzt_rung_basis`, `invzt_run_ladder`, `invzt_report_ladder` |

## Module map

### `invz_tensor/` (this directory, `invzt_` prefix; all new code)

| File | Task | Responsibility |
|---|---|---|
| `invzt_qgrid.m` | T3 | THE q-grid builder — `'halfopen'` / `'legacy_inclusive'` conventions, weights, hash (LOCKED contract) |
| `invzt_jq_tensor.m` | T3 | Cached `[12,12,nq]` sublattice×Cartesian coupling tensor + `lat` struct assembly |
| `invzt_chi_rpa.m` | T4 | Page-wise 12×12 tensor RPA propagator at one frequency |
| `invzt_gcc_lattice.m` | T4 | Weighted BZ + sublattice average of the site-diagonal cc RPA response |
| `invzt_chi0_split.m` | T5 | Transition-mask split χ0 = χ_dom + χ_rest, with the LOCKED elastic convention |
| `invzt_solve_point.m` | T6/T9/T12/T13 | A1 bridge / A2 medium / A3 self-energy point solver (mode-switched); the module's main entry point |
| `invzt_critical.m` | T7 | PM-side critical transverse field Bc via `sign(crit)` bisection |
| `invzt_tc_pm_extrap.m` | T7 | Handle-based small-B-proxy Tc estimator (two-point crit extrapolation) |
| `invzt_chi_realaxis.m` | T8 | A1 scalar-Σ real-axis continuation ONLY (not A3 — see "Open items") |
| `invzt_emt_matrix.m` | T9 | A2 direct matrix effective-medium closure |
| `invzt_kernels.m` | T11 | KMS-scaled, exprel-stable divided-difference kernels (φ / I2 / I3) |
| `invzt_vertex4.m` | T11 | Component-labelled four-point vertex engine — DENSE ordered path-sum (production) |
| `invzt_gamma4.m` | T11/T12 | Vectorized connected four-point cumulant Γ4 over the whole (n,l) grid (~100× faster than row-by-row `invzt_vertex4` stage `'Gamma'`; bit-identical to ~1e-12) |
| `invzt_vertex3.m` | T11 | Component-labelled three-point vertex (sum-rule / condition-space route) |
| `invzt_threestate.m` | T12 | Explicit three-state toy single-ion model (independent-Δ1 contract, v3) |
| `invzt_sigma_tensor.m` | T12 | A3 self-energy assembly: vertex → Dyson resummation → A2 medium → outer fixed point |
| `invzt_rung_basis.m` | T13 | Basis-content-defined A4 rung → projector + basis energies (multiplet-complete) |
| `invzt_run_ladder.m` | T13 | A4 ladder driver (data-only, budget-refusing) |
| `invzt_report_ladder.m` | T13 | Serializes a completed A4 run into an ODD-LOG-pasteable text report (writes no files) |
| `invzt_active_projector.m` | T9/T12 (`/simplify`) | Shared frequency-consistent active-subspace projector (union-of-ranges rank-reveal); called by `invzt_emt_matrix` + `invzt_sigma_tensor` |
| `invzt_str.m` | `/simplify` | Shared value→display-string coercion for error messages (replaced 4 per-file copies) |
| `README.md` | T15 | This file |

### Shared engine and oracle

| Path | Responsibility |
|---|---|
| `../invz_common/` | Branch-shared single-ion engine (16 functions, pure `git mv` out of `invz_projected/`, zero logic changes): `invz_ion`, `invz_single_ion`, `invz_chi0z`, `invz_chiperp`, `invz_check_transverse_mf`, `invz_cfrot`, `invz_field_vec`, `stevens_ops`, `invz_matsubara`, `invz_twolevel`, `invz_g`, `invz_lambdas`, `invz_sigma`, `getf`, `invz_const`, `invz_is_gamma_equiv`, plus tensor-only `invz_cache_key.m`. `invz_projected/tests` remains the acceptance suite of record for these 16; any caller must `addpath` this folder explicitly. |
| `../verify_tensor_vertex.py` | Component-labelled vertex ORACLE (T10): mpmath first-principles reference of record for the A3 vertex machinery. Writes the committed `invz_tensor/tests/fixtures/vertex_oracle.json`; MATLAB suites never call Python at run time — the committed JSON IS the oracle. |

### Tests (`invz_tensor/tests/`)

| Path | Role |
|---|---|
| `invzt_anchors.m` | Plain fixture (not a test; `runtests` on a folder does not collect it) — tensor-owned pinned digit anchors, measured independently on this tree so CORE tests never depend on `invz_projected/tests` at run time |
| `test_invzt_jq_tensor.m` (also holds the q-grid convention tests), `test_invzt_chi_rpa.m`, `test_invzt_chi0_split.m` | A0/A1 building blocks |
| `test_invzt_solve_point.m`, `test_invzt_critical.m` | A1 bridge solver + critical finder |
| `test_invzt_chi_realaxis.m` | A1 real-axis continuation |
| `test_invzt_emt_matrix.m` | A2 matrix effective medium |
| `test_invzt_vertex.m` | T10/T11 — kernels + vertex4/vertex3 vs the mpmath oracle |
| `test_invzt_a3_threestate.m` | A3 three-state model + `invzt_sigma_tensor` (the framework §11.8 emergence gate) |
| `test_invzt_a4_ladder.m` | A4 ladder structural + (INVZ_SLOW-gated) production tests |
| `interop/test_invzt_jq_parity.m`, `interop/test_invzt_rpa_parity.m`, `interop/test_invzt_solve_parity.m`, `interop/test_invzt_critical_parity.m` | Live parity tests against `invz_projected` (optional suite — see "Suite topology") |
| `fixtures/projected_baseline.json` | T2 — frozen projected reference numbers (committed generator + filtered-clean dirty flag over physics-source paths, base git `22c9ce3`) |
| `fixtures/vertex_oracle.json` | T10 — the mpmath vertex oracle (138/138 checks) |
| `exploratory/explore_tensor_blocks.m`, `exploratory/freeze_projected_baseline.m`, `exploratory/perf_vertex_scaling.m` | Read-only measurement scripts (not tests); provenance for the anchors, the frozen baseline, and the T11 performance gate |

(The q-grid convention tests — `test_qgrid_conventions` etc. — are local functions inside `test_invzt_jq_tensor.m`, per the T3 plan; there is no separate `test_invzt_qgrid.m`.)

## The mode ladder

`invzt_solve_point(ion, T, B, lat, opts)` is the module's single point-solver entry; two
independent axes select what it computes:

**`opts.mode`** — the physics layer:
- `'a1'` (default) — A1 dominant-sector bridge: `ctil = cdom/(1+Σ) + crest` propagated through
  the full tensor RPA; Anderson-accelerated outer loop on the scalar Σ.
- `'a2'` — A2 direct matrix effective medium (`invzt_emt_matrix`); same outer loop, K is now a
  `[3,3,nwn]` matrix instead of a scalar.
- `'a3'` — A3 genuine tensor self-energy: `invzt_sigma_tensor` runs its own inner fixed point on
  `Vmat` (the four-point-vertex self-energy correction), consuming the A2 medium each iteration.

**`opts.nlevels`** — the single-ion construction (routes which Hilbert space the mode above acts on):
- `'std'` (default) — the full electronuclear `invz_single_ion` (136-dim: 17 electronic × 8 nuclear I=7/2).
- `'three'` — the explicit `invzt_threestate` toy (hyperfine excluded), the A3 validation model.
- `'eN'` / `'eNxI8'`, N ∈ {3, 6, 17} — A4 ladder rungs (`invzt_rung_basis`): the field/MF single-ion
  Hamiltonian built and diagonalized IN a multiplet-complete reduced CF basis; `xI8` tensors that
  electronic basis with the complete I=7/2 nuclear space. `'e17xI8'` = 136 = `'std'`.

**`opts.dress`** (A3 only) — which channels the vertex self-energy dresses:
- `'full'` (default) — dresses ALL Cartesian channels; includes the genuine beyond-E1 transverse-
  spectator dressing Jensen's dominant-only rule leaves bare (a reported physical correction, not a defect).
- `'dominant'` — E1-matched truncation: the vertex dresses ONLY the dominant cc transition
  (spectator-normalized, constraint 3), leaving the transverse spectator and non-dominant cc paths
  bare — the framework §11.8 emergence comparator (A3-dominant ≡ A1/E1 up to the O(1/z²)
  resummation-scheme difference, constraint 8).

The A4 ladder (`invzt_run_ladder`) walks the outer product of `{a1, a3-dominant, a3-full} ×
{odd off, odd on} × {'three','e3','e6', budget-refused: 'e17','e3xI8','e6xI8','e17xI8'}` at one
validation point, seed-continuous across the whole walk.

## Flag surface

All options use `getf(opts, name, default)` (silently defaulted, never required). Selected
flags, module-wide:

| Flag | Default | Meaning |
|---|---|---|
| `odd` | `true` | `false` zeroes the Cartesian-off-diagonal entries of every sublattice block of `lat.Jt` (a copy) — the cc sector no longer sees the transverse ODD blocks. |
| `chi_rest` | `true` | `false` drops the non-dominant (`crest`) part of the local χ0 — only the dominant (ground-manifold) sector is renormalized/mediated. |
| `Esplit` | `0.4653` meV | Dominant/rest split energy (half the ~0.9306 meV Γ-doublet-to-first-excited-CF gap); passed to `invzt_chi0_split` and (mode `'a3'`, `dress='dominant'`) to the dominant-group selector. |
| `chi0_diag` | `false` | TEST HOOK. Zeroes the cross-Cartesian elements of the local tensor `ctil`/`ctil0` before use, so with `odd=false` the cc sector EXACTLY decouples (enables exact-identity gates). |
| `force_sigma0` | `false` | (`invzt_chi_realaxis` only) TEST HOOK. Forces `Sigma_w = 0` identically, isolating the bare-propagator identity gate. |
| `mtol` | `0` | (`invzt_vertex4` only) Matrix-element pruning tolerance; `0` = no pruning (exact). Dropped `\|element\| < mtol`, discarded weight bounded in `out.pruned_bound`. |
| `rank_tol` | `1e-12` | Active-subspace threshold (relative to the largest active singular value) for every rank-revealing solve: `invzt_emt_matrix`'s common projector, `invzt_sigma_tensor`'s `resum_dyson`, the Hermitian-eig criticality square root. |
| `chiperp_scale` | `1` | (`invzt_threestate` only) Rescales the matched spectator coupling ρ → `chiperp_scale·ρ` and RE-REFINES (Δ1, m0) at that ρ so the doublet always still reproduces (Δ, M²) and only χ⊥ varies (∝ scale²); `chiperp_scale = 0` gives the exact two-level limit with `|3>` decoupled. |
| `far_excited` | `false` | (`invzt_threestate` only) `true` moves the spectator gap Δ2 from 0.9306 meV (LiHoF4 CF gap) to 40 meV (a "far" spectator regime probe). |
| `impl` | `'dense'` | (`invzt_vertex4` only) `'dense'` is the sole production engine. `'factored'` is DISABLED — errors `invzt:factoredUnproven` (Task-10 oracle check recorded `factored_ok = false`; see LOCKED conventions). |
| `dress` | `'full'` | (mode `'a3'` only) `'full'` \| `'dominant'` — see "The mode ladder" above. |
| `nlevels` | `'std'` | Single-ion construction rung — see "The mode ladder" above. |

**Note on `dd_tol`:** there is no exposed `dd_tol` option. The divided-difference kernels
(`invzt_kernels.m`, `invzt_gamma4.m`) use a HARDCODED small-argument branch threshold of `1e-10`
(matching the oracle generator verbatim: φ 1e-12, φ′ 1e-10, φ″ 1e-8, divided-difference
denominators 1e-10) so that repeated-node (degenerate) rows select the correct Hermite limit.
This is a kernel implementation constant, not a caller-tunable option.

Other opts of note (outer-loop control, all `getf`-defaulted, not part of the "flag surface"
proper): `mix_outer`, `tol_outer`, `max_outer`, `anderson_depth` (A1/A2: 5; A3: 1 — Anderson was
observed to converge to a spurious A3 fixed point, see `invzt_sigma_tensor.m`), `Sigma_seed` /
`Vmat_seed` (warm-start), `Ecut`, `hyp`, `transverse_mf`, `Lmax` (A3 internal Matsubara cutoff),
`static_idx` (`invzt_emt_matrix`, opt-in static-slot Hermiticity gate).

## Suite topology

Three suites, run from repo root (path contains spaces — always quote it):

**CORE** (after every task from T3 on; must pass with `invz_projected` absent from the path):
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
```
(`runtests` on a folder does not recurse into `interop/`.)

**INTEROP** (optional; live parity tests against `invz_projected`):
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```

**PROJECTED** (regression gate — must stay identical to the pre-`invz_tensor` frozen baseline):
```
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
```

Slow gates are prefixed `INVZ_SLOW=1` (tests use
`assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only')`) — this covers the A4
production ladder run (see "Open items").

CORE test boilerplate (no `invz_projected` on path):
```matlab
function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end
```
INTEROP suite adds two lines: `addpath(...'invz_projected')` and `addpath(...'invz_projected','tests')`
(for `invz_odd_anchors`).

**Grid contract (LOCKED, T3).** All grids come from `invzt_qgrid(n, conv)`, returning
`g.qvec` (Γ-excluded), `g.w` (uniform weights summing to 1, the dropped-Γ point simply unrepresented),
`g.conv`, `g.hash`. Two conventions, never mixed inside one comparison:
- `'halfopen'` (default for all new tensor work) — `ndgrid((0:n-1)/n)` minus Γ; n distinct points
  per axis, no duplicated reciprocal-periodic face.
- `'legacy_inclusive'` (parity with projected/legacy anchors ONLY) — exactly
  `qVec_generator(...,'grid',[n n n],'range',[-0.5 0.5])` + the same Γ filter (endpoint-inclusive,
  spacing 1/(n−1) — duplicates the ±0.5 faces by design, matching the projected production
  convention bit-for-bit).

Grid conv + hash thread into every cache key, the `lat` struct's provenance, and every reported
anchor.

**Cache namespace.** `invz_tensor/cache/` (`.gitignore = *`, never shared with
`invz_projected/cache`). Keys: `invz_cache_key('jqt1', dpRng, [qhash; convcode], pkey)` with
`pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]`; the loader `isequal`-verifies the
stored `pkey` + `qvec` before trusting a hit.

**Tier-2 scope box.** The tensor branch does NOT re-implement the projected Tier-2 additive
layer. A1 cross-validates against projected *Tier-1-only* numbers; variable-moments physics
re-enters ONLY at A3, via the mixed cumulants of the exact vertex, cross-validated against the
projected Tier-2 share (≈2.8%) as a REPORT (see the A4 headline below). Never stack an A3 result
with the projected Tier-1/Tier-2 corrections.

## LOCKED conventions

These bind every function in the module; a violation is a bug, not a judgment call.

1. **Tensor index order.** Composite index sublattice-major, Cartesian-minor:
   `i = 3*(s-1) + mu`, s = 1..4, mu = 1(a), 2(b), 3(c). Slices: cc = `Jt(3:3:12,3:3:12,iq)`,
   ca = `Jt(3:3:12,1:3:12,iq)`. Local response embeds as `kron(eye(4), chi0_3x3)`.
2. **K bookkeeping (framework §11.8, "bridge medium").** K is defined against the
   RENORMALIZED propagator, never the bare one: `K = 1/Gloc − 1/G̃0`, with
   `G̃0 = −[χ_dom/(1+Σ) + χ_rest]_cc` site-local. Defining K against bare G0 (the superseded
   2026-07-16 draft's error) is wrong at O(χ_rest·Σ) and breaks the exact reduction to the
   scalar closed form when weak transitions are dropped.
3. **Symmetric bracket (framework §11.8, "matrix resummation").** Beyond two levels the
   first-order correction is carried ADDITIVELY and resummed via the manifestly Hermitian
   sandwich `G̃0 = G0·(G0+𝒱)^-1·G0` — NEVER `inv(G0)` alone (a multi-pole G0 has real zeros
   between poles, where a ratio Σ = 𝒱/G0 acquires artificial poles). In the scalar two-level
   limit 𝒱 = G0·Σ and this collapses to `G̃0 = G0/(1+Σ)`.
4. **Negative-frequency tensor transpose relation (constraint 9).** From the Lehmann form,
   `χ_ρσ(−iωn) = χ_σρ(iωn)` — the TRANSPOSE, not componentwise evenness. The same relation
   holds for K and every channel-matrix object; scalar/diagonal cases degenerate to evenness.
   **Finding (A2, LOCKED):** χ0(iωn) is genuinely NON-HERMITIAN off the static slot (measured
   anti-Hermitian fraction ~6–17%, up to 9.75% at Bx = 0.5 T — the physical gyrotropic ~B
   response). K and χ_bar inherit the transpose relation to ~1e-17 and are Hermitian ONLY at
   the static wn = 0 slot; single-frequency Hermiticity is deliberately NOT asserted off that
   slot (it would contradict this constraint). Consumed AS-IS in A3 (`pt.Kmat`), never
   symmetrized.
5. **ΔTc decomposition language.** Sequential condition/Σ-space legs only
   (`uniform_only` / `qstruct_only` + closure defect); additive "(a)/(b)/(c)" language is
   FORBIDDEN (uniform-shift theorem, framework §11.3).
6. **Elastic convention (`invzt_chi0_split`, T5).** The dominant block's elastic term
   restricts the connected degenerate-pair sum to dominant-dominant pairs and replaces the
   global-mean counterterm with the DOMINANT-GROUP mean moment `Jdom = Σ_{p∈dom} P_p·M_pp`
   (not renormalized by `Σ_{p∈dom} P_p`) — a literal transcription of the locked-box
   convention; `crest` absorbs the residual by construction. Reported as
   `mspec.elastic_conv_share`, bounded by the excited-state population share.
7. **Criticality: Hermitian eigendecomposition, never `sqrtm` (v3).** With `ctil0` the static
   renormalized 3×3 and `C12 = kron(eye(4), ctil0)` PSD: `[U,D] = eig((C12+C12')/2)`;
   `d = real(diag(D))`, clipped below `rank_tol`; `S = U*diag(sqrt(max(d,0)))*U'`;
   `M = eye(12) − S*JtGamma*S`; `crit = min(real(eig((M+M')/2)))`. Avoids `sqrtm`'s
   complex/non-Hermitian noise near criticality; `crit > 0` in the PM phase; clipped mass +
   active rank recorded. **Phase key is `sign(pt.crit)`, never `pt.converged`** — the A1
   Anderson solver converges metastable PM fixed points inside the ordered phase.
8. **Three-state independent-Δ1 contract (v3, `invzt_threestate`).** A direct low-doublet
   tunnelling term `(Delta1/2)*Sx12` splits the doublet independently of the spectator
   coupling ρ (`Ja0`/`Jb0`); ρ carries ONLY the transverse/spectator channel, so ρ→0 leaves a
   GENUINE two-level system (splitting survives as Δ1) — this is what makes the exact
   scalar-compatibility identity (item 9 below) consistent. `(Delta1, m0, rho)` solve three
   independent targets `(Delta, M2, chi_perp)` by Newton on the full `eig(H3)` — no
   overconstraint.
9. **DENSE-primary vertex.** `invzt_vertex4`'s dense ordered path-sum is THE production engine
   and correctness reference (reproduces the mpmath oracle to ≤4e-15). The O(N³) `'factored'`
   contraction is an UNPROVEN conjecture: Task-10's oracle check recorded `factored_ok = false`
   — the naive resolvent chain drops the φ/KMS boundary + initial-state structure and
   mismatches dense by 1e15–1e17 (spurious poles). `opts.impl='factored'` is disabled and
   errors `invzt:factoredUnproven`; there is no O(N³) production path.

## Headline results (docs/ODD-LOG.md §A0–§A4; report, never tuned)

- **A0 (preflight):** tensor-owned anchors bit-identical to the superseded P0.3 anchors
  (confirms the `invz_common` move is behavior-neutral); on-axis small-q ODD decay confirmed
  LINEAR in q; frozen projected baseline captured (`fixtures/projected_baseline.json`).
- **A1 (bridge):** K bookkeeping verified to reduce EXACTLY to `invz_emt_scalar` in the
  decoupled limit (frozen no-ODD baseline `dSigma0 = −4.1e-5`, gate 2e-3). **A1 proxy-Tc
  (0.05 T, 16³/dpRng-30) = 1.5599 K** vs the grid-matched projected closed form 1.5442 K
  (+0.016 K A1 enhancement — retardation ~null + transverse RPA + `chi_rest`, same grid, so
  distinct from the 16³→∞ Richardson drop to 1.509 K). Real-axis continuation: exact
  γ-uniform identity to RelTol 1e-8; no-ODD peak-parity vs `invz_chi_realaxis` 0.0000 meV.
- **A2 (matrix medium):** `invzt_emt_matrix`'s scalar reduction matches `invz_emt_scalar`'s
  `med.K` to 7.8e-14 (positive sign, no flip); a2≡a1 exact identity on diagonal input holds
  EXACTLY at the K-map level.
- **A3 (genuine tensor self-energy) — the framework §11.8 emergence, SHARPENED.** The
  original single-number "slope ≥ 2.3" gate was internally inconsistent with LOCKED
  constraint 8 (Dyson-vs-additive resummation differ at O(1/z²), so a crit-shift slope ~2 is
  inherent, not a pass/fail threshold). Emergence is validated in **two sectors** instead:
  - **SCALAR (exact):** the A3 vertex reduces to `invz_sigma` at ρ→0 to **3.24e-11**
    (≤ 1e-8 gate) — the strongest single emergence statement.
  - **ODD (matched truncation):** `dress='dominant'` (E1's dominant-only rule; the transverse
    spectator is provably bare there) reduces to A1/E1 — dominant-dress ratio **rd(1) = 1.0159**
    vs full-A3 **rf(1) = 1.1132**; the beyond-E1 excess COLLAPSES 1.113 → 1.016 (86% removed),
    confirming the transverse-spectator dressing is what full-A3 adds on top of E1. The
    residual 1.6% (`|rd−1| = 0.0159 ≤ 0.05` band) is the constraint-8 O(1/z²) resummation
    ambiguity (`resum_spread_crit = −3.93e-2`). This SUPERSEDES the earlier "slope ≥ 2.3"
    framing.
  - Perf gate (dense, 12 h budget): AFFORDABLE `{three, e3, e6}` (e6 near the boundary at
    9.86 h); REFUSE `{e17 (~196 h), every ×I8, e17×I8 = 136 (~5.9e5 h)}`.
- **A4 (state-space ladder) — the basis-content climb.** At the anchor (2.0 K, 0.5 T), the
  beyond-E1 transverse-spectator-dressing share (`rf − 1`) CONVERGES as CF content grows:
  **three 11.3% → e3 11.4% → e6 2.6%**, landing on the projected Tier-2 share (≈2.8%) — the
  small toy models OVERESTIMATE the beyond-Gaussian content; a real 6-state CF basis screens
  it down to the projected value. The virtual-completeness deficit (rung χ0 vs the full-136
  χ0, a diagnostic, not a bound) drops monotonically **0.041 → 0.002**. `crit_shift_odd`
  (+ = ODD lowers Tc, the projected ΔTc > 0 direction): three/e3 +0.028/+0.031 (1.6 K, 0.1 T),
  e6 +0.0547 (anchor point). Budget (12 h, dense): affordable three/e3 (0.62 h), e6 (9.86 h);
  refuse e17 (196 h), e3×I8 (2524 h), e6×I8 (4e4 h), e17×I8 = 136 (5.9e5 h) — recorded in
  `out.skipped_rungs`, not run.

**Cross-validation summary:** the A3 ladder's beyond-E1 correction converging to the
projected Tier-2 2.8% share, from the opposite direction (a genuine tensor self-energy
gaining CF content, vs. the projected module's quenched-Gaussian internal-field dressing), is
an independent confirmation that the two routes describe the same physics once the tensor
side's basis is complete enough.

## Open items

- **The O(N³) "factored" vertex is a dead end, not a backlog item to retry as-is.** Task-10
  disproved the naive chained-resolvent conjecture (`factored_ok = false`). A cheaper-than-dense
  engine, if wanted, needs a different starting point — e.g. a transition-Liouville-space
  (superoperator) formulation that keeps the KMS/initial-state structure — not a patch of the
  current `'factored'` stage.
- **Production ladder beyond e6 is user-left.** e6 costs ~10 min/point but a full Tc scan at
  that rung is ~10 h wall time (`INVZ_SLOW`-gated, not run as part of routine verification).
- **e17 and every ×I8 rung are budget-refused**, not attempted: dense-vertex cost scales
  ~N⁴ and the smallest refused rung (e17) already projects to ~196 h under the 12 h gate.
- **A3 real-axis continuation is deferred.** `invzt_chi_realaxis` continues only the A1
  scalar-Σ; continuing the full `Vmat(iωn)` to the real axis needs either direct real-axis
  kernel evaluation or a fitted continuation (documented in the function header, not
  implemented).
- **A3 true-zero-field Tc is deferred.** Mode `'a1'`/`'a3'` both require a symmetry-breaking
  transverse field (`invzt:a1ZeroField` below 1e-6 T); every reported A3 Tc-like number
  (including the A4 ladder's) is the small-Bx proxy (0.05 T), never a true B=0 solve. True
  zero-field physics stays with the projected closed form (`invz_odd_zero_field`).
- **Deferred-Minor triage list.** Several non-blocking Minors were logged across A0–A4
  reviews (e.g. A0's filtered-flag root whitelist omitting `ellipsoid_demagn.m`, inert at
  `demag=0`) — see the Minors noted inline in each `docs/ODD-LOG.md` §A0–§A4 entry; none
  gate correctness, but a future pass should sweep them.

## See also

- Plan: `../docs/superpowers/plans/2026-07-17-invz-tensor-full.md` (every design decision,
  hard mathematical constraints 1–9, global constraints, per-task acceptance criteria).
- Measurement record: `../docs/ODD-LOG.md` §A0–§A4 (and §P0–§V4 for the projected-side ODD
  physics this module cross-validates against).
- Framework: `../jensen_1z_framework.html` §11 (ODD within 1/z) and §11.8 in particular
  (the emergent-boundary test this module implements).
- Projected module: `../invz_projected/README.html` §1.9 (Tier-1/Tier-2 ODD, the parity
  target for A1 and the Tier-2 cross-validation target for A3/A4).
