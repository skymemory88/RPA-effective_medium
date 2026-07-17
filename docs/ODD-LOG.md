# ODD-LOG — Off-Diagonal Dipolar (ODD) extension to `invz_projected/`

Running log for the ODD main-body implementation (`odd_implementation_plan.html`,
phases P0 → T1 → T2 → T3 → V4). Each phase appends a dated entry: what was
implemented/measured, test status, headline numbers. Physics numbers are
**reported, never tuned** — the only physics gates are internal identities and
flag-off regressions against published benchmarks.

Repo: `invZ expansion/` (path has spaces — quote it). Convention: energies meV;
`C = invz_const()`; χ = −G, χ₀ in meV⁻¹, `crit` dimensionless; ferromagnetic
J > 0; 4 sublattices; cc channel 4×4 per q.

---

## P0 — Preflight (P0.1–P0.5)

**Date:** 2026-07-16 · **Branch:** `invz-1z-lihof4` · **Git:** `360dfab`
(`360dfab` differs from the fast-baseline commit `5f4ff92` only by two plan `.md`
docs — no module code). **Phase:** READ-ONLY (no module code touched).
**Scripts:** `invz_projected/tests/exploratory/explore_chiperp.m` (P0.2),
`invz_projected/tests/exploratory/explore_odd_blocks.m` (P0.3). ODD lattice sums
at **dpRng = 30** (production default). MATLAB R2025a.

All P0 anchors are pinned full-precision in `invz_projected/tests/invz_odd_anchors.m`.

### P0.1 — Full dipole tensor upstream of the cc projection

- **`MF_dipole(q, N, a, tau[, geom])`** (repo root) returns `dip` `[3,3,4,4]` =
  `dip(mu, nu, s, s')`: Cartesian `mu,nu = 1..3` (a,b,c = x,y,z) × sublattice
  `s,s' = 1..4`. Units **Å⁻³, pre-gfac**. Element assembly (line 31):
  `dip(nn,mm,nt,mt) = -exp(-i q·r) · Tf(:,nn,mm)` with the dipole tensor factor
  `Tf(:,nn,mm) = 3 r_n r_m / r^5 − δ_nm / r^3` (geometry-only, gfac-free); pair
  Hermiticity `dip(:,:,mt,nt) = conj(dip(:,:,nt,mt))` (i.e. `J_ij = J_ji*`). The
  5th-arg `geom` (reciprocal basis + per-pair shifted `r` + the 9 tensor columns)
  is built once and reused across a q-sweep — bit-identical to the 4-arg call.
  Self-sublattice pair `{1,1}` excludes the central site via the `r2 < 0.01` cut.
- **cc extraction** — `invz_jq_modes.m` line 66:
  `Jcc = -squeeze(C.gfac*dip(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:))`.
  The ODD off-diagonal blocks tap the **same** tensor:
  `J^{cα}(q)_{ss'} = -C.gfac*dip(3,α,s,s')` (α: 1 = a, 2 = b). Exchange is
  Cartesian-diagonal (`ex(3,α)=0` for α≠3) → contributes nothing off-diagonal;
  Lorentz/demag are Cartesian-diagonal → no ODD shape terms. δJ is a pure
  dipole lattice sum.
- **`info.Jcc0` assembly** — at Γ (`dip0` from the priming call):
  `Jcc0d = -C.gfac*dip0(3,3,:,:) + lorz`, `lorz = 4*pi/(3*ion.Vc)*C.gfac`;
  uniform projection `info.Jcc0_dipole = v'·Jcc0d·v` (`v = [1 1 1 1]/2`);
  `info.Jcc0 = info.Jcc0_dipole + 4*ion.J12`. Measured (dpRng 30):
  `Jcc0_dipole = 6.8244e-3` meV (published 6.821e-3), `Jcc0 = 6.4244e-3` meV
  (published 6.421e-3; `4·J12 = -0.4e-3`). ✔
- **`grep -rn "Jxx0\|Jcc0\|MF_dipole" invz_projected/*.m` — coupling flow:**
  `MF_dipole` is called **only** in `invz_jq_modes.m` (priming line 60, per-q
  line 64) and `invz_jq_path.m` (Γ-limit line 68) — a single dipole source.
  `invz_jq_modes` → `info.Jcc0`/`Jnu`/`info.Jaa0` → `invz_bz_couplings` →
  `invz_run_phase_diagram` (`J0 = info.Jcc0`, `Jxx0`, `Jf`) → solvers
  (`invz_solve_point`/`_ordered` via `opts.J0eff`/`opts.Jxx0`), critical finders,
  spectra maps/paths. `ion.J0eff`/`ion.Jxx0` are **fallback constants only**; the
  authoritative live values are the derived `info.*` (dpRng-dependent). The ODD
  extension will tap `dip(3,1)`/`dip(3,2)` from the identical `MF_dipole` call.

### P0.2 — Units + transverse susceptibility χ⊥

- **χ⊥(1.53 K, [0 0 0], hyp, elastic on)**, symmetrized (a,b)×(a,b) block of
  `invz_chi0z` at z = 0 (default `transverse_mf = 'legacy_x'`):

  | element | value (meV⁻¹) |
  |---|---|
  | χ_aa | 17.63784561529863 |
  | χ_bb | 17.637845615298673 |
  | off-diag / antisym (max abs) | 6.82e-16 (machine zero) |
  | max abs imag(block) | 1.27e-15 |

  χ_aa = χ_bb to 4.3e-14 (C4 at Bx = 0); gyrotropic/antisymmetric part is
  machine-zero. **Elastic share of χ_aa = 6.39e-4** (0.064%) → χ⊥ is Van
  Vleck-dominated (motivates computing it ONCE, not self-consistently — T1.2).
- **Band verdict:** measured **17.638 meV⁻¹**. Expectation: full-CF band
  16–17 meV⁻¹; failure modes ≈ 11 (truncation) or ×2-off ≈ 34 (convention slip).
  17.638 sits at/just above the top of the band (~4% high, the legitimate full
  electronuclear value), **far from both failure modes → PASS, proceed.**
- **Dimensional closure (E1 needs NO extra g-factors):** `C.gfac = 0.08388`
  meV·Å³ carries `μ0/4π·(gL·μB)²`; `C.gfac*4/ion.Vc = 1.16544e-3` meV = `J_D`
  (matches the `invz_const` comment). `Jcc0_dipole = 6.8244e-3` meV cross-checks
  the 6.821 μeV anchor. Therefore `J^{cα} = -gfac·dip` is in **meV**, χ₀ in
  **meV⁻¹**, and `δJ = V·χ·V` is in **meV** — no extra g-factors. ✔
- **Bx = 0:6 T sweep at 0.31 K** (elastic on, hyp; all points MF-converged,
  residual < 1e-12):

  | Bx (T) | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
  |---|---|---|---|---|---|---|---|
  | χ_aa | 17.909 | 21.027 | 10.367 | 5.106 | 3.079 | 2.117 | 1.586 |
  | χ_bb | 17.909 | 20.226 | 18.891 | 15.888 | 13.683 | 12.147 | 11.046 |

  χ_aa (along the field a-axis) **peaks at Bx = 1 T then collapses** — the known
  (0.31 K, 1–2 T) island. Every point converges tightly, so the peak/collapse is
  **physical, not a convergence artifact**.
  **⚠ ESCALATION (Task 3):** the max relative step is **0.51** (Bx 1→2 T), so the
  plan's `test_smooth_along_Bx` gate `max(|diff|/|χ|) < 0.25` is **not
  achievable** on this curve. The pinned anchor is the measured truth; the Task 3
  smoothness threshold must be recalibrated to it (or restricted to Bx ≥ 2 T).

### P0.3 — ODD symmetry check (blocks `J^{cα}(q) = -gfac·dip(3,α,:,:)`, dpRng 30)

- **(i) On-axis rays** (max |element| of the 4×4 block, meV):

  | q | maxJca [q 0 0] | maxJcb [q 0 0] | maxJca [0 0 q] | maxJcb [0 0 q] |
  |---|---|---|---|---|
  | 1e-1 | 1.7881e-3 | 2.1e-17 | 1.7e-17 | 2.4e-17 |
  | 1e-2 | 1.8565e-4 | 8.7e-18 | 1.4e-17 | 1.3e-17 |
  | 1e-3 | 1.8591e-5 | 9.9e-18 | 1.3e-17 | 1.3e-17 |

  Along **[q 0 0]**, `J^{ca}` **decays LINEARLY** in q (ratio ≈ 10 per decade →
  0 at Γ; C2-about-c), `J^{cb}` ≈ 0 (machine). Along **[0 0 q]** both blocks are
  machine-zero at every q (exact symmetry). Blocks **do decay on-axis** — the
  physical requirement holds.
  **⚠ ESCALATION (Task 2 / Global Constraints):** because the decay is **linear**,
  the residual at q = 1e-3 is `maxJca = 1.859e-5` meV = **2.894e-3 · Jcc0**, which
  is **NOT ≤ 1e-6 · Jcc0**. The plan's assumed "smallest-|q| block norms ≤
  10⁻⁶·Jcc0" and Task 2's `verifyLessThan(m(3), 1e-6*Jcc0 + 5e-9)` are
  **unachievable at q = 1e-3** (would need q ≈ 3e-7, below MF_dipole resolution).
  The pinned `A.odd_onaxis_smallq.maxca` is truth; Task 2's absolute bound must be
  recalibrated to the measured linear-decay residual (e.g. `m(3) < 3e-3·Jcc0`, or
  gate the decade ratio `m(k)/m(k+1) ≈ 10`).
- **(ii) Tilted ray q·[1 0 1]/√2** (max |J^{ca}|, meV): 1.9672e-3 (q=1e-1),
  3.1885e-4 (q=1e-2), 1.3515e-5 (q=1e-3). The **q = 0.1 value ≈ 1.97 μeV ≈ 0.54 ×
  (4π·gfac/Vc = 3.661 μeV)** is the direction-dependent **macroscopic** magnitude
  (consistent with the plan's ~1.8 μeV tilted-ray bound). The apparent decay at
  smaller q is the **finite-cutoff resolution artifact** — `MF_dipole`'s
  real-space cutoff cannot resolve |q| below ~1/(dpRng·a) (documented in
  `invz_jq_path.m`), so the sum reverts toward the spherical average. **Bound:**
  small-q decay assertions are valid **on high-symmetry on-axis rays ONLY**;
  tilted rays carry a non-vanishing ~1–2 μeV where resolved. Never special-case
  q = 0 by grid extrapolation.
- **(iii) Generic q = [0.31 0.17 0.09]:** maxJca = **4.0909e-3 meV = 0.637 · Jcc0**,
  maxJcb = 2.394e-3 meV. Off-diagonal ODD coupling is an **O(Jcc0)** effect at
  generic q (vanishing only near high-symmetry directions).
- **(iv) Smallest 16³ shells:** `[1/16 0 0]` → maxJca = 1.142e-3 meV
  (**0.178 · Jcc0**); `[0 0 1/16]` → maxJca ≈ 1.1e-17 (machine 0). On the
  production grid the near-a*-axis shell carries ~0.18·Jcc0 of off-diagonal
  block; the c*-axis shell is zero.

### P0.4 — Cache audit + pre-parfor discipline

- **`invz_jq_modes` cache key:** `jq4_<dpRng>_<hash(qvec(:))>_<hash(pkey)>.mat`,
  `pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; demag; alpha; 4]`
  (trailing **4** = schema v4, which added `info.geomD`/`geomX`). Loader
  `isequal`-verifies **both** `pkey` and `qvec` and requires the `Juni` field;
  stale/legacy entries fall through, recompute, and overwrite. `hash_vec` =
  `<numel>v_<uint32 hex of single(sum(v.*(1:n)'))>`.
- **ODD namespace decision (Global Constraints):** ODD geometric blocks get their
  **own** namespace `odd1_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat` (schema tag
  `odd1`), `pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]`; one file
  per grid stores `Vca`,`Vcb`,`Vcc`,`info`,`pkey`,`qvec` (isequal-verified load).
  Existing `jq4_` caches are never touched. The ODD `opts.odd` path in
  `invz_jq_modes` will NOT write the `jq4_` cache (its `Jnu` depend on Xp).
- **Pre-parfor discipline (P0.4):** `invz_run_phase_diagram` builds ALL lattice
  sums **before** entering `parfor` (line 40): `[Jf, info, Jxx0] =
  invz_bz_couplings(ion)` and `Tc0 = invz_critical_T0field(...)` are computed once
  up front; workers receive broadcast scalars/arrays and do **no disk I/O**. The
  Task 5 ODD wiring must mirror this: build `[Vca,Vcb,Vcc,infoB] =
  invz_odd_blocks(ion, qvec, ...)` ONCE at the same pre-parfor point and thread it
  via `opts.odd_blocks` (solvers read blocks from opts, never the cache).

### P0.5 — Frozen benchmark baseline

- **Git:** `360dfab` (branch `invz-1z-lihof4`); working tree clean apart from the
  new P0 files. Module code is byte-identical to `5f4ff92`.
- **FAST suite** (cited from the controller's run at `5f4ff92`; log
  `.superpowers/sdd/baseline_fast_5f4ff92.log`):
  **109 Passed, 0 Failed, 12 Incomplete** · 28.03 s.
- **SLOW suite** (`INVZ_SLOW=1`, this run at `360dfab`; log
  `.superpowers/sdd/baseline_slow_360dfab.log`):

  **121 Passed, 0 Failed, 0 Incomplete** · 224.68 s · `SLOW_EXIT=0`.
  (`INVZ_SLOW=1` un-skips the 12 previously-Incomplete slow tests → 109 + 12 =
  121, all pass.)

  **Benchmark console lines (slow):** the slow tests do NOT print a benchmark
  table — they assert the published values internally via `verifyEqual`, all of
  which PASSED at `360dfab`. The only stdout beyond progress dots is two expected
  warnings, both from **passing** tests:

  ```
  [Warning: invz_peak_energy: wmin = 1 excludes (nearly) the entire w grid
   (max(w) = 0.6) -- every column censored to NaN. ...]   % tested behavior,
      test_invz_spectra_qpath>test_qpath_peak_censoring (line 85)
  [Warning: Mean field not converged after 800 iterations: |dmf| = 1.26e-06 meV] % near-converged,
      test_invz_tensor_ref>test_reproducibility_of_logged_table (line 71)
  ```

  Frozen benchmark targets (asserted internally, all green):

  | benchmark | target (AbsTol/RelTol) | source test |
  |---|---|---|
  | Σc(0) Richardson(12³,24³) | 0.3004 (±0.006); controller 2·S24−S12 = 0.2980 | `test_invz_sigma_crit` L47 |
  | Σc(0) Richardson(8³,4³) | 0.3447 (±0.001) | `test_invz_sigma_crit` L16 |
  | `info.Jcc0` | 6.421e-3 meV (RelTol 0.03) | `test_invz_sigma_crit` L48 |
  | Tc(B=0) (1/z) | 1.74 K (±0.08) | `test_invz_critical` L43 |
  | `pt.Sigma0` | 0.0932 (±0.02) | `test_invz_critical` L55 |

Every later task re-verifies its flag-off path against this frozen baseline
(counts may only rise on Passed; nothing may Fail).

---

## T1.1 — Geometric ODD lattice sums, cached (`invz_odd_blocks.m`)

**Date:** 2026-07-16 · **Branch:** `invz-1z-lihof4` · **Git (pre-commit):** `cf570de`
· **Phase:** T1 (first ODD code task). MATLAB R2025a.
**Created:** `invz_projected/invz_odd_blocks.m`,
`invz_projected/tests/test_invz_odd_blocks.m`.

### What was implemented

`[Vca, Vcb, Vcc, info] = invz_odd_blocks(ion, qvec, opts)` — the geometric
foundation every later ODD task consumes. Mirrors `invz_jq_modes`' geometry
priming (single `MF_dipole`/`exchange` build at q=[0 0 0], reused across the
q-sweep), per-q loop, Γ handling and cache mechanics.

- **`Vca`/`Vcb`** `[4,4,nq]` complex = `-C.gfac*dip(3,1|2,:,:)` per q —
  **dipole-only** (exchange is Cartesian-diagonal, Lorentz/demag diagonal → no ODD
  shape terms). **NOT Hermitized** (`J^{ca}(q)` is not Hermitian); pair identity
  `J^{ca}(-q) = conj(J^{ca}(q))` is a test gate.
- **`Vcc`** `[4,4,nq]` = the SAME assembly `invz_jq_modes` eigendecomposes,
  `-C.gfac*dip(3,3,:,:) + sign(ion.J12)*ex(3,3,:,:)` + `lorz` at Γ-equivalent q,
  Hermitized per q. Verified bit-for-bit against `invz_jq_modes` by the parity
  test (eigenvalues to AbsTol 1e-12; `info.Jcc0`/`Jaa0` to RelTol 1e-12).
- **`info`**: `dpRng`, `Jcc0`, `Jaa0`, `Jcc0_dipole`, `Jaa0_dipole`, `lorz`.
- **Demag guard**: `error('invz:oddDemag', ...)` if `ion.demag ~= 0` (intrinsic-only
  layer; demag stays in `invz_jq_modes`).
- **Cache**: own namespace `odd1_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat`,
  `pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]` (25 elems), one file
  stores `Vca,Vcb,Vcc,info,pkey,qvec`; loader `isequal`-verifies **both** pkey and
  qvec. Never touches `jq4_`.

### Test status (TDD)

- **RED** (undefined function): 4 fast tests errored; `test_ds2023_geometry_sums`
  (pure-geometry unit guard, calls `MF_dipole` directly) **passed even at RED** —
  DS2023 Suppl. Table I sums land on 36.73/17.93 · a⁻⁶ (a=5.175 Å) to RelTol 1%
  with **no systematic factor** (the `{s,1}` lower-triangular indexing for s=1..4
  covers all four sublattice partners, self-site excluded by the `r2<0.01` cut).
- **GREEN**: focused file **5 Passed / 0 Failed** (fast); slow-gated Parseval
  Incomplete when unset. Full fast suite **114 Passed / 0 Failed / 13 Incomplete**
  (baseline 109/0/12 + 5 new fast + 1 slow-gated). `INVZ_SLOW=1` on this file only:
  **6 Passed / 0 Failed / 0 Incomplete**.

### T1.1(iv) Parseval (slow test) — headline numbers

n=8 uniform mesh (incl. Γ), dpRng 20, row s=1:

| quantity | value (meV²) |
|---|---|
| lhs = ⟨Σ_{s'}|Vca(1,s')|²⟩_BZ | 1.3415735084e-05 |
| rhs = Σ_s Σ (gfac·Tf(:,3,1))² (real-space) | 1.3455701810e-05 |
| abs ref = gfac²·36.73/5.175⁶ | 1.3454750450e-05 |

- **residual |lhs−rhs|/rhs = 2.97e-03** (tol 1e-2) ✔
- **|lhs−abs_ref|/abs_ref = 2.90e-03** (tol 1.5e-2) ✔ — validates gfac placement.

The ~0.3% residual is the expected n=8 superlattice-folding term (r⁻⁶-suppressed),
well inside tolerance. First block build (cache write) = 1.43 s.

### Cache evidence

`odd1_*.mat` files created under `invz_projected/cache/`:
`odd1_10_6v_4005c28f_25v_45cc976e.mat` (1.4 kB, from the cache round-trip test),
`odd1_20_1536v_49041200_25v_45cc976e.mat` (139 kB, from the Parseval mesh). Bitwise
round-trip verified; a J12 change (distinct pkey) misses and recomputes.

### Concern (escalated, not blocking)

The brief's `test_cache_roundtrip_selfverifying` asserted
`verifyFalse(isequal(V3a, V1a))` after a 5% `ion.J12` bump. This is **unsatisfiable**:
`Vca`/`Vcb` are dipole-only per the interface spec, hence J12-independent —
empirically `isequal(V3a,V1a)=1`, maxdiff = 0 (Vcc changes, maxdiff = 1e-5). The
assertion contradicts the brief's own interface. Minimal intent-preserving fix
applied (own new test file only): assert on the **cc** block, which carries the
exchange (|J12|, sign J12) and is the observable that proves the pkey-miss +
recompute. Documented inline in the test. No module code or pinned anchor touched.

---

## T1.3 — ODD-mediated coupling δJ^cc(q): E1/E4/E5 + `invz_jq_modes` `opts.odd`

**Date:** 2026-07-16 · **Branch:** `invz-1z-lihof4` · **Git (pre-commit):** `dbdf0d6`
· **Phase:** T1. MATLAB R2025a. χ⊥ = `invz_chiperp(ion, 1.53, [0 0 0])`
(diag 17.63784562 meV⁻¹, the pinned P0.2 anchor).
**Created:** `invz_projected/invz_odd_deltaJ.m`.
**Modified:** `invz_projected/invz_jq_modes.m` (strictly additive `opts.odd`
branch; `git diff` shows 116 insertions, **0 deletions** — the default path is
byte-untouched), `invz_projected/tests/test_invz_odd_blocks.m` (+4 tests).

### What was implemented

- **`[dJ, d, dinfo] = invz_odd_deltaJ(Vca, Vcb, Xp)`** — E1 contraction per q
  (`dJpre = Vca·Xp₁₁·Vca' + Vca·Xp₁₂·Vcb' + Vcb·Xp₂₁·Vca' + Vcb·Xp₂₂·Vcb'`,
  Hermitized; PSD by construction since `dJpre = [Vca Vcb](Xp⊗I₄)[Vca Vcb]'`),
  E4 self-site subtraction (per-sublattice BZ-mean removed from the **diagonal
  only**), E5 uniform reduction `d = mean_s mean_q dJpre(s,s,q)` with the S4
  assert (four per-sublattice means to 1e-10 rel, `invz:oddS4`) and a
  machine-real assert on the diagonal before `real()` (`invz:oddImag`). Header
  carries the **no-double-counting comment trail** (plan T1.3/§8): the
  subtracted on-site constant multiplies (σᶻ)² = 1 in the strict two-level
  limit — pure energy shift, no Tc content; its physical residue (internal
  field renormalizing Δ, M², n01) is Tier 2's; grid matrices carry
  post-subtraction δJ, `Jcc0` carries the explicit −d, exactly once. Caller
  contract: qvec must be a full uniform BZ mesh; Γ-less driver grids fine
  (O(1/Nq) difference); never q-paths.
- **`invz_jq_modes` `opts.odd = false (default) | struct('Xp', [2×2])`** — the
  odd branch diverts right after dpRng/cache resolution into a new local
  function `jq_modes_odd`: Xp validated ([2,2] real symmetric finite, error id
  `invz:oddXp`), demag ≠ 0 rejected (`invz:oddDemag`, intrinsic-only), blocks
  from `invz_odd_blocks` (odd1_ cache), per q `M = Vcc + dJ` (Hermitized) →
  `Jnu = sort(real(eig(M)))`, `Juni = v'Mv`. **`info.Jcc0 = infoB.Jcc0 − d`**
  (E5 — THE line carrying the DS2023 MF mechanism downstream).
  `info.Jshape_cc = 0` (demag = 0 enforced); `info.geomD/geomX` provided by a
  one-off priming rebuild (~0.16 s at dpRng 30; `invz_odd_blocks` doesn't
  export geom) so the info contract stays whole for `invz_jq_path`-style
  consumers; `info.odd = struct(d, dJ_mean_diag [4,1], dJ_max,
  uniform_residual, Xp)`. The ODD path **never reads or writes `jq4_`**
  (verified two ways: code inspection — the branch returns before the jq4_
  cache block — and empirically: 0 new `jq4_*.mat` after odd calls on both
  grids below, `cache = true` on 16³).

### Test status (TDD)

- **RED:** `test_deltaJ_identities` (undefined `invz_odd_deltaJ`),
  `test_jq_modes_odd_zero_Xp_equals_off` + `test_jq_modes_odd_on_structure`
  (no `info.odd` field) errored; `test_jq_modes_odd_off_bitwise` **passed at
  RED** as designed (it is the regression gate — `opts.odd = false` was a
  no-op before and after).
- **GREEN:** focused file **9 Passed / 0 Failed / 1 Incomplete** (slow-gated
  Parseval). Full fast suite **122 Passed / 0 Failed / 13 Incomplete** (29.8 s;
  baseline 118/0/13 + the 4 new tests, nothing fails). Bitwise regression
  `isequaln({Jnu,info,Juni})` on the default path: PASS.

### Headline numbers (report, not tuned)

χ⊥(1.53 K, 0 T), Γ-less uniform meshes `(0:n−1)/n`:

| quantity | test 6³ / dpRng 15 | production 16³ / dpRng 30 |
|---|---|---|
| **d** | **4.7352243883e-04 meV = 0.4735 μeV** | **4.7461418306e-04 meV = 0.4746 μeV** |
| d per sublattice (rel spread) | 1.1e-16 | 0 (bit-identical) |
| Jcc0 → Jcc0 − d | 6.4291e-3 → 5.9556e-3 meV | 6.4244e-3 → 5.9498e-3 meV |
| d / Jcc0 | 7.37% | 7.39% |
| dJ_max (pre-sub) | 6.189e-4 meV | 6.180e-4 meV |
| presub_min_eig | +3.9e-35 (PSD exact) | +7.9e-35 |
| postsub_diag_bzavg | 1.0e-19 | 2.5e-19 |
| uniform_residual (smallest shell, max over shell) | 1.541e-4 meV = 2.40e-2·Jcc0 (q = 1/6) | 2.633e-5 meV = 4.10e-3·Jcc0 (q = 1/16) |
| blocks build (cache warm) | 0.4 s | 32.7 s (one-off; `odd1_30_12285v_4c0ca1e0_25v_45cc976e.mat`) |

- **Expected-band verdict: PASS.** d = 0.474–0.475 μeV sits inside (at the top
  of) the 0.3–0.5 μeV band, grid- and dpRng-stable to 0.23% between the two
  grids. Cross-check: Parseval estimate d ≈ χ_aa·⟨Σ|Vca|²⟩ + χ_bb·⟨Σ|Vcb|²⟩ =
  17.638 × 2 × 1.3455e-5 = 4.747e-4 meV — matches to 0.03%. The uniform
  coupling drops 7.4%, the right neighborhood for DS2023's 5–8% MF Tc effect.
- **uniform_residual vs the plan's ≤ 1e-4·Jcc0 claim:** holds **on-axis c\***
  exactly as P0.3 predicted — the c*-shell E1 value is machine zero
  (2.5e-33 meV = 3.9e-31·Jcc0). The reported `uniform_residual` takes the max
  over the whole smallest-|q| shell, whose a*/b* members carry the
  **linear-in-q** ODD residual (P0.3 recalibration): 4.10e-3·Jcc0 at q = 1/16,
  2.40e-2·Jcc0 at q = 1/6 — grid geometry, not a symmetry violation. The plan's
  1e-4 figure is met only on the c* axis; off-axis shells follow the P0.3
  linear-decay law (≈ χ⊥·|Jca(q)|², |Jca| ∝ q).

### Finding (baseline property, not ODD-induced — flagged for T1.5)

The T1.3 acceptance-(v) ordering-margin gate `max Jnu < Jcc0 − d` **passes on
the 6³ test grid** (margin +3.665e-4 meV odd; +4.279e-4 off). On the
**production 16³ grid it is violated ALREADY WITH ODD OFF**: the top branch at
the smallest **transverse** shell q = [0 1/16 0] (= max Juni there too) is
6.45263e-3 meV > Jcc0 = 6.42444e-3 (excess +2.82e-5 meV, +0.44%), 4/4095
q-points above. With ODD on the picture is unchanged in kind: the same 4
points, excess +3.71e-5 meV vs the shifted Jcc0 − d. This is the known
dipolar q→0 transverse limit (q̂_c = 0 ⇒ no macroscopic suppression; J(q→0⊥)
approaches the Lorentz-included Γ value, and the O(q²) lattice correction is
positive along a*/b*) — the c* shell sits far below (−3.22e-3 meV). ODD
neither creates nor spreads the excess (count 4 → 4; it shrinks slightly in
relative weight against the finite-q gain). **Consequence:** never gate
`max(Jnu) < Jcc0` on production-fine Γ-less grids; the physical ordering
comparison at Γ is `info.Jcc0`/`Juni`-based, and the near-Γ transverse shell
rides the discontinuous macroscopic limit. T1.5's criticality measurements
should confirm the RPA instability still lands at the uniform Γ mode through
the 1/z machinery (which consumes `info.Jcc0`, not grid `max(eig)`).

---

## T1.4 — ODD wiring: point solvers, critical finders, T0field handles, driver switch

**Date:** 2026-07-17 · **Task:** main-body plan Task 5 (T1.4) · **Commit:** (this commit)

### What was implemented

- **`invz_solve_point` / `invz_solve_point_ordered` `opts.odd`** (false default):
  requires `opts.odd_blocks = struct(Vca,Vcb,Vcc,Jcc0-UNSHIFTED)` from
  `invz_odd_blocks` AND `Jnu_flat = []` (both guarded, `invz:oddArgs`). When on:
  `Xp = invz_chiperp(ion,T,Bx, {hyp,Jxx0,transverse_mf})` with the SAME resolved
  option values as the solver's own single-ion call (T1.2 same-converged-state);
  `[dJ,d] = invz_odd_deltaJ`; modes rebuilt per q from `Vcc + dJ` (Hermitize,
  `sort(real(eig(·)))`); **`J0eff ← J0eff − d` applied exactly once in the
  solver** (T1.3 bookkeeping: grid diagonal already carries −d via E4, J0eff the
  explicit −d via E5, no other q = 0 handling; callers keep passing UNSHIFTED
  `info.Jcc0`). In the ordered solver the shifted J0eff feeds BOTH `siopts.J0z`
  and `pt.crit` (single shift, two uses — it is the physical uniform coupling).
  `pt.odd = struct(d, Xp)` on full solve paths. Paramagnetic solver reuses
  `chiperp` `info.si` for its own si (bit-identical option set — saves one
  136-dim diagonalization; ordered solver does NOT reuse, it needs the
  ordering-MF state; χ⊥ stays paramagnetic-MF by T1.2 design).
- **`invz_crit_at`**: opts/Jf already forward verbatim → zero threading changes.
  One safety addition: `invz:odd*` errors are STRUCTURAL misuse and now rethrow
  from the classifying catch instead of being absorbed as an "ordered-phase"
  verdict (can only fire with the flag on; flag-off catch behavior unchanged).
- **`invz_critical`**: untouched (opts reach the solver unfiltered via crit_at).
- **`invz_critical_T`**: `invz:oddTc0` guard at the top of `adaptive_anchor`
  (after the explicit-`opts.Tc0` early return, before any solve) — with
  `opts.odd` on the adaptive window must not silently anchor at the no-ODD Tc0.
  Explicit `opts.window` or `opts.Tc0` paths never hit the guard.
- **`invz_critical_T0field(ion, Sc, J0eff)`**: `Sc`/`J0eff` now numeric OR
  `function_handle` of T; numerics wrapped as constant handles → identical
  arithmetic/floats/bisection path (exact-equality regression-gated).
- **`invz_ion`**: `ion.odd = 0` documented default (drivers read ion.odd,
  libraries read opts.odd; intrinsic-only, demag = 0).
- **`invz_run_phase_diagram`**: `ion.odd` branch builds the blocks ONCE
  pre-parfor on the same 16³ Γ-less grid (P0.4); finder calls get
  `Jnu_flat = []` + `odd/odd_blocks` + UNSHIFTED `J0eff`. **Anchor choice:
  inline-handles route** — Tc0 from the generalized `invz_critical_T0field`
  with script-local `J0T(T) = J0 − d(T)`, `ScT(T) = Σc(J0 − d(T), modes(Vcc + dJ(T)))`
  (the `invz_odd_zero_field` mode-'full' algebra; clearly marked SEAM to be
  replaced when T1.5 lands). Known 16³ `invz:sigmaCritExcluded` wrinkle
  (ODD-LOG §T1.3, pre-existing without ODD) is surfaced ONCE by a probe call,
  then silenced for the ~31 bisection repeats (warning state restored after).

### Test status (TDD)

- **RED:** `test_invz_odd_solve.m` — args-guard, pm-point, T0field-handles and
  oddTc0 tests failed/errored as expected; `test_solve_point_flag_off_bitwise`
  passed at RED by design (it is the regression gate: `odd = false` was an
  ignored unknown opt before, and must stay a no-op path after).
- **GREEN:** focused file **5 Passed / 0 Failed**. Full fast suite
  **127 Passed / 0 Failed / 13 Incomplete** (30.3 s; baseline 122/0/13 + 5 new,
  nothing fails). File-scoped heavy regressions (brief Step 4):
  `test_invz_solve_point` 5/0/0 and `test_invz_critical` 0/0/5
  (all slow-gated/filtered) — identical outcomes to baseline.
- Smoke (non-suite): ordered-solver odd path at (1.30 K, 0 T) converges,
  m0 5.0948 → 4.9373 (ODD weakens the moment), d = 0.475 μeV, guards fire;
  driver `checkcode` clean (one pre-existing parfor-broadcast advisory).

### Headline numbers (report, not tuned; 8³/dpRng-15 test grid, cache warm)

- **Overhead** (pm-point gate, 1.80 K, 0.1 T): ODD-on solve 0.08 s vs 0.06 s
  off — gate `t_on < 1.2·t_off + 0.5 s` passes with wide margin (the si reuse
  leaves only the χ⊥ chi0z evaluations + δJ contraction + 511 4×4 eigs as
  overhead, ~0.02 s).
- **Sign contract** at the PM gate point: crit 0.0249 → 0.0967 (increases),
  d = 0.473 μeV.
- **Plan report point (1.60 K, 0.1 T), both flags, no gate:** ODD-off does NOT
  converge (ordered side, crit = −0.2339); **ODD-on converges paramagnetic,
  crit = +0.0342** — i.e. on this test grid the ODD shift moves the (1.60 K,
  0.1 T) point across the boundary, the T1.4 Tc-suppression signal in raw form.
- **Anchor seam check** (same grid, handle route): Tc0 = 1.6736 K (off) →
  1.4351 K (odd), ΔTc = 0.238 K, d(Tc0_odd) = 0.4746 μeV. Production-grid
  zero-field numbers (12³/24³ Richardson) are T1.5's measurement — not run here.

---

## T1.5 — Zero-field ODD measurement engine (`invz_odd_zero_field.m`) + T1.6 README §1.9

**Date:** 2026-07-17 · **Branch:** `invz-1z-lihof4` · MATLAB R2025a. Production grids
**12³/24³, dpRng 30**, Σc-benchmark generator mesh (`qVec_generator` range `[-0.5 0.5]`,
Γ-excluded — `test_invz_sigma_crit.m` lines 41–42), linear O(1/n) Richardson
`X_rich = 2·X₂ − X₁` (line 46). Physics numbers **reported, never tuned**.

### What was implemented

- `invz_odd_zero_field(ion, opts)` — grids `{12,24}`, dpRng 30, cache, mode ∈
  `off | full | uniform_only | qstruct_only`. Blocks built ONCE per grid outside every
  root find (P0.4); χ⊥(T)/δJ(T)/d(T) rebuilt inside via memoized (Sc,J0) handles fed to
  the generalized `invz_critical_T0field`. Excluded-mode counts taken explicitly (`nex`),
  `invz:sigmaCritExcluded` suppressed with the error-safe house pattern.
- **SEAM replacement:** the T1.4 inline-handle block in `invz_run_phase_diagram.m` (and its
  `odd_d_at`/`odd_Sc_at` helpers) is replaced by
  `Tc0 = invz_odd_zero_field(ion, struct('mode','full','grids',{{16}},'dpRng',30,'cache',true))`
  — same mode-'full' governing algebra, on the driver's single 16³ grid (not the {12,24}
  Richardson value), so the adaptive-window anchor matches the parfor solves' mesh. Flag-off
  path byte-identical; the engine's internal warning suppression makes the wrapper moot.

### Decomposition — controller adjudication (2026-07-17), after the Task-6 BLOCKED round

The source-plan ΔTc-space split was ill-posed (measured on 12³/24³ Richardson):

- **c1≡a THEOREM (E4/E5 validation).** Literal (c) = modes(Vcc+δJ+d·I), J0=Jcc0 is a
  *simultaneous uniform +d shift* of (a)'s couplings AND J0; R2007 criticality depends only
  on J0−Jν (invariant) plus a d-order numerator shift, so **Tc_c1 ≡ Tc_a to numerical
  precision** (1.50937 both, |Δ|<1e-14). This is a closed-form theorem: it independently
  validates that the E4 self-site subtraction and the E5 −d reduction are the *same*
  constant d, applied once — the bookkeeping is consistent.
- **naive-(b) invalid regime (DS2023's naive-MF inconsistency).** Source-plan (b) =
  modes(Vcc), J0=Jcc0−d drives J0 below the **peak finite-q Vcc mode**: at 24³,
  max(eig Vcc)=0.006143 > Jcc0−d=0.005941, so 96 modes are excluded (`nex.b_naive` = 0/96
  at 12³/24³) and the "uniform mode is critical" assumption is violated (ordering would move
  to finite q). This is exactly DS2023's naive-MF limit. The full-ODD (a) avoids it because
  δJ (from the same χ⊥) **lowers the peak**, 0.006143 → 0.005727 at 24³, keeping the uniform
  mode critical. Kept as a reported diagnostic with its exclusion counts.
- **GOVERNING split = sequential condition/Σ-space factorial** (neither leg enters the
  invalid regime): (b) condition-level `(Jcc0−d)·χ0cc = 1+Sc_off` with `Sc_off` FROZEN at the
  no-ODD value (no exclusions possible); (c) Σ-level `Jcc0·χ0cc = 1+Sc_odd(T)` with
  `Sc_odd(T)=Σc(Jcc0−d, modes(Vcc+δJ))` (the full ODD Σc at its own consistent config).
  closure_defect = (a) − [(b)+(c) − off] is REPORTED, never gated.

### Headline numbers (Richardson 12³/24³, dpRng 30 — report, not tuned)

| variant | J0 | Σ source | Tc (K) | ΔTc from off | nex 12/24 |
|---|---|---|---|---|---|
| off (published route) | Jcc0 | Σc(Jcc0, Vcc) | **1.74347** | — | 0 / 0 |
| **(a) full** | Jcc0−d | Σc(Jcc0−d, Vcc+δJ) | **1.50937** | **0.2341** | 0 / 0 |
| (b) condition-level | Jcc0−d | **Sc_off (frozen)** | 1.61514 | 0.1283 | 0 / 0 |
| (c) Σ-level | Jcc0 | Sc_odd | 1.62968 | 0.1138 | 0 / 0 |
| — b_naive (diagnostic) | Jcc0−d | Σc(Jcc0−d, Vcc) | 1.28269 | 0.4608 | 0 / **96** |
| — c1_literal (theorem) | Jcc0 | Σc(Jcc0, Vcc+δJ+d·I) | 1.50937 | 0.2341 | 0 / 0 |
| — c_factorial (diagnostic) | Jcc0 | Σc(Jcc0, Vcc+δJ) | 1.79430 | −0.0508 | 0 / 0 |

- **Off route reproduces the published anchors:** Σc(0) 12³/24³ = 0.26388/0.28093 →
  Richardson **0.29798** (published 0.3004, target 0.2980, AbsTol 0.006 ✓); Tc per-grid
  1.79094/1.76720 → Richardson **1.74347** (target 1.743, AbsTol 0.006 ✓).
- **ΔTc(0) = 0.234 K (13.4%): Tc(0) 1.743 K → 1.509 K.** per-grid Tc(a) 1.55505/1.53221.
- **ΔΣc(0) = +0.0908:** Σc(0) 0.29798 (off) → 0.38880 (ODD).
- **d(Tc) = 0.483 μeV** (12³ 0.4914, 24³-at-Tc_rich 0.4832) — inside the 0.3–0.5 μeV report
  band; χ⊥-flat ⇒ d(T) nearly T-independent. 7.5% of Jcc0_dipole = 6.42 μeV.
- **GOVERNING split closes:** closure_defect = **+0.00802 K** — only 3.4% of the 0.234 K
  effect (source-plan literal split gave 0.461 K; the governing legs are additive to ~3%).
  Directionality: (a),(b),(c) all < off; (a) the largest suppression.
- **Physics context (report, models differ — plan §8):** the no-ODD 1/z baseline missed
  ≈0.21 K vs experiment (1.74 → exp 1.53 K). Tier-1 ODD alone delivers 0.234 K, landing at
  1.509 K — closing essentially the whole gap (a hair below 1.53 K). This is **more**
  suppression than DS2023's 3-state MF (5%, 2.27→2.18 K) because the 1/z fluctuation channel
  amplifies the ODD coupling's q-structure; the comparison is qualitative (their spin-½ MC
  Tc = 1.6295 K; different CF treatment, Gaussian-vs-non-Gaussian statistics). **Never** import
  their 0.805 longitudinal rescaling or perturbative hyperfine (plan §8).

### Excluded-mode reconciliation (controller ask — characterize, don't chase)

Two apparently-conflicting numbers, resolved by mesh convention:

- **T1.3 ledger:** on the **16³ ndgrid** mesh (`(0:n-1)/n`, Γ-dropped, 4095 pts), the top
  Vcc mode exceeds Jcc0 by **+2.82e-5 at 4/4095 points** (axis-aligned near-Γ shell).
- **T1.5 here:** on the **24³ Σc-benchmark generator** mesh (`linspace(-0.5,0.5,24)`, no Γ
  node, 13824 pts), **max(eig Vcc) = 0.006143 < Jcc0 = 0.006424** — off-route nex = 0.

The cc mode's q→0 limit is **direction-dependent** (non-analytic macroscopic dipolar term
∝ q̂_c q̂_c, the plan's "q=0 pitfall"). The two meshes sample different smallest-|q|
directions: the ndgrid axis-aligned shell hits a direction where the top mode marginally
overshoots Jcc0; the generator's offset (tilted, Γ-free) shell stays just below it. Neither
affects the governing legs: (b) freezes Σc at the consistent no-ODD config, and (a)/(c) use
δJ, which lowers the peak below Jcc0−d. The naive-b 96 exclusions at 24³ are the *reduced* J0
(Jcc0−d = 0.005941) dropping below the same 0.006143 top Vcc mode — a different threshold
than Jcc0, and exactly the DS2023 inconsistency the diagnostic exists to expose.

### Bc(T) boundary shifts (16³ ndgrid, dpRng 30; `invz_critical`, window [0.1 7] T)

| T (K) | Bc off (T) | Bc on (T) | shift (T) |
|---|---|---|---|
| 0.31 | 4.0552 | 3.9187 | **−0.1365** |
| 0.80 | 3.3687 | 3.1562 | **−0.2125** |
| 1.20 | 2.7812 | 2.4187 | **−0.3625** |
| 1.50 | 2.0687 | 1.2187 | **−0.8500** |

All shifts **negative** (ODD suppresses ordering — the sign contract). **Attenuation pattern:**
the shift magnitude GROWS as Bc drops below the DS2023 ≈3.5 T crossover (−0.14 → −0.21 →
−0.36 → −0.85) — i.e. ODD is attenuated toward high Bx/low T (Bc > 3.5 T) and strongest near
the boundary's high-T foot, consistent with DS2023 (qualitative). The 1.5 K row sits only 9 mK
below the ODD Tc(0)=1.509 K, so Bc_on collapses well below the 2 T default window floor
(measured 1.2187 T at 1.5 K); the finder window floor was lowered from the default 2 T to
0.1 T to locate it (gate unchanged — sign only).

### Test status (TDD)

- **RED → GREEN:** `test_invz_odd_physics.m`. Fast structural gate (8³/dpRng15) AMENDED per the
  controller adjudication: directionality only (ODD lowers Tc; b_cond, c_sig < off; d>0; the
  c1≡a theorem at AbsTol 2e-3; naive/factorial finite) — magnitudes/closure REPORTED not gated.
- **Fast:** `test_zero_field_structure_fast` passes (8³/dpRng15: off 1.816, a 1.580, b_cond 1.676,
  c_sig 1.713, closure +0.0076, c1_lit≡a).
- **Slow (INVZ_SLOW=1, 45 s warm):** `test_zero_field_off_matches_published` ✓,
  `test_odd_headline` ✓ (reproducibility gate ACTIVE against pinned anchors, RelTol 1%),
  `test_boundary_shift` ✓ — **3 slow + 1 fast = 4/4 pass**.
- **Anchors pinned** (`invz_odd_anchors.m`, full precision): `odd_Tc_rich = 1.509370677196421`,
  `odd_d_at_Tc = 0.00048311966308299265`, `odd_Sc_rich = 0.38879543801229982`. Headline slow test
  re-run after pinning → reproducibility gate passes.

### Timings

Warm-cache `odd1_` (qVec_generator 12³/24³ dpRng 30): 13.9 s / 109.8 s one-off build (24³ the long
pole). Warm mode='full' {12,24} run = 10.6 s. Full slow file (off + headline + Bc-table) = 45 s.

---

## T2 — Tier 1b: retardation sensitivity (T2.1 retarded δJ(q, iωn) + T2.2 decision)

**Date:** 2026-07-17 · **Task:** main-body plan Task 7 (T2.1 + T2.2) · **Branch:**
`invz-1z-lihof4` · **Commit:** (this commit) · MATLAB R2025a.

### What was implemented

- **`invz_emt_scalar` generalized** (T2.1 plumbing): `Jnu_flat` may now be `[nJ, nw]`
  (per-frequency mode spectra; column n pairs with the caller's Matsubara grid, wn(1) = 0
  first). Any vector input takes the ORIGINAL code path with its expressions textually
  unchanged; a matrix with identical columns is **bitwise-equal** to the vector path
  (elementwise multiply/divide only, no accumulation-order differences — regression-gated
  by `test_emt_matrix_column_constant_bitwise`). The `opts.debug` closure block got the
  same `[nJ,nw]` branch; shape-mismatched matrices error (`invz:emtJnu`).
- **`invz_solve_point` `opts.odd_retarded`** (requires `opts.odd`, else `invz:oddArgs`):
  evaluates `invz_chiperp` at z = iωn on the solve's own `(T, Ecut)` Matsubara grid,
  forms the scalar surrogate `r_n = Xpn(1,1,n)/Xpn(1,1,1)` (asserted: r_1 = 1 to 1e-12,
  finite, 0 < r_n ≤ 1 + 1e-12, monotone non-increasing — error id `invz:oddRn`), and
  builds the per-frequency spectra by **first-order eigenvalue perturbation** around the
  static-ODD modes: `Jnu_n(q,ν) = Jnu_static(q,ν) + (r_n − 1)·u_ν'·δJ(q)·u_ν`, with the
  eigvec weights captured during the static rebuild from a SEPARATE `eig(M,'vector')`
  call so the values-only static spectrum stays bitwise-identical. Flattening is the
  static path's own `(:)` column-major order (all q at ν=1, then ν=2, …) for every
  column. `pt.odd` gains `r_n`, `ab_ratio_max`, `bb_aa_dev_max`, `retarded_exact`.
- **Key linearity fact (documented in the solver header):** E1 and E4 are linear in
  χ⊥, so under the scalar surrogate `Xp_n = r_n·Xp_0` the retarded coupling matrix is
  `δJ_n = r_n·δJ_0` **exactly** (and `d_n = r_n·d`) — the surrogate is exact in the
  MATRIX; only the eigenvalue map is approximated to first order. The exact cross-check
  path `opts.odd_retarded_exact` therefore does the full per-(q,n) eig of
  `Vcc + r_n·δJ` (test/diagnostic only, wins over `odd_retarded` when both set).
- **J0eff keeps the STATIC −d**: criticality (`pt.crit`) sits in the elastic n = 0
  sector where r_1 = 1 exactly and χ⊥(0) is exact; the n-dependent d_n = r_n·d never
  enters the uniform q = 0 bookkeeping (this is also why `invz_critical_T0field`'s
  closed form — and hence Tc(0) — is intrinsically static; plan §4 mitigation).
- **`opts.odd_rn_override`** (test hook): scalar or [nwn,1] replacing the computed r_n;
  skips the χ⊥(iωn) solve. With override = 1 the perturbation term is exactly 0·w and
  the EMT receives constant columns → bitwise-static (gated).

### r_n profile (χ⊥(iωn)/χ⊥(0), full 136-state electronuclear χ0)

At the plan §4 benchmark point **(1.53 K, 0 T)**, Ecut 40 meV (nwn = 50):

| n | ωn (meV) | r_n |
|---|---|---|
| 0 | 0 | 1.000000 |
| 1 | 0.8284 | **0.679424** |
| 2 | 1.6568 | 0.425389 |
| 3 | 2.4852 | 0.303567 |
| 4 | 3.3136 | 0.233259 |

r(ω₁ = 0.828 meV) = **0.679** vs the plan §4 two-level estimate ε₁²/(ε₁²+ω₁²) = 0.558:
the measured decay is SLOWER because the full CF + hyperfine χ⊥ carries Van-Vleck weight
at higher-lying poles (ε ≫ ε₁), which retard more weakly than the single-gap estimate.
Monotone to r_49 = 0.0062; χ_aa(0) = 17.64 meV⁻¹ (T1.2 anchor band). At the PM test
fixture (1.8 K, 0.5 T): r = [1, 0.6114, 0.3721, 0.2622, …].

**Scalar-surrogate smallness logs** (justify treating χ⊥(iωn) as r_n·χ⊥(0)):
max_n |χ_ab|/χ_aa = 9.5e-16 and max_n |χ_bb − χ_aa|/χ_aa = 2.4e-15 at (1.53 K, 0 T)
(machine zero — C4 unbroken at B = 0); 2.17e-3 and 2.79e-2 at (1.8 K, 0.5 T). The
residual matrix error from the aa/bb split at finite B is folded into the
surrogate-vs-exact gate below.

### Surrogate vs exact (per-(q,n) eig) — (1.8 K, 0.5 T), 6³/dpRng 10

Σ0: static 0.253803 → retarded 0.253803 (Δ = **3.9e-9**, retarded LOWER — the expected
direction: the static form overweights the n ≠ 0 channels). Surrogate-vs-exact
|ΔΣ0| = **1.76e-11** (gate 1e-3 — first-order perturbation is more than adequate;
per-n eig stays a test-only path). The retardation itself is NOT small in the channels:
K(n=1) shifts −6.6%, K(n_max) −9.0%, max_n |ΔΣ(n)| = 1.1e-5 — it is the elastic-sector
dominance of the λ-sums that makes Σ0 nearly immune. Sum-rule residual unchanged at the
1e-6 level (T2.1 acceptance; gated).

### T2.2 decision measurement (16³ ndgrid, dpRng 30, warm `odd1_` cache)

**Fixture finding (Tc leg amended — measured, pre-existing):** with ODD on, the ordered
side of the boundary NEVER converges to a paramagnetic fixed point — probed at
Bx = 0.1/0.3/0.5/2/2.5/3 T across each boundary: converged points exist only ≥ ~50–80 mK
ABOVE the crossing, none below (at 0.1 T even the PM side within the adaptive window
fails: the near-degenerate-doublet small-B caveat in `invz_critical_T`'s own header).
`invz_critical_T`'s converged-sign-change bracket is therefore structurally unreachable
at ANY field with ODD on — `invz_critical` survives only through its `para_edge`
fallback, which the T-cut finder lacks. **Flagged for V4.1** (ODD phase-diagram T-cuts
will hit this; fixing `invz_critical_T` is outside this task's file scope). The brief's
Tc(0.1 T) leg was replaced by a deterministic PM-side estimator at the healthy 2 T
proxy: crit(T) sampled on the house 1/30 K grid (converged votes only, identical grid
for both flags), two lowest converged PM points extrapolated to crit = 0 — the
static-vs-retarded DIFFERENCE is then a pure crit-shift readout with the estimator's
own PM-side bias cancelling.

| quantity | static | retarded | shift |
|---|---|---|---|
| Tc*(2 T) PM-extrapolated (K) | 1.259129 | 1.259151 | **+0.022 mK** |
| Bc(1.2 K) (T) | 2.4187 | 2.4187 | **0** (< 0.02 T tol) |
| Tc(0) closed form (K) | 1.509 | ≡ static | **0 by construction** (elastic sector) |

**Verdict (10 mK rule): |ΔTc| = 0.022 mK ≪ 10 mK → the STATIC form is frozen as the
default.** The retarded path stays available behind `opts.odd_retarded` (exact
cross-check behind `opts.odd_retarded_exact`). T1.5's headline numbers stand — no
re-measurement needed. README §1.9 documents the decision.

### Test status (TDD)

- **RED → GREEN:** `test_invz_odd_retarded.m`. Pre-implementation: EMT matrix test
  FAILED (bitwise mismatch through the old `Jnu_flat(:)` flatten), physics test ERRORED
  (`pt.odd.r_n` absent); the r_n-unity test passes trivially pre-implementation (unknown
  opts ignored) — its value is as the post-implementation bitwise regression gate, and
  it stays green with the branch live.
- **Fast (3 tests, 0.7 s):** EMT constant-column bitwise ✓; r_n ≡ 1 override bitwise-
  static ✓ (matrix route, no hidden vector short-circuit needed); physics + surrogate
  gate ✓. Full fast suite **131 passed / 0 failed / 17 incomplete** (baseline 128/0/16 +
  3 new + 1 new slow-gated).
- **Slow (INVZ_SLOW=1):** decision measurement ✓ (4/4 in the file).

### Timings

Fast trio 0.7 s. Slow decision leg 10 s warm (2 × 10-point Tc scans + 2 × Bc finders;
the retarded solve adds one 136-dim diagonalization + a ~50-point χ0(iωn) sweep + one
eigvec pass per point — ~25% overhead at 16³). Probe sweeps (0.1–3 T convergence
census): ~7 min total, one-off.

---

## T3.1 — Tier-2 internal-field covariance C (E3): `invz_odd_fieldvar.m` (2026-07-17)

Opens Phase T3 (Dollberg's variable-moments mechanism). `[C, info] =
invz_odd_fieldvar(ion, pt, S, T, opts)` evaluates E3 at a CONVERGED 1/z ODD point
solve: the 2×2 covariance of the internal transverse field seen by the doublet,

    C_ab = (1/(4·Nq)) Σ_q tr[ V_ac(q) · Scc(q) · V_cb(q) ]        (meV²)

with `V_ac = Vca(:,:,iq)'`, `Scc(q) = Σ_ν u_ν S_ν(q) u_ν'`, and the equal-time
per-mode structure factor from the converged EMT lattice propagator (framework
eq 7, identical form to `invz_emt_scalar`): `Gq_ν = G./(1 + (Jν_odd − K).*G)`,
`S_ν(q) = −(1/β)·Σ_n wts_n·Gq_ν(iω_n)`. **Sign convention** mirrors the solver's
sum rule `sum(wts.*pt.G)/β = −si.JzJz_fluct` (χ = −G): S_ν ≥ 0 for physical
modes at a converged PM point, making every per-(q,ν) term PSD. Modes/eigvecs
from `eig(Vcc + dJ)` per q (Hermitized, sorted ascending with U permuted in
lockstep — the T2.1 `w_odd` pattern); dJ is rebuilt from `pt.odd.Xp` via
`invz_odd_deltaJ`, so C sits at the SAME converged state as the solve
(deterministic rebuild). Guards: `invz:fieldvarArgs` (pt lacks ODD fields / not
converged / malformed S), `invz:fieldvarGrid` (`numel(pt.G)` vs the (T, Ecut)
Matsubara grid — pt must be solved with the SAME Ecut, default 40).

**Basis choice (plan T3.1 "document the choice"):** contraction in the MODE
basis — with y_α(ν) = Vc_α' u_ν, cyclic invariance gives per q
tr[V_ac·Scc·V_cb] = Σ_ν S_ν · (y_b(ν)' y_a(ν)) — note the CONJUGATE ORDER: the
(a,b) element carries y_b'y_a, so C_ba(q) = conj(C_ab(q)) and each per-q term
is exactly Hermitian (PSD when S_ν ≥ 0). Costs two extra 4×4 products per q on
top of the eig that is recomputed anyway; Scc is never formed in production.
The explicit-Scc sublattice-basis route is the brute-force cross-check inside
`test_fieldvar_structure` (on-axis q AND generic q = [1/3 1/6 1/6], AbsTol
1e-14·scale — both bases agree to machine precision). The assembled C is
asserted real symmetric to 1e-10 relative and defensively symmetrized
(`invz:fieldvarHerm`); ±q imaginary cancellation measured resid ~1e-21 meV².

**Measured symmetry finding:** on the 6³ mesh the per-q E3 integrand is REAL at
EVERY q (max rel imag ~1.5e-16; inversion symmetry) — the ±q cancellation is
trivially satisfied per point, while the per-q off-diagonal (rel ~0.25 of the
diagonal) cancels to ~0 only in the q-sum. Consequence: a conjugate-order slip
would be numerically invisible on this lattice; the order is fixed by the
header algebra + the two-basis cross-check.

### Headline numbers (6³ ndgrid − Γ, dpRng 10 warm cache, Ecut 40 → nwn = 43; report, not tuned)

| fixture (T, Bx) | C_aa (meV²) | C_bb (meV²) | C_ab (meV²) | heq (T) | tail_share | S_ν min |
|---|---|---|---|---|---|---|
| (1.8 K, 0.50 T) | 4.560e-4 | 4.565e-4 | −3.1e-8 | **0.2953** | 7.7e-4 | 12.0 |
| (1.8 K, 0.05 T) | 4.574e-4 | 4.574e-4 | −3.4e-10 | **0.2956** | 7.7e-4 | 11.9 |

- heq ≈ **0.30 T** vs Dollberg's distribution width h ≈ 0.4 T (their Fig. 4) —
  same order, qualitative comparator only (log-only, never asserted). C in the
  expected O(1e-4–1e-3) meV² band (d ≈ 0.47 μeV, S_cc ~ O(12–30)).
- C4: C_aa = C_bb to 9.8e-6 relative at Bx = 0.05 T (1.0e-3 at 0.5 T); the
  off-diagonal is 4 orders below the diagonal and shrinks with Bx → 0.
- Static approx (`opts.static_approx`, gated, NEVER silent — `info.
  static_approx`): S_ν ≈ kBT·χ_ν(q,0) keeps only the n = 0 Matsubara term →
  static C is LOWER, rel diff ‖Cs−Cf‖/‖Cf‖ = **0.062** at 0.5 T (0.057 at
  0.05 T), heq 0.286 vs 0.295 T. The n ≥ 1 quantum terms carry ~6% of the
  covariance at 1.8 K — the full sum stays the default.
- Truncation tail (Ecut 40): last-frequency contribution 7.7e-4 of the total
  (ω_n^−2 decay; gate < 0.05 passes with 60× margin).

### Test status (TDD)

- **RED → GREEN:** `test_invz_odd_fieldvar.m` written first — all 3 tests
  errored (`invz_odd_fieldvar` undefined), then passed post-implementation.
  One amendment during the loop: the structure test's added "generic q has
  complex off-diagonal" sanity line was measured FALSE (integrand real at every
  q, see the symmetry finding) and was replaced by the two-point two-basis
  cross-check.
- **Fast (3 tests, 0.45 s):** structure (real symmetric, PSD, nonzero, tail
  < 0.05, two-basis cross-check at 2 q) ✓; C_aa = C_bb at Bx → 0 (RelTol 5e-2,
  measured 9.8e-6) ✓; static-approx gated + logged ✓. Full fast suite
  **134 passed / 0 failed / 17 incomplete** (baseline 131/0/17 + 3 new).

### Timings

`invz_odd_fieldvar` itself: **8–20 ms** at 6³/nwn 43 (one 4×4 eig + two 4×4
products + a [43×4] propagator per q — negligible next to the 0.07–0.2 s warm
point solve). Fast-suite addition ~0.5 s.

---

## T3.2 — Gauss–Hermite-dressed doublet: `invz_twolevel_avg.m` (2026-07-17)

`[tla, avg] = invz_twolevel_avg(ion, T, Bx, C, opts)` averages the electronic
two-level response over the quenched Gaussian transverse-field distribution
N(0, C) (C from T3.1, meV²) and fits effective parameters so the Σ machinery
keeps its two-level algebra.

> **FLAG (plan T3.2, carried verbatim-in-spirit in the function header):** this
> quenched-Gaussian, STATIC dressing of the doublet is the **least rigorous
> step of the whole plan** — thermal/quantum field fluctuations are frozen into
> a classical distribution and the averaged response is re-compressed into a
> single pole. Treat downstream numbers accordingly. (README wiring: Task 10.)

**Mechanics.** 2-D Gauss–Hermite (Golub–Welsch on the Hermite Jacobi matrix,
local subfunction, no toolbox; ngh 7 default; weights sum to 1, asserted at
1e-12), C diagonalized first; node fields enter as pure field shifts
`B_node = invz_field_vec(Bx) + [ha hb 0]/(gL·µB)` (NO single-ion code changes).
Default `opts.avg = 'response'`: average g(iωn) node-by-node on `opts.wn`, fit
`Delta_eff = sqrt(tail/gbar0)`, `n01_eff = gbar0·Delta_eff/2` — the g(0) and
ωn⁻² tail conditions match EXACTLY (tail coefficient Σw·2·n01·Δ analytic per
node, not fitted from the grid); `M2_eff = Σw·M2·g0/gbar0` (response-weighted:
M² sits outside g, and this choice makes `M2_eff·ḡ(0)` reproduce the averaged
static χ exactly). `opts.avg = 'params'` (plain (Δ, M2, n01) averages) kept for
the plan's one-shot comparison. C = 0 short-circuits to the literal
`invz_twolevel` struct (bitwise, tested). `avg.G0`/`avg.chi0cc0` (node-averaged
FULL electronuclear cc propagator, hyp = true) are OPT-IN via `opts.G0 = true`
(adjudication 2026-07-17: `opts.wn` alone no longer triggers the 49
electronuclear diagonalizations — that side effect was burning fast-suite
budget; Task 10 passes `G0 = true`).

### Monotonicity ray (1.6 K, 2 T, C = s·eye(2), ngh 7, Ecut 40) — report, not tuned

bare: Δ = 0.151653245 meV, M2 = 27.7913213, χ0_2l = M2·g0 = 183.435 meV⁻¹.

| s (meV²) | Δ_eff (resp) | Δ_eff (params) | M2_eff (resp) | χ0_2l = M2_eff·2n01_eff/Δ_eff | fit_resid |
|---|---|---|---|---|---|
| 1e-5   | 0.151766437 | 0.151713319 | 27.7918226 | 183.415 | 7.38e-6 |
| 8e-5   | 0.152554388 | 0.152134168 | 27.7952784 | 183.279 | 5.90e-5 |
| 2e-4   | 0.153887502 | 0.152856893 | 27.8009871 | 183.046 | 1.47e-4 |
| 6.4e-4 | 0.158593356 | 0.155518644 | 27.8196293 | 182.206 | 4.71e-4 |

- **Level repulsion ✓:** Δ_eff monotone UP, ≥ bare (both modes).
- **Variable-moments suppression ✓ (the gated observable):** the averaged
  static two-level susceptibility χ0_2l is monotone DOWN and below bare —
  the criticality-relevant direction the Tier-2 mechanism needs.

### Finding — plan expectation REFUTED: M2_eff RISES along the C-ray at 2 T

The plan's acceptance "M²_eff ≤ M²₀ monotonically in ‖C‖" encodes the near-B=0
intuition and is measurably FALSE at the 2 T fixture, in BOTH averaging modes.
Root cause (implementation-independent): bare M2(Bx) at 1.6 K is decreasing but
**CONVEX** at 2 T — M2''(2 T) = **+0.751 /T²** (3-pt/5-pt agree) — so a
symmetric quenched spread RAISES the average by ½·M2''·σ_B²: predicted +7.17e-4
at s = 1e-5 (a-axis), measured **+7.1658e-4** by a GH-free dense-trapezoid
Gaussian average (also validates the quadrature). Directional split: a-axis
disorder raises M2 (+1.0e-3 response-weighted, the extra ~+3e-4 being the
cov(M2, g0) > 0 term), b-axis alone lowers it (−5.1e-4, the |B|-lengthening
route) — the a-convexity wins in the isotropic sum. The direction is
plan-compliant at Bx = 1 T (dM2_a and dM2_b both < 0) and reversed at
2/3/4/4.5 T (+1.62e-3/+2.14e-3/+1.16e-3/+7.9e-4 at s = 8e-5, params mode).
Adjudicated (Task-9 round, plan amended b691b35): gate χ0_2l, REPORT M2_eff.

### Fit residual and sum rule

`fit_resid` = max relative mismatch of the three fit conditions reproduced by
the fitted form. g(0) and tail are exact by construction (~1e-16); the honest
leg is the truncated Matsubara sum rule: per node (1/β)Σ wts·g =
n01·coth(βΔ/2) = tanh·coth = **1 exactly** (analytic identity, verified in the
code comment), truncated-grid value ≈ 1 − 2·n01·Δ/(π·Ecut) (asserted within
2e-2 of 1 per node, `invz:tlavgSumrule`); the fitted pair (Δ_eff, n01_eff) no
longer satisfies n01 = tanh(βΔ/2), so the residual is SECOND order in the node
spread, i.e. linear in ‖C‖: measured band **7.4e-6 → 4.7e-4** on the ray (ratios
track s exactly). Gate adjudicated 1e-6 → **1e-3**. Sum-rule survival:
`sumrule_avg = 0.998793928` at C = 2e-4·eye(2) (truncation-level, not drift;
bare-node 0.998824 — the averaging shift is 3.0e-5, same second-order origin).

### Degenerate nodes (T3.4 groundwork)

Nodes with Δ_node < 1e-4 meV (invz:degenerateDoublet — e.g. the exact-origin
node at Bx = 0) are evaluated in their h → 0⁺ limit: g0 = β, tail = 0,
n01 = 0, m = 0, g(iωn) = β·(ωn = 0) (elastic Curie weight; per-node sum rule
= 1 exactly); other errors (incl. invz:orderedPhase) rethrown. The node M² is
the **basis-invariant** `|Mz11|² + |Mz12|² = (Mz²)(1,1)` — equal to m² in the
Jz-diagonal basis, i.e. the doublet's Curie weight (controller-approved: at
exact degeneracy `eig` returns an arbitrary doublet basis, so |Mz12|² alone
would be randomized in [0, m²]). Measured at (1.55 K, B = 0, C = 2e-4·eye(2)):
`n_degenerate = 1` of 49 (origin node, raw Δ = 1.4e-14 meV), off-origin nodes
lift the degeneracy to ≤ 0.0488 meV, dressed **Δ_eff = 5.224e-3 meV** (finite,
> 0), M2_eff = 30.28 (the m² band), fit_resid = 1.5e-8.

### Response vs params (one-shot comparison, plan T3.2) and GH convergence

- C = 2e-4·eye(2), (1.6 K, 2 T): Δ_eff **0.153888 (response) vs 0.152857
  (params)** meV — the inequivalence is at the 0.7% level of the shift.
- GH 5/7/9 at the same C: Δ_eff = 0.153887502443 / 0.153887502321 /
  0.153887502321 — |Δ(5→7)| = 1.2e-10, |Δ(7→9)| < 1e-12: machine-converged by
  ngh 7 (smooth integrand); ngh 7 default confirmed.

### Test status (TDD + adjudication round)

- **RED:** brief's tests transcribed verbatim first → 6/6 errored (undefined).
- **Round 1:** 5/6 green; `test_gh_weights_and_monotonicity` failed exactly on
  the two M2 legs + the fit_resid 1e-6 gate → escalated BLOCKED with the
  convexity attribution (task-9-report.md). **Adjudication** (plan amendment
  b691b35): M2 legs → χ0_2l suppression gate + M2 direction reported;
  fit_resid gate 1e-3; `opts.G0` opt-in; degenerate-M² deviation approved.
- **GREEN:** 7/7 (6 amended-brief tests + `test_G0_optin` covering Task 10's
  consumer path: G0 shape/finiteness, chi0cc0 = −G0(1) > 0, G0-without-wn
  errors), **0.95 s** (was 23.3 s before the G0 opt-in gating). Full fast
  suite **141 passed / 0 failed / 17 incomplete** (baseline 134/0/17 + 7 new).

### Timings

Two-level-only averaged call (ngh 7, 49 nodes, no G0): ~0.05–0.1 s. With
`G0 = true` (49 hyp-true 136-dim MF solves + chi0z): ~1.5–2.5 s per call
(2.42 s at the B = 0 fixture pre-gating). Fast-suite addition: **+0.95 s**.

---

## T3.3 — Tier-2 outer self-consistency in `invz_solve_point` + T3.4 small-Bx characterization (2026-07-17)

`opts.odd_tier2 = true` (requires `opts.odd`; `invz:oddArgs` otherwise) closes
Dollberg's variable-moments loop. After the Tier-1 solve converges: seed
`C = invz_odd_fieldvar(ion, pt, blocks, T)` from the CONVERGED Tier-1 point →
`[tla, avg] = invz_twolevel_avg(..., C, 'G0', true)` (GH-dressed doublet +
disorder-averaged FULL electronuclear cc propagator) → re-run the inner
EMT↔Σ loop with `tl = tla`, `g = invz_g(tla)`, `G0 = avg.G0` (previous Σ as
warm seed — the plan's warm-start economy) → re-evaluate C, mix with damping
0.5, repeat until `max|dC| < opts.tol_tier2` (1e-3 relative; cap
`opts.max_tier2` = 8). Non-convergence (inner OR outer₂) is masked EXACTLY
like the EMT loop: `pt.converged = false` — the ordered-phase signal the
critical finders rely on. New fields (only when the flag is on): `pt.C`
(final mixed covariance; NaN(2) when Tier 1 never converged),
`pt.tier2_iters`, `pt.tla` (final dressed params; `pt.tl == pt.tla` on a
completed Tier-2 solve), `pt.tier2_resid`.

**What Tier 2 swaps and what it does NOT (design, per plan):** the dressed
doublet enters Σ (tl, g) and the local propagator (G0), so `pt.chi0cc0` and
`pt.crit = 1 + Σ0 − J0eff·χ̄0cc(0)` carry the DRESSED χ̄0cc — that is the
variable-moments suppression channel. The LATTICE geometric side (the δJ
rebuild from Vcc + deltaJ and the −d on J0eff) stays at the Tier-1 static
χ⊥ (`pt.odd.Xp`) — χ⊥ is Van-Vleck-dominated and is NOT re-dressed (the
χ⊥-held-fixed design decision, README §1.9). `pt.si` stays the bare
single-ion state and `pt.sumrule_rel` keeps the bare `si.JzJz_fluct`
reference (documented O(C) diagnostic offset: 0.0564 → 0.0582 at the gate).

**Loop placement (implementation note):** the flag-off inner Σ loop is
textually untouched (bit-for-bit non-negotiable); the Tier-2 path re-enters
the loop through `local_sigma_loop`, a verbatim mirror local function used
ONLY by the outer₂ iteration and cross-referenced in both places.

**Retarded × Tier-2 decision:** `opts.odd_tier2` combined with
`opts.odd_retarded(_exact)` errors `invz:oddArgs` ("not yet validated") — a
real coupling subtlety, not caution theatre: `invz_odd_fieldvar` assembles
the equal-time Scc (E3) from the STATIC mode spectrum `eig(Vcc + dJ)` with
this solve's G/K, while a retarded solve's lattice propagator carries
per-frequency modes `Jnu_n(q,ν)`; E3 has no validated retarded counterpart.
The static default (T2.2, |dTc| = 0.022 mK) makes the combination unused in
practice.

### Headline numbers (report, not tuned)

**Convergence gate (1.80 K, 0.05 T; 6³/dpRng 10 warm) — the plan's (1.55 K,
0.05 T) point is REPORTED below (brief amendment: the gate must sit at a
guaranteed-PM point):**

| quantity | Tier 1 | Tier 1+2 |
|---|---|---|
| converged | 1 | 1 (tier2_iters = **2**, resid 7.6e-4) |
| crit | 0.091270 | **0.094680** (increases: suppression ✓) |
| Δ (meV) | 1.3033e-4 | **1.1789e-2** (level repulsion ✓, ×90: the internal field, not the 50 mT applied field, sets the dressed splitting) |
| Σ0 | 0.259103 | 0.255571 (falls — but the dressed χ0cc0 drop wins) |
| χ0cc0 (meV⁻¹) | 195.7178 | 194.5543 (−0.59%; J0eff·Δχ ≈ +0.0070 on crit vs ΔΣ0 = −0.0035) |
| point-solve time | 0.14 s | 5.2 s (dominated by the 49-node electronuclear Ḡ0 average, ~2.4 s/outer₂ iter) |

**C evolution across outer₂ iterations at the gate (max_tier2 sweep):** seed
(Tier-1 covariance) C_aa = 4.574e-4 → k=1: C_new rel step 1.51e-3, mixed
4.5708e-4 → k=2: rel step **7.56e-4 < 1e-3 converged**, final C_aa =
4.5690e-4 meV² (C_ab ~ −3e-10). Monotone contraction, no oscillation — the
χcc feedback is stabilizing exactly as the plan argued (dressing suppresses
Scc, which shrinks C by ~0.2% total); mixing 0.5 never stressed. Expected
2–4 iterations: measured 2 (1.80 K) and 1 (1.55 K).

**Plan point (1.55 K, 0.05 T), REPORTED:** Tier 1 converged PM crit =
+0.005426 (1.55 K sits above Tc_ODD(0) = 1.509 K — already paramagnetic at
Tier 1); Tier 1+2 converged PM, crit = **+0.008431**, tier2_iters = 1, C_aa
= 4.4449e-4 meV², Δ_eff = 0.011488 meV. The point lands on the PM side of
BOTH the Tier-1 and the Tier-1+2 boundaries.

**IR safety (1.2 K, 16³/dpRng 30; slow gate):** Bc(1.2 K) = 2.4187 T
(para-edge estimate — the `invz:orderedSideNoConverge` warning is the
EXPECTED pre-existing T2.2 behaviour: the ordered side never converges with
ODD on; the value matches T2.2's 2.4187 T exactly). Approaching the boundary
from above:

| B | crit | ‖C‖ (meV²) |
|---|---|---|
| Bc + 0.5 T | 0.0823 | 4.25e-4 |
| Bc + 0.2 T | 0.0425 | 4.25e-4 |
| Bc + 0.05 T | 0.0264 | 4.24e-4 |

crit falls 3× while ‖C‖ is FLAT to 0.2% (gate ‖C‖(+0.05)/‖C‖(+0.5) < 20
passes with margin ~0.998): the ODD blocks vanish at q → 0, so the softening
uniform mode never feeds C — the mechanism's built-in infrared safety,
measured even stronger than the "grows but saturates" expectation.

**Combined ΔTc at the 0.5 T proxy (16³/dpRng 30; PM-side crit(T)
extrapolation estimator — brief amendment: `invz_critical_T` cannot bracket
with ODD on, T2.2 finding — common grid 1.30:1/30:1.90 K; top raised from
1.80 during TDD: the off-config PM-side slowing band at 0.5 T is ~0.06 K
wide (1.7333/1.7667 non-converged, first converged PM point 1.8000, crit
+0.0253), so a 1.80 top left a single converged point):**

| config | extrapolated from T (K) | crit at those points | Tc* (K) | distance flag |
|---|---|---|---|---|
| off | [1.8000, 1.8333] | [0.02526, 0.03375] | **1.7008** | 3.0 grid steps |
| +Tier 1 | [1.6000, 1.6333] | [0.03408, 0.04408] | **1.4865** | 3.4 grid steps |
| +Tier 1+2 | [1.6000, 1.6333] | [0.03701, 0.04731] | **1.4802** | 3.6 grid steps |

- Combined suppression 0.2206 K; split **Tier 1 : Tier 2 = 97.2% : 2.8%**
  (0.2143 K : 0.0063 K). Tier 2 is a small, correctly-signed addition
  (t_t12 ≤ t_t1 ≤ t_off gates pass); consistent with Tier 1's zero-field
  0.234 K (T1.5) — the 0.5 T proxy reads slightly lower, as expected off
  axis.
- Estimator adaptation (documented in the test helper): T2.2's cross-run
  same-grid-points guard compared two nearly identical configs; here the
  three configs' Tc differ by ~0.2 K, so each extrapolates near ITS own
  boundary with per-config internal guards (≥2 converged PM points, crit >
  0, PM-side monotonicity) and the points used are REPORTED. The two ODD
  configs happen to use IDENTICAL points [1.6000, 1.6333], so the Tier-2 leg
  IS a pure dressed-crit readout; the off leg sits on its own points. All
  three distance flags fired at ~3 grid steps — the estimator's common
  PM-side extrapolation bias (the crit = 0 crossing hides inside the
  non-convergent band, T2.2 fixture finding); it cancels in the differences.

### T3.4 — small-Bx / B = 0 characterization (report)

- **Exact B = [0 0 0]:** both Tier-1 and Tier-1+2 solves throw
  `invz:degenerateDoublet` ("Doublet splitting 1.07e-14 meV too small: Bx=0
  limit needs the closed-form Sigma_c") from the BARE `invz_twolevel`
  scaffold in the Tier-1 path — the error fires before any Tier-2 machinery
  runs. **The guard is deliberate, not accidental:** the Σ machinery needs a
  finite bare doublet to seed the solve; the sanctioned zero-field route
  remains the closed form (`invz_sigma_crit` / `invz_critical_T0field` /
  `invz_odd_zero_field`), and the 0.05 T proxy is the practical small-field
  route (plan T3.3).
- **The Tier-2 dressing machinery itself is B = 0-capable** (Task 9's
  degenerate-node limit): `invz_twolevel_avg(1.55 K, B = 0, C(1.55 K))`
  returns Δ_eff = 0.011313 meV with n_degenerate = 1/49 (origin node only;
  the internal field lifts the other 48), fit_resid = 3.1e-7, n01_eff =
  0.0423. So a hypothetical B = 0 Tier-2 entry point would only need a
  dressed (not bare) scaffold — noted as a Tier-3/V4 design option, not
  implemented (YAGNI: the closed form owns B = 0).
- **Quantitative vs directional verdict: the hyperfine caveat STANDS.** The
  dressed splitting at the 0.05 T proxy is Δ_eff ≈ 11.5–11.8 μeV — only
  ~3.4× the hyperfine scale A ≈ 3.36 μeV. The quenched-Gaussian dressing
  acts on the ELECTRONIC doublet (hyperfine enters only through the full-χ0
  node average in Ḡ0), and the single-pole recompression is the plan's
  flagged least rigorous step — at scales where Δ_eff ~ A the electronuclear
  ladder structure matters. Zero-field/small-Bx Tier-2 results are therefore
  DIRECTIONAL (sign and rough magnitude of the suppression), not
  quantitative.
- **Task-9 WATCH item (n01_eff > 1 at low T/large spread): NOT triggered.**
  Along the whole characterization: n01_eff = 0.0380 (1.80 K), 0.0430
  (1.55 K), 0.0423 (1.55 K, B = 0); fit_resid 2.0e-7–3.3e-7 (fit gate 1e-3);
  n_degenerate 0/49 at Bx = 0.05 T (the applied field lifts every node),
  1/49 at exact B = 0; node_Delta spread 1.3e-4–0.105 meV at 1.55 K.

### Test status (TDD)

- **RED:** `test_invz_odd_tier2.m` written first — flag-off test FAILED on
  the structural guards (unknown opts are ignored, so the isequaln leg
  passes trivially pre-implementation; the verifyError legs are what bite),
  convergence-gate test ERRORED (`pt.tier2_iters` undefined); slow pair
  filtered.
- **GREEN (fast, 2 tests, 7.8 s):** flag-off isequaln + three `invz:oddArgs`
  guards (requires-odd, ×retarded, ×retarded_exact) ✓; convergence gate +
  suppression/level-repulsion directions + the 1.55 K report point ✓. Full
  fast suite **143 passed / 0 failed / 19 incomplete** (baseline 141/0/17 +
  2 new fast + 2 new slow-gated).
- **GREEN (slow, INVZ_SLOW=1 file-scoped): 4/4 in 56 s** (warm odd1 caches;
  the 30–60 min budget was pessimistic — non-converged grid votes fail
  fast). One amendment during the slow TDD loop, documented in the test:
  combined-grid top 1.80 → 1.90 K (off-config slowing band, measured above);
  the IR fprintf widened to report crit alongside ‖C‖.

### Timings

Tier-2 point solve: 5.2 s at 6³ (vs 0.14 s Tier 1; ~2.4 s per outer₂
iteration is the 49-node electronuclear Ḡ0 average — mesh-independent), ~8 s
at 16³ near the boundary. Slow file 56 s total (IR leg: one `invz_critical`
+ 3 point solves + 3 fieldvar; combined leg: 19-point grid × 3 configs).
Fast-suite addition +7.9 s (one Tier-2 solve at the gate + the 1.55 K
report point dominate).

---

## V4 — Headline overlay, robustness error bars, consolidation, handoff (V4.1–V4.3)

**Date:** 2026-07-17 · **Task:** main-body plan Task 11 (V4.1–V4.3) · **Branch:**
`invz-1z-lihof4` · MATLAB R2025a. The LAST main-body task: driver overlay,
robustness sweeps, suite consolidation, handoff. All numbers **reported, not tuned.**

### V4.1 — Quick ODD overlay (`invz_run_phase_diagram`, `overlay_quick` block)

The driver gained an opt-in `overlay_quick` block (default false → the standard
sweep is byte-identical). It draws the headline overlay with **B-CUT points
ONLY** (`invz_critical`), at three temperatures below Tc_ODD(0) = 1.509 K, plus
the closed-form zero-field endpoints and R2007 Fig.1 experimental anchors, with
a second Σ(0) panel. **The T-cut wall binds (ODD-LOG T2.2):** `invz_critical_T`
CANNOT bracket with ODD on (no metastable PM window below the boundary; it lacks
the `para_edge` fallback), so the near-Tc0 T-cut region of the ODD curves is
LEFT BLANK — the missing para_edge analog for the T-cut finder is flagged as a
V4-scope item (not fixed here). B-cuts survive because `invz_critical`'s
`para_edge` fallback returns the paramagnet-edge estimate when the ordered side
won't reconverge.

**Figure:** `Data/Phase_ODD_overlay_quick.fig` (saved, 18.5 kB). All five overlay
elements present + the Σ(0) panel + closed-form endpoints:

| element | present | detail |
|---|---|---|
| 1/z baseline (ODD off) | ✔ | B-cut crossing + closed-form Tc(0) endpoint |
| 1/z + Tier-1 ODD | ✔ | B-cut para_edge + closed-form Tc(0) endpoint |
| 1/z + Tier-1+2 ODD | ✔ | finite-B only (Tier-2 Tc* unavailable at B=0 — legend/caption say so) |
| mean field (Σ=0) | ✔ | crit_MF = 1 − J0·χ0cc0; MF Tc(0) = 2.259 K (README "MF 2.26 K" ✓) |
| experimental anchors | ✔ | R2007 Fig.1 (Bitko 1996 / Babkevich 2016): Tc(0)=1.53 K, Hc(T→0)≈4.95 T |
| Σ(0) 2nd panel | ✔ | critical self-energy along the boundary, off vs ODD, B=0 endpoints Richardson |

**Overlay B-cut boundary (generator-16³, dpRng 30, para_edge for ODD):**

| T (K) | Bc off (T) | Bc +T1 (T) | Bc +T1+2 (T) | Bc MF (T) | Σ0 off / T1 / T1+2 |
|---|---|---|---|---|---|
| 0.80 | 3.3257 | 3.1037 | 3.1086 | 3.8596 | 0.1703 / 0.2037 / 0.2064 |
| 1.20 | 2.7319 | 2.3280 | 2.3081 | 3.5203 | 0.2254 / 0.2884 / 0.2903 |
| 1.40 | 2.2884 | 1.6074 | 1.5466 | 3.2672 | 0.2466 / 0.3300 / 0.3321 |

- **Sign contract holds everywhere:** Tier-1 ODD lowers Bc at every T; Σ0
  INCREASES with ODD (fluctuation weight gained). MF boundary sits well ABOVE
  the 1/z boundary (it needs χ0cc0 = 1/J0 < (1+Σc)/J0). Consistent with the
  §T1.5 ndgrid Bc table (0.80 K: −0.222 T here vs −0.2125 there; mesh + para_edge
  differences account for the ~few-mT offset).
- **Tier-2 is tiny and sign-noisy at the low-T foot:** +T1→+T1+2 shifts −0.5 mT
  (0.80 K, an inversion within the finder's 0.02 T tolerance), −19.9 mT (1.20 K),
  −60.8 mT (1.40 K) — i.e. Tier-2 grows toward the high-T foot and stays ≤ ~4–5%
  of the Tier-1 shift, matching the 97.2:2.8 split at the 0.5 T proxy (§T3.3). The
  0.80 K wobble is genuine smallness, not a bug (ODD is attenuated at low T /
  high field, §T1.5). Tier-2 B-cuts finished within the 15 min budget (no drop).
- **Closed-form endpoints (Richardson 12/24):** Tc(0) MF 2.259 / off 1.74347 /
  Tier-1 ODD 1.50937 K; Σc(0) off 0.29798 → ODD 0.38880.

The **production sweep** (dense T grid, full para-edge boundary, hours, one run
per config) is documented in the driver header and LEFT TO THE USER (repo
precedent).

### V4.2 — Robustness sweeps (each ONE logged point) + error bars

**(i) Grid convergence — Tc(0), Σc(0), d at 12³/16³/24³ + Richardson (mode 'full', dpRng 30):**

| n | Tc_off (K) | Σc_off | Tc_ODD (K) | Σc_ODD | d (μeV) | ΔTc (K) |
|---|---|---|---|---|---|---|
| 12³ | 1.79094 | 0.26388 | 1.55505 | 0.34489 | 0.4912 | 0.23589 |
| 16³ | 1.77945 | 0.27208 | 1.54424 | 0.35525 | 0.4871 | 0.23521 |
| 24³ | 1.76720 | 0.28093 | 1.53221 | 0.36684 | 0.4830 | 0.23499 |
| **Rich (12,24)** | **1.74347** | **0.29798** | **1.50937** | **0.38880** | **0.4831** | **0.23410** |

- **Monotone in n** (Tc(12) > Tc(16) > Tc(24) > Tc_rich, both off and ODD) — no
  erratic behavior, no escalation.
- **ODD is NOT slower-converging than the baseline Σc.** The discretization
  errors |Tc(n) − Tc_rich| are **off [0.0475, 0.0360, 0.0237]** vs
  **ODD [0.0457, 0.0349, 0.0228]** — nearly identical, ODD marginally *faster*.
  So δJ shifting weight to finite q does NOT degrade the O(1/n) grid convergence
  of Tc(0) (the Richardson (12,24) correction remains valid for the ODD route).
- ΔTc(0) is remarkably grid-stable (0.23589 → 0.23499 → 0.23410): the two per-grid
  Tc bars (~0.023–0.024 K) are strongly correlated and **cancel in the difference**.
- (d-value definition note: this table's d column is evaluated at each grid's
  own Tc — 12³ 0.4912, 24³ 0.4830 μeV — whereas §T1.5 quoted d at Tc_rich for
  both grids (12³ 0.4914, 24³-at-Tc_rich 0.4832 μeV); same physics, different
  evaluation temperature, hence the small digit drift between the two tables.)

**(ii) dpRng sensitivity of d (16³, dpRng 20/30/40, d at 1.53 K):**

| dpRng | d (meV) | d (μeV) |
|---|---|---|
| 20 | 4.871794e-04 | 0.487179 |
| 30 | 4.871684e-04 | 0.487168 |
| 40 | 4.871722e-04 | 0.487172 |

- **Spread max−min = 1.11e-08 meV = 0.0023 % of d(30)** → d is dpRng-flat.
  This **confirms the effective r⁻⁶ short-rangedness**: contracting two r⁻³ ODD
  kernels through χ⊥ gives an r⁻⁶ mediated coupling whose uniform reduction d is
  fully captured by dpRng 20 (the dpRng error bar on d is negligible). The
  smallest ODD-block shell itself is the P0.3(iv)/§T1.3 near-a*-axis shell
  (~0.178·Jcc0 at the 16³ ndgrid [1/16 0 0]; the c*-axis shell is machine-zero);
  the contracted observable d washes out that geometry entirely.

**(iii) Gauss–Hermite node count (copied from §T3.2, machine-converged):**
Δ_eff(ngh 5/7/9) = 0.153887502443 / 0.153887502321 / 0.153887502321 meV,
|Δ(5→7)| = 1.2e-10, |Δ(7→9)| < 1e-12 → **GH error bar NEGLIGIBLE** (ngh 7 default).

**FINAL HEADLINE TABLE (each number + error bar):**

| quantity | value | error bar | source |
|---|---|---|---|
| Tc(0), 1/z baseline (off) | 1.74347 K | ± 0.024 K (grid: \|Tc(24)−rich\|) | V4.2(i) |
| Tc(0), 1/z + Tier-1 ODD | **1.50937 K** | ± 0.023 K (grid) | V4.2(i) |
| **ΔTc(0)** | **0.2341 K** (13.4 %) | ± 0.001 K (grid; the Tc bars cancel) | V4.2(i) |
| gap closed vs exp (1.53 K) | ≈ **111 %** (21 mK overshoot: 1.509 vs 1.53) | — | §T1.5 |
| **ΔΣc(0)** | **+0.091** (0.298 → 0.389) | ± 0.005 (grid) | V4.2(i) |
| **d(Tc)** | **0.483 μeV** | ± 0.0023% (2.3e-5 fractional) (dpRng); ± 0.004 μeV (grid) | V4.2(ii) |
| GH quadrature (Δ_eff) | converged | ± 1.2e-10 (negligible) | §T3.2 |

**Bc(T) shifts (16³ ndgrid, §T1.5 canonical table):** all negative
(−0.1365 / −0.2125 / −0.3625 / −0.8500 T at 0.31 / 0.80 / 1.20 / 1.50 K), worst
−0.850 T at 1.50 K; **attenuation toward high field** — the shift magnitude GROWS
as Bc drops below the DS2023 ≈ 3.5 T crossover (ODD strongest at the high-T foot,
attenuated at high Bx / low T). The V4.1 overlay B-cuts corroborate the pattern.

### V4.3 — Consolidation

- **ODD fast-test additions: 22.35 s** across the 8 ODD test files (file-scoped,
  fast): chiperp 1.65 / odd_blocks 2.30 / odd_fieldvar 0.62 / odd_physics 6.42 /
  odd_retarded 0.83 / odd_solve 1.27 / odd_tier2 8.18 / twolevel_avg 1.07 s
  (34 Passed / 0 Failed / 7 Incomplete). **Under the 30 s budget → NO INVZ_SLOW
  demotion needed** (heaviest = odd_tier2 8.2 s, odd_physics 6.4 s; both stay
  fast).
- **FAST suite: 143 Passed / 0 Failed / 19 Incomplete** (53.1 s) — matches the
  frozen baseline exactly (P0.5: 109/0/12 at `5f4ff92`; +34 ODD passes, +7
  slow-gated, **0 failed**).
- **SLOW suite (`INVZ_SLOW=1`, everything): 162 Passed / 0 Failed / 0 Incomplete**
  (191.4 s, `SLOW_EXIT=0`) — green (P0.5 baseline: 121/0/0 at `360dfab`; +41 ODD
  passes as the slow gates un-skip).
- **Flag-off parity:** 0 failed in BOTH suites is the gate; the frozen benchmark
  digits (Σc(0) Richardson 0.2980, `info.Jcc0` 6.421e-3 meV, Tc(0) 1.74 K,
  `pt.Sigma0` 0.0932; P0.5 table) are asserted internally by the pre-existing
  tests and all stayed green — the ODD extension is bit-for-bit additive vs the
  `5f4ff92` frozen baseline.

### Handoff, README, plan state

- **Handoff:** `docs/SESSION-2026-07-16-invz-odd-mainbody.md` (module map of every
  new function/flag/cache key; adjudication history; T2.2 decision; T3.4 verdict;
  open Tier-3 items + deferred Appendix-A/`invz_tensor` pointer; production runs
  left to the user).
- **README §7:** one row (Tc(0) 1/z+ODD reproducibility anchor) + a pointer to
  §1.9 and the pinned anchors / slow reproducibility test.
- **Anchors:** `invz_odd_anchors.m` already pins `odd_Tc_rich` / `odd_d_at_Tc` /
  `odd_Sc_rich` (used by the headline slow test); no new pins needed (YAGNI — the
  V4.2 error bars are reported, not asserted).
- **Plan HTML:** V4.1/V4.2/V4.3 flipped to `[x]` + manifest `done`; **all 20
  main-body tasks now show `[x]` (grep: 20 `[x]`, 0 `[ ]`, 0 `todo`).**

### Timings

V4.2 grid sweep ~1.5 min (warm 12/24, 16³-generator build 40 s); dpRng sweep
~1.2 min (dpRng-20 build 10 s, dpRng-40 build 61 s, dpRng-30 warm); overlay
~4 min (Tier-2 B-cuts the pole, within budget). Fast suite 53 s; slow suite 191 s.

---
