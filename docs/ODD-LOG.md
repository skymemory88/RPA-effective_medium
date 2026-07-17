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
below the ODD Tc(0)=1.509 K, so Bc_on collapses toward B≈0; the finder window floor was lowered
from the default 2 T to 0.1 T to locate it (gate unchanged — sign only).

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
