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
