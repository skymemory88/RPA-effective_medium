# Session Handoff — 1/z Effective-Medium Method for LiHoF4 (`invz/`)

Date: 2026-07-06 → 2026-07-07. Branch: `invz-1z-lihof4` (pushed to
`https://github.com/skymemory88/RPA-effective_medium`, PR can be opened at
`.../pull/new/invz-1z-lihof4`). This note is the single re-entry point for
debugging or further development.

## 1. What was built

A self-contained MATLAB module `invz/` implementing Jensen's high-density 1/z
expansion (J. Jensen, PRB **49**, 11833 (1994)) for LiHoF4, following the
LiHoF4-specific application of Rønnow et al., PRB **75**, 054426 (2007)
("R2007"). Scope: **paramagnetic phase only** (thermal + quantum paramagnet,
transverse field B‖a). It computes the 1/z-renormalized susceptibility
χ(q,ω) one order beyond MF-RPA.

Reference documents (all in this repo):
- `jensen_1z_framework.html` — full derivation of the theory with equation
  numbers ("HTML eq N" ↔ "J 2.x" of Jensen's paper). Verified line-by-line
  against the published PDF this session (two errata in the *published paper*
  are documented inside: J 2.24 prints m²n₀₁, should be m·n₀₁).
- `docs/superpowers/plans/2026-07-06-invz-lihof4-susceptibility.md` — the
  implementation plan (11 tasks). Note two stale points: qVec_generator's
  parameter is `'grid'` not `'size'`, and `invz_chi_realaxis` returns
  `chi_cc_q [nJsel×nw]` keyed by coupling values, not q-points.
- Three paper PDFs in the repo root (untracked): Jensen 1994 (HoF3),
  R2007 (LiHoF4 1/z), Kovacevic et al. PRB **94**, 214433 (2016)
  (electronuclear model + GHz χ″ data — the observable target).

## 2. Module map (call chain top-down)

```
invz_ion.m                LiHoF4 parameters (CF from R2007 Table I; A=3.36e-3 meV;
                          J12=-0.1e-3 meV; J0eff=6.421e-3 meV; Jxx0=3.512e-3 meV)
invz_single_ion.m         H_CF + A I·J + Zeeman + transverse MF fixed point;
                          exact diag, 17 (hyp=false) / 136 (hyp=true) states.
                          si.Mx/My/Mz are ELECTRONIC J ⊗ 1_nuclear in eigenbasis.
invz_chi0z.m              Kubo spectral χ0(z), 3×3, arbitrary complex z;
                          elastic (degenerate+diagonal, β-weighted) term added
                          only at |z|<ztol. Real axis: opts.elastic=false.
invz_matsubara.m          bosonic grid wn=2πn/β, n≥0, weights [1;2;2;...].
invz_twolevel.m           electronic doublet params Δ, M², m, n01, g0
                          (errors on Δ<1e-4: invz:degenerateDoublet; |m|>1e-3:
                          invz:orderedPhase).
invz_g.m                  g(z) = 2 n01 Δ/(Δ²−z²).
invz_jq_modes.m           4×4 cc sublattice coupling (recycles repo-root
                          MF_dipole.m + exchange.m), eigenvalue branches
                          Jnu [nq×4] meV; self-verifying cache in invz/cache/.
invz_sigma_crit.m         closed-form Σc(0) = mean(J/(J0−J))  (R2007 eq 10).
invz_emt_scalar.m         K(iωn) fixed point: G=G0/(1+Σ+KG0); Gq=G/(1+(J−K)G);
                          K=⟨J·Gq⟩/G  (R2007 eqs 7-9).
invz_lambdas.m            λp = (1/β)Σn wn K gᵖ           (HTML eq 21 / J 2.19)
invz_sigma.m              α (HTML 27/J 2.17), γ (HTML 23/J 2.18), Σ=α+γg.
invz_solve_point.m        outer loop at (T,Bx): full 136-state G0cc(iωn) +
                          two-level Σ ↔ EMT K; returns pt.{Sigma0,alpha,lambda,
                          K,G,Sigma,tl,chi0cc0,crit,sumrule_rel,converged}.
invz_critical.m           bisection of pt.crit over Bx; NON-FINITE crit ⇒
                          ordered side (para EMT fixed point diverges there).
invz_critical_T0field.m   zero-field Tc from 1+Σc = J0·χ0cc(0;T) (elastic-
                          dominated; uses the closed-form Σc, not α+γg).
invz_chi_realaxis.m       Σ(ω) = α + (M²/n01²)[λ1 −(1−n01²)K(ω)]·g(ω+iη);
                          χ̃0 = χ0cc/(1+Σ); mode-RPA χ = χ̃0/(1−Jν χ̃0).
invz_run_phase_diagram.m  driver: boundary Bc(T) + zero-field Tc (~1-2 h full).
invz_run_spectra.m        driver: χ″(ω) RPA vs 1/z at several fields, 0.31 K.
```

Key design decisions (from R2007, section II):
- MF/RPA layer is "just a bigger Hilbert space": the full electronuclear χ0
  feeds G0 and K, so hyperfine enhancement enters K non-perturbatively.
- Σ is computed from the ELECTRONIC two-level (split doublet) parameters and
  divides the whole cc component (doublet↔higher-CF contributions to cc are
  <1%, R2007's own justification).
- The 1/z lattice sums run over the 4 eigenvalue branches of the 4×4 cc
  sublattice coupling matrix (scalar theory); the tensor 4×4 RPA observable
  layer is future work.

## 3. Validation status (all in the test suite)

| Benchmark | Published | This module | Test |
|---|---|---|---|
| fcc Watson Σc | 0.3447 | 0.34465 (Richardson 40³/80³) | test_invz_sigma_crit (fast) |
| LiHoF4 Σc(0) | 0.3004 | 0.2980 (Richardson 12³/24³) | test_invz_sigma_crit (slow) |
| 𝒥_D·D_aa/cc(0) | 3.912 / 6.821 μeV | 3.910 / 6.824 | test_invz_jq_modes (fast) |
| Tc(B=0), 1/z | 1.74 K | 1.74 K (MF: 2.26 K) | test_invz_critical (slow) |
| Hc(0.31 K) | 42.4–43 kOe | in [4.0,4.6] T | test_invz_critical (slow) |
| Σ(0) at that boundary | 0.0932 | 0.0932 ± 0.02 | test_invz_critical (slow) |
| Soft-mode minimum | ≈0.19 meV (calc) | 0.22 meV | test_invz_chi_observable (slow) |
| CF level scheme | Table II: 0*,11,32,72,84 K | 0*,10.8,32.2,72.0,83.6 K | test_invz_single_ion |
| g∥ | 13.78 | 13.78 | test_invz_single_ion |
| Jensen α/γ formulas | — | 1e-9 match vs old src/ implementation | test_invz_sigma (auto-skips, see §5) |

Run commands (MATLAB is NOT on PATH on this machine):
```bash
cd "<repo root>"    # path contains spaces — always quote
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
# fast: ~2-3 s of test time, expect 26 passed / 0 failed / 5 filtered
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "setenv('INVZ_SLOW','1'); results = runtests('invz/tests'); assertSuccess(results)"
# slow: ~1-2 min with warm invz/cache/, ~15 min cold (dipole grids recompute)
```

## 4. Numerical landmines discovered this session (do not rediscover)

1. **BZ averages of J/(J0−J) are O(1/n)-biased** by the integrable ~1/q²
   singularity at Γ (raw 16³ gives 0.272 for a converged 0.298). Always
   Richardson-extrapolate grid pairs: Sc ≈ 2·S(2n) − S(n).
2. **Σ(0) is NOT monotonic on cooling at fixed field.** It rises from its
   high-T decay, peaks near k_B·T ~ Δ, and saturates to a finite zero-point
   value (α → M²[λ2 − (g0/2)λ1] as T→0). Don't write tests or convergence
   heuristics that assume monotonicity.
3. **Inside the ordered phase the paramagnetic EMT fixed point does not
   exist** — invz_solve_point returns non-finite values with converged=false.
   invz_critical uses this as the ordered-side signal. Always check
   pt.converged when calling solve_point directly.
4. **Matsubara truncation tails**: (1/β)Σg = 1 converges like
   (n01·Δ·β/π²)/nmax — grid identity tests need Ecut ≳ 2000·Δ for 1e-3
   accuracy. λ2-type sums (ωn⁻⁴) converge much faster.
5. **Jq cache**: keyed by a hash but content-verified (pkey+qvec isequal on
   load). A demonstrated single-precision hash collision (5% J12 retune →
   same key) is caught by the verification. If you change lattice/ion
   parameters and see something odd anyway, `rm invz/cache/*.mat`.
6. `qVec_generator` parameter is `'grid'` (not `'size'`), and it fprintf's
   diagnostics — wrap in `evalc` inside tests.
7. MF_dipole returns MINUS the standard dipole sum ([3,3,4,4], Å⁻³); invz_jq_modes
   flips sign and adds the Lorentz term as a scalar broadcast at Γ-equivalent q
   only; Γ diagnostics use the uniform-mode Rayleigh quotient v=[1 1 1 1]/2
   (max(eig) is WRONG for the aa channel). Full convention block with candidate
   table at the top of invz_jq_modes.m.
8. B64s sign: +0.98e-5 meV works (both signs physical per R2007; at B=0 the
   flip is a basis change — level scheme test pins it).

## 5. Known gaps / deferred items

- `test_matches_existing_src_formula` auto-skips (assumeTrue) because the old
  `src/` EMT solver was removed from the tree (see §6). It re-arms if src/ is
  restored from history.
- `invz_chi_realaxis` opts.Jfull/npass multi-pass branch is implemented but
  untested and never exercised (single-shot default is what all tests/drivers
  use); it doesn't check med.converged. Harden before first real use.
- README soft-mode row: computed 0.22 meV vs published calc 0.19 — R2007
  itself needed a ×1.15 energy scale factor vs experiment; treat sub-20%
  agreement here as expected, not a bug.
- Small-Bx caveat (R2007's own): when the doublet splitting ≲ hyperfine width
  (~1.5 K), the two-level determination of Σ overestimates the Tc
  suppression (their comparison vs QMC). The Bx=0 endpoint uses the
  closed-form Σc route instead.
- sumrule_rel ~ 0.05 is EXPECTED (the Dyson-resummed G obeys the sum rule
  only approximately — Jensen's own remark after J 2.23). Don't chase it to 0.
- Minor code polish never done (harmless): unused gL in invz_const; inner
  per-z loop in invz_chi0z unvectorized; a few tautological test assertions
  documented in `.superpowers/sdd/progress.md` (local only, gitignored).

## 6. Repo state / history

- Branch `invz-1z-lihof4`, 23 commits over `main` (023597e), all pushed.
- The old real-frequency "Track-A" EMT stack (effective_medium.m, src/,
  tests/, plot_comparison.m, sum_rule_check.m, SUM_RULE_IMPLEMENTATION.md)
  was committed at **bee3ca2** and then REMOVED at 1e8dc7d — recover any of it
  with `git show bee3ca2:src/emt_solve_point.m` etc. Its solver had a real
  pre-existing defect: the noninteracting-limit test failed (0/2 points
  converged, residual 0.015 vs 1e-7 tol).
- Removed earlier, backed up only transiently (recreate from bee3ca2^ history
  if ever needed): effective_medium_realw.m, "effective_medium-Yikai's
  MacBook Pro.m".
- Kept as the RPA baseline / tensor-observable layer: MF_RPA_Yikai.m (loads
  eigen-data from the external LiReF4_MF_Yikai pipeline), RPA_line.m,
  ellipsoid_demagn.m.
- `invz/cache/` and `.superpowers/` are gitignored (cache is warm on this
  machine only). `.superpowers/sdd/progress.md` holds the full task-by-task
  development ledger with every review finding.

## 7. Where to go next (in rough order of value)

1. **Physics production runs**: `invz_run_phase_diagram` (full boundary,
   ~1-2 h) and `invz_run_spectra`; compare against R2007 Figs 1-2 and
   Kovacevic Fig 3(d). The GHz-frequency χ″ maps need ω grids in the
   1.7–5.6 GHz range (× C.Gh2mV to meV) and the electronuclear elastic-sector
   physics is already in χ0.
2. **Ordered phase** (the big one): implement the m≠0 elastic sector —
   HTML eqs 37–40 (α_m with λ3, the ξ tanh resummation) and the modified
   mean field H_MF via HTML eqs 41–43 (J 2.31–2.33). Entry points:
   invz_twolevel currently ERRORS on |m|>1e-3 — that guard marks exactly
   where the extension goes; invz_solve_point's non-finite behavior in the
   ordered phase is the other boundary.
3. **Tensor observables**: feed χ̃0 (cc ÷ (1+Σ)) into a 4-sublattice 3×3 RPA
   (MF_RPA_Yikai's machinery is the template) for neutron cross-sections /
   transverse components.
4. **Thermodynamics**: δU and heat capacity via HTML eqs 34–35 (J 2.21–2.22)
   and the free-energy consistency test eq 44 (J 2.34) — note δU(T→0) ≠ 0
   (zero-point shift; Jensen quotes −9.06 μeV/ion for HoF3).
5. **Performance**: Ewald summation to replace the brute-force dipole sums if
   denser q-grids are needed; parfor over q in invz_jq_modes.
