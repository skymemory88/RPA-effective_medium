# Session note — 2026-07-15: EMT closed-form solve + Codex efficiency-review triage

Branch: `invz-1z-lihof4`. This session evaluated `Code_efficiency_review_by_codex.md`
(2026-07-14 Codex code-efficiency review of `invz/`, 10 findings) and implemented its top
item. That standalone review file is **deleted** in this session's commit; its findings are
preserved below so nothing is lost — findings #2–#10 remain valid, unaddressed backlog.

## What changed (finding #1 — closed-form EMT)

`invz/invz_emt_scalar.m`: replaced the 35-iteration effective-medium fixed-point loop with the
exact closed form. The CPA closure `K = (1/N)Σ_q J·Gq/G` with `Gq = G/(1+(J−K)G)` and
`G = G0/(1+Σ+K·G0)` is exactly solvable because **K cancels out of Gq**. With `D = 1+Σ`:

```
Gq = G0/(D + J·G0)                 (K-free)
A  = (1/N) Σ_q  J/(D + J·G0)       (K-independent)
K  = A·D/(1 − A·G0)                (linear in K → exact)
G  = G0/(D + K·G0)
```

Re-derived independently and confirmed: at this fixed point R eq 8 also gives `mean_q Gq = G`
exactly, so the closure diagnostic is machine-zero when a solution exists.

Two improvements beyond speed:
- **Blockwise frequency evaluation** (`opts.freq_block`, default 4096) so the `[nJ × nw]`
  implicit-expansion temporary never materialises at full width — this, not the closed form
  alone, is what delivers the review's low-T memory win (185 MiB at nw=… collapses to a block).
- **Explicit finiteness guard** on `med.converged = all(isfinite(G)) && all(isfinite(K))`.
  Fixes a **latent bug in the old code**: at a singular medium (`1 − A·G0 → 0`, i.e. a vanishing
  RPA denominator) it reported `converged = true` because MATLAB `max()` silently ignores the
  `NaN` residual. `invz_solve_point`/`invz_critical` rely on the singular case reading as
  "no paramagnetic solution" (the ordered-phase signal), so the guard is load-bearing, not cosmetic.

`opts.tol/max_iter/mix/K0` are still accepted (no caller break) but are now no-ops.

## Verification

- **New unit tests** in `invz/tests/test_invz_emt_scalar.m` (3 added):
  `test_solution_is_exact_fixed_point` (eq-8 global-relative residual + closure < 1e-12),
  `test_matches_reference_iteration` (agrees with an independent tight iteration across coupling
  factors 0.5/0.95/1.6 incl. super-critical), `test_nonfinite_denominator_not_converged`
  (singular `G0=−1/J` ⇒ not converged). The exactness and non-finite tests were proven RED
  against the old implementation and GREEN against the new (git-stash A/B); the other two are
  characterisation locks.
- Exactness measured: eq-8 residual `8e-16`, closure `1e-15` (machine precision). The old
  per-entry-floored residual metric was replaced by a global-relative one because it is
  meaningless where `K → 1e-7` at high Matsubara frequency (test metric fixed, not the impl).
- **Fast suite: 49 passed, 0 failed, 9 filtered. Full `INVZ_SLOW=1` suite: 58 passed, 0 failed**
  (596 s). The slow set includes the deep near-critical integration tests
  (`test_invz_critical` Tc0 + Bc@310mK, `test_ordered_phase`, `sigma_crit/lihof4`, `soft_mode`,
  demag transverse) that the original review explicitly did **not** run — closed form + guard
  reproduce downstream physics exactly.
- Kernel speedup (16³ grid, nw=240): **~84×** (35 iters → 1), old/new agree to `5e-12`.
  End-to-end is smaller: EMT was ~0.68 s of a ~1.1 s point solve, so this alone takes the solve
  to ~0.44 s (~2.5×); the larger point-solve wins come from findings #2/#3 below.
- `checkcode` clean on `invz_emt_scalar.m` (removed now-unused `nJ`).

## Codex efficiency-review backlog (findings #2–#10, NOT yet done)

Preserved verbatim-in-substance from the deleted `Code_efficiency_review_by_codex.md`. Spot-checked
against source this session: #2, #3, #4, #10(gL) confirmed accurate; #1's algebra proven exact.

- **#2 (High) — DONE (2026-07-15).** Paramagnetic spectra repeated expensive single-ion work.
  `invz_solve_point` now returns `pt.si` (matching the ordered solver). `invz_chi_realaxis`
  accepts `opts.chi0cc_w` (a precomputed bare cc), so the RPA-overlay and 1/z calls of one field
  point share a single `chi0cc` instead of each rebuilding the single ion + `invz_chi0z`.
  `invz_spectra_map/one_field` reuses `pt.tl`/`pt.si` on the converged paramagnet and threads the
  overlay's `chi0cc_w` into the 1/z call (both branches). Output is bit-identical (regression:
  full `INVZ_SLOW=1` suite 60/60). Tests: `test_returns_single_ion_state`,
  `test_chi0cc_passthrough_is_consumed` (perturbation proves consumption). NOT done: the phase-
  continuation sub-item (skip the ordered-branch test at clearly-PM points) — folded into #6.
- **#3 (High) — WAIVED (2026-07-15).** `invz_chi0z`'s full 3×3 tensor is a required substrate for
  the planned ODD work (`odd_implementation_plan.html`: Tier-1 `invz_chiperp` extracts the `(a,b)`
  transverse block; the full-tensor RPA parity layer needs all 9 Cartesian blocks; the plan is
  written against the current source). An opt-in `comp` fast path was prototyped and reverted per
  user direction — leave the full-tensor kernel alone. If the ~9× production `cc` speedup is later
  wanted, an *additive* `opts.comp` (default = all 9) is ODD-compatible.
- **#4 (High)** Cold q-space construction rebuilds q-invariant geometry per q-point. `MF_dipole`
  rebuilds lattice grid, displacements, `r^-3/r^-5`, cutoff masks, and all 9 Cartesian components
  every call; only `exp(-iq·r)` and the final contraction are q-dependent. ~0.114 s/call ⇒ ~8 min
  for a cold 4096-pt grid. Fix: split geometry prep from Fourier eval; batch q-points; scalar `cc`
  path (+ `aa` only at Γ). Files: `invz_jq_modes.m:116-130`, `MF_dipole.m:25-53`, `exchange.m:30-51`.
- **#5 (Medium)** Reciprocal grid has duplicate periodic boundary faces: inclusive
  `linspace(-0.5,0.5,N)` double-counts the ±0.5 faces (16³ evaluates 4096 pts, only 3375 distinct;
  17.6% redundant). Fix: half-open grid (integer mesh / N); re-confirm quadrature normalisation
  against convergence tests; longer-term irreducible-BZ set with weights. Files:
  `qVec_generator.m:188-197`, `invz_spectra_map.m:50-53`, `invz_spectra_qpath.m:68-71`,
  `invz_run_phase_diagram.m:53-56`.
- **#6 (Medium)** Critical-boundary searches don't warm-start from nearby solutions; every (T,B)
  restarts `Sigma`/`K`/mean fields from defaults, and `invz_critical_T` scans a full coarse grid
  before picking the highest-T crossing. Fix: accept initial state in `invz_solve_point`; carry last
  converged state through scans; tolerance-based termination in `invz_critical_T0field` (fixed 60
  iters → bisect to tol); consider Anderson/quasi-Newton near critical slowing. Files:
  `invz_critical.m:31-68`, `invz_critical_T.m:57-79`, `invz_solve_point.m:22-33`,
  `invz_critical_T0field.m:7-11`.
- **#7 (Medium)** Config knobs scattered (scripts require source edits; grid/cutoff hardcoded above
  the "knobs" section; nonempty `qpath` silently unuses `fields`/`w`; two sources for uniform
  couplings — pasted `ion.J0eff`/`ion.Jxx0` vs derived `info.Jcc0`/`info.Jaa0`). Fix: make drivers
  callable functions returning results; one validated cfg struct (`invz_defaults()`); explicit
  `cfg.mode`; treat uniform couplings as derived; `arguments` blocks. Files: `invz_ion.m:3-24`,
  `invz_run_spectra.m:24-60`, `invz_run_phase_diagram.m:53-79`, `invz_jq_modes.m:134-138`.
- **#8 (Medium)** Duplicated code to consolidate (do AFTER perf fixes so numerics stay pinnable):
  `invz_solve_point[_ordered]` (option parsing / Matsubara / EMT-Σ loop / packing);
  `invz_spectra_map`+`invz_spectra_qpath`+`invz_run_phase_diagram` (q-grid prep / zero-filter /
  `invz_jq_modes` / `Jaa0` fallback); the two `invz_plot_spectra_*` (log scaling / percentile /
  colorbar); duplicated local `getf`; `invz_twolevel[_ordered]`.
- **#9 (Low)** Comments carry heavy historical narrative (`invz_jq_modes.m` ~56% comment,
  `invz_run_phase_diagram.m` ~53%, `invz_critical_T.m` ~37%, `invz_spectra_qpath.m` ~42%).
  Move derivations / rejected alternatives / benchmark tables / bug history into README/design note;
  keep headers to I/O, invariants, numerical limits, non-obvious conventions.
- **#10 (Low)** Legacy/minor: `invz_const.m:6` unused `gL` (Code Analyzer); `invz_plot_spectra_qpath.m:17-20`
  pre-`S.x` compat branch (decide support / isolate in loader); `invz/cache/` `jq_`/`jq2_`/`jq3_`
  schema generations not auto-cleaned; `qVec_generator.m` unconditional diagnostics + unused
  `display_qvectors`; phase codes 0/1/2 want named constants; terse names (`si`,`tl`,`sg`,`med`,
  `Jf`,`c0`) fine in kernels but hurt driver readability.

Recommended order (per review): #1 ✓ → #2 → #6 (warm start) → #4 + #5 → #7 → #8 → #9/#10.

## Gotchas

- **OneDrive + git** (recurring, see 2026-07-14 note): `.git` lives inside OneDrive; watch for
  conflict-copy branch refs.
- Removing the standalone `Code_efficiency_review_by_codex.md`: its content now lives only here.
