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
- **#4 (High) — DONE (2026-07-15).** `MF_dipole`/`exchange` gained an optional precomputed-geometry
  arg (5th/5th) and return the geometry struct; `invz_jq_modes` builds the q-independent lattice
  geometry (reciprocal basis, meshgrid/`r0`, per-pair cutoff `r`, `r^-3`/`r^-5`, the 9 tensor-factor
  columns, `nN`) ONCE and reuses it across the whole q-loop and the Γ `dip0`. Cold-grid speedup
  without touching the cache key or any numerics. BIT-IDENTICAL verified against a golden from the
  unmodified code: max abs diff 0 for `dip`, `exchange`, `Jnu`, `Juni`, `Jcc0`, full `info`
  `isequal`. The full 3×3 Cartesian tensor output is UNCHANGED (kept for ODD, like #3 — the review's
  "scalar cc dipole path" sub-item was NOT taken, as ODD needs all blocks). External callers
  (`MF_RPA_Yikai`, `RPA_line`, `invz_jq_path`) use the 4-arg form and are untouched. Test:
  `test_invz_dipole_geometry_reuse` (4-arg == 5-arg-reused-geom over 7 q). The `exchange` cutoff
  quirk (`max(vecnorm(a,1).^2)` vs `max(vecnorm(a,1)).^2`, equal by coincidence for this lattice)
  is preserved verbatim.
- **#5 (Medium) — OPT-IN PROVIDED (2026-07-15).** `qVec_generator` gained an `endpoint` option
  (default `true` = the historical inclusive `linspace`, byte-identical). `endpoint=false` gives the
  half-open grid (`lo + (0:N-1)/N*(hi-lo)`) with no duplicate ±0.5 face. NOT flipped by default:
  the LiHoF4 benchmarks (Tc0≈1.7795 K on the 16³ grid, Σ_c, Bc) are calibrated on the inclusive
  grid, so changing the quadrature shifts every grid-dependent result and requires re-validation
  against R2007 — a physics decision left to the user. Test `test_qVec_generator` pins both the
  unchanged default and the half-open point set. Irreducible-BZ-with-weights left as future work.
- **#6 (Medium) — PARTIAL (2026-07-15).** DONE: `invz_solve_point` accepts `opts.Sigma_seed`/
  `opts.K_seed` (same length as `wn`; only changes the iteration path, not the converged fixed
  point — `test_warm_start_same_fixed_point` proves same result + strictly fewer outer iters);
  `invz_critical_T0field` bisection now stops at bracket width 1e-9 K (~31 vs the old fixed 60
  136-state solves). DELIBERATELY NOT DONE: wiring warm-start into `invz_critical`/`invz_critical_T`
  boundary classification — near critical slowing-down a warm seed can flip a point from (cold)
  non-converged to (warm) converged, shifting the classified `Bc` and breaking the field-cut vs
  temperature-cut mirror consistency (`test_tc_at_fixed_field_crossing`); the ordered-phase
  non-convergence signal these searches depend on must stay seed-independent. Anderson/quasi-Newton
  and the multiple-crossing-audit option (bullet 6) not pursued (the full-grid re-entrance scan in
  `invz_critical_T` is intentional).
- **#7 (Medium) — PARTIAL (2026-07-15).** DONE: the pasted `ion.J0eff`/`ion.Jxx0` are relabelled
  as explicit REFERENCE/FALLBACK values with the authoritative live `info.Jcc0`/`info.Jaa0`
  (from `invz_jq_modes`, dpRng-dependent) called out as the source of truth. DELIBERATELY NOT DONE:
  converting the drivers (`invz_run_spectra`, `invz_run_phase_diagram`) into callable functions with
  a validated `invz_defaults()` cfg struct — the house style is hand-edited script variables (see
  2026-07-14 note), and an `invz_defaults()` nobody wires in would be dead code. `arguments`-block
  option normalisation not pursued (would touch every solver signature; deferred with #8).
- **#8 (Medium) — PARTIAL (2026-07-15).** DONE: the duplicated local `getf` in `invz_spectra_map`
  and `invz_spectra_qpath` is centralised to `invz/getf.m` (both local copies removed;
  `test_getf_shared_helper`). DEFERRED (review's own guidance — do these as a dedicated refactor
  now the numerics are pinned): the `invz_solve_point[_ordered]` prepared-core extraction, the
  `invz_spectra_*`/`invz_run_phase_diagram` shared lattice-prep, the two `invz_plot_spectra_*`
  helpers, and `invz_twolevel[_ordered]` dedup — each carries numerical-drift risk that outweighs
  the readability gain in this pass.
- **#9 (Low) — WAIVED (2026-07-15).** Not trimming the in-code physics narrative. This is research
  code whose long comments (the `invz_jq_modes` sign-resolution + candidate table, demag semantics,
  critical-classifier history) are load-bearing documentation tightly coupled to subtle lines;
  moving them to `README.html` reduces discoverability for whoever maintains the kernels. The review
  rates this "Low"; per code-review discipline this is a reasoned decline, not an omission.
- **#10 (Low) — MOSTLY DONE (2026-07-15).** DONE: removed unused `gL` in `invz_const.m` (Analyzer
  clean); `qVec_generator.m` gained `verbose` (default true), gated all diagnostics, removed the dead
  `display_qvectors` helper + its commented call, and the three production callers
  (`invz_spectra_map`/`invz_spectra_qpath`/`invz_run_phase_diagram`) now pass `verbose=false`;
  deleted orphaned `invz/cache/jq_*`/`jq2_*` (kept current `jq3_`). NOT DONE (deliberate): the
  `invz_plot_spectra_qpath` pre-`S.x` compat branch is KEPT — it is covered by
  `test_qpath_plot_coordinate_fallback`, so removing it would drop tested behaviour. Phase-code
  named constants and the terse-name renames (`si`/`tl`/`sg`/…) were skipped as high-churn,
  low-value edits with regression risk and no functional benefit.

### Codex efficiency-review — final disposition (2026-07-15)

All 10 findings triaged and actioned in this session (commits on `invz-1z-lihof4`):

| # | Sev | Status | Note |
|---|-----|--------|------|
| 1 | Critical | DONE | closed-form EMT (~84× kernel), non-finite guard |
| 2 | High | DONE | share si/tl/chi0cc per field point (bit-identical) |
| 3 | High | WAIVED | full chi0z tensor is required ODD substrate |
| 4 | High | DONE | dipole/exchange geometry reuse (bit-identical) |
| 5 | Medium | OPT-IN | half-open grid available (`endpoint=false`), default unchanged |
| 6 | Medium | PARTIAL | seed capability + T0field early-stop; not wired into boundary search (physics) |
| 7 | Medium | PARTIAL | couplings relabelled fallback; driver→function declined (house style) |
| 8 | Medium | PARTIAL | `getf` centralised; large refactors deferred (numerics now pinned) |
| 9 | Low | WAIVED | in-code physics narrative kept (research-code asset) |
| 10 | Low | MOSTLY | gL, qVec verbose, dead code, cache cleanup; compat branch kept (tested) |

"WAIVED/PARTIAL/OPT-IN" items each carry a documented rationale above; none is a silent omission.
Every code change was gated behind the full `INVZ_SLOW=1` suite (0 failures throughout).

## Gotchas

- **OneDrive + git** (recurring, see 2026-07-14 note): `.git` lives inside OneDrive; watch for
  conflict-copy branch refs.
- Removing the standalone `Code_efficiency_review_by_codex.md`: its content now lives only here.
