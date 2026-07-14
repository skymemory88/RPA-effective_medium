# Session note — 2026-07-14: demag knob rework + exploratory q-path dispersion

Branch: `invz-1z-lihof4`, commits `8b0854e`…`b435377` (+ final cleanup commit after this note).
Plan of record: `docs/superpowers/plans/2026-07-14-demag-knob-qpath-spectra.md` (rev 2, incorporating
`q-spectrum_plan_review_by_Codex.md`). Executed task-by-task with per-task spec + quality reviews.

## What changed (physics)

**Demagnetization is now handled per Rønnow et al. PRB 75, 054426 (2007)** — fixing
`Code_review_byCodex.md` Finding #1 (High), where the draft folded the sample shape into
`info.Jcc0` and thereby (wrongly) shifted Bc/Tc through the critical condition. The canonical
semantics, used verbatim in all comments/README:

1. `info.Jcc0`, `Jnu`, and the ordering-channel contribution to criticality are **demag-invariant**
   (R2007: the demag field cancels from the critical condition; ordering at q→0⁺, not strict q=0 —
   stated by R2007 *against* a suggestion of Chakraborty et al., so cite "per R2007").
2. `Tc(B=0)` is **exactly** demag-invariant (transverse moment vanishes there).
3. `Bc(T)` vs **applied** transverse field **can shift** through the demag-aware `info.Jaa0`
   (internal-vs-applied transverse field relation).
4. q-path spectra omit the strict-uniform `Jshape_cc` transform (finite-q probe = intrinsic) but
   still see demag through `info.Jaa0`.

Mechanics: `invz_jq_modes` exports `info.Jshape_cc = 4*dm_cc` (applied once in `invz_chi_realaxis`
as `chi_meas = chi/(1+Jshape_cc*chi)` — algebraically identical to the old shifted pole, so the
observable is unchanged; only criticality/Σ/mean-field moved) and demag-aware
`info.Jaa0 = Jaa0_dipole + 4*J12`, threaded as `opts.Jxx0` through every single-ion touchpoint
(`invz_solve_point[_ordered]`, `invz_twolevel[_ordered]`, `invz_chi_realaxis`). The Lorentz cavity
term stays unconditionally on (mandatory dipole-sum-split term, not a physical toggle — decided,
no knob). `ion.Jxx0` remains the fallback default; drivers now hoist the live `info.Jaa0`
(≈3.5104e-3 vs hardcoded 3.512e-3 meV at dpRng 30 — intended ~0.05% transverse drift, absorbed by
test tolerances).

**User knob** (both drivers, hand-edited script variables per house style):
`ion.demag = 1; ion.alpha = 0.5;` right after `ion = invz_ion()` (0=off/intrinsic default;
alpha: 1 sphere, 0 c-needle, Inf disk).

## What changed (new capability)

**Exploratory q-path dispersion** (R2007 Fig. 3 *trends*, not quantitative reproduction):
- `invz_jq_path` — path couplings with a direction-aware Γ-limit guard: within trust radius
  `2.5·2π/(dpRng·min‖aᵢ‖)` of a Γ-equivalent point, the truncated `MF_dipole` sum (which collapses
  approaching e.g. (2,0,0) — measured 0.0016 vs correct 0.0064 meV at h=1.999, dpRng 30) is
  replaced by `eig(J_reg(Γ) + gfac·(4π/Vc)·(1/3−k̂_z²))`. In-plane approach ⇒ uniform-mode Lorentz
  value = `info.Jcc0` (continuous endpoint by construction). Guard covers Γ-equivalent points ONLY.
- `invz_spectra_qpath` — one 1/z medium solve per (T,B) via new `invz_solve_auto` (ordered-first,
  `invz:*`-only catches with diagnostics, rethrows programming errors), then the vectorized `Jsel`
  pole formula along the path. Censored, parabolic-refined peak picker (`peak_wmin=0.05` meV skips
  the hyperfine pole; boundary maxima → NaN, never reported).
- `invz_plot_spectra_qpath` + third view in `invz_run_spectra.m`: knobs `qpath`, `Bq` (scalar →
  colormaps; vector → E_peak overlay, Fig. 3 style), `wq` (to 0.85 meV — Fig. 3 reaches ~0.75),
  `dispScale` (R2007 scale theory by 1.15).
- Output is a **branch susceptibility** (sorted-eigenvalue index, no crossing tracking), NOT
  neutron intensity; inherits Codex #2 (BZ quadrature), #3 (real-axis positivity), #4 (bare-MF
  FM/PM handoff) — all documented in headers/README.

## Verification

- Fast suite: 42 passed, 0 failed, 9 filtered. Full `INVZ_SLOW=1` suite: **51 passed, 0 failed**
  (843 s wall, run twice independently).
- New regression locks: demag invariance of Jnu/Jcc0/Σ_c/Tc0 (bit-identical); pinned-`Jxx0` crit
  invariance + unpinned transverse shift (SLOW); Γ-approach cutoff regression (dpRng 20 vs 40,
  h=1.90…2.0); anisotropic k̂_z² limit (out-of-plane vs in-plane closed-form shift); peak
  censoring; observable-rescale algebraic identity; `invz:*`-only catch/rethrow policy.
- `invz_jq_modes` cache schema bumped `jq_` → `jq2_`: old `invz/cache/jq_*.mat` are orphaned and
  safe to delete.

## Deferred / out of scope (unchanged from plan rev 2)

Ewald summation (guard is the sanctioned deferral), off-diagonal dipolar (ODD/full-tensor) terms
(cf. Dollberg et al. PRB 105, L180413 — invz's cc/aa scalar decoupling misses these; legacy
`MF_RPA_Yikai.m` keeps the full tensor), aa-channel observable, eigenvector-overlap branch
tracking, Codex #2/#3/#4.

## Gotchas discovered this session

- **OneDrive + git**: OneDrive created a conflict-copy branch ref
  (`invz-1z-lihof4-Yikai's MacBook Pro`, pointing at our own Task 6 commit). Deleted via
  `git branch -D`; harmless, but expect recurrence while `.git` lives inside OneDrive.
- `ion.renorm` in `MF_RPA_Yikai.m` is the manual MF spin-length renormalization from R2007 —
  deliberately absent from invz (the 1/z expansion supplies it).
