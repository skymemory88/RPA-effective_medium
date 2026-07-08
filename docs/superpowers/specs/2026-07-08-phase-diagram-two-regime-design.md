# Two-regime phase-boundary search for `invz_run_phase_diagram`

Date: 2026-07-08. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-08); AMENDED the same day during implementation —
see Amendments below, which govern wherever later sections conflict.

## Amendments (2026-07-08, user-approved during implementation)

Diagnostics during Task 1 (the crossing-consistency test passed on the first
run; convergence scans at B = 0.2/0.3/0.5 T) established:

1. At small fields the paramagnetic solve develops non-convergence patches
   near the boundary; the non-finite ⇒ ordered classifier reads them as
   ordered and biases Tc(B) upward (+0.04–0.05 K at 0.2–0.3 T; gentler
   `max_outer`/`mix_outer` does NOT cure them). At B ≥ 0.5 T the method is
   clean: Tc(0.5 T) = 1.777 K, just below the closed-form Tc0 = 1.7795 K
   evaluated on the same 16³ q-grid. (The spec's 1.74 K baseline below is
   the Richardson-extrapolated value — the wrong comparison for a 16³
   computation; the undershoot expectation was right once compared
   grid-consistently.)
2. `Bs` is therefore floored at 0.5 T (default `[0.5 0.75 1.0 1.25 1.5]`).
   The 0 < B < 0.5 T boundary segment spans only ~4 mK in temperature and
   is represented by the closed-form Tc0 endpoint on the plot.
3. The `Tsplit` knob is REMOVED (user decision: trimming `Ts` is redundant —
   one can simply shorten `Ts`). The Tc(B) bisection window is fixed inside
   the driver as commented constants `Tlo = 1.0` K / `Tmax = 2.0` K, with
   the documented constraint `max(Bs) < Bc(Tlo) ≈ 2.8 T`. The script never
   modifies `Ts`; the default `Ts` now ends at 1.6 K.
4. The near-zero-field test becomes `test_tc_small_field`: B = 0.5 T,
   window [1.0 2.0] K, bounds (1.70, 1.79) K.

## Problem

`invz_run_phase_diagram.m` traces the LiHoF4 paramagnetic-side phase boundary
by bisecting the critical field Bc(T) at each temperature in `Ts`
(`invz_critical`: sign change of `pt.crit` from `invz_solve_point` over Bx).
Near the classical critical point (T → Tc0 ≈ 1.74 K, B → 0) the boundary is
nearly parallel to the field axis, so a fixed-T (vertical) cut crosses it at a
glancing angle: tiny temperature errors become huge Bc errors, brackets become
hard to establish, and the last few `Ts` points (the 1.64–1.8 K cluster) are
ill-conditioned or fail outright.

## Solution overview

Split the boundary trace into two regimes at a user-facing knob `Tsplit`:

- **Low-T regime (T ≤ Tsplit)**: unchanged — bisect Bc(T) at fixed T with
  `invz_critical`.
- **High-T regime (boundary above Tsplit)**: invert the roles — at each fixed
  field in a new `Bs` list, bisect the critical temperature Tc(B) over T with
  a new mirror function `invz_critical_T`. A fixed-B (horizontal) cut crosses
  the near-vertical boundary transversally, so it is well-conditioned exactly
  where the fixed-T cut is not.

The ordered/paramagnet classifier is identical in both directions: `pt.crit`
non-finite or ≤ 0 ⇒ ordered side (the paramagnetic EMT fixed point does not
exist inside the ordered phase). Orientation is also the same — bracket low
edge ordered, high edge paramagnetic — so the bisection loop mirrors
`invz_critical` line for line.

## Components

### 1. New function `invz/invz_critical_T.m`

```
Tc = invz_critical_T(ion, Bx, Jnu_flat, opts)
```

- Bisects `pt.crit` (via `invz_solve_point(ion, T, Bx, Jnu_flat, opts)`) over
  temperature at fixed transverse field Bx.
- `opts.window` = `[Tlo Thi]` temperature bracket, default `[1.0 2.0]` K.
- `opts.tol` = bisection resolution in K, default `0.01`.
- Bracket assert: requires ordered at `Tlo` (crit non-finite or ≤ 0) and
  paramagnetic at `Thi` (crit finite and > 0). The error message states the
  likely cause: Bx exceeds Bc(Tlo) — lower `max(Bs)` or `Tlo`.
- All other `opts` fields pass through to `invz_solve_point` (`J0eff`, `hyp`,
  `Ecut`, …), exactly as `invz_critical` does today.

`invz_critical.m` is NOT modified (it is validated by the slow test suite).

### 2. Driver changes `invz/invz_run_phase_diagram.m`

Knobs at the top of the script:

```matlab
Tsplit = 1.5;   % regime split: boundary below Tsplit via Bc(T), above via Tc(B)
Ts = [0.05 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.64 1.68 1.72 1.74 1.76 1.78 1.8];
                % master temperature list; auto-trimmed to Ts(Ts <= Tsplit)
Bs = [0.1 0.25 0.5 0.75 1.0 1.25 1.5];  % fields for the high-T Tc(B) regime
Tmax = 2.0;     % paramagnetic upper bracket edge for Tc(B) (above Tc0 = 1.74 K)
```

- The master `Ts` list is kept in the script and trimmed with
  `Ts = Ts(Ts <= Tsplit)`; the current 1.64–1.8 K cluster (the ill-conditioned
  points) thereby moves to the new regime automatically.
- One flat `parfor` over all `numel(Ts) + numel(Bs)` jobs — every point in
  both regimes is an independent bisection. Jobs `1..nT` run the unchanged
  Bc(T) path (window `[0.5 7]` T as today); jobs `nT+1..end` run
  `invz_critical_T` with window `[Tsplit, Tmax]`. Each job keeps the existing
  try/catch-to-NaN pattern so one failure never kills the sweep.
- Bracket validity is geometric: at `T = Tsplit` a point is ordered exactly
  when `B < Bc(Tsplit)`, so every `Bs` entry below `Bc(Tsplit)` brackets
  cleanly; an entry above it fails its bracket assert cleanly to NaN without
  affecting other jobs.
- The zero-field endpoint is unchanged: closed-form
  `invz_critical_T0field(ion, invz_sigma_crit(...), J0)`.
- Header comment documents why `Bs` starts at 0.1 T:
  1. `invz_twolevel` raises `invz:degenerateDoublet` when the field-induced
     doublet splitting Δ < 1e-4 meV (Bx ≈ 0);
  2. R2007's small-Bx caveat — the two-level Σ overestimates the Tc
     suppression when the doublet splitting ≲ the hyperfine width, so
     Tc(B → 0) slightly undershoots the closed-form Tc0. Expected physics,
     not a bug.
- Plot: both branches on one axes — `plot(Ts, Bc*10, 'o-')` for the low-T
  branch, `plot(TcB, Bs*10, 's-')` for the high-T branch — plus the Tc0
  black-square marker, as today. Additionally assemble a merged, T-sorted
  `phase_boundary = [T, B]` array in the workspace for downstream use (named
  to avoid shadowing MATLAB's built-in `boundary`).

### 3. Cost

Per-point cost is essentially unchanged: tol 0.01 K over a ~1 K window is
~9–10 EMT solves per point, comparable to tol 0.02 T over the 6.5 T field
window. Total runtime scales with `numel(Ts) + numel(Bs)`; parallel speedup
is near-linear up to the job count as before.

## Error handling

- Degenerate doublet at tiny Bx: guarded by the 0.1 T minimum in the default
  `Bs` (documented in the driver header).
- `Bs` entry above `Bc(Tsplit)`: clean bracket-assert failure → NaN for that
  point only, with a hint in the message.
- Non-convergent EMT solves inside the ordered phase are the classifier
  signal, not an error (unchanged behaviour inherited from `invz_critical`).

## Testing

New test points appended to `invz/tests/test_invz_critical.m` (reusing its
`lihof4_couplings` helper; both SLOW-gated with `INVZ_SLOW=1` like the
existing critical-point tests):

1. **Crossing consistency** — at T* = 1.4 K (mid-slope, where both cut
   directions are well-conditioned): compute `B* = invz_critical(ion, 1.4,
   Jf, ...)`, then `invz_critical_T(ion, B*, Jf, struct('window',[1.0 2.0],
   ...))` must return 1.4 K within 0.05 K. Validates the mirror against the
   already-validated function with no new reference data.
2. **Near-zero-field sanity** — `invz_critical_T` at Bx = 0.2 T returns Tc in
   (1.55, 1.76): at or below Tc0 = 1.74 K, pinning the expected small-Bx
   undershoot direction.

## Success criteria

- Fast test suite still passes with no failures (26 passed; filtered count
  grows from 5 to 7 because the two new tests are SLOW-gated).
- Both new slow tests pass with `INVZ_SLOW=1`.
- Driver's high-T branch returns finite Tc(B) for every default `Bs` entry,
  monotonically decreasing with B.
- Merged boundary joins smoothly at (Tsplit, Bc(Tsplit)) and reproduces the
  R2007 Fig 1 shape, including the near-vertical drop to (1.74 K, 0).

## Out of scope

- Slanted/normal cuts or predictor–corrector boundary continuation (rejected
  alternative — serial, more code, overkill).
- Generalising `invz_critical` into a line-segment bisection (rejected — it
  would refactor a validated function for no current need).
- Any change to the solver layer (`invz_solve_point` and below).
