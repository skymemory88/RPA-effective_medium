# Session handoff — ordered-side (FM) a1 + dominant-vertex a3d for LiHoF4

**Date:** 2026-07-19/20 · **Branch:** `invz-1z-lihof4` · **Plan:** `docs/superpowers/plans/2026-07-19-invzt-ordered-side.md`
**Status:** Tasks 1–6 and 7A–7D complete (11 of the 12 plan tasks); this note + `invz_tensor/README.html` + `docs/ODD-LOG.md` §A-ordered are Task 8. Final whole-branch review is next.

This note is the human handoff for the ordered-side plan. The module reference (mode ladder, module map, LOCKED conventions, headline results, open items) lives in `invz_tensor/README.html`; the measured physics log is `docs/ODD-LOG.md` §A-ordered. This note adds the *why*, the execution story (the diagnosis that motivated the plan and the four amendments made while executing it), and the operational lessons.

---

## What was built

Two things, mirroring the projected branch's ordered stack onto the tensor lattice engine (`invzt_gcc_lattice` 12×12 RPA):

1. **Ordered-phase `a1`, end to end:** a spontaneous-moment point solver (`invzt_solve_point_ordered`), a real-axis continuation (ordered branch of `invzt_chi_realaxis`), a **stability-based** phase dispatcher (`invzt_solve_auto`), and a rewired `invzt_run_spectra.m` that now spans the quantum phase transition instead of stopping at it.
2. **`a3d`:** the **full-response, fixed-rank field-adapted dominant-vertex ordered a3d** — the full 136-state whole-cc ordered medium with the four-point vertex computed only in a 16-state field-adapted dominant basis, Matsubara-only, compact `cc;cc` storage, hard performance budget guards (Tasks 7A–7D).

Full-dress 136-state ordered `a3` remains permanently out of scope (budget-refused, `invzt:orderedMode`) — it is a different name from `a3d` on purpose (round-2 review minor: never conflate the two).

### 2026-07-20 electro-nuclear spectra correction

The historical driver diagnosis below concerns the broad predominantly electronic branch. The current `invzt_run_spectra.m` intentionally uses `[0, 0.018]` meV to resolve the experimentally measured electro-nuclear soft mode. Two independent near-QCP defects were subsequently fixed:

1. `invzt_solve_auto` now supplies the ordered a1 leg with the boundary-matched linearized Jensen modified field, so its moment vanishes at the PM instability instead of retaining the large bare-MF moment to ~5 T.
2. The PM dominant manifold now uses the fixed-rank lowest 16 electronuclear states. The former `E < 0.4653` meV cut changed rank 11→10→9→8 at 4.65/4.76/4.88 T even though the state-16→17 gap remained ~2.4–2.5 meV; those bookkeeping changes were the source of the repeated unphysical spectral jumps.

Measured after the correction (0.1 K, 16³ halfopen/dpRng-30): ordered through 4.64 T, PM from 4.65 T, zero masked; electro-nuclear peak 0.000113 meV at 4.64 T and 0.000754 meV at 4.65 T; PM-side hardening is smooth to 0.005609 meV at 5.05 T. Tensor CORE: 98 passed / 0 failed / 5 expected `INVZ_SLOW` incompletes.

---

## The driver-knob diagnosis (three measured failures, 2026-07-19)

Before Task 6 rewired `invzt_run_spectra.m`, the driver's PM-only knobs — ported from the projected driver's GHz-unit window — were measured to fail in three independent ways at the ordered-side production anchor (T = 0.1 K):

1. **Window too narrow.** `w` topped out at 0.025 meV while the soft mode sits at 0.26–0.34 meV near B<sub>c</sub> — **97% of the spectral weight fell above the window**.
2. **Aliasing broadening.** `eta = 2e-5` meV sat below the frequency-grid step (1.25e-4 meV): a Lorentzian narrower than the sampling aliases between grid points, measured as a **10× peak-height error**.
3. **Masked boundary.** Fields 4.0–4.65 T were masked by the PM-only solver — the driver simply had no ordered leg to fall back on.

Fixes shipped in Task 6: window widened to `[0, 0.6]` meV (601 points), `eta` raised to 2e-3 meV with a hard guard (`invzt_run_spectra:etaStep` errors if `eta` < the w step), `peak_wmin = 0.05` meV added so `Epeak` tracks the soft mode instead of the 0.003–0.009 meV hyperfine line, and the field sweep now uses `invzt_solve_auto` so it spans B<sub>c</sub> instead of stopping at the old PM-only mask.

---

## The four execution amendments

**1. Whole-cc ordered medium (Task 3, user-approved 2026-07-20).** The plan's draft dominant/rest split (`E < Esplit`) was implemented, tested, and measured to carry only **47.6%** of the cc weight at 0.1 K/3 T ordered points. The projected-parity gate (5e-3 tolerance) was structurally unsatisfiable at the resulting `dSigma0 = 8.6e-3`, `dalpha_m = 5.8e-2` — while `chi0cc0` stayed bit-identical and `dm0 = 0`, proving the divergence was a pure medium-convention mismatch, not a physics bug. Root cause: the projected ordered reference itself renormalizes the FULL electronuclear `G0`, no split — the dominant/rest split is a PM/E1-only prescription. Fix: replace the split with `c0 = invz_chi0z(si, T, iwn, elastic=true)` renormalized whole, and make `Esplit`/`chi_rest` on this solver error `invzt:orderedSplitKnobs` (fail loud, not silently ignored).

**2. Stability-based phase dispatcher (already in the plan, from the pre-execution review — P0-1).** Not discovered mid-execution but worth restating because it shapes everything downstream: the round-1 draft's "ordered-first" pattern was CONFIRMED broken before Task 5 was implemented — the bare mean-field ordered moment (`m0`) was reproduced bit-identically as 1.6553 / 1.5109 / 1.1717 at 4.65 / 4.70 / 4.80 T, i.e. it persists **well past** the true QPT (crit sign change at 4.65–4.70 T). An ordered-first rule would therefore misassign the whole `[Bc, ~5.0 T]` window as FM. Task 5 was redesigned around a PM-first stability rule before any code was written: the PM leg's `crit > 0` (the tensor's own QPT criterion) decides first; the ordered leg is consulted only when the PM sample is invalid.

**3. Basis falsification → fixed-rank field-adapted basis (Task 7B, user-approved 2026-07-20).** The plan's round-2-recommended fixed **zero-field** `e2xI8` content basis was built and measured: `chi_share` = 0.00000–0.00013 across 3.0–4.65 T (essentially zero — the zero-field Ising doublet has ⟨+|J<sub>z</sub>|−⟩ = 0, and the true low-energy manifold rotates far from the zero-field content under the transverse ordering field). The moving `E < Esplit` energy cut was independently re-measured on this basis question too: ndom drifts 13→8 across the same field range. **Scorecard: Esplit-cut 48% / zero-field-content 0% / fixed-rank-adapted 97.7%.** The shipped basis fixes the RANK (16, constant by construction) and adapts the CONTENT (rebuilt at every field from the ordered mean-field Hamiltonian's 16 lowest eigenstates) — resolving the round-2 P0-3 count-drift objection without reintroducing the zero-field failure. A same-day refinement (Codex T7B feedback) reframed the object as a fixed-**rank** spectral subspace and added `var_share` (0.665–0.838) as an explicit caveat: high static `chi_share` (~0.98) does **not** certify vertex convergence.

**4. Complete hybrid map (Task 7D, user-required rework).** The as-built a3d fixed point iterated V on the **isolated** 16-state system and assembled the hybrid `chi_til` only afterward. This was REJECTED: the returned `Vcc` was not generated by the returned `Kmat` (an invariant violation), and the omitted 119-state medium feedback modifies K at the *same* vertex order — not automatically negligible, since `var_share` at 3 T is only 0.665. Required rework: iterate V on the COMPLETE hybrid map (dress → `chi_H = chi_full + (chi_dom_til − chi_dom)` → EMT → contract), via a `chi_base` seam added to `invzt_sigma_tensor` (default zero = bit-identical to the legacy dominant-only path, verified). **The vindication:** self-generation residual dropped from RED 1.503 to **0.000e+00** on the complete map, and the isolated-vs-complete comparison measured **dVcc(0) = 28.5%** — a large, same-order effect, confirming the rejection was correct rather than merely cautious (dcrit shifted −0.014 in the same comparison).

---

## Scope decisions (with dates)

- **2026-07-19 (plan header, user):** transverse field only for the ordered branch — `invzt:orderedLongitudinal` on any `|Bz| > bz_tol`; no `forced_moment` port.
- **2026-07-19 (plan header, user, round 2 — superseding a round-1 deferral):** ordered `a3` ships as **`a3d`**, the full-response dominant-vertex hybrid, not the full 136-state calculation (permanently budget-refused).
- **2026-07-19 (plan header):** documentation must say "full-response, fixed-rank field-adapted dominant-vertex ordered a3d" — never "full tensor ordered a3" (round-2 minor, so it can never be mistaken for the refused calculation).
- **2026-07-19 (plan header):** the ordered Jensen moment-form bridge is an approximation diagnostic for a3d (static/n=0 gate only; finite-frequency and βΔ=3 rows are REPORT), never an exact correctness gate — the dense-vertex oracle rows are the load-bearing correctness check.
- **2026-07-20 (Task 3 amendment, user-approved):** whole-cc ordered medium (amendment 1 above).
- **2026-07-20 (Task 7B amendment + refinement, user-approved / Codex-verified):** fixed-rank field-adapted vertex basis (amendment 3 above).
- **2026-07-20 (Task 7D rework, user-required):** complete hybrid map (amendment 4 above).

---

## Operational notes

- **The `pgrep -x MATLAB` lesson.** On this machine MATLAB's process name is **`MATLAB_maca64`**, not `MATLAB` — `pgrep -x MATLAB` **never matches**, so every "is the run still alive" check built on it was silently wrong. The correct liveness check is `pgrep -f MATLAB_maca64`. This was discovered mid-plan when a 0.1 K convergence-anchor build (`conv2`) that looked dead was actually alive the whole time (a silent `Nv = 24` build); the controller briefly duplicated it via a second `nohup` launch, and the duplicate had to be killed. The separate "subagent long MATLAB calls get auto-backgrounded and die silently" problem (below) is real and distinct — some of those apparently-dead runs may in fact have completed unobserved once the correct `pgrep` pattern is used to re-check.
- **Auto-backgrounded runs.** Subagent-issued long MATLAB calls get auto-backgrounded by the harness and can die silently mid-run (observed window: 50 min–2 h), even when the controller is nominally monitoring them. From Task 4 onward the controller ran verification suites directly (foreground or explicitly monitored background) rather than trusting an implementer subagent's own long-running MATLAB call; implementers were restricted to edit+commit once a task's fast tests were green.
- **Production a3d build: user-launched `nohup` ONLY.** The rank-16 production point (0.1 K, 3 T, 16³ halfopen/dpRng-30) measured **7.43 h** — inside the ≤12 h benchmark gate, but far longer than any interactive or CI-supervised MATLAB session survives in this environment. It is gated behind `INVZ_SLOW` (test `test_a3d_production_point` in `invz_tensor/tests/test_invzt_solve_point_ordered.m`) and MUST be launched independently of any agent session, e.g.:

  ```bash
  cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion"
  nohup env INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch \
    "cd('invz_tensor/tests'); r = runtests('test_invzt_solve_point_ordered'); disp(table(r))" \
    > invzt_a3d_production_$(date +%Y%m%d_%H%M).log 2>&1 &
  ```

  Check liveness with `pgrep -f MATLAB_maca64` (never `pgrep -x MATLAB`, see above), and tail the log file for the `a3d PRODUCTION @ 0.1 K/3 T: ...` line the test prints on success.

---

## Suite status

**CORE Incomplete allowlist (5 names)** — the printed `{r([r.Incomplete]).Name}` from `runtests('invz_tensor/tests')` must equal exactly:

1. `test_ladder_production_slow` (A4 production ladder, pre-existing baseline)
2. `test_tcut_matches_field_cut_slow` (T-cut round-trip, pre-existing baseline)
3. `test_tcut_odd_on_slow` (T-cut odd-on, pre-existing baseline)
4. `test_tcut_adaptive_slide_slow` (T-cut adaptive-slide, pre-existing baseline)
5. `test_a3d_production_point` (**new** — the a3d production stub added by this plan)

All five use `assumeTrue(tc, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only')` — the repo's slow-gate convention (a skip that *registers* as Incomplete rather than silently passing or erroring). Tensor CORE progressed 63/0/4 → 67/0/4 → 75/0/4 → 78/0/4 → 81/0/4 → 89/0/4 → 96/0/4 across Tasks 3 through 7D (pre-fix); the fix commit that added `test_a3d_production_point`'s `assumeTrue` gating is what brings the Incomplete allowlist to 5. PROJECTED (`invz_projected/tests`) stayed at 143/0/19 throughout — never disturbed.

**Slow-suite deferral.** Running the CORE suite with `INVZ_SLOW=1` set will trigger all five gated tests, including the 7.4 h a3d production point above — this is *not* something to run inside an interactive or CI-style session. Defer it to the user-launched `nohup` command in "Operational notes" above; the four T-cut/A4 slow tests are comparatively cheap and can be run under `INVZ_SLOW=1` in a normal session if only those are wanted (they will simply also re-run the a3d point if left in the same `runtests` call, so isolate it if that's not desired).

---

## Provenance

- **Plan of record:** `docs/superpowers/plans/2026-07-19-invzt-ordered-side.md` (two Codex review rounds, both disposition tables at the end of the file — every finding CONFIRMED on re-verification before adoption).
- **Commit trail:** `c0c723e` (T1) → `7834593` (T2) → `ef54a52` (T3) → `456c924` (T4) → `e584555` (T5) → `d2c4bae` (T6) → `500033a` (7A) → `67054f1` + `b5305ee` (7B + refinement) → `e5e3e22` (7C) → `a25a6cf` (7D) → `886d105` (7D rework) → `1a94f51` (7D production-stub fix) → this doc + README + ODD-LOG (Task 8).
- **Physics log:** `docs/ODD-LOG.md` §A-ordered (headline numbers, full provenance).
- **Module reference:** `invz_tensor/README.html` (mode ladder, module map, LOCKED conventions, headline results, open items — all updated for a1-ordered and a3d).
- **Sibling handoff note:** `invz_tensor/SESSION-2026-07-18-invz-tensor-full.md` (the PM-side A0–A4 build this plan extends).
