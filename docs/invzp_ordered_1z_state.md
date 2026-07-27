# Projected-spin 1/z ordered leg → Ewald lattice fix — consolidated project state (READ THIS FIRST)

**Updated 2026-07-27.** Single self-contained "what a clean session needs to know" doc for the
projected-spin (`invz_projected`) ordered-1/z / BZ-lattice work. It **maps and summarizes** the
frozen pre-registrations and committed reports; it does not supersede them. Supersedes the older
`invzp_z^-1_spectra_fix_attempt_Claude.md` (deleted) and folds in its still-load-bearing content.

## TL;DR (current state)

`invz_run_spectra` masks the 1/z ordered panel (`phase_1z=0`, `Sigma0=NaN`) at every field: the
projected "jensen" ordered leg does not converge on the **real** BZ couplings. **The symptom is still
open.** Two candidate root causes have been tested to conclusion:

1. **Lattice (tested; a real defect, but NOT the masking cause).** `MF_dipole`'s brute-force spherical
   real-space dipolar sum is a **conditionally-convergent lattice sum** whose truncated-sphere shape
   term is offset/cutoff-dependent (Phase-1 audit). This is genuine and was repaired by an Ewald
   summation, now built, gate-validated, and integrated opt-in through all six production surfaces
   (Steps 4–5 COMPLETE, suite 394/0/23, default still `bruteforce`). **It does not un-mask.** Measured
   2026-07-25: Ewald and brute-force couplings differ by `‖ew−bf‖₂/‖bf‖₂ = 1.2e-3`, and
   `invz_hmf_ordered` fails **identically on both backends** at 1 T and 3 T, at every grid size
   `N = 8, 12, 16, 24`. See §"The 2026-07-25 pivot" below.
2. **Ordered static medium (tested; FALSIFIED as specified).** The resummed static medium carries a
   q-denominator pole that dies ~7 % below the uniform instability. A strict-`O(1/z)` one-shot
   truncation `K0 = Jbar − mu2·Gref` was specified, preregistered, implemented additively, and
   **FAILED its own Gate 0** on 2026-07-27 (clauses a, c, e). Full analysis:
   `invzp_convg_diagnosis_Claude.md` §9.1–§9.4; verdict:
   `docs/invzp_strict_medium_gate0_report.md`.

**Ewald Steps 6–7 (Gate-D rerun, Gate-E physics anchors, default flip) were never executed** — they
were sequenced behind a masking fix that has not arrived. The production default remains
`bruteforce` everywhere.

## Reading map (load-bearing docs, in order)

0. `docs/INVZ-DEVELOPMENT-RECORD.md` — the consolidated execution history for the whole repo, and the
   index of where every other durable document lives. `docs/INVZ-DESIGN-RATIONALE.md` — the
   consolidated design decisions. Both replace the retired `docs/superpowers/` + `.superpowers/`
   directories.
1. **This doc** — the map + take-homes.
2. `docs/invzp_ewald_design.md` — the FROZEN Ewald design (architecture, scope, the "adds 0 at Γ"
   decision, the Step-5 gate structure).
3. `docs/invzp_ewald_derivation.md` — the math (screened Hessian, k=0/Γ) + the numerical Verification
   A–E that closed it.
4. `docs/invzp_ewald_prereg.md` — the FROZEN acceptance pre-registration (params, metrics, Gates A–F,
   the jensen/HMF target) + the §12 Errata.
5. `docs/invzp_ewald_integration_map.md` — the Step-5 integration surfaces (`Jpath_base_cc`/`Jgamma_cc`,
   the six consumers, the wrapper API).
6. The executed Step-4 plan (`docs/superpowers/plans/2026-07-24-ewald-dipolar-primitive.md`) was retired 2026-07-27; recover with `git show <removal-commit>^:<path>`. See `docs/INVZ-DEVELOPMENT-RECORD.md` §1.
7. `docs/invzp_phase1_report.md` + `docs/invzp_phase1_quadrature_prereg.md` — the audit that motivated
   Ewald (committed, frozen).
8. `docs/invzp_task2_report.md` + `docs/invzp_task2_prereg.md` — the LATTICE/MESH-UNRESOLVED verdict.
9. `invzp_convg_diagnosis_Claude.md` (root) — the consolidated ordered-leg/QCP diagnosis, including
   the independent QCP replication and surviving Stage 1/2 conclusions in §9.5.
10. Memory `invzp-jensen-realcoupling-nonconvergence` (`~/.claude/.../memory/`) — the sharp cross-session
    statement of the discovery + the physics correction.

## History (condensed; commit ranges recoverable from `git log`)

- **Stage 1 (DONE):** split the RPA/auto overlay from the 1/z self-energy-dressed curve — they have
  DIFFERENT critical fields and must not be forced to close (`S.phase_1z`/`S.crit_pm`/`S.suspect`,
  `invz_boundary_field.m`).
- **Stage 2 (DONE, 7 tasks):** full nonlinear ordered 1/z thermodynamics on the **projected** path
  (Jensen J 2.28–2.29 elastic static sector + J 2.31–2.33 H_MF self-consistency + J 2.34 free-energy
  validation) — this introduced the `'jensen'` ordered mode that later fails on real couplings. A
  ~13.7% same-order static-elastic residual in the closed 2×2 free-energy check is a KNOWN, tested
  approximation limitation (pinned by `test_invz_deltaF`), not a defect.
- **The discovery (2026-07-23):** every Stage-2 test used a 24-point SYNTHETIC `Jnu` fixture; NONE ran
  the real 16384-point `invz_bz_couplings` grid through the jensen leg. On the real grid it masks at
  every field.
- **Stage 2c (DONE, diagnose-first):** fixed two real code defects — P0-3 (failed-node seed carriers
  contaminate later seeds) and P0-4 (node acceptance used pre-mix step, not the full coupled residual)
  — and verified they are NOT the blocker. The Task-2 causal-discriminator matrix → verdict
  **LATTICE/MESH-UNRESOLVED**; §10 traced the endpoint-inclusive BZ grid's duplicated boundary faces
  (16³ → only 15³=3375 distinct periodic points) over-weighting the unweighted-mean EMT.
- **Phase 1 (DONE, committed `086d102`):** coupling-only audit — the half-open uniform grid removes the
  duplicate faces (item 1 passes) but the multiset stays **offset-sensitive (~143× tol, non-shrinking)
  + cutoff-drifting** → frozen stop rule → escalate to Ewald.
- **Ewald Step 4 (DONE 2026-07-24, `086d102..fcb3031`):** see below.
- **Ewald Step 5 (DONE 2026-07-24/25, `fcb3031..2ee310b`, 10 tasks):** opt-in backend integrated
  through all six surfaces (`jq_modes`/`bz_couplings`/`jq_path`/`spectra_map`/`spectra_qpath`/
  `run_spectra`), v5 backend-separated cache with exact-`isequaln` validation, symmetric provenance
  forwarding + conflict rejection, Γ-metadata double-Lorentz-free. Suite **394/0/23**, whole-branch
  review 0 Critical / 0 Important, **default stays `bruteforce`**.
- **The 2026-07-25 pivot — the lattice hypothesis is closed out.** With the backend finally reachable
  from production, three probes were run *before* committing to Steps 6–7:
  (i) Ewald ≈ brute-force couplings — `Jcc0` reldiff `4.3e-4`, sorted-multiset
  `‖ew−bf‖₂/‖bf‖₂ = 1.2e-3`, `J_max`/`J_min` identical to ~`1e-4`;
  (ii) `invz_hmf_ordered` at 1 T and 3 T fails **identically on both backends** — `hstar=NaN`,
  `status='node_failed'`, `slope0=NaN` (even the `h=0` predictor fails), `D_uni<0` at 30–79 % of 33
  nodes, `Dq_neg_count` → the full ~16384;
  (iii) grid sweep at 1 T: `node_failed` at `N = 8, 12, 16, 24` — grid-independent, not a
  quadrature-density artifact.
  **Conclusion: the masking is solver-side, in the ordered static medium.** Ewald fixed a separate,
  real coupling-correctness defect. It is not the masking fix.
- **Strict-order static medium (2026-07-25 → 07-27, `2ee310b..6937e9d`, 19 tasks):** the next
  hypothesis — truncate the `ω_n=0` medium to strict `O(1/z)`, removing the resummed q-denominator
  feedback. Built additively (15 new production files, default unchanged and bit-identical),
  preregistered before the first strict run, and **FAILED Gate 0**: the reference denominator leaves
  its domain at *low* field, the omitted-order ratio exceeds its frozen limit on every measured row
  including those that solve cleanly, and both PM controls converge to a *negative* mass. Six gate
  tests are red **by design** pending a dated preregistration amendment. Read the consolidated
  `invzp_convg_diagnosis_Claude.md` and implementation plan `invzp_convg_fix_Claude.md` before
  proposing a successor.

## The Ewald reformulation — decisions (authority: design/prereg)

- **Scope:** projected non-ODD path only; additive. `MF_dipole`, the tensor path, legacy RPA, and the
  ODD path untouched (later, separate adopters — all would benefit).
- **Convention (frozen):** `T=+∂∂(1/r)`; `dip=−Σ'T e^{−iq·(R+d)}` (total-displacement), Å⁻³,
  `[3,3,ntau,ntau]`. Three parts: erfc-screened real + Gaussian-screened reciprocal (`e^{+iG·d}`,
  `k=q+G`) + self; exact `k=0` omitted (`conducting_k0_omitted`/tinfoil). Half-open q-reduction + gauge
  restore.
- **Γ crux — resolved:** the Ewald branch **adds 0 at Γ** — the isotropic Lorentz is already in
  `dip_reg = dip_sphere − (4π/Vc)(1/3)δ` — so the four-bullet demag semantics are preserved.
- **Acceptance:** primitive convergence is not enough — the frozen Phase-1 BZ grid/offset gates are
  rerun on Ewald (Gate D), and physics anchors **including a jensen/HMF target** (Gate E) must pass
  before the default flips.

## Step-4 build + gate roll-up (all PASS; `086d102..fcb3031`, 7 files, all additions, suite 269/0/23)

Files: `invz_dipole_ewald.m` (root) + `invz_projected/tests/{invz_ewald_fixtures, invz_ewald_metrics,
invz_scalar_ewald_ref, test_invz_dipole_ewald, test_invz_dipole_ewald_ref, test_invz_dipole_ewald_gammaC}.m`.

| Gate | Worst frozen metric |
|---|---|
| A1 screened-Hessian FD (+ prod-`g_ab` bridge) | analytic 6.9e-9 / bridge `M_id` −1.8e-15 |
| A2 small-shell real signs | `M_T` −2.4e-9 |
| A3 recip-gauge covariance + periodic spectrum | `M_id` −4.1e-14 / coupling −1.0e-10 |
| A4 half-open reduction + completeness | `M_id` −5.0e-14; all `G` present |
| A5 structural identities incl. `U'·Dcell·U` | `M_id` −4.1e-14 |
| A6 self term | `M_T` −3.1e-11 |
| A7 boundary guard + no 4th term | `M_id` −4.1e-14 |
| A8 counts vs enumeration + 3 pre-alloc caps | 227=227; caps fire correctly |
| A9 α + separate-axis cutoff independence + `M_J` | default-vs-joint `M_T` 6.5e-5 / `M_J` 2.6e-5 |
| B independent scalar-Coulomb oracle (separate path) | off-diag `M_FD` −2.0e-8 (never = `M_T`) |
| C1/C2 filter omits only `G=0` at Γ; `dip_reg(0)` | `M_id` −4.1e-14 |
| C3 isolated `(4π/Vc)q̂q̂ᵀe^{−|q|²/4α²}` projector | `M_id` −2.0e-14 |
| C4 even/odd remainder | `E_even∝s²`, `E_odd∝s` |
| C5 uniform-mode support (isolated `G=0`) | `v'Δv` margin 0.0 |

Built with a fresh implementer plus an independent review per task and a final whole-branch review —
all clean. The detailed per-task ledger used to live in `.superpowers/sdd/` (git-ignored, now
retired); its durable content is consolidated into `docs/INVZ-DEVELOPMENT-RECORD.md`.

**Known erratum (E1):** prereg §5 Gate-C check-5 prose omits the `exp(−|q|²/4α²)` factor check 3
requires; the **code is correct** (see prereg §12 Errata). Relevant when Step 5 wires caller-level
Gate-C.

## Roadmap — revised 2026-07-27

**Ewald track (paused, not abandoned):**

1. ~~**Step 5 — Integration**~~ **DONE** (`fcb3031..2ee310b`). Default stays `bruteforce`.
2. **Step 6 — Gate-D** (NOT STARTED): rerun Phase-1 BZ grid/offset gates on Ewald, both Γ policies.
3. **Step 7 — Gate-E/F + adopt** (NOT STARTED): physics anchors; select Γ (P-complete vs P-drop);
   flip default; recompute `Bc_PM`.

Steps 6–7 were sequenced behind "the masking symptom is fixed", and it is not. Gate E as written
includes a **jensen/HMF target** that no current scheme can meet, so Gate E cannot pass until the
ordered leg solves. Either the default flip is decoupled from the masking fix — Ewald's coupling
correctness stands on its own evidence and does not need the ordered leg to work — or Steps 6–7 wait.
**That is an open decision, not an oversight.** The unexecuted plan is
`docs/invzp_ewald_step5-7_plan.md`.

**Ordered-solver track (open):**

The strict-`O(1/z)` static medium was the tested successor hypothesis and it failed Gate 0. Per the
frozen preregistration §3, any of: carrying another moment, changing `Gref`, or truncating other
Matsubara sectors is a **new theory candidate requiring its own spec and fresh preregistration** —
never an in-run fallback. The consolidated diagnosis §9.4/§10.2 records what a successor must answer:
the PM boundary has now been measured, leaving the reference construction at large moment and whether
the `ω_n=0` scheme jump propagates past the primitive level into `invz_lambdas` as sharper open
questions.

The original ordered exact-`h` lattice audit ("Phase 3") and the underived-unstable-branch paths
(3A/3B) remain hard-stopped.

## Durable gotchas (carry these forward)

- **Convergence ≠ a physical state.** Since `Gstat=−chi<0` and `max(Jnu)<J0eff`,
  `Dq = D_uni + (J0eff−J(q))·chi ≥ D_uni`, so any `Dq≤0` ⇒ `D_uni≤0` — the trial state is on the
  UNSTABLE side of the uniform static mode (expected: H_MF integrates outward from the unstable h=0 PM
  state below Bc). A Picard/EMT loop reporting `converged` proves nothing physical without checking
  `D_uni`/`min(Dq)`.
- **Never regularize a pole / widen a tolerance / weaken a gate to make a solve converge or a plot fill
  in.** Broadening a static/real response is a THEORY decision, not a numerical patch.
- **Test against real coupling DENSITY, not a synthetic proxy.** This trap has now fired **twice**.
  (i) The 24-point fixture masked a production failure through 7 Stage-2 tasks and full green suites.
  (ii) The same fixture then reported the strict medium's omitted-order ratio as `0.0350`, comfortably
  inside its `0.10` limit; on the real 16384-point multiset the same quantity is `1.20`–`16.4`, a
  factor of 34–470 worse. Any diagnostic whose value grows with `|Gref|` **must** be validated on the
  real density before it is treated as evidence of control. A synthetic fixture can carry the right
  cardinality and still miss the distribution's tail — indeed the Stage-2c discriminator matrix showed
  exactly that: matching real *cardinality* alone converged far better than reality, while matching
  cardinality **and** distribution shape reproduced reality almost exactly.
- **Cite Jensen's own J-numbers**, never `jensen_1z_framework.html`'s equation numbers (its Section-9
  numbering has drifted +2 vs `HEAD`).
- **`MF_dipole` cost ~ O(N_grid³·dpRng³):** N=20,dpRng=50 ≈ 466 s for ONE (convention,offset). A
  brute-force lattice/offset/dpRng sweep is multi-HOUR — evidence for Ewald over more brute-force rungs.
- **Frozen pre-registrations stay frozen:** any tolerance/target change needs a NEW dated
  re-registration (see prereg §12), never an in-place edit.
- **Git safety (a prior run corrupted the index):** stage by explicit named path, NEW commits only,
  never touch the user's unrelated uncommitted edits (`invz_tensor/*`, `README.html`,
  `invz_run_spectra.m`, `jensen_1z_framework.html`, `docs/SESSION-*.md`) or root-level
  `*_Codex.md`/plan/review notes.

## How to verify independently

- Read the frozen authority docs (list above).
- `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "addpath('.'); addpath('invz_common'); runtests('invz_projected/tests')"` → 269/0/23; per-gate margins print in the output.
- Inspect `invz_dipole_ewald.m` against the frozen convention (prereg §1).
