# invZ expansion — development record

**Consolidated 2026-07-27.** This file replaces the execution ledger that lived at
`.superpowers/sdd/progress.md` (git-ignored, therefore not recoverable from history) and indexes the
implementation plans under `docs/superpowers/plans/`. It exists so those directories can be removed
without losing the record of what was built, why, and what was decided along the way.

Two things this file is **not**: it is not the physics (see the README pages), and it is not the
design rationale (see `docs/INVZ-DESIGN-RATIONALE.md`).

---

## 1. Where knowledge lives after the superpowers directories go away

| topic | authoritative file | status |
|---|---|---|
| Physics, data flow, function reference, limitations — projected branch | `invz_projected/README.html` | user-maintained, current |
| Same, tensor branch | `invz_tensor/README.html` | user-maintained, current |
| Jensen 1/z framework (theory reference) | `jensen_1z_framework.html` | reference, do not edit |
| ODD implementation | `odd_implementation_plan.html` | reference, do not edit |
| Three-level extension | `three_level_1z_extension.html` | reference, do not edit |
| Design decisions and rejected alternatives | `docs/INVZ-DESIGN-RATIONALE.md` | consolidated from the specs |
| **What was built, when, and what it measured** | **this file** | consolidated from the ledger |
| Ewald: derivation / design / preregistration / integration map | `docs/invzp_ewald_*.md` | frozen |
| Strict static medium: preregistration | `docs/invzp_strict_medium_prereg.md` | **frozen — never edit in place** |
| Strict static medium: Gate-0 verdict | `docs/invzp_strict_medium_gate0_report.md` | final |
| Strict static medium: failure analysis | `invzp_convg_fix_attmpt_claude.md` | final |
| Ordered residual contract | `docs/invz_ordered_residual_contract.md` | binding |
| Stage-2c discriminator matrix | `docs/invzp_task2_prereg.md`, `docs/invzp_task2_report.md` | frozen + final |
| BZ-quadrature/Γ audit | `docs/invzp_phase1_quadrature_prereg.md`, `docs/invzp_phase1_report.md` | frozen + final |
| ODD numerical log | `docs/ODD-LOG.md` | append-only |
| Ordered 1/z consolidated state | `docs/invzp_ordered_1z_state.md` | superseded in part by this file |
| Per-session narratives | `docs/SESSION-*.md` | historical |

| Archived one-off diagnostic scripts | `docs/diagnostics/` (+ its README) | unmaintained, off-path |
| Unexecuted Ewald Steps 5–7 plan | `docs/invzp_ewald_step5-7_plan.md` | **live forward work** |

### Recovering the retired plans and specs

`docs/superpowers/` was removed on 2026-07-27. **Four of its files had never been committed** — the
Stage-2 ordered-thermodynamics plan, the Ewald Step-5 plan, the strict-static-medium plan (4983
lines), and the approved strict-medium **spec** (869 lines) — so they were committed first,
specifically so removal would not destroy them. All 22 plans and specs are therefore recoverable:

```
git log --diff-filter=D --name-only -- 'docs/superpowers/*'   # find the removal commit
git show <commit>^:docs/superpowers/specs/<file>              # read one
git checkout <commit>^ -- docs/superpowers                    # restore the tree
```

The unexecuted Ewald Step-5/6/7 plan was **not** deleted — it describes forward work, so it was moved
to `docs/invzp_ewald_step5-7_plan.md` and stays in the working tree.

`.superpowers/` (git-ignored scratch, 146 MB) was deleted outright and is **not** in history. Its
durable content is this file plus `docs/diagnostics/`; see that directory's README for what was
dropped and why each dropped item is regenerable.

---

## 2. Phase inventory

Branch `invzp-stage2c-diagnostic`, 420 commits total, 98 ahead of `main`. Nothing merged or pushed.
Depth below is proportional to what is independently verifiable; early phases are summarised from
their commit ranges, recent ones from measurements reproduced by the controller.

### 2025-11 — origin
RPA susceptibility, first effective-medium implementation, Anderson-mixing circular-buffer fixes.
Pre-dates the structured plan era.

### 2026-07-06 → 07-08 — projected branch foundations
`2026-07-06-invz-lihof4-susceptibility.md` (11 tasks) and
`2026-07-08-phase-diagram-two-regime.md` (3 tasks). Established the projected 1/z solver, the
two-regime phase-diagram driver, and the spectra energy-unit convention.

### 2026-07-14 — demagnetisation knob, q-path, spectra
`2026-07-14-demag-knob-qpath-spectra.md` (10 tasks). Produced the four-bullet demagnetisation
semantics still in force.

**Durable finding.** The q-path dispersion is the **uniform FM projection** `v'*Jcc*v` (Jensen
`P.Juni`), *not* `max(eig)`. Using `max(eig)` produced a spurious mirror-about-`h=1.5` artifact.

### 2026-07-16 → 07-17 — field angle, in-plane rotation, ODD, tensor branch
Four plans: `invz-field-angle` (5 tasks), `invz-inplane-rotation`, `invz-odd-mainbody`,
`invz-tensor-odd` → `invz-tensor-full` (15 tasks). Built the `invz_tensor` branch: q-grid, 12×12
tensor `J(q)`, RPA susceptibility, `chi0` split, A1/A2/A3/A4 solver modes, the four-point vertex and
its independent Python oracle, and the reproducibility pinning.

**Pinned convention:** `cfRot r ⇔ phi_ab = +r` (same sign); production angle `phi_ab = −11°`.

**Notable:** the A4 ladder shipped with an `eye(3)` climb-breaker fix that was *proven* at `e6` and
crash-reproduced on revert; the vertex oracle (`verify_tensor_vertex.py`, 138/138 checks) is an
independent reimplementation, not a replay of stored answers.

### 2026-07-18 → 07-19 — critical-T finder, run drivers, ordered side
`invzt-critical-t` (T-cut finder, two-regime driver, spec rev 5 recording execution findings E1/E2 —
branch-tracked rolling seed, cross-finder path hysteresis in the shallow-`crit` odd-on corner),
`invzt-run-drivers` (5 tasks), `invzt-ordered-side` (a1 point solver, moment-form real-axis
continuation, stability-based phase dispatcher, a3d full-response hybrid).

**Decision of record here, and it recurs:** the phase dispatcher is **stability-based** — the PM
`crit` decides first, and a bare-MF moment never assigns the phase.

### 2026-07-21 — QCP stage-1 split
`invzp-qcp-stage1-split-overlays.md`. The 1/z leg got its own stability-gated phase in
`invz_spectra_map`; the automatic overlay was honestly renamed `Bc_auto` and restricted to valid-PM
anchors; phase-0 and suspect columns are masked from both overlays before peak extraction.
Diagnostics `S.phase_1z`, `S.crit_pm`, `S.suspect` added. Diagnosis: `invzp_QCP_diagnosis.md`.

### 2026-07-22 — Stage 2, ordered thermodynamics
`2026-07-22-invzp-stage2-ordered-thermodynamics.md`, Tasks 0–5 plus 6/6b.

Built: the ordered elastic static single-site function (J 2.28–2.29) with sign-anchor identity
gates; the static-sector EMT closure for the ordered elastic propagator; applied-field/`H_MF`
self-consistency and the spontaneous root (J 2.31–2.33) on a geometric grid with a
grid-convergence gate; the `'jensen'` ordered mode; and the free-energy consistency route (J 2.34).

**Result — closure achieved:** `Bc_1z = 3.025 T` versus the PM `crit`-zero at `3.033 T`,
`|ΔB| = 0.008 T` (~0.3 %). `D_ord` monotone, `m` ratio `0.28`, pole shift `0.0025 meV`, reversal
exact.

Two things to carry forward:
- The re-densify amendment in Task 3 cut the near-boundary root error from **11.16 % → 0.34 %**.
- `invz_deltaF_ordered.m` is by its own docstring a **partial hybrid diagnostic**, explicitly *not*
  `dF(m=0)`, meaningful only at a common cutoff, with a documented ~13.7 % same-retained-order
  static-elastic residual. It is validation-only and must never be promoted into a "two routes
  agree" gate.
- A `max_outer` raise (PM 60→200, ordered 80→200) was the user's own uncommitted worktree edit,
  caught bundled into a task commit by review and re-committed with explicit attribution
  (`a1f69f2`). **Check for user edits before staging.**

### 2026-07-23 — Stage 2c, the discriminator matrix
The masking symptom below `Bc_1z` was now the blocker. A pre-registered 40-cell causal-discriminator
matrix (`docs/invzp_task2_prereg.md` → `docs/invzp_task2_report.md`, 992 lines) tested lattice
against solver.

**Raw counts (integrity check, not the verdict):** every swept cell `hmf_status = node_failed`,
`n_nodes = 34` everywhere. Stable accepted nodes exist at all four fields (6–10 of 34).
`isolated_cold == isolated_seed2 == swept` identically at every field → seed and continuation
independence. Downsampling ds2/ds4/ds8 ≈ unchanged.

**The decisive cells:**
- **G4** (real cardinality, *synthetic* distribution shape) converges **much better** than real
  (2.85 T: 33 accepted / 8 stable / 25 unstable, versus real 9/6/3).
- **G5** (real cardinality **and** real distribution shape) **reproduces the real lattice almost
  exactly** (2.85 T: 9/6/3 vs real 9/6/3; 1.1732 T: 10/10/0 vs real 10/10/0).

So the failure tracks the *distribution shape* of the coupling multiset, not the mesh. Verdict:
lattice/mesh-unresolved. A follow-on BZ-quadrature/Γ audit (`docs/invzp_phase1_*.md`) found the
legacy BZ quadrature had duplicate faces — 15³ distinct of 16³, 17.6 % redundant — and pointed at
the conditionally-convergent brute-force dipolar sum in `MF_dipole`.

Also delivered: the pre-declared ordered residual/scaling contract with a non-mutating complete-
residual checker (`docs/invz_ordered_residual_contract.md`), and the shared pure seed-safe nested
ordered node solver with checker-gated acceptance.

**Bug worth remembering:** `struct('cell_filter', {})` builds a **0×0 empty struct array**, which
silently broke every downstream `getf` on the matrix driver's default path — a path that had no
test. Fixed by building `run_opts` field-by-field, plus a regression test.

### 2026-07-24 — Ewald dipolar primitive and opt-in integration
Two plans: `2026-07-24-ewald-dipolar-primitive.md` (Step 4, the additive primitive plus Gates A/B/C)
and `2026-07-24-ewald-step5-integration.md` (Tasks 1–10, production integration).

Step 4 built the primitive with a frozen preregistration and gate ladder: screened-Hessian FD,
small-shell signs/self-term/boundary/counts, gauge covariance, periodic spectrum, structural
identities and cell-phase covariance, α/cutoff-ladder independence, an independent scalar-Coulomb
oracle, and primitive-level projector/remainder/uniform-mode checks.

Step 5 carried the opt-in backend to all six surfaces (`jq_modes`, `bz_couplings`, `jq_path`,
`spectra_map`, `spectra_qpath`, `run_spectra`) with a v5 backend-separated cache, exact-`isequaln`
cache validation, symmetric provenance forwarding and conflict rejection, and Γ-metadata that is
double-Lorentz-free. **Final suite 394 passed / 0 failed / 23 incomplete. Whole-branch review:
0 Critical / 0 Important. Production default stays `bruteforce`.**

**Two test-integrity failures caught here, both instructive:**
1. The Gate-C7 implementer initially used `dpRng` (a real `MF_dipole` cutoff, `(2·dpRng+1)³` mesh)
   as a large per-test cache tag. Result: multi-GB `geomD`, `"Variable 'info' was not saved"`
   warnings, and therefore **silently corrupted cache saves** — every "warm" read was a miss
   masquerading as a hit, defeating the gate's entire purpose. Self-caught and fixed.
2. Review then found that **no assertion proved a genuine cache read at all**: every cold→warm
   `isequaln` passes identically whether the warm call hit or silently recomputed, because a
   recompute overwrites the same deterministic filename. Fixed by a payload-poison test — corrupt
   the on-disk value, keep the cache metadata valid, and assert the warm call *returns the poison*.

A second review finding: the schema-corruption test hard-deleted all pre-existing caches matching a
prefix before snapshotting, which would destroy a future production cache. Replaced with
stash-and-restore under `onCleanup`. **Rule: never delete unrelated caches.**

**`local_qhash` requires the JVM** on the default path; `-nojvm` is unsupported. Documented rather
than worked around, because this is a plotting-oriented codebase.

### 2026-07-25 — the pivot: Ewald is not the masking fix
Three probes, run before committing to any further lattice work:

| probe | result |
|---|---|
| coupling comparison, 16³ | `Jcc0` reldiff `4.3e-4`; sorted-multiset `‖ew−bf‖₂/‖bf‖₂ = 1.2e-3` — essentially the same input |
| `invz_hmf_ordered` at 1 T and 3 T, **both backends** | identical failure: `hstar=NaN`, `node_failed`, `slope0=NaN`, `D_uni<0` at 30–79 % of nodes, `Dq_neg_count` → the full ~16384 |
| grid sweep, Ewald, 1 T | `node_failed` at `N = 8, 12, 16, 24` — grid-independent |

**The masking is solver-side, in the ordered static medium.** Ewald repairs a real and separate
offset-sensitivity defect in the couplings; it does not un-mask. Ewald Steps 6–7 (physics anchors,
default flip) were never executed.

### 2026-07-25 → 07-27 — strict-order static medium
`2026-07-25-invzp-strict-static-medium.md`, 19 tasks / 128 steps. 37 commits, 53 files,
+6744 / −320, 15 new production files, all additive, default unchanged and bit-identical.

Hypothesis: truncate the `ω_n=0` medium to strict `O(1/z)` — `K0 = Jbar − mu2·Gref`, one-shot, no
denominator feedback — removing the resummed q-denominator pole that dies ~7 % below the uniform
instability.

**Gate 0: FAIL.** Clauses (a), (c) and (e) of the frozen five-clause predicate fired. The full
analysis is `invzp_convg_fix_attmpt_claude.md`; the verdict document is
`docs/invzp_strict_medium_gate0_report.md`. In one line each:

- **(a)** the reference denominator leaves its domain at **low** field — `B ∈ {0.05, 0.25, 0.5, 1} T`
  fail; `B ∈ {2, 2.5, 2.9, 3.0} T` all return `ok`. The failure is deep in the ordered phase, not
  near the boundary, which inverts the expectation.
- **(c)** `max(omit_max)` exceeds the frozen `0.10` on **all 24 rows**, including every row that
  solves cleanly (smallest `0.17782`). The synthetic fixture had said `0.0350`.
- **(e)** the same 12 rows, and independently **both PM controls**, which converge, report
  `medium_status='ok'`, and return a **negative** mass.

The run stopped at diagnosis, as its preregistration required. No tolerance widened, no
regularisation, no scheme switch, no default flip.

Six gate tests remain **red by design** — three gate *definitions* (G1, G2, G13) proved defective,
and every relevant constant is frozen, so correcting them needs a dated preregistration amendment
(`Task-17_prereg_amendment_DRAFT.md`, unsigned). Suite at HEAD: `P=576 F=6 I=29 T=611`.

---

## 3. Standing rules and conventions

Accumulated across the project; all still in force.

**Physics / notation**
- `G = −chi`, units meV⁻¹. Ferromagnetic couplings are **positive** `J`. State this verbatim in new
  physics docstrings.
- q-path dispersion = uniform FM projection `v'*Jcc*v`, never `max(eig)`.
- Population moments: `mu_n = mean((J−Jbar).^n)`. MATLAB's `var()` uses `N−1` and is **wrong** here.
- Phase labelling is stability-based: PM `crit` decides first; a bare-MF moment never assigns a phase.
- Cite Jensen J-numbers (J 2.31, J 2.34, …), not volatile HTML equation numbers.

**Numerical discipline**
- Never broaden or regularise a static response, add a pole regulariser, flip a sign, or widen a
  tolerance as a convergence patch. If a tolerance appears wrong, escalate — do not tune it.
- A frozen preregistration is amended by appending a **new dated section**, never by editing in place.
- Do not select a threshold after seeing the value it will be applied to.
- Report ratios in consistent units. (An absolute displacement over a relative residual is not a
  ratio; that mistake shipped once and was caught in review.)
- "Converged" is not "physical". Compute stability diagnostics for every node; gate acceptance
  separately.

**Testing**
- `INVZ_SLOW=1` gates the only tests that reach the adaptive-grid paths (extension, redensification,
  near-boundary root). A green default suite does **not** certify those paths. Run the blast-radius
  trio under `INVZ_SLOW=1` after any refactor of `invz_hmf_ordered`, `invz_ordered_node_solve`, or
  `invz_ordered_residual`.
- Run MATLAB in the **foreground** with the result inline. A backgrounded-and-polled suite run lost a
  whole task's report once and left a log truncated mid-run.
- MATLAB is not on `PATH`: use `/Applications/MATLAB_R2025a.app/bin/matlab`. Scripts run via
  `run(...)` execute with the script's own directory context, so `addpath` inside them needs
  absolute paths.
- Never delete unrelated cache files in a test. Stash and restore under `onCleanup`.
- A cold→warm cache assertion proves nothing unless a poisoned payload makes silent recompute fail.

**Repository hygiene**
- The working tree carries ~30 unrelated user edits and 11 stashes. **Never** `git stash`,
  `git restore`, `git checkout -- <path>`, or `git clean`. Stage exact paths only, never `git add -A`.
- No authorship trailer on this branch. Verify commit bodies with `%B`.
- Error identifiers are `invz:*` / `invzt:*` only.

**MATLAB traps that have actually bitten**
| trap | consequence |
|---|---|
| `all([])` is `true` | vacuous guards pass |
| `isfield([], 'x')` is `false` | absent-vs-empty confusion |
| `max` ignores `NaN` | a missing value reads as "not the maximum" |
| `nargout` is `−1` for anonymous handles | output-count guards cannot fire; `@() error(...)` throws `MATLAB:maxlhs` |
| `struct('f', {cell})` | yields a **scalar** struct, not an array |
| `struct('f', {})` | yields a **0×0 empty struct array** — every downstream field read breaks |
| subfunctions are not externally callable | a testable helper needs its own file |
| `for x = col(:)` | iterates once over the whole column, not element-wise |

---

## 4. Open items

1. **Sign or reject the strict-medium preregistration amendment.** Until then the branch cannot be
   green. Draft: `Task-17_prereg_amendment_DRAFT.md`.
2. **Merge or park the strict-medium branch.** Gate 0 failed, but the primitives are independently
   useful (domain scanner, coverage accounting, two-tier verdict, omitted-order diagnostics) and the
   default path is untouched.
3. **Review backlog, ~12 items**, mostly cases where the plan mandated something the review rubric
   treats as a defect: tests that assert nothing; a test whose name promises legacy-drift detection
   it cannot deliver because it compares two runs of the *new* code; a third verbatim copy of the
   backaction block; an unvalidated `hmax_abs`; a bare `opts.J0eff` read.
4. **Spec vocabulary question:** `phase_1z_reason` is declared a closed 11-value list; the
   implementation now has a twelfth, `response_failed`. Needs ratification or a rename.
5. **G5 tension:** the Task-18 brief says Gate 0 cannot pass with a missing finest integral, but
   `rep.pass` is specified as exactly `~(a||b||c||d||e)` with no sixth term. G5 is measured, fails on
   every field, and is not a Boolean. Recorded, not silently reconciled.
6. **Ewald Steps 6–7** (physics anchors, default flip) are planned but unexecuted. The plan is
   `docs/invzp_ewald_step5-7_plan.md`, kept in the working tree for exactly that reason.
7. **Deferred Ewald minors** (M-t1/t2/t4/t5/t7abc/t9-1/t9-2) plus the preregistration §5
   exponential-factor erratum, which needs a dated §12 E2 re-registration.
8. **The strict-medium failure leaves live physics questions** — see §11 of
   `invzp_convg_fix_attmpt_claude.md`. The most consequential: the strict candidate's own PM phase
   boundary was never determined, yet it reports a negative PM mass at both frozen controls; and the
   `ω_n=0` scheme jump was measured only at primitive level, so its propagation into `invz_lambdas` is
   unknown and cheap to measure.
9. **G1 and G5 share a tolerance-structure defect.** Both compare O(h²)-accurate discrete quantities
   against absolute thresholds chosen independently of the grid they run on; both fail while their
   underlying sequences converge cleanly at second order (G5's eight refinement ratios all lie in
   `[4.0118, 4.0175]`, against the exact second-order value of 4).
   Any successor preregistration must check each threshold's feasibility against the frozen grid
   *before* freezing both, or register convergence-order criteria instead of absolute floors.
