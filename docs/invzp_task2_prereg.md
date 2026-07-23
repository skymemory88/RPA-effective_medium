# Task 2 — Pre-registered thresholds, stop rules, and experiment matrix (FROZEN)

**Status: FROZEN pre-registration.** Set by the user 2026-07-23 and committed BEFORE any
discriminator run. These definitions are not to be tuned after seeing the matrix. Any implementation
that must depart from a definition here STOPS and reports; it does not silently reinterpret.

**Notation.** `D_uni = 1 + (J0eff − K0)·Gstat` (uniform static pole observable). `Dq = 1 + (J(q) −
K0)·Gstat` per flat mode. "Accepted" / "checker-accepted" = passes the Task-1 checker
`invz_ordered_residual` (all four blocks A–D + finite). Baseline physical case: T = 0.1 K, real
couplings `invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30))`.

---

## A. Per-node stability classification (stability scalar)
Define the **stability scalar** `s = min(D_uni, min(Dq))`. Classification:

| class | rule |
|---|---|
| **stable** | `s > +1e-3` |
| **marginal** | `|s| ≤ 1e-3` |
| **unstable** | `s < −1e-3` |
| **unconverged** | separate classification (checker NOT accepted, or non-finite) |

**Only checker-accepted exported states are classified stable/marginal/unstable.** A node whose
solve is not checker-accepted is **unconverged**, full stop — it is NOT labelled "unstable" no matter
how negative its transient `Dq`/`D_uni` got during iteration. **Negative transient iterates are
diagnostic evidence, not proof of an unstable solution.** (Consequence: the current physical
failure — 25/34 nodes hitting max_iter — is `unconverged`, not yet `unstable`; proving genuine
instability requires an *accepted* state with `s < −1e-3`, per F.)

## B. Convergence tolerances (inherited from Task 1 — not re-opened)
The committed checker contract (`docs/invz_ordered_residual_contract.md`): blockA/C outer-Sigma
`tol_outer = 1e-8`; blockB static abs `1e-10` + rel `1e-10·|Gstat|`; blockD dynamic `K(2:end)`
identity at `tol_outer`. Task 2 uses these verbatim.

## C. Agreement across seeds and isolated-vs-swept (reproducibility of an accepted state)
Compare the **complete exported state**, not merely `Sigma0`/`D_uni`:
`Sigma(:)`, `K(:)`, `lambda(:)`, `K0`, `Gstat`, `D_uni`, and the endpoint/root quantities
(`hstar`, `m_star`, `r_star`). Agreement test per quantity is a **combined absolute + relative**
tolerance — relative tolerance ALONE is invalid near `D_uni = 0`:
`|a − b| ≤ AbsTol_q + RelTol·max(|a|,|b|)`, with **RelTol = 1e-6** and `AbsTol_q` a
**dimensionally scaled** floor per quantity: `AbsTol_q = 1e-8 · scale_q`, where
`scale_q` = 1 for the dimensionless Σ/`D_uni`, `|J0eff|` for the energy quantities `K/lambda/K0`,
a characteristic `|Gstat|` for the inverse-energy `Gstat`, and the field/moment scale for
`hstar`/`m_star`. (The implementer fixes each `scale_q` in code from the quantity's units and
documents it; no post-hoc adjustment.)
- **Cold vs multi-start** (same node, different seeds) and **isolated vs swept** (a node solved
  alone vs reached by continuation): the accepted state must agree under the above test.
  Disagreement ⇒ the accepted state is seed/continuation dependent (a defect to report).

## D. Mesh size / offset / dpRng convergence — exact tolerance + identical classification
Across the accepted lattice ladder (grid size, grid offset, `dpRng`), require **BOTH**:
1. **Identical stability classification** (the discrete class of A is the SAME at each field across
   every rung of the ladder), and
2. numeric convergence of `s` (and its components `D_uni`, `min(Dq)`) under the combined test
   `|Δ| ≤ 1e-6 + 1e-4·max(|·|,|·|)` (AbsTol 1e-6, RelTol 1e-4 on the dimensionless `s`).
- **Grid-offset semantics (dipolar Γ):** a shifted q-grid must use a deliberate, tested offset that
  preserves the dipolar-Γ convention — a naive translation can change whether a Γ-equivalent point
  receives the Lorentz term. Two offsets are run: the unshifted grid, and a half-step
  body-centred shift with explicit, tested Γ handling. q/branch provenance is retained for the
  closest modes at every rung.
- **dpRng ladder:** `dpRng ∈ {30 (baseline), 40}`.
- **Verdict on failure:** if the stability class **disagrees across the ladder** (e.g. sign flip) OR
  `s` **fails to converge under refinement**, the result is **LATTICE/MESH-UNRESOLVED** — a distinct
  hard-stop classification (see F), NOT evidence for 3A or 3B.

## E. Existence bar and frozen physical field locations
- **Field derivation (before any ordered run):** compute `Bc_PM`, the paramagnetic critical field,
  from the **independent PM-mass calculation** (the PM-leg mass `→ 0`, i.e. the field where
  `1 − J·chi0_PM` closes — implemented from the existing PM path, documented). Freeze the ordered-run
  fields as: **`{0.25, 0.55, 0.80}·Bc_PM`** plus the existing **`2.85 T` defect anchor** (four
  physical fields, transverse, real 16³, T = 0.1).
- **Existence bar (defensible ordered continuation):** at least **3 of the 4 physical fields** must
  yield a **stable (class A) checker-accepted** ordered state, each agreeing across **≥ 2 seeds and
  ≥ 2 mesh offsets** (per C/D). **Synthetic fixtures do NOT count** toward this physical
  existence bar (they are lattice-diagnostic only, per G).

## F. Trigger / stop rules (classification → path)
The committed report classifies the blocker and selects the next path. **Lattice/mesh-unresolved is
adjudicated FIRST and is a hard stop before 3A/3B:**

| verdict | trigger |
|---|---|
| **LATTICE/MESH-UNRESOLVED** (hard stop, precedes 3A/3B) | Per D: stability class disagrees across the size/offset/dpRng ladder, or `s` fails to converge under refinement. The lattice construction/normalization must be resolved before any solver or theory path is entertained. |
| **Task 3A** (solver conditioning repair) | A **checker-accepted, stable** solution exists **at the same node where the production iteration fails**, **reproducible across seeds AND lattice realizations** (isolated/cold/downsampled reaches a class-A accepted state the dense swept run does not). Failure is conditioning, not physics. |
| **Task 3B** (theory / reformulation) | A **residual-clean unstable INTERVAL** along the required HMF path — a run of accepted states with `s < −1e-3` over a contiguous field/`h` interval, **not one node and not a transient iterate** — **persists under HMF-grid refinement, seeds, and offsets.** Genuine unstable-branch traversal in the static approximation. |
| **3B then 3A** | Unstable interval is intrinsic on part of the required range while a stable pole-free branch exists elsewhere. |
| **UNSUPPORTED** | No configuration yields a defensible stable continuation (existence bar E unmet) AND no residual-clean unstable interval is established AND lattice is resolved — i.e. **unsupported by the current derivation and the completed search.** This is NOT a proof that no mathematical solution exists; it records that the current static-Jensen construction + this search found none. |

- **Seed/state defect** is treated as already largely EXCLUDED (Task 1b-ii-A: physical still masks
  with the verified P0-3 rollback). Task 2 CONFIRMS this quantitatively via C (cold-vs-continued full
  -state agreement), rather than re-litigating it.
- No path is selected because it produces fuller plots.

## G. Experiment matrix (NOT trimmed — lattice construction is a live causal possibility)
Each cell records: checker-accept?, the class of A (with `s`, `D_uni`, `min(Dq)`), the full state for
C, seed provenance, and q/branch of the closest pos/neg modes (Task-0 trace). **Serial first**;
parallel only if it does not weaken determinism or provenance.
1. **Isolated vs swept** at the four E fields (node converges alone but not in continuation?).
2. **Cold vs continued** seeds at the four E fields (confirm seed-defect exclusion via full-state C).
3. **Physical downsampling:** deterministic subsets of the 16³ q-set at {full, 1/2, 1/4, 1/8}
   (isolates coupling density from distribution).
4. **Cardinality-controlled synthetic:** duplicate the 24-pt synthetic distribution up to the real
   cardinality (isolates cardinality from distribution shape). *(Lattice-diagnostic; not an
   existence endpoint.)*
5. **Histogram-matched fixture:** synthetic couplings matched to the real distribution shape at real
   cardinality. *(Lattice-diagnostic; not an existence endpoint.)*
6. **Grid offsets** (D) and **dpRng ladder** (D) at the four E fields.

## H. Deliverable and gate
A committed report (`docs/…`) with the COMPLETE matrix (every cell, not merely the verdict), the A
class counts, the C/D agreement outcomes, and the F classification with the triggered verdict. This
is the Task-2 decision gate: **HARD STOP** — the user decides 3A / 3B / 3B-then-3A / lattice-first /
unsupported. Tasks 3A/3B/4 are NOT started here.
