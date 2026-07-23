# Phase 1 — BZ-quadrature/Γ coupling audit: FROZEN pre-registration

**Status: FROZEN** (user-approved 2026-07-23: D1 approved, D2 approved with the amendments below,
D3 confirmed). These rules are fixed BEFORE any coupling audit result is inspected. Phase 1 is
**coupling-only and additive**: NO ordered / critical-field / Tc solves, NO change to production
behaviour or the frozen Task-2 pre-registration (`docs/invzp_task2_prereg.md`). Any departure from a
rule here STOPS and reports rather than reinterpreting.

## Scope and the field-independence fact
Phase 1 runs **checks 1–6 only** (point uniqueness; reciprocal periodicity; cardinality + Γ count;
weight normalization; coupling-multiset summaries + refinement; offset sensitivity). The physical
benchmark gate (§Phase-2) requires model solves and is therefore **Phase 2, not a Phase-1
deliverable**.

**The projected BZ couplings are field-independent.** `invz_bz_couplings(ion, opts)` takes no field
argument (verified: `invz_bz_couplings.m:1,14,16`); the couplings depend only on
`(ion, grid, offset, dpRng, grid-convention, Γ policy)`. Phase 1 therefore has **no field axis at
all** — statistics are NOT evaluated "at a field." (The four physical fields are downstream of
`Bc_PM`, which is Phase 2/3.)

## Frozen reference scale
`J_ref = |Jcc0|` of the full-Γ (P-complete) convention at the baseline lattice (N=16, dpRng=30):
**`J_ref = 0.006424435656 meV`** (computed once, frozen here before any audit result is viewed;
`invz_bz_couplings(invz_ion(), struct('grid',[16 16 16],'dpRng',30))`, `info.Jcc0`).

## Conventions under test
- **Grid:** half-open periodic (`qVec_generator(..., 'endpoint', false)`: `q_a = lo + (0:N-1)/N`),
  no duplicate ±0.5 face. The legacy endpoint-inclusive grid is retained ONLY as a labeled
  **legacy-parity diagnostic** (it is expected to fail item 1; that is the documented baseline, not
  a stop).
- **Offsets:** the eight `{0,½}³` phase offsets (**eight *including* the `000` baseline**), built by
  constructing the refined `2N` grid **once** and partitioning it into the eight sub-offsets.
  Shifted coordinates wrapped into one BZ (`mod(q+0.5,1)−0.5`).
- **Grid sizes:** `N ∈ {12, 16, 20}`. An optional `N=24` confirmation rung may be added ONLY if
  pre-registered before its result is inspected.
- **Cutoffs:** `dpRng ∈ {30, 40, 50}`. Optional `dpRng=60` confirmation rung, same pre-registration
  rule.
- **Γ policy — BOTH reported symmetrically (D3); neither pre-selected:**
  - **P-complete:** keep Γ; uniform weight over all kept points; the uniform critical pole is
    handled separately in the `D_uni`/`Dq` diagnostic (Phase 3), not by dropping it here.
  - **P-drop:** remove the exact-Γ row; renormalize the remaining weights consistently.
  Weights and Γ treatment are defined **globally per convention**, never inferred from row counts.
  **Under P-drop the eight independently-renormalized offset subsets no longer have exact
  refined-`2N`-grid partition equivalence** — results are reported separately per Γ policy and the
  partition identity is not used after filtering without qualification.
- **Weights:** UNIFORM `1/Npts` per kept point. (The ordered EMT is an unweighted `mean` over
  `Jnu_flat` with no weight argument — `invz_emt_scalar.m:43,48`, `invz_emt_static_ordered.m:48,50`
  — so a uniform point set needs no solver change; nonuniform weights are OUT OF SCOPE for this
  route. `invz_tensor/invzt_qgrid.m` already implements this "uniform weight over kept points"
  contract.)

## Frozen pass/fail rules

**1. Point uniqueness.** Two reduced-q points are identical iff `max|Δq_wrapped| < tol_uniq`,
**`tol_uniq = 1e-12`**. PASS: the half-open grid has exactly `N³` distinct periodic points (0
duplicates) at every N and offset. The legacy inclusive grid is expected to FAIL (report `N³`
nominal vs `(N−1)³` distinct — for N=16, 4096 vs 3375); that failure is the documented baseline.

**2. Reciprocal periodicity.** For a sample of integer reciprocal vectors `G`, the **sorted
eigenvalue-branch spectrum** (or equivalent coupling multiset) at `q+G` equals that at `q` within
`|ΔJ| ≤ AbsTol_J + RelTol_J·|J(q)|`, **`AbsTol_J = 1e-10 meV`, `RelTol_J = 1e-8`**. Test the sorted
branch spectrum / coupling multiset — **NOT** elementwise equality of raw sublattice matrices (which
may differ by a reciprocal-space gauge transformation). PASS required for every convention.

**3. Cardinality + Γ count.** Per convention/offset/N: report nominal vs distinct cardinality, the
exact-Γ row count, the individual weight rule, and the total weight. PASS: P-complete has cardinality
`N³` (half-open) with Γ present only in offsets that contain it; P-drop has cardinality
`N³ − (Γ count)` with Γ count 0 after filtering. The one-row cardinality difference between a
Γ-containing and a Γ-free offset under P-drop is expected and handled by the global weight
definition.

**4. Weight normalization.** `|Σ weights − 1| ≤ 1e-12` for every convention/offset/N.

**5. Coupling-multiset summaries + refinement (field-free).**
- **Summaries:** normalized **mean, variance, min, max**, and the quantiles
  **{0.05, 0.25, 0.5, 0.75, 0.95}** of the flat coupling multiset; plus the energy scalars
  **`J0eff`, `Jcc0`, `max(Jnu)`**.
- **Normalization (units-correct):** mean, min, max, quantiles → `s = statistic / J_ref`; variance
  → `s = variance / J_ref²`.
- **Normalized-shape gate:** `|s₂ − s₁| ≤ 1e-6 + 1e-3·max(|s₁|,|s₂|)`.
- **Energy-scalar gate** (`J0eff`, `Jcc0`, `max(Jnu)`): `|J₂ − J₁| ≤ 1e-6·J_ref + 1e-4·max(|J₁|,|J₂|)`.
- **Refinement:** report BOTH grid steps `12→16` and `16→20`, and BOTH cutoff steps `30→40` and
  `40→50`. **Gate convergence on the finest comparison only** — grid gate `N=16→20`, cutoff gate
  `dpRng=40→50` — AND additionally require the spread **not to grow** relative to the preceding step
  (final-step spread ≤ preceding-step spread). (Requiring every step to individually meet the final
  tolerance would wrongly demand the coarsest rung already be converged.)

**6. Offset sensitivity.** Report the spread of the item-5 summaries across the eight `{0,½}³`
offsets at **every** grid size; apply the pass/fail agreement gate (item-5 tolerances) at the finest
rung **`N=20, dpRng=50`**. The final offset spread must not exceed the earlier spread. If the finest
rung fails, report the failed statistic and its direction before extending the ladder.

## Phase 2 (frozen now; NOT run in Phase 1) — physical benchmark gate
A candidate convention is physically accepted only if the following anchors PASS under it, each at
its **own documented physical tolerance**:
- LiHoF4 `Σ_c` and `Jcc0` anchor — `test_invz_sigma_crit.m`;
- zero-field `Tc` case — `test_invz_critical.m`;
- 310 mK critical-field case — `test_invz_critical.m`;
- qualitative ordered-state and soft-mode case — `test_invz_ordered_phase.m`.

**Existing regression tests encode the legacy construction and remain UNCHANGED.** Add
**parameterized benchmark wrappers** that inject the candidate couplings when evaluating the new
convention — separating (i) legacy-parity regression from (ii) physical acceptance of the candidate
quadrature. The exact grid-dependent values in `invz_odd_anchors.m` are legacy/ODD diagnostics ONLY
and are **not** convention-selection gates (matching old-grid numbers cannot establish that a new
quadrature is physically correct). Because these are critical/ordered solves, they are a Phase-2
selection gate, not a Phase-1 deliverable.

## Selection / escalation / stop rules (frozen)
- **Convention selection (Phase 2):** choose the convention that (i) passes items 1–4 identically at
  every rung, (ii) passes item-5 convergence and item-6 offset agreement, and (iii) passes the
  Phase-2 physical benchmark set. If both Γ policies qualify, the production choice is the one the
  **theoretical derivation** supports (theory input required at Phase 2, not now).
- **Escalation:** if the corrected half-open uniform construction still fails item-5 `dpRng`
  convergence (couplings keep moving as the real-space cutoff grows), that indicates a
  conditionally-convergent real-space dipolar sum → escalate to an **Ewald / convergence-accelerated
  dipolar summation** (Lorentz + demagnetization separated analytically) BEFORE Phase 3.
- **Stop / hard-fail:** if NO convention passes items 1–6, Phase 1 reports "no coherent BZ quadrature
  at this construction level" and the lattice question escalates to the dipolar-sum method as the
  operative blocker — still short of any 3A/3B path.

## Execution order (frozen)
1. Freeze this pre-registration (done — this file). 2. Implement and run Phase-1 checks 1–6 only
(coupling-only, additive). 3. Review the Phase-1 evidence **without selecting Γ post hoc**. 4. Phase
2: resolve Γ theoretically and run the exact physical benchmark set through parameterized
candidate-coupling wrappers. 5. Select + freeze the quadrature, weights, and Γ policy only if BOTH
the numerical and physical gates pass. 6. Recompute `Bc_PM` and re-freeze the audit fields if the
selected convention changes the production couplings. 7. Only then resume the ordered exact-`h`
lattice audit and the subsequent `nH` refinement.

**The frozen Task-2 §D observable tolerance remains unchanged throughout.**
