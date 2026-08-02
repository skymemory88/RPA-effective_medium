# Multilevel 1/z theory and implementation evidence map

- Status: evidence index; not an implementation plan
- Active plan: [`../../invzt_multilevel_1z_next_steps.md`](../../invzt_multilevel_1z_next_steps.md)

## Evidence classes

- **Current verified fact:** checked against the present
  `invzt-multilevel-1z` tree.
- **Retained historical evidence:** recorded in an identified Git object but
  not yet reproduced on this branch.
- **Inference:** a consequence of verified or retained evidence that still
  needs a direct branch-local check.
- **Unresolved:** a required theory or implementation decision with no accepted
  answer yet.

## Current verified sources

1. `jensen_1z_framework.html`, especially sections 2--6 and 11.8, is the local
   two-level Jensen derivation and establishes the reduction target. Section
   11.8 records the additive multilevel correction and symmetric matrix
   resummation; its banner now makes clear that it is reference material.
2. `invz_tensor/invzt_sigma_tensor.m` already carries an additive matrix
   correction `Vmat`, contracts a connected four-point object with the medium,
   and uses the symmetric-bracket map. It is a structural reference, not the
   new production implementation.
3. `invz_tensor/invzt_vertex4.m` is the dense component-labelled path-sum
   reference. Its own contract disables the attempted factored implementation
   because the factorization lost KMS/initial-state structure.
4. `invz_tensor/invzt_chi0_split.m` implements the retired A1 dominant/rest
   partition. The multilevel producer must not depend on that partition as its
   physical closure.
5. `docs/execution/invzt_ewald_certification.md` and
   `invz_tensor/tests/test_invzt_jq_tensor_ewald.m` certify the full-tensor
   lattice convention and default. Multilevel work consumes that lattice
   object; it does not create a second dipolar convention.

## Retained historical sources pending reproduction

| Git object | Retained evidence | Required treatment |
|---|---|---|
| `93f20d9:three_level_1z_extension.html` | Reviewed finite-level derivation: connected spectral cumulants, repeated-node limits, KMS boundary layers, ordered three-point/tadpole terms, and two-level reduction. | Recover equations selectively into the theory contract; independently check signs and normalization. |
| `09f0d3a:invz_functional/invzf_multilevel_cumulant.m` | Schema-v2 exact retained-rank C2/C3/C4 oracle using block-triangular matrix exponentials; guards multiplet cuts and grades reality. | Recover only the local oracle into a tensor-owned namespace; do not revive the abandoned functional package wholesale. |
| `09f0d3a:invzp_functional_prototype_record.md` | Two-level state-path agreement to `1.8e-15`; full 136-state static C4 fixture; rank-48/64 dynamic C4 drifts of order `1.1e-4`--`1.6e-4`; `functional_use_authorized=false`. | Treat numbers as historical until reproduced by branch-local tests. |
| `fcfd75f:invz_functional/invzf_electronuclear_cumulant.m` | Wiring from the 136-state single-ion solution to the oracle and a measured dense-backend failure around large `beta*energy_span` for nonzero-frequency C2. | Preserve the fail-closed status check and replace the unsafe automatic backend threshold before use. |
| `fcfd75f:docs/execution/invzp_exec_wp1_gate.m` | Dense/action, Hermiticity, static derivatives, permutation, KMS/reality, high-frequency, and rank-ladder gates. | Translate useful gates into tensor-owned tests; do not import projected execution policy. |

Historical measurements are leads, not branch-local certification. The first
work package reproduces the independent identities before the oracle can feed a
lattice map.

## Theory constraints carried into the active plan

1. The multilevel object is an additive connected-vertex correction
   \(\mathcal V\), not a scalar ratio \(\Sigma=\mathcal V/G_0\). A multipole
   \(G_0\) has zeros between poles, so the ratio creates artificial poles.
2. The local C4 must retain all time orderings, KMS endpoint terms, disconnected
   subtractions, and Hermite limits at repeated energies/frequencies.
3. Rank truncation is controlled by observable ladders and complete multiplets.
   Discarded Boltzmann weight alone does not bound virtual intermediate states.
4. In an ordered/source-biased state, connected C3 and one-point/tadpole terms
   are generally nonzero. A zero-source PM formula cannot be reused unchanged.
5. A fixed point is not a thermodynamic selector. Until a common functional or
   another derived selector is verified, coexisting ordered roots remain
   accepted-but-unranked.
6. Ewald provenance, lattice grid, local-state rank, Matsubara cutoffs, and
   cumulant backend are part of every saved result's identity.

## Explicitly rejected routes

- Dividing an additive multilevel correction by `G0` to force it into the A1
  scalar self-energy form.
- Treating each transition as an independent two-level system and summing its
  self-energy; shared states generate connected cross-transition cumulants.
- Restoring the dense \(O(N^4)\) A3 production route without new exact
  compression evidence and a passed cost gate.
- Enabling the old unproven `invzt_vertex4` factored path.
- Reviving the incomplete Gaussian/local-trace functional skeleton as a phase
  selector.
- Using convergence, residual size, seed order, or continuation direction as a
  free-energy ranking rule.
- Replacing nonfinite results with floors, clipping, or accepted NaNs.

## Open questions owned by the plan

- Whether a longitudinal multilevel vertex embedded in the full tensor medium
  is sufficiently accurate, or mixed Cartesian cumulants are required.
- Which finite-level contraction can meet production cost while preserving the
  exact two-level and KMS identities.
- Whether an ordered multilevel common functional exists for the chosen
  approximation, and if not, how far PM-only output remains useful.
- How to construct the real-axis vertex without a numerically ill-conditioned
  analytic continuation or a frozen-static shortcut.
