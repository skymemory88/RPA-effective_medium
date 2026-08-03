# Full-tensor multilevel 1/z execution journal

- Plan: [`../../invzt_multilevel_1z_next_steps.md`](../../invzt_multilevel_1z_next_steps.md)
- Documentation index: [`../invzt_documentation_index.md`](../invzt_documentation_index.md)

## Recording policy

This is the append-only execution record for the multilevel Jensen-style
\(1/z\) program. Each attempt records the question, exact change or command,
decisive observation, evidence location, verdict, rollback status, and next
action. Large logs and binary artifacts belong under `docs/diagnostics/` and are
referenced rather than pasted here.

Evidence labels are:

- **Verified:** reproduced on the current branch by a named check.
- **Historical:** retained at a named Git object but not reproduced here.
- **Inference:** logically supported but not directly measured.
- **Hypothesis:** proposed explanation awaiting a discriminating test.
- **Unresolved:** proof obligation not yet satisfied.

No checkpoint may upgrade historical evidence or inference to verified without
a recorded branch-local check.

## Standing constraints

1. Work in `invz_all_levels/`, the peer multilevel module in the host checkout.
   The retired `invzt-multilevel-1z` worktree is retained only by archival tag
   `archive/invzt-multilevel-1z-20260802`.
2. Do not edit `invz_projected`. Treat `invz_common` and `invz_tensor` as
   shared inputs unless a shared promotion is explicitly authorized and tested.
3. The certified Ewald implementation remains the tensor lattice default.
   Brute force is an explicit diagnostic comparator only.
4. The A1/A2/A3 programs are retired on this branch. Historical source may be
   read or selectively recovered, but production dispatch is not revived by
   implication.
5. Never conceal a failed gate with floors, clipping, NaNs, broader tolerances,
   seed order, or an unrecorded fallback.
6. Preserve exact rollback boundaries. Every code work package lands separately
   after its focused tests pass.
7. Five consecutive variants that fail the same proof obligation trigger a
   return to the governing equation/representation before a sixth is attempted.
8. Pause for human direction before changing material conventions, accepting a
   non-derived phase selector, relaxing a preregistered production budget, or
   adopting a labelled approximation as an exact full-tensor result.

## Baseline

- **Verified:** branch `invzt-multilevel-1z` began this program at
  `425143d8633a32457b21c1aca0024b821064adf3`
  (`feat(invzt): certify Ewald as tensor default`).
- **Verified:** `invz_tensor/invzt_jq_tensor.m`, spectra, phase-diagram, and
  ladder drivers select Ewald by default; the focused certification is
  `invz_tensor/tests/test_invzt_jq_tensor_ewald.m` and the retained record is
  `docs/execution/invzt_ewald_certification.md`.
- **Verified:** the current tensor tree contains the old A1/A2/A3 code as
  reference/legacy implementation, but it has no active multilevel local
  cumulant package and no multilevel production mode.
- **Verified:** before this documentation cleanup, the tensor worktree was
  clean at the baseline commit.
- **Verified:** this pass performs no production-code edit.
- **Protection audit:** at the start of the documentation edit, the projected
  worktree was on `invzp-exec-convg-plan` at `29d4117e55ba0a2e25543f5a2a19dea195169741`;
  its porcelain-status SHA-256 was
  `3bf0fab6131abc56cf3a75386c695b3669bd36312b4ce2ddac2d127f50bd505f`.
  The same digest must be observed at completion of this pass.

## Retained historical evidence, not yet reproduced

- `09f0d3a:invz_functional/invzf_multilevel_cumulant.m` contains a finite-level
  connected C2/C3/C4 oracle using block matrix exponentials.
- `09f0d3a:invzp_functional_prototype_record.md` reports two-level oracle
  agreement, one full-rank static C4 fixture, and useful rank-ladder results,
  while explicitly leaving functional use unauthorized.
- `fcfd75f:invz_functional/invzf_electronuclear_cumulant.m` records a
  nonfinite dense-backend regime at large `beta*energy_span` for a nonzero
  Matsubara C2; exponential action was the observed workaround.
- `93f20d9:three_level_1z_extension.html` contains the reviewed multilevel
  derivation and ordered three-point/tadpole requirements.

These facts justify a recovery-and-audit work package, not direct production
reuse.

## Decision ledger

| Date | Decision | Basis | Revisit condition |
|---|---|---|---|
| 2026-08-02 | Preserve `invzt_a_tier_guide.html` as clearly labeled historical reference. | Explicit user request; it remains pedagogically useful. | User requests removal or replacement. |
| 2026-08-02 | Retire A1/A3 diagnoses, plans, and journals from this branch. | They describe a stopped convergence repair and dense-A3 route, not the chosen multilevel program. | New evidence makes one a necessary source; recover only the needed section with provenance. |
| 2026-08-02 | Keep `jensen_1z_framework.html` as two-level reduction reference. | It fixes conventions and the exact Jensen limit. | Superseded by a reviewed derivation with equal or stronger coverage. |
| 2026-08-02 | Do not duplicate the 2.06 MiB historical three-level HTML into the active tree. | Git preserves the exact artifact; a compact evidence map prevents a broken or competing plan. | Offline distribution requires a self-contained copy. |
| 2026-08-02 | Recover only the finite-level cumulant oracle, not the abandoned functional package wholesale. | The oracle has independent identities; the broader Gaussian/local functional route had unresolved structural defects. | A new derivation proves a common functional and its exact correspondence. |
| 2026-08-02 | Stage longitudinal multilevel dressing before optional mixed-Cartesian promotion. | It directly tests the Jensen generalization with one local operator and avoids the known dense component-labelled A3 cost. It remains explicitly an approximation until its omission error is bounded. | Mixed-channel discriminant fails or full component route becomes affordable. |

## Rejected-route ledger

| Route | Rejection evidence | Status |
|---|---|---|
| Scalar `Sigma = V/G0` for a multipole propagator | Multipole zeros create artificial ratio poles; additive `V` avoids the division. | Rejected representation. |
| Independent two-level corrections for each transition | Shared levels produce connected cross-transition cumulants. | Rejected physics approximation unless separately derived and labeled. |
| Existing dense component-labelled A3 as production engine | `invzt_vertex4` is \(O(N^4)\); earlier cost gate stopped the route and its factored version is unproven. | Retained only as small-system reference. |
| Old Gaussian/local-trace functional skeleton | Historical work did not establish the required nonlinear common functional/physical sheet. | Do not revive without a fresh derivation. |
| Fixed-point convergence as phase selector | Convergence does not rank coexisting thermodynamic states. | Permanently rejected selector. |
| Brute-force dipolar sum as default | Certified Ewald convention is present and tested. | Diagnostic comparison only. |

## Checkpoint 0 — 2026-08-02 — Documentation reset

**Question.** Can the tensor branch be reduced to one unambiguous multilevel
implementation program without touching the working projected branch?

**Changes.** Removed projected-spin convergence plans, journals, proposals, and
diagnostic reports from this tensor branch while retaining its module README as
reference-only; imported the user-requested A-tier guide with a historical
banner; labeled the Jensen document as reference; added the documentation
index, theory evidence map, this journal, and the multilevel next-steps plan.

**Outcome.** Passed. The active set contains eight Markdown/HTML documents and
one authoritative next-steps plan. YAML parsed as
`invzt_multilevel_1z_plan/v1` with all 11 gates present; every local document
link resolved; all three retained HTML documents passed internal-anchor checks;
and `git diff --check` passed. The nine focused Ewald MATLAB tests passed in
7.3 s. Historical Git objects for the reviewed derivation, local oracle,
electronuclear wrapper, identity gate, and removed projected documents were
confirmed recoverable. The reviewed three-level HTML is 2,159,385 bytes
(2.06 MiB), which supports indexing rather than duplicating it.

The projected worktree remained on `invzp-exec-convg-plan` at `29d4117`; its
final porcelain-status SHA-256 was identical to the protection-audit digest,
`3bf0fab6131abc56cf3a75386c695b3669bd36312b4ce2ddac2d127f50bd505f`.

**Rollback.** All documentation changes are relative to baseline `425143d` and
contain no production-code edits. Revert the eventual documentation commit to
return exactly to that baseline documentation state.

**Next action.** After checkpoint 0 is committed, begin ML0 task 1: recover the
schema-v2 local cumulant oracle into tensor-owned names without production
wiring, then establish the independent two-level and backend tests before any
lattice-map work.

## Attempt template

Copy this section for each new attempt.

### Checkpoint N — YYYY-MM-DD — short title

**Proof obligation.** The single bounded claim being tested.

**Pre-state.** Commit, options, fixtures, cache state, and relevant artifact
hashes.

**Change or command.** Exact files/symbols changed and reproducible command.

**Observation.** Decisive numbers, failures, warnings, and evidence paths.

**Classification.** Verified fact, inference, hypothesis, or unresolved
contradiction.

**Verdict.** `pass`, `fail`, `stop_cost`, `stop_theory`, `diagnostic_only`, or
`needs_human_decision`.

**Rollback.** Commit/revert boundary and whether rollback was exercised.

**Next action.** One bounded next proof obligation.
