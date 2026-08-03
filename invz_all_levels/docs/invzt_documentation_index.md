# Full-tensor documentation index

- Status date: 2026-08-02
- Module: `invz_all_levels`
- Retired branch archive: `archive/invzt-multilevel-1z-20260802`
- Baseline: `425143d8633a32457b21c1aca0024b821064adf3`

## Purpose

This index prevents historical A1/A2/A3 material from being mistaken for an
active work order. The only active implementation program in this module is the
multilevel Jensen-style \(1/z\) program. `invz_projected` and `invz_tensor` are
peer modules; their current code is not changed by this migration.

## Active documents

| Document | Authority |
|---|---|
| [`../invzt_multilevel_1z_next_steps.md`](../invzt_multilevel_1z_next_steps.md) | Machine-readable implementation plan, gates, budgets, and completion criteria. |
| [`execution/invzt_multilevel_1z_journal.md`](execution/invzt_multilevel_1z_journal.md) | Append-only attempt, evidence, decision, and rollback record. |
| [`execution/invzt_ewald_certification.md`](execution/invzt_ewald_certification.md) | Certified Ewald convention, validation evidence, and production-default record. |
| [`reference/invzt_multilevel_theory_sources.md`](reference/invzt_multilevel_theory_sources.md) | Evidence map for current and historical theory/code sources that the plan may recover or test. |

If these documents disagree, the next-steps plan controls implementation and
the journal records what was actually done. A discrepancy must be resolved in
both before production authorization.

## Reference-only documents

| Document | Status |
|---|---|
| [`../../jensen_1z_framework.html`](../../jensen_1z_framework.html) | Two-level analytical derivation and reduction reference. Its historical implementation prose is not an active plan. |
| [`../invzt_a_tier_guide.html`](../invzt_a_tier_guide.html) | Historical A1/A2/A3 guide retained at the user's request. It explains the retired vocabulary and dense-A3 cost stop. |
| [`../../invz_projected/README.html`](../../invz_projected/README.html) | Reference-only module documentation for the user-certified projected implementation. |

The reviewed multilevel derivation formerly stored as
`three_level_1z_extension.html` is not duplicated into the active tree. Its
source-of-record is Git object
`93f20d9:three_level_1z_extension.html`; decisive contracts are indexed in
[`reference/invzt_multilevel_theory_sources.md`](reference/invzt_multilevel_theory_sources.md).
This keeps the active documentation small while preserving exact recovery.

## Retired or excluded material

The following A-tier documents are not active on this branch and must not be
used as implementation instructions:

- `invzt_converg_diagnosis.md`
- `invzt_convergence_next_steps.md`
- `docs/execution/invzt_convergence_journal.md`
- `docs/execution/invzt_projected_convergence_port_audit.md`
- `invzt_a3_next_steps.md`
- `docs/execution/invzt_a3_journal.md`
- `docs/execution/invzp_convergence_journal.md`
- `converg_diagnosis.md`
- `convg_solution_suggestion.md`
- `invzp_convergence_dead_ends.md`
- `invzp_convergence_next_steps.md`
- `long_term_code_design.md`
- `thermodynamic_check_theory.html`
- `docs/invz_projected_field_orientation_recovery.md`
- `docs/diagnostics/invzp_approximation_wp6/README.md`

The first six tensor A-tier files were deliberately not imported into the clean
full-tensor branch. The projected reports listed afterward were removed here
because they are not tensor-branch dependencies. Their exact state remains
recoverable at projected commit `29d4117`. These removals do not alter or
downgrade the projected-spin implementation; that separate branch is
user-certified as fully working.

## Maintenance rule

New implementation decisions go into the active plan and new measurements go
into the journal. Do not create parallel `diagnosis`, `next_steps`, or tier
journals for the same program. A superseded document must be reclassified here
in the same commit that introduces its replacement.
