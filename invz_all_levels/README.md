# invz_all_levels

`invz_all_levels` is the peer implementation module for the finite-level,
full-tensor Jensen-style \(1/z\) calculation.  Its target is the retained
136-state Ho3+ local spectrum coupled to the complete tensor dipole-dipole
medium; it is not a checkout of the whole repository and it does not own the
projected or legacy tensor implementations.

## Module boundary

This module sits alongside `invz_projected` and `invz_tensor`.

- New multilevel cumulants, contexts, maps, solvers, response code, drivers,
  tests, and diagnostics belong here.
- `../invz_common` supplies shared, stable local-ion and numerical utilities.
- `../invz_tensor` is the current shared full-tensor/Ewald baseline.  Do not
  overwrite it from historical all-level work; promote a shared change only
  deliberately and with its own validation.
- `../invz_projected` is a separate projected implementation and is not a
  dependency of the multilevel closure.

## Active material

The recovered roadmap is [invzt_multilevel_1z_next_steps.md](invzt_multilevel_1z_next_steps.md), with its execution journal and theory evidence
under `docs/`.  References in those recovered documents to the retired
`invzt-multilevel-1z` branch/worktree are historical.  Active work now occurs
in this module inside the host checkout.

The Ewald certification record in
`docs/execution/invzt_ewald_certification.md` is retained evidence from the
retired worktree.  The current `../invz_tensor` implementation is the shared
baseline requested for future work; rerun the focused Ewald gate before using
the old certificate as evidence for a changed implementation.

## Archive

The former `invzt-multilevel-1z` worktree and branch were retired during this
migration.  Its final commit is retained by the non-development tag
`archive/invzt-multilevel-1z-20260802` for provenance and recovery.
