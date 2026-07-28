# Diagnostic numerical oracles

This directory contains only reusable, current diagnostic machinery. It is
not production code and must not be added recursively to the MATLAB path.

## Retained

| directory | purpose |
|---|---|
| `invzp_solver_stability_2026-07-27/` | Diagnostic-only ordered-node equations, Jacobian, and Newton oracle used to isolate the noncontractive nested iteration. |
| `biased_smooth_r_2026-07-28/` | Guarded root enumeration, pseudo-arclength/fixed-coordinate continuation, event evidence, graph assembly, and smooth-`r` selection oracles for the documented low-field backup route. |
| `invzp_qcp_grid_2026-07-28/` | Coupling-only and state-only `12^3--24^3` QCP grid gate, compact edge-pair trace, and Jensen area-rule oracle. |

Each retained directory has its own scope and acceptance notes. Neither is
wired into the production spectra dispatcher.

## Retired historical probes

The raw `diag*.m/.log`, `task*.m`, and synthetic `probe_qcp.m` files were
removed on 2026-07-28. Their measured conclusions had already been
consolidated into the durable reports, several scripts depended on deleted
scratch `.mat` files, and all hard-coded a workstation path. Keeping them in
the working tree therefore suggested a reproducibility they no longer had.

They remain recoverable from Git:

```text
git show 31a7fd0:docs/diagnostics/<path>
git restore --source=31a7fd0 -- docs/diagnostics/<path>
```

The authoritative current documents are:

- `invzp_convg_diagnosis.md`
- `invzp_convg_fix.md`
- `non-convg_mechanism_diagnosis_Codex.md`
- `biased_convergence_solution.md`
- `docs/invzp_strict_medium_gate0_report.md`
- `docs/invz_ordered_residual_contract.md`

Fresh test or probe scripts should be created only for a specific gate,
removed after their durable result is recorded, and must not rely on
untracked cached inputs without recording how those inputs are rebuilt.
