# Diagnostic and verification scripts

One-off MATLAB probes written during development. They are **not** production code and **not** on the
MATLAB path.

**Archived 2026-07-27** from the retired `.superpowers/sdd/` scratch directory (git-ignored, so not
recoverable from history). **Purged 2026-07-27**, same day: 71 of the 88 archived scripts were removed
as superseded process scaffolding for closed tasks. What remains is (a) the scripts that the surviving
committed documents cite as provenance, and (b) the current open investigation.

Everything removed was committed before the purge and is recoverable:

```
git log --oneline -- docs/diagnostics            # find the pre-purge commit
git checkout <commit> -- docs/diagnostics/<file>
```

## Path safety

Do **not** `addpath(genpath(<repo root>))` — that would put this directory on the path. Every script
here sets its own absolute `addpath` calls and is meant to be run directly with `run(...)`. The scripts
are unmaintained and may drift from the APIs they call.

## Current investigation

| item | what it is |
|---|---|
| `claude_convg_2026-07-27/` | the ordered-leg non-convergence diagnosis written up in `../../invzp_convg_diagnosis_Claude.md`: `diag1`–`diag11` plus their logs, on the real 16³ multiset at T = 0.31, 0.10 and 1.00 K. `diag7`/`diag8`/`diag9` are the ones that bear on the Gate-0 verdict (head-to-head resummed-vs-strict, B_c(T), and the Gate-0 field set re-run at the temperature it fits). `diag10` verifies the Stieltjes reduction of both medium sectors and the zero-parameter `x` vs `supp(ρ)` computability criterion; `diag11` re-measures the independent QCP claims now consolidated in §9.5 (the response evaluator is healthy given a state; the exact Γ-exclusion gap) |
| `invzp_solver_stability_2026-07-27/` | solver-only follow-up on the unchanged resummed equations: positive-control Jacobian oracle, traversal-direction test, local Picard-to-residual-Newton repairs, and the moderate-field branch gap. The retained Newton kernel is diagnostic-only and is not wired into HMF quadrature |
| `biased_smooth_r_2026-07-28/` | pure implementation oracle for the documented non-default smooth-`r(h)` selector. It ranks only pre-enumerated, fully certified paths and performs no continuation, smoothing, endpoint solving, or production dispatch |

## Provenance for numbers in committed documents

Each entry names the document it backs. Scripts whose backing documents were deleted on 2026-07-27
were removed with them.

| script | backs |
|---|---|
| `task18_g7.m` | the G7 scheme-jump table in `../invzp_strict_medium_gate0_report.md` §4. Note it calls `invz_emt_scalar` directly with `Sigma = zeros(size(wn))` — an undressed, primitive-level comparison — and its denominator `abs(a.K(2)-a.K(1))` is Matsubara-frequency dispersion, **not** q-dispersion |
| `task18_merge_report.m`, `task18_verify_g5g17.m`, `task18_nestcheck.m` | the Gate-0 aggregation, the G5/G17 measurements, and the grid-nesting check (`../invzp_strict_medium_gate0_report.md`) |
| `task12_verify.m` | the independent controller reproduction of the `+55 %` root shift on the synthetic fixture; this is a synthetic wiring result, not a real-density validity claim |
| `task17_g11_detail.m`, `task17_g11_bareref.m`, `task17_characterize.m`, `task17_g3_mcheck.m` | the Task-17 G11 anchor, three-scheme comparator, and G3 characterisation recorded in `../../Task-17_failure.md` |
| `task14_g9.m`, `task14_g9_compare.m`, `task15_g9.m`, `task15r_g9_cmp.m`, `task15r_g9_run.m` | the G9 default-path bit-identity comparisons |
| `probe_qcp.m` | the phase-dispatch probe behind the original QCP investigation, now consolidated in `../../invzp_convg_diagnosis_Claude.md` §9.5. **Caveat:** it runs on a 24-point synthetic `Jnu = linspace(...)` lattice — exactly the fixture the record flags as having masked two separate production failures. Do not read its numbers as real-lattice results |

`suite.m` — the runner behind the authoritative pass/fail counts quoted in the committed reports — was
removed on 2026-07-27 together with the suites it ran (see below). Its recorded counts stand as
historical measurements; they are no longer reproducible from this tree.

## What was not archived (deleted with the scratch directory, all regenerable)

- **~140 MB of `.mat` run data** (`task18_chunk_*.mat`, `task18_final_rep.mat`, `task18_couplings.mat`,
  `task2_matrix_results.mat`, `phase1_checkpoint.mat`, …). The coupling fixture is digest-pinned
  (`ddb9532d…`, verified on every invocation) and the drivers are committed (`invz_gate0_report.m`,
  `invz_static_domain_scan.m`); every measured value is quoted at full precision in the reports.
- **65 `review-*.diff` review packages** — pure `git diff A..B` output.
- **118 task briefs, implementer reports and review write-ups** — process scaffolding, consolidated
  into `../INVZ-DEVELOPMENT-RECORD.md`.
- **The execution ledger** `progress.md` — consolidated into the same file.

## Production-tree script cleanup, 2026-07-27

Separately from the diagnostics purge above, the production trees were stripped to production code
only, ahead of the ordered-leg fix attempt (`../../invzp_convg_fix_Claude.md`). 144 files removed, all
committed beforehand and recoverable with `git checkout HEAD -- <path>`:

| removed | count | reason |
|---|---|---|
| `invz_projected/tests/` | 93 | the MATLAB regression suites, their fixture builders, and `exploratory/`. Removed by request — the fix attempt builds its own temporary tests as it goes |
| `invz_tensor/tests/` | 29 | same, incl. `interop/` parity tests and the `fixtures/*.json` oracles. `vertex_oracle.json` is regenerable from `../../verify_tensor_vertex.py` |
| `invz_projected/invz_phase1_*.m` (8 of 9) | 8 | the BZ-quadrature/Γ audit harness. Unreachable from any driver; its prereg and report were deleted the same day. **`invz_phase1_qgrid.m` was kept** — it is live production, called by `invz_bz_couplings.m:62` for the offset/convention/gammaPolicy grid route |
| `invz_projected/invz_task2_*.m` | 13 | the stage-2c discriminator matrix harness. Whole chain unreachable; its prereg and report were deleted the same day |
| `docs/diagnostics/suite.m` | 1 | ran `runtests('invz_projected/tests')`; inoperable once the suites went |

Verified after removal: no surviving file makes a non-comment reference to anything deleted.

`invz_ordered_trace.m` and `invz_ordered_trace_resolve.m` lost their only caller
(`invz_task2_run_config`) in this cleanup and were **deliberately kept** — they wrap the
behaviour-neutral `opts.trace` hook in `invz_hmf_ordered`, which Phase 0.3 of the fix plan builds on.
They are now standalone entry points. Also kept as documented entry points with no internal caller:
`invz_cfrot.m`, `invz_deltaF_ordered.m`, `invzt_vertex3.m`, `invzt_vertex4.m`, `invzt_report_ladder.m`.

22 comment lines across 17 surviving production files still name deleted test files as the
regression that pinned a given identity (e.g. `invz_common/invz_gstat_ordered.m:17`). These were left
in place: they record *which invariant* the code was pinned to, which outlives the file that checked it.
