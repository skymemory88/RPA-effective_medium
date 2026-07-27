# Diagnostic and verification scripts (archive)

**Archived 2026-07-27** from the retired `.superpowers/sdd/` scratch directory, which was git-ignored
and therefore not recoverable from history. These are one-off MATLAB probes written during
development — controller verifications, red/green captures, gate measurements, and reproductions of
published numbers. They are **not** part of the production code and **not** on the MATLAB path.

They are kept because several are the provenance for numbers quoted in committed documents.

## Path safety

Do **not** `addpath(genpath(<repo root>))` — that would put this directory on the path. Every script
here sets its own absolute `addpath` calls and is meant to be run directly with `run(...)`. Filenames
were checked for collisions against the rest of the repo at archive time; there were none, but the
scripts are unmaintained and may drift from the APIs they call.

## The ones that back committed claims

| script | backs |
|---|---|
| `task18_g7.m` | the G7 scheme-jump table in `invzp_strict_medium_gate0_report.md` §4 and `invzp_convg_fix_attmpt_claude.md` §5.5. Note it calls `invz_emt_scalar` directly with `Sigma = zeros(size(wn))` — an undressed, primitive-level comparison — and its denominator `abs(a.K(2)-a.K(1))` is Matsubara-frequency dispersion, not q-dispersion |
| `task18_merge_report.m`, `task18_verify_g5g17.m`, `task18_nestcheck.m` | the Gate-0 aggregation, G5/G17 measurements, and the grid-nesting check |
| `task12_verify.m` | the independent controller reproduction of the `+55 %` root shift (`invzp_convg_fix_attmpt_claude.md` §4.1) |
| `task17_g11_detail.m`, `task17_g11_bareref.m`, `task17_characterize.m`, `task17_g3_mcheck.m` | the Task-17 G11 anchor, the three-scheme comparator, and the G3 characterisation |
| `task14_g9.m`, `task15_g9.m`, `task14_g9_compare.m`, `task15r_g9_cmp.m` | the G9 default-path bit-identity comparisons |
| `suite.m` | the suite runner used for the authoritative counts |

## What was NOT archived

Deleted with the scratch directory, all regenerable:

- **~140 MB of `.mat` run data** (`task18_chunk_*.mat`, `task18_final_rep.mat`, `task18_couplings.mat`,
  `task2_matrix_results.mat`, `phase1_checkpoint.mat`, …). Regenerable: the coupling fixture is
  digest-pinned (`ddb9532d…`, verified on every invocation) and the drivers are committed
  (`invz_gate0_report.m`, `invz_static_domain_scan.m`). Every measured value is quoted at full
  precision in the committed reports.
- **65 `review-*.diff` review packages** — pure `git diff A..B` output, regenerable verbatim.
- **118 task briefs, implementer reports and review write-ups** — process scaffolding. Their durable
  content is consolidated into `../INVZ-DEVELOPMENT-RECORD.md`.
- **The execution ledger** `progress.md` — consolidated into the same file.
