# Execution diary — `blind_convg_plan.md` Part I + `invzp_convg_fix.md`

Branch: `invzp-exec-convg-plan` (forked from `invzp-stage2c-diagnostic` @ `0eeee64`).
Started 2026-07-29.

Purpose: one append-only record of *what was tried, what it produced, and what it
taught*. Each entry is a checkpoint boundary: the repo is committed at the end of
every entry so the next attempt can be regressed cleanly with
`git revert`/`git reset --hard <checkpoint>`.

## Governing constraints (do not relax without amending the plan)

- **Promotion boundary** (plan §4): no legacy-equation diagnostic — DOS
  prescription, defactored continuation, discovered branch, TI score, polished
  root — may alter production phase labels until the stationary-functional,
  completeness, stability, and regression gates pass. Every change lands
  default-off / opt-in.
- **Repo standards**: no pole floors, no clipped denominators, no
  masked-to-bare substitution, no seed-dependent phase labels, no
  smoothness-as-thermodynamics.
- **Claim labels**: every recorded result is tagged *verified fact* /
  *inference* / *assumption* / *hypothesis* / *unresolved contradiction*.
- **Priority** (plan §6): fix-doc §12 (WP0–WP3, exact-cluster coefficient
  harness) is the critical path. Everything in the S1–S5 band is a bounded
  side-packet that must not displace it.

## Order of work

| # | Packet | Kind | Status |
|---|---|---|---|
| S1 | Classify the 3.825 T node failures by proximate cause | measurement | see entries |
| S2 | Full-closure defactoring (plan §2.6) | cheap fix attempt | see entries |
| S3 | Newton-polish the 4.05 T twins (blind WP1c) | measurement | see entries |
| S4 | Closed-loop path-dependence test | measurement | see entries |
| S5 | Document the `'last'`-crossing convention | docs | see entries |
| S6 | WP0 freeze the microscopic contract | critical path | see entries |
| S7 | WP1 exact local source oracle | critical path | see entries |
| S8 | WP2/WP3 diagram manifest + exact-cluster harness | critical path | see entries |

Rationale for starting at S1 rather than at the critical path: the presenting
symptom is a mask, and the plan permits days-scale legacy diagnostics to run
alongside the weeks-scale WP0–WP3. A classification of *why* the 17/34 nodes
fail costs one solver run and decides whether any cheap fix exists at all. Per
the operating policy, a measurement precedes machinery.

---

## Entries

<!-- newest last; each entry: attempt / outcome / learned / checkpoint -->
