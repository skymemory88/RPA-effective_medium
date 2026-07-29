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

### E1 — S1 measurement: what actually fails at 3.825 T

**Attempted.** Reproduce the documented 3.825 T failure and classify every h-node by
proximate cause before proposing any fix. Harness:
`docs/execution/invzp_exec_s1_failure_census.m` (per-node census through
`invz_ordered_trace` → `invz_hmf_ordered`, production `jensen`/`resummed` path),
`invzp_exec_s1_iter_anatomy.m` (outer-iteration history of one node),
`invzp_exec_s1_config_sweep.m` (the same census over damping/budget/acceleration
configurations). Fixed conditions throughout: T = 0.1 K, 16³ grid, dpRng 30,
`transverse_mf = legacy_x`, `hyp = true`, cold predictor + warm node-to-node seeding.
Logs and `.mat` under `docs/execution/out/`.

**Outcome — verified facts.**

1. *The documented failure count is a property of the iteration configuration, not of
   the field.* At Bx = 3.825 T the census gives:

   | `mix_outer` | `max_outer` | accel | status | failed | h\* |
   |---|---|---|---|---|---|
   | 0.20 | 1000 | none | `node_failed` | 11/34 | — |
   | **0.30** | **1000** | none | **`ok`** | **0/43** | **0.0172215** |
   | 0.40 | 1000 | none | `node_failed` | 1/34 | — |
   | 0.50 | 1000 | none | `node_failed` | 2/34 | — |
   | 0.70 | 200 | none | `node_failed` | 22/34 | — |
   | 0.70 | 1000 | none | `node_failed` | 17/34 | — |
   | 0.40 | 5000 | none | `node_failed` | 1/34 | — |
   | 0.40 | 1000 | `signed_aitken1` | `node_failed` | 1/34 | — |

   The `mix_outer = 0.70, max_outer = 1000` row reproduces `invzp_convg_diagnosis.md`
   §5's "17 of 34 sampled nodes fail" exactly, which identifies the configuration that
   measurement was taken under.

2. *The 3.825 T column is not unsolvable by Picard.* At `mix_outer = 0.30` every node
   converges, the full A–D contract is satisfied at each consumed node, the bisection
   and root phases run (43 node evaluations rather than 34), and the column returns an
   accepted ordered root h\* = 0.0172215 with a positive stability classifier. The
   **hard core — nodes failing under every configuration tried — is empty.**

3. *Every failure is `max_iter`; none is a domain rejection.* Across all 34 nodes in
   every configuration, `medium_status` is `not_applicable` (the normal resummed
   value) and `term_reason` is `converged` or `max_iter`. No
   `medium_out_of_domain`, no `degenerate_doublet`, no non-finite state. The static
   closure itself is essentially always converged at the point of failure
   (`resid_B ~ 1e-11` while `resid_A ~ 3.7e-3`): what fails is the **outer Σ↔K map**,
   not the static EMT sub-solve.

4. *The binding node's failure mode is an aperiodic wander through the Gstat pole.*
   For node 23 (h = 0.00241) at `mix_outer = 0.40`, over 1000 outer iterations:
   `resid_map` reached 5.40 × 10⁻⁵ at iteration 803 and then **grew again**
   (log₁₀ slope +1.64 × 10⁻³ per iteration over the last 200); ΔK₀ alternated sign
   139 times in 199 steps with lag-1 autocorrelation −0.521 but **no dominant period**
   (lag-2…8 all |ρ| < 0.07); `gstat_local_denom` ranged over **[−2.91, +2.97]**, i.e.
   the iterate crossed the local Gstat pole repeatedly; `Gstat` swung over
   [−2413, +5253]; `Dq_neg_count` wandered between 96 and 16248 of 16384 modes and
   `Dq_abs_min` fell to 2.4 × 10⁻⁵.

**Learned.**

- The presenting symptom at 3.825 T is **iteration dynamics**
  (`invzp_convg_diagnosis.md` §7.3), not the resummed-denominator domain boundary
  (§7.2). §7.2's "many deep-ordered masks are raised here" does not describe this
  field. This is a correction to the diagnosis's framing of *this* data point, and
  it is recorded as such rather than folded away.
- §5's "Both cold and warm field seeding leave the visible sliver unconverged" is
  true at `mix_outer = 0.70` and **false at 0.30**. Seeding was never the free
  variable that mattered; damping was.
- `signed_aitken1` changes nothing here, which is consistent with — and independently
  confirms — its narrow contract: the failure has no clean alternating mode for it to
  extrapolate.
- Because the mask/no-mask outcome moves with `mix_outer`, a phase label read off this
  column is iteration-configuration dependent. That is exactly the class of result the
  promotion boundary exists for: it is evidence that the mask was numerical, **not**
  evidence that h\* = 0.0172215 is the equilibrium root. Nothing here ranks roots or
  proves completeness, and no production default is changed on the strength of it.

**Open question carried forward.** Do the configurations that *do* close this column
close it to the *same* root? If different `mix_outer` values return different h\*, the
column is a live demonstration of iteration-dependent phase output. Measured in E3.

**Checkpoint.** `ckpt-01-s1-census`.

### E2 — S5: the `'last'`-crossing convention is now documented

**Attempted.** Plan §4's "latent defect, documented not patched": record the
`find(..., 1, 'last')` root-selection rule at `invz_hmf_ordered.m:362` in the code
itself. Verified first that it was genuinely undocumented — the only nearby comment
(`invz_hmf_ordered.m:413`) justifies the *increasing*-crossing sign pattern, and
nothing anywhere justified taking the **largest-h** crossing when several exist.

**Outcome.** A comment block above the selection line now separates the two choices
made there: the sign pattern (derived — Landau minimum, not maximum) and `'last'`
(a convention with nothing behind it), cites the multi-root evidence that makes the
convention load-bearing (diagnosis §6.1's seven h = 0 roots at 1.5 T, §6.2's two
states that both pass A–D at 4.05 T), and records that it must not be swapped for an
argmin-Φ score and is replaced only by `invzp_convg_fix.md` WP7. **No behaviour
change**; the executable line is byte-identical.

**Learned.** The convention is not merely undocumented, it is *silently* load-bearing:
`prof.status` reaches `'ok'` and a root is published without any record that a choice
among admissible crossings was made. Any future certificate vocabulary (plan §2.7)
needs a `multiple-crossings-present` flag emitted from this site, otherwise the
interim reports cannot distinguish `unique-accepted-within-budget` from
`multiple-accepted-unranked` at all. Logged as an open item; not implemented here,
because emitting it changes the trace schema.

**Checkpoint.** `ckpt-01-s1-census` (same commit).
