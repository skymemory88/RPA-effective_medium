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

**Learned (E2).** The convention is not merely undocumented, it is *silently* load-bearing:
`prof.status` reaches `'ok'` and a root is published without any record that a choice
among admissible crossings was made. Any future certificate vocabulary (plan §2.7)
needs a `multiple-crossings-present` flag emitted from this site, otherwise the
interim reports cannot distinguish `unique-accepted-within-budget` from
`multiple-accepted-unranked` at all. Logged as an open item; not implemented here,
because emitting it changes the trace schema.

**Checkpoint.** `ckpt-01-s1-census` (same commit).

### E3 — S2 cheap fix attempt: full-closure defactoring. **Correct, gated, and it does not help.**

**Attempted.** Plan §2.6 / review §4: the `stable_form` reassociation fixes `Gtil0`/`r`
at the local `Gstat` pole but the *other* consumer,
`Gq = Gstat./(1 + (J(q) − K0).*Gstat)`, still evaluates Inf/Inf = NaN there. Hypothesis:
E1 showed the binding node's outer iteration crosses `gstat_local_denom = 0` repeatedly,
so removing that NaN might let the map walk through the pole and converge.

Implementation (all opt-in, default-off):
- `invz_common/invz_reciprocal_static_closure.m` — **one** definition of the q-average in
  the reciprocal coordinate, `Gq = 1/(1/Gstat + J(q) − K0)`, with `1/Gstat` built from
  `Gstat`'s own parts (`z = d0/(G0inel0 + xi*G0el0*d0)`) so no infinity is ever formed.
  Exact algebraic reassociation: no floor, no broadening, no added tolerance, no sign change.
- `invz_emt_static_ordered.m` — new `opts.closure_coordinate` (`'factored'` default |
  `'defactored'`). The **closure condition and `out.resid` are unchanged**
  (`|mean_q Gq − Gstat| < resid_tol` in both), so acceptance semantics do not move.
- `invz_ordered_residual.m` — its private `local_blockB_defactored` now calls the shared
  function instead of carrying a second copy, so audit and iteration cannot drift apart.
- `invz_gstat_ordered.m` — records `G0inel0`/`G0el0` (values only; `G0bare` is their sum
  and the split is not recoverable from it).

**Outcome.**

*Gates pass* (`invzp_exec_s2_defactor_gate.m`, fixture = the real T = 0.1 K / 3.825 T /
h = 0.00241 node, 16384 modes):
- G1 equivalence over the K₀ range the failing node visits: worst relative disagreement
  7.8 × 10⁻¹³ (typical ~10⁻¹⁵), and — the load-bearing part of the claim — the agreement
  is **best near the pole** (6.7 × 10⁻¹⁵ for |d₀| < 0.15) and worst away from it, i.e. the
  error does not track the pole. A tighter ulp-scaled budget was tried first and rejected
  as the wrong bar: both arrangements inherit a common-mode amplified rounding error from
  `d0 = 1 + Σ₀ + K₀·G0inel0`, which is itself formed by cancellation.
- G2 at the pole: the factored form is non-finite exactly at d₀ = 0, the reciprocal form is
  finite there and at every offset tested, and it reproduces the analytic limit
  `Gq → 1/(J(q) − K0)` to 2.2 × 10⁻¹⁵.
- G3 behaviour-neutrality: the default path is **bit-identical** to the explicit
  `'factored'` option in (K₀, Gstat, resid, iters), and the two coordinates reach the same
  fixed point to 1.0 × 10⁻¹⁵.

*Column effect: none.* `invzp_exec_s2_column_effect.m` at 3.825 T:

| `mix_outer` | `max_outer` | factored | defactored |
|---|---|---|---|
| 0.40 | 1000 | 1 failed, ids {23} | 1 failed, ids {23} |
| 0.70 | 1000 | 17 failed, same ids | 17 failed, same ids |
| 0.70 | 200 | 22 failed, same ids | 22 failed, same ids |
| 0.30 | 1000 | `ok`, h\* = 0.017221493 | `ok`, h\* = 0.017268978 |

The failure sets are **identical, node for node**, in every masked configuration. This
also re-verifies behaviour neutrality at the column level: the factored rows reproduce
E1's pre-change counts exactly.

**Learned — why the hypothesis was wrong.** The static closure iterates on K₀ at *fixed*
(λ, Σ₀); within one static solve K₀ stays in a well-conditioned region and the pole is
never actually reached. The pole crossings E1 measured are **outer**-loop events, where
`d0` moves because **Σ₀** moves between outer iterations. So the Inf/Inf site the review
identified is real but is not on the path that fails. Defactoring the inner q-average
therefore cannot address the outer map's non-contractivity, and it did not.

**Disposition.** Kept, default-off. It removes a demonstrated NaN failure mode (G2), it
de-duplicates the reciprocal-coordinate arithmetic, and it costs nothing when unused. It
is **not** presented as a convergence fix, and per the promotion boundary it may not
change a production phase label. Rejected as a fix for the 3.825 T masks — recorded in
the rejected-branch list below so it is not retried.

**Checkpoint.** `ckpt-02-s2-defactor`.

### E4 — Root and state sensitivity: two A–D-accepted states at one h, from a 10⁻¹² perturbation

**Attempted.** E1's open question — do the configurations that close the 3.825 T column
close it to the *same* root? `invzp_exec_e3_root_sensitivity.m` perturbs a closing
configuration along two axes that cannot change the fixed-point set: the arithmetic
coordinate (factored vs defactored, disagreeing by ≤ 8 × 10⁻¹³ per G1) and `mix_outer`
(damping changes iteration dynamics only).

**Outcome — verified facts.**

1. **Two distinct A–D-accepted states at the same h.** Comparing the two arithmetic
   variants node by node at `mix_outer = 0.30`, the difference is *not* a smooth
   accumulation. Nodes 1–4 differ only at the 10⁻¹⁰–10⁻¹² level (rounding). Then:

   | node | h | rel Δr | rel ΔK₀ | iters (fac / def) |
   |---|---|---|---|---|
   | 23 | 0.00241 | 1.38 × 10⁻² | **7.8 × 10⁻¹** | 70 / 35 |
   | 24 | 0.00299 | 1.74 × 10⁻² | **7.2 × 10⁻¹** | 26 / 92 |

   Both runs report **zero failed nodes**, so both states at each of these h passed the
   full A–D contract at `tol_outer = 1e-8`. A K₀ differing by 78% is not a tolerance
   artifact: these are two different accepted fixed points at one operating point,
   selected by a 10⁻¹²-level arithmetic difference upstream.

2. **The published root moves by about three root-tolerance units, not more.**
   h\*(factored) = 0.0172214927186, h\*(defactored) = 0.0172689777947 → 2.76 × 10⁻³
   relative. Across the damping window {0.28, 0.30, 0.32, 0.34, 0.36} — all of which
   close — h\* takes exactly four distinct values, {0.0172215, 0.0172373, 0.0172531,
   0.0172690}, **evenly spaced by 1.58 × 10⁻⁵**. That spacing is the bisection's final
   bracket width under `tol_root = 1e-3` (≈ 1.7 × 10⁻⁵ at h ≈ 0.0172), so the h\* spread
   is at the root finder's own reported precision. `mix_outer = 0.26` fails 2 nodes.

**Learned.**

- Fact 1 extends the diagnosis's coexisting-accepted-roots evidence
  (`invzp_convg_diagnosis.md` §6.2, recorded at 4.05 T) to 3.825 T, and strengthens it:
  there, two roots were found by *seeding*; here they are reached by an arithmetic
  perturbation smaller than the acceptance tolerance by four orders of magnitude. The
  A–D contract is not selective enough to distinguish them — exactly the gap
  `invzp_convg_fix.md` §7 exists to close.
- Fact 2 is the honest limit on how far fact 1 goes *at this field*: the coexistence
  happens at nodes far below the crossing (h = 0.0024, 0.0030 vs h\* = 0.0172), where the
  contribution to h₀ = ∫r dh′ is small, so the published root is not visibly destabilised
  here. It would be an overclaim to call the root chaotic on this evidence. What is
  established is that the *state path* is not determined, and that the current reporting
  cannot tell the user which of two accepted states each node took.
- This makes the missing certificate flag from E2 concrete: the interim vocabulary needs
  `multiple-accepted-unranked` to be *emittable per node*, not only per column.

**Checkpoint.** `ckpt-02-s2-defactor` (same commit).

---

## Rejected branches (do not retry without new contradicting evidence)

| Branch | Rejecting evidence |
|---|---|
| Defactoring the static q-average as a fix for the 3.825 T masks | identical failure sets node-for-node in all three masked configurations (E3); the pole crossings are outer-loop, not inner-closure, events |
| `signed_aitken1` as a fix for the 3.825 T masks | no change at `mix_outer = 0.40` (E1); the binding node has no clean alternating mode (lag-2…8 autocorrelation \|ρ\| < 0.07) |
| Reading the 3.825 T mask as a resummed-denominator domain rejection | `medium_status` is `not_applicable` at every node in every configuration; every failure is `max_iter` with `resid_B ~ 1e-11` (E1) |
| Larger iteration budget alone | `max_outer` 1000 → 5000 at `mix_outer = 0.40` changes nothing (E1) |
