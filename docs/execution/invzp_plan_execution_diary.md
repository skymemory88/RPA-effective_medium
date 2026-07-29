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

### E5 — Critical path: the exact-cluster coefficient harness (fix doc §12) is built and gated

**Attempted.** `invzp_convg_fix.md` §12's named first deliverable: "a small exact-cluster
coefficient harness tied to an automatically generated diagram manifest", implementing
§6's four steps for the two-site X–X bond.

A capability inventory of `invz_functional/` was run first, to avoid rebuilding what the
repo has. Result: the exact cluster oracles exist (`invzf_cluster_exact` for 1–8 sites,
`invzf_two_site_exact`, `invzf_centered_*`) but **only produce exact numbers** — no file
extracts coupling-order coefficients from them and none compares against a candidate;
exact rank-2/3/4 local cumulants exist (`invzf_multilevel_cumulant`) but are wired
end-to-end only for the two-level model, not the 136-state one; and **no diagram
enumeration or provenance machinery exists anywhere in the repo**.

**Delivered.**
- `invz_functional/invzf_taylor_extract.m` — certified Taylor extraction from a handle
  J ↦ F(J): Chebyshev-node polynomial fits repeated over a halving sequence of radii, so
  the **spread across radii is a measured error bar per coefficient** rather than an
  assumed accuracy. A comparator whose own error is unknown cannot grade anything.
- `invz_functional/invzf_pair_local_manifest.m` — the independent side. Predicts the pair
  free-energy coefficients from **local data only** (m and Cₙ at J = 0) via the
  linked-cluster cumulant expansion, with the derivation recorded in the header and each
  retained monomial carrying its topology and symmetry factor:
  a₁ = −m² (one line, C1 on each site); a₂ = −m²C₀ (two lines, C2 on one site and two
  disconnected C1 on the other, factor 1) − (1/2β)ΣₙCₙ² (two lines, C2 on each site,
  factor ½ — the single ring). Order 3 is derived **only at h = 0**, where the local
  X-flip unitary maps J → −J so every odd coefficient vanishes; order 3 at h ≠ 0 and
  order 4 are returned `derived = false` and are reported **ungraded, never passed**.
- `invz_functional/invzf_cluster_coeff_harness.m` — the comparator and regression fixture.
- `docs/execution/invzp_exec_wp3_harness_gate.m` — the gate.

**Outcome — all four gates pass.**
- **W1 extractor self-test** on exp(0.7J): orders 0–2 recovered to ≤ 8.3 × 10⁻¹²
  relative; the reported error bar bounds the true error at every order (worst ratio
  1.00). Order 4 is only ~2 × 10⁻⁶ accurate at r₀ = 0.05 — and the harness *reports* that,
  which is the property being tested.
- **W2 parity fixture** (h = 0): a₁ = −5.4 × 10⁻¹⁸, a₃ = 1.0 × 10⁻¹³, a₅ = −4.9 × 10⁻¹⁰,
  each inside its own error bar, against a₂ = −1.708. The extractor can certify a
  coefficient as zero, which is the harder half of a coefficient gate.
- **W3 manifest vs exact cluster**: match at orders 0, 1, 2 at both h = 0 and h = 0.05,
  worst |diff| 7.4 × 10⁻¹²; order 3 matches at h = 0 (parity). Orders reported `ungraded`
  are exactly those the manifest does not derive.
- **W4 negative control**: a self-consistent mean-field pair functional reproduces orders
  0 and 1 and is **rejected at order 2**, and its shortfall equals the predicted missing
  single-ring term −ΣₙCₙ²/(2β) to **1.26 × 10⁻⁷ relative**.

**Learned.** W4 is the load-bearing capability: the harness does not merely say a
candidate is wrong, it says **which diagram is missing and by how much**. That is what
makes it a cheap falsifier for a functional candidate, as §12 argues, and it is now
demonstrated rather than asserted. Two flaws in my own first verdict logic were found and
fixed by running it: a relative tolerance keyed to the coefficient's own magnitude
declared a correctly-vanishing coefficient "ungradeable" (fixed by keying tolerances to
the expansion scale), and a self-test threshold that ignored the reported error bar
counted an honestly-reported 10⁻⁶ accuracy as a failure.

**What remains open on the critical path, stated plainly.**
- WP2's *automatic* diagram generator does not exist and was not built. The manifest here
  is hand-derived and machine-checked, covering orders 0–2 at any source and order 3 at
  zero source. Orders 3 (h ≠ 0) and 4 are **frozen as exact numbers** in the regression
  fixture so whoever derives them has an immediate target: at (Δ, M, h, β) =
  (0.30, 1.00, 0.05, 20.0), a₃ = +0.448263449489 and a₄ = +6.16588863563.
- WP1 is complete for the two-level model but **not** for the 136 electronuclear states:
  `invzf_multilevel_cumulant` is a general exact rank-2/3/4 engine that is simply not
  wired to `invz_single_ion`'s spectrum and operator. That wiring, plus KMS/permutation/
  high-frequency identity checks, is the next bounded packet.
- The harness is candidate-agnostic (it takes a handle), so it will grade a 2PI candidate
  the day one exists, with no change.

**Checkpoint.** `ckpt-03-wp3-harness`.

### E6 — WP1: the 136-state local source oracle is wired, and its own gate **fails at nonzero frequency**

**Attempted.** Close the WP1 gap the inventory identified: `invzf_multilevel_cumulant` is
an exact rank-2/3/4 cumulant engine for any (E, O, β), and `invz_single_ion` already
exports the 136-state eigenvalues and Jz in that eigenbasis, but nothing connected them —
so C3/C4 for the electronuclear problem existed only as finite-difference h-derivatives
of C2. Delivered `invz_functional/invzf_electronuclear_cumulant.m` (the wiring, with the
source imposed through `hz_fixed` and `transverse_mf = 'none'` so the generator is local
and not mean-field-dressed) and `docs/execution/invzp_exec_wp1_gate.m` (the identities
WP1 names: KMS, permutation, Hermiticity, static derivative, high-frequency, plus a
backend cross-check and a rank ladder).

**Outcome — four of seven pass; two fail, and a third defect was found on the way.**

| gate | result | number |
|---|---|---|
| P0 backend cross-check at T = 0.6 K | PASS | dense vs exponential-action agree to 3.9 × 10⁻¹⁵ |
| P1 Hermiticity of Jz | PASS | 1.2 × 10⁻¹⁶ |
| P2 static derivative C₂(0) = dm/dh_z | PASS | 2.16 × 10⁻⁹ relative |
| P3 second derivative C₃(0,0,0) = d²m/dh_z² | PASS | 1.98 × 10⁻⁸ relative |
| P4 permutation symmetry, C₃ (6 orderings) and C₄ (24) | PASS | 1.9 × 10⁻¹⁵ / 3.3 × 10⁻¹⁵ |
| **P5 KMS/reality** | **FAIL** | C₂ perfect (5.7 × 10⁻¹⁶, Im = 0); **C₃(2,−5,3) vs conj C₃(−2,5,−3) differ by 2.43 × 10⁻³** |
| **P6 high-frequency tail** | **FAIL** | ωₙ²C₂(n) at n = 40, 80, 160, 320 is 30.67, 52.13, 107.15, 191.44 — **growing**, ratios 1.70, 2.06, 1.79, where it must approach a constant |

**Defect found while building the gate (independent of P5/P6).**
`invzf_multilevel_cumulant`'s default `dense_beta_span_limit = 2000` selects the dense
block-exponential backend up to β·span = 2000, but that backend returns
`status = 'nonfinite'` — a **NaN cumulant** — for **rank 2 at nonzero Matsubara labels**
once β·span exceeds roughly 10³. Reproduction: LiHoF4 ion, B = [4 0 0] T, h_z = 0.02 meV,
labels [3 −3]; `ok` at T ≥ 0.6 K (β·span 813) and `nonfinite` at T = 0.3 K (β·span 1627).
The 17-state electronic-only problem fails identically at the same span, so it is a span
effect, not a dimension effect. Rank 2 at zero frequency and rank 3 at nonzero frequency
are unaffected. Selecting the exponential-action backend gives a finite,
limit-independent value (37.25855779 at both `dense_beta_span_limit` = 1000 and 100), and
P0 confirms the two backends agree to 4 × 10⁻¹⁵ where both run.
`invzf_electronuclear_cumulant` now **errors** rather than returning the NaN, because a
silent NaN would propagate straight into a functional coefficient.

**Learned.** The three failures share one location: **nonzero Matsubara frequency.**
Everything static is excellent — the static-derivative identities close to 10⁻⁸–10⁻⁹ on
all 136 states, which is a genuinely new capability. Everything at nonzero frequency is
suspect: NaN under the default backend, a 2.4 × 10⁻³ conjugation error at rank 3, and a
rank-2 tail that decays like ~ωₙ⁻¹·¹ instead of ωₙ⁻². For a Hermitian operator the
leading tail must be ⟨[[X,H],X]⟩/ωₙ², so ωₙ²C₂ must plateau; it does not.

**Consequence for the critical path, stated plainly.** WP2 needs exactly these
frequency-resolved cumulants to build the diagram manifest's coefficients, and its
Matsubara sums converge only as fast as the tail. **The local source oracle is therefore
not yet fit for WP2**, and this is a blocker, not a nuisance. Resolving it (is the tail
anomaly the engine's large-frequency accuracy or a genuine property? is the rank-3
conjugation error the same root cause?) is the next packet on the critical path, ahead of
any manifest work. Note that P2/P3 passing means the fault is *not* in the wiring, the
source convention, or the spectrum.

**Two errors of mine that the gate caught, recorded because they are the reason to run
gates at all.** (i) The P3 normalisation: I divided C₃(0,0,0) by β and it failed by a
factor of exactly β = 38.68 — with the engine's convention the fully static rank-3
cumulant *is* d²m/dh_z² directly. (ii) The first `invzf_electronuclear_cumulant` read
`cum.value`, a field the engine does not have; it returns `.connected`.

**Checkpoint.** `ckpt-04-wp1-oracle`.

### E7 — WP0: contract frozen; the two independent coupling constructors disagree at finite q

**Attempted.** `invzp_convg_fix.md` WP0 and its gate, "two independent constructors
reproduce all frozen microscopic anchors and the same serialized coupling set". Delivered
`invz_functional/invzf_wp0_contract.m` (versioned freeze of the ion anchors, unit and
Boltzmann constants read from `invz_const`, the Matsubara convention read from
`invz_matsubara`, the q-grid/Fourier/self-site conventions, and exact-byte digests of the
coupling set under each dipolar backend) and `docs/execution/invzp_exec_wp0_gate.m`.

**Outcome.**
- **G0a reproducibility PASS** — repeated construction is bit-identical for each backend
  separately (digests `499922e6c9f8…` brute-force, `87489dd5d8a2…` Ewald).
- **G0b independent constructors at finite q — FAIL.** Over all 16384 (q, branch) pairs:
  max absolute difference 1.136 × 10⁻⁴ meV, max relative difference 12.68%.
- **G0c Γ-point difference PASS** — the backends differ at Γ by the documented Lorentz
  term to within 0.4% of it (lorz = 1.220 × 10⁻³ meV; Jcc0 differs by 2.77 × 10⁻⁶ meV).
- **G0d limit order RECORDED** (not tested): stationary equations at nonzero
  symmetry-breaking source → thermodynamic limit at fixed source and declared dipolar
  shape convention → source to zero; q → 0⁺ and ω → 0 sit inside the finite-volume
  sequence and never substitute for it.

**G0b resolved (dpRng ladder, against the single Ewald reference).**

| dpRng | max abs diff (meV) | max rel diff | RMS diff (meV) |
|---|---|---|---|
| 30 (production default) | 1.136 x 10^-4 | 0.1268 | 3.32 x 10^-6 |
| 45 | 2.113 x 10^-5 | 0.0874 | 9.82 x 10^-7 |
| 60 | 2.894 x 10^-5 | 0.0614 | 8.81 x 10^-7 |
| 90 | 1.273 x 10^-5 | 0.0254 | 4.00 x 10^-7 |

Verdict **(a): convergence toward Ewald**, i.e. brute-force real-space truncation, not a
convention mismatch. RMS and max-relative fall monotonically (8.3x and 5.0x over a 3x
cutoff increase); max-absolute falls overall by 8.9x but non-monotonically at 45 -> 60,
which is expected of a max taken over 16384 modes as the truncation shell moves between
them. Fitted decay max|diff| ~ dpRng^(-1.80).

So the two constructors do agree in the converged limit, and the correct statement of G0b
is: **the production dpRng = 30 coupling set carries a real-space truncation error of
1.14 x 10^-4 meV max / 3.32 x 10^-6 meV RMS against the convergent Ewald sum.** WP0's gate
as literally written -- "the same serialized coupling set" -- is therefore **not passed at
the production default**; it is passed only as a limit statement. That error had not been
measured before, and it is now a declared entry in the frozen contract rather than an
unknown.

**Superseded reading.** G0b's failure was not interpretable at first: the brute-force backend is a
real-space lattice sum truncated at `dpRng = 30` while Ewald is convergent, so the
disagreement may be brute-force truncation error rather than a convention mismatch. A
`dpRng` ladder (30/45/60/90) against the single Ewald reference is running to decide
between monotone convergence and a plateau. The ladder above settled it in favour of the
former.

**Learned.** The production coupling set is *not* an
exactly-defined object at the 10⁻⁴ meV level today — the two constructors the repo
provides differ on it by ~1e-4 meV at the production cutoff, and nothing previously
measured that. That is precisely the
sort of unfrozen anchor WP0 exists to expose before a functional is built on top of it.

**Checkpoint.** `ckpt-04-wp1-oracle` (same commit).

### E8 — S3: the 4.05 T twins are genuinely distinct, and there are far more than two

**Attempted.** Plan §4's first bounded diagnostic (blind WP1c): Newton-polish the 4.05 T
coexisting roots to tight backward error with conditioning estimates, to decide whether
`invzp_convg_diagnosis.md` §6.2's two A–D-accepted states are two distinct high-accuracy
zeros or a pseudo-root pair that only passes because the acceptance tolerance is loose.
Harness: `docs/execution/invzp_exec_s3_twin_polish.m`, using the existing node context /
Newton / enumeration kernels unmodified.

**Outcome — verified facts.**

1. **They are distinct.** Both roots polish cleanly:

   | root | r | polished residual_inf | rcond | pole margin |
   |---|---|---|---|---|
   | A | 1.0876581 | 1.94 × 10⁻¹⁵ (at tol_outer 1e-14) | 5.19 × 10⁻⁶ | 0.008 |
   | B | 0.90653911 | 6.71 × 10⁻¹³ (tightest accepted, 1e-12) | 1.42 × 10⁻⁷ | 0.0016 |

   Separation ‖u_A − u_B‖∞ = 2.40 × 10⁻², |r_A − r_B| = 0.181. The separation exceeds the
   worst polished backward error by **3.57 × 10¹⁰** and the A-block acceptance tolerance by
   2.40 × 10⁶. Verdict: *distinct*, in the required form.

2. **The multiplicity is far larger than two.** The coexistence scan (26 h values × 27
   seeds) finds **6 to 11 distinct accepted root clusters at essentially every h**, e.g.
   11 clusters at h = 0.0055 with r spanning 0.891–1.092, 9 at h = 0.0035 spanning
   0.822–1.143. The diagnosis's "two states pass the A–D contract" was a two-seed view of
   a much larger set.

**Learned.**

- The A–D contract does not select. It admits a whole family, and the family is dense
  enough in r that no local smoothness, proximity, or continuity criterion can separate
  it — which retires the last version of the "choose the smoothest root" idea more
  decisively than §6.2's two-root counterexample did.
- Polishing to 10⁻¹⁵ backward error is **not** certification. No interval-Newton /
  Krawczyk enclosure was computed, so the correct statement is that these are two
  numerically distinct high-accuracy zeros of the residual, not that the set of zeros is
  known or complete. `invzp_convg_fix.md` §7.2 remains the only route to a completeness
  claim.
- Caveat on the count: "distinct cluster" is whatever `invzp_enumerate_ordered_roots`'s
  clustering metric says, so some clusters may be near-duplicates. The *polished pair* is
  not subject to that caveat — its separation is six orders above the acceptance
  tolerance.
- The r values found (1.088 / 0.907) are **not** the pair recorded in §6.2 (0.768 /
  0.822); r ≈ 0.822 does appear in the scan at h = 0.0035. The recorded pair was one
  sample of the family, not a distinguished feature.

**Checkpoint.** `ckpt-05-s3-twins`.

---

## Rejected branches (do not retry without new contradicting evidence)

| Branch | Rejecting evidence |
|---|---|
| Defactoring the static q-average as a fix for the 3.825 T masks | identical failure sets node-for-node in all three masked configurations (E3); the pole crossings are outer-loop, not inner-closure, events |
| `signed_aitken1` as a fix for the 3.825 T masks | no change at `mix_outer = 0.40` (E1); the binding node has no clean alternating mode (lag-2…8 autocorrelation \|ρ\| < 0.07) |
| Reading the 3.825 T mask as a resummed-denominator domain rejection | `medium_status` is `not_applicable` at every node in every configuration; every failure is `max_iter` with `resid_B ~ 1e-11` (E1) |
| Larger iteration budget alone | `max_outer` 1000 → 5000 at `mix_outer = 0.40` changes nothing (E1) |

---

## S4 — closed-loop path-dependence test: SPECIFICATION, not yet executed

Design decisions taken by the user 2026-07-29. Not started; this section is the brief.

**Plane.** Rectangular loop in (h, Bx) at fixed T = 0.1 K:
(h1,B1) -> (h2,B1) -> (h2,B2) -> (h1,B2) -> (h1,B1), traversed in both directions.

**Location.** Bx in the certified 4.60-4.90 T window (revised by the user 2026-07-29 after
the 'not constructible' risk was raised; the earlier choice was 3.0-3.8 T). A loop requires
EVERY point on all four legs to be an accepted node, and E1 showed many nodes in the
3.0-3.8 T band fail under some configurations. Step 1 below is still reconnaissance -- the
rectangle is chosen from the converged sub-region rather than assumed -- but in this window
it should succeed. Consequence to state in the writeup: a loop that closes HERE says
nothing about the hard band; only a loop that FAILS here would be decisive.

**Endpoint test.** ||delta u||_inf < 1e-6 on the full u = [Sigma; K0]. This is six orders
below the 4.05 T inter-root separation (2.4e-2, E8), so it cannot mistake two distinct
roots for one. A loop failing this is reported as `branch transition`, NEVER as evidence
of non-integrability -- that distinction is the whole point of the test.

**Budget.** ~25 points per leg, ~100 node solves per loop direction.

**Steps.**
1. Reconnaissance: at T = 0.1 K, for Bx in {3.0, 3.2, 3.4, 3.6, 3.8} and a coarse h grid,
   record which (Bx, h) nodes reach an accepted A-D state, and at which mix_outer. Use
   invzp_exec_s1_failure_census as the pattern. Output: the converged region's boundary.
2. Choose (h1, h2, B1, B2) as the largest rectangle wholly inside that region. If no
   rectangle exists, REPORT THAT and stop -- do not shrink the loop until it fits by
   excluding the hard points, which would test only the easy region while claiming to test
   the hard one.
3. Walk each leg with warm continuation, recording at every point the complete state
   (m, Sigma, K0), r, the accepted/failed verdict, and any change in which root cluster
   the state belongs to (a branch transition).
4. Compute the loop integral of m dh along the h legs plus the corresponding Bx-leg terms,
   in both traversal directions, with the quadrature error estimated by halving the point
   count.
5. Verdict, in exactly one of these forms:
   - `loop closes`: |loop integral| is within its quadrature error AND the endpoint test
     passes. Proves NOTHING about global integrability -- it is a necessary condition only,
     and cross-sheet constants remain unfixed.
   - `loop does not close`: |loop integral| exceeds its quadrature error with the endpoint
     test passing. Falsifies a common potential in these coordinates.
   - `branch transition`: the endpoint test fails, or a branch transition was recorded on
     any leg. The loop integral is then uninterpretable and must not be quoted as evidence
     either way.
   - `not constructible`: no rectangle of converged nodes exists in the window.

**Do not** loosen the endpoint tolerance, substitute a bare state for a masked one, or
report a loop integral computed across a branch transition.

---

## PRIORITY REVISION 2026-07-29 (user): all-field susceptibility is the goal

The user's primary objective is to **compute chi at every magnetic field**, not to settle
the epistemics first. Strict zero field may be skipped if the double degeneracy is
problematic there. This reorders the remaining work.

**New top packet -- S9: field-wide damping-ladder census.** E1 established that the
3.825 T mask is outer-map iteration dynamics, not a domain failure, and that
mix_outer = 0.30 closes that column completely (0/43 node evaluations fail, accepted A-D
at every node -- 43 not 34 because bisection and the root phase run to completion only
there; the E1 table above is the authority for this count). It is NOT established that
this generalises. The decisive measurement:

  For Bx over the full range of interest (e.g. 3.0-6.0 T at the driver's 0.025 T step,
  plus a coarse low-field pass), run the column at a LADDER of mix_outer values
  {0.7, 0.5, 0.4, 0.3, 0.2, 0.15} with max_outer 1000, and record for each field the
  first ladder rung that closes the column, or that none does.

  Output: the fraction of the field axis recoverable by damping alone, and the residual
  set of fields that no rung closes. Reuse invzp_exec_s1_failure_census /
  invzp_exec_s1_config_sweep unchanged -- this is a measurement, not new machinery.

**If most columns close:** wire an opt-in adaptive-damping mode (try rungs until the
column closes, record which rung was used) into the ordered path, default-off, and expose
it in invz_run_spectra.m. That gives all-field chi spectra now.

**The honest label that must ride with it.** A column closed by a damping ladder yields a
SPECTRUM at an accepted state; it does not yield a certified equilibrium phase label. E4
showed a 1e-12 arithmetic perturbation moving nodes to different A-D-accepted states, and
E8 found 6-11 accepted states at one h. So such columns must be reported as
`converged-under-damping-ladder, branch-not-certified` -- which is exactly the interim
vocabulary of plan SS2.7/SS4, and does not breach the promotion boundary as long as the
phase LABEL is not claimed. Producing chi where it was previously masked is a legitimate
deliverable; claiming it is the equilibrium branch is not.

**Resulting order:** S9 census -> (conditional) adaptive-damping mode -> S4 loop test ->
WP1 nonzero-frequency blocker -> WP2 manifest. S4 drops below S9 because it is an
epistemic check, while S9 is on the path to the user's actual deliverable.

## E9 (2026-07-29) -- S10: deep-ordered root cause is a finite-q SOFT MODE, not a solver defect

Trigger: S9 found 1.0 T failing every rung of the damping ladder (962.7 s), with the
failure count FLAT in mix_outer (10-12 of 34 over 0.15-0.50) where 3.825 T went 22 -> 0.

**Verified facts.**
1. No local linear regime exists. The tail ratio of the raw map step never settles
   (IQR 0.9-2.5 on all ten commonly-failing nodes), so the dominant Jacobian eigenvalue
   is not measurable and neither "lambda > 1 unstable" nor "lambda < 1 slow" applies.
   (invzp_exec_s10_outer_eigenvalue.m; out/s10_B1_summary.mat)
2. The LOCAL Gstat pole is on the PATH, not the solution. Accepted nodes sit at
   |d0| in [0.476, 1.61] with none below 0.1, identically across all four alpha, while
   the path reaches |d0| ~ 1.8e5 and changes sign 912-1410 times per 1000 iterations.
   (invzp_exec_s10_pole_locality.m; out/s10_pole_locality.mat)
3. **The failures correspond ONE-TO-ONE with a near-vanishing q-resolved denominator.**
   At 1.0 T, 11 of 34 nodes have a FINITE Dq_absmin and 11 fail -- the same 11. No node
   with finite Dq_absmin is accepted; no failing node lacks it. Six have Dq_min < 0, i.e.
   D(q) = 1 + (J(q) - K0)*Gstat CROSSES ZERO in the Brillouin zone (min |D(q)| down to
   8.19e-4). Every node of the converging 3.825 T column has Dq_absmin = NaN.
4. 1.0 T is NOT a globally failing column: 24 of 34 nodes converge in 50-83 outer
   iterations. Failures are two populations -- a contiguous block h in [0, 1.28e-4]
   (including h = 0 exactly, the double degeneracy), and isolated nodes 19/22/23/24
   sitting between converging neighbours (the E8 multiplicity signature). The column
   returns hstar = NaN because the root bracket needs the small-h end that is missing,
   not because the field is unsolvable.

**Interpretation (inference, not yet independently confirmed).** At low transverse field
and small longitudinal source the 1/z effective medium develops a divergent static
response at a finite wavevector. The static closure q-averages 1/D(q); where D(q) changes
sign that average has a genuine pole and no admissible fixed point exists. The solver is
not failing to find a state -- it is correctly declining to cross a pole, which is the
deliberate behaviour described in invzp_convg_diagnosis.md Sec 7.2.

**REJECTED BRANCHES (both refuted by pre-stated criteria, before any solver was written).**

| branch | prediction | measured | verdict |
|---|---|---|---|
| pole-regular arithmetic (`closure_coordinate='defactored'`, `stable_form`) clears the mask | fewer failures at 1.0 T | mix=0.30: **11/34 -> 11/34, identical node-for-node** [1 2 3 4 5 6 7 19 22 23 24]; mix=0.15: 10/34 both | REFUTED. It regularises the LOCAL d0 pole, but only ~2% of iterations sit within |d0| < 0.1 -- the iterate leaps across that pole rather than landing on it, and the binding singularity is in D(q), a different denominator. |
| backtracking line search / step limiter would prevent the excursion | outward steps preferentially worsen the residual | frac worsening = 0.469 overall vs **0.467** restricted to steps landing at |d0| > 10; corr(d|d0|, d|resid|) = **-0.002** | REFUTED. No outward preference and no correlation: a random wander, not a systematic overshoot. A line search would reject ~47% of steps at random. |

**Lesson.** Both rejected fixes were algorithmic responses to what is a structural feature.
The `Dq_absmin` column was present in the S1 census table from the start and would have
identified this immediately; the damping-ladder framing (S9) directed attention to
iteration dynamics and delayed reading it. When a failure count is FLAT in the solver's
own knob, that is early evidence the obstacle is not in the solver.

**Consequence for the plan.** An adaptive-damping production mode is NOT justified and is
withdrawn as a deliverable: it cannot close columns whose obstruction is a sign change in
D(q). The open question is now physical -- whether the finite-q instability is (a) a real
competing ordering wavevector at low field, (b) an artefact of the 1/z truncation's
q-dependence, or (c) the documented 1.5 T fold region -- and it must be settled before any
low-field chi is trusted, however cleanly a solver might be made to converge there.
