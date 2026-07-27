# The 1/z ordered-solver convergence fix — attempt record

**Author:** Claude · **Date:** 2026-07-27 (rev. 2 — corrected after external review, see §12)
**Branch:** `invzp-stage2c-diagnostic`
**Scope:** from the closure of the Ewald opt-in integration (`2ee310b`, 2026-07-25) to the Gate-0
verdict (`6937e9d`, 2026-07-27). Everything before that — the Stage-2c lattice diagnostic and the
Ewald primitive itself — is out of scope and recorded in `docs/INVZ-DEVELOPMENT-RECORD.md`.

**Outcome: the preregistered candidate FAILED Gate 0.** The run stopped at diagnosis, as its own
preregistration required. Nothing was regularised, no tolerance was widened, no scheme was switched,
and the production default was never flipped.

---

## 0. One-paragraph summary

The 1/z ordered spectra panel masks below `Bc_1z ≈ 3.025 T`: the Jensen H_MF ordered solve returns no
root, so the phase label falls through to "solver failed". After the Ewald work established that this
is **not** a lattice/quadrature defect, the hypothesis under test became: the *resummed* static medium
carries a q-denominator pole that dies just below the uniform instability, and truncating that medium
to strict `O(1/z)` — `K0 = Jbar - mu2*Gref`, one-shot, no denominator feedback — removes the pole and
lets the ordered root exist. A 19-task plan built that closure as an additive, opt-in scheme (default
unchanged), preregistered a five-clause pass/fail predicate before any strict run, and then measured
it. On a synthetic coupling multiset it worked and moved the spontaneous root by +55%. On the real
16384-point BZ multiset it failed three of the five clauses. The truncation is not controlled where it
matters, and the failure is at **low** field, deep in the ordered phase — not near the boundary where
the resummed scheme dies.

---

## 1. Where this started: what Ewald fixed, and what it did not

Ewald was pursued because the brute-force dipolar lattice sum (`MF_dipole`) is conditionally
convergent, making the couplings offset-sensitive. That is a real defect and the Ewald primitive
repairs it. It is **not** the masking cause. Three probes settled this (2026-07-25):

| probe | measurement | conclusion |
|---|---|---|
| 1 — coupling comparison, 16³ | `Jcc0` reldiff `4.3e-4`; sorted-multiset `‖ew−bf‖₂/‖bf‖₂ = 1.2e-3`; `J_max/J_min` identical to `~1e-4` | Ewald and brute-force couplings are essentially the same input |
| 2 — `invz_hmf_ordered` at 1 T and 3 T, **both** backends | `hstar = NaN`, `status = node_failed`, `slope0 = NaN` (the `h=0` predictor itself fails), `D_uni < 0` at 30–79 % of 33 nodes (min `−1.3 … −47`), `Dq_neg_count →` the full ~16384 | identical failure on both backends; whole-BZ static `D_q ≤ 0` crossing |
| 3 — grid sweep, Ewald, 1 T | `node_failed` at `N = 8, 12, 16, 24` | grid-independent; not a quadrature-density artifact |

The outer Picard map limit-cycles (`term_reason ∈ {converged, max_iter}` alternating) on
unstable-branch states. **The masking is solver-side, in the ordered static medium.** Ewald was
completed as opt-in integration (Step 5, suite 394/0/23, default stays `bruteforce`) and the effort
pivoted here.

*Note for provenance:* everything below ran on the **legacy brute-force** couplings, not Ewald. The
Ewald backend is integrated but not default, and the Step-6/7 physics-anchor and default-flip work was
never executed. That plan (`docs/superpowers/plans/2026-07-24-ewald-step5-integration.md`) is written
but unexecuted beyond Step 5.

---

## 2. The hypothesis, and why it was plausible

Full derivation: `docs/superpowers/specs/2026-07-25-invzp-ordered-solver-static-medium-design.md` §A–§B.

**Jensen's condition is kept verbatim.** Only the `ω_n = 0` medium changes. With
`h0(hmf) = ∫₀^hmf r dh'`, `r = G0bare/Gtil0`, `F = h0 − J0eff·m`, and J 2.31 (`dm/dh = −G0bare`):

```
F'(h) = r(h) + J0eff·G0bare(h) =: crit(h)        [dimensionless mass]
F(h)  = ∫₀^h crit dh'
```

so `crit(h*) > 0` is exactly the Landau local-minimum condition, and `crit(0)` is the existing
`prof.slope0`. This gave the ordered leg a stability criterion it previously lacked.

**Both legs truncate identically.** Expanding the ordered closure (`Gq = Gs/(1+(J−K0)Gs)`,
`mean_q Gq = Gs`) and the PM closed form (`invz_emt_scalar`) independently in the local dressed
propagator gives the *same function*:

```
K0 = Jbar − mu2·G + mu3·G² + (2·mu2² − mu4)·G³ + (mu5 − 5·mu2·mu3)·G⁴ + O(G⁵)
```

verified term-by-term to `O(G⁴)`, difference exactly 0. No q-denominator survives formal expansion.
Under `mu2 ~ 1/z`, the strict `O(1/z)` medium is `K0 = Jbar − mu2·Gref` with `Gref` **not** fed back.
Using the same truncation on both sides preserves the exact `m → 0` PM/ordered identity — which
matters because an `O(1/z²)` mismatch at a continuous boundary, where the target mass is exactly zero,
can become the *leading* residual.

**The a-priori breakdown ladder** (`|G0|` = bare local static weight, `D ≈ 1`):

| `\|G0\|` | event | kind |
|---|---|---|
| **155.7** | `D_uni = 1 + (J0eff − K0)G ≤ 0` — physical uniform instability | real physics, keep as observable |
| **167.1** | first q-average pole, `D + J_max·G0 = 0` — **the resummed closure dies** | arithmetic |
| **208.2** | self-consistent `mu2` branch point | arithmetic |
| **416.2** | one-shot `mu2` pole, `1 + Jbar·G0 − mu2·G0² = 0` | arithmetic |

The uniform instability precedes the resummed pole by ~7 %. That 7 % gap is the whole reason
"below `Bc` ⇒ node dies", and the hierarchy is monotone in *how much beyond-order denominator feedback
is retained*. Removing the feedback should buy a factor `2.49` of domain. **That was the argument, and
it is still correct as far as it goes** — it is a statement about where the *arithmetic pole* sits, not
about whether the truncated series is accurate there.

**Two conventions were carried:** `strict_1z_dyson_ref` (`Gref = G0bare/(1+Sigma0)`, the frozen
selection) and `strict_1z_bare_ref` (`Gref = G0bare`, a systematic comparator only).

**Preregistration.** `docs/invzp_strict_medium_prereg.md`, frozen 2026-07-25 **before any strict run**,
fixing: the coupling fixture and its digest, the field set, `nH = [33 65 129]`, every tolerance
(`crit_tol = 1e-6`, `omit_promote = 0.10`, `ref_margin = 1e-6`, `pole_cont_tol = 1e-3`,
`K_atol = 1e-14`, `K_rtol = 1e-12`, `I_atol = 1e-10`), the five-clause Gate-0 predicate, and — §9 — a
*blind* Stage-4 predicate written before Stage 4 could be planned. Amendment rule: a new dated
section, never an edit in place.

---

## 3. What was built

37 commits, 53 files, +6744 / −320. Fifteen new production files, all additive:

```
invz_common/   invz_coupling_moments.m        mu_n = mean((J-Jbar).^n)  [NOT var(): N-1 is wrong here]
               invz_medium_moment_closure.m   K0 = Jbar - mu2*Gref, + omitted-order diagnostics
               invz_static_medium_reference.m Gref construction, both conventions
               invz_check_static_medium.m     one scheme-resolution authority, stamped into both legs
               invz_pm_verdict.m              three-way dispatcher: PM / ordered / genuinely-failed
               invz_boundary_interval.m
               invz_is_recoverable_solver_error.m  whitelist = {orderedPhase, degenerateDoublet} exactly
               invz_try_solver_call.m
               invz_exact_numeric_digest.m
invz_projected/invz_hmf_grid.m                geometric h-grid, hgrid = hmax*ratio.^((nH-1):-1:0)
               invz_hmf_status.m              status precedence (degenerate binds over everything)
               invz_static_domain_scan.m      prospective domain scanner (Gate 0 stage a)
               invz_pole_continuity.m         actual-path pole-jump measurement (G17)
               invz_gate0_report.m            the Gate-0 driver
               invz_gate0_aggregate.m         the (a)-(e) predicate, independently forceable
```

Constraints honoured throughout: the production default stayed `static_medium = 'resummed'` and
**bit-identical** (gate G9); the strict path is reachable only by explicit opt-in; `invz:*` error
identifiers only; the recoverable whitelist was never widened; the spec and the frozen preregistration
were never edited.

Two structural additions worth naming, because they outlived the hypothesis:

- **The two-tier verdict.** `res.stability` is computed for *every* node but gates acceptance for
  **none**. Diagnosis and acceptance were deliberately decoupled, which is why the failure below is
  fully characterised rather than just "it didn't converge".
- **Coverage accounting.** Every node lands in exactly one of four buckets plus `unrecognized`, and
  `n_accounted == numel(trc.nodes)` is asserted. This is clause (b), and it is the one clause that
  never fired anywhere — the accounting itself is sound.

---

## 4. What was measured, in the order it was measured

### 4.1 The synthetic fixture said yes (Task 12, mid-plan)

Same fixture, both schemes, reproduced independently by the controller:

| | `resummed` | `strict_1z_dyson_ref` |
|---|---|---|
| status | ok | ok |
| `hstar` | `0.0085473576` | `0.0132482595` **(+55 %)** |
| `slope0` | `−0.0633068` | `−0.1634994` |
| `crit_star` | `0.0868790` | `0.1460789` (both > 0) |
| `∫Sigma0 dh` | `−1.67716e-3` | `−1.98539e-3` |
| `∫(r−1) dh` | `−1.73252e-4` | `−4.93252e-4` |
| ratio `(r−1)/Sigma0` | `0.1033` | `0.2484` |
| `max(omit_max)` | n/a | **`0.0350`** — comfortably under `omit_promote = 0.10` |

This looked like success. It also produced one durable, correct result: **using `∫Sigma0 dh` as a
proxy for `∫(r−1) dh` overstates the correction by ~10× (resummed) / ~4× (strict).** The spec's
insistence on measuring both integrals was load-bearing, not defensive.

**The lesson, in hindsight:** the `0.0350` was an artifact of the 24-point synthetic multiset. It was
read as evidence the truncation was controlled. It was evidence about the fixture.

### 4.2 The real multiset said no (Tasks 17–18)

Fixture, verified by digest on every invocation before any solve:

```
grid [16 16 16], dpRng 30, dipole 'bruteforce', cache false     (55.9 s build)
digest = ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17   MATCH vs frozen
n = 16384 (4096 q × 4 branches)
J0eff = 6.424435656e-3 meV        Jbar = 1.207664433e-4
mu2   = 5.482637653e-6            mu4  = 2.3894·mu2²        mu3 = −3.42228e-11 (skew −0.0027)
J_max = 5.985139e-3               J_min = −6.763100e-3
```

**Gate 0 was run at full preregistered coverage:** 8/8 ordered fields `[0.05 0.25 0.5 1 2 2.5 2.9 3.0] T`
× 3 `nH` = 24/24 rows, both PM controls `[3.1 3.5] T`, plus the exact `B = 0` hard-domain control.
Nothing was stopped early. (Task 17's earlier G11 covered 1 of 8 ordered fields and 0 PM controls; the
gap is stated here because an unqualified "Gate 0 failed" from Task 17 alone would have implied
coverage it did not have.)

```
fail_a = TRUE    12 of 24 rows — all three nH at B ∈ {0.05, 0.25, 0.5, 1} T
fail_b = FALSE   coverage identity held on every row; n_unrec = 0 everywhere
fail_c = TRUE    max(omit_max) exceeds 0.10 on ALL 24 rows; smallest = 0.17782 (3 T, nH=129)
fail_d = FALSE   the one genuine local Gstat crossing (2 T) resolves 'ok' with a decreasing jump
fail_e = TRUE    the same 12 rows, AND independently both PM controls

rep.pass = ~(a || b || c || d || e) = FALSE
```

## VERDICT: FAIL — clauses (a), (c), (e) fired.

---

## 5. Five adverse findings

Three of these are Gate-0 clause failures; G5 and G7 are separate numerical-quality measurements with
no pass/fail role in `rep.pass` (G7 has no preregistered threshold at all). They are grouped here
because they are all evidence about the candidate, not because they carry equal weight.

### 5.1 The reference denominator leaves its domain — at LOW field

`ref_denom_nonpositive` at `B ∈ {0.05, 0.25, 0.5, 1} T`, all three `nH`. At `B ∈ {2, 2.5, 2.9, 3.0} T`
the solve returns `status = 'ok'` with a finite root at every `nH`.

**This inverts the expected shape.** The hypothesis predicted trouble near the instability, where the
resummed pole sits. Instead strict is *algebraically solvable and in-domain at the sampled
boundary-adjacent fields*, and breaks *deep in the ordered phase*, where the order parameter is large.
(That is a statement about four sampled fields, not a continuity claim about the interval.) At 1 T the
affected region is a bounded interior window — four consecutive nodes with solved `ok` nodes on both
sides:

```
 idx   h              m         medium_status           term_reason
  20   0.0021219622   4.16509   ok                      converged
  21   0.0026332230   4.54081   ref_denom_nonpositive   medium_out_of_domain
  22   0.0032676658   4.81953   ref_denom_nonpositive   medium_out_of_domain
  23   0.0040549699   5.00369   ref_denom_nonpositive   medium_out_of_domain
  24   0.0050319653   5.11630   ref_denom_nonpositive   medium_out_of_domain
  25   0.0062443557   5.18629   ok                      converged
```

Nodes 1–20 rise smoothly `m: 0 → 4.17`; nodes 25–34 run `m: 5.19 → 5.44`; all `ok`.

**Precisely what this establishes** (this wording matters and was corrected after external review):
the current warm-started **and** cold-retry Picard attempts *encounter* `1 + Sigma0 ≤ 0` in a bounded
intermediate-`h` window. It does **not** establish that `1 + Sigma0` passes through zero as a property
of a converged solution, nor that no positive-denominator fixed point exists there. The flagged nodes
are not converged states: `invz_ordered_node_solve.m:226-228` breaks the moment the medium status
leaves `{not_applicable, ok}`, *before* `invz_lambdas` at `:232`; and `:101-105` / `:165-171` retry
once from a cold start. Two independent seeds were tried and both halted. Proving a branch crossing
would need a separate continuation study — legitimate diagnosis toward a *future* candidate, never a
rescue of this one.

### 5.2 The omitted order is not small on the real density

`max(omit_max)` exceeds the frozen `0.10` on **all 24 rows**, including every row that solves cleanly.
The smallest value anywhere is `0.17782` — 78 % over the limit, at `B = 3 T`, `status = 'ok'`, finite
root. Clause (c) is therefore independent of clause (a): the truncation is uncontrolled even where
nothing looks wrong.

| profile | `max(omit_max)` |
|---|---|
| synthetic fixture, strict (Task 12) | `0.0350` |
| real, `strict_1z_dyson_ref`, 1 T | `1.20077` (profile) / `1.201054` (full ledger incl. predictor) |
| real, `strict_1z_bare_ref`, 1 T | `16.4358` |

**Mechanism.** `omit_cubic = |2·mu2² − mu4|·Gref²/mu2`. With `mu4 = 2.3894·mu2²` this is
`0.3894·mu2·Gref²` — it grows exactly where `|Gref|` is large, i.e. right next to the near-singular
nodes. The first omitted coefficient being "less than one" was never a small-error guarantee, and the
spec said so; the synthetic fixture simply never probed `|Gref|` large enough to show it.

Same conclusion from the path integrals, reported like-for-like rather than as one headline number:

| profile | `∫(r−1)dh / ∫Sigma0 dh` | vs real bare-ref |
|---|---|---|
| bare-ref, real couplings | `25.4842` | — |
| strict, synthetic fixture | `0.248441` | real is **102.6×** |
| resummed, synthetic fixture | `0.103301` | real is **246.7×** |

### 5.3 The PM leg converges to a negative mass

```
B = 3.1 T   converged = 1   medium_status = 'ok'   crit = −0.884482
B = 3.5 T   converged = 1   medium_status = 'ok'   crit = −0.528219
```

Both are frozen PM control fields, and prereg §3(e) requires each to return "a converged finite
**positive**-mass PM state". Both violate that control. Cross-checked by calling `invz_solve_point`
directly, outside all Gate-0 machinery — not a driver artifact.

**What this does and does not say.** The defensible statement is that the strict candidate reports a
converged, in-domain, **negative** PM mass at both control fields, implying it predicts the PM state
remains unstable through 3.5 T. It is *not* established that the sign is "wrong relative to the strict
candidate's own phase boundary": `Bc_1z = 3.025 T` is the legacy/resummed boundary, and the strict
candidate's own boundary was never independently determined. For reference only, `resummed` at 3.5 T
returns `crit = +0.374153` — but it does **not converge** there, so that number cannot serve as an
independent validation of either sign.

This is a **second, independent** failure route. It has nothing to do with the ordered-leg denominator
window in §5.1, and it fires clause (e) on its own.

### 5.4 The path integrals miss their frozen tolerance (G5)

G5 is a separate numerical-quality gate, not one of the (a)–(e) clauses, and is reported rather than
folded into the Boolean. **It fails its preregistered 65→129 acceptance tolerance on every field**,
including all four that return `status = 'ok'`.

**But the refinement sequence is convergent.** On the four fields with valid integrals the successive
differences fall by almost exactly the second-order factor:

| B (T) | `d33→65` / `d65→129`, `∫Sigma0` | same, `∫(r−1)` |
|---|---|---|
| 2 | `1.55618e-05 / 3.87346e-06` = **4.018** | `2.94724e-05 / 7.34317e-06` = **4.014** |
| 2.5 | **4.017** | **4.015** |
| 2.9 | **4.014** | **4.013** |
| 3 | **4.012** | **4.012** |

All eight ratios lie in `[4.0118, 4.0175]` against the exact second-order value of 4.

So G5 is an adverse result about *how tight the frozen tolerance is relative to `nH = 129`* — the same
structural problem as G1 in §7 — and **not** evidence that the Jensen path integral fails to converge.
An earlier revision of this record said the integrals "do not converge under grid refinement". That
was wrong and is withdrawn.

The low-field `∫(r−1)dh = NaN` at all three `nH` is likewise **a consequence of the reference-domain
failure in §5.1** — there is no valid integrand to integrate — not an observed quadrature result. It
must not be read as one.

Related: prereg §4's binding caution — that the ~0.3 % PM boundary shift (`Bc_1z 3.025` vs bare
`3.033`) bounds *neither* path integral deep in the ordered phase — stands. At `B = 0.05 T`,
`∫Sigma0 dh = −10.75`, and at low field `∫(r−1)dh` is not even finite. (An earlier revision called
that integral "~3500× larger than a 0.3 % shift could bound". That divides a dimensional integral by a
dimensionless percentage and is undefined; withdrawn. The correct statement is simply that a
boundary-shift observation supplies **no** bound on either ordered-path integral, which is what the
prereg said and what the measurements are consistent with.)

### 5.5 The `ω_n = 0` scheme jump is not small (G7) — at primitive level

The strict truncation applies at `ω_n = 0` only, so the static slot differs between schemes while
`ω_1 → 0` as `T` falls — a scheme discontinuity across adjacent Matsubara grid points. Measured at
`B = 6 T`, `Ecut = 40 meV`:

```
T=0.05 K   |K1s−K1r| = 1.492e-4   |K2−K1| = 2.524e-4   ratio 0.591
T=0.10 K   |K1s−K1r| = 1.307e-4   |K2−K1| = 2.460e-4   ratio 0.531
T=0.31 K   |K1s−K1r| = 8.617e-5   |K2−K1| = 1.905e-4   ratio 0.452
T=1.00 K   |K1s−K1r| = 5.955e-5   |K2−K1| = 3.434e-4   ratio 0.173
```

**Read this precisely.** The measurement (`.superpowers/sdd/task18_g7.m`) calls `invz_emt_scalar`
directly with `Sigma = zeros(size(wn))` — an *undressed, primitive-level* comparison, not a converged
solve. The denominator `|K(2) − K(1)|` is `K` at Matsubara index 2 versus index 1 **on the resummed
scheme**, i.e. **Matsubara-frequency dispersion of `K`** — it is *not* the coupling's q-dispersion.

So: the static-only truncation produces a **17–59 % jump relative to the measured Matsubara-frequency
dispersion in this primitive-level comparison.** Propagation of that jump into `invz_lambdas` and
downstream spectral quantities was **not measured** and remains open. "Assume it is negligible" is
nevertheless no longer available as a default.

> **Erratum affecting a committed document.** `docs/invzp_strict_medium_gate0_report.md` §4 describes
> this same ratio as "not negligible relative to the coupling's own q-dependence". That label is
> wrong for the reason just given — the denominator is a frequency difference, not a q-difference.
> G7 is explicitly non-gating (prereg §6 sets no threshold), so the Gate-0 verdict is unaffected, but
> the sentence should be corrected.

---

## 6. What held, and is worth keeping

- **Coverage accounting (clause b) never fired anywhere.** Every node accounted for, every predictor
  phase present, every finite `hstar` backed by a value-agreeing `root` ledger entry, no dropped chunk,
  on all 24 rows.

  **Scope limit, and a latent defect found in review.** That completeness is verified for the
  **ordered ledger only**. `invz_gate0_aggregate(ordered_rows, pm_rows, expected_fields, expected_nH)`
  drives its coverage loop from `expected_fields × expected_nH` against `ordered_rows`; it accepts
  **no expected-PM-field list** and simply iterates whatever `pm_rows` it is handed. An empty
  `pm_rows` yields `any(false(1,0)) == false`, so the PM half of clause (e) would be **silently
  skipped** — the `all([])`/`any([])` trap, in my own aggregator. **This run is unaffected**: both PM
  controls were passed in and both fired (e), as recorded in the verdict document §3.2. But before the
  aggregator is reused, it needs an explicit expected-PM-field check that fails closed on missing or
  duplicate controls. Until then, describe the accounting as complete for the ordered ledger, not for
  all Gate-0 inputs.
- **The removable-singularity exemption (clause d) worked as designed.** Spec §1 argued that a *local*
  `Gstat` denominator crossing is not fatal because `r` and `crit` are continuous through it, whereas
  the *reference* denominator `1+Sigma0` has no such cancellation and keeps a hard rule. Measured: the
  only genuine crossing (`B = 2 T`) shows `jump_exceeded` at `nH=33` and resolves `ok` at both 65 and
  129 with the jump *decreasing* (`8.346e-4 → 1.594e-4` for `r`; `4.729e-4 → 5.04e-5` for `crit`),
  well inside `pole_cont_tol = 1e-3`. Distinguishing the two denominator classes was correct.
- **The `m → 0` cross-leg closure identity is exact.** Three of four G2 comparisons are *exactly zero*,
  including `K_ordered = K_PM` **bitwise** — stronger than §7.3 promised, which explicitly disclaimed
  bitwise agreement.
- **G3:** `r = 1 + Sigma0` holds at `RelTol 1e-12` for every `K0 ∈ {0, 1e-3, 5e-3, 0.05}`.
  Characterisation identified `G0el0 = 0` as the algebraic condition producing it — the relevant
  *consequence* of the physical `m → 0` limit, which §7.3 independently names as where the two `Gref`
  expressions coincide.
- **G15 / G15b:** fatal identifiers escape every layer; a domain outcome keeps its category
  (`ref_margin = 1e9` drove all 33 nodes through `ref_denom_small` and produced
  `status = medium_out_of_domain`, confirming the intended path fires).
- **G9 bit-identity:** the default `resummed` path is unchanged.
- **The `B = 0` control behaved exactly as specified**: `degenerate_doublet`, forced by
  `invz_hmf_status`'s binding precedence at the `h=0` predictor node, and excluded from clauses (a)/(e)
  **by construction** — never passed into the aggregator — rather than by the numbers cooperating.

---

## 7. Three gates that were wrong — a separate failure class

Seven of Task 17's eleven gate tests failed. **Three of those failures are defects in the gate
definitions, not in the code under test.** Distinguishing them from the real Gate-0 result took
independent measurement, and getting it wrong in either direction would have been serious: dismissing
a real failure, or "fixing" code that was correct.

### G1 — the identity holds; the frozen gate is incompatible with the frozen grid

All three G1 forms compare **secant/trapezoidal** quantities against **differential** identities, so
the residual is O(h²) panel error, not an identity violation:

| `nH` | G1a `dm/dh + G0bare` | G1b `dF/dh − crit_avg` | G1c `dF/dm − crit/chi_path` |
|---:|---:|---:|---:|
| 33 | `1.8216829803e-02` | `5.8431222186e-03` | `2.2609487299e-04` |
| 65 | `4.8713648324e-03` | `1.4726351883e-03` | `6.4840455639e-05` |
| 129 | `1.2528393091e-03` | `3.6828822520e-04` | `1.7383399960e-05` |

Reduction ratios `3.9678, 3.9986` → empirical order `1.9883, 1.9995`. Second order, confirmed. The
gate's own monotonicity assertions passed; only the absolute `≤ 1e-6` check failed. A pure
`trapz(exp)` control fixture passed, independently confirming the quadrature machinery.

**`1e-6` is not mathematically impossible — it is incompatible with the frozen `nH = 129`.** At second
order it needs roughly `sqrt(368) × 129 ≈ 2500` nodes. Both the threshold (§5) and the grid (§8) are
frozen, so **no dial can be turned without a dated amendment.** That deadlock is self-inflicted and is
the clearest process lesson of the whole run.

**All three G1 forms are provably the same test.** `invz_hmf_ordered.m:397-401` builds
`h0 = cumtrapz(hgrid, r)` and `F = h0 − J0eff·m`; `:539` sets `crit = r + J0eff·G0bare`. Because `h0`
is a cumulative trapezoid of `r`, the secant `diff(h0)./diff(h)` is *exactly* the midpoint average
`r_avg`. Write `a = Δm/Δh`, `g = midpoint-avg(G0bare)`, `ρ = midpoint-avg(r)`, `J = J0eff`. Then
`ΔF/Δh = ρ − J·a` and `crit_avg = ρ + J·g`, so

```
G1a  =  a + g
G1b  =  ΔF/Δh − crit_avg              =  −J·(a + g)                    =  −J · G1a
G1c  =  ΔF/Δm − crit_avg/chi_path     =  (ρ−J·a)/a + (ρ+J·g)/g
                                      =  ρ·(a + g)/(a·g)               =  ρ/(a·g) · G1a
```

Since `chi_path = −G0bare > 0` gives `g ≠ 0`, `a > 0` along the profile, and `ρ ≈ 1` near onset,
**G1b and G1c both vanish exactly when G1a vanishes.** They are rescalings of one identity, not
independent checks. The differing magnitudes and slightly differing reduction ratios come from the
`−J` and `ρ/(a·g)` prefactors varying along the grid, so the max-norm is attained at different nodes.

The `−J·(a+g)` relation was confirmed numerically at `5.7e-15`, but the derivation above shows both
relations are algebra, not measurement. **Prereg §5 registers all three forms as independent; all
three are one gate.** The registered coverage of §5 is therefore substantially narrower than it
appears — a genuinely independent check would have to be a local/thermodynamic test, not another
rearrangement of the same secant identity.

> An earlier revision of this record, and the current amendment draft, claimed only that G1b is
> dependent and proposed retaining G1a **and G1c** as two independent checks. That is wrong: the
> algebra above shows G1c is dependent too. External review supplied the G1c relation; I verified it
> by re-deriving it. The amendment draft has been corrected accordingly.

### G2 — the closure identity is exact; the failing assertion is out of scope

The one failing comparison is `|p.K0_pm0 − ptp.K(1)| = 1.021112e-12` against a gate of `1.6e-14`. But
`K_atol/K_rtol` is the **block-B recomputation** gate, whose own preregistered text reads: *"A one-shot
formula recomputed from the exported state should agree to floating-point reassociation only; this gate
catches mis-wiring, not physics."* It is scoped to one exported state. `p.K0_pm0` and `ptp.K(1)` are
outputs of **two independently terminated outer iterations**. `1.02e-12` absolute is `3.4e-10`
relative — cross-solver iteration residue, which this gate was never designed to bound.

### G13 — no PM-slot leak, but the asserted invariant is not one the code promises

`state.lam` is computed *inside* the outer loop (`invz_ordered_node_solve.m:232`); `state.K(1)` is
overwritten *after* the loop by the post-loop static refresh (`:257`). The test compares a pre-refresh
quantity against a post-refresh one, so the `4.705e-10` residual is **structural, not noise**. The
sentinel's actual purpose is satisfied: the ordered `K(1)` differs from the PM slot by `4.71e-6`
(relative `1.73e-3`), and the counterfactual leaked lambdas sit `2.974e-5` away — about `1.7×10⁶` times
the observed absolute residual `1.772e-11`.

> **A correction I had to make.** I first reported that margin as "63,000×", dividing an **absolute**
> lambda displacement by a **relative** residual. Dimensionally invalid. Caught by external review. The
> corrected figure is larger and strengthens the conclusion — but the original arithmetic was wrong,
> and a wrong number in the right direction is still a wrong number.

**All three require a dated preregistration amendment, not a test edit**, because every relevant
constant is frozen. Draft: `Task-17_prereg_amendment_DRAFT.md`, written to be appended verbatim as a
new §11, leading with what it does *not* change. **Unsigned. Until it is signed, the six tests stay
red and the branch stays red — deliberately.**

---

## 8. Methodological lessons

1. **Synthetic fixtures flatter.** `max(omit_max)` was `0.0350` on 24 synthetic points and `1.20`–`16.4`
   on 16384 real ones — a factor of 34–470. The synthetic fixture had the right *cardinality logic* but
   not the real distribution's tail. Any diagnostic whose value depends on `|Gref|` growing large must
   be validated on the real density before it is treated as evidence of control.
2. **Freezing two dials that trade against each other can deadlock you.** §5 froze `1e-6`; §8 froze
   `nH = 129`. They are incompatible for an O(h²) residual, and nothing in the preregistration process
   caught it because each looked reasonable alone. Preregistration should include a feasibility check
   of each threshold *against the grid it will be evaluated on*.
3. **"Converged" is not "physical" — twice.** First the resummed limit cycle (`term_reason` alternating
   `converged`/`max_iter` on `D_uni < 0` states); then the strict PM leg reporting `converged = 1`,
   `medium_status = 'ok'`, and a negative mass. A solver's own convergence flag is not a physics gate.
   The two-tier verdict (`res.stability` computed everywhere, gating nothing) was the right structure.
4. **Blindness is not health.** `strict_1z_bare_ref` was the only scheme returning `status = 'ok'` at
   1 T — with `omit_max = 16.4358`, 164× the limit, and a root of `6.67e-6`. Its inability to report a
   reference-denominator failure is an *expected limitation of the convention* (`denom ≡ 1`, nothing to
   screen), not a clean bill of health. A comparator that cannot fail a check is not passing it.
5. **Report ratios in consistent units.** This one bit **twice**, and both times in my favour, which is
   exactly why it is dangerous: the G13 "63,000×" margin divided an absolute displacement by a relative
   residual, and the G5 "3500×" divided a dimensional integral by a dimensionless percentage. Neither
   ratio existed. Both were caught in external review, not by me. A number that flatters the argument
   deserves *more* scrutiny than one that does not — check the units before quoting the factor.
6. **Do not background a long MATLAB run and poll its log.** The Task 17 implementer did exactly that,
   died mid-poll, wrote no report, and left a suite log truncated mid-run. Everything had to be
   recovered from leftover scratch files and re-measured. Foreground only, result inline.
7. **Verify claims; trust measurements.** Across this run, every *quoted number* an agent reported was
   correct. Every *assertion* of the form "this branch is unreachable / this is inert / this is safe"
   that I checked turned out wrong. One implementer declared the domain branch "unreachable at every
   fixture I have"; a one-line `ref_margin = 1e9` reached it immediately.
8. **Pre-dispatch brief verification was the highest-yield activity.** Checking each task brief against
   the actual source before dispatching caught, among others, that `nargout` returns `−1` for anonymous
   handles — so the brief's `invz_try_solver_call(@() error(...))` throws `MATLAB:maxlhs` instead of
   returning the identifier. Cheaper to catch before dispatch than after.

---

## 9. What is NOT established

Stated explicitly, because several of these are easy to over-read from the above:

- **Not established:** that no valid positive-denominator ordered fixed point exists at low field. Only
  that two Picard seeds both halt on the domain boundary.
- **Not established:** that `1 + Sigma0` crosses zero as a property of a converged solution.
- **Not established:** that the strict scheme *extends* the ordered domain relative to resummed at
  2–3 T. `resummed` was measured on the real multiset only at 1 T (`node_failed`) and at the two PM
  controls (non-convergent). Its status at 2, 2.5, 2.9 and 3.0 T was never measured head-to-head, so the
  comparison at those fields is inference, not data.
- **Not established:** that the masking symptom is fixed, or fixable, by any static-medium change.
- **Not established:** that the moment series is a *controlled* truncation anywhere near the critical
  path. §B was always explicit that it is a *formal order statement*, and clause (c) is now direct
  evidence that the first omitted term is not small on this density.
- **Not established:** the strict candidate's own PM phase boundary. §5.3's negative masses are
  measured against the *frozen controls*, not against a strict-scheme boundary, which was never
  determined.
- **Not established:** that the `ω_n = 0` scheme jump propagates into `invz_lambdas` or the spectra.
  §5.5 is a primitive-level measurement with `Sigma = 0`.
- **Not established** — and specifically **withdrawn** — that the path integrals fail to converge under
  grid refinement. They converge at second order (§5.4); they miss an absolute tolerance.
- **Not refuted:** the underlying `crit(h) = F'(h)` stability formulation of §A. It is independent of
  the medium truncation and remains available to any future candidate.

---

## 10. State of the branch

- **HEAD `6937e9d`**, 37 commits ahead of the Ewald closure `2ee310b`, 98 ahead of `main`. Nothing
  merged, nothing pushed.
- **Suite: `P=576 F=6 I=29 T=611`** (controller-verified independently of the implementer's report).
  The six failures are exactly the Task-17 gate tests described in §7, and nothing else:

  ```
  test_invz_strict_identities/test_G1a_dm_dh_equals_minus_G0bare
  test_invz_strict_identities/test_G1b_F_prime_equals_crit
  test_invz_strict_identities/test_G1c_dF_dm_equals_crit_over_chi_path
  test_invz_strict_identities/test_G1d_second_order_convergence
  test_invz_strict_identities/test_G2_onset_coincidence_at_m_zero
  test_invz_strict_contracts/test_G13_pm_slot_does_not_leak_into_ordered_lambdas
  ```

  That is **four G1 tests, G2, and G13** — an earlier revision miscounted them as "five G1 forms plus
  G13", dropping G2 from the inventory. `test_G1d_nested_trapezoid_is_second_order`, the synthetic
  quadrature control, **passes**. All six are red **by design**, pending the amendment.
- Gate-0 verdict document committed at `docs/invzp_strict_medium_gate0_report.md` (one non-gating
  erratum, §5.5 above).
- Untracked working records: `Task-17_failure.md`, `Task-17_prereg_amendment_DRAFT.md`.
- The frozen preregistration and the approved spec were **never edited**, verified by
  `git log a923691..HEAD` on those paths returning zero commits.
- ~30 unrelated dirty user paths and 11 stashes were preserved untouched throughout.

**Open decisions:**

1. Sign or reject the prereg amendment (§7). Until then the branch cannot be green.
2. Merge the branch as a diagnostic capability, or park it as a recorded negative result. The primitives
   in §3 are independently useful — the domain scanner, the coverage accounting, the two-tier verdict,
   the omitted-order diagnostics — and the default path is untouched.
3. A backlog of ~12 review findings the plan itself mandated (tests that assert nothing; a test whose
   name promises legacy-drift detection it cannot deliver, because it compares two runs of the new
   code; a third verbatim copy of the backaction block), plus one spec question: `phase_1z_reason` is
   declared a closed 11-value list and the implementation now has a twelfth, `response_failed`.
4. G5 tension: the Task-18 brief says Gate 0 cannot pass with a missing finest integral, but `rep.pass`
   is specified as exactly `~(a||b||c||d||e)` with no sixth term. G5 is measured, fails on every field,
   and is not a Boolean. Recorded, not silently reconciled.

---

## 11. What a next candidate would have to answer

Per prereg §3, any of the following is a **new theory candidate** requiring its own spec and a fresh
preregistration — never an in-run fallback from this one. Listing them is diagnosis, not a proposal:

1. **Where does the low-field breakdown actually live?** The failure is at large order parameter, not
   near the instability. The `1+Sigma0` reference was chosen for its `m → 0` cross-leg identity; that
   choice is untested at large `m`, and §5.1's window sits exactly there.
2. **Any candidate must survive clause (c) on the real density.** Since `omit_cubic ∝ mu2·Gref²`, a
   scheme is only controlled where `|Gref|` stays bounded. Retaining `mu3`, or `mu4`, does not obviously
   help — the ladder in §2 shows the domain gain is monotone in *removed* feedback, while the truncation
   error grows with retained order times `Gref` powers. Whether any truncation order has a window where
   both are acceptable on this multiset is an open, and answerable, question.
3. **The PM negative-mass result (§5.3) needs its own explanation** regardless of the ordered leg. It is
   the cleanest single anomaly in the run: converged, in-domain, negative, at both frozen PM controls,
   reproducible outside the driver. A necessary first step is to determine the strict candidate's *own*
   PM boundary, which was never established — the result implies it predicts instability through 3.5 T.
4. **The `ω_n = 0`-only restriction needs its propagation measured (§5.5).** The 17–59 % figure is a
   primitive-level, undressed comparison against Matsubara-frequency dispersion. Whether that jump
   survives into `invz_lambdas` and the spectra is **unmeasured** — and it is a cheap measurement.
   Either the truncation extends to the low-`ω_n` sectors consistently, or the jump has to be shown not
   to propagate.
5. **The G1 and G5 tolerance structure needs rebuilding, not just re-registering.** Both gates compare
   O(h²)-accurate discrete quantities against absolute thresholds chosen independently of the grid they
   are evaluated on, and both consequently fail while the underlying sequences converge cleanly at
   second order. Any successor preregistration should either target Richardson-extrapolated values or
   register convergence-order criteria — and must check each threshold's feasibility against the frozen
   grid *before* freezing both.

---

## 12. Review-verification log (rev. 2)

`invzp_fix_attmpt_review_Codex.md` (Codex, 2026-07-27) raised eight required corrections against
rev. 1. Each was checked against source, the recorded run, or arithmetic before being adopted. The
review's disposition — Gate-0 FAIL accepted, do not proceed to Stage 4 — is unchanged by any of them.

| # | claim | how verified | outcome |
|---|---|---|---|
| 1 | All three G1 forms are algebraically dependent, not two | re-derived `G1c = ρ(a+g)/(a·g)` independently; checked `g ≠ 0` from `chi_path > 0` | **ADOPTED** — §7 rewritten; amendment draft corrected |
| 2 | G5 misses its tolerance but the sequence *is* second order | computed all eight refinement ratios from the recorded table: `4.012–4.017` | **ADOPTED** — §5.4 rewritten, overclaim withdrawn |
| 3 | The "3500×" comparison is dimensionally invalid | units inspection | **ADOPTED** — removed, §5.4 |
| 4 | Suite inventory is four G1 + G2 + G13, not five G1 + G13 | the six `FAILEDNAME` lines from my own controller run | **ADOPTED** — §10 corrected |
| 5 | "Wrong sign above `Bc_1z`" overstates; `Bc_1z` is the legacy boundary | prereg §3(e) wording; strict boundary never determined | **ADOPTED** — §5.3 requalified |
| 6 | G7 is primitive-level (`Sigma = 0`), and `\|K(2)−K(1)\|` is frequency not q dispersion; propagation unmeasured | read `.superpowers/sdd/task18_g7.m` — line 20 `S = zeros(size(wn))`, line 24 `abs(a.K(2)-a.K(1))` on the resummed struct, no `invz_lambdas` call | **ADOPTED** — §5.5 requalified; committed-report erratum raised |
| 7 | `invz_gate0_aggregate` has no expected-PM-field list | read `invz_gate0_aggregate.m` — signature line 1, ordered coverage loop `:86-90`, PM loop `:98-103` iterates `numel(pm_rows)` only | **ADOPTED** — §6 qualified; latent defect recorded |
| 8 | `ratio >= 3` is post-hoc in timing even if structural in reasoning | timeline: the criterion was drafted after the failing run | **ADOPTED** — amendment draft relabelled |

All four wording changes were adopted: "five adverse findings" (G5/G7 are non-gating), the
boundary-adjacent phrasing in §5.1, and both retained distinctions (transient-vs-converged in §5.1,
removable-`Gstat`-vs-reference-denominator in §6).

**One item extends beyond the review.** Correction 6 also falsifies a sentence in the *committed*
verdict document (`docs/invzp_strict_medium_gate0_report.md` §4, "the coupling's own q-dependence").
That is flagged in §5.5 and awaits a decision — it is non-gating, but it is wrong.

**Net effect on the verdict: none.** Corrections 1, 3, 4, 6 and 8 concern gate design, arithmetic
presentation, or inventory. Corrections 2, 5 and 7 weaken specific *interpretations* — the G5
non-convergence claim is withdrawn outright, the PM sign claim is narrowed to a control violation, and
the coverage claim is narrowed to the ordered ledger. Clauses (a), (c) and (e) rest on none of them.

---

*Every measurement quoted here was produced by a preregistered run against a digest-verified coupling
fixture, or independently re-measured by the controller. No tolerance in this document was widened,
re-derived, or selected after seeing the value it was applied to. Where a claim in rev. 1 proved
wrong, it is withdrawn in place and named as withdrawn rather than quietly edited out.*
