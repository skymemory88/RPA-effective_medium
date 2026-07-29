# Blind second-opinion convergence plan and its adjudication

**This file merges two former documents** (both deleted 2026-07-29 to avoid
confusion): `blind_plan_adjudication.md` (now Part I) and
`blind_second_opinion_convg_plan.md` (now Part II). **Part I is the only operative
record. Part II is historical evidence, not an implementation specification, and
must not be executed as a plan**; its body is preserved verbatim, while its
title, status banner, and cross-references have been updated. The reviewed critique remains in
[`blind_plan_review_codex.md`](blind_plan_review_codex.md).

---

# Part I — Adjudication of the Codex review of the blind second-opinion plan

**Date:** 2026-07-29
**Inputs:** the blind plan (Part II below),
[`blind_plan_review_codex.md`](blind_plan_review_codex.md) (review),
[`invzp_convg_diagnosis.md`](invzp_convg_diagnosis.md),
[`invzp_convg_fix.md`](invzp_convg_fix.md),
[`jensen_1z_framework.html`](jensen_1z_framework.html) (targeted verification reads).
**Method:** each load-bearing review objection was evaluated; the two objections that
turn on what the framework document actually establishes were verified directly
against `jensen_1z_framework.html` rather than adopted on the reviewer's authority.

## 1. Verdict

The review is **substantially correct on every load-bearing objection**. Its bottom
line is adopted: the thermodynamic-integration (TI) selector must not become a
production phase-selection rule; the source-dependent stationary-functional program
in `invzp_convg_fix.md` remains the load-bearing root fix; the blind plan's
contributions are retained as bounded diagnostics, representation studies, and
solver infrastructure feeding that program.

**Re-evaluation result:** the revised program is accepted with the corrections
incorporated below. Three residual overclaims have been removed: route disagreement
is evidence of inconsistency but does not quantify non-integrability; a continuum
DOS result classifies the legacy equations rather than establishing physical
correctness; and a discovered-set report cannot use the word “minimum” until a
derived functional supplies the ordering. The historical Part II remains
non-operative.

**Independent re-verification (2026-07-29, second evaluator):** all corrections of
the re-evaluation were checked and accepted; burden-of-proof phrasing was
additionally tightened in §2.1/§3.1 (certification requires a positive derivation —
absence of counterevidence never suffices); the relation to the canonical program
is recorded in §6.

## 2. Point-by-point disposition

### 2.1 Φ = ∫F dm is not an exact selector (review §2.2) — **adopted, verified**

Verified directly: framework §9.4 (`jensen_1z_framework.html:516-521`) derives a
*difference* between the 1/z and pure-MF theories **at the same single-ion state**,
integrated from the saturation anchor (eq. 46), and explicitly designates agreement
of the (46) and (37) routes as "a stringent global check on a numerical
implementation" — i.e., not automatic. That check is measured to fail at ~13.7% in
the diagnosed ordered regime (diagnosis §7.4). This is evidence of route
inconsistency, but it does **not** by itself quantify non-integrability: the two
routes have different approximation scopes in the current hybrid, and their
disagreement is not a mixed-derivative or closed-loop measurement. The failed check
and the (coordinate-limited) non-gradient residual measurement corroborate the
refusal, but the operative bar is burden of proof: a line integral becomes a
certified potential only by derivation from a stationary functional, and no such
derivation exists — even perfect route agreement would not supply one (§2.3).
[verified fact + bounded inference]

### 2.2 The m = 0 anchor is not common (review §2.3) — **adopted**

Seven distinct h = 0 roots at 1.5 T (diagnosis §6.1) means per-sheet anchoring
Φ_sheet(0) = 0 silently equates unknown additive offsets — exactly the quantities a
cross-sheet ranking needs. The supplementary point is also correct: continuation in
T or Bx requires the −S dT and transverse-work terms of the full differential;
integrating H dm alone on such paths is wrong even for a consistent theory.
[verified fact (root count); elementary thermodynamics]

### 2.3 Route discrepancy is not an error bound (review §2.4) — **adopted**

|e₁−e₂| bounds neither |e₁| nor |e₂|; shared systematic content (omitted diagrams,
common normalization) defeats it. Corrected status: **disagreement falsifies the
claim that these two legacy routes evaluate the same hybrid potential; agreement is
necessary but never sufficient for that particular claim.** Because the second route
is not exact for the ordered 136-state hybrid, its agreement is neither necessary nor
sufficient for a future, separately derived 136-state functional unless that
functional reproduces both routes as controlled limits. ε_route is retained as a
*discrepancy map*, not an error bar. [logical fact]

### 2.4 The internal-energy route is not exact for the hybrid (review §2.5) — **adopted, verified**

Verified directly: eq. (35)/(36) are derived via the two-level equation of motion
(Ẋ governed by Δa₁₁ alone) with `G_X = G/M² (m=0)` stated in the derivation
(`jensen_1z_framework.html:423-429`); "exact given G" (line 429) is a statement
about the closed two-level paramagnet. Eq. (37) is tagged "population term is
Jensen's resummed form" (line 432). Use on ordered sheets of the 136-state
projected hybrid is an additional approximation, so the blind plan's "independent
exact route" claim is retracted. [verified fact]

### 2.5 DOS quadrature needs a limit-prescription contract (review §3) — **adopted**

The distinct-objects list (ordinary real integral / principal value / i0±
boundary values / balanced grid sequences) is correct; the framework's DOS remark
(§10 step 3) selects none of them. Adopted requirements for any DOS work:
(i) a declared prescription and order of limits; (ii) direct adaptive-q
cross-checks; (iii) separate convergence of band edge, B_c, accepted width, and
response; (iv) an explicit Γ/shape/uniform-term contract; (v) scope declaration for
field/T/frequency-dependent coupling rebuilds (ODD, retarded paths); and
(vi) **an in-band static operating point is a stability warning to be classified,
never a state to be admitted by choosing a principal value.** The DOS study keeps
its diagnostic value either way: it determines whether the 3.825 T-class masks
converge to finite-grid artifacts or to well-defined continuum-limit conditions of
the **legacy equations**. It does not establish that those equations or conditions
are the correct physical theory. [adopted; consistent with diagnosis §7.1]

### 2.6 Stable reassociation is incomplete (review §4) — **adopted**

`stable_form` defactors G̃₀ and r, but `Gstat` itself still diverges and
`invz_emt_static_ordered` evaluates `Gq = Gstat./(1+(Jf-K0).*Gstat)` — an Inf/Inf
form at the exact crossing. Any experiment must defactor the complete closure and
its residual (limiting G_q, G_bar, J_loc), starting from the defactored research
coordinate in `invz_ordered_residual.m`. [verified against the two independent
code reads; the specific line was cited identically by the blind agent and the
reviewer]

### 2.7 Budgeted enumeration ≠ completeness (review §5.2) — **adopted**

Certificate vocabulary is stage-dependent:

- before a stationary functional exists: `unique-accepted-within-budget`,
  `multiple-accepted-unranked`, `masked(cause)`, or `incomplete`;
- after a derived functional exists but before root completeness is certified:
  `lowest-Gamma-among-discovered`, never “equilibrium”; and
- `certified-equilibrium` is permitted only after validated root counting (interval
  Newton/Krawczyk or an equivalent certified reduction), functional-value error
  control, and relaxed-Hessian classification, i.e. fix-doc WP6–WP7/§7.2–7.4.

“Failure to find another branch is not proof none exists.” [logical fact]

### 2.8 Endpoint tier ≠ relaxed Hessian (review §5.3) — **adopted**

`crit_star`/`D_uni`/`Dq_min` do not test relaxed variations of Σ and K; the blind
plan's "read as Φ-curvature in m" is retracted. Thermodynamic stability
classification requires the Schur-complemented Hessian of the derived functional
(fix doc §7.3). [inference, uncontested]

### 2.9 Near-QCP language (review §5.4) — **adopted**

The 4.60–4.90 T fixture is a **"verified finite-grid state/response window,"** not a
"certified equilibrium solver" on any domain: no enumeration or functional ranking
was performed there either. All handoff language updated accordingly.

### 2.10 Evaluator-mode functional use (review §6) — **adopted with a recorded nuance**

Adopted: a nominal O(1/z²) equation-mismatch bound at selected points does not by
itself bound functional-*value* differences nor preserve the ordering of nearly
degenerate phases — and near-degeneracy is precisely the regime where selection
matters (first-order boundaries). Nuance recorded for the log: evaluating a true
functional Γ on an off-shell root incurs an error that is *second order* in the
state displacement only when a nearby stationary point of that same smooth Γ
exists, the displacement is bounded in the functional's natural norm, and the
relevant Hessian is bounded. This local Taylor fact does not bound truncation,
normalization, additive-constant, or missing-diagram errors. Evaluator mode is
therefore not categorically void, but it is admissible only with all of those
bounds compared against the actual phase gap—which in the near-degenerate regime
normally forces replacement mode. The exact-cluster coefficient gate (fix doc §6,
WP3) remains the safer first gate.

## 3. Precisions recorded so the log does not overcorrect

1. **The blind plan's premise was framework-textual, not invented.** The framework
   asserts "the thermodynamic relation dF/N = −gμ_B⟨J_x⟩ dH holds with the MF
   moment" as a consequence of tadpole cancellation
   (`jensen_1z_framework.html:517`). What fails is the *promotion*: a
   one-parameter, per-sheet antiderivative cannot be elevated to a globally
   consistent cross-sheet potential while the theory's own §9.4 two-route check
   fails at 13.7% — and, per §2.1/§2.3, could not be so elevated without a
   derivation even if that check passed. The identity may still hold to retained order on-shell for a
   properly derived functional; that is exactly what fix-doc WP4/WP5 must decide.
2. **WP1b is motivated by the framework's own test.** The two-route comparison is
   the §9.4 “stringent global check” generalized beyond its derived two-level,
   \(m=0\) scope. On the ordered 136-state hybrid it is an explicitly approximate
   discrepancy diagnostic, not a framework-mandated identity or an error estimate.

## 4. Resulting merged program

The canonical program of `invzp_convg_fix.md` (WP0–WP10) is load-bearing and its
sequencing (§12: begin with WP0–WP3, the exact-cluster coefficient harness) is
unchanged. The blind plan's items slot in as follows (superseding Part II §3/§6):

- **Immediate bounded diagnostics (days, existing machinery):**
  - Newton-polish the 4.05 T twins to tight backward error with conditioning
    estimates (blind WP1c; review §5.1 concurs) — determines whether the two
    observations converge to distinct high-accuracy numerical zeros. A local
    interval/Krawczyk enclosure is still required before calling either zero
    mathematically certified.
  - State-resolved two-route (46)-vs-(37) discrepancy map (blind WP1b,
    reinterpreted per §2.3 above) — approximation-consistency evidence, not an error
    bar or a quantitative integrability measure. Add a direct path-dependence
    measurement (two-path / closed-loop ∫m dH comparison at fixed T, Bx) as a
    necessary-condition test: a nonzero loop falsifies a common potential, while a
    numerically closed loop does not prove global integrability or fix cross-sheet
    constants. The compared paths must begin and end at the same complete
    \((m,\Sigma,K)\) state and must record every branch transition; matching only
    \(m\) or \(H\) is not a closed-loop test.
  - Φ-profile vs. 'last'-crossing comparison (blind WP1a) as a **diagnostic
    only**, tracking one verified sheet and carrying unknown integration
    constants explicitly.
- **Representation studies:** DOS/adaptive-q convergence study under the §2.5
  contract; results inform the lattice-representation decisions of fix-doc WP5
  and the convergence axes of WP8. Full-closure defactoring per §2.6. These studies
  classify the continuum limit of the legacy equations only; they cannot promote
  the legacy theory.
- **Solver infrastructure:** coupled-u Newton–Krylov, deflation, pseudo-arclength
  (blind WP3) is built as *discovery* infrastructure and reused directly inside
  fix-doc WP6; outputs are conditional discovered sets with declared budgets.
- **Reporting:** the certificate data model is adopted for interim reporting with
  the stage-qualified class names of §2.7/§2.9. Before Γ exists, the current map may
  use only `unique-accepted-within-budget`, `multiple-accepted-unranked`,
  `masked(cause)`, and `incomplete`; it may not report a minimum or equilibrium.
- **Latent defect, documented not patched:** the `find(..., 1, 'last')` crossing
  selection in `invz_hmf_ordered.m:362` is an undocumented convention. It is now
  recorded as such. It is **not** to be replaced by argmin Φ (also uncertified);
  it is replaced only by the derived functional's selection when WP7 lands.
- **Promotion boundary:** every legacy-equation diagnostic above remains
  default-off or outside the production path. No DOS prescription, defactored
  continuation, discovered branch, TI score, or polished root may alter production
  phase labels until the corresponding stationary-functional, completeness,
  stability, and regression gates pass.

## 5. Rejected branches added to the decision record

| Rejected branch | Rejecting evidence |
|---|---|
| TI/Φ = ∫F dm as production equilibrium selector | §9.4 scope (framework:516-521); 13.7% failure of the framework's own two-route check (diagnosis §7.4); non-common h = 0 anchor (diagnosis §6.1); review §2.2-2.4 |
| ε_route (route disagreement) as a certification error bar | measures \|e₁−e₂\| only; shared-systematic failure mode (review §2.4) |
| Eq. (36)/(37) as an exact independent route for the ordered 136-state hybrid | two-level, m = 0 derivation scope (framework:423-432); review §2.5 |
| PV-regularized acceptance of in-band static operating points | prescription is an analytic choice, not a quadrature refinement; in-band zero of inverse response is a stability warning (review §3.2) |
| `stable_form`-only enabling as the resummed-path fix | `Gq` Inf/Inf at crossing survives (review §4) |
| "certified-equilibrium" labels from budgeted enumeration | completeness unprovable without validated root counting (review §5.2) |
| Endpoint tier as thermodynamic Hessian | internal (Σ, K) directions untested (review §5.3) |
| Immediate argmin-Φ replacement of the 'last'-crossing rule | swaps one uncertified convention for another; selection waits for the derived functional (this adjudication §4) |

## 6. Relation to the canonical program (`invzp_convg_fix.md`)

**Independent convergence (mutual corroboration).** Produced blind, Part II
nevertheless converged on the canonical program's central claims: the root problem
is the missing equilibrium selector, not insufficient iteration (fix doc §1 ↔
Part II §0); the only sound selector is a single source-derived stationary
functional (fix doc §2.1 ↔ Part II §2.1 end-state layer); discovery machinery must
never decide phases (fix doc §7.1 ↔ Part II §2.2C); honest
`ambiguous`/`unranked`/`incomplete` outcomes are mandatory (fix doc §2.4 ↔ Part II
§2.2D); and a cheap falsification gate precedes expensive machinery (fix doc §12 ↔
Part II WP1). Because Part II was produced without access to the fix doc, this
agreement is independent evidence for the shared skeleton (review §9 concurs).

**Resolved divergence.** The one load-bearing divergence — Part II's attempt to
bridge to production with a TI/Φ selector before the functional exists — is
rejected (§2.1–2.3, §5). The canonical selection rule (stationary Γ value plus
relaxed Hessian, fix doc §7.3–7.4) stands unmodified.

**Where the canonical program is strictly stronger.** (i) Variables: fix doc §4.1
forbids the scalar/doublet projection until a controlled-limit identity or
measured error bound justifies it; everything in this document operates inside the
legacy scalar projection and therefore diagnoses the *legacy* theory only.
(ii) Certification: fix doc §7.2's interval/Krawczyk root counting is the only
route to completeness claims (§2.7 above). (iii) Functional origin: the fix doc
derives Γ from the exact source construction with finite-cluster coefficient
oracles (fix doc §4–§6, including the already-falsified varied-D candidate);
Part II's Ω sketch is a WP4/WP5 candidate input, not a derivation (review §6).

**What this document adds to the canonical program.**
1. The M1/M2/R decomposition — in particular the representation axis R with its
   concrete evidence (interval-rank saturation at the band top; ~N^(−1.1) ladder
   exponents) — as an input to fix-doc WP5 (lattice representation) and WP8
   (independent convergence axes).
2. Days-scale legacy diagnostics (§4 above) that can run alongside the weeks-scale
   fix-doc WP0–WP3 without touching production: twin polishing (also a validity
   check on the 4.05 T coexistence evidence cited in fix-doc §11 item 7), the
   state-resolved two-route discrepancy map, and the closed-loop integrability
   test with complete-state endpoints.
3. Concrete code-level findings: the `'last'`-crossing convention
   (`invz_hmf_ordered.m:362`), the incomplete `stable_form` defactoring, and the
   pseudo-root caution for near-tolerance accepted twins.
4. The staged certificate vocabulary and promotion boundary (§2.7, §4), which
   operationalize fix-doc §2.4's "completeness or explicit ambiguity" during the
   interim period that the fix doc's end-state gates do not themselves regulate.
5. The S3 triangular-Jacobian reading (Part II §1), which sharpens — and is
   consistent with — fix-doc §3's bounded interpretation of the asymmetric-Jacobian
   evidence.

**Priority rule.** Fix-doc §12 remains the critical path: WP0–WP3 (exact-cluster
coefficient harness) before any large lattice implementation. The §4 diagnostics
here are bounded side-packets; none may displace that priority or alter production
phase labels (§4 promotion boundary).

## 7. Document status

- **Part II below** — retained verbatim as a **historical, non-operative artifact
  that must not be implemented**. Its §2 selection
  principle, ε_route certification rule, "certified-equilibrium" class,
  Φ-curvature reading, and WP5 TI selector are **superseded** by this
  adjudication. Its M1/M2/R decomposition, WP1 experiments (reinterpreted), DOS
  study (contracted), solver infrastructure, honest-outcome classes, S3
  (triangular-Jacobian caution), and the `'last'`-crossing discovery stand.
- `invzp_convg_fix.md` — unchanged; remains canonical handoff document 2 of 2 and
  the load-bearing program.

---

# Part II — Historical blind artifact — DO NOT IMPLEMENT

> **Status (2026-07-29): superseded and non-operative.** This part is preserved
> verbatim only as the blind artifact and provenance record. Do not implement it
> or treat its work-package ordering as current. Its production selection rule
> (TI/Φ selector, ε_route
> certification, "certified-equilibrium" classes, Φ-curvature stability reading)
> is **superseded** by Part I above, following the review in
> [`blind_plan_review_codex.md`](blind_plan_review_codex.md). Its diagnostic and
> solver-infrastructure content stands as adjudicated there.

**Provenance:** produced 2026-07-29 by a subagent that was barred from reading
`invzp_convg_fix.md`, `backup_plan_biased_convergence_solution.md`,
`convergent_suggestions_claude.md`, and the deleted
`biased_convergence_solution.md` (including via git). Its inputs were
`invzp_convg_diagnosis.md`, `jensen_1z_framework.html`,
`invz_projected/README.html`, and direct code reading. The main-session
assistant had already had `invzp_convg_fix.md` auto-loaded into context via an
`@`-mention and therefore did not author this document; it is relayed verbatim
below. Blindness declaration is at the end of the plan.

---

# Independent fix plan — ordered-state non-convergence and equilibrium selection in the projected 1/z code

**Scope:** `invz_projected/` + `invz_common/` (Jensen 1/z, LiHoF4, 136-state electronuclear, scalar cc projection).
**Baseline:** `c015333`, as described in `invzp_convg_diagnosis.md` (read in full).
**Independence:** produced without reading `invzp_convg_fix.md`, `backup_plan_biased_convergence_solution.md`, `convergent_suggestions_claude.md`, or the deleted `biased_convergence_solution.md` (declaration at end).

---

## 0. Summary of the recommendation

The ordered-phase problem is **two distinct multiplicity problems stacked on one fragile numerical representation**:

- **M1 (node level):** at fixed molecular field h, the coupled state u = [Σ(iωₙ); K₀] can have several admissible solutions (the 4.05 T twins, diagnosis §6.2).
- **M2 (moment level):** at fixed (T, Bx), the moment equation F(h) = h₀(h) − J₀m(h) can have several roots/folds per solution sheet (1.5 T, diagnosis §6.1).
- **R (representation):** the EMT q-average is computed as a mean over a finite discrete coupling multiset (16³×4 = 16384 rational terms), which is meromorphic in the static medium with grid-specific poles exactly where the ordered static state operates (top of the cc band, Γ-adjacent; diagnosis §3.2, §7.1).

My recommendation, in one sentence: **fix R by evaluating the EMT average through a converged density-of-states (DOS) quadrature instead of a raw multiset sum; solve M2 immediately with the theory's own, already-latent Landau potential Φ(m) = ∫₀^m F dm′ (exact at retained order via the tadpole-cancellation identity dF = −m dH); solve M1 by sheet enumeration (deflated Newton + arclength) ranked by thermodynamic integration along connected sheets, cross-checked by a second route with an explicit error band; and build the order-consistent stationary functional as the end-state consistency layer — but only after a cheap two-route measurement at the known counterexample points has confirmed (or falsified) that ranking-by-integration is decisive.** The cheapest falsification experiment (WP1) costs days and uses only existing machinery.

---

## 1. Root-cause framing adopted, and where I sharpen the existing diagnosis

I accept the diagnosis's four-part mechanism (§7.1–7.4) as correct. Framing sharpened in four places:

**S1 — The four causes are not co-equal; separate representation from theory.** [inference] Causes 1–2 (grid meromorphy; resummed static denominator) are *representation/estimator* problems for the same N→∞ object the theory defines; causes 3–4 (non-contractive map; missing selector) are *solver/theory* problems. The evidence that a large share of deep-ordered masks are representation artifacts is concrete:
- The static-sector operating point y = K₀ − 1/G_stat sits in the top discrete interval of the coupling multiset (interval ranks 16384/16376 at the failed 4.3 T node, diagnosis §4; |Jf| = 16³·4 = 16384).
- The excluded Γ-adjacent gap shrinks like N^(−1.103) and the accepted ordered width like N^(−1.076) (diagnosis §3.2) — the signature of discrete gap structure at the band top, where the dipolar J_cc(q→0) is direction-nonanalytic and a finite grid samples it worst.
- The framework's own recipe prescribes a DOS representation for these averages (`jensen_1z_framework.html` §10, step 3: "the averages are efficiently done once per lattice via a density-of-states representation of J(q)").

[hypothesis, testable in WP2 gate] A substantial fraction of the 3.825 T-class masks and the interval-rank chatter will disappear under a converged DOS quadrature, without any change of theory and without any floor.

**S2 — The moment-level selector already exists inside the theory and is nearly free to compute.** [verified fact + inference] The code's F(h) is, unit-for-unit, the applied longitudinal field: h₀ = ∫₀^h r dh′ with r = G₀(0)/G̃₀(0) implements J 2.33 (framework eq. 45; `invz_projected/invz_hmf_ordered.m:443`), and F = h₀ − J₀eff·m (`invz_hmf_ordered.m:447`) is gμ_B·H_applied in meV. Because the tadpole-cancelling split makes ⟨J_z⟩ exact (J 2.24; framework §9.1, §9.4), dF_free/N = −m·d(gμ_B H) is an **exact identity of the truncated theory**, so along one connected solution sheet

  **Φ(m) ≡ ∫₀^m F dm′  (per site, meV), with dΦ/dm = F,**

is the free-energy difference from the m = 0 anchor, and the equilibrium among that sheet's stationary points is **argmin Φ** — a trapz over the already-computed profile. This also exposes a latent defect: production currently selects the **last** increasing crossing of F (`invz_hmf_ordered.m:362`, `find(..., 1, 'last')`) — an embedded, unjustified selection heuristic that should be replaced by argmin Φ. Note the crossing test is monotone-in-h; through folds (dm/dh < 0) "increasing in h" and "increasing in m" differ, so the stability read must be taken in m along the sheet.

**S3 — The Jacobian-asymmetry evidence is weaker than it looks, and the 4.05 T "true multiplicity" claim needs one cheap confirmation.** [inference] J₂₁ = 0 exactly with J₁₂ ≈ −1.66×10⁻³ (diagnosis §7.4) means the tested coordinate pair is structurally *triangular* — one tested equation simply does not depend on the other tested variable. That proves the residual is not a gradient *in those coordinates*; it is nearly no evidence about whether a functional exists in conjugate variables (e.g., (G, K) with the Matsubara measure), because a triangular Jacobian is exactly what a non-conjugate slice of a gradient system looks like. Do not fit metrics blindly; derive the functional and read off its natural coordinates (WP4). Separately: two "accepted" roots separated by 1.04×10⁻³ under an A-block tolerance of 1e-8 on a **non-contractive** map can in principle include a pseudo-root (residual below tol without a genuine zero nearby, because near-unit Jacobian eigenvalues widen the acceptance band in state space). Newton-polish both 4.05 T states to machine precision before treating M1 multiplicity as established (WP1c; the repair kernels under `docs/diagnostics/invzp_solver_stability_2026-07-27/` already exist for this).

**S4 — The theory, as assembled, is probably not exactly integrable — plan for a bounded inconsistency, not for perfection.** [inference from documented facts] Jensen's ordered elastic sector contains uncontrolled resummations: the tanh in ξ (J 2.29; `invz_common/invz_gstat_ordered.m:44-45`), the static replacement g(ω_m ± ω_n) → g(ω_m) (framework §9.2), and the resummed static closure with its own local denominator 1 + Σ₀ + K₀G₀inel₀ (`invz_gstat_ordered.m:46`). These are same-retained-order choices, so the state equations need not be the exact Euler–Lagrange system of any single functional; the 13.7% two-route ΔF disagreement (diagnosis §7.4; `invz_projected/invz_deltaF_ordered.m:2-10`) is the measured size of that inconsistency in one regime. Consequence: **any selection rule must carry an explicit error band**, and "degenerate within theory error" must be a legal, reportable outcome. A selection rule that always returns a unique winner would be claiming a precision the truncation does not have.

---

## 2. The fix

### 2.1 Theory level — the selection principle and its mathematical objects

**Variables (per (T, Bx) column):** the single-ion state at molecular field h (fixed-field, no ordering update — P0-1); the node state u = [Σ(iωₙ); K₀] (the acceptance-contract state, diagnosis §2.4); the moment m(h) = ⟨J_z⟩; the auxiliary applied longitudinal field H (continuation axis); sheets 𝒮 = maximal connected solution manifolds of the coupled equations over (h; T, Bx, H).

**Equations (unchanged production theory):** dynamic EMT (framework eqs. 16–17, 31); ordered Σ (eq. 39–40; `invz_common/invz_sigma_ordered.m`); elastic static closure (eqs. 41–42; `invz_common/invz_gstat_ordered.m` + `invz_projected/invz_emt_static_ordered.m:58-71`); H_MF↔H relation (eq. 45; `invz_hmf_ordered.m`).

**Selection principle (production rule):**

1. **Admissibility** = the existing A–D contract + finiteness + domain validity (unchanged; diagnosis §2.4). Admissibility never selects.
2. **Within one sheet (M2):** equilibrium = argmin over the sheet's stationary points (m = 0 and F = 0 crossings) of the sheet potential Φ(m) = ∫₀^m F dm′, evaluated along the sheet with quadrature error tracked. This is exact at retained order because dF_free = −m dH is an identity of the tadpole-cancelling split (framework §9.4). Stability of the winner = the existing endpoint tier (crit_star > 0, D_uni > 0, Dq_min > 0), now read as Φ-curvature in m.
3. **Across sheets (M1):** rank sheets by thermodynamic integration (TI) of the same exact differential along connected paths to a **common reference where the solution is unique** — the m = 0 PM state at h → 0 where sheets merge with the PM predictor, or, for sheets not reaching h = 0, continuation in the auxiliary applied field H (or in T) to a merge/fold connection; Φ is continuous through folds because it is a line integral along the connected curve. Two sheets anchored to the same reference are ranked by their Φ values.
4. **Error band:** every ranking is accompanied by a second, independent route — the internal-energy route: U from framework eq. 36 (exact given G), F by temperature integration of eq. 37 along the sheet from a PM-connected point. The certification requirement is |ΔΦ(sheet A, sheet B)| > ε_route, where ε_route is the measured two-route disagreement for those states (the honest generalization of the 13.7% number). Otherwise the point is classified **degenerate-within-theory-error**, reported as such.
5. **Sheets that cannot be connected to any common reference within the declared continuation budget** are reported **unranked** (with the budget stated). "Not found/not connected within budget" is never converted into "absent/metastable."

**End-state consistency layer (the stationary functional).** The right long-term object is the order-truncated linked-cluster/Baym–Kadanoff-type functional for a single site in the retarded cavity field K with the EMT embedding term (the scalar-EMT analog of the DMFT/Soven functional):

  βΩ[m or H_MF; Σ(iωₙ); K(iωₙ), K₀] = βF_ion(H_MF) + Legendre terms (h−H)m − J₀m²/2 + medium-embedding trace terms Σₙ wₙ ⟨ln[1 + (J(q)−Kₙ)Gₙ]⟩_q + the first-order single-site term built from the one-K-line cumulant including the S₄ semi-invariant (framework eqs. 19–29),

with stationarity required to reproduce: δΩ/δK = 0 ⇒ eq. 16; δΩ/δΣ = 0 ⇒ eqs. 24/39; ∂Ω/∂H_MF = 0 ⇒ eq. 45; and Hessian positivity ⇒ the stability tier. Substantial prototype machinery already exists and should be audited and reused, not rebuilt: `invz_functional/invzf_ring_scalar.m` (ring free energy (1/2β)Σₙwₙ⟨ln(1−J_qCₙ)+J_qCₙ⟩_q with analytic derivatives and a fail-closed pole status), `invzf_scalar_functional.m` (f = f₀(h) + (h−H)m − J₀m²/2 + f_ring with gradient/Hessian), and `invzf_stationary_scalar.m` (root enumeration + Schur-curvature stability + lowest-functional selection) — currently two-level research oracles, not production (diagnosis §2.1). The known obstruction (S4 above) must be confronted head-on: the derivation must state exactly where Jensen's ξ-tanh/static-g elastic sector deviates from the functional's Euler–Lagrange equations and bound the deviation as O(1/z²); the outcome decides whether the functional is used (i) as an evaluator on Jensen solutions (if the deviation is within retained-order error) or (ii) as a replacement scheme promoted through a frozen Gate-0-style predicate.

### 2.2 Algorithm level

**A. DOS-quadrature EMT (fixes R; not a regularizer).** Replace the discrete-multiset averages ⟨·⟩_q in both the dynamic slot (`invz_projected/invz_emt_scalar.m` call sites) and the static ordered closure (`invz_emt_static_ordered.m:62-65`, `Gq = Gstat./(1+(Jf-K0).*Gstat)`) with integrals ∫dJ ρ(J)(...) against a once-per-lattice converged density of states ρ(J) of the cc coupling branches, evaluated by piecewise-analytic panel integration (log/rational antiderivatives), with the Γ-point directional nonanalyticity handled by explicit angular averaging of the q→0 limit. Inside-band static operating points then meet the correct branch-cut/PV structure — a *reportable physical* condition (overlap with the coupling continuum), not a spurious inter-grid-point pole. This is the framework's own §10-step-3 prescription and is a refinement of the estimator of the same N→∞ object — it is not a pole floor, not clipping, and changes no equation. Convergence must be demonstrated by a DOS-resolution ladder replacing the N-ladder of diagnosis §3.2.

**B. Exact stable reassociation on the resummed path.** `invz_gstat_ordered.m:49-66` already contains the exact (≤ ~1 ulp) reassociation that makes the G_stat local-denominator crossing removable in G̃₀ and r (limits G̃₀ → −1/K₀, r → −G₀bare·K₀), currently wired only for strict mode; the historical arithmetic manufactures Inf/(−Inf) = NaN at the crossing (spec G17 note in the file). Enable it for the resummed production path behind its own mini-gate (bit-comparison away from the crossing per the G9 contract). This is exact algebra, not a floor.

**C. Simultaneous coupled solve + enumeration.** Reformulate the node problem as one rootfind in u = [Σ; K₀] (the residual checker's 'coupled' formulation already defines the residual; `invz_projected/invz_ordered_residual.m:506-531`), solved by Newton–Krylov with trust region, seeded by short Picard; **deflated Newton** for bounded enumeration of coexisting roots at fixed h; **pseudo-arclength continuation** in h (and in the auxiliary H and T axes) to build sheets through folds. The nested static Picard (`invz_emt_static_ordered.m:60-71`) survives only inside the Picard initializer. The diagnostic Newton kernels (`docs/diagnostics/invzp_solver_stability_2026-07-27/`) and branch-topology scripts (`docs/diagnostics/biased_smooth_r_2026-07-28/`) are the seed code. Enumeration output is a declared-budget solution *set*, never a single seed-dependent state.

**D. Certificate data model.** Each (T, Bx) column exports one certificate: {enumerated sheet set + budget; per-sheet Φ profile + quadrature error; TI anchors used; two-route ε_route; winner or class ∈ {certified-equilibrium, degenerate-within-error, unranked-disconnected, masked(cause)}; stability tier of the winner}. Certificates are functions of the enumerated set only — no seed-dependent labels by construction.

---

## 3. Work packages

**WP0 — Freeze and record (0.5–1 day).** Decision record for this plan; freeze coupling digest provenance before any new gate (the §9 caveat: current digest `499922e6…` vs the stale Gate-0 constant `ddb9532d…` in `invz_projected/invz_gate0_report.m` — never silently update; freeze constructor+options+ordering+digest together). Gate: record merged; digests reconciled or explicitly quarantined.

**WP1 — The cheap falsification measurements (2–4 days; existing machinery only). This runs first.**
- **1a (M2 selector oracle):** at 1.5 T (folds) and one mid-window field, compute Φ(m) = ∫F dm′ by trapz on existing accepted profiles; rank all F = 0 crossings and m = 0; compare against the current 'last'-crossing choice; quantify quadrature error by nH-refinement. Extend/compare with the retained `docs/diagnostics/invzp_qcp_grid_2026-07-28/invzp_area_rule_oracle.m`.
- **1b (two-route error model):** for one state pair with both routes computable, evaluate the Φ/TI route and the U(36)/(37)-temperature route between the same two states; measure ε_route in the ordered window (turn "13.7%" into a state-resolved error model).
- **1c (4.05 T twins):** Newton-polish both accepted states to machine precision (repair kernels); confirm genuine M1 multiplicity vs pseudo-root; if genuine, attempt to connect each to h → 0 or to a fold via arclength continuation and rank by Φ with the 1b error band.
- **Gate (decision point):** if the two routes agree on the winner with |ΔΦ| > ε_route → TI ranking is viable as the production rule and WP4 becomes a consistency upgrade; if routes disagree or the margin is smaller than ε_route at representative points → WP4 (functional) is promoted to load-bearing and blocks WP5. Either outcome is recorded progress.

**WP2 — DOS-quadrature EMT (1–2 weeks).** Build converged ρ(J) (with Γ-limit angular treatment); implement panel-analytic averages for dynamic and static sectors; enable the exact reassociation (2.2B) behind its mini-gate. Gates: (i) PM certified window (diagnosis §3.1: 19 ordered + 42 PM columns, B_c bracket 4.690–4.695 T) reproduced within stated tolerance; (ii) DOS-resolution ladder shows a converging accepted-width and B_c (replacing the non-converged N-ladder of §3.2); (iii) re-run the 3.825 T sliver and the 4.3 T failed node: every previously masked node either solves or is reclassified with a *physical* cause (in-band continuum overlap / genuine fold / no solution) — count the artifact fraction (tests hypothesis S1).

**WP3 — Coupled Newton + enumeration + sheets (1–2 weeks, parallel with WP2 after its interface freezes).** Productionize 2.2C. Gates: recover both 4.05 T roots from neutral deflation seeds; recover the 3.6 T Newton-repairable nodes; zero unexplained `max_iter` terminations on the test fields — every failure classified.

**WP4 — Stationary functional derivation (2–4 weeks, parallel track; blocking only if WP1 gate fails).** Derive Ω (2.1); audit/reuse `invz_functional/` prototypes; hard cheap anchor first: verify on a small Matsubara grid that δΩ reproduces the PM-sector equations to machine precision (m = 0, where no elastic-resummation ambiguity exists), and that ∂(βΩ)/∂β reproduces eq. 37 and the two framework §9.4 consistency tests. Then the ordered sector with the explicit remainder terms quantifying the ξ-tanh deviation. Gate: a written derivation with the deviation bounded and a numeric identity check; a scheme-promotion decision (evaluator vs replacement) recorded with evidence.

**WP5 — Production selection rule + certificates (1–2 weeks).** Implement 2.2D on top of WP2+WP3, ranking by Φ/TI with ε_route bands (or by Ω if WP4 promoted a replacement). Replace the 'last'-crossing heuristic with argmin Φ. Gates: certificates for the full 3.5–5 T scan at T = 0.10 K; the 4.05 T column resolves to certified or degenerate-within-error (never silent); no seed-dependence under seed scrambling of the enumeration.

**WP6 — Regression, certified-domain statement, docs (≈1 week).** Publish the certified domain map; freeze new fixtures and gates (with WP0 digest discipline); update `invz_projected/README.html` operational status.

---

## 4. Disposition of existing machinery

| Machinery | Disposition | Reason |
|---|---|---|
| A–D acceptance contract (`invz_ordered_residual.m`) | **Keep** — as the admissibility layer only | Necessary state-validity check; never a selector (diagnosis §2.4) |
| Safeguarded Aitken (`invz_ordered_node_solve.m:346-455`) | **Keep, default-off**, inside the Picard initializer only; retire from any "convergence fix" narrative | Correctly safeguarded but deliberately narrow (diagnosis §4); superseded by Newton (WP3) |
| Warm continuation (`qcp_down`) | **Repurpose** as enumeration/continuation seeding inside WP3 | Useful discovery tool; not a selector (diagnosis §2.3, §8.1) |
| Newton repair kernels (`docs/diagnostics/invzp_solver_stability_2026-07-27/`) | **Promote/repurpose** into WP3's production Newton–Krylov | Proven to recover 3.6 T nodes (diagnosis §5) |
| Branch-topology tools (`docs/diagnostics/biased_smooth_r_2026-07-28/`) | **Keep** as WP3 enumeration infrastructure; smooth-r objective stays **retired** | Topology probes valid; smoothness-as-thermodynamics banned (diagnosis §8.3) |
| ΔF routes (`invz_deltaF_ordered.m` + U/T route) | **Repurpose** as WP1's measurement instruments; upgrade to common-reference TI with error bands; never an absolute F | Two-route disagreement is the error model, not a defect to hide (diagnosis §7.4) |
| Strict-medium candidate (`invz_medium_moment_closure.m`, `invz_static_medium_reference.m`) | **Keep** as diagnostic comparator; do not promote | Failed Gate-0 clauses (a),(c),(e) (diagnosis §8.2); DOS work is orthogonal to truncation order |
| Gate-0 harness (`invz_gate0_aggregate.m`/`invz_gate0_report.m`) | **Keep** the pattern for all new promotions; fix digest provenance first (WP0) | §9 provenance caveat is binding |

---

## 5. Rejected alternatives, with rejecting evidence

1. **More Picard iterations / mixing / damping tuning** — rejected: map is non-contractive and basin-sensitive (diagnosis §7.3); 3.825 T fails on both cold and warm seeding (§5); iteration dynamics do not alter the fixed-point set (§8.1, binding).
2. **General Aitken/extrapolation acceleration** — rejected as a fix: the safeguard itself correctly refuses outside its alternating-tail signature (4.3 T node, §4); §8.1 binding.
3. **Warm continuation as equilibrium selector** — rejected: continuation-connected component at 4.05 T has no in-domain closing endpoint while a second accepted root exists (§6.2); §8.1 binding.
4. **Smoothest-r / nearest-previous-root selection** — rejected: §6.2 is the decisive counterexample; retired in §8.3 (binding).
5. **Pole floors, clipped denominators, bare-state substitution, masking-and-interpolating** — rejected: repo standards + §8.1 (binding). The DOS quadrature (WP2) and exact reassociation (2.2B) are categorically different: converged estimators/exact algebra of the same defined object, gated by reproduction of certified results.
6. **Promoting the strict one-shot static medium** — rejected: failed Gate-0 clauses (a),(c),(e); ordered state persisted to ~4.6 T with no state 4.7–4.85 T (§8.2). Absence of the legacy pole ≠ complete thermodynamics.
7. **Brute-force coupling-grid refinement (bigger N) instead of DOS** — rejected: measured ladder N = 12→24 shows B_c still moving (+0.021 T) and accepted width shrinking ~N^(−1.08) (§3.2); refinement densifies poles rather than removing the class (§7.1).
8. **Ranking by the existing partial ΔF as-is** — rejected: routes disagree by ~13.7% and the function is explicitly outside its validated domain as an absolute selector (§7.4; `invz_deltaF_ordered.m:2-10`). Usable only after WP1's common-reference + error-band upgrade.
9. **Forcing gradient structure by metric-fitting the current residual coordinates** — rejected: J₂₁ = 0 exact (structural triangularity) cannot be repaired by any diagonal metric (S3); conjugate variables must come from the derived functional, not from a fit.
10. **Real-axis broadening** — rejected: failure occurs during imaginary-axis state construction (§8.1, binding).

---

## 6. Cost/risk and the minimum path to an honest production answer

| WP | Relative cost | Main risk |
|---|---|---|
| WP0 | trivial | none |
| WP1 | small (days) | routes may disagree → decision, not loss |
| WP2 | moderate | in-band PV/edge handling must be honest; mitigated by artifact-count gate |
| WP3 | moderate | enumeration completeness unprovable → mitigated by declared budgets + honest classes |
| WP4 | large | elastic-sector deviation may be O(1/z) rather than O(1/z²) → forces scheme decision; the PM-sector anchor makes failure visible early |
| WP5 | moderate | ε_route may be large in deep-ordered regime → "degenerate" class absorbs it honestly |
| WP6 | small | none |

**Minimum honest path:** WP0 → WP1. If WP1's gate passes: WP2 + WP3 → WP5-lite (TI ranking with error bands) → WP6, with WP4 as a parallel consistency upgrade. If WP1's gate fails: WP0 → WP1 → WP4 (blocking) → WP2/WP3 → WP5.

**Interim deliverable (available immediately after WP0+WP1, before any large implementation):** re-issue the current production map with per-column certificate classes — *certified* (unique accepted state, stable endpoint, within the verified near-QCP window of diagnosis §3.1), *multiple-accepted-unranked*, *masked (cause)* — plus the explicit statement: "the current theory is certified as an equilibrium solver only on the near-QCP window (4.60–4.90 T fixture at T = 0.10 K, finite-16³ regression sense, B_c^{1/z} = 4.6925 T bracketed; mode shape robust to 0.38–0.53% under grid alignment, absolute field coordinate not grid-converged); everywhere else states are locally admissible research output." That statement is already fully supported by the diagnosis (§3, §11) and requires no new computation.

---

## 7. Repo-standards compliance

- No pole floors / clipping: WP2 is a converged estimator of the defined N→∞ average; 2.2B is exact algebra with a removable-singularity limit; both gated on reproducing certified results.
- No masked-to-bare substitution: certificate classes preserve masks with causes.
- No seed-dependent phase labels: certificates are functions of the enumerated set with declared budgets; seed-scrambling gate in WP5.
- No smoothness-as-thermodynamics: smooth-r stays retired; ranking is by Φ/TI/Ω only.
- No unperformed-verification claims: every gate above is a measurement; "degenerate-within-error" and "unranked" are legal outcomes; ε_route is measured, not assumed.

---

## Evidence actually consulted

- `invzp_convg_diagnosis.md` — in full (428 lines).
- `jensen_1z_framework.html` — headings; §3 (bare propagator), §4 (power counting/RPA), §5 (effective medium), §6 (1/z series, S₄, Σ, sum-rule α), §7 (sum rule), §8 (energy/heat capacity), §9.1–9.4 (ordered phase, elastic Σ, H_MF relation eq. 45, free-energy consistency), §10 (computational recipe, DOS remark, boundary search), §11 intro.
- `invz_projected/README.html` — §1.5–1.9.1 (ordered phase, shape/ODD context), heading map.
- Code read in full: `invz_projected/invz_ordered_residual.m`, `invz_ordered_node_solve.m`, `invz_solve_point_ordered.m`, `invz_hmf_ordered.m`, `invz_deltaF_ordered.m`, `invz_emt_static_ordered.m`; `invz_common/invz_medium_moment_closure.m`, `invz_static_medium_reference.m`, `invz_sigma_ordered.m`, `invz_gstat_ordered.m`.
- Headers only (as evidence of existing prototypes): `invz_functional/invzf_scalar_functional.m`, `invzf_stationary_scalar.m`, `invzf_two_site_exact.m`, `invzf_ring_scalar.m`.
- Directory listings: repo root, `invz_projected/`, `invz_common/`, `invz_functional/`, `docs/diagnostics/`, `References/`.
- Reported-not-rederived (flagged as such where used): the 13.7% two-route number, J₁₂/J₂₁ values, 4.05 T twin-root measurements, N-ladder numbers, 3.825 T node census — all from the diagnosis; no MATLAB was executed (read-only constraint).

## Blindness declaration

The blindness constraint was maintained. I did not open, grep, or retrieve (via git or otherwise) any content of `invzp_convg_fix.md`, `backup_plan_biased_convergence_solution.md`, `convergent_suggestions_claude.md`, or the deleted `biased_convergence_solution.md`. My only exposure to prior proposal content is the four-line summary of the fix-document's program in `invzp_convg_diagnosis.md` §11 (unavoidable per the task brief) and incidental one-line pointers to that document inside code comments/docstrings (e.g., tolerance-constant summaries citing "invzp_convg_diagnosis.md §8.2"). The plan above — in particular the Landau-potential/TI selection rule with a measured two-route error band, the DOS-quadrature EMT, the deflation/arclength sheet enumeration, and the ordering of the falsification experiment before machinery — was developed independently from the theory framework, the diagnosis's evidence sections, and the code itself.
