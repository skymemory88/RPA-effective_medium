# invzp ordered-leg non-convergence — implementation plan

**Written 2026-07-27.** Derived from the consolidated diagnosis in
`invzp_convg_diagnosis_Claude.md` (source-first, then reconciled against `docs/`, an independent QCP
analysis now preserved in its §9.5, and external review).

**Independently audited 2026-07-27 against the current MATLAB source.** The upstream diagnosis is
substantially supported, but the original Phase 1 was not implementable as written:
`invz_ordered_residual` is an acceptance checker that returns scalar norms and pass flags, not a vector
residual for Newton, and its Block A calls the nested fixed-point map that the experiment is meant to
replace. The continuation experiment and its decision fork have therefore been rewritten below. The
audit also corrects the proposed node census, trace placement, h-grid work, and the scope of `omit_max`.
A fresh MATLAB R2025a production-path run gave
`hstar=NaN`, `status=node_failed`, 33 profile nodes/6 accepted, and 34 ledger nodes/6 accepted with the
predictor failed. A separate R2025a semantics check confirmed that bare `max` omits NaN by default.
Subsequent review feedback was folded in only where it survived an independent source/algebra audit:
the Phase-1 Jacobian structure is explicit, the static residual is defactored to remove its
`Gstat = 0` normalization degeneracy, and the four-variable reduction remains only a
branch-conditioned analysis tool. The final pass replaces the coordinate-singular `Q` event by the
physical `mean(Gq)` event, evaluates the diagnostic in the exact reciprocal chart `z = 1/Gstat` with
`stable_form = true`, and defers structured-linear-algebra optimization until a profile or a
low-temperature run justifies it.

**Execution update, 2026-07-27.** First-hand inspection of the current `invz_run_spectra` output added
an important constraint: convergence is non-uniform across the ordered field range; usable columns
can occur near the apparent QCP while moderate-field columns (for example 1.5 T) are masked. A direct
solver-only pilot on the exact legacy 16³/brute-force production coupling confirmed that pattern.
At `T = 0.10 K`, Picard accepted 30/33 profile nodes at 3.6 T but only 11/33 at 1.5 T. Reversing the
3.6 T node traversal left exactly the same three failed nodes, whereas a defactored residual Newton
corrector repaired only those nodes in 4/4/8 iterations; the unchanged HMF equations then gave
33/33 nodes, `hmf = 0.01949623263696515 meV`, and a stable endpoint. At 1.5 T the same method repaired
the predictor and two nodes but left the path incomplete (13/33). Thus a local numerical-stability
defect is established at 3.6 T, while the moderate-field failure needs continuation rather than a
one-node solver swap. A subsequent 300-node natural-`h` trace at 1.5 T remained A–D accepted down to
`h = 0.003414300661 meV` without a pole or mean-cancellation event, but it did not have a sufficient
branch-identity gate. Repeating the trace with scaled variables and an analytic tangent predictor
instead stopped at `h = 0.005243548157786061 meV`, where `rcond(J) = 1.027e-11`,
`|deta/ds| = 1.771e-6`, and `sigma_min(J) = 9.672e-6`, while the augmented tangent block remained
regular with `sigma_min([J R_eta]) = 0.4581`. Scaled pseudo-arclength crossed this fold with all A–D
audits accepted and changed `deta/ds` from `-1.771e-6` to `+8.358e-4`. Thus the clean high-`h`
branch turns back toward increasing `h`; the earlier lower-`h` secant roots were reached by a branch
jump and are not evidence of a continuous Jensen path. Repeating the local segment at one-quarter
the initial arc step resolved `min(h) = 0.005243548147122800 meV`, bracketed the zero of `deta/ds`,
and kept successive oriented-tangent overlap above `0.9999999999`.

An independent low-`h` test then used three declared deterministic cold seeds, none derived from the
high-`h` trace. All converged to the same `h = 0` root within a scaled-state diameter of `4.278e-12`.
Upward tangent continuation reached a second fixed-`h` rank loss at only
`h = 1.130400753628218e-5 meV`: `sigma_min(J) = 1.558e-5` while
`sigma_min([J R_eta]) = 0.5629`. Pseudo-arclength crossed it with all A–D audits accepted and changed
`deta/ds` from `+2.677e-7` to `-5.243e-7`, resolving a local maximum
`h = 1.130402502867847e-5 meV`. The directly traced low- and high-`h` neighborhoods are disjoint and
turn away from the intervening range. Roots returned by the earlier secant experiment between the two
folds must lie on at least one additional branch segment. More remote folds could in principle join
segments globally, but the two local folds already establish multiple real fixed-`h` roots and hence
the branch-selection obstruction.

That result changes the near-term priority, not the rigor threshold. The retained Newton kernel lives
under `docs/diagnostics/invzp_solver_stability_2026-07-27/` and is not wired into production HMF:
an A–D-accepted fixed-`h` root is not automatically the continuous branch required by the Jensen
integral. The low-`h` overlap test is now complete and has reached the second row of the Phase-1 fork:
multiple reproducible real branches/components. No fixed-node Newton fallback is eligible for
production, and a wider field map or expensive grid/cutoff funnel would not resolve which component
belongs in the Jensen integral. The next implementation requires a separately justified and
preregistered branch-selection prescription (or the common-functional formulation). The intermediate
secant-root component may be traced only if its topology is needed to choose or falsify such a
prescription.

**Retained continuation oracle, 2026-07-28.** The defactored residual/Jacobian and exact HMF
fixed-`h` node construction are now exposed only under diagnostics, together with a scaled bordered
pseudo-arclength tracer. On the frozen 0.10 K, 1.5 T legacy 16³ fixture, 149 full-state roots traced
the clean high-`h` component through its regular fold; every independent A--D audit passed. A local
Hermite tangent reconstruction gives
`h_fold=0.0052435482911986821 meV`, only `1.44e-10 meV` from the earlier one-off quarter-step result.
This validates the retained continuation primitive, not a Jensen path: disconnected root discovery,
single-valued section construction, endpoints, reverse/cold-seed agreement, and refinement gates
remain outstanding. The dense oracle cost 1292.9 s, making the already derived
diagonal-plus-low-rank bordered solve the next measured implementation requirement.

**Theory-route update, 2026-07-28.** The explicit-prescription backup is now specified in
`biased_convergence_solution.md`, with global smoothness of `r(h)` as its primary declared bias and
state/QCP continuity as binding gates. The preferred-route audit is recorded in
`invzp_functional_derivation_audit.md`: Jensen's field equation supplies a one-dimensional potential
only after an internal branch has been selected, and the scalar EMT closure has an isolated log
potential, but neither is a common functional for the coupled resummed equations. The current Newton
residual is not a direct or diagonally rescaled gradient. The strict scalar two-level WP0 derivation
is now closed at the unicyclic `O(1/z)` scope in `invzp_functional_wp0_spec.md` and
`invzp_functional_wp0_ring_derivation.md`. The isolated `invz_functional/` prototype implements the
ring common functional and exact two-site oracle without modifying or dispatching from the production
solver. This remains the preferred staged route.

**First production-input functional pilot, 2026-07-28.** The isolated strict functional now has a
read-only adapter to the actual scalar BZ coupling multiset and transverse electronic doublet. At
`T=0.1 K, Bx=1.5 T` it finds two stable, degenerate ordered minima on the legacy `16^3` couplings;
the positive state is `h=0.0329775225244 meV`, `m=5.13313920359`, with minimum ring denominator
`0.9322`. Complete root re-enumeration is stable over Matsubara cutoffs `32,64,128` and BZ grids
`8,12,16`. This is independent evidence that the moderate-field Jensen masking is not simply the
absence of an ordered state. It is not yet a spectral fix: the fixed electronic-doublet mode is
about `85 GHz` even after the uniform RPA shift, whereas the reported window is `0--6 GHz`, so the
electronuclear manifold is essential. The strict electronic pilot also develops a Gaussian-pole
no-state interval around `3.2--3.5 T`, between its ordered and paramagnetic solutions. Production
integration is therefore still prohibited; the next isolated gate is a source-biased
electronuclear local reference, followed—only if the no-state interval remains—by a newly specified
stationary skeleton/2PI extension.

**Electronuclear gate result, 2026-07-28.** A full 136-state source-biased local oracle now supplies
`f0`, `m`, connected `C2(iwn)`, and nested-stencil source/beta derivatives from one Hamiltonian.
It agrees with the independent `invz_chi0z` path to `1.9e-13` and closes the local Maxwell identities
to roughly `1e-6` relative or better. To preserve a single generator, this first oracle uses
`transverse_mf='none'`; the production transverse molecular field would require its own conjugate
functional variables. On the actual `16^3` couplings it finds a cutoff-stable `1.5 T` minimum at
`h=0.0340929334913 meV`, `m=5.30675931175`. Hyperfine structure restores ordered minima at
`3.6,4.0,4.5,4.6 T` and a stable symmetric state at `4.9 T`, but no strict-ring stationary root at
`4.7,4.8,4.85 T`. The pole/no-state interval is narrower and on the experimental field scale, but
it survives. The next candidate must therefore be a stationary nonlocal-return skeleton containing
the complete one-`C4` and `C3-C3` 2PI cores; a zero-frequency-only mass patch is disallowed.
That candidate and its stop gates are now pre-registered in
`invzp_functional_wp4_skeleton_spec.md`.

**WP4 entry-gate result, 2026-07-28.** The smallest skeleton has not yet been run, but its immutable
two-level inputs are now complete.  Exact frequency-labelled `C2/C3/C4` and 1PI
`gamma2/gamma3/gamma4` use ordered-simplex matrix exponentials and pass the source-derivative,
permutation, zero-source beta-delta, and near-degenerate gates.  A centred pair fixes the leading
`C3-C3` coefficient; a centred three-site chain fixes the nonzero-source mixed-`C4` coefficient that
must cancel reducible `C3*C2*C3` content after 1PI conversion.  Optional q-resolved Hermitian
coupling pages are exposed for both brute-force and Ewald backends without changing ordinary
production outputs.

**WP4 gate verdict, 2026-07-28.** The preregistered varied-covariance ansatz was then implemented
temporarily and failed that mixed-chain gate already at zero source, where `C3=Phi33=0`.  The failure
is analytic: varying a Gaussian local covariance trace plus the exact local 1PI `gamma4` generates a
spurious local-return denominator and changes the retained `a^2b^2` coefficient.  The temporary
solver was removed before commit.  A corrected functional would need the exact local bilocal/2PI
Legendre kernel, not another symmetry-factor or pole patch.  The preferred route has therefore
reached a genuine new-derivation boundary; production spectra remain blocked and the documented
smooth-`r(h)` prescription is now the only implementation-ready backup route.

Line numbers below are as of the 2026-07-27 working tree; confirm each before editing.

---

## 0. What this plan is, and the one thing it is not

The failure is **thermodynamic state construction**, upstream of the response: `invz_chi_realaxis` has no
lattice loop on the production path (no caller sets `opts.Jfull`; `invz_chi_realaxis.m:57`) and returned
60/60 finite samples at each of the four measured valid-state anchors. The ordered leg fails because the Jensen H_MF
quadrature must evaluate the ω = 0 medium from `h ≈ 0` upward, and on the real coupling multiset that
medium's q-average is taken across denominators that change sign, so the outer Picard map is
meromorphic, pole-sensitive, and non-contractive rather than merely slowly converging.

**The plan does not contain a step that "fixes the physics", because the evidence does not yet identify
which fix is needed.** Everything measured to date shows that the present *damped Picard* method is not
a usable ordered-path solver over the required domain. Nobody has run a branch-tracked solve of a
properly exposed vector residual. That experiment can distinguish numerical non-contraction from
algebraic branch failure, but even a unique finite-grid branch does not by itself prove a physical
off-shell continuation. Phase 1 is therefore a discriminator, not a solver-versus-theory binary fork.

### Invariants every phase must respect

| # | invariant | why |
|---|---|---|
| **I1** | all pre-existing numeric/status outputs of the default `static_medium = 'resummed'` path stay **bit-identical**, except a separately identified non-finite-acceptance defect | additive fields necessarily change whole-struct equality; frozen gate G9 checks the established outputs |
| **I2** | the frozen Gate-0 preregistration is not amended in-run, and no rejected candidate is promoted | prereg §3 forbids in-run fallback |
| **I3** | acceptance rules are not relaxed; a masked column stays masked | §5/D2 — the all-node gate is severe but presently correct |
| **I4** | no pole regularisation (broadening, `ε` in a denominator, a floor on `Dq`) | the pole is a statement, not a numerical nuisance |
| **I5** | every numerical-method change is additive or flag-gated, default off | so I1 holds by construction |

Any task that cannot satisfy I1 by construction carries an explicit G9 comparison in its acceptance
criteria.

---

## Phase 0 — instrumentation and a computability map

**Goal:** convert "masked column, no information" into "this node, this reason", and produce the one
artifact with immediate scientific value. **Touches no theory or preregistration; pre-existing finite
default outputs remain unchanged.** The sole intended behavioural correction is fail-closed handling of
a partially non-finite residual in 0.3. Estimated ~1–2 days.

### 0.1 Failure ledger and node census (D2 / D3 / R6) — one implementation task

* **Files:** `invz_projected/invz_solve_point_ordered.m` (jensen early return),
  `invz_projected/invz_hmf_ordered.m`, and `invz_projected/invz_hmf_status.m`.
* **Change:** attach `pt.hmf_prof = hprof` on the `hstar = NaN` return and change the reducer to
  `[status,status_detail]` while preserving its first output exactly. This is one patch, not separate
  "export" and "census" work.
* Extend the existing fixed node record with `h`, outer-iteration count, raw static residual, and the
  four scalar A–D residuals, populated from the already-available `info`; do not create another node
  schema. On a non-`ok` exit the reducer
  filters `[predictor,nodes]` into a compact `status_detail` containing only the fields needed to explain
  the verdict: `h`, accepted, outer iterations, A–D residuals, raw static residual, `medium_status`,
  `Dq_min`, `min(abs(Dq))`, `D_uni`, `ref_denom`, predictor flag, and bracket membership. Return
  `status_detail = []` on `ok`, so successful untraced calls allocate no ledger.
* Per-outer-iteration records remain behind `opts.trace`; extend the existing `append_iter` record in
  place in 0.2 rather than creating another tracer.
* Define a deterministic **binding node** using the same precedence as the status reducer:
  first degenerate doublet, then first medium-domain event, then first non-finite failed residual,
  otherwise the largest normalized failed residual (stable evaluation-order tie break). On the early
  return, copy that node's `medium_status`, `ref_denom`, and `ref_margin` to the existing point-level
  fields. Do not use an undefined "worst node" rule.
* **Acceptance/I1:** at `T = 0.31 K`, `B = 1 T`, jensen, resummed, the point remains
  `'node_failed'`; `prof` has 33 profile-node records, while `status_detail` includes the separate
  `h = 0` predictor and reports 34 total records, 6 accepted, and `predictor_failed = true`. Its
  failing-node list reproduces `docs/diagnostics/claude_convg_2026-07-27/diag2_B1.log`. Additive fields
  do not change the verdict or gate. Confirm that successful untraced calls allocate no failure ledger,
  and measure failed-call memory/time separately from G9 (which checks results, not performance).
* **Explicitly not in scope:** relaxing the gate, skipping nodes, or interpolating `r(h)` across them.
  That was R3 and it is withdrawn — the failed nodes sit exactly where the path crosses poles or may
  change branch, so interpolation manufactures the Jensen integral; and the `r(0) = 1 + Σ(0)`
  predictor identity needs the very `Σ(0)` that failed to converge.
### 0.2 Pole-proximity instrumentation (Step 1a)

* **Primary file:** `invz_projected/invz_ordered_node_solve.m`. The node solver already has `G0`,
  `Sigma`, `K0`, `Gstat`, and the coupling vector at each outer iteration, so the EMT leaf interfaces
  need not be expanded just to duplicate these diagnostics.
* **Change:** extend the existing `append_iter` record; do not build a parallel diagnostic path. Compute
  each denominator vector once per traced iteration and record
  * dynamic sector — `x = −(1+Σ(1))/G0(1)`, `x − max_q J_ν`, `min(abs(1+Σ(1)+J_νG0(1)))`,
    and the count of nonpositive denominators;
  * static sector — `y = K0 − 1/Gstat`, `min_q Dq`, `min_q |Dq|` (reuse the existing closest-positive
    and closest-nonpositive values), `D_uni`, the count of
    nonpositive `Dq`, the interval/rank of `y` in the sorted coupling multiset, `Gstat`,
    `gstat_local_denom`, and `xi`;
  * whether the inner static closure closed while the outer residual did not (the discriminator
    preserved in the consolidated diagnosis §9.3/§9.5; measured at 3.6 T: static 4.9e−11, outer
    4.0e−2).
* **Store only under `opts.trace`**; the untraced path remains untouched. Cache coupling extrema/order
  once per traced node rather than sorting inside the outer loop.
* **Acceptance:** when `Gstat < 0`, reproduce diag11's identity
  `min_q Dq − D_uni = (J0eff − max_q J_ν)·|Gstat|` to ≤ 1e−12 relative at the 4.0 T jensen endpoint
  (measured 0.0854417 vs 0.0854417). The sign qualification matters: for `Gstat > 0`, the minimum
  is controlled by `min_q J_ν`, not `max_q J_ν`. Here `J0eff` is the separately supplied effective
  uniform coupling, not a `q = 0` row of the sampled `J_ν` array.
* **I1:** untraced path allocation-free and numerically untouched.

### 0.3 NaN-propagating residuals (D5)

* **Files:** there are **four**, not three, projected-path copies:
  `invz_ordered_node_solve.m:234`, `invz_solve_point.m:282`,
  `invz_solve_point.m:389` (Tier-2 helper), and `invz_solve_point_ordered.m:328`.
* MATLAB R2025a was checked directly: `max([NaN 2])` returns `2`, while
  `max([NaN 2],[],'includemissing')` returns `NaN`. Replace all four copies with one shared finite
  max-absolute-difference helper. Move the checker's local `robust_max_abs` to that helper too, so the
  solver gates and checker cannot drift. The helper should reject both NaN and Inf.
* **Risk, and how to handle it:** this *can* change acceptance — but only where a partially non-finite Σ
  is currently being declared converged, which is the defect. **If the G9 comparison shows any
  default-path result changing, stop and report the differing point as a finding** rather than
  proceeding; a difference means a state was previously accepted with non-finite Σ and that is a result
  in its own right.
* **Acceptance:** a focused helper guard covers finite, one-NaN, one-Inf, and shape-mismatch inputs;
  established finite G9 anchors remain bit-identical. If an old point changes only because its iterate
  was partially non-finite, report that point as the defect this task was intended to expose.

### 0.4 h-grid re-grading (R2) — deferred until after Phase 1

Do **not** implement the proposed `'root_scaled'` grid in Phase 0. `F(h)` contains
`∫_0^h r(h')dh'`, so a coarse root bracket cannot be obtained without first computing the very
low-`h` path the proposal hoped to skip. Refining only inside a bracket would also leave the integral's
quadrature error uncontrolled. The continuation experiment should choose steps from conditioning, not
from the production quadrature grid.

After a branch prescription exists, design an error-controlled adaptive quadrature over the full
`[0,h*]` interval and compare its integral/root against the legacy geometric grid. Keep the current
three-argument `invz_hmf_grid` untouched until then.

### 0.5 Deliverable — the `(T, B)` computability map, staged

* **New diagnostic driver** under `docs/diagnostics/` only if the map is intended as a retained
  scientific artifact. Ad-hoc guards belong under `/tmp`, in keeping with the cleaned worktree.
* **Stage A, on the critical path:** emit the four 0.31 K anchors already required by G9 plus the
  1.00 K/2.9 T positive control. This validates the schema without a long sweep.
* **Stage B, after the Phase-1 pilot:** map
  `T ∈ {0.10, 0.31, 0.60, 0.90, 1.20, 1.50} K`. Start on a 0.4 T field grid from 0.5 T to above
  `B_c(T)` and refine every status-transition or non-monotone interval to 0.2 T. Do not assume the
  availability set is monotone. Use the measured 4.709 / 4.021 / 3.562 / 3.192 / 2.717 / 2.001 T
  values (§9.4 vii) as placement guides, labelled as interpolated estimates unless this run refines
  them. This gives the same declared transition resolution without paying for a full Cartesian sweep
  in featureless regions.
* **Record per point:** PM convergence/stability, `x`, `x − max_q J_ν`, Jensen status, accepted-node
  census, pole-distance metrics, and the same under `strict_1z_dyson_ref` as a labelled comparator.
  `path_omit_max` is meaningful only for the strict moment closure and is expected to be NaN on the
  resummed path.
* **Output:** a **scheme/grid-specific numerical availability map**, not a physical phase map:
  {stable PM available} / {ordered Jensen state available} / {ordered Jensen state unavailable} /
  {no accepted state}, plus the raw table at full precision.
* **Why this matters now:** §9.4(viii) already shows the ordered leg solving at 2.9–3.0 T at
  `T = 1.00 K` while solving nowhere on that field set at `T = 0.10 K`. That boundary has never been
  mapped. It tells you where this implementation can currently construct an ordered state; it does not
  validate that state as a controlled ordered 1/z thermodynamics.

---

## Phase 1 — branch-tracked full-residual experiment

**Goal:** determine whether the existing Jensen equations have a reproducible real continuation along the
H_MF path at fields where Picard fails, as specified in the consolidated diagnosis §10.2. It can
establish the algebraic branch structure of the discretized equations; it cannot, by numerical
convergence alone, establish the physical off-shell prescription.

**Entirely a diagnostic script — no production file is modified, so there is no G9 or preregistration
exposure.** Build the first guard in `/tmp`; retain it under `docs/diagnostics/` only if its outputs
become part of the scientific record.

The first guards did become part of the scientific record and are retained in
`docs/diagnostics/invzp_solver_stability_2026-07-27/`. The implementation deliberately stops at a
fixed-node corrector. One-off scaled-tangent and pseudo-arclength drivers remained under `/tmp`; they
established the 1.5 T fold described above but are not a production solver. Temporary end-to-end
wiring was used once to establish the 3.6 T numerical result above, then removed because it had no
branch-identity gate.

### 1.1 Expose a square vector residual

At fixed `h`, use unknowns `u = [Σ(1..n_w); K0]`. Do **not** pass
`invz_ordered_residual` to Newton: it returns four scalar norms/pass flags, not a vector with
`n_w+1` components, and Block A internally invokes the nested static iteration.

This experiment is specifically for the existing **resummed** equations. Retain the fixed-`h`
single-ion/two-level domain screen before constructing the residual; strict one-shot media have a
different static equation and are not silently routed through this residual. Under resummed,
`medium_status = 'not_applicable'`, so bypassing `invz_emt_static_ordered` does not drop a hidden
reference-denominator gate.

Define one explicit, non-iterating residual evaluation:

1. Compute the dynamic medium `Kdyn = invz_emt_scalar(G0, Σ, Jnu, eopts).K` and set
   `K = [K0; Kdyn(2:end)]`. Use the existing helper initially; add a dynamic-only variant only if a
   profile shows the discarded slot is material.
2. Derive `λ = invz_lambdas(K, g, wts, beta, [1 2 3])`.
3. Derive `Σmap = invz_sigma_ordered(tl, λ, K, g, beta).Sigma`.
4. Evaluate `[Gstat, go] = invz_gstat_ordered(..., struct('stable_form',true))` **once**, without
   calling `invz_emt_static_ordered`. This is diagnostic-only exact reassociation, not regularisation.
   Use the finite reciprocal coordinate
   `z = 1/Gstat = d0/(G0inel0 + go.xi*G0el0*d0)`, where
   `d0 = go.gstat_local_denom`; this expression remains finite at the removable `d0 = 0`,
   `Gstat = ±Inf` crossing. For the load-bearing diagnostic values, evaluate
   `Gtil0 = 1/(z-K0)` and `r = go.G0bare*(z-K0)` directly from this same `z`; away from singular
   limits they must agree with `go.Gtil0`/`go.r` to floating-point tolerance.
5. With `E_q = z + Jnu-K0`, choose a common scale
   `s = max(abs(z),Jscale)`, `Jscale = max(abs(Jnu(:)))`, and weights `w_q = s/E_q`.
   Then
   `Gbar = mean(w_q)/s`,
   `Jloc = mean(Jnu.*w_q)/mean(w_q)`, and
   `RK = (K0-Jloc)/Jscale`.
   Use the analytic `z -> ±Inf` limit `Jloc -> mean(Jnu)`, `Gbar -> 0`; require finite positive
   `Jscale`. Return `RΣ = Σmap-Σ` and `RK`.

This is the original defactored EMT equation
`K0 = mean(Jnu.*Gq)/mean(Gq)`, evaluated with
`Gq = 1/(z+Jnu-K0)`. It is algebraically identical to the previous `P/Q` form away from limits, but
has one stable chart through `Gstat = ±Inf` and an explicit limit at `Gstat = 0`. At `Gstat -> 0`,
`RK -> (K0-mean(Jnu))/Jscale` instead of vanishing for every `K0`. Do not use `Q-1` or
`mean(Gq)-Gstat` as the load-bearing residual: both contain the spurious `Gstat = 0` factor.

Preregister two dimensionless event measures:

* `pole_margin = min(abs(E_q))/s` approaching zero — a genuine lattice denominator pole;
* `mean_margin = abs(Gbar)*Jscale` approaching zero — at finite `z`, a genuine cancellation in the
  denominator of the defactored K closure; as `abs(z)/Jscale -> Inf`, the separately identifiable
  `Gstat -> 0` boundary with an unbounded Jensen integrand. The analytic `Jloc -> mean(Jnu)` limit
  keeps the residual evaluable in the latter case, but does not make that thermodynamic path point
  admissible.

Keep `Q` only as a reported coordinate. The useful precision monitor is
`q_cancel = abs(mean(w_q))/mean(abs(w_q))`, which is invariant to the common scale and does not
misclassify the traversable `z = 0` crossing. If `Q` is reconstructed, use `Q = z*Gbar` with its
analytic `z -> ±Inf` limit rather than forming the original `Dq`. Require finite stable `Gtil0` and
`r`; allow `Gstat = ±Inf` when `z`, `Gtil0`, `r`, and the residual remain finite. Thresholds classify
and halt unresolved events; they are never inserted into an equation.

After a corrector claims convergence, construct the complete state and run
`invz_ordered_residual` as an **independent acceptance audit**. Its A–D blocks must all pass. This uses
the existing checker for the job it was designed to do without making its scalar norms into Newton
equations. The `mean(Gq)` event gate, finite-`r` gate, and defactored `RK` remain independently binding
because the current resummed Block B uses the factored `mean(Gq)-Gstat` residual and can accept
`Gstat = 0`.

At an exact `d0 = 0` crossing the production checker may produce NaN because its resummed arithmetic
does not use `stable_form`. Do not modify it for this experiment and do not waive A–D. Instead require
accepted correctors on both sides whose residuals and stable `r` converge to the same limit; the
coordinate-singular point itself is recorded as a traversed limit, not counted as an accepted node.

**Jacobian implementation.** On a static coupling vector, each dynamic slot is
`K_n = kappa_n(Sigma_n)`. Its derivative is analytic: the existing q-loop gives
`A'_n = -mean(Jnu./(D_n + Jnu*G0_n).^2)`, after which `kappa'_n` follows by the quotient rule. Since
the only cross-frequency quantities are `lambda_1`, `lambda_2`, and `lambda_3`, the `RΣ` block is a
diagonal matrix plus rank at most three. The `RK` row is **dense**, not sparse:
`Gstat` depends on `lambda_1` and `lambda_2` through `xi`, so every dynamic slot contributes through
those two sums. Even so, with the dense `K0` column and bordered `RK` row included, the complete
Jacobian is diagonal plus a rank-at-most-five correction. Differentiate the stable static row:
differentiate the closed quotient for `z`, then
`dE_q = dz-dK0`,
`dGbar = -mean(dE_q./E_q.^2)`,
`dN = -mean(Jnu.*dE_q./E_q.^2)`, and
`d(K0-N/Gbar) = dK0-(Gbar*dN-N*dGbar)/Gbar^2`, where
`N = mean(Jnu./E_q)`. These formulas stay finite at `z = 0`.

Stage the linear algebra to get an informative result before optimizing:

1. Bring up the residual and corrector on the small `T = 1.00 K`, 2.9 T positive control
   (`n_w ≈ 75`) with a centred-finite-difference dense Jacobian. It is slow per step but cheap enough
   here, and separates residual/driver errors from analytic-derivative errors.
2. Assemble the analytic Jacobian densely from the formulas above, compare its entries and Newton step
   with the finite-difference oracle, then use MATLAB's dense solve for the first 0.31 K edge run
   (`n_w ≈ 240`).
3. Only if profiling says the dense linear solve is material, or before extending the experiment to
   `T = 0.10 K` (`n_w ≈ 740`), expose those same terms as diagonal-plus-low-rank factors and add the
   Woodbury/bordered solve. Keep the dense path as the oracle/fallback; do not write another residual.

The optional four-variable reduction is **not** on the initial implementation path. Add it only if the
full solve exposes multiple branches that need per-frequency scalar-root enumeration; it is
branch-conditioned and fails at implicit scalar folds.

### 1.2 Continuation driver

1. Start at the largest `h` node where Picard converges cleanly (measured at 1 T: nodes 28–33, 13–20
   iterations, `D_uni`/`Dq_min` ≈ 0.8–0.99).
2. Implement adaptive natural-parameter continuation downward in `h` first, using one damped Newton
   corrector and step-retry policy. Add pseudo-arclength only if the fixed-`h` Jacobian loses rank while
   the augmented tangent remains regular or the computed tangent has `dh/ds` approaching zero. Then
   augment the same residual on `[u;h]` with one bordered arclength equation; a tangent predictor plus
   a fixed-`h` corrector is not pseudo-arclength and cannot pass a fold.
3. Record `x`, stable `z` and `y = K0-z`, `pole_margin`, signed `Dq_min`, `D_uni`, the coupling interval
   containing `y`, `Gbar`, `mean_margin`, `q_cancel`, reported `Q`, `Gstat`,
   `gstat_local_denom`, `xi`, `Gtil0`, `r`, residual norms, step size, corrector iterations, and a
   smallest-singular-value/condition estimate. Compute legacy `Dq` metrics only where their arithmetic
   is finite; `pole_margin` is the load-bearing pole measure.
4. Halt on `pole_margin` reaching its declared threshold; `mean_margin` reaching its declared threshold
   (classified as finite-`z` cancellation or `Gstat -> 0`/unbounded `r`); rank loss of an already
   augmented system; corrector step collapse after retries; or a branch-identity jump not resolved by
   step reduction. Do **not** halt merely because `Q -> 0`, `Gstat -> ±Inf`, or `z -> 0`: the exact
   reciprocal chart is designed to traverse that local-denominator crossing.
5. Use the staged Jacobian policy in §1.1. Finite differences are a positive-control oracle only;
   analytic dense solves are acceptable for the initial 0.31 K campaign, and Woodbury is a profiled
   optimization or low-temperature prerequisite rather than an up-front deliverable.

### 1.3 Reproducibility protocol, applied as a funnel

The positive-control and first edge pilots use the registered 16³/`Ecut` configuration. They validate
implementation only and support no domain claim. Once a pilot produces a candidate branch, promote
only that branch and its event/end points through:

* decreasing-`h` and increasing-`h` traces, with the reverse trace seeded independently where possible;
  merely retracing from the final forward solution tests determinism, not branch independence;
* cold correctors at registered checkpoints versus warm continuation;
* at least **three** coupling grids (12³, 16³, and 20³) for a convergence trend. Two grids can show a
  difference, not separate branch structure from the Γ-exclusion artifact;
* the registered `Ecut` plus one Matsubara-cutoff refinement;
* tolerances declared **before** the run for the scaled defactored vector residual, A–D audit, `m0`,
  `Σ(0)`, `K0`, `hstar`, the Jensen integral, `pole_margin`, `mean_margin`, `q_cancel`, and branch
  metrics.

Run that protocol at production size before making a continuation-domain claim. Do not take the
Cartesian product of every exploratory field, direction, grid, and cutoff: refine the representative
branch plus each distinct fold/pole/cancellation end point, and expand only if different pilot fields
show different branch topology.

### 1.4 Field funnel

1. `T = 1.00 K`, 2.9 T — positive control where resummed Picard already solves.
2. `T = 0.31 K`, 3.6 T — first failed edge case and cheapest informative branch test.
3. `T = 0.31 K`, 1.0 T — deep-ordered stress test only after the edge corrector is validated.
4. Expand to 0.31 K at 2.0/3.0 T and 1.00 K at 2.5 T only to map a branch-topology change or support a
   promoted domain claim. All four 0.31 K fields are `node_failed` today; 2.9 T is the positive control.

### 1.5 The fork

| outcome | reading | next |
|---|---|---|
| defactored-residual- and A–D-accepted, finite-Jensen-integrand, direction-independent branch with a grid/cutoff limit | strong evidence that Picard non-contraction is a numerical defect of the current solver | keep continuation diagnostic/flagged; then justify the off-shell branch physically before any default change |
| folds or multiple reproducible real branches | the algebraic equations need a branch-selection prescription | common-functional work, or a separately preregistered continuation prescription |
| no branch within the declared search, with reproducible pole/rank-loss events | evidence that the chosen **real finite-grid continuation** fails | report its domain; this does not prove that every analytic/PV/complex continuation or the H_MF representation is impossible |
| the 2.9 T positive control is not reproduced | the continuation implementation or residual is invalid | fix the experiment before drawing a physics conclusion |

**Phase 2′ (only after the first outcome and an explicit branch prescription).** Add
`opts.ordered_solver = 'picard' | 'continuation'`, default `'picard'` (I1), initially as an experimental
path. Establish equivalence at every point where Picard succeeds, full `[0,h*]` quadrature convergence,
and grid/cutoff convergence where Picard fails. A default flip is a separate reviewed and preregistered
production decision; finite-grid uniqueness alone is not sufficient.

---

## Phase 2 — theory candidates (separate research tracks)

These candidates are not logically gated only on a non-unique Phase-1 branch. They change either the
closure or the coupling measure and therefore remain theory candidates even if the finite-grid
continuation is unique. **One** fully specified candidate may be preregistered and evaluated at a time;
I2 forbids switching between them as an in-flight fallback.

### Candidate A — replace the strict reference `Gref` (motivated, not yet specified)

§9.4(viii) isolates the one temperature-robust Gate-0 failure: clause (a), `ref_denom_nonpositive`,
fires at 0.05–1 T at `T = 0.10 K` and 0.05–2.0 T at `T = 1.00 K`. Clause (c) is inside its frozen gate
at T ≈ 1 K (omit 0.058–0.075) and clause (e)'s PM half fired on controls placed 1.2–1.6 T *inside* the
ordered phase. Thus one important failure is localized to `Gref = G0bare/(1+Σ(0))`, chosen
for its `m → 0` cross-leg identity and **untested at large `m`** — which is the record's own §11 item 1.

"Build the reference from the ordered local propagator" is not yet an implementable candidate. Before
code, specify the formula, order counting, domain, and whether it depends on `K0`, `λ`, or `xi`.
Dependence on those quantities can make the supposedly one-shot closure self-consistent again and
reintroduce the same denominator feedback. Derive the candidate from the ordered expansion, prove the
`m → 0` cross-leg identity, and define its omitted-order estimator before preregistration.

### Candidate B — density / continuation prescription (Step 1b)

Replace the empirical measure by a smooth `ρ`, evaluated inside the support by singularity subtraction,
`PV∫ρ(J)/(J−x)dJ = ∫[ρ(J)−ρ(x)]/(J−x)dJ + ρ(x)·ln|(J_max−x)/(x−J_min)|`, with `S(x ∓ i0) = PV ∓ iπρ(x)`
as the one-sided alternative.

**This is a new approximation, not a stabler evaluation.** For the real 16384-value multiset `S_N` is a
rational function with 16384 simple poles — meromorphic, not cut — so smoothing changes the medium by an
amount nobody has bounded. It must not be sold as an implementation change, and `−πρ(x)` must not be
described as a decay rate without a dynamical derivation (`ρ` is a density of couplings over **q**).

### Preregistration requirements for whichever is chosen

1. Register a **`(T, B)` domain**, not a single cut. The frozen fixture is the worst case: §9.4(vii)
   shows its field set is a `T ≈ 1 K` set evaluated at `T = 0.10 K`.
2. Place ordered fields and PM controls against **measured `B_c(T)`** (§9.4 vii), not against an
   inherited constant.
3. For a candidate that retains the **same coupling-moment expansion**, gate its own omitted terms at
   the registered temperature and every required path node. The present `omit_max` is tied to
   `K0 = Jbar − μ2Gref` and its specific `μ3/μ4` remainder; it is not automatically the right gate for
   a density prescription or a differently derived functional.
4. Keep the two-tier verdict (`converged` vs `stable_1z`); do not collapse them.
5. State the m → 0 cross-leg identity as a pass/fail gate.

**Prerequisite decision (yours):** whether to amend the frozen fixture. That is a dated amendment, not
something to fold into a candidate spec.

---

## Phase 3 — the durable formulation

The consolidated diagnosis §10.2's scientifically durable route, regardless of the Phase-1 outcome:
derive the ordered thermodynamics **and** the response from one truncated free-energy/effective-action
functional, with χ from the stationary Hessian. The ordered state then comes from stationarity rather
than from integrating a separately resummed susceptibility across an undefined unstable branch.
Multiple stationary solutions and spinodals may still require numerical continuation, but the
functional supplies the missing thermodynamic comparison/selection principle.

The 2026-07-28 focused audit in `invzp_functional_derivation_audit.md` found no already implied common
functional or thermodynamic branch rule in the current hybrid. It did identify the branch-conditioned
`Phi_path` and an isolated scalar-EMT log potential; both are useful checks for WP0, but neither
authorizes selecting the present multiple roots.

Their checklist, adopted unchanged: strict retained-order vacuum / one-point / two-point diagrams;
validation against a small exact finite-cluster oracle; analytic handling of the elastic zero-frequency
sector; resummation only where it follows from stationarity of the same functional; PM/FM boundary
closure on the real coupling density. Entry point: `rigorous_z^1_extension_Codex.md`.

One constraint to carry into the spec: every approximation needs an omitted-order/error diagnostic
derived from **its own** expansion. The current `omit_max` remains valid only for the existing
coupling-moment closure and reference convention; it cannot be inherited by a different functional
merely because both formulas contain powers of a propagator. A functional derivation supplies a
consistent construction and branch-selection criterion, not by itself a lattice sum that converges
where the present one does not.

---

## Interim production posture (in force throughout)

* Keep `jensen` columns masked where the path is incomplete. Do not relax (I3).
* Keep `ordered_mode = 'bare'` available as an **explicitly labelled** bare-H_MF moment-form 1/z
  response — never substituted into a column requested as `jensen`. It converged across the measured
  real-field anchors because it never evaluates the Jensen medium along the unstable path, and its
  moment onsets at the mean-field boundary, so it is a labelled product and not an answer.
* Report per ordered column: dynamic support margin, `min|Dq|`, static branch interval,
  scheme-appropriate omitted-order diagnostics, `Σ(0)`, and the node census.
* Do not promote either strict scheme (I2).

---

## Explicitly excluded, with the reason

| excluded | why |
|---|---|
| raising `max_outer`, softening `mix_outer`, widening tolerances as the production fix | the pole-sensitive map has O(10²–10³) sign-indefinite gain; the tested `mix_outer = 0.02` reaches a deeply unstable PM root (`crit = −3.669`) and still returns `node_failed` on the ordered path after 10,354 iterations |
| skipping failed H_MF nodes or interpolating `r(h)` across them | the failures sit at pole/branch events — interpolation manufactures the integral (withdrawn R3) |
| making the `h = 0` predictor non-fatal via `r(0) = 1 + Σ(0)` | needs the `Σ(0)` that failed to converge |
| regularising the pole | I4 |
| switching dipolar backends | Ewald vs brute-force differ by 1.168e−3 in `norm(sort J)`; both return `node_failed` at 1 T. Ewald is right on its own merits, unrelated to this failure |
| promoting `strict_1z_dyson_ref` on convergence evidence | failed a preregistered gate; refining *why* is not overturning it |
| treating D4 (damping Σ but not λ/K0s) as a defect to fix | block-Picard with an inner converged solve is legitimate; no evidence separate damping helps. It belongs in Phase 1's coupled solve, not a maintenance list |

---

## Verification checklist (applies to every phase)

1. **G9 bit-identity of pre-existing outputs** on the finite default resummed path — capture a fresh
   baseline and compare before/after each numerical Phase 0 change set. Implement 0.1/0.2 together
   because they extend the same node/trace records; audit the behaviour-changing 0.3 separately. Added
   diagnostic fields are excluded from whole-struct equality.
2. **Coupling provenance:** verify fixture digest `ddb9532d…` on the registered 16³ baseline
   (already enforced); record separate exact digests for 12³/20³ rather than expecting the 16³ digest.
3. **Anchors**, not only the 24-point synthetic fixture: a deep ordered point (1 T), the failure edge
   (3.6 T), a near-boundary accepted point (3.8 T), and a PM point (4.6 T), all at `T = 0.31 K`.
4. **Defactored square vector residual plus the four-block A–D audit, finite-`r` gate,
   `pole_margin`, and `mean_margin`**, not inner static closure or small `Q` alone — a small
   `resid_static` is not evidence of a physical medium (measured 6e−11 with `Dq_min = −155`).
5. **Jacobian oracle:** first reproduce the positive control with centred finite differences; then
   compare the analytic dense Jacobian and Newton step against that oracle. If Woodbury is implemented,
   compare its bordered step with the same dense solve before using it.
6. **Forward/reverse and cold/warm agreement** for anything involving continuation.
7. **At least three coupling grids and one `Ecut` refinement** for any continuation-domain claim.
8. **`(T, B)` reporting**, never a single cut, for any domain claim.

## Practical notes

* MATLAB process name is `MATLAB_maca64` (`pgrep`); long runs must be launched with `nohup … &`, never
  through a pipe — piping `matlab -batch` buffers all output and the run looks dead.
* Write progress to a log file line by line (`fopen`/`fprintf`/`fclose` per line) so a long run can be
  followed; poll with `until grep -q "^done" <log>`.
* Keep one-off fresh guards in `/tmp` and remove them after use. Only retained scientific diagnostics
  belong in `docs/diagnostics/`; they are **not** on the MATLAB path and set explicit `addpath` calls.
  Do not `addpath(genpath(<repo root>))`.

## Effort and sequencing

| phase | scope | effort | gated on |
|---|---|---|---|
| 0 | two instrumentation change sets + staged numerical-availability map; h-grid change deferred | ~1–2 days | nothing — start here |
| 1 | defactored reciprocal-chart residual + staged continuation experiment | estimate the campaign after the 2.9 T positive-control and 3.6 T edge pilots | Phase 0.1–0.3 |
| 2′ | experimental continuation solver behind a flag | estimate only after Phase 1 | grid-limited, A–D-accepted branch plus an explicit prescription |
| 2 | one fully specified theory candidate + preregistration | weeks | derivation and candidate selection, not merely a Phase-1 outcome |
| 3 | common-functional formulation | research programme | durable track regardless of Phase 1 |

## Recommended execution order

1. Implement 0.1 and 0.2 as one traced-diagnostics patch, capture Stage-A output, and run one G9
   comparison. Implement 0.3 as a separate fail-closed patch and run its focused guard plus G9.
2. Build the Phase-1 residual and finite-difference positive control at 1.00 K/2.9 T. Add and validate
   the analytic dense Jacobian, then run the 0.31 K/3.6 T edge pilot. Do not build pseudo-arclength,
   Woodbury, the four-variable reduction, the full field sweep, or the three-grid campaign until its
   corresponding trigger above occurs.
3. Before pseudo-arclength, use the retained corrector only at isolated Picard-failed nodes and require
   step-reduction branch continuity. This already repairs the 0.10 K/3.6 T path, but not 1.5 T.
4. The scaled tangent predictor has now confirmed the fixed-`h` rank/tangent trigger at 0.10 K/1.5 T,
   and pseudo-arclength has confirmed that the high-`h` branch turns toward increasing `h`. The
   independent low-`h` trace also terminates in a regular fold, far below the high-branch fold, and
   turns toward decreasing `h`. This enters the second row of the Phase-1 fork. Do not run the
   expensive reproducibility funnel, wire a fallback into production, or consider Phase 2′ until an
   explicit branch prescription exists.

Two scientific decisions remain external to this implementation sequence: Candidate A still needs a
derived formula before it can be compared with Candidate B, and any amendment of the frozen Gate-0
fixture must be dated explicitly. The evidence in consolidated diagnosis §9.4(vii)/(viii) informs but
does not make either decision; the existing FAIL verdict remains in force.
