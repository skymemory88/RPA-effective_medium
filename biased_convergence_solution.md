# Biased convergence solution — explicit smooth-`r(h)` branch prescription

**Recorded:** 2026-07-28

**Status:** backup authorized for diagnostic/oracle development after the first common-functional
skeleton candidate failed its immutable weak-coupling gate; non-default and not production-integrated

**Production trigger:** not yet satisfied. Production use requires both a separately preregistered
branch prescription and a dated, reviewable finding that the preferred common-functional route is
impractical within a prospectively declared scope, effort budget, and deadline. That finding must
record the attempted exact-local-bilocal/2PI deliverables, the blocking evidence, and why no retained
alternative can meet the declared accuracy and domain gates. Rejecting one skeleton candidate
authorizes oracle work only; it does not exhaust that route.

## 1. Purpose and declared bias

The resummed ordered equations can have several finite, A–D-accepted roots at the same auxiliary
Jensen field `h`. At `T = 0.10 K`, `B = 1.5 T`, branch-tracked continuation has already resolved two
regular folds and an intervening root segment. A nonlinear solver can enumerate these roots but cannot
decide which one belongs in the Jensen integral.

This backup makes that missing decision explicit:

> **Among admissible root paths, select the path whose Jensen integrand `r(h)` is globally smoothest,
> subject first to continuity of the state and to the physical endpoint/QCP constraints.**

The motivation is the continuous quantum phase transition: the physical state and its observables are
expected to approach the paramagnetic state continuously as the physical field approaches the QCP.
A branch prescription that creates avoidable discontinuities in `r`, `Sigma`, `K0`, or the moment is
therefore disfavoured.

This is deliberately called a **biased** solution. Continuity across the QCP is continuity with
respect to the physical field or temperature. It does not by itself prove smoothness along the
off-shell auxiliary `h` path used by J 2.31–2.33. Smooth `r(h)` is consequently a transparent modelling
assumption and selection principle, not a theorem derived from the current resummation.

## 2. Non-negotiable admissibility gates

Smoothness never rescues an invalid root. Every candidate state at every `h` must first satisfy:

1. the defactored reciprocal residual tolerance;
2. the independent ordered-state A–D audit;
3. finite `r`, `Gtil0`, `Sigma`, `K0`, and required response quantities;
4. the declared pole and `mean(Gq)` cancellation margins;
5. the two-level/local-domain checks;
6. the selected grid, Matsubara cutoff, and coupling provenance;
7. local branch-identity control: a corrector remains inside its preregistered predictor tube.

Failed nodes are not skipped, interpolated, broadened, or assigned a large but finite surrogate cost.
If no admissible continuous path covers the complete Jensen interval, this backup returns
`branch_unresolved` and the spectrum remains masked.

## 3. Selection hierarchy

The following order is binding. A lower item cannot override a higher one.

### 3.1 Physical endpoint and QCP continuity

Require the selected ordered solution to:

- have its own finite increasing root
  `F_path(hstar_path) = h0_path(hstar_path)-J0*m(hstar_path) = 0`, with
  `crit(hstar_path)>0`;
- approach the accepted paramagnetic solution as `B -> Bc(T)^-`, with preregistered tolerances on
  `m`, `Sigma(0)`, `K0`, `r`, the static mass, and the soft-mode observables;
- preserve the `m -> 0` ordered/paramagnetic identity;
- retain the expected continuous onset of the order parameter.

`hstar_path` is an output of each candidate `r(h)` path, not a common input. Freezing the current
production `hstar` would be circular because `h0_path` is the integral of the branch-dependent `r`.
The QCP comparison is an external acceptance gate applied after each candidate endpoint is solved.

The branch selector must use no measured spectral peak, fitted critical field, or other downstream
empirical agreement as an input. Those remain external validation data.

### 3.2 State continuity

Candidate roots at neighbouring `h` values may be joined only when a branch-tracked corrector or a
pseudo-arclength segment demonstrates local continuity. Define a scaled state distance, for example

\[
d_u^2 =
\frac{1}{n_\omega}\sum_n
\left|\frac{\Delta\Sigma_n}{S_{\Sigma,n}}\right|^2
+\left|\frac{\Delta K_0}{J_{\rm scale}}\right|^2
+\left|\frac{\Delta m}{S_m}\right|^2 ,
\]

with every scale and threshold frozen in the preregistration. An apparent crossing in `r` alone does
not authorize switching between distant states.

At a fold, continue in arclength to enumerate both fixed-`h` roots. A branch switch at a common `h` is
allowed only if the state-distance and tangent-line gates demonstrate a continuous connection.
Otherwise it is a discontinuous jump and is inadmissible, however smooth an interpolation of `r`
would look.

Pseudo-arclength is an enumeration tool, not the Jensen integration coordinate. A viable candidate
must be a single-valued monotone-in-`h` section `u(h)` containing exactly one selected root at every
`h` from zero to its own `hstar_path`. When an arclength curve folds, its repeated `h` values are
distinct branch labels; the loop itself cannot be inserted into `int r(h) dh`. No switch across the
fold is allowed unless a separate accepted connection satisfies the state/tangent gates at common
`h`.

### 3.3 Primary selector: global smoothness of `r(h)`

Enumerate all root segments first on a common preregistered physical-`h` window. Construct every
admissible single-valued section, solve that section's own `hstar_path`, and then let
`x=h/hstar_path` for comparing complete candidates on a common adaptively refined dimensionless mesh.
Selection is global rather than greedy: choosing the nearest or smoothest root at one node must not
pre-empt a smoother complete path.

Raw derivatives on `x` cannot be compared directly. In particular,

\[
\int_0^1 r_{xx}^2\,dx
=(h^\star_{\rm path})^3
\int_0^{h^\star_{\rm path}} r_{hh}^2\,dh ,
\]

so the unnormalised dimensionless curvature systematically favours a shorter candidate; both it and
the raw slope-change metric also scale with the amplitude of `r`. Instead define

\[
D_r=\int_0^1 r_x^2\,dx,\qquad
S_r=\frac{\max_i|\Delta r_x|}{\sqrt{D_r}},\qquad
Q_r=\frac{\int_0^1 r_{xx}^2\,dx}{D_r}.
\]

For nonconstant paths these shape scores are invariant under `r -> a*r+b`, `a != 0`, and do not
prefer one `hstar_path` over another when the candidates have the same normalized shape. This is a
deliberate comparison of shapes on their complete Jensen intervals, not an assertion that `Q_r`
equals the physical-`h` bending energy. The binding lexicographic criterion is:

1. reject every path with a finite jump in `r`;
2. among the remaining paths, minimize `S_r`;
3. then minimize `Q_r`.

The preregistration must freeze numerical floors for `D_r` and the raw curvature integral. If both
are below their compatible constant-path floors, assign `S_r=Q_r=0`. If `D_r` is below its floor but
the curvature test is not, the derivative representation is unresolved and the path cannot be
selected.

For discrete evaluation, use derivatives of the branch-tracked representation, not a smoothing spline
that invents values between unresolved roots. The finite-difference stencil, endpoint treatment,
adaptive refinement rule, and numerical tie tolerance must be preregistered. Report both the
normalized scores, the raw `x` and physical-`h` diagnostics, and their mesh-refinement drift. Raw
metrics are report-only and cannot break a normalized-score tie. If the normalized winner changes
within the refinement uncertainty, return `branch_ambiguous`.

### 3.4 Tie breakers

These tie breakers apply only after every candidate has passed the mandatory reconstruction failures
in section 6; they cannot rescue missing or disagreeing forward/reverse/cold-seed evidence. If two
complete paths are indistinguishable in normalized `r` smoothness within the frozen numerical
uncertainty, apply, in order:

1. smaller maximum scaled state curvature;
2. closer continuous approach to the paramagnetic endpoint over the registered near-QCP field set;
3. otherwise return `branch_ambiguous`.

Do not break a tie using agreement with the desired spectrum.

### 3.5 Existing one-dimensional potential: binding diagnostic, not primary selector

Once an internal branch has already been chosen, J 2.31–2.33 imply

\[
F(h)=h_0(h)-J_0m(h),\qquad
\Phi_{\rm path}(m)=\int_0^m F(m')\,dm',
\]

with `F'(h) = crit(h)`. An accepted spontaneous root must have `crit(hstar_path) > 0`, the
local-minimum condition on that chosen path. Compute and report `Phi_path` for every complete
smooth-`r` candidate.

At present it does not override the primary selector because its normalization has not been closed
against J 2.34, and constructing it already assumes a branch for `r(h)`. It therefore cannot by itself
normalize disconnected internal root components. If a later derivation supplies the missing common
normalization and demonstrates that `Phi_path` compares all candidates thermodynamically, that
physical free-energy comparison supersedes this biased prescription.

Even on a selected branch, `Phi_path` requires `h -> m` to be single-valued and invertible over the
integrated segment (`dm/dh=-G0bare>0` in the present convention) and every compared path to use the
same physical `m=0` reference. If either condition fails, `Phi_path` is reported as undefined rather
than used as a comparison.

## 4. Root enumeration and global path construction

For each preregistered `(T,B)`, first freeze a common search window `[0,h_search_max]` that does not
depend on a candidate endpoint:

1. Construct independent roots from declared cold seed families at `h = 0`, `h = h_search_max`, and
   fixed interior checkpoints.
2. Trace every distinct accepted root in both arclength directions. Stop only on a declared event,
   augmented rank loss, domain exit, or completed connection to an already identified segment.
3. Cluster roots using the scaled independent state `[Sigma;K0]`, not `r` alone. The exported `K`
   and `lambda` are derived coordinates and do not receive a second vote.
4. At common `h`, retain all distinct roots and their tangent lines.
5. Build a graph whose vertices are accepted roots and whose edges are demonstrated local
   continuation segments. Split folds into distinct branch labels and enumerate only
   single-valued-in-`h` sections.
6. On every section covering `h=0`, construct `h0_path`, locate all increasing
   `F_path(h)=0` roots inside the search window, and form one candidate for each admissible
   `hstar_path`. A complete candidate covers `[0,hstar_path]` with exactly one root per `h`; visual
   proximity is not an edge.
7. Evaluate the selection hierarchy on all complete paths. Dynamic programming may be used after the
   curvature cost is written in a local augmented-state form; exhaustive enumeration remains the
   small-fixture oracle.

The result record must preserve every enumerated branch and cost, not only the selected one.

## 5. Preregistration and validation funnel

Freeze before the first production-spectrum comparison:

- the `(T,B)` pilot and validation domains;
- the common `h_search_max` rule, independent of candidate `hstar_path`;
- coupling grids and Γ policy;
- Matsubara cutoffs;
- residual, event, state-distance, tangent, and ambiguity thresholds;
- cold seed families and root-clustering tolerances;
- the `r`-smoothness discretization and tie uncertainty;
- quadrature tolerances for the resulting Jensen integral;
- near-QCP continuity metrics and their field approach sequence.

Run the funnel in this order:

1. **Implementation oracle:** synthetic small coupling set with deliberately constructed folds and a
   known smooth path.
2. **Measured pilot:** exact legacy 16³ coupling at `T = 0.10 K`, `B = 1.5 T`, retaining the already
   observed low-, intermediate-, and high-`h` segments.
3. **Direction and seed audit:** forward/reverse arclength plus independent cold checkpoints.
4. **Numerical limit:** at least 12³, 16³, and 20³ coupling grids and one Matsubara-cutoff refinement.
5. **Physical-field continuity:** a dense field sequence approaching `Bc(T)` from below and matching
   paramagnetic anchors above it.
6. **External validation:** only after selection is frozen, compare spectra and thermodynamic
   observables with experiment.

The pure decision layer of step 1 is frozen in
`docs/diagnostics/biased_smooth_r_2026-07-28/invzp_select_smooth_r_oracle.m`. It accepts only
pre-enumerated complete paths on an identical dimensionless mesh and never differentiates,
interpolates, or repairs them. Synthetic in-memory fixtures pin the ranking and ambiguity semantics.
The separate generic pseudo-arclength primitive crosses both folds of the analytic cubic oracle, and
its invZ adapter retains complete frequency-resolved states. On the frozen 0.10 K, 1.5 T, legacy 16³
fixture it traced 149 accepted roots from the clean high-field endpoint through the known regular
fold; every A--D audit passed. Its fold-coordinate interpolations were later rejected as precision
claims because one is inconsistent with an accepted extremal node; only the tangent sign change and
regular bordered topology remain evidence. A retained explicit-seed
enumerator now records every fixed-`h` solve and clusters the complete independent coordinates
`[Sigma;K0]` with an uncertainty band; 25 seeds at `h=0` produced seven distinct A--D-accepted roots.
A retained fixed-`h` helper also
stages cheap regular segments around local arclength traces. The returning low-fold leg has reached
`h=6.890625e-9 meV`, but its connection to any enumerated `h=0` root remains unresolved, so the
single-valued path-construction fixture is still incomplete. The dense high-fold trace cost
1292.9 s. Profiling showed that local node and Richardson field-derivative evaluation, not the dense
bordered solve, dominates this cost; exact last-point caching is now active, while an exact low-rank
solve is retained only as an oracle. The tracer now also caps attempts per accepted step and equation
evaluations per corrector attempt; hitting either envelope is a retained `budget_exhausted` failure,
not a branch endpoint or an accepted extrapolation.

The first graph layer is now frozen as
`docs/diagnostics/biased_smooth_r_2026-07-28/invzp_assemble_audited_graph.m`. It normalizes only
explicit accepted-root and local-continuation certificates, requires tangent evidence for any
equal-`h` adjacency, permits resolved termination only for preregistered successful statuses, and
retains missing event/termination evidence as nonconnections. It never
matches roots by proximity, infers an endpoint, splits a fold, or claims a complete Jensen
candidate. Consequently it can expose the currently demonstrated components without changing the
unresolved endpoint verdict.

The fixed-`h` trace now has a fail-closed adapter into that graph schema:
`docs/diagnostics/biased_smooth_r_2026-07-28/invzp_fixed_h_trace_graph_inputs.m`. New Newton records
embed their normalized thresholds; legacy records are labelled as caller-declared unless strict
embedded provenance is required. The adapter accepts no bare event Boolean: a passing certificate
must name a preregistered certifier, bind itself to the exact trace and record pair, and provide
positive full-edge lower bounds for every registered signed margin. Applied without such certificates
to the 33-node root-6 opposite-leg trace, it emits 33 vertices but deliberately leaves all 32
adjacencies blocked. Continuous event certification, not schema translation, is now the next missing
edge gate.

The unresolved low endpoint has a diagnostic positive-side coordinate
`q=(h/h_reference)^2`, with
`h_reference=6.890625e-9 meV`, the last accepted h-coordinate handoff, so that handoff is exactly
`q=1`. This reference is only a coordinate scale. The adapter preserves the original equations
pointwise for every `q>0` and has advanced
the returning leg from `q=1` to `q=0.775` with A--D acceptance, whereas the h-coordinate corrector
stagnates because its field-Jacobian error is large relative to the bordered condition. The attempted
`q=0.7625` point exposed a machine-resolution floor solely in the static closure; row scaling did
not remove it. A default-off one-shot scalar `K0` polish, bounded to 4096 physical-`K0` ULPs and
subject to every unchanged full gate, removed that floor and carried the same component to
`q=0.062004970571964718`. The trace then stopped on q-Jacobian refinement disagreement (about 6%
between two individually drift-accepted stencils) with bordered `rcond≈1.05e-12`. This remains only
a conditioning experiment for the positive approach: `q=0` is an explicit unresolved event, because
no finite one-sided stencil proves the absence of a `sqrt(q)` term. Neither the polish nor the
positive trace can supply or select an h=0 endpoint.

A derivative-free fixed-`h` check then removes the q-Jacobian from this conclusion. Row
equilibration raises the handoff Jacobian condition estimate from raw `8.27e-13` to `6.10e-9`
without changing an equation or raw acceptance gate, and carries the component to
`h=3.4035747258282251e-10 meV`. Refinement stops there; a direct zero-field solve fails at residual
`5.077e-2`, and the last positive state is at least `0.312` away from every root in the frozen
seven-root `h=0` census. The backup selector must therefore treat this endpoint as unresolved. It
must not extrapolate the positive component to zero or substitute the closest census root.

This is component-specific, not a universal failure of all ordered roots. A later simple-first
screen found another strict component (root 6) with a regular fold near
`h≈9.54114e-6 meV`; its raw turn estimate has no certified coordinate enclosure. Its returning leg
crosses `h=0`; a local refined bracket and one frozen
fixed-field correction identify census root 6 with cluster distance `1.42e-14`, after a predictor
correction of `7.43e-12`. Both pass the registered `1e-9` same-root threshold. Root 4 also crosses a
separate fold near `5.56814e-6 meV` under the exact defactored audit coordinate. These results add
certified vertices and local fold evidence, but they do not complete the candidate graph: the
root-7 boundary layer remains unresolved, roots 1--5 still have untraced legs, and continuous
signed-event edge certificates are still missing.

A derivative-free reverse check has since removed the root-6 opposite-leg endpoint ambiguity. All
33 decreasing fixed-`h` nodes from the common handoff to zero pass A--D without an event failure;
the largest frozen predictor distance is `0.0247702`, and the endpoint matches census root 6 at
distance `7.92e-14`. This is endpoint and predictor-tube evidence, not a substitute for the
continuous signed-event edge certificate. The root-7 boundary layer, other untraced root legs, and
the complete single-valued selector therefore remain open.

A fixed-`w`, `w=z-K0`, replay of the 32 positive root-6 targets confirms that the opposite leg is
numerically much better conditioned in this coordinate: every target passed the coupled residual,
bordered-rank, event, and independent A--D point gates. This does not satisfy the graph's edge gate.
Finite accepted samples do not prove existence, uniqueness, nonsingularity, or event clearance over
the intervening continuum. The scalar interval calculation is therefore retained only as
conditional clearance: it emits `event_bracket_ok=false` for every adjacency and supplies no
arbitrary caller-chosen bound for the normalized `|z|` pole margin. Production authorization still
requires a validated continuation enclosure (or the separately declared candidate prescription);
the current fixed-`w` result is a solver-coordinate diagnostic.

The implementation policy is consequently staged and simple-first: use fixed-`h` solves on regular
segments, invoke arclength only near a measured rank/tangent turn, use the defactored closure only as
the registered diagnostic continuation coordinate, retain both its provenance and the amplified raw
closure, and checkpoint each root before any summary work. This does not authorize replacing the
frozen production/raw A--D certification. Do not add a general h-coordinate `K0` polish: the existing
4096-ULP proposal was tested at the root-4/root-5 frontiers and correctly declined both as outside
its envelope, while ordinary step shrink accepted the next points.

## 6. Failure conditions

Keep the column masked if any of the following occurs:

- no single-valued continuous accepted section covers `[0,hstar_path]`;
- a candidate has no unique admissible increasing `F_path=0` endpoint;
- the selected path changes under grid/cutoff refinement beyond its error budget;
- the smoothness winner is tied within uncertainty;
- the selected state is discontinuous as `B -> Bc^-`;
- the Jensen integral is not quadrature-converged;
- forward/reverse/cold-seed reconstruction is unavailable or selects different paths;
- an event is crossed by interpolation or regularisation rather than by an accepted continuation;
- the selected path violates a thermodynamic identity required by the current approximation.

## 7. Possible implementation shape

If the prescription passes its preregistration, add it only behind an explicit experimental option,
for example

```matlab
opts.ordered_solver = 'picard' | 'biased_smooth_r';
```

with `picard` remaining the default. Export:

- all branch records and connectivity;
- the selected branch identifiers;
- `D_r`, `S_r`, `Q_r`, raw physical/dimensionless diagnostics, state-continuity metrics, tie
  margins, and refinement drift;
- a machine-readable selection status and failure reason;
- the exact prescription/version digest.

A future default change requires a separate review and decision. The common-functional formulation
remains preferred because it would compare stationary states thermodynamically rather than choosing
one by geometric smoothness.

## 8. Relationship to the preferred route

The preferred staged recommendation remains:

1. audit the current derivation for an existing thermodynamic functional, analytic-continuation rule,
   or stability criterion;
2. implement that criterion if it is genuinely implied by the approximation;
3. otherwise derive the common free-energy/effective-action functional;
4. use this smooth-`r(h)` prescription only as the documented backup.

Stage 1 is complete: the focused audit found useful branch-conditioned and fixed-variable
potentials, but no common functional or binding physical selector. Stage 2 is therefore inapplicable
to the current equations. Stage 3 produced a strict scalar common functional and then a
pre-registered nonlocal-return skeleton candidate. The strict ring theory passed its gates but
retained a Gaussian no-state interval; the stationary skeleton failed its exact mixed-chain
coefficient because it lacked the exact local bilocal Legendre kernel. That candidate was rejected
rather than tuned, but the exact local bilocal/2PI route remains open. The rejected skeleton
authorizes no production fallback and does not establish impracticality of the preferred route. The
backup remains authorized at oracle level; production activation awaits the separately recorded
prospective impracticality finding and every existing preregistration, branch-identity, refinement,
and production-stop gate above.

Relevant current records:

- `invzp_convg_fix_Claude.md`
- `invzp_convg_diagnosis_Claude.md`
- `docs/diagnostics/invzp_solver_stability_2026-07-27/README.md`
- `rigorous_z^1_extension_Codex.md`
