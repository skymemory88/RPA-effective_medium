# Biased convergence solution — explicit smooth-`r(h)` branch prescription

**Recorded:** 2026-07-28

**Status:** decision and full-state continuation oracles added after the first common-functional
skeleton candidate failed its immutable weak-coupling gate; non-default and not production-integrated

**Trigger:** use only if the common-functional route proves impractical on the required timescale and
the prescription below is separately preregistered before spectra are inspected

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

Use the lexicographic smoothness criterion:

1. reject every path with a finite jump in `r`;
2. among the remaining paths, minimize the maximum resolved slope change;
3. then minimize the integrated curvature

\[
\mathcal C_r =
\int_0^1
\left(\frac{d^2r}{dx^2}\right)^2 dx .
\]

For discrete evaluation, use derivatives of the branch-tracked representation, not a smoothing spline
that invents values between unresolved roots. The finite-difference stencil, endpoint treatment,
adaptive refinement rule, and numerical tie tolerance must be preregistered. Report both the
continuum-scaled value and its mesh-refinement drift.

Because different candidates can have different `hstar_path`, also report the same diagnostics in
physical `h` units. The dimensionless `x` criterion remains the preregistered selector, but if its
winner changes under an admissible normalization or within the refinement uncertainty, return
`branch_ambiguous` rather than exploiting the candidate-dependent interval length.

### 3.4 Tie breakers

If two complete paths are indistinguishable in `r` smoothness within the frozen numerical uncertainty,
apply, in order:

1. smaller maximum scaled state curvature;
2. closer continuous approach to the paramagnetic endpoint over the registered near-QCP field set;
3. agreement of forward, reverse, and independently cold-seeded traces;
4. otherwise return `branch_ambiguous`.

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
3. Cluster roots using the scaled full state, not `r` alone.
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
fold and reproduced its field to `1.44e-10 meV`; every A--D audit passed. The global cold-seed root
enumerator and single-valued path-construction fixture remain to be built before step 1 as a whole is
complete. The dense low-temperature trace cost 1292.9 s. Profiling showed that local node and
Richardson field-derivative evaluation, not the dense bordered solve, dominates this cost; exact
last-point caching is now active, while an exact low-rank solve is retained only as an oracle.

## 6. Failure conditions

Keep the column masked if any of the following occurs:

- no single-valued continuous accepted section covers `[0,hstar_path]`;
- a candidate has no unique admissible increasing `F_path=0` endpoint;
- the selected path changes under grid/cutoff refinement beyond its error budget;
- the smoothness winner is tied within uncertainty;
- the selected state is discontinuous as `B -> Bc^-`;
- the Jensen integral is not quadrature-converged;
- the rule selects different paths under forward and reverse reconstruction;
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
- `C_r`, state-continuity metrics, tie margins, and refinement drift;
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
pre-registered nonlocal-return skeleton candidate.  The strict ring theory passed its gates but
retained a Gaussian no-state interval; the stationary skeleton failed its exact mixed-chain
coefficient because it lacked the exact local bilocal Legendre kernel.  That candidate was rejected
rather than tuned.  The trigger for this backup is therefore satisfied at the implementation-oracle
level, while all non-default, branch-identity, refinement, and production-stop rules above remain
binding.

Relevant current records:

- `invzp_convg_fix_Claude.md`
- `invzp_convg_diagnosis_Claude.md`
- `docs/diagnostics/invzp_solver_stability_2026-07-27/README.md`
- `rigorous_z^1_extension_Codex.md`
