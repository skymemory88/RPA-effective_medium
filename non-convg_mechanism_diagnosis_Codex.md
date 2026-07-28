# Why ordered-state invZ convergence can look stochastic

**Recorded:** 2026-07-28
**Scope:** the thermodynamic state construction upstream of
`invz_projected/invz_run_spectra.m`, with particular attention to the
`T=0.10 K`, `q=[0 0 0]`, `B=0--9 T` visual-spectra use case.

## Executive diagnosis

The failure is not in the real-frequency susceptibility calculation. A field
column is masked before the susceptibility is evaluated because the ordered
Jensen solver cannot construct every state on the auxiliary molecular-field
path needed for its free-energy integral.

The apparently stochastic pattern is the combined effect of four distinct
problems:

1. **Numerical instability of the production iteration.** The damped Picard
   map contains Brillouin-zone averages with many nearby finite-grid poles.
   In the ordered region it is generally not a contraction. Small changes in
   its input can cause a large, sign-changing output.
2. **A bad continuation coordinate.** Even when the coupled algebraic
   equations have regular solutions, their branch can fold with respect to
   the imposed molecular field `h`. A fixed-`h` Jacobian then becomes nearly
   singular and a natural-`h` solver stalls or jumps branches.
3. **A branch-selection problem, not just a root-finding problem.** Several
   fully residual-accepted real roots coexist at the same `h`. Newton or
   pseudo-arclength continuation can find roots, but the equations alone do
   not say which root supplies the physical, continuous Jensen integrand.
4. **A finite-grid coupling edge.** Excluding Γ leaves
   `Jmax<J(0)`. The resulting ordered-availability sliver shrinks
   approximately as `1/N` and the inferred critical field moves with grid
   size. Near-QCP convergence on one grid is therefore not evidence of a
   continuum ordered solution.

The first two layers admit numerical improvements. The third prevents a
numerically converged patch from immediately becoming the default physics,
and the fourth prevents a single-grid preview from becoming a quantitative
QCP claim. The current finite-`16^3` spectrum is useful for visual
inspection; a final result additionally needs branch and grid convergence.

## 1. What actually fails

The production chain is

```text
invz_run_spectra
  -> invz_spectra_map
     -> invz_solve_point_ordered
        -> invz_hmf_ordered
           -> many ordered node solves along h
```

`invz_hmf_ordered` needs the full path from `h=0` to a nonzero endpoint
because the Jensen condition contains a path integral. If any required node
is not accepted, `invz_hmf_status` returns `node_failed`; the ordered state is
unavailable and the corresponding spectrum column is masked.

This all-node rule is severe but mathematically necessary. Filling a failed
interval by interpolation would assign values to the Jensen integral without
knowing whether the state stayed on the same branch, passed a pole, or folded
back. That would manufacture a free energy rather than evaluate one.

When a valid ordered or paramagnetic state is supplied,
`invz_chi_realaxis` produces finite spectra in the measured tests. Thus the
phrase “the susceptibility does not converge” is shorthand for:

> the thermodynamic state required by the susceptibility was not constructed.

## 2. The local numerical mechanism

### 2.1 Pole-sensitive effective-medium averages

The dynamic medium contains denominators of the form

```math
D_{q,n}=1+\Sigma_n+J_qG_{0,n}.
```

The ordered static closure can be written

```math
G_q=\frac{G_{\rm stat}}
          {1+(J_q-K_0)G_{\rm stat}}.
```

Define

```math
z=\frac{1}{G_{\rm stat}},\qquad w=z-K_0.
```

Then the static lattice problem becomes the scalar Stieltjes-type average

```math
G_q=\frac{1}{w+J_q},\qquad
\bar G(w)=\frac{1}{N_q}\sum_q\frac{1}{w+J_q},
```

with

```math
\frac{d\bar G}{dw}
=-\frac{1}{N_q}\sum_q\frac{1}{(w+J_q)^2}.
```

There is a pole at every `w=-Jq`. Between poles `Gbar` is monotone, but close
to a pole its slope is arbitrarily large. When denominators of opposite sign
occur in the same lattice average, large positive and negative terms nearly
cancel. The local coupling

```math
J_{\rm loc}(w)
=\frac{\sum_q J_q/(w+J_q)}
       {\sum_q 1/(w+J_q)}
```

is additionally singular when its denominator is near zero.

This is not ordinary “critical slowing down.” It is a meromorphic numerical
map whose gain and sign can change abruptly as an iterate crosses pole and
mean-cancellation surfaces.

### 2.2 The self-energy amplifies the static-medium sensitivity

The ordered node update feeds the lattice medium into
`invz_sigma_ordered`. In a measured `T=0.31 K`, `B=1 T` example, the local
linear response was approximately

```math
\frac{\partial\Sigma_0^{\rm map}}{\partial K_0}
\simeq 7.8\times10^2\ {\rm meV}.
```

Consequently, a small change in the pole-sensitive `K0` can cause an
order-one or much larger change in `Sigma(0)`. The production outer update

```math
\Sigma^{(k+1)}
=(1-\alpha)\Sigma^{(k)}
 +\alpha\,\Sigma^{\rm map}(\Sigma^{(k)},K_0^{(k)})
```

is then not locally contractive. Smaller damping can slow the excursions and
occasionally reach a root, but it cannot make a pole-crossing map uniformly
smooth or select which root should be used.

This explains the observed alternating signs, large residual excursions, and
iteration counts that hit their budgets. Raising the iteration limit treats
the symptom, not the unstable map.

### 2.3 Why the result can depend on tiny perturbations

No retained solver uses random numbers. Repeating the same serial calculation
with the same inputs should be deterministic. Nevertheless, the output can
*appear* stochastic because:

- several root basins coexist;
- basin boundaries are close to the warm-start trajectory;
- the map has very large derivatives near pole/cancellation surfaces;
- finite line searches and iteration budgets impose sharp accept/reject
  boundaries;
- row imbalance can turn a solvable Newton equation into an apparently
  rank-deficient one;
- last-bit changes from arithmetic order or parallel execution can be
  amplified enough to cross a basin boundary.

This is deterministic sensitive dependence, not evidence of a stochastic
physical state.

There is also no evidence that the order of Matsubara frequencies is the
source. The production Picard iteration updates the complete `Sigma` vector,
and the diagnostic Newton solve treats `[Sigma(:);K0]` as one coupled system.
Frequency blocking in `invz_emt_scalar` is a memory implementation detail,
not frequency-by-frequency continuation.

### 2.4 Direct tests of the two simplest proposed fixes

Two default-off experiments were run on 2026-07-28 and then removed from the
production path after they failed.

First, the ordered static closure was rewritten in the algebraically
equivalent reciprocal/defactored `K0` coordinate. At `T=0.10 K`, `B=1.5 T`,
this accepted `10/33` profile nodes with `max_outer=200`. Raising the limit
to `1000` accepted the same `10/33` nodes, with the same termination reasons;
every accepted node's `r`, `Sigma(0)`, `K0`, `D_uni`, and `Dq_min` was
bit-identical. The extra iterations therefore supplied no missing state.

Second, the full accepted Matsubara `[Sigma;K0]` state at every HMF profile
node was temporarily exported and used only as the seed for the nearest
log-`h` node at the next physical field. No state was interpolated, and all
ordinary A--D acceptance gates remained binding.

- A cold `4.04 T` control accepted `33/33` nodes and returned
  `hstar=0.0148737 meV`.
- A complete `4.05 T` profile also accepted `33/33` nodes and returned
  `hstar=0.0147337 meV`.
- Seeding `4.04 T` from that complete `4.05 T` profile, with a maximum
  matched log-`h` distance of only `0.00559` (about `0.56%`), reduced the
  result to `31/33`; six warm attempts needed the existing cold retry and
  two still failed.
- Seeding `4.00 T` from `4.05 T` likewise left it at `32/33`. Seeding only
  the `h=0` predictor did not repair it either.

Thus the lack of cross-field seeding in the production `parfor` is not the
primary defect. A nearby accepted self-energy can be an excellent geometric
guess yet belong to the wrong attraction basin of the next field's Picard
map. Cross-field warm starting can make convergence worse, so it cannot be
promoted as a general fix or used to bridge an incomplete field. The
temporary ledger interface was removed; only this result is retained.

## 3. Why moderate field can fail while points near the QCP converge

The observation is counter-intuitive only if one equates the nonlinear
iteration with physical relaxation. They are different:

- physical critical slowing concerns the dynamics or response of the system
  near a continuous transition;
- Picard/Newton convergence concerns the geometry and conditioning of the
  chosen algebraic map and coordinates.

Several effects make the near-QCP ordered calculation accidentally easier:

1. The ordered moment and `hstar` are small, so the Jensen path is short.
2. The Brillouin-zone sum omits the exact Gamma point. Therefore the largest
   retained `Jq` is below `J(0)`, leaving a finite-grid sliver below the
   uniform instability in which every sampled denominator can remain on the
   same side of its pole.
3. The warm start may remain inside the attraction basin continuously
   connected to the paramagnetic solution.

At moderate field the endpoint can itself be regular, but the longer Jensen
path must still begin at `h=0` and pass through off-shell states. Along that
path it can encounter pole-sensitive averages, multiple roots, and folds.
Thus “near the QCP works, 1.5 T fails” is evidence against a simple
critical-slowing explanation and in favour of path geometry plus numerical
stability.

The near-QCP success is not, by itself, a continuum guarantee. The protective
Gamma-exclusion gap shrinks as the Brillouin-zone grid is refined.

## 4. From unstable iteration to multiple branches

### 4.1 A converged residual is not a branch label

The diagnostic formulation solves the coupled residual

```math
R(\Sigma,K_0;h)=0
```

directly. This avoids the noncontractive production Picard map and repairs
some failed nodes rapidly. At `T=0.10 K`:

- at `3.6 T`, Picard accepted `30/33` profile nodes and the coupled Newton
  corrector repaired the three missing nodes;
- at `1.5 T`, Picard accepted only `11/33`, and isolated fixed-node repair was
  insufficient.

The difference shows that part of the problem is a simple solver defect, but
not all of it.

At `1.5 T`, a frozen 25-seed census at `h=0` produced 16 accepted attempts
clustered into **seven distinct full `[Sigma;K0]` roots**. All satisfy the
same residual and independent A--D audits. Therefore:

```text
small residual != unique state != physical branch
```

A Newton fallback without branch history can silently replace a failed Picard
node with a different algebraic root.

### 4.2 Fixed-`h` folds

On a solution curve, the fixed-`h` Jacobian

```math
J_u=\frac{\partial R}{\partial(\Sigma,K_0)}
```

can become singular even while the full augmented derivative

```math
\left[J_u\ \ \frac{\partial R}{\partial h}\right]
```

remains regular. That is a fold in the projection onto `h`, not necessarily
an end of the solution curve.

The retained `1.5 T` traces measured this behaviour at both high and low
`h`: fixed-`h` conditioning collapsed, whereas bordered
pseudo-arclength continuation crossed the fold with accepted residual and
event audits. The tangent component `dh/ds` changed sign.

This creates a mathematical obstruction for the Jensen construction. A
pseudo-arclength curve may be perfectly continuous in its own parameter `s`
while containing two or more states at the same `h`. The Jensen integral
requires a single-valued state section over `h`, so tracing through a fold
does not by itself say which leg to integrate.

### 4.3 Coordinate conditioning

Some apparent singularity is caused by coordinates rather than the root
manifold:

- the static residual row can be millions of times larger than the self-energy
  rows;
- row equilibration changed a measured raw/equilibrated reciprocal condition
  estimate from about `8.3e-13` to `6.1e-9`;
- a small, fully re-audited `K0` polish repaired a last-bit static-closure
  stall;
- reparameterizing by `w=z-K0` gave well-conditioned bordered solves along a
  retained root-6 leg.

The latest discrete fixed-`w` replay matched 33 fixed-`h` states, corrected
all 32 positive targets in at most four Newton iterations, and reached small
residuals with a positive bordered condition estimate. This is strong
evidence that `w` is a useful numerical coordinate.

It is **not** yet a continuous-edge certificate. Sampled accepted points do
not prove that a connected solution exists throughout every intervening
`w` interval, that the bordered Jacobian stays nonsingular there, or that an
event was not crossed and re-entered. The current work therefore retains
fixed-`w` as a conditioning/solver result and fails closed on graph-edge
promotion.

## 5. Why the obvious simple fixes are inadequate

| Proposed shortcut | What it can help | Why it is not a final fix |
|---|---|---|
| More Picard iterations | A weakly contractive case | Does not remove poles, folds, or multiple basins |
| Smaller mixing | Can turn jumps into slow drift | No uniform contraction guarantee; can take hundreds or thousands of iterations and still select the wrong root |
| Looser tolerances | Makes more nodes print “converged” | Converts numerical uncertainty into hidden physics error |
| Newton at every failed node | Repairs locally regular nodes | Can converge to any coexisting root and jump branch |
| Pseudo-arclength | Crosses fixed-`h` folds | Produces a curve, not necessarily a single-valued Jensen section |
| Interpolate failed nodes | Makes quadrature numerically complete | Invents the integrand across the region where branch identity is least known |
| Use the bare ordered endpoint | Deep ordered endpoints are easy to solve | Changes the thermodynamics and order-parameter onset |
| Remove the resummed denominator | Removes the numerical poles | Changes the approximation; the tested strict candidate failed its registered validity gate |
| Accept any smooth-looking spectrum | Useful exploratory visualization | A visually smooth mode does not validate the state used to compute it |

The central issue is that **root finding and state selection are different
tasks**. Numerical conditioning can answer “can this root be followed?” It
cannot alone answer “is this the state the truncated theory prescribes?”

## 6. What would constitute a mathematically clean solution

There are two principled routes.

### Route A — a common thermodynamic functional

If the coupled resummed equations can be derived as stationarity conditions
of one off-shell functional, competing roots can be compared through the same
thermodynamic object, stability can be graded consistently, and Maxwell/
field-derivative identities provide independent checks.

The present derivation work has not yet established such a functional for
the complete retained approximation. An exact local bilocal/2PI-kernel route
remains open, so this route is not declared exhausted.

### Route B — an explicit branch prescription

The documented backup in `biased_convergence_solution.md` reconstructs all
admissible candidates and selects the root path with the smoothest
dimensionless `r(h)` history, subject first to:

- complete residual and A--D acceptance;
- no unregistered pole/mean event;
- forward/reverse/cold-seed agreement;
- a unique physical Jensen endpoint;
- continuity toward the QCP.

The selection principle is consistent with a continuous phase transition:
the state should not jump discontinuously across the QCP. It is still a
prescription that must be declared and tested, not something a Newton solver
may decide implicitly.

## 7. Practical path to a visual preview

The preview and the final solver should have different, explicit claims.

### Preview claim

An opt-in preview may display a column when:

1. the coupled residual and independent A--D audit pass at every used node;
2. a deterministic guarded continuation follows one recorded branch;
3. predictor distance, condition estimate, pole margin, and mean margin stay
   inside frozen thresholds;
4. any ambiguity, event, fold that cannot be sectioned, or provenance
   mismatch masks the column rather than substituting another root;
5. output metadata says `diagnostic_preview`, including the branch/seed ID.

This is sufficient to reveal whether spectra improve or whether new
artefacts appear. It is not sufficient to make the branch the production
default or claim a complete physical phase diagram.

### Final claim

Default production use additionally needs either:

- the common-functional derivation and its thermodynamic gates; or
- complete candidate enumeration and the preregistered smooth-`r`/QCP
  branch prescription, with grid, Matsubara, direction, and quadrature
  refinement.

For a rigorous continuous fixed-`w` edge, finite samples must be replaced by
an interval Newton/Krawczyk-style cover (or an equivalent validated
continuation theorem) proving existence, local uniqueness, nonsingularity,
and event clearance over every slab. That proof is valuable for the final
claim, but it need not block an honestly labelled visual preview.

### 7.1 Exact audit-level failure mechanism found on 2026-07-28

The QCP-anchored 4.05 T experiment isolated a simpler defect inside the diagnostic acceptance
layer. The simultaneous equations use independent coordinates `[Sigma(:);K0]`, but the old
Block-A audit replayed the production map, including a second nested static Picard solve. In the
multi-root regime that inner solve could select a different `K0` branch and report the difference as
if the supplied simultaneous root had failed its equations.

At seven representative corrected roots between `h=0.0041142866` and `0.010569793 meV`, coupled
Block A equalled the Sigma part of `invz_ordered_node_equations` exactly in reported floating
arithmetic, and defactored Block B equalled the absolute static component exactly. Blocks A/C/D
were at `0` to `7.3e-15` and Block B stayed below `7.72e-10`. The old nested replay on those same
accepted states falsely reported `3.40e-3` to `3.12e-2`.

The optional coupled audit now mirrors the simultaneous equations directly, independently rebuilds
canonical `K`/`lambda`, and contains no nested static solve unless an explicitly report-only debug
flag is enabled. Exact wiring checks reject stale `lambda` or `K(1)~=K0s`. The default nested
production audit remains unchanged.

With that correction, 126 consecutive QCP-anchored arclength states extended the smooth branch from
`h=0.0041142866` to `0.0117280521 meV`, all below `9.76e-10` coupled residual. Equation convergence
and branch selection nevertheless remain distinct: at the identical
`h=0.01171732294388717 meV`, two A--D-accepted roots have `r=0.768127507` and `r=0.822169537`,
separated by only `0.00103702` in the frozen scaled state metric. Natural-`h` Newton selected the
latter; QCP-anchored arclength selected the former smoothly.

### 7.2 QCP-first susceptibility result on 2026-07-28

The research-priority observable does not require the unresolved low-field
component. A direct production run at `T=0.10 K`, `q=[0 0 0]`, `B=4.60--4.90
T`, and `0--6 GHz` completed all 61 field columns:

- 19 Jensen-ordered columns and 42 stable paramagnetic columns;
- no masked, suspect, or non-finite peak columns;
- `Bc_1z = 4.6925 T`, with the finite-sweep bracket
  `[4.690,4.695] T`;
- a continuous V-shaped 1/z mode, falling from about `1.00 GHz` at
  `4.60 T` to `0.174 GHz` at `4.690 T`, then rising from `0.110 GHz`
  at `4.695 T` to about `1.02 GHz` at `4.90 T`.

The nonzero sampled minimum is not evidence for a finite critical gap: the
field mesh does not contain the inferred critical field, the frequency mesh
is `0.01 GHz`, and the real-axis HWHM is `5e-5 meV` (about `0.012 GHz`).
Resolving the asymptotic gap more closely requires a local field/frequency
refinement and a controlled broadening check, not a different ordered
branch.

This result separates two questions that had become entangled. The
near-QCP susceptibility and poles are now available from the accepted
production equations. Robust coverage of the low-field ordered phase
remains an algorithmic branch-selection problem, but it is not a blocker
for comparison with experiments near the QCP.

### 7.3 Finite-grid and edge-pair result on 2026-07-28

The last sentence above must be read as a visual-availability statement,
not a quantitative comparison authorization. A state-only grid ladder found:

| grid | `Bc_1z` (T) | contiguous accepted ordered width (T) |
|---:|---:|---:|
| 12³ | 4.682284546 | 0.382284546 |
| 16³ | 4.692758179 | 0.267758179 |
| 20³ | 4.699093628 | 0.224093628 |
| 24³ | 4.702957153 | 0.177957153 |

The width scales approximately as `N^-1.076`, while the excluded-Γ coupling
gap scales as `N^-1.103`. The QCP-side accepted band is therefore a
finite-grid computability sliver. The `0.02067 T` movement of `Bc` across
the ladder also falsifies the prior expectation that it would stay within
`0.01 T`.

A phase-aligned response gate gives a more favorable and complementary
result. At ten matched offsets `B-Bc(N)` between `-0.080` and `+0.080 T`,
all 40 susceptibility columns are finite and correctly phased. The
`12^3--24^3` soft-mode peak curves differ by only `0.38--0.53%`, and a
halved frequency step reproduces the extreme-grid spreads. The mode shape
is therefore robust after alignment; the dominant unresolved grid effect
is where the curve sits on the absolute field axis.

The representative `16^3` rejected/accepted pair at 4.400/4.425 T does not
change static pole interval: every iterate remains on the rightmost
`y>Jmax` root. The rejected predictor already has a static closure residual
near `5e-11`; its outer Σ residual contracts only by about `0.916` per
iteration and misses the 200-step budget. Raising the cap to 1000 accepts
4.400 T but still fails 4.300 T. This is direct evidence that the mask edge
is controlled by outer-map conditioning near the finite-grid pole, not by
critical slowing or a universally wrong inner root.

The exact area identity
`h0-J0*m=integral(crit dh)` improves the conditioning audit but does not
remove this obstacle. Doubling the HMF grid approximately quarters the
two-route quadrature discrepancy, yet the finer grid samples a new failed
node. The reproducible packet is
`docs/diagnostics/invzp_qcp_grid_2026-07-28/`.

## 8. Current work estimate

For a research-priority QCP visual inspection, no further solver work is
required: the verified finite-`16^3`, `4.60--4.90 T`, `0--6 GHz` subset
displays both the susceptibility map and extracted peak-energy curve. The
interactive driver currently spans 3--6 T, so masked lower-field columns are
expected.

For a quantitative QCP comparison, work remains: extend the coupling/state
grid ladder until `S_N(J0)`, solver-grade `Bc`, and accepted-state support
have a defensible continuum limit. The phase-aligned susceptibility peak
curve is already a sub-percent four-grid regression; spectral weight and an
independent analytic pole solve remain ungraded.

The following results still govern any later attempt to make the full
`0--9 T` ordered map complete:

1. a coupled Newton correction repairs the small `3.6 T` node deficit but
   does not complete the `1.5 T` path;
2. the defactored static coordinate can repair an isolated static solve but
   does not complete the production outer map, even with 1000 iterations;
3. matched field-to-field profile seeding can reduce, rather than improve,
   the accepted-node count.

The retained fixed-`w` samples still do not certify continuous graph edges.
Accordingly, low-field production coverage still requires either the
explicit Route-B branch prescription with ambiguity masking or a successful
common-functional Route-A selector. That longer task is now secondary and
can be resumed without holding up QCP spectra.

## 9. Bottom line

The user's visual observation is consistent with the measured mathematics:
the nonconvergence is nonuniform and substantially numerical. Near-QCP
success and moderate-field failure do not contradict a continuous physical
transition because the failing iteration is not physical dynamics.

The simple diagnostic improvement is to replace a noncontractive, badly
scaled natural-`h` iteration by a coupled, guarded solve in a better
coordinate. The non-simple parts are preventing that stronger solver from
quietly choosing among several valid roots and establishing a coupling-grid
limit. A fixed-G rightmost-root bisection is not a drop-in cure: the
implemented inner equation recomputes `Gstat(K0)`, and the measured edge
failure occurs in the outer Σ block after the inner block is already closed.
