# Ordered \(1/z\) root-search problem

## Purpose

This document is a self-contained numerical problem statement to be read
together with `jensen_1z_framework.html`, especially Sections 9.2--9.4. It is
intended to support requests for mathematical, physical, and numerical advice.
It distinguishes verified observations from hypotheses and does not assume
that all nonconvergence has one cause.

The ultimate project objective is a certified projected-spin \(1/z\)
susceptibility calculation over transverse fields

\[
0 < B_x \leq 9\ {\rm T},
\]

possibly treating strict \(B_x=0\) separately. The diagnostics summarized here
use \(T=0.1\) K unless stated otherwise.

## Executive summary

Equation (45) is the global ordered-state equation, but the present code
normally fails before its final root search. Its integrand
\(r(h)=G_0(h)/\widetilde G_0(h)\) requires a physical coupled
effective-medium solution at every intermediate \(h\). At present:

- the separately evaluated \(h=0\) endpoint has no admissible static root
  along the production outer trajectory, so the integral cannot be seeded;
- failed low-\(h\) nodes can commit a partially updated self-energy and make a
  later node depend on grid history;
- coarse continuation can miss a nearby solution component; and
- one certified high-\(h\) component at 4 T terminates at a physical static
  pole at finite positive \(h\), rather than extending to zero.

Thus the evidence does not presently indict the derivation of equation (45).
It shows that the approximate self-consistent closure used to construct its
integrand does not yet provide a single-valued, admissible, selected branch on
the required interval.

## Coordinate and notation convention

The framework document writes the ordered moment and its conjugate molecular
field along \(x\). In the production calculation the applied transverse field
is called \(B_x\), while the ordered moment is along \(z\). These are related by
an axis relabelling; the mathematical structure of Section 9.3 is unchanged.

The implementation absorbs \(g\mu_B\) into its longitudinal molecular-field
energy and denotes it by \(h\) or `hmf`. Thus:

- \(B_x\) is a fixed external transverse-field parameter;
- \(h\) is the longitudinal \(H_{\rm MF}\) integration/root variable at that
  fixed \(B_x\);
- \(h_0\) is the energy-unit version of \(H_0\);
- \(J_{0,\rm eff}\) is the uniform/order-mode coupling used in the spontaneous
  longitudinal self-consistency equation.

## Governing outer equation

Section 9.3 derives

\[
\langle J_z\rangle(h)
 =-g\mu_B\int_0^h G_0(0;h')\,dh',
\tag{43, relabelled}
\]

and

\[
H_0(H_{\rm MF})
 =\int_0^{H_{\rm MF}}
   \frac{G_0(0;H')}{\widetilde G_0(0;H')}\,dH'.
\tag{45}
\]

For zero applied longitudinal field,

\[
h_0=J_{0,\rm eff}m(h),
\]

where \(m(h)=\langle J_z\rangle\) is obtained directly from the single-ion
calculation. The production root problem is therefore

\[
\boxed{
F(h)=
\int_0^h r(h')\,dh'
-J_{0,\rm eff}m(h)=0,
\qquad
r(h)=\frac{G_0(0;h)}{\widetilde G_0(0;h)} .
}
\]

A nonzero root gives the spontaneous ordered molecular field. Equation (45)
is consequently the global equation in which the problem appears, but its
integrand is the output of another coupled nonlinear problem at every \(h\).

## Nested numerical problem at fixed \(h\)

For each candidate \(h\), the code must:

1. diagonalize the full single-ion electronuclear Hamiltonian at fixed
   \((B_x,h)\), producing \(m(h)\), the bare response \(G_0(i\omega_n;h)\),
   and its elastic/inelastic static pieces;
2. solve the dynamic effective-medium equations for \(K_{n>0}\) from the
   current self-energy \(\Sigma_n\);
3. solve the ordered static effective-medium closure for \(K_0\) and
   \(G_{\rm stat}\);
4. compute the fluctuation sums \(\lambda_p\);
5. update the ordered self-energy \(\Sigma_n\);
6. repeat steps 2--5 until the complete outer residual closes.

The corrected deterministic diagnostic regards this as a map

\[
\Sigma_{\rm new}=\mathcal F_h[\Sigma].
\]

Production currently evaluates the same ingredients through a damped Picard
loop rather than a simultaneous nonlinear root solve.

### Three distinct roots

The word “root” refers to three nested problems that must not be conflated:

| Level | Unknown | Equation | Meaning of failure |
|---|---|---|---|
| Static closure | \(x=\widetilde G_0(0)\) | \(\widehat R(x;\Sigma,h)=0\) | At the current \((\Sigma,h)\), no root passes the physical static gates |
| Coupled node | \(\Sigma\) and dependent \(K,\lambda,x\) | \(\mathcal F_h[\Sigma]-\Sigma=0\) | The chosen nonlinear method did not certify a complete fixed point; another basin/root may still exist |
| Section 9.3 outer root | \(h=H_{\rm MF}\) | \(F(h)=\int_0^h r(h')dh'-J_{0,\rm eff}m(h)=0\) | No ordered molecular field can be certified from the selected integrand profile |

The production status `no_admissible_static_root` belongs to the first level
at one iterate of the second level. It is not, by itself, a nonexistence proof
for either the complete coupled-node root or the final Section 9.3 root.

### Bounded static closure

The static problem is written in the physical coordinate

\[
x\equiv\widetilde G_0(0)
 =\frac{G_{\rm stat}}{1-K_0G_{\rm stat}},
\qquad
-\frac{1}{J_{\rm sup}}<x<0,
\]

where

\[
J_{\rm sup}
=\max\!\left[J_{0,\rm eff},\sup_{q,\nu}J_{q\nu}\right].
\]

For the frozen production-equivalent \(16^3\) coupling mesh,

\[
J_{\rm sup}=J_{0,\rm eff}=0.00642166180942\ {\rm meV},
\qquad
\max_{q,\nu}J_{q\nu}=0.00637172573616\ {\rm meV}.
\]

The separately weighted uniform/directional mode, rather than the largest
sampled mesh value, therefore sets the physical endpoint.

For a given \(x\),

\[
\Phi(x)=\left\langle\frac{x}{1+J_qx}\right\rangle_q,
\qquad
K_0(x)=\frac{1}{\Phi(x)}-\frac{1}{x},
\]

and the scalar equation is

\[
\widehat R(x)
=\Phi(x)
-G_{\rm stat}\!\left(
K_0(x),\Sigma_0,\lambda_1,\lambda_2,\ldots
\right)=0.
\]

The entire configured interval is searched. A candidate is accepted only if,
in addition to satisfying the scalar residual, it has:

- positive supremum mass \(1+J_{\rm sup}x\);
- positive uniform mass
  \(1+(J_{0,\rm eff}-K_0)G_{\rm stat}\);
- positive lattice denominators in both the \(x\) and medium coordinates;
- physical and finite elastic-resummation variable \(\xi\);
- finite static quantities and a negative \(G_{\rm stat}\);
- root and equivalent closure residuals within \(10^{-10}\).

Zero admissible roots makes \(\mathcal F_h[\Sigma]\) undefined at that
\(\Sigma\). More than one admissible root makes it multivalued unless an
additional branch selector is supplied. The static solver deliberately ignores
its legacy `K0_seed` and enumerates the bounded roots, so its mathematical root
set is seed-independent at fixed input data.

## Current production algorithm

The bounded static root enumerator is wired into the production ordered
profile. The deterministic outer-map evaluator and adaptive continuation used
in the later diagnostics are **not** wired into production.

At each fixed transverse field \(B_x\), production:

1. evaluates a separate \(h=0\) node to obtain
   \(r(0)=G_0(0;0)/\widetilde G_0(0;0)\) and the linear ordering predictor;
2. constructs a geometric grid, normally 33 nodes, from a small positive
   \(h_{\min}\) to a bare-mean-field-derived \(h_{\max}\);
3. evaluates nodes from low \(h\) to high \(h\), passing the returned
   \(\Sigma\) and static carrier to the next node as warm starts;
4. approximates equation (45) by a cumulative trapezoidal rule seeded with
   \(r(0)\);
5. forms \(F(h)=h_0(h)-J_{0,\rm eff}m(h)\);
6. brackets the last negative-to-positive nonzero root and refines it with
   direct node evaluations.

The current strict contract rejects the whole ordered point if:

- the \(h=0\) predictor node fails;
- any profile node fails;
- a bisection or final-root node fails; or
- a predicted nonzero root cannot be bracketed/refined.

This strict rejection is intentional: finite values returned by a failed
static or outer solve are not treated as physical susceptibility data.

## Observed failure classes

### 1. The \(h=0\) endpoint makes equation (45) undefined in production

In the 1, 2, and 3 T resolution census, every \(h=0\) predictor returned

```text
no_admissible_static_root
```

and `predictor_converged=false`. Production therefore has no valid \(r(0)\).
It sets the cumulative equation-(45) integral to `NaN` and ultimately masks
the ordered point.

This explains why the present production run masks the whole ordered regime
while the paramagnetic regime remains normal: the PM path does not require
this ordered \(H_{\rm MF}\) profile and integral.

Changing the number of positive-\(h\) nodes cannot affect the independently
evaluated \(h=0\) node.

This status does **not** prove that the full coupled
\((\Sigma,K,\lambda,x)\) equations possess no root at \(h=0\). It proves only
that the production outer trajectory reached a \(\Sigma\) for which the full
bounded static scan found no admissible \(x\). A disconnected coupled solution
or a different outer basin has not been excluded.

### 2. Increasing equation-(45) node resolution is not a solver repair

Only `nH` was changed in a production-path census at 1, 2, and 3 T. All other
options were frozen:

- `nH = 33, 65, 129`;
- `mix_outer = 0.35`;
- `max_outer = 200`;
- `tol_outer = 1e-8`;
- `Ecut = 40` meV;
- static `scan_points = 4097`;
- static `endpoint_margin = 1e-10`;
- static `resid_tol = 1e-10`.

The result was:

| \(B_x\) | converged profile nodes, 33 | 65 | 129 | first converged \(h\), 33/65/129 |
|---:|---:|---:|---:|---:|
| 1 T | 10/33 | 17/65 | 31/129 | 0.00624162 / 0.00774547 / 0.00862825 |
| 2 T | 8/33 | 15/65 | 30/129 | 0.00891140 / 0.00891140 / 0.00844323 |
| 3 T | 6/33 | 11/65 | 22/129 | 0.0115294 / 0.0115294 / 0.0109237 |

All nine final points remained masked. At every resolution and transverse
field, the failed nodes form one contiguous low-\(h\) block and the successful
nodes form the complementary high-\(h\) block. All profile failures are
`no_admissible_static_root`; no multiple-root, unresolved-search, or pure
outer-iteration failure was reported.

The grids are exactly nested: every 33-node value occurs on the 65-node grid
and every 65-node value occurs on the 129-node grid, to about \(10^{-17}\) in
the internal field coordinate.

At 2 and 3 T, every shared node has exactly the same verdict. The slightly
lower first-success value at 129 nodes is one newly inserted geometric
half-step, a factor 0.9474635 below the coarser boundary. This is finer
localization of the same transition, not restoration of a lost branch.

### 3. Failed nodes contaminate the production warm-start path

At 1 T, the verdict changes at identical shared nodes:

\[
\begin{array}{c|c|c}
h & \text{coarser path} & \text{denser path}\\
\hline
0.00624162311 &
\text{success},\ \Sigma_0=-0.292673 &
\text{failure},\ \Sigma_0=-0.821288\\
0.00774546581 &
\text{success},\ \Sigma_0=-0.176826 &
\text{failure},\ \Sigma_0=-0.861066
\end{array}
\]

These differences cannot be quadrature or node-location errors because the
field values are shared exactly.

The production node evaluator mutates `Sigma` after every successful inner
static iteration. A later iteration at the same node may encounter
`no_admissible_static_root` and reject the node, but the partially updated
`Sigma` is still returned to the profile sweep. The next node therefore starts
from a state produced by a failed predecessor.

The evidence strongly supports failed-node state contamination as the source
of the 1 T shared-node contradictions. A direct transactional rollback test is
still required before calling it uniquely causal.

The consequence for equation (45) is fundamental: production does not yet
evaluate a grid-independent function \(r(h)\). It can instead evaluate

\[
r(h;\text{history of previous failed nodes}).
\]

Adding more quadrature nodes can therefore make the result worse by adding
more failed state transitions.

### 4. Coarse continuation can miss a nearby coupled component

At 4 T, a deterministic high-to-low traversal of the original 33-node grid
certified only the five highest-\(h\) nodes. A direct transfer from node 29 to
node 28 had no admissible static root.

Adaptive downward parameter steps, while preserving the physical acceptance
gates, subsequently reached node 28. This established that:

- absence of a root after one coarse state transfer does not prove absence of
  a nearby coupled branch;
- adaptive continuation can be necessary even when the final fixed point is
  locally contractive;
- reversing the coarse grid is not sufficient.

A separate generic root-polishing defect was found during this work. Sign
brackets formerly stopped at the first point satisfying
\(|\widehat R|\le10^{-10}\), while the algebraically equivalent closure
residual could remain barely above \(10^{-10}\). This made admissibility depend
on scan-grid phase. Sign roots are now polished to the coordinate tolerance,
and the affected false failures disappeared. That numerical defect is fixed
and should not be conflated with the remaining blanket production mask.

### 5. The continued 4 T component reaches a physical-domain endpoint

Continuing downward from the certified 4 T node-28 root did not reach node 27.
Instead, accepted steps accumulated at

\[
h_c\simeq0.0080428632,
\]

where the uniform and supremum masses approached their common zero:

\[
1+J_{\rm sup}x\longrightarrow0,
\qquad
1+(J_{0,\rm eff}-K_0)G_{\rm stat}\longrightarrow0.
\]

The last accepted roots were stable over a 3-by-3 static
scan-density/endpoint-margin ladder. A retained proposal below the edge had
zero mathematical and zero admissible static roots over the same ladder.
Meanwhile:

- the outer-map dominant eigenvalue remained approximately \(-0.20\), hence
  locally contractive;
- nonuniform medium and dynamic denominator margins remained positive.

The stopping mechanism is therefore not simply insufficient Picard damping or
field resolution. The followed high-\(h\) component reaches the boundary of
the physical static domain at finite positive \(h\).

This component cannot supply a continuous integrand for equation (45) from
\(h=0\) to the ordered root. The result does not rule out a disconnected
admissible component below \(h_c\), nor does it establish which component is
thermodynamically stable.

### 6. At least one 1 T node is locally noncontractive

At the previously studied 1 T node 22, admissible undamped and damped Picard
trajectories leave the static domain. The last admissible state's dominant
outer-map eigenvalue is approximately

\[
\lambda_{\rm dom}=1.35957.
\]

No positive scalar damping can make that local mode contractive. This is a
different obstruction from the contractive 4 T component endpoint and is
evidence against a universal mixing-parameter repair.

It remains possible that another admissible coupled root exists at the 1 T
node; the Picard result is not an all-roots nonexistence proof.

## Central interpretation

Equation (45) assumes that, along the integration path:

1. a self-consistent physical \(\widetilde G_0(0;h)\) exists at every \(h\);
2. a unique physical branch, or an independently justified branch selector,
   defines the integrand;
3. the selected branch can be followed consistently from \(h=0\) to the
   desired ordered root;
4. numerical evaluation of a node is independent of the sequence of failed
   trial nodes used to reach it.

The current implementation satisfies none of these prerequisites globally.
There is no evidence that the algebraic derivation of equation (45) is itself
incorrect. The problem is that its current ordered effective-medium integrand
is not yet a certified single-valued function over the required interval.

At least three separable causes are present:

1. **Algorithmic:** failed nodes commit partially updated continuation state.
2. **Nonlinear/root-search:** coarse transfers and locally noncontractive
   branches defeat damped Picard even when another nearby solution may exist.
3. **Mathematical/physical:** an admissible high-\(h\) component can terminate
   at the static susceptibility pole before reaching \(h=0\).

There may also be a representation issue: production combines a full
136-state electronuclear \(G_0\) with a two-level electronic ordered
self-energy/vertex. This mismatch is verified structurally but has not yet
been shown to cause the component endpoint or the \(h=0\) failure.

## What is established and what remains open

### Established

- The blanket ordered mask is presently triggered by the failed \(h=0\)
  integrand endpoint and strict rejection logic.
- Increasing `nH` does not repair that endpoint or the final point.
- The 1 T production sweep is path-dependent because failed-node state is
  retained.
- Some coarse transfer failures can be crossed with smaller continuation
  steps.
- The followed 4 T high-\(h\) component nevertheless terminates at a
  resolution-stable physical-domain pole.
- A scalar damping change cannot repair every observed class.
- The corrected static root enumerator no longer has the earlier
  residual-polishing/grid-phase defect.

### Not established

- That the full coupled equations have no admissible solution at \(h=0\).
- That no disconnected low-\(h\) component exists.
- That the high-\(h\) continued component is the equilibrium component.
- That components may be joined across the pole.
- That equation (45) permits a discontinuous branch switch or piecewise
  integration without an additional thermodynamic construction.
- That the present electronic/electronuclear hybrid is a consistent
  realization of Jensen's finite-moment ordered theory.
- That equation (45), rather than the approximate closure used to evaluate its
  integrand, must be modified.
- That any susceptibility computed from failed or pseudo-root nodes is valid.

## Advice sought

The most useful external suggestions would address the following bounded
questions.

### Mathematical/physical questions

1. Does Jensen's derivation of equation (45) require one differentiable
   equilibrium component to connect \(h=0\) to the finite ordered solution?
2. If the self-consistent static component terminates at
   \(1+J_{\rm sup}x=0\), should this be interpreted as a spinodal/collective
   instability, a failure of the ordered effective-medium approximation, or a
   signal to change the reference phase?
3. If disconnected admissible components exist, may equation (45) be evaluated
   piecewise with a branch switch? If so, what thermodynamic matching or
   integration constant is required?
4. Should the free-energy relation in Section 9.4, equation (46), be used as
   the branch selector before attempting equation (45)?
5. Is it consistent for the ordered branch used in equation (45) to approach
   \(h=0\) through the paramagnetic solution, or must the ordered elastic
   closure remain valid all the way to the endpoint?
6. Could the full-electronuclear-\(G_0\)/two-level-self-energy mismatch
   invalidate the finite-\(h\) component topology even if the equations are
   numerically solved exactly?

### Numerical questions

1. What is the most appropriate simultaneous formulation of the coupled
   residual in \((\Sigma,K_0,x,\lambda)\), avoiding a nested Picard iteration
   that becomes undefined when the inner root set changes?
2. Would pseudo-arclength continuation in \((h,\Sigma,x)\) correctly classify
   folds separately from the hard pole \(1+J_{\rm sup}x=0\)?
3. How should disconnected coupled roots be enumerated in a high-dimensional
   Matsubara self-energy problem: deflation, homotopy, interval methods,
   reduced coordinates, or another method?
4. How should a continuation algorithm behave when a trial node leaves the
   physical static domain? The immediate proposal is transactional rollback
   to the last certified state plus a separate fresh-start attempt.
5. What numerical evidence would be sufficient to distinguish genuine
   nonexistence of a coupled root from failure of the selected outer solver?
6. Once a physical branch is selected, what error-controlled quadrature is
   appropriate for equation (45), particularly near narrow electronuclear
   response features or a component endpoint?

## Proposed next bounded diagnostic

No production integration is proposed yet. The least expensive discriminating
test is:

1. make the diagnostic node transition transactional;
2. on rejection, restore the last certified \(\Sigma\) and static carrier;
3. separately evaluate the same target from a documented fresh start;
4. retest only the two contradictory shared 1 T nodes;
5. verify whether their result becomes independent of the 33/65/129 path.

This can confirm or reject the state-contamination hypothesis. It cannot prove
that a physical branch reaches \(h=0\), cross the certified 4 T pole, or select
an equilibrium component.

If transactional evaluation removes the shared-node contradiction, the next
mathematical task is a bounded search for disconnected coupled components
below the 4 T endpoint, followed by thermodynamic selection. Equation (45)
should be integrated only after those steps establish a valid branch.

## Implementation and evidence index

- Framework derivation:
  [`jensen_1z_framework.html`](jensen_1z_framework.html), Sections 9.2--9.4.
- Production ordered profile and equation-(45) root:
  [`invz_projected/invz_solve_point_ordered.m`](invz_projected/invz_solve_point_ordered.m),
  nested function `invz_hmf_ordered`.
- Bounded physical static closure:
  [`invz_projected/invz_emt_static_ordered.m`](invz_projected/invz_emt_static_ordered.m).
- Deterministic coupled outer map:
  [`invz_projected/invz_ordered_outer_map.m`](invz_projected/invz_ordered_outer_map.m).
- Execution journal and exact checkpoint claims:
  [`docs/execution/invzp_convergence_journal.md`](docs/execution/invzp_convergence_journal.md),
  especially Checkpoints 3--6.
- Resolution-census artifact:
  `docs/diagnostics/invzp_outer_wp2/wp2_hmf_node_resolution_census.mat`.
  Its checkpointed one-off generator is recoverable from commit `02bc7d3`.
- 4 T adaptive continuation and endpoint artifacts:
  `docs/diagnostics/invzp_outer_wp2/wp2_4t_adaptive_boundary_continuation.mat`,
  `wp2_4t_adaptive_target_audit.mat`,
  `wp2_4t_adaptive_node28_to27.mat`, and
  `wp2_4t_adaptive_component_edge_audit.mat`.
- Rejected explanations and conditions for reopening them:
  [`invzp_convergence_dead_ends.md`](invzp_convergence_dead_ends.md).

The machine-readable artifacts contain the full numerical states. This
document intentionally retains only the decisive values needed to reproduce
or challenge the conclusions.
