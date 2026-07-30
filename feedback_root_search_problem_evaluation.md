# Evaluation of `feedback_root_Search_problem.md`

## Scope and verdict

This note evaluates the feedback against:

- Sections 9.3--9.4 of `jensen_1z_framework.html`;
- the current ordered static and outer implementations;
- the checkpointed 1--4 T diagnostics; and
- the exact elastic-sector identities already established in
  `invzp_convergence_next_steps.md`.

The feedback contains several high-value ideas, but its central conclusion is
stronger than the available evidence. The followed 4 T high-\(h\) component
does terminate at a uniform-mode stability boundary while the reduced outer
map remains contractive. That is good evidence for a physical-domain endpoint
of that component, not an ordinary Picard failure. It does **not** prove that
every failed low-\(h\) node at every transverse field is one connected
spinodal branch, that no disconnected admissible component exists, or that
there is exactly one nonzero root of equation (45).

The most useful recommendations are:

1. derive and test a saturation-anchored form using the exact integrand
   \(r-1\), not generally \(\Sigma_0\);
2. make node evaluation pure, or equivalently commit continuation state only
   after acceptance;
3. derive a reduced simultaneous residual before claiming that the Matsubara
   problem is exactly four-dimensional; and
4. use bordered pseudo-arclength continuation and a full residual-Jacobian
   diagnostic to distinguish folds from physical-domain endpoints.

The proposed removal or reinterpretation of the uniform-mode gate is not
supported by the implementation or the physics represented by that gate.

## Claim classification

| Feedback claim | Evaluation |
|---|---|
| \(F'(h)=r(h)[1+J_{0,\rm eff}\widetilde G_0(h)]\) | **Correct** on a differentiable certified component. |
| A zero uniform mass is a stationary point of \(F\) | **Correct** on that component, subject to the same differentiability assumptions. |
| The certified 4 T component ends at a stability boundary rather than loss of Picard contractivity | **Strongly supported**, but the smallest singular value of the full simultaneous residual Jacobian has not yet been measured. |
| Every low-\(h\) failure is the same unstable middle branch | **Not established**. The 1 T result is continuation-history dependent, and disconnected roots have not been excluded. |
| Exactly one nonzero equation-(45) root is possible | **Not established** without connectedness, monotonicity, and single-crossing assumptions that have not been proved. |
| The old class-A separator is algebraically the exact uniform mass because \(\xi\) is small | **Incorrect as a general statement** for the implemented ordered elastic closure. |
| Re-anchor equation (45) at saturation | **Promising but conditional**. The exact tail integrand is \(r-1\), and the boundary condition must be derived and numerically validated. |
| Use \(\Sigma_0\) in the saturation-tail integral | **Only conditionally correct** where \(r-1=\Sigma_0\); it is not the general ordered-state identity. |
| The \(16^3\) lattice average includes the exact-\(\Gamma\) channel as a discrete weight | **False for the current code**. No such weight is added to \(\Phi\). |
| The hard endpoint is caused by turning a zero-measure \(\Gamma\) cusp into a pole in \(\Phi\) | **Contradicted by the implementation premise**. The lattice average and the separately tested uniform susceptibility are distinct. |
| The complete coupled problem has exactly four scalar unknowns | **Plausible reduction to derive**, not yet an established equivalence. The dynamic per-frequency implicit equations must first be eliminated uniquely. |
| A sigmoid coordinate removes spurious domain exits | **Potentially useful for simultaneous Newton**, but it cannot create a missing physical root and does not explain the current bounded enumerator's `no_admissible_static_root`. |
| Pure node evaluation | **Correct and high priority** for eliminating verified failed-state contamination. |
| Bordered pseudo-arclength continuation | **Correct next-level diagnostic** after a well-defined simultaneous residual exists. |
| Adaptive integration with \(h=h_c+t^2\) | **Plausible after branch selection**; it cannot bridge an interval on which no selected admissible component has been constructed. |
| The hybrid representation must shift \(h_c\) appreciably but cannot remove the spinodal | **Speculative**. The mismatch is real, but neither quantitative conclusion has been demonstrated. |

## 1. Exact derivative identity and what it proves

Write

\[
F(h)=H_0(h)-J_{0,\rm eff}m(h),\qquad
H_0'(h)=r(h)=\frac{G_0(h)}{\widetilde G_0(h)}.
\]

With the project's absorbed \(g\mu_B\) convention, equation (43) gives

\[
m'(h)=-G_0(h).
\]

Therefore, wherever a differentiable certified solution component exists,

\[
\begin{aligned}
F'(h)
 &=r(h)+J_{0,\rm eff}G_0(h)\\
 &=r(h)\left[1+J_{0,\rm eff}\widetilde G_0(h)\right].
\end{aligned}
\]

For an accepted static state, \(G_0<0\) and
\(\widetilde G_0=x<0\), so \(r=G_0/x>0\). Hence the sign of \(F'\) is the sign
of the uniform mass

\[
D_x=1+J_{0,\rm eff}x.
\]

The code records the equivalent medium-coordinate mass

\[
D_{\rm uni}=1+(J_{0,\rm eff}-K_0)G_{\rm stat}.
\]

Since

\[
x=\frac{G_{\rm stat}}{1-K_0G_{\rm stat}},
\qquad
1+J_{0,\rm eff}x
=\frac{D_{\rm uni}}{1-K_0G_{\rm stat}},
\]

and the accepted coordinate transformation has
\(1-K_0G_{\rm stat}=G_{\rm stat}/x>0\), the two masses have the same sign.

This proves that a uniform-mass zero is also a stationary point of \(F\) on a
single differentiable component. It does not supply values of \(F'\) where no
certified component has been found, prove that the high-\(h\) component
continues below its endpoint, or rule out a disconnected component.

### Consequence for the 4 T result

At 4 T, adaptive continuation followed one admissible high-\(h\) component to
\(h_c\simeq0.0080428632\), where both the physical endpoint and uniform-mode
mass approach zero while the reduced outer fixed-point map remains
contractive. This is strong evidence that the followed component reaches a
physical-domain endpoint rather than failing because Picard's local spectral
radius reaches one.

The feedback's fold/endpoint distinction is therefore useful, but it should
not yet be phrased as a full Jacobian proof. The existing diagnostic measured
the outer-map spectrum after eliminating the static solve; it did not measure
the smallest singular value of a simultaneous residual in
\((h,\Sigma,K,\lambda,x)\).

## 2. Why the global spinodal conclusion is premature

Four missing proof obligations prevent the stronger conclusion:

1. **Connectedness is unproved.** A bounded static root enumeration at one
   frozen outer state does not enumerate all roots of the coupled outer
   problem. Disconnected coupled components below \(h_c\) have not been
   excluded.
2. **The 1 T topology is not path independent.** Shared nodes on the
   33/65/129 grids changed classification because a rejected node's partially
   updated self-energy was carried forward. A contaminated contiguous failure
   block is not by itself a topological proof.
3. **The low-\(h\) derivative sign is not observed.** The algebraic identity
   for \(F'\) applies on an existing certified branch. It cannot be used to
   assign \(F'<0\) throughout an interval where no branch values have been
   computed.
4. **Single crossing is unproved.** “Exactly one nonzero root” additionally
   requires a selected connected component and sufficient monotonicity or a
   proof that its uniform mass changes sign only once.

The correct current statement is narrower: the certified 4 T high-\(h\)
component ends at its uniform-mode stability boundary. It is reasonable to
call that endpoint spinodal-like. Applying the same interpretation to every
field is a hypothesis for subsequent continuation and branch-search work.

### The old class-A separator is not exact

The feedback replaces the production static quantity by

\[
\widetilde G_0\simeq\frac{G_0^{\rm inel}}{1+\Sigma_0}
\]

and argues that small \(\xi\) makes this exact enough. In the implemented
ordered elastic closure,

\[
G_{\rm stat}=U+V,\qquad
U=\frac{G_0^{\rm inel}}
        {1+\Sigma_0+K_0G_0^{\rm inel}},
\qquad
V=\xi G_0^{\rm el},
\]

and the exact identity is

\[
\frac{1}{\widetilde G_0}
=\frac{1+\Sigma_0}{G_0^{\rm inel}}
-\frac{V}{U(U+V)}.
\]

A small multiplier \(\xi\) does not bound \(V\) without also controlling the
elastic bare response and the denominators. The exact correction changed the
accepted 1 T pole ratios by about 3--10%, changed one stored anomalous node by
about 53%, and reclassified stored 1 T node 9. The reported 42--60% values are
dominant two-level shares of the full bare response, not elastic shares.
Therefore the historical separator was useful diagnostically but is not a
zero-exception algebraic stability theorem.

## 3. Correct form of a saturation anchor

Equation (45) supplies

\[
\frac{dH_0}{dh}=r(h).
\]

Define \(\delta(h)=H_0(h)-h\). Then the exact differential relation is

\[
\delta'(h)=r(h)-1.
\]

If a derivation establishes \(\delta(\infty)=0\) on the relevant certified
component, the corresponding anchor is

\[
\boxed{
H_0(h)=h-\int_h^\infty [r(h')-1]\,dh'
}
\]

and the zero-applied-longitudinal-field residual becomes

\[
F_\infty(h)
=h-\int_h^\infty [r(h')-1]\,dh'
-J_{0,\rm eff}m(h).
\]

The feedback instead uses \(\Sigma_0\). That substitution is exact in the
special limit where \(\widetilde G_0=G_0/(1+\Sigma_0)\), including the
implemented \(m=0\) limit, but `invz_gstat_ordered.m` explicitly retains the
ordered elastic correction for finite moment. The general tail test must
therefore measure \(r-1\), not merely \(\Sigma_0\).

This is still a serious proposal because it could determine the integration
constant using only the certified high-field component. Before adopting it:

1. derive, rather than assume, the applicable asymptotic condition on
   \(\delta H\);
2. verify that a single certified component reaches sufficiently large \(h\);
3. measure the decay and convergence of \(\int_h^\infty(r-1)\,dh\);
4. verify the result against the original lower anchor in a regime where both
   paths are admissible; and
5. check whether the re-anchored relation is the same physical \(H_0\), rather
   than a different integration constant selected across a component change.

Section 9.4 anchors the **free-energy difference** at saturation. The local
framework text says fluctuations are quenched and sets
\(\delta\Phi=0\); it does not, by itself, explicitly prove the
\(\delta H(\infty)=0\) condition needed above. That condition is physically
plausible but remains a proof obligation.

The boundary-exact linear anchor near the continuous transition is also
valid for its stated near-critical purpose. It cannot by itself meet the
project objective of a certified finite-moment result throughout
\(0<B_x\leq9\) T.

## 4. Lattice average versus the uniform channel

The feedback asks whether the exact-\(\Gamma\) channel is added to
\(\Phi(x)\) with a discrete weight. The answer is:

\[
\boxed{\text{No discrete exact-\(\Gamma\) weight is added to \(\Phi\).}}
\]

In `invz_emt_static_ordered.m`, the lattice average is formed only from the
provided Gamma-dropped mesh:

\[
\Phi(x)=x\,\operatorname{mean}_{q,\nu}
\frac{1}{1+J_{q\nu}x}.
\]

The separately verified uniform/directional value
\(J_{0,\rm eff}\) is used to define the physical search interval and is tested
through the independent uniform-mode mass. It is not inserted into the
Brillouin-zone mean.

Numerically,

\[
J_{0,\rm eff}=0.006421661809416939,
\qquad
J_{\rm mesh,max}=0.006371725736157789\ {\rm meV}.
\]

The relative gap is about \(0.778\%\) when normalized by \(J_{0,\rm eff}\),
not \(0.08\%\). At \(x=-1/J_{0,\rm eff}\), the mesh denominators remain finite.
Thus a three-dimensional band-edge cusp in the thermodynamic lattice average
and a zero of the separately selected uniform susceptibility denominator are
different mathematical objects.

A density-of-states or analytic band-edge quadrature could improve the
accuracy of \(\Phi\), and is appropriate under work package 4. It cannot by
itself remove the uniform susceptibility instability or justify deleting its
gate. The load-bearing implementation premise of recommendation (c) is
therefore false.

## 5. Evaluation of the algorithmic proposals

### Reduced four-scalar residual

Equations (39)--(42) make the updated self-energy depend on a small collection
of fluctuation sums and static variables, so a substantial reduction is
plausible. It is not yet demonstrated that
\((\lambda_1,\lambda_2,\lambda_3,x)\) is an exactly closed four-dimensional
system. For every nonzero Matsubara frequency, \(K_n\), \(\Sigma_n\), and
\(g_n\) still obey an implicit scalar relation. A valid reduction must:

1. eliminate each frequency equation on a unique physical branch;
2. retain all dynamic-denominator gates;
3. prove residual equivalence with the full Matsubara formulation; and
4. reproduce known healthy 4 T states before multistart or deflation is used.

Globalized Newton would remove dependence on Picard contractivity as a local
convergence requirement. It would not, by itself, prove root completeness or
select the thermodynamic branch.

### Sigmoid coordinate

Writing \(x=-\sigma(u)/J_{\rm sup}\) can enforce the open interval during a
simultaneous unconstrained Newton solve. The current bounded static enumerator
already evaluates only inside that physical interval, so its
`no_admissible_static_root` status is not caused by a trial iterate stepping
outside it. Reparameterization cannot create a root that does not exist in the
interval; it can also become poorly scaled as the endpoint is approached.

### Pure node evaluation

This recommendation directly addresses a verified defect. The node evaluator
should return a candidate result without mutating continuation state, and the
caller should commit `Sigma` and the static carrier only after all acceptance
gates pass. Transactional rollback is the same required state-transition
semantics implemented at the caller boundary; a pure interface makes the
contract harder to violate. This is the safest bounded algorithmic change.

### Pseudo-arclength and bordered Newton

This is an appropriate method for following a coupled component through an
ordinary fold and for comparing:

- the smallest singular value of the unbordered simultaneous residual
  Jacobian;
- the uniform/supremum mass;
- the lattice denominator margin; and
- the bordered-system conditioning.

It should follow derivation and validation of a simultaneous residual. A hard
physical boundary can itself make a Jacobian poorly conditioned, so the
classification must use both residual geometry and physical margins rather
than one scalar signature.

### Adaptive integration

Adaptive quadrature or embedded ODE integration is useful only after a
physical branch and integration constant have been selected. If the endpoint
has a verified square-root approach, \(h=h_c+t^2\) is a sensible local
coordinate. It does not solve the current branch-existence/matching problem.

## 6. Representation mismatch

The use of a full 136-state electronuclear bare response with a two-level
ordered self-energy/vertex is a verified structural mismatch and should be
audited before assigning physical meaning to a calculated endpoint. The
42--60% dominant-sector shares show that the omitted sector is not obviously
negligible.

They do not quantify the shift in \(h_c\), prove that the block boundary moves
appreciably, or prove that the same spinodal topology survives a consistent
multilevel closure. Those are hypotheses to test, not conclusions.

## Recommended ordering of follow-up work

1. **Purity/transaction test:** eliminate rejected-state carryover and retest
   only the contradictory shared 1 T nodes. This closes a known algorithmic
   ambiguity without making a physical assumption.
2. **Saturation-anchor derivation and cheap tail feasibility test:** derive
   the asymptotic condition, then measure \(r-1\) on a certified high-\(h\)
   component. Do not substitute \(\Sigma_0\) without checking the elastic
   correction.
3. **Reduced-residual derivation on healthy fixtures:** determine whether the
   per-frequency equations can be uniquely eliminated and verify exact
   equivalence at known solutions.
4. **Bordered continuation:** use the validated simultaneous residual to
   classify folds, domain endpoints, and possible disconnected components.
5. **Representation and lattice convergence audits:** quantify, rather than
   assume, how the hybrid local representation and continuum band-edge
   treatment move the component boundary.

Do not, on the present evidence:

- delete or relax the uniform-mode stability condition;
- declare all low-\(h\) failures one connected spinodal branch;
- infer uniqueness of the nonzero equation-(45) root;
- use unconverged PM endpoint iterates as certified anchors; or
- replace \(r-1\) by \(\Sigma_0\) in a finite-moment ordered tail without an
  explicit equivalence check.
