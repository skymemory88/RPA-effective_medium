# Ordered \(1/z\) convergence diagnosis

## Scope

This is the current technical diagnosis of the projected-spin ordered-state
failure. Read it with `jensen_1z_framework.html`, especially Sections
9.2--9.4. The objective is a certified susceptibility over

\[
0<B_x\leq 9\ {\rm T},
\]

with strict zero field treated separately if necessary. Numerical values below
are for \(T=0.1\) K unless stated otherwise.

The present blanket mask is a solver/path failure upstream of the real-axis
susceptibility. It is not evidence that the physical susceptibility vanishes,
nor is it currently evidence that Jensen's equation (45) is algebraically
wrong.

## Governing problem

After relabelling the ordered axis from \(x\) in the framework to \(z\) in the
code, equation (45) gives

\[
H_0(h)=\int_0^h r(s)\,ds,\qquad
r(h)=\frac{G_0(0;h)}{\widetilde G_0(0;h)} ,
\]

where \(h\equiv H_{\rm MF}\) is stored in energy units. At zero applied
longitudinal field, the ordered root is

\[
\boxed{F(h)=\int_0^h r(s)\,ds-J_{0,\rm eff}m(h)=0.}
\]

The quadrature is therefore the outermost of three different nonlinear
problems:

| Level | Unknown | Required condition |
|---|---|---|
| Static closure | \(x=\widetilde G_0(0)\) | an admissible root of \(\widehat R(x;\Sigma,h)=0\) |
| Coupled node | \(\Sigma\), with dependent \(K,\lambda,x\) | \(\mathcal F_h(\Sigma)-\Sigma=0\) |
| Ordered field | \(h\) | \(F(h)=0\) on one selected node-solution component |

`no_admissible_static_root` is a statement about the first problem at one
trial \(\Sigma\). It is not a proof that the complete coupled node has no
other root, and neither node failure nor node success selects the equilibrium
ordered-field root.

### Physical static domain

The bounded static solver uses

\[
x=\frac{G_{\rm stat}}{1-K_0G_{\rm stat}},\qquad
-\frac{1}{J_{\rm sup}}<x<0,
\]

\[
\Phi(x)=\left\langle\frac{x}{1+J_qx}\right\rangle_q,\qquad
K_0(x)=\frac{1}{\Phi(x)}-\frac{1}{x},
\]

and enumerates the configured interval for roots of

\[
\widehat R(x)=\Phi(x)-G_{\rm stat}
  (K_0(x),\Sigma_0,\lambda_1,\lambda_2,\ldots).
\]

For the production-equivalent \(16^3\) lattice,

\[
J_{\rm sup}=J_{0,\rm eff}=0.00642166180942\ {\rm meV},
\qquad
\max_{q,\nu}J_{q\nu}=0.00637172573616\ {\rm meV}.
\]

An exported root must have positive supremum, uniform, lattice, and dynamic
masses; finite physical \(\xi\); negative finite \(G_{\rm stat}\); and
equivalent closure residuals below \(10^{-10}\). The static solver is
seed-independent at fixed inputs because it enumerates bounded roots rather
than iterating from `K0_seed`.

At a fixed \(h\), production still solves the coupled node with damped Picard
iteration. A candidate continuation state \((\Sigma,K_0)\) is now committed
only after all node gates pass. The strict `full_profile` default uses an
independent \(h=0\) predictor and 33 positive geometric nodes, and rejects the
ordered point if the predictor, any profile node, or final refinement fails.
The two visual integration modes remain implemented but are not selected by
the production driver.

## Current diagnosis

### 1. The immediate strict-mask trigger is the missing lower path

Across the controlled 0.5--2.2 T census, every independent \(h=0\) ordered
predictor returns `no_admissible_static_root`. Each 33-node positive profile
has one transition: 22--26 rejected low-\(h\) nodes followed by 7--11 accepted
high-\(h\) nodes. The same predictor failure occurs in the strict 1, 4, and
4.68 T probes.

Consequently, strict production has neither a valid \(r(0)\) nor a complete
set of certified nodes joining \(h=0\) to the high-\(h\) component. The
equation-(45) integral is undefined under the current acceptance contract and
the ordered point is masked. The PM calculation is unaffected because it
does not use this ordered \(h\)-profile.

Increasing `nH` cannot repair the separately evaluated \(h=0\) node. It can
only localize, or alter transfers near, the boundary between the rejected
low-\(h\) region and the certified high-\(h\) component.

### 2. Rejected-state history was a real defect, but is now removed

Before transactional rollback, failed nodes returned partially updated
\(\Sigma\) and contaminated the next warm start. On nested 1 T grids, the
same \(h\) could therefore pass on a coarse path and fail after extra rejected
nodes on a denser path.

After accepted-state-only commit:

- all shared-node verdicts agree on the 33/65/129 grids;
- \(h=0.0062416231096\) meV is accepted at all resolutions with
  \(\Sigma_0=-0.29267261\) to \(-0.29267263\);
- \(h=0.0077454658051\) meV is accepted at all resolutions with
  \(\Sigma_0=-0.176825909\); and
- the maximum jointly accepted differences are
  \(2.28\times10^{-8}\) in \(\Sigma_0\) and
  \(1.00\times10^{-11}\) meV in \(K_0\).

The strict profiles nevertheless remain `node_failed`: only 10/33, 20/65,
and 39/129 positive nodes certify, and the \(h=0\) predictor still fails.
Rollback restores determinism but does not create the missing component.

### 3. More than one nonlinear obstruction is present

At 4 T, adaptive high-to-low continuation crosses a coarse transfer failure
and reaches the next grid node, proving that one failed transfer is not a
root-nonexistence result. Continuing further reaches

\[
h_c\simeq0.0080428632\ {\rm meV},
\]

where

\[
1+J_{\rm sup}x\rightarrow0,\qquad
1+(J_{0,\rm eff}-K_0)G_{\rm stat}\rightarrow0.
\]

The endpoint is stable over the 2049/4097/8193 static scan ladder and three
endpoint margins. A below-edge frozen proposal has no mathematical or
admissible static root over that ladder. The coupled outer map remains locally
contractive, with dominant eigenvalue about \(-0.20\), and nonuniform/dynamic
margins remain positive. The followed high-\(h\) component therefore meets a
physical static pole at finite \(h\); smaller continuation steps do not join
it to zero.

At a previously examined 1 T node, by contrast, the last admissible outer-map
state has

\[
\lambda_{\rm dom}\simeq1.35957.
\]

No positive scalar damping can make that local Picard mode contractive.
Another coupled root or component has not been excluded. These distinct
contractive-endpoint and noncontractive-iteration cases rule out a universal
mixing or iteration-budget repair.

On any differentiable certified component,

\[
F'(h)=r(h)\,[1+J_{0,\rm eff}\widetilde G_0(h)].
\]

Thus a uniform-mass zero is a stationary point of \(F\) on that component.
This supports a spinodal-like interpretation of the 4 T edge, but it must not
be extended through uncertified intervals or applied to every low-\(h\)
failure.

### 4. The small low-field susceptibility is related, but is not the mask cause

Phase acceptance occurs before `invz_chi_realaxis`; a tiny final response
cannot directly trigger the current phase mask. The measured electronic
two-level matrix element also has the opposite node ordering from the proposed
failure mechanism:

- failed \(h=0\) predictors have \(M^2=27.55\)--30.04;
- failed positive nodes extend up to \(M^2=27.55\)--30.00; and
- certified high-\(h\) nodes include the smallest measured values, from
  \(M^2=0.0220\) at 0.5 T to \(4.235\) at 2.2 T.

The apparently dangerous ordered self-energy factors have the exact
reassociation

\[
A=\lambda_2-\tfrac12[g_0+\beta(1-n_{01}^2)]\lambda_1,\qquad
B(z)=\lambda_1-(1-n_{01}^2)K(z),
\]

\[
\Sigma(z)=-\alpha_m
-\frac{2m^2}{n_{01}^2}B(0)g(z)
+\frac{M^2}{n_{01}^2}\{A+B(z)g(z)\}.
\]

If the other quantities remain finite and \(n_{01}\ne0\), the \(M^2\to0\)
limit is finite. The largest sampled \(m^2/M^2\) is
\(1.37\times10^3\), yet the cancelling product agrees with its stable form to
\(2.23\times10^{-16}\), and the complete self-energy agrees to
\(4.45\times10^{-15}\) while finite. Direct-ratio overflow appears only near
\(M^2=4.89\times10^{-308}\) in the frozen test, far below the measured
minimum.

The low-field near-zero susceptibility and the solver problem can therefore
share a physical origin—strong polarization and redistributed multilevel
spectral weight—without the small number causing numerical breakdown.
Accepted high-\(h\) states can have tiny response; the mask instead arises
because the required anchor/path includes rejected states.

### 5. The representation mismatch is quantitative, not yet the topological cause

Production combines a full 136-state electronuclear bare response with an
electronic two-level ordered self-energy/vertex. This is not internally closed.
The projected two-level static weight divided by the full electronuclear
inelastic weight is 7.03--9.19 at the first accepted 0.5--1 T nodes and reaches
88.2 at one failed node, so it cannot be interpreted as a bounded sector
share.

A field-adapted lowest-16 electronuclear subspace is spectrally smooth: its
16/17 gap is 0.618--1.235 meV, retained thermal mass is one, and the minimum
overlap across a failed/accepted boundary is 0.999996. It is nevertheless not
a controlled vertex reduction. At 0.5 T its first-accepted static share is
98.0%, but its connected-variance share is only 9.86%; at the high endpoint
the shares are 16.1% and 0.521%. On the tested
16/24/32/48/64/96/136 rank ladder, every non-predictor sample needs all 136
states to exceed 90% connected-variance coverage.

An internally closed electronic two-level control also fails every \(h=0\)
predictor and has a shorter certified high-\(h\) component than the hybrid:

| \(B_x\) (T) | closed-model certified nodes, including \(h=0\) | lowest certified \(h\) (meV) |
|---:|---:|---:|
| 0.5 | 9/34 | 0.00785085 |
| 1.0 | 7/34 | 0.01192745 |
| 1.8 | 5/34 | 0.01740279 |
| 2.2 | 5/34 | 0.01659925 |

The hybrid remains an important susceptibility systematic, but neither cheap
representation replacement restores the lower integration path.

### 6. A good quadrature estimator is sufficient only after anchoring and branch selection

The visual experiments show that equation (45) need not be integrated to
machine precision for the final root to look physically sensible. Dense
finite-node filtering produced plausible morphology over much of the ordered
range. This supports using an error-controlled estimator once \(r(h)\) is a
certified selected function.

It does not solve the present problem:

- the two-endpoint rule gives finite coarse roots at 3.5--4.5 T only by
  joining a PM lower endpoint with negative uniform mass to a positive-mass
  ordered component across an unresolved pole;
- the filtered 256-node low-field roots use unconverged PM anchors and bridge
  the missing lower interval;
- at 4.68 T, 31 certified retained nodes still give
  `filtered_no_bracket`; and
- discarding `NaN`/`Inf` nodes removes evidence of the missing path but cannot
  determine the integration constant or justify a component switch.

An upper/saturation anchor also fails for the current hybrid. Writing
\(\delta h=H_0-h\) gives

\[
\delta h'(h)=r(h)-1,\qquad
\delta h(h)=-\int_h^\infty[r(s)-1]\,ds
\]

only if \(\delta h(\infty)=0\), \(r-1\) is integrable, and the certified model
reaches that asymptotic regime. On accepted hybrid nodes from 8 to 128 times
the ordinary profile ceiling, fitted \(|r-1|\) exponents are 0.593 at 1 T and
0.811 at 4 T, both below the integrability threshold. Full electronuclear
connected-variance exponents are about 0.17 and the variance is still
\(O(1)\). The matched closed two-level control gives integrable exponents
1.733 and 1.568. Driving the finite 136-state hybrid to numerical saturation
requires fields where omitted higher multiplets invalidate the retained-space
model. A finite upper cutoff would therefore impose an arbitrary additive
constant.

## Crux of the unresolved problem

Equation (45) requires more than finite endpoint values or a good estimator.
It requires:

1. a coupled node solution at each integration coordinate;
2. a branch/component selection rule;
3. a justified integration constant; and
4. a path on which \(r(h)\) is well defined and sufficiently regular.

The present certified positive-\(h\) solution is separated from \(h=0\) by a
region where nested Picard/static evaluation is undefined or leaves the
positive-mass domain. It is unresolved whether:

- a disconnected admissible coupled component exists below the observed edge;
- a simultaneous solver can find roots missed by nested Picard;
- the constrained states used in deriving equation (45) must themselves obey
  the current positive-mass stability gates; or
- the formal construction requires an analytic continuation or a different
  thermodynamic reference through that interval.

Principal-value poles, unstable roots, and discontinuous component switches
must not be admitted without a derivation that supplies the corresponding
thermodynamic branch and integration constant.

## Next discriminating work

The next step is not another node-count or damping trial.

1. Derive an exactly equivalent reduced simultaneous residual, preferably
   eliminating uniquely solvable nonzero-frequency \(K_n\) variables and
   retaining all static, uniform, lattice, and dynamic gates.
2. Reproduce healthy 4 T fixtures and existing component points before using
   globalized Newton, multistart, deflation, or bordered pseudo-arclength
   continuation to search for disconnected roots.
3. In parallel at the mathematical level, re-derive whether the equation-(45)
   integration path is an equilibrium/stable path or a constrained path that
   may cross the present positive-mass boundary. This decision controls
   whether the current domain is a physical obstruction or an over-restrictive
   numerical gate.
4. Build error-controlled quadrature only after the component, branch
   selector, and anchor are defined.

## Evidence index

- Framework: `jensen_1z_framework.html`, Sections 9.2--9.4.
- Production profile: `invz_projected/invz_solve_point_ordered.m`.
- Bounded static closure: `invz_projected/invz_emt_static_ordered.m`.
- Deterministic outer map: `invz_projected/invz_ordered_outer_map.m`.
- Detailed execution record:
  `docs/execution/invzp_convergence_journal.md`, checkpoints 16--19.
- Transactional control:
  `docs/diagnostics/invzp_outer_wp2/wp2_hmf_node_transaction_census.mat`.
- Low-field \(M^2\) controls:
  `docs/diagnostics/invzp_outer_wp2/wp2_low_field_m2_census.mat` and
  `docs/diagnostics/invzp_outer_wp2/wp2_m2_asymptotic_check.mat`.
- Representation controls:
  `docs/diagnostics/invzp_representation_wp3/wp3_dominant_manifold_census.mat`
  and
  `docs/diagnostics/invzp_representation_wp3/wp3_closed_twolevel_landmarks.mat`.
- Saturation-tail control:
  `docs/diagnostics/invzp_integral_wp5/wp5_saturation_tail_census.mat`.
- 4 T continuation and endpoint:
  `docs/diagnostics/invzp_outer_wp2/wp2_4t_adaptive_node28_to27.mat` and
  `docs/diagnostics/invzp_outer_wp2/wp2_4t_adaptive_component_edge_audit.mat`.
- Rejected explanations and reopening conditions:
  `invzp_convergence_dead_ends.md`.

`diag_rev3_check.mat` is retained as a legacy frozen-state fixture, not as the
current solver result. Two maintained static-domain regression tests and four
reproducibility diagnostics still load its 1 T/3 T node payload. It may be
moved to a scoped fixture directory in a separate migration, but deleting it
alone would break those checks.
