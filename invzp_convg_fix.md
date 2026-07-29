# A formal replacement for the invZ projected branch-selection algorithm

**Canonical handoff document 2 of 2**

**Status:** recommended research program; not yet implemented

**Priority:** theoretical soundness and algorithmic completeness over compatibility,
runtime, or preservation of the current scalar closure

The current implementation and failure evidence are summarized in
[`invzp_convg_diagnosis.md`](invzp_convg_diagnosis.md). The protected Jensen
foundation remains
[`jensen_1z_framework.html`](jensen_1z_framework.html).

## 1. Decision

The formal fix is to replace the present chain of independently motivated
self-consistency maps with a **single source-dependent stationary thermodynamic
functional** for the full electronuclear lattice problem. The order parameter,
propagator/medium, self-energy, free energy, stability matrix, and response must
all be derivatives of the same truncated functional.

The replacement must:

1. retain the exact 136-state local Ho electronuclear Hamiltonian;
2. retain component, sublattice, Matsubara-frequency, and momentum labels until a
   reduction is proven;
3. use a skeleton-complete expansion at a declared order;
4. find all stationary branches in a certified search domain;
5. select equilibrium by the functional value, with metastability classified by
   the reduced Hessian; and
6. derive real-axis response from the same irreducible kernels.

This is more expensive than repairing the projected Picard solver. That cost is
intentional. A globally convergent nonlinear solver for the current equations
would still leave the central physical question unanswered: which of several
valid roots is the equilibrium state?

## 2. Non-negotiable standard

A candidate theory is not accepted because it produces a smooth spectrum or
fills the masked fields. It is accepted only when the following proof obligations
are met.

### 2.1 Common origin

There is one explicitly written functional \(\Gamma\), with stated independent
variables, normalization, source convention, double-counting terms, and
truncation order. Every state equation used by the code is obtained by
differentiating that \(\Gamma\). No equation may be inserted because it gives a
desirable resummation.

### 2.2 Skeleton completeness

At the retained order, all diagrams allowed by the same counting rule and
symmetries are present. In particular, a selected \(C_4\) contribution cannot be
kept while an equally ordered \(C_3\)-\(C_3\) contribution or its source
derivative is omitted.

### 2.3 Thermodynamic and response consistency

The free energy used to rank roots is the stationary value of \(\Gamma\). The
static susceptibility is the derivative of the same stationary solution with
respect to its physical source, including relaxation of every internal
variational variable. Dynamic response uses the corresponding irreducible kernel
and obeys causality, Kramers–Kronig, spectral-moment, and equal-time sum rules.

### 2.4 Completeness or explicit ambiguity

The algorithm either certifies that every stationary solution in a declared
bounded domain has been found, or returns `ambiguous/incomplete`. It may not
silently select the root produced by one seed.

### 2.5 Controlled limits

The theory must reduce to exact single-ion, noninteracting, symmetric-phase,
mean-field, Gaussian/RPA, and finite-cluster limits where those limits apply. The
coupling grid, Matsubara rank, local-state rank, momentum grid, real-axis
broadening, and skeleton order are independent convergence axes.

## 3. Why the current residual cannot simply be integrated

The present accepted-state contract is valuable as a numerical audit, but it was
assembled from a nested sequence of EMT, static closure, lambda, and self-energy
maps. In the tested positive-control coordinates its finite-difference Jacobian
is not symmetric; one measured pair is

\[
J_{12}\simeq-1.6649\times10^{-3},\qquad J_{21}=0.
\]

A continuously differentiable Euclidean gradient would have symmetric mixed
partials. Therefore the existing residual is not already
\(\nabla\Gamma\) in those coordinates.

This result has a limited but decisive interpretation:

- it **does not** prove that the underlying physical approximation has no
  variational representation in better variables;
- it **does** rule out declaring a line integral of the current residual to be a
  path-independent free energy; and
- it requires deriving new state equations from a functional rather than
  reverse-engineering a scalar score from the old maps.

The current `invz_deltaF_ordered` route discrepancy of about 13.7% is consistent
with this limitation. Thermodynamic integration may remain an independent audit,
but not the primary selector until path independence and the relevant derivative
identities have been proved.

## 4. Exact source construction

### 4.1 Microscopic variables

Start from the full local Hamiltonian already represented by the 136-state
electronuclear basis:

\[
H_{\rm loc}
=H_{\rm CF}+H_{\rm Z}+H_{\rm hf}
\text{any explicitly declared local terms}.
\]

Let \(O^a_i(\tau)\) denote the physical magnetic operators, with the compound
index \(a\) retaining Cartesian component, sublattice, and any additional channel
required by the chosen exact representation. The intersite interaction is a
fully indexed \(J^{ab}_{ij}\), including dipolar, exchange, and—when promoted—
off-diagonal-dipolar terms in one convention.

Do not project to a scalar doublet or a single uniform eigenmode at this stage.
Such a reduction is allowed only after a controlled-limit identity or a measured
error bound.

### 4.2 Generating functional

With all physical fields included in the microscopic action \(S\), introduce
auxiliary linear and symmetric bilocal sources:

\[
Z[\eta,Q]
=\operatorname{Tr}\,\mathcal T_\tau
\exp\!\left[
-S+\int\eta_a O^a
+\frac12\iint O^a Q_{ab}O^b
\right],
\qquad W[\eta,Q]=\log Z[\eta,Q].
\]

Using one frozen sign and normalization convention, define

\[
m_a=\frac{\delta W}{\delta\eta_a},\qquad
G_{ab}=2\frac{\delta W}{\delta Q_{ab}}-m_am_b.
\]

The exact double Legendre transform

\[
\Gamma[m,G]
=\eta\!\cdot\!m
+\frac12 Q:(G+mm)-W[\eta,Q]
\]

is stationary at zero auxiliary sources. This definition fixes what “free
energy,” “self-consistency,” and “response” mean before any approximation is
made.

With the dimensionless convention \(W=\log Z\) above, the stationary Helmholtz
free energy is \(\beta^{-1}\Gamma\) (or the corresponding grand potential after
including its physical sources). The implementation must not compare a
dimensionless \(\Gamma\) from one convention with an energy-valued functional
from another.

The transform exists only on its representable source domain. Use the
Legendre–Fenchel transform where a global transform is required and a branchwise
Legendre transform only where the source-response Jacobian is demonstrably
nonsingular. WP4 must identify the domain, singular boundary, and any loss of
one-to-one source inversion; a solver may not extrapolate \(\Gamma\) through an
unproved branch.

There is also an essential order of limits. A finite system at zero longitudinal
source has a symmetry-restored convex Gibbs functional, not two spontaneously
ordered equilibrium states. The broken-symmetry branches must be defined by
solving at nonzero conjugate source, taking the thermodynamic limit at fixed
source and declared dipolar shape/boundary convention, and only then taking the
source to zero. Exact finite clusters are coefficient and source-identity
oracles; they are not expected to exhibit a true spontaneous transition.
Metastable continuations require a stated constrained/coarse-grained potential
or branchwise analytic construction and must never be confused with the exact
convex equilibrium functional.

The implementation should derive signs and factors from automated source
derivatives of small exact systems; they must not be copied by analogy from a
fermionic formula.

### 4.3 Exact local generator and vertices

Construct the isolated-site \(W_{\rm loc}[\eta,Q]\) by exact diagonalization of
the full local Hamiltonian. Its connected time-ordered derivatives generate
\(C_1,C_2,C_3,C_4,\ldots\), including:

- arbitrary source field and ordered moment;
- bosonic Matsubara frequencies with exact conservation;
- all retained component indices;
- KMS/cyclic identities;
- coincident-frequency and degenerate-level limits; and
- analytic derivatives with respect to the physical field.

\(C_1\)–\(C_4\) are the minimum for the first non-Gaussian construction.
Derivatives of \(C_3\) and \(C_4\) that enter response may require \(C_5\) and
\(C_6\). Those vertices must be generated when the derivative counting requires
them; treating lower vertices as source-independent would break the common-origin
contract.

## 5. Controlled approximation

### 5.1 Expand the functional, not selected equations

Introduce an explicit bookkeeping parameter—coordination order, interaction
order, or another stated small parameter—and perform the linked-cluster/2PI
expansion of \(\Gamma[m,G]\). Enumerate vacuum skeletons algorithmically, reduce
them by symmetry, and attach:

- topology and automorphism factor;
- component and sublattice contractions;
- momentum/frequency conservation;
- retained-order label; and
- the source derivatives that generate state and response kernels.

The diagram manifest is part of the theory specification. A code path is rejected
if a term in its residual lacks a parent vacuum diagram or if differentiating the
manifest produces a term absent from the implementation.

Skeleton completeness and a common functional are necessary, not sufficient.
A finite 2PI truncation can still violate convexity, produce extra stationary
points, or fail spectral positivity/causality. The cluster, Hessian, controlled-
limit, and real-axis gates below remain independent rejection tests;
“\(\Phi\)-derivable” is not itself a promotion verdict.

### 5.2 Exact local reference, direct lattice propagation

Use the exact local problem as the reference and keep the lattice propagator
\(G^{ab}(\mathbf q,i\omega_n)\) explicitly momentum resolved in the first formal
implementation. This avoids imposing a scalar effective-medium closure before a
corresponding log functional and double-counting correction have been derived.

An EMT/DMFT-like compression can be reconsidered later, but only if:

1. its functional is stated;
2. differentiating it reproduces the impurity and lattice Dyson equations;
3. its local/lattice double counting is explicit; and
4. its finite-grid poles are distinguished from physical singularities.

### 5.3 Required same-order terms

For the first non-Gaussian order, the manifest must test explicitly whether the
following occur at the same declared order:

- a local \(C_4\) vertex closed by two lattice lines;
- two \(C_3\) vertices linked by three lattice lines;
- tadpole/counterterm structures induced by a nonzero moment;
- source derivatives of every retained vertex;
- sublattice- and component-exchange partners; and
- the vacuum terms whose derivatives supply the state equations.

The list is a classification obligation, not permission to assume the
coefficients. Coefficients and signs come from the exact source construction and
finite-cluster matching.

## 6. Independent coefficient oracles

Before solving the infinite lattice, build exact two-, three-, and where needed
four-site clusters using the same local Hamiltonian and interaction convention.
Expand their exact \(\log Z\), moments, connected correlators, and susceptibilities
in the coupling parameter.

For every retained monomial:

1. compute its coefficient by exact finite-cluster differentiation;
2. compute it independently from the diagram generator;
3. compare before any resummation or numerical fitting; and
4. freeze the matching coefficient and convention in a regression fixture.

This gate has already exposed a material failure in an earlier varied-\(D\)
skeleton attempt: at a symmetric fixture an exact static mixed-chain coefficient
was \(-0.0944822\), while the candidate gave \(-0.366694\). Because this is an
immutable series coefficient, the candidate was rejected rather than tuned.

The same standard applies to every proposed resummation. Agreement of a final
critical field cannot compensate for a wrong low-order coefficient.

## 7. Stationarity and phase-selection algorithm

### 7.1 Discovery

Discretize only after the functional and derivatives are fixed. Solve

\[
\frac{\delta\Gamma}{\delta m}=0,\qquad
\frac{\delta\Gamma}{\delta G}=0
\]

with analytic or automatic-differentiation Jacobian-vector products. Use
globalized Newton/trust-region methods for local convergence, pseudo-arclength
continuation through folds, and deflation or symmetry-aware multi-start searches
to discover disconnected stationary components.

These methods are discovery tools. None decides phase stability.

### 7.2 Certification

Partition a physically bounded domain into continuation slabs. Use interval
Newton/Krawczyk operators, validated linear solves, or an equivalent certified
root-counting method to:

- enclose every nonsingular stationary point;
- detect unresolved boxes near singular bifurcations;
- certify empty boxes;
- preserve symmetry-related multiplicities; and
- report an incomplete domain rather than choose from an unproved list.

For the large discretized system, a practical hierarchy is acceptable:
certify reduced critical eigenspaces and tail bounds analytically, then apply
interval methods to the reduced border system. The reduction itself needs a
residual and spectral-gap bound.

### 7.3 Stability

At a stationary point, eliminate internal propagator variations before
classifying order-parameter stability. If the Hessian is blocked as

\[
H=
\begin{pmatrix}
\Gamma_{mm}&\Gamma_{mG}\\
\Gamma_{Gm}&\Gamma_{GG}
\end{pmatrix},
\]

then the relaxed inverse susceptibility is the Schur complement

\[
\chi^{-1}_{\rm relaxed}
=\Gamma_{mm}
-\Gamma_{mG}\Gamma_{GG}^{-1}\Gamma_{Gm},
\]

subject to the declared sign convention. Testing only \(\Gamma_{mm}\) while
holding \(G\) fixed would misclassify a coupled instability.

After removing exact symmetry zero modes:

- a positive reduced Hessian identifies a local minimum/metastable phase;
- a negative direction identifies a saddle or maximum;
- a zero eigenvalue triggers bifurcation analysis and finer certification.

### 7.4 Equilibrium

Evaluate the **same** \(\Gamma\), including every constant and
double-counting term, at all certified local minima. The global minimum is the
equilibrium phase; other minima are metastable. Equal minima locate a first-order
boundary, while a vanishing relaxed Hessian eigenvalue locates a continuous
boundary.

If the stationary search or free-energy error bars cannot order two candidates,
the result is `ambiguous`, not a seeded guess.

## 8. Response from the same theory

Differentiate the stationary equations with respect to the physical magnetic
source. This automatically includes vertex motion and the relaxation of \(G\).
The Matsubara susceptibility must agree with the inverse reduced Hessian in the
static limit.

For real frequency:

1. continue the irreducible kernels or solve their spectral representation;
2. rebuild the retarded response with the same skeleton content;
3. verify analyticity in the upper half-plane and nonnegative physical spectral
   weight where required;
4. check Kramers–Kronig and equal-time moment sum rules; and
5. show convergence independently in Matsubara cutoff, spectral quadrature,
   broadening, local rank, and momentum grid.

Real-axis broadening is a resolution parameter, never a repair for a missing
imaginary-axis state.

## 9. Work packages and gates

### WP0 — Freeze the microscopic contract

- Freeze the Hamiltonian, energy units, Fourier convention, self-site
  subtraction, Ewald/shape convention, component basis, and coupling-grid digest.
- Freeze the finite-volume sequence and the order of the thermodynamic,
  symmetry-breaking-source, and zero-frequency/momentum limits.
- Record which ODD, exchange, demagnetization, and hyperfine terms are included.
- Keep `jensen_1z_framework.html` immutable as a comparison derivation; do not
  assume every hybrid resummation in the legacy code must survive.

**Gate:** two independent constructors reproduce all frozen microscopic anchors
and the same serialized coupling set.

### WP1 — Exact local source oracle

- Generate \(W_{\rm loc}\), \(C_1\)–\(C_4\), and derivative-required higher
  vertices from all 136 states.
- Implement degeneracy-safe divided differences or contour formulas.
- Verify KMS, permutation, Hermiticity, static derivative, and high-frequency
  identities.

**Gate:** analytic/source derivatives match high-precision finite differences and
direct Lehmann sums on random and degenerate fixtures.

### WP2 — Strict linked-cluster functional

- Enumerate the vacuum diagrams at the declared first non-Gaussian order.
- Derive moment, propagator, and response equations symbolically from that
  manifest.
- Produce a machine-readable provenance map from each residual term to its parent
  diagram.

**Gate:** independent symbolic and numerical differentiation agree term by term.

### WP3 — Exact finite-cluster matching

- Match every low-order coefficient against exact clusters.
- Include ordered, symmetric, zero-temperature, finite-temperature, and
  component-mixing fixtures.

**Gate:** exact coefficients agree before fitting; any immutable mismatch rejects
the functional.

### WP4 — Nonlinear partial Legendre transform

- Construct and validate the bilocal effective action for the non-Gaussian local
  reference.
- Prove convexity/domain statements that are actually used by the solver.
- Resolve the currently missing nonlinear relation between source, connected
  propagator, and local irreducible vertices.

**Gate:** forward source solve followed by the Legendre inversion returns the
source, \(m\), \(G\), and functional value within certified errors.

### WP5 — Full tensor 2PI assembly

- Assemble the q-resolved, full-component functional with all same-order
  skeletons and double-counting terms.
- Differentiate it to produce stationarity, Hessian, and source-response kernels.

**Gate:** mixed derivatives agree, exact limits reduce correctly, and no residual
term lacks diagram provenance.

### WP6 — Global stationary solver

- Add globalized Newton/trust-region iteration, pseudo-arclength continuation,
  deflation, and certified slab/root counting.
- Preserve every stationary branch and bifurcation event in a versioned trace.

**Gate:** synthetic fold/cusp fixtures and exact finite clusters have complete
root counts; seed order does not change the certified set.

### WP7 — Equilibrium selection

- Classify the relaxed Hessian.
- Compare stationary functional values with numerical error bounds.
- Return equilibrium, metastable, unstable, or ambiguous status.

**Gate:** exact finite-cluster phase rankings and susceptibilities are recovered;
coexisting minima are ordered independently of continuation direction.

### WP8 — Independent convergence

- Converge local rank, Matsubara rank, q grid/Ewald range, diagram order, nonlinear
  tolerances, and certification bounds separately.
- Perform continuum extrapolation only when its asymptotic regime is demonstrated.

**Gate:** the claimed precision has an error budget with no cancellation between
unconverged axes.

### WP9 — Real-axis response

- Derive and implement the retarded kernel from the same functional.
- Test causality, positivity, Kramers–Kronig, equal-time, and spectral moments.

**Gate:** static response agrees with the relaxed functional Hessian and all
response identities close within the declared discretization error.

### WP10 — Production promotion

- Introduce the formal solver behind a new explicit theory flag.
- Compare with Jensen, strict-ring, exact-cluster, and experimental observables
  without tuning to fill masks.
- Change the default only after every preceding gate passes.

**Gate:** no undocumented fallback, seed-dependent phase label, pole floor, or
masked-to-bare substitution remains.

## 10. Disposition of existing approaches

### Safeguarded Aitken

Retain as a default-off accelerator for the legacy `resummed` solver. It has a
valid narrow contract and useful regression evidence. It is not part of the
formal phase-selection proof.

### Warm continuation, Newton repair, and branch graphs

Retain as discovery and diagnostic tools. They provide initial guesses and reveal
folds, disconnected components, and local solver weaknesses. Their output must
feed the certified stationary search; it must not directly assign equilibrium.

### Smooth-\(r\) objective

Retire as a proposed final selector. It has no source-functional derivation and
can choose a smooth metastable or disconnected component.

### Strict ring/local functional

Retain its exact local generator, algebraic checks, cluster oracles, and
high-frequency machinery. Do not promote the current strict-ring lattice theory:
its Gaussian pole/no-state interval shows that removing the legacy resummed pole
is insufficient.

### Earlier varied-\(D\) and quadratic bilocal prototypes

Retain only as falsification evidence and oracle infrastructure. The varied-\(D\)
candidate has a failed immutable coefficient. The quadratic bilocal prototype is
valid only in bounded \(C_3=0\), finite-cutoff fixtures and is not the required
nonlinear Legendre construction.

### Potthoff/self-energy-functional shortcuts

Do not use a reference-system self-energy functional merely because it supplies
a scalar score. Its universality and representability assumptions must be proved
for this non-Gaussian spin/electronuclear action and chosen variables. It may be
considered only after it passes the same exact-cluster and mixed-derivative gates.

## 11. Definition of done

The root problem is formally fixed only when all of the following are true:

1. one versioned functional specification generates state, free energy,
   stability, and response;
2. all retained diagrams and coefficients have machine-auditable provenance;
3. exact local and finite-cluster coefficient gates pass;
4. mixed derivatives and controlled limits pass;
5. the solver certifies all stationary states in its declared domain or returns
   `ambiguous`;
6. phase selection uses the stationary functional and relaxed Hessian, not seed
   history;
7. the 1.5 T folds, 3.825 T failure region, 4.05 T coexistence, and QCP window are
   all resolved without pole floors or branch heuristics;
8. the static response equals the appropriate functional derivative;
9. real-axis response passes causality and sum-rule gates; and
10. every numerical and truncation axis has an independent error budget; and
11. the symmetry-breaking and dipolar thermodynamic limits are declared and
    numerically stable under the prescribed order.

If a skeleton-complete candidate still has a no-state interval, that is a failure
of the candidate approximation at that order. The correct response is to derive
the next complete order or change the controlled expansion—not to clip a
denominator, inject a preferred root, or hide the interval.

## 12. Recommended first action

Begin with WP0–WP3, not with another production solver modification. The most
valuable immediate deliverable is a small exact-cluster coefficient harness tied
to an automatically generated diagram manifest. It can falsify an unsound
functional cheaply relative to a full lattice implementation and prevents a
large computation from concealing a wrong sign, symmetry factor, or omitted
same-order topology.

Only after those coefficients close should effort move to the nonlinear bilocal
Legendre transform and the full q-resolved tensor 2PI solver.

## 13. Theoretical orientation

The recommended construction follows the stationary-functional logic of
Luttinger–Ward/Baym–Kadanoff and the two-particle-irreducible effective action,
adapted explicitly—not by formula substitution—to a finite-dimensional quantum
spin/electronuclear reference:

- J. M. Luttinger and J. C. Ward, *Phys. Rev.* **118**, 1417 (1960),
  doi:10.1103/PhysRev.118.1417.
- G. Baym, *Phys. Rev.* **127**, 1391 (1962),
  doi:10.1103/PhysRev.127.1391.
- J. M. Cornwall, R. Jackiw, and E. Tomboulis, *Phys. Rev. D* **10**, 2428
  (1974), doi:10.1103/PhysRevD.10.2428.

These references provide architecture, not a ready-made LiHoF4 formula. The
actual functional, signs, symmetry factors, local vertices, and double counting
must be derived and verified by the source and cluster program above.

The independently supplied McKenzie–Stamp spin-bath paper is a particularly
useful WP1/WP2 normalization oracle. Its Hubbard–Stratonovich construction for
its stated quantum-Ising/spin-bath Hamiltonian organizes the action in exact
single-site connected cumulants and recovers RPA at quadratic order, so its
source signs, beta factors, and \(C_2/C_3/C_4\) organization should be
cross-audited against the new local generator. It does **not** furnish the
nonlinear bilocal/2PI Legendre functional for the current
\((\Sigma,K_0)\) equations or a value that can be evaluated on the existing
roots; it cannot authorize a retrospective free-energy score.

- R. D. McKenzie and P. C. E. Stamp, *Phys. Rev. B* **97**, 214430 (2018),
  doi:10.1103/PhysRevB.97.214430. A copy is retained under `References/`.
