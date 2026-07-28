# Rigorous and system-general 1/z extension

**Status:** active isolated scalar prototype; not connected to production
**Recorded:** 2026-07-22  
**Motivation:** the Task 6b closed-model J 2.34 comparison leaves a stable 13.65% route mismatch. The independent assessment (adjudicated against the published Jensen paper) is folded into §7b of the Stage-2 plan ([docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md](docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md)) and preserved in the consolidated [ordered-leg diagnosis](invzp_convg_diagnosis_Claude.md); the standalone review and QCP diagnosis notes were retired after consolidation.

**Feasibility audit, 2026-07-28:** `invzp_functional_derivation_audit.md` tested whether the current
resummed/Jensen implementation already contains the needed functional. It found a
branch-conditioned one-dimensional `Phi_path` and an isolated scalar-EMT log potential, but no common
functional for the coupled ordered equations and no thermodynamic selector for their multiple roots.
The present Newton residual also fails the symmetry/zero-pattern test for a direct or diagonally
rescaled gradient. WP0 therefore remains the correct entry point; production equations are not to be
modified during the derivation.

WP0 is now closed at the strict unicyclic `O(1/z)` scope by
`invzp_functional_wp0_spec.md` and `invzp_functional_wp0_ring_derivation.md`. The complete retained
vacuum is the `C2` ring family; its source derivatives generate the ordered `C3` one-point term and
the symmetric `C4+C3^2` Hessian. Exact nonzero-source two-site coefficients and the leading Jensen PM
dynamic coefficient pass independent checks. An explicit `C4` vacuum belongs to the deferred
`O(1/z^2)` bicyclic scope, where all `C3-C3` classes and local `C5/C6` derivatives must also be kept.
The isolated `invz_functional/` prototype implements the closed scope without EMT or production
dispatch.

A first read-only production-input pilot now maps the actual scalar coupling multiset and transverse
electronic doublet into that strict functional. At `T=0.1 K, Bx=1.5 T` it finds stable ordered minima
with converged Matsubara and lattice-grid root enumeration, but its electronic mode is about `85 GHz`
and cannot represent the observed `0--6 GHz` hyperfine spectrum. It also has a Gaussian-pole
no-stationary-state interval around `3.2--3.5 T`. These are stop signals for production wiring, not
failures to be hidden: WP1 must next supply a source-biased electronuclear local reference, after
which survival of the pole interval would trigger WP4 rather than an ad hoc regularization.

That WP1 local gate is now implemented for the full 136-state manifold with
`transverse_mf='none'`, the only scalar choice that keeps the local free energy and response in one
Hamiltonian without adding transverse conjugate variables. It matches the independent local
susceptibility path to `1.9e-13`, gives a cutoff-stable `1.5 T` minimum on the production scalar
couplings, and moves the ordered solutions to `3.6--4.6 T`. A narrower no-state interval remains at
`4.7--4.85 T` before the stable symmetric `4.9 T` state. WP4 is therefore activated at the
specification level in `invzp_functional_wp4_skeleton_spec.md`; no static-only pole patch is
authorized.

The WP4 two-level entry fixtures are now complete.  Exact dynamic `C2/C3/C4` and their 1PI
amputations pass analytic source, permutation, beta-delta, and degeneracy gates.  Centred two- and
three-site clusters freeze the leading `C3-C3` coefficient and the nonzero-source mixed-`C4`
cancellation coefficient.  Both dipolar backends can now expose the exact q-resolved Hermitian
coupling pages behind an optional output while leaving ordinary production calls unchanged.  The
next authorized task is the smallest isolated varied-covariance skeleton, not production wiring.

## 1. Decision summary

The present Jensen ordered machinery should remain available as an explicitly labelled, inexpensive approximation. It should not become the universal thermodynamic foundation of the code.

A future general extension is possible, but it requires a new derivation and solver whose self-energy, order-parameter equation, susceptibility, internal energy, and free energy all come from one truncated generating functional. The project should retain the useful, systematic part of Jensen's `1/z` idea while removing the HoF3-motivated static-elastic shortcuts.

The recommended hierarchy is:

1. Build a strict, non-resummed first-order reference from an explicit free-energy/effective-action functional.
2. Validate it on exactly solvable finite clusters and against known two-level limits.
3. Only then introduce a self-consistent skeleton resummation derived from the same functional.
4. Generalize the local state from two levels to the full electronuclear manifold through numerically generated connected vertices.

This is an extension project rather than a repair to Stage 2. Until it exists, `jensen`/`jensen_static` results must carry their current approximation scope, and J 2.34 must not be represented as a universal percentage-accuracy guarantee.

## 2. What is general in Jensen, and what is not

Jensen's paper treats a material deliberately suited to the approximation: HoF3 has two low singlets well separated from the remaining levels and predominantly dipolar intersite coupling. The paper computes the first correction in a high-density `1/z` expansion. See [Jensen's published paper](<References/Z^-1 renormalization of the mean field behavior of the dipole-coupled singlet-singlet system HoF3.pdf>) and the [author's preprint](https://arxiv.org/abs/cond-mat/9312060).

| Layer | Generality assessment |
|---|---|
| Exact local diagonalization and Lehmann response | General for any finite local Hilbert space and operator set |
| Expansion around a single-site/mean-field reference in powers of `1/z` | Systematic when the coordination/high-density power counting is controlled |
| Connected local cumulants and lattice contractions | General; local vertices can be generated from energies and operator matrices |
| Effective-medium closure | Potentially reusable, but its functional derivability must be proved rather than assumed |
| Reduction to `Delta`, `M²`, `m`, and `n01` | Two-level/singlet-singlet specific |
| J 2.26–2.27 moment-form ordered self-energy | Specific to the two-level algebra and the approximations made in its derivation |
| J 2.28–2.29 zero-frequency elastic closure | Not universal; it neglects dynamical elastic structure and elastic-peak broadening |
| The `tanh` elastic-weight resummation | A higher-order stabilization chosen for the HoF3 low-temperature problem, not a strict first-order result |
| The repository's full-weight/two-level-`xi` hybrid | A repository modelling choice that preserves the PM boundary; it is not Jensen's original closed theory |

The future core must therefore contain no required dependence on the two-level fields `tl.Delta`, `tl.M2`, `tl.m`, or `tl.n01`. A two-level adapter may remain as an analytic benchmark.

## 3. Required standard of rigor

“Rigorous” here should mean **controlled and thermodynamically consistent at a declared truncation**, not an exact solution of an infinite quantum lattice.

The target is a source-dependent functional, schematically

\[
\Gamma[m,G;h,T,J] = \Gamma_0[m,G;h,T] + \Phi_{1/z}[m,G;J],
\]

with the retained diagrams and their `1/z` power counting listed explicitly. The physical state is stationary:

\[
\frac{\delta\Gamma}{\delta m}=0,\qquad
\frac{\delta\Gamma}{\delta G}=0.
\]

The exact signs and numerical prefactors must be derived in the repository's `G = -chi` and ferromagnetic-positive-`J` conventions; the equations above are architectural, not a transcription specification.

At the stationary point, all reported observables must be obtained from this same object:

\[
m=-\frac{\partial F}{\partial h},\qquad
\chi=-\frac{\partial^2F}{\partial h^2},\qquad
U=\frac{\partial(\beta F)}{\partial\beta}.
\]

This is the relevant principle behind Baym–Kadanoff/Φ-derivable conserving approximations: a common functional ties static thermodynamics to the self-energy instead of constructing them through separate approximations. See [Baym, *Phys. Rev.* 127, 1391 (1962)](https://doi.org/10.1103/PhysRev.127.1391) and [Potthoff's functional review](https://arxiv.org/abs/cond-mat/0406671).

For spin/Hubbard operators, the electronic Luttinger–Ward formulas must not be copied mechanically. The functional should be derived from the exact local source-generating functional or an equivalent Hubbard–Stratonovich/linked-cumulant construction.

### Non-negotiable consistency rules

- Every retained self-energy diagram has its corresponding vacuum/free-energy diagram.
- Field and temperature derivatives include the explicit source dependence of the local vertices.
- Elastic and inelastic contributions use one Lehmann/KMS representation; no zero-frequency replacement is made only in one observable.
- Tadpoles are removed through a declared centred operator `delta J = J - <J>` and a consistent mean-field subtraction.
- A resummation is permitted only if it follows from stationarity of a declared truncated functional.
- Static response and real-frequency spectra must be continuations of the same Matsubara solution.
- Numerical projection of the local Hilbert space is an independently converged approximation, not silently combined with the `1/z` truncation error.

## 4. Recommended mathematical route

### 4.1 Exact local generating data

For a local Hamiltonian in a source field,

\[
H_{\mathrm{loc}}(h)=H_{\mathrm{CF}}+H_{\mathrm{hf}}-\sum_\mu h_\mu J_\mu,
\]

diagonalize the local problem and construct connected Matsubara vertices of the centred operators:

- `C1`: the ordered moment;
- `C2`: the full dynamic susceptibility, including the degenerate/static contribution;
- `C3`: the source/condition-space vertex needed when differentiating through an ordered state;
- `C4`: the connected four-point vertex entering the leading fluctuation self-energy and vacuum diagrams.

Use Lehmann divided differences with explicit KMS limits so degeneracies and `omega_n = 0` are handled analytically. Do not split out an “elastic pole” and later approximate its dynamics.

### 4.2 Strict first-order reference

Before any Dyson-style self-consistency, evaluate the complete set of diagrams at the declared first order using the same reference propagator. This gives a transparent power-counting oracle. It may be less numerically smooth near criticality, but it is the correct object against which later resummations should be judged.

The theory note for this stage must contain:

1. the Hamiltonian and sign conventions;
2. the definition of `z` or the high-density scaling for long-range dipolar couplings;
3. every retained vacuum, two-point, and one-point diagram with its symmetry factor;
4. the relationship between the vacuum functional and the self-energy/field vertex;
5. the first omitted diagram class;
6. the expected coupling/`1/z` power of the physical error.

No numerical acceptance exponent should be chosen before this diagram audit is complete. The zero-mean synthetic fixture happens to start at `J²`; that does not prove that every next-order residual must be `J³` rather than `J⁴` or another symmetry-selected power.

### 4.3 Functional self-consistency

After the strict-order implementation passes its exact references, construct an optional self-consistent skeleton version by replacing bare internal lines with stationary dressed lines **inside the same Φ functional**. This may improve behaviour near the transition while preserving thermodynamic identities.

The current mixture—Dyson-resummed dynamics, a separately resummed static elastic weight, and a first-order energy expression—is specifically what the extension must avoid.

### 4.4 Thermodynamics and response

Evaluate the stationary functional directly to obtain `F`. Use the envelope theorem when differentiating a converged stationary solution: implicit derivatives of stationary variables cancel, leaving the explicit field or temperature derivative, provided the stationarity residual is below the derivative error budget.

Implement two independent numerical audits:

- finite differences or complex-step derivatives of the evaluated `F`, where analyticity permits;
- analytic moment, susceptibility, and internal-energy formulas derived from the same diagrams.

Their agreement is a numerical correctness gate. The difference between the truncated functional and exact finite-cluster results is the physical truncation error. These are distinct error budgets and must not be combined into one percentage.

## 5. Repository reuse and isolation plan

The project should begin as a new isolated module—provisionally `invz_functional/`—and consume stable lower-level code without changing the production solvers.

### Reuse directly

- `invz_common/invz_single_ion.m`: exact crystal-field, Zeeman, and optional 136-state electronuclear diagonalization.
- `invz_common/invz_chi0z.m`: full multilevel Lehmann two-point response, including its explicit elastic convention.
- `invz_common/invz_matsubara.m`: frequency grid and weights, after confirming that signed-frequency needs are supported by the new contractions.
- `invz_tensor/invzt_kernels.m`: divided-difference/KMS kernels and exact degeneracy limits.
- `invz_tensor/invzt_vertex3.m`: existing connected three-point engine.
- `invz_tensor/invzt_vertex4.m`: dense connected four-point reference, already checked against an independent high-precision oracle.
- `invz_tensor/invzt_emt_matrix.m`: candidate lattice/medium engine, subject to a functional-derivability audit.
- `invz_tensor/invzt_ordered_vertex_basis.m`: candidate field-adapted truncation mechanism and its subspace diagnostics.

### Reuse only as comparison backends

- `invz_common/invz_sigma.m` and `invz_sigma_ordered.m`;
- `invz_common/invz_gstat_ordered.m`;
- `invz_projected/invz_emt_static_ordered.m`;
- `invz_projected/invz_hmf_ordered.m`;
- `invz_projected/invz_deltaF_ordered.m`.

These encode the Jensen/moment-form route being tested; they must not define the new functional.

### Important existing limitation

`invz_tensor/invzt_sigma_tensor.m` already builds a genuine dynamic connected-four-point self-energy and is the strongest starting point. It is **not yet the requested solution**, because it supplies no matching vacuum/free-energy functional. Its symmetric-bracket Dyson resummation and hybrid `chi_base` map must be derived from a common functional before being admitted into the rigorous backend.

The existing `a3d` rank-16 production point costs approximately 7.4 hours. A new functional implementation must therefore start with small local Hilbert spaces and cached vertex data. Full 136-state dense `C4` is not a reasonable first target.

### Suggested new interfaces

Names are provisional; the data contracts matter more than filenames.

| Function | Responsibility |
|---|---|
| `invzf_local_state` | Wrap exact local energies, populations, centred operators, sources, and basis provenance |
| `invzf_vertices` | Produce/cache `C2`, `C3`, and `C4` with symmetry and truncation diagnostics |
| `invzf_diagrams_o1z` | Evaluate the enumerated first-order vacuum and self-energy diagrams |
| `invzf_functional` | Evaluate `Gamma`, `F`, and individual diagram contributions |
| `invzf_stationary` | Solve the common `m,G` stationarity equations and report residuals |
| `invzf_observables` | Derive `m`, `chi`, and `U` from the stationary functional |
| `invzf_cluster_exact` | Exact finite-cluster partition function and response oracle for small models |
| `invzf_continue` | Optional real-axis continuation, added only after Matsubara thermodynamics closes |

Every returned solution should carry at least:

- the Hamiltonian/source and lattice hashes;
- local basis definition and captured population/variance;
- Matsubara and lattice grids;
- functional and stationarity residuals;
- separate frequency, lattice, local-rank, and solver error estimates;
- diagram-resolved contributions to `F`, `U`, and `Sigma`;
- the declared theory mode: strict-order or functional-self-consistent.

## 6. Work packages and exit gates

### WP0 — derivation and feasibility, no production edits

Deliver a theory specification containing the diagram list, power counting, signs, symmetry factors, and functional derivatives. Prove or disprove that the current scalar and tensor EMT closures can be obtained from the proposed functional.

**Exit gate:** an independent hand calculation on the closed two-level model reproduces the PM Jensen first-order self-energy and free-energy correction without using J 2.28–2.29.

**2026-07-28 status:** passed at the strict unicyclic scope. The exact pair free-energy coefficient,
nonzero-source derivatives, and local Matsubara `j^2` coefficient reproduce the ring/Jensen expansion.
The full resummed EMT primitive is explicitly not claimed as the strict vacuum functional.

### WP1 — exact local vertices

Generalize the existing tensor vertex engine into a stable local-vertex service. First support a two-level model, then three-level fixtures, then a field-adapted multilevel basis.

**Exit gates:**

- basis-unitary invariance;
- KMS and signed-frequency transpose identities;
- exact degeneracy/Hermite limits;
- connected-cumulant subtraction and equal-time sum rules;
- agreement with the existing JSON high-precision oracle;
- rank-ladder convergence of both static susceptibility and connected variance.

**2026-07-28 status:** the scalar two-level service now includes analytic `C2`, exact
frequency-labelled `C3/C4`, all static/source identities, connected beta-delta subtraction, and
dynamic 1PI `gamma2/gamma3/gamma4`.  A full source-biased electronuclear `C2` oracle with verified
nested-stencil source derivatives is also implemented.  Component-labelled electronuclear `C3/C4`
and multilevel rank gates remain open.

### WP2 — strict functional on a scalar closed model

Implement the vacuum diagrams and their derivatives without self-consistent resummation.

**Exit gates:**

- `J = 0` gives the exact local partition function;
- differentiating the functional reproduces the separately implemented strict-order self-energy and moment vertex;
- field and temperature thermodynamic identities close to the measured quadrature/derivative error;
- expansion coefficients match an exact finite-cluster calculation as the coupling scale tends to zero.

**2026-07-28 status:** initial scalar implementation exists in `invz_functional/`. It contains the
finite-mode ring, analytic source derivatives, common `(m,h)` functional/Hessian, and a stationary
root enumerator that compares stable roots by functional value. The stationary envelope identity
`U=d(beta F)/d beta` and the independent field identities for `m` and `chi` pass. The scalar
Matsubara free-energy tail has a rigorous bound plus a Richardson audit. Complete stationary-root
re-enumeration across Matsubara cutoffs and production BZ grids is implemented, as is a
single-generator electronuclear local reference.  The non-Gaussian skeleton is deliberately a
separate WP4 resummation rather than a silent addition to this strict ring theory.

### WP3 — exact finite-cluster oracle

For the two-level Hamiltonian, build small clusters with an explicit coupling matrix satisfying the same no-self-site convention. Compute

\[
F=-T\log\operatorname{Tr}e^{-\beta H},\quad m,\quad U,\quad \chi
\]

directly. Extract weak-coupling coefficients by symmetric coupling scans or analytic derivatives where possible.

Start with dimensions small enough for dense diagonalization. Increase site count or use sparse methods only after the oracle itself is independently verified. A full 136-state multi-site cluster is not an initial deliverable.

**Exit gate:** the functional reproduces the exact cluster's coefficient at every retained order; the residual follows the analytically predicted first-omitted power.

**2026-07-28 status:** dense scalar cluster oracles through eight sites are implemented. The pair
`j^2` coefficient passes, and a four-site ring fixture gives exact-minus-functional powers
`epsilon^4` at the symmetric point and asymptotically `epsilon^3` at nonzero source, matching the
first omitted `C4` and `C3-C3` classes.  Fixed-reference centred clusters additionally isolate the
leading nonzero-source `C3-C3` coefficient and the mixed-chain connected-`C4` coefficient required
to grade the 1PI cancellation.

### WP4 — Φ-derived self-consistent solver

Introduce stationary dressed lines and compare strict-order versus self-consistent results.

**2026-07-28 specification:** the first candidate varies the full covariance and builds its skeleton
from the subtracted nonlocal return `R_ij=D_ij-delta_ij*C2`. The one-`C4` double-return and the
same-order `C3-C3` three-line core are inseparable. The latter is momentum/sublattice dependent and
cannot be reconstructed from `Jnu_flat`; q-resolved coupling matrices or an equivalent real-space
transform are required.  Those matrices and all two-level entry fixtures are now available, so the
smallest isolated implementation is authorized.  See `invzp_functional_wp4_skeleton_spec.md`.

**Exit gates:**

- the stationarity norm converges under a safeguarded solver (a 2PI functional may be a saddle, so monotonic decrease is not assumed without a separate proof);
- results are independent of warm/cold start within a tracked physical basin;
- `m = -dF/dh`, `chi = dm/dh`, and `U = d(beta F)/d beta` agree within the numerical budget;
- no route-specific static substitution or partial resummation remains.

### WP5 — multilevel and tensor integration

Move from the scalar closed model to the field-adapted multilevel/tensor machinery. Reuse `invzt_vertex3/4`, but run an explicit rank ladder rather than accepting rank 16 from static `chi_share` alone. The existing record shows that rank-16 `chi_share` can be about 0.98 while connected-variance coverage is only 0.665–0.838.

**Exit gates:** all reported observables have separate, converged local-rank, Matsubara, and lattice-grid uncertainties; the result is invariant under rotations within the retained spectral subspace.

### WP6 — production comparison and continuation

Compare the new backend with:

- the exact closed-cluster oracle;
- the current Jensen-static backend;
- Jensen's published HoF3 low-temperature 2–3% Eq. 2.34 check;
- the current LiHoF4 projected and tensor phase-boundary/spectral anchors.

Only after Matsubara thermodynamics closes should real-axis continuation and production sweeps be connected.

## 7. Validation hierarchy

The new project should use several distinct gates rather than one percentage comparison.

1. **Algebraic gates:** KMS, connectedness, symmetry, basis invariance, zero-coupling limits.
2. **Stationarity gates:** residuals of every varied object, with no stale tuple reuse.
3. **Thermodynamic identity gates:** derivatives of the same computed functional agree to the numerical error budget. These should approach solver precision, not a 10–15% physics tolerance.
4. **Order gates:** exact-cluster residuals scale with the first omitted order derived in WP0.
5. **Discretization gates:** independent Matsubara cutoff/density, lattice grid, temperature grid, and field grid convergence.
6. **Projection gates:** local-rank ladders track `F`, `m`, `U`, `chi`, and vertex norms—not `chi_share` alone.
7. **External benchmark gates:** Jensen/HoF3 and existing LiHoF4 anchors, each with its own physical and numerical error budget.

The two J 2.34 routes can remain as an audit, but the primary new comparison should be between analytic derivatives of the common functional and independent numerical derivatives. If those disagree materially, it is an implementation defect by construction.

## 8. Open theory decisions that must be resolved before coding

1. **Expansion variable:** precise `1/z`/high-density scaling for a long-range dipolar, multi-sublattice lattice.
2. **Functional representation:** strict linked-cluster effective action versus a 2PI/Φ functional; the strict implementation should come first either way.
3. **EMT status:** whether `invz_emt_scalar`/`invzt_emt_matrix` is stationary under an explicit medium functional, and what vacuum term generates its `K` relation.
4. **Mean-field subtraction:** exact handling of `J(0)`, `J(ii)=0`, demagnetization, ODD-mediated couplings, and double counting.
5. **Ordered zero mode:** treatment of spontaneous-symmetry limits and diagonal/degenerate matrix elements without an ad hoc static pole split.
6. **Local truncation:** fixed-rank field-adapted basis, adaptive rank, or another compressed representation; required convergence evidence for each.
7. **Resummation policy:** strict order only, or a skeleton self-consistent mode in addition; never mix them silently.
8. **Analytic continuation:** whether spectra are derived directly from the same diagrams or continued numerically from Matsubara data.
9. **Computational target:** projected scalar production, full tensor production, or a controlled benchmark engine only.

## 9. Principal risks and mitigations

| Risk | Mitigation |
|---|---|
| A self-energy is implemented without the matching vacuum functional | WP0 blocks coding until diagram/functional derivatives are explicit |
| Dense multilevel `C4` is computationally prohibitive | small-model-first design, caching, symmetry, rank ladders, and hard allocation/work budgets |
| A fast factorization changes KMS structure | keep the dense oracle authoritative; the existing unproved factorized vertex path remains disabled |
| Local projection looks converged in `chi` but not in vertices | gate connected variance, `C3/C4` norms, and final observables over rank |
| Self-consistency converges to a spurious basin | cold/warm/multistart checks plus stationarity/free-energy comparison |
| Criticality amplifies numerical derivatives | analytic functional derivatives, envelope-theorem residual checks, adaptive finite-difference windows |
| ODD/transverse feedback is double counted | include every interaction in the WP0 Hamiltonian and diagram inventory before contraction |
| A numerically fitted tolerance is mistaken for a derivation | tie gates to measured numerical error or a derived first-omitted order |

## 10. Recommended first continuation session

Do not begin by editing `invz_sigma_ordered` or widening Task 6b's gate. The first future session should:

1. Create a theory specification for a scalar two-level Hamiltonian with an explicit finite coupling matrix.
2. Enumerate the first-order vacuum and two-point diagrams and take the functional derivative on paper.
3. Map each local `C2/C3/C4` object to `invzt_vertex3/4` and identify any missing frequency/source indices.
4. Determine whether the current EMT closure has a generating functional; if not, initially use the direct lattice propagator without EMT.
5. Design a two- or four-site exact-diagonalization oracle and freeze exact `F`, `U`, `m`, and `chi` values at several coupling scales.
6. Only after the derivation is independently reviewed, write a task-by-task implementation plan for `invz_functional/`.

At the time of writing, the repository is a shared, dirty worktree with unrelated user changes. A future implementation should start on a dedicated branch/worktree, inventory `git status`, preserve all existing edits, and avoid modifying the Jensen production path until the isolated backend passes its own references.

Routine regression commands from the repository root are:

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"
```

Do not enable the whole tensor `INVZ_SLOW` suite casually: it includes a roughly 7.4-hour `a3d` production point. New functional tests should live in their own test folder and keep production-scale gates separately selectable.

## 11. Definition of done

The extension is scientifically complete only when all of the following hold:

- its core equations contain no mandatory two-level moment parametrization;
- every self-energy/field vertex is the derivative of a documented common functional;
- `F`, `U`, `m`, and `chi` are produced by that functional and satisfy the corresponding identities within numerical error;
- strict-order coefficients match exact finite-cluster coefficients through the retained order;
- the first omitted error scales as derived, rather than as fitted after seeing one fixture;
- local-rank, Matsubara, lattice, and nonlinear-solver uncertainties are reported separately;
- the PM and ordered solutions meet at the same continuous boundary without a boundary-specific patch;
- Jensen-static results can be recovered as an optional comparison, but are not needed by the new backend;
- external HoF3 and LiHoF4 benchmarks are reproduced within declared physical and numerical scopes;
- documentation distinguishes thermodynamic consistency of the truncation from accuracy relative to the exact infinite lattice.

At that point, the current 13.65% issue becomes a useful historical benchmark: the new functional route should either close it to numerical precision within the truncated theory or identify, diagram by diagram, the same-order content omitted by Jensen's static ordered construction.

## 12. Side tasks carried from the 2026-07-22 Stage-2 closeout

Two concrete follow-ups were recorded when Stage 2 ([docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md](docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md)) completed and merged to `main` at `0edc0ab`. Both are smaller and more immediate than the WP0–WP6 programme above and can proceed independently of it; neither is a prerequisite for the extension project, but 12.1 shares its external-benchmark gate (§7 item 7; WP6).

### 12.1 HoF3 published-benchmark discriminator — scientific clearance for J 2.34

**Why.** The closed-model J 2.34 comparison leaves the stable ~13.65% two-route mismatch that motivates this whole note. It has been *adjudicated* as a probable same-retained-order static-elastic approximation residual — implementation not shown defective — but it is not scientifically *cleared*. The most informative single discriminator is a faithful reproduction of Jensen's own full-machinery low-temperature HoF3 Eq. 2.34 check, whose published result is 2–3%.

**Task.** Using the existing Jensen-static ordered machinery (`invz_gstat_ordered`, `invz_emt_static_ordered`, `invz_hmf_ordered`, `invz_deltaF_ordered`), reproduce Jensen's HoF3 physical lattice and parameters and run the low-temperature field-route-vs-internal-energy-route Eq. 2.34 comparison, budgeting numerical quadrature error separately from the physical residual.

- Reproduces 2–3% while the synthetic zero-mean fixture stays near 13–14% → the synthetic residual is convincingly fixture-dependent static-approximation error, and the current implementation is cleared for its declared scope.
- Faithful HoF3 reproduction also lands near 13–14% → a same-order implementation defect, or a material mismatch with Jensen's procedure, remains open, and the §7b diagnostic stance must hold.

Until this runs, `jensen`/`jensen_static` J 2.34 stays a labelled diagnostic, never a percentage-accuracy guarantee (§1).

### 12.2 Independent RPA evaluator on 1 − J·χ0(0) — Stage-1 diagnosis regression 4

**Why.** Stage 1 ([docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md](docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md)) did not discharge "regression 4": a critical-field evaluator built independently of the current solve path, from the bare RPA condition `1 − J·χ0(0) = 0`. A dead phase-0 `chirpa` computation remains in place pending it.

**Task.** Implement an independent RPA dispatcher evaluating `1 − J·χ0(0)` (bare — no ordered/Σ dressing) as a cross-check on the critical field, then remove the dead phase-0 `chirpa` computation it supersedes. Self-contained and unrelated to the J 2.34 static-elastic question; recorded here only because both were the open verification items at Stage-2 closeout.
