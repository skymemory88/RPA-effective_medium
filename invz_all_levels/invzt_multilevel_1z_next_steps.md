---
schema: invzt_multilevel_1z_plan/v1
plan_id: invzt-ml1z-2026-08-02
status: ready_for_ml0
module_path: invz_all_levels
host_worktree: invZ expansion
host_git_branch: invzp-exec-convg-plan
retired_git_branch: invzt-multilevel-1z
baseline_commit: 425143d8633a32457b21c1aca0024b821064adf3
created: 2026-08-02
active_theory_mode: multilevel_1z
retired_theory_modes: [a1, a2, a3, a3d]
protected_modules: [invz_projected]
shared_modules: [invz_common, invz_tensor]
shared_module_policy: explicit_promotions_only
primary_temperature_K: 0.1
transfer_temperatures_K: [0.3, 1.0]
diagnostic_fields_T: [0.0, 0.09, 0.36, 0.99, 1.5, 3.0, 4.6, 4.68, 4.9, 6.0]
production_field_grid_T:
  start: 0.0
  stop: 9.0
  step: 0.09
local_rank_nominal_ladder: [32, 48, 64, 96, 136]
local_rank_policy: promote_each_nominal_cut_to_complete_multiplets
frequency_cutoff_meV_ladder: [10, 20, 40]
lattice_ngrid_ladder: [8, 12, 16, 20]
budgets:
  production_sweep_wall_clock_hours: 24
  forecast_safety_factor: 2.0
  peak_rss_bytes_per_worker: 4294967296
  local_cache_bytes_per_worker: 1073741824
  retained_sweep_bytes: 536870912
  production_workers: 12
tolerances:
  two_level_oracle_relative: 1.0e-11
  static_c2_derivative_relative: 1.0e-7
  static_c3_derivative_relative: 1.0e-4
  permutation_relative: 1.0e-9
  kms_reality_relative: 1.0e-10
  backend_crosscheck_relative: 1.0e-10
  local_rank_observable_relative: 5.0e-4
  fixed_point_residual_relative: 1.0e-8
  returned_costate_relative: 1.0e-8
  lattice_grid_observable_relative: 1.0e-3
  matsubara_observable_relative: 1.0e-3
  unlisted_dimensionless_observable_relative: 1.0e-3
  equal_time_sum_rule_relative: 5.0e-3
observable_tolerances:
  static_susceptibility_relative: 1.0e-3
  pm_mass_absolute: 1.0e-4
  peak_energy_GHz_absolute: 5.0e-3
  integrated_spectral_weight_relative: 5.0e-3
  phase_boundary_field_T_absolute: 1.0e-2
  phase_boundary_temperature_K_absolute: 1.0e-2
gates:
  documentation_reset: complete
  ml0_local_oracle: pending
  ml0_5_cost_architecture: blocked_by_ml0
  ml1_theory_contract: blocked_by_ml0
  ml2_local_producer: blocked_by_ml0_5_and_ml1
  ml3_paramagnetic_map: blocked_by_ml2
  ml4_paramagnetic_certification: blocked_by_ml3
  ml5_ordered_theory: blocked_by_ml4
  ml6_branch_selection: blocked_by_ml5
  ml7_real_axis: blocked_by_ml6
  ml8_production_integration: blocked_by_ml7
production_authorized: false
---

# Full-tensor multilevel 1/z implementation plan

- Journal: [`docs/execution/invzt_multilevel_1z_journal.md`](docs/execution/invzt_multilevel_1z_journal.md)
- Theory evidence: [`docs/reference/invzt_multilevel_theory_sources.md`](docs/reference/invzt_multilevel_theory_sources.md)
- Ewald baseline: [`docs/execution/invzt_ewald_certification.md`](docs/execution/invzt_ewald_certification.md)

## 1. Objective and completion claim

Implement and certify a finite-level Jensen-style \(1/z\) calculation for the
full-tensor LiHoF4 response without compressing the local ion to one two-level
transition. The local connected cumulants must be evaluated from the retained
electronuclear spectrum, contracted with the effective medium as an additive
vertex correction, and propagated through the certified full-tensor Ewald
lattice response.

The program is complete only when the 101-field production spectra sweep at
\(T=0.1\) K is finite, phase-consistent, rank/cutoff converged within the frozen
tolerances, reproducible in serial and parallel, and within the frozen resource
budgets. Until then, output is diagnostic and `production_authorized` remains
false.

This plan does **not** promise that replacing the two-level closure will by
itself resolve every old nonconvergence. It tests that hypothesis while
separating three possible outcomes:

1. a representation defect is removed and the multilevel map converges;
2. a valid map still has multiple or unstable states, requiring ordered
   thermodynamics and branch selection; or
3. the exact retained-level contraction is computationally infeasible and the
   route stops before production refactoring.

## 2. Authority and evidence policy

Implementation authority, in descending order, is:

1. the equations and invariants accepted in ML1 of this document;
2. branch-local tests and measurements recorded in the journal;
3. the two-level reduction in `jensen_1z_framework.html`;
4. current structural references in `invz_tensor/`;
5. historical Git objects listed in the theory evidence map.

Historical code is evidence to audit, not code to trust by age. A copied
function remains diagnostic until its identities are reproduced on this branch.
If a current test contradicts a retained result, record the contradiction and
stop the dependent package rather than adjusting a tolerance after seeing it.

For arrays, every declared relative comparison uses

\[
\operatorname{relerr}(A,B;S)=
\frac{\max|A-B|}{\max\{\max|A|,\max|B|,S\}},
\]

with the nonzero scale \(S\) declared by the owning test before measurement.
For vertex and co-state residuals, \(S=4096\,\epsilon\max|G_0|\) on the active
subspace. Listed production observables use the per-observable YAML tolerances;
the generic unlisted tolerance may be used only for a dimensionless observable
whose definition and scale are preregistered in the journal. Near a zero or
phase boundary, use the listed absolute tolerance rather than a relative ratio.

## 3. Staged physics architecture

### 3.1 ML1Z-L: first production candidate

The first candidate evaluates exact retained-level connected cumulants of the
longitudinal local operator \(J_z\), produces an additive longitudinal vertex
\(\mathcal V_{zz}(i\omega_n)\), and embeds that vertex in the full \(3\times3\)
local response and full tensor Ewald effective medium. The medium contraction is
through \(K_{zz}\); transverse and mixed local vertices are not renormalized at
this stage.

This is a **longitudinal-vertex/full-medium approximation**, not exact
component-complete tensor \(1/z\). It is the simplest direct multilevel
generalization of Jensen's fluctuation channel and avoids the known dense
component-labelled A3 cost. Its omitted mixed-vertex error must be measured by
the ML0.5/ML4 discriminants.

### 3.2 ML1Z-T: conditional component-complete promotion

If ML1Z-L fails a mixed-channel control or changes production observables beyond
their listed tolerances, generalize the local oracle to ordered operator lists
and evaluate
\(C^{(4),c}_{\mu\nu\rho\sigma}\). This route is authorized only if its
identities and full-sweep cost forecast pass before production wiring.

The current `invzt_vertex4` dense path sum is a small-system reference for this
gate, not the production algorithm. Its disabled factored path remains disabled
unless a new exact equivalence proof and numerical cross-check pass.

### 3.3 Ordered and real-axis layers are downstream

PM fixed points are implemented first. Ordered work starts only after the local
and PM maps are certified. Real-axis spectra start only after the state selector
is certified, so the spectra driver cannot evaluate a different phase or theory
mode from the selected Matsubara state.

## 4. Candidate mathematical contract to freeze in ML1

These equations define the object to be checked, not permission to skip their
derivation.

Use bosonic labels \(n\in\mathbb Z\),
\(\omega_n=2\pi n/\beta\), and the repository Fourier convention

\[
X(\tau)=\beta^{-1}\sum_n X_n e^{-i\omega_n\tau}.
\]

For a source-biased local Hamiltonian and centered local operators, define the
connected four-point cumulant with external labels
\((n,-n,\ell,-\ell)\). All imaginary-time orderings, KMS endpoint layers,
repeated-node Hermite limits, and allowed disconnected pairings are retained.

The candidate first-order medium contraction is

\[
\mathcal V_{\mu\nu}(n)
=\frac{1}{2\beta}\sum_{\ell}\sum_{\rho\sigma}
K_{\rho\sigma}(\ell)
C^{(4),c}_{\mu\nu\rho\sigma}(n,-n,\ell,-\ell).
\]

For ML1Z-L, only \(\mu=\nu=\rho=\sigma=z\) is retained. No calculation may
form \(\mathcal V/G_0\). On the active local-response subspace, the candidate
symmetric resummation is

\[
\widetilde G_0
=G_0\,[G_0+\mathcal V]^{-1}G_0,
\qquad \widetilde\chi_0=-\widetilde G_0,
\]

implemented by solves, never literal matrix inversion. ML1 must confirm the
sign, the \(1/(2\beta)\) normalization, the active-subspace rule, and the exact
two-level reduction to Jensen's \(G_0/(1+\Sigma)\).

The full tensor effective medium continues to satisfy

\[
\overline\chi
=\langle(\widetilde\chi_0^{-1}-J(\mathbf q))^{-1}\rangle_{\mathbf q},
\qquad
K=\widetilde\chi_0^{-1}-\overline\chi^{-1},
\]

through the active-subspace solve contract already used by
`invz_tensor/invzt_emt_matrix.m`. The Ewald lattice object and its exact
provenance are immutable inputs to this map.

## 5. Non-negotiable invariants

1. **Module isolation.** New multilevel code belongs under `invz_all_levels`.
   Do not edit `invz_projected`; promote a shared change in `invz_common` or
   `invz_tensor` only deliberately, with its own validation.
2. **Ewald identity.** Ewald remains the default. Every point stores backend,
   Ewald parameters, grid, boundary convention, and cache identity. Brute force
   is never an implicit fallback.
3. **Additive vertex.** Never divide by a multipole propagator and never hide a
   zero with a floor.
4. **Exact local definition.** Cumulants keep all orderings, KMS terms,
   disconnected subtractions, and repeated-node limits.
5. **Controlled rank.** A nominal rank is promoted to a complete multiplet and
   the realized rank is returned. Thermal discarded weight is not a virtual
   truncation bound; observable rank ladders decide acceptance.
6. **Immutable map.** Context construction, one map evaluation, iteration, and
   state selection are separate. A returned point identifies both the input
   co-state used to evaluate \(\chi,K\) and the mapped output vertex.
7. **Residual consistency.** A converged return satisfies both vertex
   self-generation and response dressing at the reported co-state. A capped
   one-step diagnostic is labeled `map_evaluation`, not `converged=false` state
   ambiguity.
8. **No silent options.** Requested and realized rank, cutoffs, backend, source,
   theory mode, phase, and retry path are returned and tested.
9. **Fail-closed resources.** Preflight peak temporary memory, byte-bounded
   cache residency, persisted payload, and sweep wall time. Defaults are not
   raised after a failed forecast without an explicit plan revision.
10. **Thermodynamic honesty.** Fixed-point convergence does not rank roots. No
    ordered point is called production-stable until a derived selector passes.
11. **Reproducible continuation.** Warm seeds and ordered-boundary retries are
    deterministic recovery tools, not admissibility evidence. Cold/warm,
    forward/reverse, and serial/parallel classifications must agree.
12. **No inherited scalar repairs.** The A1 `missingAreaApproximation`, scalar
    phase selector, phase-boundary rescaling, and dominant/rest split are not
    reused unless separately re-derived for the multilevel object.
13. **No bulk cache payloads.** Saved point/sweep artifacts exclude local C4
    tables and cache internals unless a diagnostic explicitly requests them.
14. **Temperature scope.** Certification is primary at 0.1 K. Checks at 0.3 K
    and 1.0 K are transferability diagnostics, not separate production claims.

## 6. Required point provenance

Every diagnostic point must expose at least:

- schema version, Git commit, MATLAB version, timestamp, and theory mode;
- temperature, applied field vector, longitudinal source/mean field, phase
  label, and selector status;
- local spectrum digest, operator digest, full/realized rank, multiplet-cut
  gap, and discarded thermal weight;
- external/internal Matsubara labels and energy cutoff;
- cumulant backend, backend status, cross-check status, and cache hit/miss;
- Ewald parameters, grid, boundary convention, and lattice digest;
- input vertex, mapped vertex, scaled/component residual maxima, evaluated
  \(\widetilde\chi_0\), \(K\), and local lattice response;
- iteration/retry/continuation history and cold/warm provenance;
- PM/ordered admissibility diagnostics, phase-selector value/status, and
  rejection reason;
- peak RSS, retained-cache bytes, elapsed context-build time, elapsed map time,
  and persisted-payload estimate.

Large vertices/cumulants may be referenced by a content-addressed diagnostic
artifact rather than embedded in each point.

## 7. Work packages

### ML0 — Recover and certify the isolated finite-level cumulant oracle

**Purpose.** Establish a trusted local C2/C3/C4 producer before changing a
lattice or production caller.

**Tasks.**

1. Recover and rename only the useful local machinery from
   `09f0d3a:invz_functional/invzf_multilevel_cumulant.m` into tensor-owned
   symbols, provisionally `invzt_ml_cumulant.m` and `invzt_ml_local.m`.
2. Preserve a source note linking every recovered block to its Git object. Do
   not import the abandoned functional assembly or its production dispatch.
3. Replace the historical unsafe span-only backend selection with explicit
   `dense`, `action`, and `auto` modes. `auto` must use a conservative domain,
   detect nonfinite/failed results, retry with action, and report the route. It
   must never silently return a NaN.
4. Implement byte/work preflight for block matrices, right-hand sides, label
   count, and cache. Use byte-based eviction.
5. Reject a rank cut through a degenerate/near-degenerate multiplet. Return both
   nominal and realized ranks.
6. Add an independent two-level state-path oracle in the test file, not by
   calling the recovered implementation through a second wrapper.
7. Test rank-2/3/4 zero and nonzero Matsubara fixtures, repeated/degenerate
   nodes, disconnected C4 subtraction, operator Hermiticity, permutation, KMS,
   and expected-reality projection.
8. Test static source identities: C2 against \(dm/dh\) and C3 against
   \(d^2m/dh^2\), with finite-difference step ladders and a plateau criterion.
9. Reproduce the historical 136-state static C4 fixture and dynamic rank ladder
   as branch-local measurements. If conventions differ, explain the mapping
   before judging the number.
10. Reproduce the known large-`beta*span`, nonzero-frequency dense failure and
    prove `auto` returns a finite action result with explicit fallback status.

**Acceptance.**

- Two-level C2/C3/C4 agrees with the independent oracle within
  `two_level_oracle_relative`.
- Static C2/C3 derivative identities pass their frozen tolerances on a stable
  step plateau.
- Dense/action overlap agrees within `backend_crosscheck_relative`.
- Permutation and KMS/reality gates pass their frozen tolerances.
- Full-rank output is finite; truncated output is never marked rank-certified
  without an observable ladder.
- No production solver or driver calls the new oracle.

**Stop.** Any unresolved normalization, KMS, repeated-node, or backend
contradiction sets `ml0_local_oracle: stop_theory` and blocks every later package.

**Artifact.** `docs/diagnostics/invzt_multilevel_ml0/` plus a journal entry.

### ML0.5 — Measure cost and choose the staged architecture

**Purpose.** Decide feasibility before building a new self-consistent map.

**Tasks.**

1. Write a disposable but versioned cold/warm benchmark for the actual
   \((n,-n,\ell,-\ell)\) label set at ranks promoted around 32, 48, 64, 96,
   and 136; use Ecut 10, 20, and 40 meV.
2. Record separately: local diagonalization, per-label C4, label-table build,
   map contraction, cache lookup, peak RSS, and retained bytes.
3. Fit and report the observed complexity with uncertainty. Do not assume the
   dense A3 \(N^4\) constant applies to the block-exponential oracle.
4. Count distinct local states generated per field by source/ordered iterations.
   Forecast at least the PM 101-field sweep and a conservative ordered inner
   solve; apply the frozen factor-two safety margin.
5. Benchmark symmetry reuse only after an identity test proves the reused label
   relation. Compare every optimized result with the uncompressed oracle at a
   small rank to `1e-11` relative.
6. Forecast ML1Z-L first. Separately forecast a minimal mixed-Cartesian
   discriminant and the component-complete ML1Z-T route.

**Decision.** Record exactly one:

- `go_longitudinal`: ML1Z-L fits all frozen budgets with the safety factor;
- `conditional_compression`: an exact, identity-tested batching/tiled/action
  design has a quantified upper forecast within budget;
- `stop_cost`: no exact route fits, so production map work stops;
- `needs_human_decision`: only a newly labeled approximation can fit.

ML1Z-T is not authorized merely because ML1Z-L passes. The wall-clock ceiling is
the primary feasibility gate; a low call count alone is not evidence.

**Artifact.** `docs/diagnostics/invzt_multilevel_cost/` with raw timings and a
compact forecast report.

### ML1 — Freeze the multilevel theory contract

**Purpose.** Resolve the equations, conventions, and approximation boundary
before wiring the oracle into EMT.

**Tasks.**

1. Re-derive the connected C2/C3/C4 spectral/block-exponential object from the
   local generating functional, including centering, disconnected subtractions,
   KMS endpoints, and repeated-node limits.
2. Derive the medium contraction's sign, \(1/(2\beta)\) factor, frequency
   assignment, transpose/conjugation relations, and source dependence.
3. Derive the symmetric-bracket resummation on a rank-deficient matrix response
   using an active subspace and solves.
4. Prove the exact two-state reduction to Jensen's scalar equation through the
   retained \(1/z\) order. Use an independent symbolic check (Wolfram is
   appropriate) and numerical fixtures; retain both artifacts if used.
5. Define ML1Z-L precisely and list every omitted mixed/transverse vertex class.
   Design one cheap discriminating observable at 3 T and 6 T plus one
   near-boundary field.
6. Define the source convention. Local derivatives vary only the conjugate
   longitudinal source; transverse mean-field feedback is held fixed unless a
   separate total-derivative object is explicitly requested.
7. State which quantities are strict first order and which are resummations
   differing at \(O(1/z^2)\). Preserve the additive-versus-resummed spread as a
   method diagnostic.
8. Derive the equal-time sum rule and high-frequency moments used by ML4/ML7.

**Acceptance.** A compact theory note maps each equation to the local producer,
map producer, residual, and test. No sign or normalization is left to fixture
tuning. ML1Z-L is visibly labeled as such in every schema.

**Stop.** Failure of the two-level reduction or an unresolved source derivative
sets `ml1_theory_contract: stop_theory`.

### ML2 — Build the immutable local multilevel producer

**Purpose.** Turn the certified oracle into a reusable, bounded context for one
local single-ion state.

**Tasks.**

1. Build the source-biased 136-state single-ion state from `invz_single_ion` and
   expose digests of energy, operator, beta, source convention, and rank.
2. Produce C2, required C3, and only the C4 label table authorized by ML1.
3. Promote nominal rank cuts to complete multiplets and run the observable rank
   ladder. Never infer virtual convergence from thermal weight.
4. Use a byte-bounded cache keyed by the complete local-state and label-table
   signature. Provide explicit cache reset/stats for tests and workers.
5. Separate cached compact tables from diagnostic payload. No point object owns
   the cache or a full C4 table by default.
6. Return backend provenance and all ML0 quality gates with each context.
7. If ML1Z-T was authorized, generalize to ordered operator lists only after
   single-operator regression remains bit-identical/tolerance-identical.

**Acceptance.** Context replay is deterministic; a one-bit option change that
affects physics changes the digest; rank/backend requests appear in output; cache
residency remains below the frozen byte cap.

### ML3 — Implement one immutable PM map and its solver

**Purpose.** Make the high-dimensional fixed point inspectable rather than
burying it inside retries or the production dispatcher.

**Provisional interfaces.**

- `ctx = invzt_ml_build_context(ion,T,B,lat,opts)`
- `ev = invzt_ml_map(ctx,Vin)`
- `pt = invzt_ml_solve_pm(ctx,seed,opts)`

Names may change once, in ML1, but the separation may not.

**Tasks.**

1. Embed ML1Z-L \(\mathcal V_{zz}\) in a \(3\times3\) additive vertex and apply
   the accepted matrix resummation on the active response subspace.
2. Evaluate full tensor EMT with the certified Ewald lattice. Contract the local
   C4 table with the evaluated \(K_{zz}\) to produce `Vmapped`.
3. Return `Vin`, `Vmapped`, evaluated `chi_til`, `K`, `chi_bar`, component-wise
   scaled residuals, active rank, and domain diagnostics from one map call.
4. Add damped iteration first. Anderson/Newton-Krylov is conditional on a
   branch-local Jacobian/residual check; it changes the solver, not the map.
5. On every exit, regenerate the reported response at the reported vertex. A
   one-step reevaluation remains explicitly a map tuple even though it is not a
   converged fixed point.
6. Define admissibility separately from convergence: finite matrices, local
   spectral positivity/reality, active-subspace conditioning, EMT domain, PM
   uniform mass, and rank/cutoff status.
7. Add zero-field symmetry handling derived from the operator basis. Do not
   throw the old `invzt:a1ZeroField` path from multilevel mode.

**Acceptance.** Synthetic and small-ion maps pass co-state/residual contracts;
fixed points are seed-independent where the measured Jacobian implies local
uniqueness; all failures carry explicit status and provenance.

### ML4 — Certify the PM map and test the original nonconvergence hypothesis

**Purpose.** Determine whether finite-level additive physics removes the old
scalar-A1 failure before ordered complexity is introduced.

**Fixtures.** Primary \(T=0.1\) K fields are the frozen YAML list. The 4.68 T
neighborhood tests the expected continuous boundary; 3 T and 6 T are accepted
anchors; 0, 0.09, 0.36, and 0.99 T deliberately probe low-field PM instability
or rejection. Transfer checks use 0.3 K and 1.0 K.

**Tasks.**

1. Run cold and warm solves at all fixtures with rank, Ecut, and lattice ladders.
2. Record convergence and admissibility separately. A converged negative-mass PM
   point is an unstable PM solution, not a solver failure and not an ordered
   state.
3. Compare ML1Z-L with the exact two-level reduction fixture. Comparing against
   current A1 production at 136 states is diagnostic only because the physical
   closures differ.
4. Measure the mixed-channel discriminant at the ML1 fields. Promote to ML1Z-T
   only if the omission changes a production-relevant observable by more than
   its listed YAML tolerance and the ML0.5 full-route gate passes.
5. Check forward/reverse field order, serial/parallel, cache cold/warm, and
   deterministic retry classification.
6. Grade equal-time sum rule, high-frequency moments, rank, Matsubara, lattice,
   additive-versus-resummed spread, and local/map residuals.
7. Produce a PM-only diagnostic plot after the gates pass. Clearly mask fields
   where PM is physically inadmissible; this is the first visual-inspection
   checkpoint, not the production spectrum.

**Decision.** Record whether old failures were representation failures, physical
PM instabilities, residual solver failures, or unresolved. Do not state “A3
resolves convergence”; A3 is not this program.

**Stop.** If no rank-converged PM result fits the budget, stop production work.
If high-field PM anchors cannot be certified, do not begin ordered theory.

### ML5 — Derive and implement the ordered/source-biased multilevel theory

**Purpose.** Supply the physics that retries could not invent in scalar A1.

**Tasks.**

1. Derive the ordered local generator with scalar longitudinal source \(h\),
   its moment, C2, C3, and C4. Include the three-point/tadpole and one-point terms
   required away from the symmetric state.
2. Derive the \(1/z\) correction to the equation of state and the correct field
   differential. Distinguish fixed-vertex partial derivatives from total
   derivatives through the self-consistent map.
3. Reduce the outer ordered problem to a declared scalar coordinate (preferably
   \(h\) or \(m\)): for each candidate coordinate, solve the inner multilevel
   vertex map to its stated branch and evaluate a scalar equation-of-state
   residual. Do not claim enumeration in the full Matsubara-vector space.
4. Derive a common thermodynamic functional whose stationarity reproduces both
   the inner map and equation of state, or prove and record that the chosen
   resummation lacks one.
5. If a common functional exists, verify its derivative numerically and use it
   to rank coexisting states. If it does not, label states
   `accepted_ordered_unranked`; they cannot fill the strict production class.
6. Verify the \(m\to0\) ordered limit against the PM map at a continuous
   boundary, including vertex, medium, susceptibility, mass, and functional
   difference.
7. Reconsider `missingAreaApproximation` only through a new derivation. The
   projected scalar implementation is not evidence that it applies here.

**Acceptance.** Static derivative identities, equation-of-state residual,
co-state consistency, boundary matching, rank/cutoff gates, and (if claimed)
functional-gradient/selector tests all pass.

**Likely limitation.** A common functional for the chosen partial resummation may
not exist. That is an acceptable scientific result; it limits output to PM and
unranked ordered branches rather than licensing a heuristic selector.

### ML6 — Enumerate scalar ordered roots and track branches

**Purpose.** Make “all relevant roots” operational and make continuation a
reproducible accelerator rather than a selector.

**Tasks.**

1. Freeze the physical outer interval in \(h\) or \(m\) from operator bounds
   and lattice/source conventions.
2. Establish whether the inner vertex map is locally unique on each outer
   interval using the measured Jacobian/spectral radius. If not, assign stable
   inner branch IDs and treat each as a separate scalar residual.
3. Adaptively subdivide every continuous scalar residual interval, identify
   poles/discontinuities, bracket every sign-changing root, and test tangencies
   with derivative/stationary-point searches.
4. Define “all relevant roots” as all isolated roots found in the frozen bounded
   outer interval on every certified inner branch at the declared search
   resolution and tangency tolerance. If inner branches cannot be bounded or
   enumerated, the completeness claim fails explicitly.
5. Add pseudo-arclength/continuation only after cold point solves pass. Run
   forward/reverse field sweeps and retain branch IDs across folds.
6. Add warm-seed retry and ordered-boundary retry in a fixed order. The final
   state must be identical in classification to cold enumeration; retry history
   is provenance, not evidence of lower free energy.
7. Require serial/parallel order independence.

**Acceptance.** Complete operational root tables, selector values/status,
boundary continuity, retry determinism, and forward/reverse agreement at all
frozen fixtures.

**Stop.** If ML5 has no selector, ML6 may publish branch tables but cannot
populate `accepted_multilevel_ordered_strict`.

### ML7 — Derive and certify the real-axis multilevel response

**Purpose.** Produce spectra from the same accepted state and theory, without
freezing a static correction or numerically continuing noisy Matsubara data.

**Tasks.**

1. Derive a direct Lehmann/complex-frequency representation of the required C2
   and contracted C4/vertex, including repeated poles and source-biased terms.
2. Use the same Ewald lattice provenance and accepted local/ordered state as the
   Matsubara point. The phase selector runs once upstream; spectra consume its
   immutable state ID.
3. Add a contract test that the spectra driver requested and realized
   `theory_mode='multilevel_1z'`, correct phase, source, rank, and selector ID.
   This directly guards the historical wrong-phase-selector failure class.
4. Verify upper-half-plane analyticity, conjugation, nonnegative physical
   spectral weight under the repository sign convention, Kramers--Kronig/static
   agreement, high-frequency moments, and equal-time sum rule.
5. Run broadening and frequency-grid ladders. Report unresolved poles rather
   than clipping or replacing them.

**Acceptance.** Matsubara-to-static, causality, sum-rule, grid/broadening, phase
identity, and rank/cutoff gates pass on PM and every strict ordered fixture.

### ML8 — Integrate and certify production

**Purpose.** Expose one explicit multilevel mode to the drivers and produce the
first production-quality visual output.

**Tasks.**

1. Add an explicit `theory='multilevel_1z'` dispatch. Do not overload `a1`,
   `a2`, `a3`, or infer the mode from fields in an options struct.
2. Keep legacy modes non-default and clearly labeled. Switch the production
   spectra default only after every prerequisite gate is complete and a clean
   session test passes.
3. Require certified Ewald input. An explicit brute diagnostic must carry a
   visibly different provenance and cannot seed a production result.
4. Run the 101-field \(0\)--9 T sweep at 0.09 T spacing, first serial on a small
   subset and then with 12 process workers. Compare overlapping point digests
   and classifications.
5. Strip C4 tables, block-exponential workspaces, and cache internals from saved
   points. Enforce the retained-sweep byte cap before saving.
6. Save a machine-readable manifest with every gate, rejection reason, resource
   measurement, branch/commit, and artifact digest.
7. Generate the spectrum/phase visual for user inspection with masked points
   distinguished as physical rejection, numerical failure, cost stop, or
   unranked ordered state.
8. Re-run from a fresh MATLAB session with caches cold and require the same
   accepted/rejected classifications and observables within tolerance.
9. Run transfer-temperature diagnostics at 0.3 and 1.0 K and record scope; do
   not broaden the production claim beyond 0.1 K without a plan revision.

**Production gate.** Set `production_authorized: true` only when every strict
field is selected, all numerical/physics/resource gates pass, the manifest is
complete, and the visual output corresponds exactly to the saved accepted
states. If ordered selection remains unavailable, the result is PM-only plus
unranked branches, not a completed production sweep.

## 8. Adversarial review checklist

Before completing each major gate, perform a separate review for:

- sign or \(\beta\)-normalization disagreement between derivation, oracle,
  contraction, residual, and tests;
- co-state/time-level mismatch between returned vertex, susceptibility, medium,
  moment, and selector;
- a requested rank/backend/cutoff silently ignored by a caller;
- a cache key that omits source, operator, spectrum, temperature, labels, rank,
  Ewald provenance, or theory schema;
- a degeneracy cut whose apparent convergence is basis-dependent;
- rank convergence in local C4 but not in the final observable;
- a retry or continuation order changing the selected phase;
- an exact Gamma/directional-Ewald convention lost between Matsubara and real
  axis;
- a phase selector evaluated from A1/A3 state rather than the multilevel state;
- a sum rule satisfied only after clipping, missing-area replacement, or a
  broadened tolerance;
- cumulative worker memory and saved payload omitted from single-call preflight;
- test fixtures that only exercise zero frequency, zero source, or a two-state
  model and therefore miss the production failure surface.

## 9. Fresh-session handoff

The next executable action is **ML0 task 1**: recover the schema-v2 local
cumulant oracle from Git history into tensor-owned names without wiring it into
production. Before editing, verify:

1. `invz_all_levels` is the active module and the host worktree's pre-existing
   status is recorded without modifying unrelated paths;
2. the shared `invz_projected` and `invz_tensor` paths are treated as external
   inputs unless an explicit shared promotion is authorized;
3. Ewald certification tests still pass;
4. the ML0 proof obligations and frozen tolerances above are copied into the
   focused test file.

Land ML0 as its own reversible commit only after the independent two-level,
static-derivative, backend, KMS, and nonfinite-fallback gates pass. Then update
the journal and YAML gate state before starting ML0.5.
