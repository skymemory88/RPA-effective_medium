# Projected-spin convergence execution journal

## Scope and recording policy

This journal tracks execution of `invzp_convergence_next_steps.md`. It records
checkpoints, decisive evidence, exact checks, assumptions, and unresolved risks.
Raw numerical logs and large diagnostic tables belong in referenced artifacts,
not here.

## Standing project constraints

These constraints apply across checkpoints, context compactions, and future
sessions:

1. The ultimate target is a certified projected-spin \(1/z\) susceptibility
   calculation across \(0<B_x\le 9\) T. The strict zero-field point may be
   excluded if its special symmetry/domain requires separate treatment.
2. Nonconvergence need not have a single cause. Classify failures by layer and
   evidence; do not force a universal repair across physically or
   algorithmically distinct regimes.
3. After five consecutive failed trial methods, stop proposing variants and
   return to the innermost governing mathematical, physical, and algorithmic
   assumptions. Record the five trials and their rejecting observables before
   any sixth method is considered.
4. Pause for human judgment whenever the next action requires an unvalidated
   physical choice, thermodynamic branch selection, a material scope change, or
   risk disproportionate to the available evidence. Steady certified progress
   is preferred to aggressive unsupervised advancement.

## 2026-07-30 — Checkpoint 0: baseline and first work packet

**Repository state.** Branch `invzp-exec-convg-plan`, starting commit `693c518`
(`chore(invzp): checkpoint reduced convergence baseline`), clean worktree.
MATLAB R2025a is available at `/Applications/MATLAB_R2025a.app/bin/matlab`.

**Governing sources.**

- `invzp_convergence_next_steps.md`: work package 1 is the immediate task.
- `invzp_convergence_dead_ends.md`: prior residual-only, strict-medium,
  principal-value, and unbounded root approaches remain rejected.
- Production ordered closure:
  `invz_projected/invz_solve_point_ordered.m` and
  `invz_common/invz_gstat_ordered.m`.
- Lattice construction:
  `invz_projected/invz_bz_couplings.m`,
  `invz_projected/invz_jq_modes.m`, and
  `invz_projected/invz_phase1_qgrid.m`.

**Proof obligations for checkpoint 1.**

1. Numerically and structurally audit
   \(J_{\rm sup}=\max[J_{0,\rm eff},\sup_{q,\nu}J_{q\nu}]\), retaining the
   uniform/directional-\(\Gamma\) contribution even though production drops
   exact \(\Gamma\) from the Brillouin-zone average.
2. Implement the exact production elastic closure as a bounded scalar problem
   in \(x=\widetilde G_0(0)\), with explicit pole, finite-value, \(\xi\), and
   uniform-mass gates.
3. Search the entire configured physical interval for every isolated root
   found under documented adaptive scan tolerances; never accept a sign change
   across a detected discontinuity.
4. Add focused MATLAB tests for algebra, domain rejection, discontinuity
   rejection, multiple-root enumeration, seed independence, and the roadmap's
   production fixtures.

**Initial code fact.** The live solver still uses a nested damped iteration on
\(K_0\) and accepts on closure residual alone. The historical one-shot
“strict-medium” branch was removed in `693c518` and is rejected by the current
ledger; it will not be restored.

**Current assumption.** The first implementation will preserve the production
full-electronuclear \(G_0^{\rm inel}/G_0^{\rm el}\) split and two-level
\(\xi\) exactly. Representation changes are deferred to work package 3.

**Lattice-supremum audit.** Using the frozen diagnostic configuration
(`16^3`, half-open, unshifted, `P_drop`, Ewald
`alpha=0.3,r_cut=16,g_cut=3,boundary=conducting_k0_omitted`):

- `info.Jcc0 = 0.0064216618094169392` meV;
- `max(Jnu_flat) = 0.0063717257361577892` meV, at mesh
  \(q=(0,0.0625,0)\);
- the mesh-to-uniform gap is \(4.993607325915004\times10^{-5}\) meV;
- a 10,001-point scan of the exact directional nonanalytic form over
  \(k_z^2\in[0,1]\) has its maximum at \(k_z^2=0\), equal to
  `info.Jcc0` within \(9\times10^{-19}\) meV.

The directional matrix is `info.Jpath_base_cc` plus a positive-semidefinite
all-sublattice rank-one term whose coefficient decreases with \(k_z^2\).
Therefore the in-plane endpoint is the directional top edge, not merely the
largest sampled direction. For this lattice provenance,
\(J_{\rm sup}=J_{0,\rm eff}\). The exact \(\Gamma\) point remains absent from
the Brillouin-zone average but is retained as the domain boundary and uniform
mass gate.

**Check command.** MATLAB batch construction through `invz_bz_couplings`,
followed by the directional form exported by `invz_jq_modes`; assertions
required `directional_max == Jcc0` and `mesh_max <= Jcc0` to \(5\times10^{-13}\)
meV. Both passed.

**Healthy-state regression fixture.** Before changing the legacy inner solve,
the full ordered solver was run at \(T=0.1\) K, \(B_x=4\) T with the same
couplings, `mix_outer=0.3`, `max_outer=500`, `tol_outer=1e-8`, and `Ecut=40`.
It converged with
`hmf=0.01577947695662214`, `K0=0.0007415288139413191`,
`Gstat=-118.0307001345370`, `Sigma0=0.01376361306413017`,
`D_uni=0.3295699256867276`, and outer residual \(5.27\times10^{-9}\).
The complete fixture and provenance are retained at
`docs/diagnostics/invzp_static_wp1/legacy_4T_fixture.mat`. This is a regression
target, not an admissibility certificate.

**Next action.** Define the bounded solver interface and diagnostics around
this verified `Jsup`, then add synthetic algebra/domain tests before running
the production fixtures.

## 2026-07-30 — Checkpoint 1: bounded physical static solver

**Implementation.**

- Added `invz_projected/invz_emt_static_ordered.m`, solving in normalized
  \(s=-J_{\rm sup}x\in(0,1)\). It exports a root table, explicit search
  resolution, discontinuity brackets, lattice provenance, and the status
  vocabulary `ok`, `no_admissible_static_root`,
  `multiple_admissible_static_roots`, or `static_search_unresolved`.
- Added `invz_common/invz_bounded_roots.m`. It enumerates all sampled
  sign-changing intervals and sampled minima of \(|\widehat R|\); a sign
  bracket is accepted only if a finite point reaches the residual tolerance.
  This prevents a pole sign change from masquerading as a root.
- Extracted the former local static solver from
  `invz_solve_point_ordered.m`. The compatibility `K0_seed` is ignored; more
  than one admissible root is reported but not arbitrarily selected.
- Made `invz_gstat_ordered.m` array-capable so the scan calls the production
  elastic closure itself rather than a replicated formula. Added exact
  `xi_numer`, `xi_denom`, `U`, and `V` diagnostics.
- Ordered outer loops now stop before feeding an uncertified static state into
  \(\lambda\) or \(\Sigma\). Fixed-\(H_{\rm MF}\) profiles retain each node's
  static status separately from `outer_iteration_failed`.

**Physical gates.** A candidate must have \(x<0\),
\(1+J_{\rm sup}x>0\), finite \(K_0,U,V,G_{\rm stat}\), nonnegative finite
\(\xi\) with positive \(\xi\) denominator, negative \(G_{\rm stat}\), closure
and root residuals at or below \(10^{-10}\), positive uniform mass, and
positive mesh denominators in both \(x\) and medium coordinates. No unproved
upper bound on \(\xi\) was introduced.

**Focused tests.**

- `test_invzp_static_domain`: 12/12 checks passed. It covers two-root
  enumeration, an even-multiplicity root, pole/discontinuity rejection,
  explicit rejection (but diagnostic retention) of a nonphysical-\(\xi\)
  mathematical root, \(m\to0\) equivalence with `invz_emt_scalar`, the
  production `J_sup` audit,
  bitwise preservation of the legacy scalar `invz_gstat_ordered` arithmetic,
  healthy 4 T reproduction, exact seed independence, stored 1 T node-9
  rejection, stored 3 T node-21--27 rejection, and positive accepted margins.
- `test_invzp_static_domain_resolution`: classifications are unchanged for
  `(scan_points,endpoint_margin)` equal to `(2049,1e-8)`, `(4097,1e-10)`,
  and `(8193,1e-12)`. The healthy-root \(x\) spread is
  \(6.94\times10^{-11}\ {\rm meV}^{-1}\).
- MATLAB Code Analyzer reports zero messages for both new solvers,
  `invz_gstat_ordered.m`, and both new tests. `git diff --check` passes.

Machine-readable results are
`docs/diagnostics/invzp_static_wp1/wp1_static_gate.mat` and
`wp1_static_resolution_gate.mat`.

**Acceptance results.**

- Healthy frozen 4 T state: one admissible root; relative to the pre-change
  fixture, \(\Delta K_0=-4.02\times10^{-15}\) meV and
  \(\Delta G_{\rm stat}=5.67\times10^{-11}\ {\rm meV}^{-1}\). The minimum
  signed accepted denominator/mass across the toy PM and production 4 T roots
  is \(0.30304635\).
- Stored 1 T node 9: unchanged state rejected and no alternative admissible
  static root found for that frozen outer state.
- Stored 3 T nodes 21--27: unchanged states rejected and no admissible static
  root found for any of those seven frozen outer states.
- \(m=0\) toy limit: agreement with the ordinary scalar medium is
  \(6.47\times10^{-12}\) in \(K_0\) and \(1.09\times10^{-13}\) in \(G\).

**Integration finding.** Re-running the complete legacy 4 T
\(H_{\rm MF}\)-integral construction no longer exports an ordered point. The
predictor and 28/33 profile nodes report `no_admissible_static_root`; only the
five largest-\(H_{\rm MF}\) nodes report `ok`. This is not a contradiction with
the healthy frozen endpoint: it confirms that the old field integral crossed
uncertified frozen outer states. The full result is retained at
`docs/diagnostics/invzp_static_wp1/new_4T_fullsolve_failed.mat`.

**Interpretation boundary.** A finite-resolution scan cannot prove the absence
of an arbitrarily narrow root. The scan size, maximum gap, endpoint margin,
root tolerance, and discontinuity brackets are therefore part of every
diagnostic. The three-level resolution ladder above is the current numerical
evidence, not an interval-arithmetic existence proof.

**Next action.** Complete an adversarial WP1 review. If no blocking defect is
found, update the roadmap checkpoint and begin work package 2 by constructing a
frozen-node outer-residual evaluator; do not try damping/iteration variants
before measuring the admissible outer domain and local Jacobian.

### Adversarial review disposition

The bounded frozen-inner solver passed the algebra, domain, resolution, and
fixture review. The review found a blocking **promotion** issue, not a false
inner root: the legacy Picard caller starts many profile nodes at an outer state
with no admissible static root, so the newly strict caller exits immediately.
Continuing through those states would require using a pseudo-root, NaN, or
substituted value and is forbidden. Therefore:

- work package 1 is accepted as a frozen-inner checkpoint;
- the complete ordered production solver is explicitly not promoted;
- the 4 T end-to-end regression is retained as the entry condition for work
  package 2; and
- no damping, iteration-budget, or warm-start trial is allowed until a
  deterministic admissible outer residual has been defined.

## 2026-07-30 — Checkpoint 2 entry: define the outer residual

At a fixed full self-energy vector \(\Sigma\), the dynamic
\(K(i\omega_n)\), \(n>0\), is algebraic through `invz_emt_scalar`. The static
candidate \(K_0(x)\) changes every
\(\lambda_p=\beta^{-1}\sum_n w_nK_ng_n^p\). Hence

\[
\lambda_p(K_0)=\lambda_{p,\mathrm{dyn}}
+ \frac{w_0g_0^p}{\beta}K_0
\]

is affine in \(K_0\). Treating the previous Picard iterate's \(\lambda\) as a
fixed hidden seed would make \(\mathcal F[\Sigma]\) non-deterministic. The next
implementation must instead evaluate this affine \(\lambda(K_0(x))\) inside
the bounded static residual. For each admissible static root it then constructs
the full \(K\), all three lambdas, and
`invz_sigma_ordered`, yielding \(\mathcal F[\Sigma]\). If the static solve has
zero or multiple admissible roots, the outer map is respectively undefined or
multivalued and must report that status rather than select a branch.

**Next proof obligations.**

1. The affine lambda construction must equal `invz_lambdas` to roundoff for
   arbitrary test vectors and every static candidate.
2. Fixed-lambda WP1 behavior must remain unchanged and all existing gates must
   continue to pass.
3. A unique-root outer evaluation must reproduce a healthy converged state's
   self-energy residual.
4. Former failures must be classified by static-domain status before any
   Jacobian or accelerator is applied.

### Deterministic outer-map implementation and first measurements

Added `invz_projected/invz_ordered_outer_map.m`. Given \(\Sigma\), it:

1. computes the dynamic medium algebraically;
2. forms the exact affine coefficients of all three lambdas;
3. passes \(\lambda_{1,2}(K_0)\) into the bounded static search;
4. exports no map for zero/multiple admissible static roots; and
5. for a unique root, reconstructs \(K,\lambda,\Sigma_{\rm map}\) and
   \(\mathcal R_\Sigma=\Sigma-\Sigma_{\rm map}\).

`test_invzp_outer_map` passes. At the retained healthy 4 T state:

- the deterministic residual is \(5.24748067\times10^{-9}\), versus the legacy
  final residual \(5.26544204\times10^{-9}\);
- \(|\Delta K_0|=1.27\times10^{-14}\) meV and the largest lambda change is
  \(1.02\times10^{-11}\);
- the affine-lambda identity agrees with a fresh `invz_lambdas` call to
  \(8.67\times10^{-19}\);
- the dynamic denominator minimum is \(0.460035161\), with no nonpositive
  dynamic denominator; and
- at the same node, the map is also uniquely defined at \(\Sigma=0\), whereas
  the 4 T \(h_z=0\) predictor is `no_admissible_static_root` at \(\Sigma=0\).

Added the diagnostic-only
`invz_projected/invz_outer_dominant_eigen.m`. A known diagonal linear map
returns its dominant eigenvalue 2. At the healthy 4 T state, central-difference
steps \(10^{-5},3\times10^{-6},10^{-6}\) give
\(\lambda_{\rm dom}=-0.00877249\), spread \(1.1\times10^{-8}\), with
eigen-residual about \(5.2\times10^{-6}\). The healthy endpoint is strongly
contractive; it does not diagnose the low-field-profile failures.

The zero-self-energy domain census
`docs/diagnostics/invzp_outer_wp2/invzp_outer_zero_census.m` evaluates the 34
retained profile fields at each of 1 T and 3 T:

- 1 T: the map is defined at 13/34 nodes (node 20 and nodes 22--33);
- 3 T: the map is defined at 6/34 nodes (nodes 28--33);
- every other zero-state evaluation reports
  `no_admissible_static_root`; and
- among defined zero-state maps, the residual is
  \(0.00271\)--\(0.393\) at 1 T and \(0.0207\)--\(0.0332\) at 3 T.

This census is only a map-domain slice. Undefined at \(\Sigma=0\) does not
prove no coupled root; defined with a nonzero residual does not prove a fixed
point. Results are in `docs/diagnostics/invzp_outer_wp2/`.

All new outer-map, Jacobian, test, and census files have zero MATLAB Code
Analyzer messages; `git diff --check` passes.

**Next action.** Define a bounded coupled-root search on the admissible outer
domain for a very small fixture set. Before choosing Newton, Anderson, or
continuation, determine whether admissible-domain perturbations exist around
each fixture and whether the local Jacobian is measurable. Do not extrapolate
the \(\Sigma=0\) census into a coupled-existence claim.

### Small-fixture Jacobians and existence probes

The \(\Sigma=0\) Jacobian measurements at four defined fixtures give:

- 1 T node 20: power iteration does not settle within 40 iterations
  (\(\lambda\) estimate \(-0.242\), eigen-residual 0.057); no scalar
  contractivity conclusion.
- 1 T node 22: \(\lambda_{\rm dom}=0.24833\), eigen-residual
  \(5.66\times10^{-9}\).
- 3 T node 28: \(\lambda_{\rm dom}=-0.01632\), eigen-residual
  \(3.44\times10^{-6}\).
- 3 T node 33: \(\lambda_{\rm dom}=0.02601\), eigen-residual
  \(3.39\times10^{-7}\).

Undamped admissible Picard was used only on the three fixtures with a converged
contractive zero-state Jacobian:

- 3 T node 28: admissible coupled root found in 5 iterations, residual
  \(1.47\times10^{-9}\), final \(\lambda_{\rm dom}=-0.01495\),
  \(D_{\rm uni}=0.12155\), dynamic minimum 0.34485.
- 3 T node 33: admissible coupled root found in 6 iterations, residual
  \(4.15\times10^{-10}\), final \(\lambda_{\rm dom}=0.02748\),
  \(D_{\rm uni}=0.79519\), dynamic minimum 0.80930.
- 1 T node 22: failed method 1. Residual \(0.393\to0.194\), then iteration 3
  leaves the admissible static domain.

For 1 T node 22, a 41-point line scan along the failed second update shows that
fractions through 0.525 are admissible and 0.55 is the first sampled
inadmissible fraction. A targeted mix of 0.5 was therefore tested:

- failed method 2. The residual sequence is
  \(0.393,0.261,0.203,0.201\), then iteration 5 leaves the domain.
- At the last admissible state,
  \(\lambda_{\rm dom}=1.35957\) with eigen-residual
  \(3.34\times10^{-9}\).

For a real \(\lambda>1\), scalar damping has local factor
\(1-a+a\lambda>1\) for every \(a>0\). No further damping values will be tried
on this branch. This new rejecting observable is recorded in
`invzp_convergence_dead_ends.md`. It does not prove that no different coupled
root exists.

Evidence:
`wp2_outer_boundary_jacobians.mat`,
`wp2_outer_contracting_probes.mat`, and
`wp2_node22_step_domain.mat` under
`docs/diagnostics/invzp_outer_wp2/`.

For the two contractive 3 T roots, a halfway-to-root start with mix 1 and a
zero start with mix 0.5 both converge to the cold/mix-1 result within the
\(10^{-8}\) outer tolerance (maximum component differences
\(7.45\times10^{-9}\) and \(6.17\times10^{-9}\)). Mix 0.5 takes 22--24
iterations instead of 5--6 but does not change the exported branch on these
fixtures.

**Trial counter for 1 T node 22:** 2 consecutive failed methods
(undamped Picard; evidence-selected mix 0.5). The next action must not be
another scalar damping value. A Newton/continuation formulation would be a
different algorithmic question and still requires explicit domain handling
and root-completeness limits.

## 2026-07-30 — Checkpoint 2 pause / handoff

**Verified conclusions.**

1. The frozen static solver now enforces the physical \(x\) domain and passes
   all WP1 fixture, limit, discontinuity, seed, and resolution gates.
2. The outer map is deterministic in \(\Sigma\); lambda is no longer a hidden
   prior-iterate state.
3. Two certified coupled roots exist on the tested 3 T high-\(h\) fixtures
   (nodes 28 and 33), and their cold/warm and mix-1/mix-0.5 branches agree
   within tolerance.
4. The tested 1 T node-22 trajectory encounters a genuine admissible-domain
   boundary and becomes locally noncontractive before reaching a root. Scalar
   damping is ruled out for that measured branch, but another coupled root is
   not ruled out.
5. The complete 0--9 T susceptibility objective remains unresolved. In
   particular, the low-\(H_{\rm MF}\) field-integral nodes and the broader
   1 T/3 T coupled-root census are not certified.

**Adversarial limits / unresolved risks.**

- Root enumeration is finite-resolution, with explicit scan/margin ladders but
  no interval-arithmetic proof against arbitrarily narrow roots.
- \(J_{\rm sup}=J_{0,\rm eff}\) is verified for the frozen standard Ewald
  production configuration. Field-dependent ODD or other coupling
  configurations require their own directional-supremum audit.
- The matrix-free dominant eigenvalue does not characterize every
  non-normal/subdominant direction. Node 20 explicitly failed to yield a
  settled single mode.
- No free-energy or thermodynamic selector exists for multiple admissible
  coupled roots.
- A domain-aware Newton/continuation method may find a node-22 root, hit a
  fold/boundary with no root in that component, or find a different component.
  Choosing and interpreting that search is the next judgment-sensitive step.

**Checks at pause.**

- `test_invzp_static_domain`: 12/12 passed.
- `test_invzp_static_domain_resolution`: passed.
- `test_invzp_outer_map`: passed.
- MATLAB Code Analyzer: zero messages on every new kernel, test, and diagnostic
  script. The modified legacy ordered driver retains its existing
  variable-growth advisories in the adaptive profile-extension block.
- `git diff --check`: passed.
- At this pause, the checkpoint changes were intentionally uncommitted pending
  review.

**Recommended next action for human review.** Review the bounded-static
acceptance contract and the node-22 evidence, then choose whether the next
bounded packet is (a) a domain-aware Newton--Krylov search at node 22,
(b) field continuation from the certified 3 T node-28 root toward nodes
27--21, or (c) a representation audit before further nonlinear solving. Do not
run another scalar damping trial.

## 2026-07-30 — Checkpoint 2 review and selected diagnostic

**Review outcome.** A visual rerun showed every ordered-state output masked,
including formerly converged points. This is consistent with the audited
control flow rather than evidence that every ordered solution disappeared:
`invz_hmf_ordered` treats the failed \(h_z=0\) predictor as fatal, sweeps the
profile from low to high \(H_{\rm MF}\), and its legacy per-node iteration
starts without a certified high-field coupled state. The deterministic outer
map and the certified high-field roots are not yet used by that production
driver. Therefore the blanket mask is currently classified as a solver/
continuation-path symptom; it is not a physical no-order conclusion.

**Authorized bounded experiment.** After locally checkpointing the accepted
WP1 and partial-WP2 work, run one diagnostic-only 4 T continuation from the
largest profile \(H_{\rm MF}\) downward. Use the deterministic domain-gated
outer map, seed only from the immediately preceding certified node, and stop
at the first undefined map, unresolved/noncontractive Jacobian, or failed
Picard solve. Do not alter the production driver, reconstruct the
\(H_{\rm MF}\) integral, choose a thermodynamic branch, or add another damping
trial in this packet.

**Pre-experiment validation.** The three focused MATLAB tests pass, the WP1
resolution classification is unchanged, Code Analyzer reports no messages on
the new kernels/tests/diagnostics, and `git diff --check` passes.

## 2026-07-30 — Checkpoint 3: coarse-grid 4 T reverse continuation

**Implementation boundary.** Starting from clean commit `1e238a3`, the
diagnostic `invzp_4t_reverse_continuation.m` traverses the retained 33-node
4 T profile from largest to smallest \(H_{\rm MF}\). Each node is seeded only
by the immediately higher certified node. A node must pass a unique bounded
static map, a resolved leading-mode screen with spectral radius below one,
undamped domain-gated Picard, outer residual \(10^{-8}\), positive uniform/
mesh/supremum masses, and positive dynamic denominators. The production
driver and \(H_{\rm MF}\) integral are untouched.

**Linear-diagnostic correction.** At node 32, real power iteration cycled and
did not settle although its Rayleigh estimates remained small. A matrix-free
nonsymmetric Arnoldi diagnostic was added and first tested on a known real
rotation map: it recovered a complex-pair spectral radius of 0.8 with maximum
active-mode eigen-residual \(1.5\times10^{-10}\). At node 32, finite-difference
steps \(10^{-5},3\times10^{-6},10^{-6}\) all resolve
\[
\lambda=0.0013664\mathbin{\pm}0.0038831i,\qquad
\rho=0.00411649 .
\]
Four additional returned modes are numerical null modes of order \(10^{-12}\);
they are excluded from relative-residual acceptance because a relative
eigen-residual is ill-conditioned at zero. This corrected a diagnostic false
stop, not a nonlinear solver failure.

**Continuation result.**

- Nodes 33 through 29 are certified in 4--6 undamped Picard iterations, with
  final residuals \(3.40\times10^{-10}\) to \(4.77\times10^{-9}\).
- These are exactly the same five high-\(H_{\rm MF}\) nodes that the strict
  legacy sweep already reached; the reverse sweep newly certifies zero coarse
  nodes.
- Across the five roots, the minimum uniform mass is 0.09770, minimum
  supremum mass 0.08162, minimum mesh-\(x\) mass 0.08877, minimum
  medium-coordinate mesh mass 0.10624, and minimum dynamic absolute
  denominator 0.34805. No nonpositive dynamic denominator occurs.
- Transferring the certified node-29 \(\Sigma\) directly to node 28 makes the
  bounded static map `no_admissible_static_root`, so the diagnostic stops as
  specified.

The node-28 frozen-seed classification is unchanged at
`(scan_points,endpoint_margin)` equal to `(2049,1e-8)`, `(4097,1e-10)`, and
`(8193,1e-12)`: zero mathematical roots, zero admissible roots, zero unresolved
minima, zero discontinuity brackets, and minimum sampled
\(|\widehat R|=3.2582\). Evidence is in
`wp2_4t_reverse_continuation.mat` and
`wp2_4t_reverse_boundary_audit.mat` under
`docs/diagnostics/invzp_outer_wp2/`.

**Conclusion and limit.** Reversing the unchanged coarse profile grid is not
a sufficient solver repair. The blanket production mask remains a
control-flow/predictor symptom, but this experiment also shows that a single
coarse transfer does not carry the certified component below node 29. It does
not establish whether the coupled branch ends between nodes 29 and 28, crosses
a fold/domain boundary, or can be followed with smaller field steps.

**Trial accounting.** The 4 T coarse reverse-continuation counter is 1 failed
method (failed to extend coverage). The 1 T node-22 counter remains 2 failed
methods and is not incremented by this different fixture. The rejected power
iteration was a corrected measurement defect, not a nonlinear trial.

**Checks at checkpoint.** All three focused MATLAB tests pass; the Arnoldi
test includes a complex pair with a four-dimensional nullspace and a separate
synthetic domain-boundary stop. Code Analyzer reports zero messages on every
new/modified kernel, test, and diagnostic script. `git diff --check` passes.

**Next bounded option.** If another unsupervised experiment is authorized,
bracket the 4 T boundary with adaptive downward field steps from node 29 while
keeping the same deterministic map and undamped acceptance contract. Stop at
the first reproducible loss of the coupled component; do not infer
thermodynamic equilibrium from continuation.

## 2026-07-30 — Checkpoint 4: adaptive 4 T node-29 to node-28 continuation

**Production status.** None of the deterministic-map or continuation
diagnostics is wired into `invz_solve_point_ordered`. The blanket ordered mask
seen in a production rerun is therefore expected: the fatal \(h_z=0\)
predictor and low-to-high legacy profile iteration remain unchanged. PM output
is unaffected.

**Initial adaptive result and contradiction.** A diagnostic controller was
started from the certified 4 T node-29 solution. It halves rejected
\(H_{\rm MF}\) proposals, grows accepted steps by at most 1.5, and otherwise
keeps the unique bounded static map, leading-mode gate, undamped Picard, outer
residual, and denominator acceptance contract fixed. The first run advanced
through nine accepted roots to \(h=0.00836206475\), then repeatedly reported
Arnoldi domain boundaries even for \(3\times10^{-9}\) self-energy
perturbations. This contradicted the positive explicit margins and a fully
defined first Picard segment, so nonlinear continuation was paused.

**Root mathematical cause.** The failing perturbation always contained the
same interior mathematical static root, with supremum mass 0.01255 and uniform
mass 0.01608. Some scan configurations exported it and others rejected it only
because:

1. sign-bracket bisection stopped at the first point satisfying
   \(|\widehat R|\le10^{-10}\);
2. the independently recomputed but algebraically equivalent closure residual
   has different local conditioning and landed at
   \(1.01\)--\(1.08\times10^{-10}\); and
3. both quantities were gated at the same \(10^{-10}\) tolerance.

Grid phase therefore changed the first below-tolerance bisection point and
produced false `no_admissible_static_root` statuses. The generic root
enumerator now continues a sign bracket to `x_tol` before applying residual
acceptance. It also does not launch a tangency minimization in a cell already
covered by an adjacent sign bracket, preventing a less accurately polished
duplicate of a simple root. A focused steep-root regression was added.

The retained failing direction then becomes uniquely admissible over the full
Cartesian product of 2049/4097/8193 scan points and endpoint margins
\(10^{-8}/10^{-10}/10^{-12}\), with root/closure residuals of order
\(10^{-12}\)--\(10^{-11}\). Arnoldi no longer encounters a map-domain
boundary. Very small finite-difference steps eventually become
roundoff-limited, as expected, but the physical classification is stable.

**Adaptive rerun.** Restarting from node 29 after the polishing fix reaches
node 28 in 11 proposals and six accepted downward steps. Direct node-29 to
node-28 transfers still have no admissible static root, so adaptive parameter
steps are genuinely required; polishing alone is not a substitute for
continuation.

At the node-28 target \(h=0.0080712748574\):

- all nine scan/endpoint combinations contain exactly one admissible root;
- maximum outer residual is \(4.07\times10^{-9}\), maximum root residual
  \(9.72\times10^{-12}\), and maximum closure residual
  \(2.82\times10^{-11}\);
- minimum supremum, uniform, mesh-\(x\), and medium-mesh masses are
  0.0010506, 0.0013877, 0.0088186, and 0.0116482;
- minimum dynamic absolute denominator is 0.31507, with zero nonpositive
  dynamic denominators;
- the dominant eigenvalue is \(-0.17726\), stable for finite-difference steps
  \(10^{-5},3\times10^{-6},10^{-6}\); and
- an undamped halfway-to-root start returns to the continued root within
  \(1.75\times10^{-9}\).

Evidence:
`wp2_4t_adaptive_boundary_continuation.mat` and
`wp2_4t_adaptive_target_audit.mat` under
`docs/diagnostics/invzp_outer_wp2/`.

**Interpretation.** This packet isolates two distinct algorithmic causes:
insufficient root polishing created false static-domain failures, while the
coarse \(H_{\rm MF}\) grid still crosses outside the continued component's
basin/domain. The result certifies one additional coarse node, not the complete
profile, the \(H_{\rm MF}\) integral, thermodynamic branch selection, or
0--9 T production susceptibility.

**Trial accounting.** The 4 T coarse reverse-grid method remains one failed
method. Adaptive continuation is a distinct successful method for the
node-29-to-node-28 gap. Rejected internal step proposals are controller
observations, not separate methods. The 1 T node-22 failed-method counter
remains 2.

**Checks at checkpoint.** `test_invzp_static_domain` passes 13/13 checks;
the 2049/4097/8193 static resolution ladder passes with healthy-root \(x\)
spread \(5.47\times10^{-12}\); `test_invzp_outer_map` passes; the target audit
is resolution-stable and branch-consistent; Code Analyzer reports zero
messages on all new/modified kernels, tests, and diagnostics; and
`git diff --check` passes.

**Next bounded option.** After validation and a clean checkpoint, extend the
same adaptive controller across only the next 4 T coarse gap (node 28 to
node 27). Reassess root count, masses, Jacobian structure, and step demand
before attempting a whole-profile controller or production integration.

## 2026-07-30 — Checkpoint 5: 4 T component endpoint below node 28

**Bounded experiment.** The same corrected adaptive controller was started
from the certified node-28 root and targeted node 27. The controller logic was
factored only enough to select the adjacent certified-start artifact; no
acceptance rule, nonlinear map, damping, or production code changed.

The controller accepted eight progressively smaller downward steps but reached
only \(h=0.008042919679\), far above node 27 at
\(h=0.00650417379\). The minimum-step gate then stopped the run. Unlike the
node-29-to-node-28 packet, step demand collapsed systematically:

- the uniform mass decreased from \(7.87\times10^{-4}\) at the first accepted
  step to \(2.74\times10^{-6}\) at the last;
- the equivalent supremum mass decreased from \(5.95\times10^{-4}\) to
  \(2.06\times10^{-6}\);
- the medium-mesh margin remained above 0.0103 and the dynamic minimum above
  0.3147; and
- the last-root dominant eigenvalue is \(-0.20142\), stable over
  \(10^{-5},3\times10^{-6},10^{-6}\), so outer noncontractivity is not the
  stopping cause.

**Component-edge audit.** The last root is uniquely admissible over the full
3-by-3 scan-density/endpoint-margin grid. Its maximum outer, root, and closure
residuals are \(2.33\times10^{-9}\), \(8.73\times10^{-12}\), and
\(3.04\times10^{-11}\). A halfway-to-root undamped start returns within
\(6.76\times10^{-9}\).

Linear fits over the last 3, 4, and 5 accepted roots give:

\[
h_c(D_{\rm uni})=0.00804286313\text{--}0.00804286319,\qquad
h_c(D_{\rm sup})=0.00804286314\text{--}0.00804286321 ,
\]

with displayed \(R^2=1\) at artifact precision. The two masses are the same
uniform/supremum pole condition expressed in the medium and physical-\(x\)
coordinates, so their common zero is expected rather than independent
evidence. The last accepted root lies only \(5.65\times10^{-8}\) above this
local extrapolated edge.

The last retained below-edge proposal at \(h=0.008042766469\) has zero
mathematical and zero admissible static roots for all nine scan configurations;
its minimum sampled scalar residual is \(1.44\times10^{-3}\). This is no
longer a marginal residual-polishing classification.

Evidence:
`wp2_4t_adaptive_node28_to27.mat` and
`wp2_4t_adaptive_component_edge_audit.mat` under
`docs/diagnostics/invzp_outer_wp2/`.

**Interpretation boundary.** The continued high-\(h\) component reaches the
physical static-domain endpoint \(1+J_{\rm sup}x=0\) at finite positive
\(H_{\rm MF}\). It cannot be step-shrunk through that pole while retaining the
current physical acceptance contract. This does not prove that no disconnected
admissible coupled component exists below \(h_c\), and it does not determine
which component is thermodynamically selected.

The production \(h_z=0\) predictor failure is therefore not merely an
initialization inconvenience for this component: the certified high-\(h\)
branch itself does not connect to \(h=0\). A production controller that simply
continues this branch cannot supply the full \(0\)-to-\(H_{\rm MF}\) integral.

**Trial accounting.** Adaptive continuation successfully crossed node 29 to
node 28, but one application of the same method failed to reach node 27 and
instead certified the component edge. The node-28-to-node-27 failed-method
counter is 1. The 1 T node-22 counter remains 2. No step-size variant will be
tried across this pole.

**Checks at checkpoint.** The generalized controller reproduces the committed
node-29-to-node-28 path and target audit. The component-edge audit passes its
3-by-3 resolution gate and halfway-start check. All three focused MATLAB tests
pass (including 13/13 static-domain checks); Code Analyzer reports zero
messages on the modified kernels/tests and all adaptive diagnostic sources;
and `git diff --check` passes.

**Required review before a next solver method.** Pause at this checkpoint and
review the underlying ordered construction:

1. determine whether a disconnected admissible coupled component exists below
   \(h_c\), using a bounded root search that does not assume continuity from
   the high-\(h\) component;
2. establish a thermodynamic selector before joining or switching components;
   and
3. reconsider the \(H_{\rm MF}\)-integral construction if its derivation
   requires one physical component to extend continuously to \(h=0\).

Do not wire adaptive continuation into production yet.

## 2026-07-30 — Checkpoint 6: production \(H_{\rm MF}\)-node resolution at 1--3 T

**Bounded experiment.** The production ordered-point solver was run directly
at \(B_x=1,2,3\) T with only `nH` changed from 33 to 65 and 129. Temperature,
couplings, the geometric field range, outer mix 0.35, 200-iteration budget,
\(10^{-8}\) outer tolerance, 40 meV cutoff, and the bounded static acceptance
contract were held fixed. No production default or control flow changed.

All nine points remain nonconverged and masked. In every case the independent
\(h=0\) predictor is `no_admissible_static_root`; changing `nH` cannot affect
that separately evaluated node. Profile-node coverage is:

| \(B_x\) (T) | 33 nodes | 65 nodes | 129 nodes | first converged \(h\), 33/65/129 |
|---:|---:|---:|---:|---:|
| 1 | 10/33 | 17/65 | 31/129 | 0.00624162 / 0.00774547 / 0.00862825 |
| 2 | 8/33 | 15/65 | 30/129 | 0.00891140 / 0.00891140 / 0.00844323 |
| 3 | 6/33 | 11/65 | 22/129 | 0.0115294 / 0.0115294 / 0.0109237 |

At all nine settings, every failed profile node is in one contiguous low-\(h\)
block and every converged node is in the complementary high-\(h\) block. All
failed static classifications are `no_admissible_static_root`; there are no
multiple-root, unresolved-search, or outer-iteration-only statuses.

**Nested-grid adversarial check.** The 33-node grid is the odd-index subset of
the 65-node grid, and the 65-node grid is the odd-index subset of the 129-node
grid, to at most \(1.04\times10^{-17}\) in \(h\). At 2 and 3 T, every shared
node retains exactly the same convergence verdict. The 129-node first-success
value is lower by the factor 0.9474635, exactly one inserted geometric
half-step. This is finer localization of the same sampled transition, not
evidence that increased resolution restores a missing component.

At 1 T, increased resolution is actively history-dependent:

- \(h=0.00624162311\) succeeds on the 33-node path with
  \(\Sigma_0=-0.292673\), but fails on the 65-node path with
  \(\Sigma_0=-0.821288\);
- \(h=0.00774546581\) succeeds on the 65-node path with
  \(\Sigma_0=-0.176826\), but fails on the 129-node path with
  \(\Sigma_0=-0.861066\).

These are identical shared field nodes, not interpolation differences. The
bounded static solve itself ignores `K0_seed`, while production `eval_node`
threads its mutable `Sigma` carrier across every low-to-high node. A node can
update `Sigma` during early outer iterations and later exit with
`no_admissible_static_root`; that partially updated state is still returned
and becomes the next node's warm start. The shared-node mismatch is therefore
strong evidence of failed-node state contamination. A direct transactional
rollback comparison remains to be performed before calling this mechanism
uniquely causal.

**Conclusion.** Increasing `nH` is not a convergence repair. It leaves the
fatal predictor unchanged, yields only one-step boundary localization at 2
and 3 T, and worsens the observed 1 T boundary through path-history
dependence. This reinforces the user's solver-level diagnosis and identifies
failure-state commit semantics as a separate algorithmic defect from the
physical static-component endpoint already certified at 4 T.

**Trial accounting.** The production node-resolution ladder is one failed
profile-level method. It does not increment the earlier 1 T node-22
nonlinear-method counter because it did not test a new coupled-root solver.
No further `nH` values will be tried without new evidence.

**Evidence and checks.**
`docs/diagnostics/invzp_outer_wp2/wp2_hmf_node_resolution_census.mat` retains
the exact settings, nine summary rows, and per-node states. Its one-off
generator was removed after checkpoint closure and remains recoverable from
commit `02bc7d3`. The generator passed MATLAB Code Analyzer and
`git diff --check`. Artifact assertions pass for all nine rows, the contiguous
failure/success topology, exact nested-grid correspondence, invariant shared
nodes at 2/3 T, and the two shared-node contradictions at 1 T.

**Next bounded option for human review.** Before changing production, make the
diagnostic node transition transactional: on any rejected node, restore the
last certified `Sigma`/static carrier (and separately compare a fresh start),
then re-evaluate only the two 1 T shared-node contradictions. This would test
the contamination hypothesis without claiming that rollback can cross a
genuine physical component boundary or select a thermodynamic branch.

## 2026-07-30 — Checkpoint 7: two-endpoint equation-(45) trapezoid

**Cleanup before the packet.** The completed node-resolution generator was
removed after its compact artifact and documentation were retained. The
maintained static-domain and outer-map regression tests were preserved, as
were older diagnostic generators that remain the reproducible definitions of
referenced artifacts. Cleanup commit: `cbba65f`.

**Proposed approximation.** For every candidate upper endpoint \(h\), replace
the complete Section-9.3 integral by

\[
h_0^{(2)}(h)=\frac{h}{2}\,[r_0+r(h)],\qquad
F^{(2)}(h)=h_0^{(2)}(h)-J_{0,\rm eff}m(h).
\]

This linearly interpolates \(r\) across the entire integration window and
uses no interior integrand values. The census used \(B_x=0.5,1.0,\ldots,9.0\)
T, the production 33-node upper-endpoint profile, outer mix 0.35, 200 ordered
iterations, 40 meV cutoff, and the unchanged bounded static gates. PM
continuation was followed transactionally from 9 T downward with a
1000-iteration diagnostic budget.

Two lower-endpoint definitions were kept distinct:

1. the strict ordered \(h=0\) endpoint currently required by production; and
2. the exact \(m=0\) PM-limit identity
   \(r_0=1+\Sigma_0\),
   \(\widetilde G_0(0)=-\chi_0/(1+\Sigma_0)\), used only when the PM fixed
   point itself converged.

An unconverged but finite PM last iterate was retained as an explicitly
uncontrolled comparator and never accepted.

**Coverage result.**

- At every bare-ordered field from 0.5 through 4.5 T, the strict ordered lower
  endpoint remains `no_admissible_static_root`.
- At 0.5 T, zero of 33 upper endpoints is admissible.
- At 1 T, ten upper endpoints are admissible but the two-endpoint residual
  remains negative through the whole valid block, placing any root above the
  sampled profile.
- From 1.5 through 3 T, finite brackets appear only if an unconverged PM last
  iterate is substituted for \(r_0\).
- The PM-limit lower endpoint converges at 3.5, 4.0, and 4.5 T and produces
  linear-interpolated coarse roots:

| \(B_x\) (T) | \(h_*^{(2)}\) | \(r_0\) | \(\widetilde G_0(0)\) | lower uniform mass |
|---:|---:|---:|---:|---:|
| 3.5 | 0.02125179 | 1.115808 | -229.9934 | -0.476940 |
| 4.0 | 0.01639873 | 1.080825 | -192.2376 | -0.234485 |
| 4.5 | 0.00870614 | 1.067688 | -164.8072 | -0.058336 |

- From 5 through 9 T, the converged PM endpoint has positive uniform mass and
  is the accepted phase; no ordered equation-(45) integral is needed.

**Endpoint audit.** At 3.5 T the two upper bracket nodes have positive uniform
masses 0.4490 and 0.5748; at 4 T, 0.3174 and 0.4331; at 4.5 T, 0.1184 and
0.1866. Their static closures, \(r\), and \(\widetilde G_0\) are finite and
admissible. The PM-limit lower endpoints are also finite, but their uniform
masses in the table above are negative.

The finite trapezoid therefore connects an unstable PM endpoint to a distinct
positive-mass ordered component while skipping the intervening zero of
\(1+J_{0,\rm eff}\widetilde G_0\). At 4 T this is consistent with the
independently certified high-\(h\) component endpoint near
\(h_c=0.0080428632\); the coarse root lies above the endpoint, but its nominal
integral crosses the unavailable interval below it.

**Conclusion.** The two-endpoint rule does “spit out” finite numbers at
3.5--4.5 T, but not across the full ordered field range. More importantly,
those numbers are finite because the approximation jumps over the exact
static-domain obstruction rather than resolving it. They are useful
diagnostics of scale only, not certified Jensen roots or susceptibility input.
No production code or default was changed.

**Trial accounting.** Endpoint-only linearization is one failed
complete-profile method. Together with the immediately preceding `nH`
increase, the post-endpoint complete-profile workaround counter is 2
consecutive failed methods. Unconverged PM-last-iterate brackets are
observations within this method, not additional trials.

**Evidence and checks.**
`docs/diagnostics/invzp_outer_wp2/wp2_endpoint_trapezoid_census.mat` retains
the 18-field table, endpoint states, residual brackets, and frozen options.
Artifact assertions pass for the strict-lower failures, PM convergence
partition, upper-endpoint topology, bracket coverage, and endpoint mass signs.
The one-off generator passed MATLAB Code Analyzer and `git diff --check`; it
was removed after checkpoint closure and remains recoverable from commit
`540ebbe`.

**Interpretation boundary.** Do not promote this quadrature rule to
production. Interior smoothness is not the decisive issue: a selected
admissible component spanning the integration interval is absent. Further
progress requires either finding and thermodynamically selecting a low-\(h\)
coupled component or deriving how equation (45) should be matched across a
phase/component change.

## 2026-07-30 — Checkpoint 8: temporary production visual-inspection mode

**Authorized temporary wiring.** The two-endpoint approximation is available
through the explicit ordered-solver option
`hmf_integral_mode='endpoint_trapezoid_visual'`. The normal solver default
remains `full_profile`. `invz_run_spectra.m` temporarily enables the visual
mode and permits a finite unconverged PM last iterate as \(r(0)\), with a
1000-iteration PM endpoint budget.

The visual branch:

1. independently evaluates the PM lower endpoint
   \(r_0=1+\Sigma_0^{\rm PM}\);
2. uses \(h_0^{(2)}(h)=h[r_0+r(h)]/2\);
3. allows failed interior profile nodes but requires two adjacent converged
   upper nodes to bracket the coarse root;
4. applies the same static gates during direct root refinement and final
   fixed-\(h\) reconstruction; and
5. exports profile status `ok_endpoint_visual`, result metadata
   `hmf_integral_mode`, and `visual_only=true`.

The driver emits a runtime warning and labels the plot
“VISUAL endpoint-Eq.-45 approximation.” Removing the three `hmf_*` fields
from `solve_opts` restores the strict production path.

**Focused end-to-end result.** With independent per-field PM endpoints, the
temporary solver returns finite ordered states at 1, 2, 3, 3.5, 4, and 4.5 T:

| \(B_x\) (T) | visual \(h_{\rm MF}\) | PM endpoint source | final outer residual |
|---:|---:|---|---:|
| 1.0 | 0.0146384 | unconverged last iterate | \(4.62\times10^{-9}\) |
| 2.0 | 0.0231908 | unconverged last iterate | \(4.51\times10^{-9}\) |
| 3.0 | 0.0252772 | unconverged last iterate | \(5.39\times10^{-9}\) |
| 3.5 | 0.0214252 | converged PM limit | \(4.66\times10^{-9}\) |
| 4.0 | 0.0165488 | converged PM limit | \(6.34\times10^{-9}\) |
| 4.5 | 0.00876178 | converged PM limit | \(5.50\times10^{-9}\) |

Every reconstructed upper state has `static_status='ok'`. At 5 T the ordered
bare bracket is unavailable and the dispatcher accepts the stable PM state.
A focused spectra-map check at 1, 4, and 5 T returns phases
`[ordered, ordered, PM]`, finite susceptibility columns, and
`visual_only=true`.

**Critical caveat.** The 1--3 T roots differ from the preceding endpoint
census because their unconverged PM last iterates were reached by independent
per-field runs rather than high-to-low PM continuation. This is direct
evidence of path dependence, not improved certification. These spectra are
for visual morphology only and must not vote on phase boundaries, equilibrium
selection, or quantitative susceptibility.

**Default neutrality and checks.** With no visual option, a direct 4 T solve
still returns the strict `node_failed` mask,
`integral_mode='full_profile'`, and the original
`no_admissible_static_root` predictor status. The maintained static-domain
suite passes 13/13 checks, the 2049/4097/8193 static resolution ladder passes,
and the deterministic outer-map regression passes. Code Analyzer reports no
messages on all three modified files; `git diff --check` passes. Regenerated
test artifacts were restored after the checks.

## 2026-07-30 — Checkpoint 9: external root-problem feedback review

**Input and retained evaluation.** User-supplied feedback is retained
unchanged in `feedback_root_Search_problem.md`. The claim-by-claim assessment,
derivations, implementation correspondence, and prioritized follow-ups are in
`feedback_root_search_problem_evaluation.md`.

**Verified useful result.** On any differentiable certified component,

\[
F'(h)=r(h)\,[1+J_{0,\rm eff}\widetilde G_0(h)].
\]

Thus the uniform-mass zero is a stationary point of the Section-9.3 residual
on that component. Combined with the 4 T continuation result, this strongly
supports a physical-domain endpoint of the followed high-\(h\) component
rather than loss of Picard contractivity.

**Interpretation boundary.** This identity does not prove that every
low-\(h\) failure is one connected unstable branch, exclude disconnected
coupled roots, establish the sign of \(F'\) where no branch was certified, or
prove uniqueness of a nonzero equation-(45) root. The 1 T shared-node
contradictions remain contaminated by rejected-state carryover.

**Two material corrections to the feedback.**

1. A saturation anchor, if its asymptotic condition is derived, has exact tail
   integrand \(r-1\):
   \(H_0(h)=h-\int_h^\infty[r(h')-1]\,dh'\).
   Replacing this by \(\Sigma_0\) is not generally valid in the finite-moment
   ordered elastic sector.
2. The exact-\(\Gamma\) channel is not added to the code's lattice
   \(\Phi(x)\) average with any discrete weight. The Gamma-dropped mesh defines
   \(\Phi\); \(J_{0,\rm eff}\) separately defines the uniform stability
   condition. The mesh/uniform gap is about 0.778%, not 0.08%. Refining the
   lattice DOS may improve \(\Phi\), but cannot by itself justify removing the
   uniform-mode gate.

**Next actions held for human review.** The lowest-risk algorithmic packet is
pure/transactional node evaluation at the two contradictory 1 T shared nodes.
The highest-value physical proposal is a derivation plus cheap high-\(h\)
feasibility test of the corrected saturation anchor. A four-scalar residual
and bordered pseudo-arclength continuation are promising only after exact
elimination/equivalence is derived and reproduced on healthy fixtures.
Feedback review is not counted as another trial method; the consecutive
complete-profile workaround counter remains 2.

## 2026-07-30 — Checkpoint 10: filtered multi-node visual integral

**Plan update.** `invzp_convergence_next_steps.md` now records the agreed
implementation order: pure accepted-state node transitions; a corrected
saturation-anchor feasibility derivation using \(r-1\); exact reduction of a
simultaneous residual if possible; bordered pseudo-arclength continuation; and
then representation/lattice audits. None relaxes the strict acceptance
contract.

**Rejected literal anchor.** The first implementation removed rejected or
nonfinite 65-node profile values and assigned \(H_0=0\) at the first retained
positive-\(h\) node. It returned `filtered_no_bracket` at 1, 2, 3, 3.5, 4,
and 4.5 T. The first retained node lay at
\(h=0.00764\)--\(0.01153\) meV and every retained residual remained negative.
At 4 T, permitting finite values from rejected nodes added no samples because
their \(r\) values were nonfinite. This dead anchor was removed from
production and recorded in `invzp_convergence_dead_ends.md`.

**Active visual construction.** The retained temporary mode is
`hmf_integral_mode='filtered_profile_visual'`. The driver sets `nH=65`,
evaluates the same independent PM endpoint used by the two-node experiment,
removes nonfinite positive-\(h\) nodes, and integrates the retained sequence.
A finite uncertified last iterate is permitted to enter only behind
`hmf_filtered_include_unconverged=true`; its count is exported. Root
refinement and the final susceptibility state must still converge. The
profile records the used-node mask, PM provenance, lower bridge width,
interior bridge count and widths, root-bracket indices, and whether the
bracket itself spans omitted nodes.

**Coverage census.** With independent per-field PM endpoints:

| \(B_x\) (T) | visual \(h_{\rm MF}\) | used positive-\(h\) nodes | PM endpoint |
|---:|---:|---:|---|
| 0.5 | no bracket | 2 | unconverged |
| 1.0 | 0.025673 | 17 | unconverged |
| 1.5 | 0.033570 | 19 | unconverged |
| 2.0 | 0.029772 | 15 | unconverged |
| 2.5 | 0.028069 | 13 | unconverged |
| 3.0 | 0.025625 | 11 | unconverged |
| 3.5 | 0.021095 | 10 | converged |
| 4.0 | 0.016021 | 9 | converged |
| 4.5 | 0.0086876 | 8 | converged |

All positive-\(h\) nodes used at these sampled fields were certified:
`filtered_unconverged_used_count=0`. Each retained block was contiguous, each
root bracket used adjacent retained nodes, and there were no interior bridge
panels. For the successful 1--4.5 T roots, the sole uncontrolled panel
connects the PM \(h=0\) endpoint to the first ordered node, over widths
\(0.00607\)--\(0.01153\) meV. Final accepted outer residuals are
\(4.56\times10^{-9}\)--\(6.23\times10^{-9}\).
At 0.5 T, the independent PM last iterate gives the pathological
\(r_0=-13.908\) and no bracket.

**Resolution behavior.**

| \(B_x\) | `nH=33` | `nH=65` | `nH=129` |
|---:|---:|---:|---:|
| 1 T | 0.028252 | 0.025673 | 0.024213 |
| 4 T | 0.016026 | 0.016021 | 0.015978 |

The 4 T visual root is resolution-stable on this ladder. The material 1 T
drift, accompanied by an upward-moving first retained node, is consistent with
the known failed-state continuation contamination and prevents a low-field
resolution claim.

**End-to-end and default checks.** A focused spectra map at 1, 4, and 5 T
returns phases `[ordered, ordered, PM]`, finite susceptibility columns,
`hmf_integral_mode='filtered_profile_visual'`, and `visual_only=true`. A
strict 4 T solve still returns `node_failed`, a NaN molecular field,
`full_profile`, `visual_only=false`, and
`predictor_static_status='no_admissible_static_root'`. The accepted 5 T PM
dispatch followed a warning from its separate bare ordered probe that mean
field had not met its tolerance after 800 iterations
(\(|dmf|=1.28\times10^{-7}\) meV).

**Checks and retained evidence.** The pure filtering/bridge helper has a
focused synthetic regression. MATLAB Code Analyzer reports zero messages on
all five touched MATLAB files. The bounded-static suite passes 13/13 checks,
the 2049/4097/8193 resolution suite passes, the deterministic outer-map test
passes, and `git diff --check` passes. Regenerated maintained-test MAT files
were restored. The compact census, resolution ladder, spectra-map result, and
strict-default comparator are in
`docs/diagnostics/invzp_outer_wp2/wp2_filtered_profile_visual_census.mat`.

**Interpretation and trial accounting.** This mode is suitable only for the
requested visual morphology check. It has more information than the two-node
rule, but its lower panel still crosses an unresolved stability/component
interval and its 1--3 T anchor remains iteration-path dependent. It is the
third consecutive complete-profile workaround that fails the certification
objective; the project-wide five-failure review trigger has not yet fired.

## 2026-07-30 — Checkpoint 11: 128-node visual production escalation

**Bounded change.** Following encouraging visual inspection of the 65-node
filtered-profile mode, `invz_run_spectra.m` now sets `nH=128`. No solver
equation, filtering rule, PM-anchor policy, acceptance gate, or strict default
changed. The plot label explicitly reports 128 nodes.

**Exact focused checks.**

| \(B_x\) | \(h_{\rm MF}\) | used positive-\(h\) nodes | first retained \(h\) | PM endpoint | final residual |
|---:|---:|---:|---:|---|---:|
| 1 T | 0.024394 | 31 | 0.008519 | unconverged | \(6.01\times10^{-9}\) |
| 4 T | 0.015982 | 18 | 0.009421 | converged | \(4.86\times10^{-9}\) |

Both final ordered states converge, both use zero uncertified positive-\(h\)
quadrature nodes, and both agree closely with the prior `nH=129` probes
(0.024213 and 0.015978). This verifies the driver setting but does not cure the
documented 1 T resolution/path dependence or certify the lower PM-to-ordered
bridge.

**Cost decision.** The full driver contains 101 transverse-field points.
Moving from 65 to 128 nodes approximately doubles its ordered-profile work.
A 256-node full visual sweep is deferred until human inspection shows that the
128-node morphology adds useful information; it would approximately double
the new cost again.

**Evidence and checks.** The exact 128-node table was appended to
`docs/diagnostics/invzp_outer_wp2/wp2_filtered_profile_visual_census.mat`.
MATLAB Code Analyzer and `git diff --check` pass. This is a resolution change
within the existing third visual workaround, not a fourth trial method.
