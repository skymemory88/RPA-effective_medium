# invZ projected: code state and non-convergence diagnosis

**Canonical handoff document 1 of 2**

**Implementation baseline before this documentation-only handoff:** `c015333`
(`feat(invzp): safeguard cold edge acceleration`)

**Status date:** 2026-07-29

**Scope:** operational state supplementary to
[`invz_projected/README.html`](invz_projected/README.html), followed by the current
diagnosis of the ordered-state non-convergence problem.

The companion design document is
[`invzp_convg_fix.md`](invzp_convg_fix.md). The theoretical foundation that must not
be edited or discarded is
[`jensen_1z_framework.html`](jensen_1z_framework.html).

## 1. Executive handoff

The projected scalar code remains usable for its verified near-QCP window and for
exploratory visual spectra. Run
[`invz_projected/invz_run_spectra.m`](invz_projected/invz_run_spectra.m). The
driver's current uncommitted settings deliberately scan a wider field range than
the certified regression window; masked columns in that wider scan are expected
provenance, not missing plotting data.

The current production theory is **not globally well posed as a phase-selection
algorithm**. Its local nonlinear equations can have several accepted roots, folds,
disconnected branches, and intervals in which the resummed static map meets a
pole or no admissible fixed point is found. Cold starts, warm field continuation,
larger iteration limits, stronger damping, and the new safeguarded Aitken proposal
can change which local basin is reached, but none supplies the missing physical
rule for selecting among coexisting roots. In particular, the narrow data sliver
near **3.825 T fails on both cold- and warm-seeded paths**; the retained diagnostics
show that this is representative of a broad path-level failure, not just a slow
alternating edge mode.

The root issue is therefore not “insufficient iteration.” It is the combination
of:

1. a meromorphic finite-grid coupling transform;
2. poles in the legacy resummed static-medium closure;
3. a non-contractive, basin-sensitive full-\(\Sigma\) map; and
4. no stationary thermodynamic functional that ranks every admissible solution.

The recently added `signed_aitken1` option safely shortens one diagnosed
alternating cold-start tail and is default-off. It must not be presented as the
general convergence fix.

## 2. Repository state

### 2.1 Active modules

| Path | Role | Handoff status |
|---|---|---|
| `invz_common/` | Shared 136-state single-ion, lattice, coupling, and response utilities | Active shared engine |
| `invz_projected/` | Scalar projected \(1/z\) production and diagnostic workflows | Main runnable implementation |
| `invz_tensor/` | Tensor/component research path through the implemented A0–A4/a3d stages | Research code; not a full real-axis replacement for the projected path |
| `invz_functional/` | Isolated strict-functional and exact-oracle experiments | Non-production research oracles; no production dispatch |
| `docs/diagnostics/` | Reproducible convergence, QCP-grid, and branch-topology scripts and retained results | Evidence archive; not a unit-test suite |

There is no current in-tree `tests/` directory. Old documentation that quoted
`runtests(...)` commands described historical trees and has been retired. Current
verification consists of focused MATLAB scripts, embedded assertions, saved
diagnostic artifacts, and direct source-contract review.

### 2.2 Preserved user worktree state

This documentation handoff deliberately does not normalize the user's
uncommitted interactive setup. At the time of handoff,
[`invz_run_spectra.m`](invz_projected/invz_run_spectra.m) contains:

- `fieldContinuation = 'none'`;
- `useParallel = true`;
- `outerMix = 0.4`, `outerMax = 5000`;
- `fields = linspace(3.5,5,61)`; and
- `showPeaks = false`.

These are exploratory runner choices, not newly certified production defaults.
The untracked cold-/warm-seeding `.fig` files under `Data/` and the reference PDF
under `References/` are also intentionally preserved and are not part of this
documentation change.

### 2.3 Production defaults and opt-in paths

The high-level defaults in
[`invz_spectra_map.m`](invz_projected/invz_spectra_map.m) are:

- `ordered_1z = 'jensen'`;
- `static_medium = 'resummed'`;
- `field_continuation = 'none'`;
- `cold_acceleration = 'none'` unless explicitly forwarded as
  `'signed_aitken1'`;
- a finite \(16^3\) coupling grid with the usual `dpRng = 30` construction unless
  the caller overrides it.

`field_continuation = 'qcp_down'` is an opt-in sequential warm-seeding experiment.
It is useful for distinguishing a false local cold-start failure from a missing
branch, but it is not an equilibrium branch selector.

The strict reference media `strict_1z_dyson_ref` and
`strict_1z_bare_ref` remain selectable diagnostics. The legacy `resummed` medium
remains the production default because the strict candidate did not pass its
promotion gate; see §8.2.

At the lower-level point solver,
[`invz_solve_point_ordered.m`](invz_projected/invz_solve_point_ordered.m), the
compatibility default for `ordered_mode` is still `bare`. The map-level driver
selects Jensen explicitly. Callers that bypass the map must therefore set their
theory choices deliberately rather than infer the production configuration from
the low-level fallback.

### 2.4 Acceptance contract

For one ordered node, the independent state used by the coupled diagnostics is

\[
u = [\Sigma(i\omega_n);K_0],
\]

while dynamic \(K(i\omega_n)\) and
\(\lambda_{1,2,3}\) are derived from \(u\). An accepted state must satisfy all four
independent checks implemented by
[`invz_ordered_residual.m`](invz_projected/invz_ordered_residual.m):

- **A — outer self-energy map:** the freshly rebuilt
  \(\Sigma_{\mathrm{map}}\) agrees with \(\Sigma\);
- **B — static EMT closure:** the \(q\)-averaged elastic closure agrees with the
  supplied \(K_0\);
- **C — derived chain:** recomputed \(\lambda_{1,2,3}\) and ordered
  \(\Sigma\) agree with the exported state;
- **D — dynamic EMT identity:** every nonzero Matsubara slot agrees with a fresh
  dynamic-medium reconstruction.

Finite values, residual tolerances, reference-domain checks, and positive
stability margins are all mandatory. Every \(h\)-node consumed by the moment
profile must pass. The outer moment solve then uses

\[
F(h)=h-J_0m(h),
\]

accepting a nonzero increasing crossing only when the stability classifier is
positive. A mask means that this contract could not be certified; it must never
be silently replaced by a bare state or by a numerically convenient root.

### 2.5 Trace formats

- [`invz_hmf_ordered.m`](invz_projected/invz_hmf_ordered.m) emits raw ordered
  traces with `schema_version = 4`; version 4 adds cold-acceleration summaries and
  proposal records.
- [`invz_ordered_trace.m`](invz_projected/invz_ordered_trace.m) emits the compact
  wrapper trace with `schema_version = 3`.

Downstream diagnostic readers must check the corresponding schema rather than
assume the two version numbers are interchangeable.

## 3. What is verified

### 3.1 Finite-\(16^3\) visual regression

At \(T=0.10\) K, \(\mathbf q=0\), fields 4.60–4.90 T on 61 points, and a
0–6 GHz response window:

- 19 columns are Jensen ordered and 42 are stable paramagnetic;
- no susceptibility peak column is masked or non-finite;
- \(B_c^{1/z}=4.6925\) T is bracketed by 4.690–4.695 T.

This is a **finite-grid visual/regression fixture**, not a grid-converged critical
field.

### 3.2 Coupling-grid ladder

For \(N=12,16,20,24\), the measured critical fields are respectively

\[
4.6822845,\ 4.6927582,\ 4.6990936,\ 4.7029572\ {\rm T}.
\]

The contiguous accepted ordered widths are

\[
0.38228,\ 0.26776,\ 0.22409,\ 0.17796\ {\rm T},
\]

approximately \(N^{-1.076}\). The excluded \(\Gamma\)-adjacent coupling gaps are

\[
0.60124,\ 0.43930,\ 0.33866,\ 0.28143\ \mu{\rm eV},
\]

approximately \(N^{-1.103}\). When fields are aligned by each grid's own
\(B_c(N)\), the peak curve varies by only 0.38–0.53%; however, the absolute
critical-field alignment moves by 0.02067 T, corresponding to roughly
0.26–0.28 GHz. Thus the near-QCP mode shape is relatively robust while its
absolute field coordinate and accepted-width boundary are not converged.

## 4. Safeguarded cold-start acceleration

Commit `c015333` adds a resummed-only, default-off
`cold_acceleration = 'signed_aitken1'` path in
[`invz_ordered_node_solve.m`](invz_projected/invz_ordered_node_solve.m). It:

1. fits one signed factor to four successive full-\(\Sigma\) increments;
2. requires a tightly alternating factor in \([-0.99,-0.50]\), small ratio spread,
   and a small single-mode fit error;
3. constructs one vector Aitken proposal;
4. evaluates fresh full-\(\Sigma\) residuals at the ordinary and proposed states;
5. accepts only a finite, domain-valid, materially improving proposal; and
6. leaves the ordinary Picard iterate untouched otherwise.

At the diagnosed 4.4 T cold edge, proposal factors
\(-0.90861,-0.91183,-0.92129\) were observed at outer iterations
143, 147, and 152, and the solve closed at iteration 153. The accepted state was
bit-identical to the ordinary `mix_outer = 0.7`, `max_outer = 1000` result and
agreed with the warm solve to about \(2.1\times10^{-10}\) in \(\Sigma\).

At 4.3 T, the moment predictor accepted 13 nodes, but the failed node
\(h=0.002811855\ldots\) had interval ranks 16384/16376 and did not satisfy the
acceleration signature. No proposal was made. This is correct safeguarding:
the accelerator is deliberately narrow.

## 5. The 3.825 T failure and its relatives

At 3.825 T, 17 of 34 sampled nodes fail. The failures include:

- a predictor endpoint with \(A\)-block residual of order unity and interval rank
  zero;
- numerous interior static-medium rank/domain events;
- a later node at which one acceleration proposal is accepted but the
  \(A\)-block residual is still \(1.82\times10^{-5}\) after 1000 iterations.

Both cold and warm field seeding leave the visible sliver unconverged. This
directly rules out the hypothesis that all such masks are caused by one cold
alternating eigenmode.

At 3.6 T, a diagnostic Newton repair can recover three failed nodes and can make
one temporary end-to-end path close. That result establishes that some Picard
failures are local solver failures, but it does **not** establish uniqueness,
thermodynamic stability, or the correct global branch. The Newton repair is
therefore retained as a diagnostic kernel, not promoted to production.

## 6. Branch topology evidence

### 6.1 Multiple folds at 1.5 T

At 1.5 T, 25 initial states yield seven distinct \(h=0\) roots. Continuation
reveals several folds, including a high-\(h\) fold near
\(0.005243548\) meV and a low-\(h\) fold near
\(1.1304\times10^{-5}\) meV. These are topology-localizing estimates, not
certified fold coordinates. Both observed legs of root 6 return to the same
zero-field root within the registered cluster metric. Local continuity or
proximity in \(r=J_0m/h\) still cannot decide which connected component is the
equilibrium phase.

### 6.2 Coexisting accepted roots at 4.05 T

At one common \(h\) near 4.05 T, two states pass the numerical A–D contract, with
\(r\approx0.768\) and \(r\approx0.822\), separated by a state-space distance of
about \(1.037\times10^{-3}\). A 1369-state trace of the QCP-connected lower-\(r\)
component keeps \(F(h)<0\) through the permitted \(h\)-ceiling and therefore has
no closing endpoint in that domain.

This is the decisive counterexample to “choose the smoothest root.” Smoothness
can identify a connected numerical component; it cannot establish which
stationary state minimizes the thermodynamic potential, nor prove that every
competing component has been enumerated.

## 7. Mechanistic diagnosis

### 7.1 Meromorphic discrete transform

On a finite coupling grid, the scalar EMT relation evaluates a discrete sum of
rational terms. Its dependence on the static medium contains grid-specific poles.
Refinement moves and densifies these features; it does not turn a finite-grid
pole encounter into evidence that the physical state is absent.

### 7.2 Legacy resummed static denominator

The production static closure has an additional resummed denominator. Near its
reference boundary, interval/domain checks intentionally reject a state rather
than cross a pole or conceal it with a numerical floor. Many deep-ordered masks
are raised here.

### 7.3 Non-contractive coupled map

The outer update is a full Matsubara-vector problem. Dynamic \(K\), the three
\(\lambda\) moments, and the static closure feed back into every
\(\Sigma(i\omega_n)\). The map can alternate, stall, or enter another basin.
Changing damping and iteration count alters iteration dynamics but does not alter
the fixed-point set or define equilibrium.

### 7.4 Missing global thermodynamic selector

[`invz_deltaF_ordered.m`](invz_projected/invz_deltaF_ordered.m) is only a partial
free-energy difference. Its two retained-order routes disagree by about 13.7% in
the diagnosed regime, and it is explicitly outside its validated domain as an
absolute selector. It cannot rank the folds and disconnected components above.

The existing residual also cannot simply be integrated to make a scalar
potential. In the tested positive-control coordinates its numerical Jacobian is
not symmetric—for example \(J_{12}\approx-1.6649\times10^{-3}\) while
\(J_{21}=0\). This does not prove that no functional exists in better physical
variables; it proves that the present residual is not already a gradient in its
current coordinates.

## 8. What has been ruled out

### 8.1 Numerical workarounds

The following may be useful diagnostics but are not root fixes:

- more Picard iterations or different scalar mixing;
- cold Aitken acceleration outside its validated alternating-tail signature;
- warm continuation from the QCP;
- local Newton or trust-region repair without global branch enumeration;
- selecting the root nearest the previous field or the smoothest \(r(h)\);
- pole floors, clipped denominators, or substituting the bare state;
- masking and interpolating across a failed ordered column;
- broadening the real-axis spectrum.

The final item is especially irrelevant: the failure occurs during imaginary-axis
state construction. Once a state is accepted, real-axis response evaluation is
healthy in the verified window.

### 8.2 Strict-medium candidate

The strict reference medium removes the particular legacy resummed pole but
failed its Gate-0 promotion clauses, including the required reference
correspondence, phase coverage, and regression conditions. In the tested
strict-ring implementation the ordered state persisted through roughly 4.6 T,
no state was obtained from 4.7–4.85 T, and the symmetric phase returned near
4.9 T. Absence of the legacy pole therefore does not by itself produce a complete
thermodynamic theory.

The frozen Gate-0 constants are `crit_tol = 1e-6`,
`omit_promote = 0.10`, `pole_cont_tol = 1e-3`, and
`ref_margin = 1e-6`. The aggregate predicate in
[`invz_gate0_aggregate.m`](invz_projected/invz_gate0_aggregate.m) fails when:

- **(a)** a solved-path node has an invalid reference-denominator status;
- **(b)** node/predictor/final-root ledger coverage is incomplete;
- **(c)** the first-omitted-term estimate is non-finite or exceeds 0.10 at an
  otherwise valid node;
- **(d)** the \(r\) or critical-mass observable is unresolved/discontinuous at a
  local-denominator crossing or worsens from 65 to 129 profile nodes; or
- **(e)** an ordered endpoint or required PM control is absent, non-finite, or
  unstable.

The measured strict candidate failed clauses (a), (c), and (e). The old frozen
coupling digest used by this driver also conflicts with the current reproducible
digest; see §9. These are promotion failures, not tolerances to relax after seeing
the output.

### 8.3 Empirical smooth-\(r\) selection

The retained branch-graph tools are valuable topology probes. The proposed
smooth-\(r\) objective is not a formal equilibrium principle and is retired as a
production solution. It may not be used to turn numerical continuity into a
claim of thermodynamic stability.

## 9. Provenance caveat

The current cached and freshly rebuilt \(16^3\) coupling sets agree bit-for-bit
with digest

`499922e6c9f8c44d51b5db06486aac345b6226d1f8096713d20916ca78612cb5`.

An older frozen Gate-0 constant in
[`invz_gate0_report.m`](invz_projected/invz_gate0_report.m) is

`ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17`.

Current extrema and physical anchors agree, but the digest provenance has not
been reconstructed. Any future quantitative gate must freeze the constructor,
options, ordering, and digest together; do not silently update the old constant.

## 10. Retained diagnostic artifacts

The explanatory README files formerly beside these artifacts have been folded
into this document. The executable and data artifacts remain.

### Cold acceleration

- `docs/diagnostics/invzp_cold_accel_2026-07-29/invzp_cold_acceleration_gate.m`
- `docs/diagnostics/invzp_cold_accel_2026-07-29/cold_acceleration_results.mat`

### Near-QCP grid and response

- `docs/diagnostics/invzp_qcp_grid_2026-07-28/invzp_qcp_coupling_scan.m`
- `docs/diagnostics/invzp_qcp_grid_2026-07-28/invzp_qcp_state_grid_gate.m`
- `docs/diagnostics/invzp_qcp_grid_2026-07-28/invzp_qcp_peak_grid_gate.m`
- `docs/diagnostics/invzp_qcp_grid_2026-07-28/invzp_area_rule_oracle.m`
- retained `coupling_scan.mat`, `state_grid_gate.mat`, `peak_grid_gate.mat/.png`,
  `area_rule_grid.mat`, and `edge_pair_trace.mat` in the same directory

### Solver repair kernels

- context, node-construction, equation, Jacobian, and Newton helpers under
  `docs/diagnostics/invzp_solver_stability_2026-07-27/`

### Branch topology

- root enumeration, pseudo-arclength, fixed-\(h\), fixed-\(w\), clustering, graph,
  and event-oracle scripts under
  `docs/diagnostics/biased_smooth_r_2026-07-28/`

These scripts are research diagnostics. Unless an individual script says
otherwise, a successful run establishes its bounded numerical claim only; it is
not an acceptance gate for the production theory.

## 11. Handoff decision

The local solver can still be strengthened for discovery, but the production
problem cannot be closed by another seeding or acceleration heuristic. The next
substantive work should implement the stationary-functional program in
[`invzp_convg_fix.md`](invzp_convg_fix.md):

1. derive state equations and free energy from one source-dependent functional;
2. enumerate and solve all stationary branches with globalized methods;
3. select stable phases by that same functional and its Hessian; and
4. derive the response from the same retained diagrams.

Until those gates pass, the projected Jensen code should be handed over as a
well-instrumented but locally valid research implementation: trustworthy in the
explicitly verified windows, transparent when it masks a state, and not certified
as a global equilibrium solver in the moderate/deep ordered regime.
