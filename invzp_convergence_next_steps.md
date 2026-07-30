# Projected-spin ordered-state convergence: next-step plan

## Objective

Determine whether the projected-spin \(1/z\) equations possess a physically
admissible ordered solution at each \(H_{\mathrm{MF}}\) node, and distinguish
that question from failure of the present fixed-point iteration. Only after the
node solutions are certified should the \(H_{\mathrm{MF}}\) integral and the
real-axis susceptibility be rebuilt.

This plan does not relax convergence requirements or use finite outputs from
failed nodes.

## Verified starting point

The following conclusions are supported by production-equivalent node traces
at \(T=0.1\) K:

1. The static sector iterates \(K_0\) alone; \(G_{\mathrm{stat}}\) is an
   algebraic function of \(K_0\). There is no missing independent
   \(G_{\mathrm{stat}}\) equation or degenerate \((G_{\mathrm{stat}},K_0)\)
   manifold.
2. The current residual-only gate can accept inadmissible static solutions.
   At \(B_x=3\) T, seven nominally converged nodes had
   \(\widetilde G_0=-269\) to \(-165\ {\rm meV}^{-1}\), with
   \(1+J_q\widetilde G_0\) changing sign across the lattice spectrum.
   Cancellation across the sampled poles allowed the scalar closure residual
   to fall below \(10^{-10}\).
3. The proposed small-\(n_{01}\) explanation is false for the measured nodes.
   The electronic two-level gap is \(0.0493\)--\(0.473\) meV at 1 T and
   \(0.266\)--\(0.439\) meV at 3 T; \(n_{01}\ge0.9935\). The measured
   cancellation metric is normally \(O(1)\), not \(10^3\)--\(10^5\), and does
   not separate failed from converged nodes.
4. The code combines an electronic two-level self-energy with a full
   136-state electronuclear \(G_0\). This is a real representation mismatch to
   audit, but it has not yet been shown to cause the node failures.
5. The fixed 33-node rule is not the underlying cause. Permissive integration
   merely propagates failed and pseudo-root values into \(H_0\).
6. The approximately \(0.059\) ordered-state sum-rule residual is not a clean
   solver diagnostic at `Ecut=40` meV. At the healthy 4 T endpoint, the bare
   full-electronuclear Matsubara sum alone has a \(7.42\%\) cutoff deficit at
   40 meV, falling to \(4.25\%\) at 80 meV and \(2.22\%\) at 160 meV.
7. The inelastic-channel cancellation

   \[
   \widetilde G_0=\frac{G_0^{\rm inel}}{1+\Sigma_0}
   \]

   is exact only when the dressed elastic contribution vanishes. In the
   production hybrid,

   \[
   U=\frac{G_0^{\rm inel}}
           {1+\Sigma_0+K_0G_0^{\rm inel}},\qquad
   V=\xi G_0^{\rm el},\qquad G_{\rm stat}=U+V,
   \]

   and the exact relation is

   \[
   \frac{1}{\widetilde G_0}
   =\frac{1+\Sigma_0}{G_0^{\rm inel}}
    -\frac{V}{U(U+V)}.
   \]

   Therefore \(J_{0,\rm eff}|G_0^{\rm inel}|/(1+\Sigma_0)\)
   is an inelastic-only diagnostic, not an unconditional admissibility
   precheck.
8. Evaluating the exact relation on `diag_rev3_check.mat` confirms that the
   stored 3 T nodes 21--27 are inadmissible:
   \(J_{0,\rm eff}|\widetilde G_0|=1.058\)--1.726. It also moves the stored
   1 T node 9 from an inelastic-only ratio of 0.981 to an exact ratio of
   1.00136, with negative uniform mass. Node 9 is therefore not a certified
   class-B/outer-map example and must be reclassified by the new solver.
9. The measured 42%/60% quantities are two-level dominant-sector shares, not
   elastic/inelastic shares. They cannot be used to bound the elastic
   correction. On accepted 1 T nodes 10--17 the exact correction changes the
   pole ratio by approximately 3--10%; on the anomalous stored node 21 it is
   approximately 53%.
10. Node 18 is not a hard-partition discontinuity. `invz_chi0z` uses a fixed
    \(10^{-8}\) meV degeneracy cutoff; its \(G_0^{\rm inel}\) at node 18 is
    unchanged for cutoffs from \(10^{-12}\) to \(10^{-4}\) meV. A finer
    \(h_z\) scan resolves a narrow, continuous electronuclear susceptibility
    peak. The geometric 33-node grid undersamples that feature.

These are the baseline claims. They should not be reopened without new
contradictory evidence.

Known rejected explanations and attempted fixes are summarized in
`invzp_convergence_dead_ends.md`. A rejected branch may be reopened only when
new evidence directly invalidates its recorded rejection condition.

## Work package 1: make the static physical domain explicit

Define

\[
x\equiv\widetilde G_0(0)
  =\frac{G_{\mathrm{stat}}}{1-K_0G_{\mathrm{stat}}},
\qquad
\Phi(x)=\left\langle\frac{x}{1+J_qx}\right\rangle_q .
\]

The physical interval is

\[
-\frac{1}{J_{\mathrm{sup}}}<x<0,\qquad
J_{\mathrm{sup}}=
\max\!\left[J_{0,\mathrm{eff}},\sup_{q,\nu}J_{q\nu}\right].
\]

For the present Ewald convention, verify rather than assume that
\(J_{\mathrm{sup}}=J_{0,\mathrm{eff}}\). The directional \(q\to0\) limit and
the uniform ordering mode must be included in this verification even though
the exact \(\Gamma\) point has zero weight in the Brillouin-zone average.

Using the closure \(G_{\mathrm{stat}}=\Phi(x)\), compute

\[
K_0(x)=\frac{1}{\Phi(x)}-\frac{1}{x},\qquad
\widehat R(x)=
\Phi(x)-G_{\mathrm{stat}}\!\left(K_0(x);
  \Sigma_0,\lambda_1,\lambda_2\right).
\]

Implement a bounded scalar solver on the physical \(x\) interval. It must:

- never evaluate across a lattice or uniform pole;
- reject \(x\ge0\), \(1+J_{\mathrm{sup}}x\le0\), nonphysical \(\xi\), or a
  nonpositive uniform mass;
- distinguish a zero of \(\widehat R\) from a discontinuity;
- search for all roots on the interval rather than assuming uniqueness; and
- report `no_admissible_static_root` separately from iteration failure.

The solver and its diagnostic table must use the exact production elastic
closure. At every candidate record

\[
G_0^{\rm inel},\ G_0^{\rm el},\ \xi,\ U,\ V,\ G_{\rm stat},\
K_0,\ \Sigma_0,\ x,\ 1+J_{\rm sup}x,\ D_{\rm uni},
\min_{q,\nu}|1+J_{q\nu}x|.
\]

Also record the inelastic-only value
\(x_{\rm inel}=G_0^{\rm inel}/(1+\Sigma_0)\) and its discrepancy from \(x\),
but do not use it as a hard gate unless \(V=0\) exactly or a rigorous bound on
the elastic correction proves that the classification cannot change. The bare
column \(J_{0,\rm eff}|G_0^{\rm inel}|\) is likewise descriptive only: it does
not furnish a universal threshold.

Monotonicity of \(\Phi\) alone does not prove uniqueness of
\(\widehat R\), because \(G_{\mathrm{stat}}(K_0)\) contains rational and
nonlinear elastic terms.

### Acceptance tests

1. Reproduce the healthy 4 T ordered static state within the existing
   numerical tolerance.
2. Reject the previously accepted 3 T pseudo-roots at nodes 21--27.
3. Reject the stored 1 T node-9 state if evaluated unchanged; then report
   whether a different admissible static root exists for the same frozen outer
   state.
4. In the \(m\to0\) limit, reproduce the ordinary paramagnetic static medium.
5. Show that every accepted root satisfies a positive denominator margin on
   the mesh and in the uniform mode.
6. Demonstrate seed independence of the inner static result.

### Execution checkpoint (2026-07-30)

The bounded physical-\(x\) inner solver is implemented in
`invz_projected/invz_emt_static_ordered.m`, with generic finite-resolution root
enumeration in `invz_common/invz_bounded_roots.m`. It searches every configured
sign-change interval and sampled \(|\widehat R|\) minimum, reports
discontinuities separately, records all requested exact elastic diagnostics,
and does not export an arbitrary branch when more than one admissible root
survives. `K0_seed` is compatibility-only and has no effect.

The focused gates in `invz_projected/tests/` pass:

- the verified production lattice has
  \(J_{\rm sup}=J_{0,\rm eff}=0.006421661809416939\) meV, while the
  \(16^3\) Gamma-dropped mesh maximum is
  \(0.006371725736157789\) meV;
- the frozen healthy 4 T state has one admissible root and reproduces the
  former state to \(|\Delta K_0|=4.02\times10^{-15}\) meV and
  \(|\Delta G_{\rm stat}|=5.67\times10^{-11}\ {\rm meV}^{-1}\);
- the unchanged stored 1 T node-9 state and 3 T nodes 21--27 are rejected;
  no alternative admissible static root is found at any of those frozen outer
  states;
- the \(m\to0\) result agrees with the ordinary paramagnetic scalar medium;
- accepted roots have positive mesh and uniform margins and are exactly seed
  independent; and
- classifications are unchanged over 2049, 4097, and 8193 scan points while
  the normalized endpoint margin is reduced from \(10^{-8}\) to \(10^{-12}\).

Machine-readable results and exact provenance are retained under
`docs/diagnostics/invzp_static_wp1/`; the compact execution record is
`docs/execution/invzp_convergence_journal.md`.

This completes the frozen-inner acceptance contract, but it does **not**
promote the full ordered susceptibility solver. With failed static values no
longer propagated, the legacy full 4 T field-profile run has no admissible
predictor root and only 5/33 admissible sampled nodes. Whether other coupled
outer roots exist is now a work-package-2 question. Until that classification
and the later field-integral rebuild are complete, the affected ordered output
must remain masked rather than restored through the old iteration.

## Work package 2: classify the outer self-energy problem

After the static map is physical and deterministic, rerun the 1 T and 3 T node
censuses. At fixed node parameters define the full outer residual

\[
\mathcal R_\Sigma[\Sigma]
  =\Sigma-\mathcal F[\Sigma],
\]

where every evaluation of \(\mathcal F\) uses the bracketed physical static
solve.

For each formerly failed node, record one of:

- admissible coupled root found;
- admissible static root exists, but the outer map is noncontractive;
- no admissible static root for the frozen outer state; or
- no coupled root found within a documented search domain.

Do not carry the old class-A/class-B labels forward as facts. They were assigned
using a mixture of residual-only acceptance and the inelastic-only
approximation. In particular, node 9 is provisional until this census, and
rejecting the stored 3 T pseudo-root does not establish that no alternative
coupled root exists.

Measure the local Jacobian or dominant eigenvalue of the converged outer map
before choosing Anderson, Newton--Krylov, or continuation. A larger iteration
budget is not evidence of existence.

### Acceptance tests

- Cold and warm starts converge to the same certified solution when only one
  admissible solution exists.
- Damping changes convergence speed, not the exported branch.
- All outer roots pass the static-domain, uniform-mass, dynamic-denominator,
  fixed-point-residual, and finite-value gates.

### Execution checkpoint (2026-07-30, partial)

`invz_projected/invz_ordered_outer_map.m` now defines
\(\mathcal F[\Sigma]\) without a hidden previous-iterate lambda. At fixed
\(\Sigma\), the dynamic \(K_{n>0}\) is algebraic and
\(\lambda_p(K_0)=\lambda_{p,\rm dyn}+w_0g_0^pK_0/\beta\) is evaluated inside
the bounded static search. Zero or multiple admissible static roots make the
outer map undefined or multivalued; no residual is exported in those cases.

Initial results are deliberately limited:

- the healthy 4 T state reproduces the former outer residual and has
  \(\lambda_{\rm dom}=-0.00877249\), stable over a finite-difference-step
  ladder;
- at \(\Sigma=0\), the map is defined at 13/34 retained 1 T profile nodes and
  6/34 retained 3 T nodes. This is a domain slice, not a coupled-root census;
- admissible coupled roots are found at 3 T nodes 28 and 33. They pass static,
  uniform, dynamic, finite-value, and outer-residual gates. Cold/halfway starts
  and mixes 1/0.5 agree within the \(10^{-8}\) outer tolerance;
- 1 T node 22 does not converge under admissible Picard. Undamped and
  evidence-selected mix 0.5 trajectories leave the static domain. At the last
  admissible mixed state, \(\lambda_{\rm dom}=1.35957\), so no positive scalar
  damping can make that local branch contractive.

Node 22 is therefore a measured noncontractive/domain-boundary case, but no
claim is made that another admissible coupled root does not exist. The broader
1 T/3 T coupled-root census remains open. Evidence and trial counters are in
`docs/execution/invzp_convergence_journal.md` and
`docs/diagnostics/invzp_outer_wp2/`.

A diagnostic-only 4 T reverse sweep on the unchanged 33-node profile certifies
nodes 33--29, exactly the five high-\(H_{\rm MF}\) nodes already reached by the
strict legacy per-node sweep. At node 28, the certified node-29 self-energy
gives `no_admissible_static_root`; this classification is unchanged over the
WP1 2049/4097/8193-point and endpoint-margin ladder. Thus reversing the coarse
node order is not a sufficient repair. The result does not rule out a coupled
branch between nodes 29 and 28; resolving that question requires bounded
adaptive field continuation, not another scalar damping value. Machine-readable
evidence is `wp2_4t_reverse_continuation.mat` and
`wp2_4t_reverse_boundary_audit.mat`.

Adaptive downward \(H_{\rm MF}\) continuation subsequently reaches node 28 in
six accepted steps. This exposed and corrected a separate numerical defect in
`invz_bounded_roots`: sign-bracket refinement stopped at the first residual
below \(10^{-10}\), while the equivalent closure residual could remain just
above \(10^{-10}\), making admissibility depend on grid phase. Sign roots are
now polished to `x_tol`, and adjacent sign brackets suppress duplicate
tangency polishing. The node-28 root is stable over the full 3-by-3
scan-density/endpoint-margin grid, is contractive, and passes every static and
dynamic margin gate. The direct coarse node-29-to-node-28 transfer still
fails, so both the polishing fix and adaptive parameter continuation are
load-bearing. Evidence is `wp2_4t_adaptive_boundary_continuation.mat` and
`wp2_4t_adaptive_target_audit.mat`.

This remains diagnostic-only. The next bounded packet is one additional 4 T
coarse gap, node 28 to node 27; whole-profile and production integration remain
premature.

That next gap reaches a certified component endpoint rather than node 27.
Adaptive continuation advances only to \(h=0.008042919679\); uniform and
supremum masses extrapolate to their common zero at
\(h_c\approx0.0080428632\), while the outer Jacobian remains contractive and
nonuniform/dynamic margins remain healthy. A below-edge frozen-seed proposal
has zero roots across the full 3-by-3 static resolution grid. Thus the
continued high-\(h\) component terminates at the physical
\(1+J_{\rm sup}x=0\) endpoint and cannot provide the profile down to \(h=0\).

The next WP2 question is no longer step-size tuning. It is whether a
disconnected admissible coupled component exists below the edge, followed by
thermodynamic selection and the WP5 integral review. Do not promote the
adaptive controller to production before those questions are resolved.
Evidence is `wp2_4t_adaptive_node28_to27.mat` and
`wp2_4t_adaptive_component_edge_audit.mat`.

A production-path resolution ladder at 1, 2, and 3 T subsequently changed only
`nH=33,65,129`. All nine final points remain masked because the separately
evaluated \(h=0\) predictor has no admissible static root. At 2 and 3 T, every
shared nested-grid node has the same verdict; 129 nodes merely locate the
success boundary one inserted half-step lower. At 1 T, the same shared field
node can succeed on a coarser path and fail on a denser path. The retained
self-energy shows that failed low-to-high nodes partially mutate the warm-start
carrier before a later static failure. Thus `nH` is not a repair, and the
production profile has an additional transactional-state defect. A rollback/
fresh-start comparison at only the contradictory shared nodes is the next
bounded diagnostic; it must not be confused with proof that a physical branch
extends to \(h=0\). Evidence is `wp2_hmf_node_resolution_census.mat`, with
interpretation and exact settings in the execution journal.

The accepted-state rollback is now implemented at every profile transition:
the \(h=0\) predictor, ordinary sweep, lower extension, and root refinement
compute candidate \((\Sigma,K_0)\) states but commit them only after the node
passes all existing gates. The exact 1 T ladder now has zero verdict
mismatches at every shared node. In particular,
\(h=0.0062416231096\) and \(0.0077454658051\) meV are accepted at all three
resolutions with matching \(\Sigma_0\) and \(K_0\); before rollback, each
failed on the denser path after inheriting a rejected last iterate. This
verifies the transactional defect and removes it as a source of subsequent
low-field evidence. It does not make the strict profile complete: the
predictor still has no admissible root and only 10/33, 20/65, and 39/129
positive-\(h\) nodes certify. Evidence is
`wp2_hmf_node_transaction_census.mat`.

A subsequent endpoint-only test replaced equation (45) by
\(h_0^{(2)}=h[r(0)+r(h)]/2\). The strict ordered \(h=0\) endpoint still fails
at every bare-ordered field in a 0.5 T census. Substituting a converged PM-limit
endpoint produces finite coarse roots only at 3.5, 4.0, and 4.5 T. In all
three cases the lower PM endpoint has negative uniform mass while both upper
bracket endpoints have positive mass. The trapezoid therefore jumps across a
static pole/component boundary; below 3.5 T even the PM lower fixed point is
uncertified. Endpoint linearization is not an all-field repair and must not be
wired into production. Evidence is `wp2_endpoint_trapezoid_census.mat`; exact
values and acceptance limits are recorded in the execution journal.

The two visual integration modes remain available as explicit diagnostic
options, but they are temporarily disabled in `invz_run_spectra.m` for a
strict-output comparison. The production driver now supplies only
`mix_outer=0.15` and `max_outer=3000`, so the ordered solver uses its
`full_profile` default with `nH=33`. No PM anchor, nonfinite-node filtering, or
bridged panel enters equation (45). Exact strict probes at 1, 4, and 4.68 T
all return `node_failed`; their \(h=0\) predictors fail and only 10, 5, and 4
of 33 positive-\(h\) profile nodes certify, respectively. Thus
`filtered_no_bracket` did not cause the general mask: it was a fail-closed
status inside a workaround that bypassed these broader strict-profile
failures. The strict comparison is expected to mask the ordered output.

### Feedback-informed implementation order

The external assessment in `feedback_root_Search_problem.md` and its checked
evaluation in `feedback_root_search_problem_evaluation.md` refine the
remaining implementation order as follows.

On every differentiable certified component, retain the exact diagnostic

\[
F'(h)=r(h)\,[1+J_{0,\mathrm{eff}}\widetilde G_0(h)].
\]

Thus a uniform-mass zero is a stationary point of \(F\) on that component.
The 4 T high-\(h\) component is strongly supported to end at such a
spinodal-like stability boundary while the reduced outer map remains
contractive. Do not extend that conclusion through an uncertified interval,
assume one uniform-mass crossing, or infer that every low-\(h\) failure belongs
to the same component. Electronuclear nonmonotonicity and disconnected coupled
components remain open.

1. **Pure accepted-state node transitions.** Refactor a fixed-\((B_x,h)\)
   evaluation into a pure function of its physical inputs and explicit seed.
   It must return a candidate state without mutating caller-owned continuation
   data. Commit `Sigma` and the static carrier only after every acceptance gate
   passes. First retest only the contradictory shared 1 T nodes from the
   33/65/129 census, with a documented fresh-start comparator. Success means
   that a shared node's verdict and accepted state are independent of rejected
   nodes earlier in the sweep. **Completed:** all shared-node verdicts now
   agree, including both previous contradictions; the strict-profile failure
   remains and is therefore not attributable solely to rejected-state
   contamination.
2. **Low-field small-\(M^2\) asymptotic packet.** Only after item 1 removes
   rejected-state history, run a bounded \(0.5\)--\(2.2\) T node census. If
   \(M^2=|\langle0|J_z|1\rangle|^2\) is the vanishing spectral weight, record
   \(M^2\), \(m^2/M^2\), \(M^2/n_{01}^2\), \(G_0^{\rm inel}\),
   \(G_0^{\rm el}\), the static \(U/V\) split, uniform mass, root status, and
   final response weight. The phase mask is decided before real-axis response
   evaluation, so a small susceptibility is not itself a valid mask reason.
   Check the exact removable product

   \[
   \frac{2m^2}{M^2}\gamma_0
   =\frac{2m^2}{n_{01}^2}
     [\lambda_1-(1-n_{01}^2)K_0]
   \]

   and derive the complete \(M^2\to0\) limit before altering arithmetic.
   Compare the present full-response/two-level-vertex hybrid with a
   dominant-sector-plus-bare-remainder construction. An implementation change
   is justified only if it preserves healthy fixtures and demonstrably removes
   a conditioning or representation failure; do not replace masked columns by
   zero susceptibility. **Conditioning checkpoint:** the 0.5--2.2 T census
   falsifies direct small-\(M^2\) arithmetic as the sampled mask cause. Failed
   predictors have \(M^2=27.55\)--30.04, whereas accepted high-\(h\) nodes
   include the measured minimum \(M^2=0.0220\). The exact prefactor
   reassociation changes measured-node arithmetic only at roundoff and the
   current form becomes nonfinite only at underflow scale. The next part of
   this packet is therefore the representation comparison, not a production
   arithmetic rewrite. **Representation checkpoint:** the fixed lowest-16
   electronuclear sector is smooth and spectrally isolated but captures only
   0.5--73% of the connected variance on the sampled accepted nodes; every
   non-predictor sample requires all 136 states to exceed 90% on the tested
   rank ladder. A closed electronic two-level control also retains a detached
   high-\(h\) component and fails at every \(h=0\) predictor. Neither cheap
   representation route supplies the missing equation-(45) path.
3. **Corrected saturation-anchor feasibility.** Starting from
   \(dH_0/dh=r(h)\), derive the asymptotic condition on
   \(\delta H=H_0-h\). If \(\delta H(\infty)=0\) is justified, test

   \[
   H_0(h)=h-\int_h^\infty[r(h')-1]\,dh'.
   \]

   Measure the tail decay on a certified high-\(h\) component and compare the
   upper and lower anchors where both are admissible. Do not replace \(r-1\)
   by \(\Sigma_0\) at finite ordered moment unless the elastic correction is
   explicitly proved absent. Transporting or extrapolating an integration
   constant from the paramagnetic side may be tested as a falsifiable
   diagnostic, but smoothness in \(B_x\) would not by itself certify that the
   same physical constant crosses the transition. **Completed feasibility
   test:** the closed two-level control has an integrable tail, but the
   production hybrid does not on a certified 1/4 T window extending to
   128 times the profile ceiling. Its fitted exponents are 0.593 and 0.811,
   while the full electronuclear connected variance is still far from
   saturation. Reaching the finite-Hilbert asymptote would leave the retained
   multiplet's physical validity. No upper-anchored production root is
   authorized.
4. **Reduced simultaneous residual.** Determine whether each nonzero-frequency
   implicit \(K_n\)--\(\Sigma_n\) relation can be eliminated uniquely, leaving
   an exactly equivalent residual in
   \((\lambda_1,\lambda_2,\lambda_3,x)\). The reduction must retain all
   physical denominator gates and reproduce healthy 4 T fixtures before
   globalized Newton, multistart, or deflation is attempted.
5. **Bordered pseudo-arclength continuation.** Once the simultaneous residual
   is verified, continue in \((h,\lambda,x)\) while recording the smallest
   singular value of the unbordered Jacobian, the uniform/supremum masses, and
   lattice/dynamic margins. Use these together to distinguish an ordinary
   fold from a physical-domain endpoint and to search for disconnected
   components. Record every uniform-mass crossing on each certified component;
   do not assume the current “last crossing” is unique.
6. **Representation and lattice audits.** Quantify how the
   full-electronuclear-\(G_0\)/two-level-vertex hybrid and controlled
   band-edge/DOS quadrature move the coupled component. The exact directional
   Gamma value is not a discrete member of the current \(\Phi\) average; its
   uniform susceptibility gate must not be removed as a mesh correction.
   Resolve the anisotropic \(q\to0\) density of states and measure its edge
   exponent rather than assuming a generic square-root cusp.

The first item is the lowest-risk algorithmic packet. The second is the
targeted entry point for the observed sub-1.8 T regime. The third is the
highest-value alternative construction of equation (45). Items 4--6 require
their preceding equivalence and branch-selection evidence; none is authorized
to relax the strict production acceptance contract.

## Work package 3: audit the electronic/electronuclear hybrid

Compare three internally defined calculations at representative healthy and
failing nodes:

1. a fully electronic two-level calculation;
2. a dominant electronuclear manifold with the remainder kept bare; and
3. the present hybrid, retained only as a diagnostic comparator.

Do not replace the dominant sector by the lowest electronuclear pair without
deriving its spectral weight. The full static response contains multiple
electronuclear and crystal-field transitions, and a collective susceptibility
peak is not a bare single-ion gap.

For each representation, compare:

- static susceptibility and uniform mass;
- \(\Sigma_0\), \(K_0\), and the admissible-root count;
- dominant and remainder contributions to the equal-time sum rule; and
- reduction to the established paramagnetic result near the continuous
  boundary.

Keep two decompositions distinct:

- elastic versus inelastic, defined by the static spectral/degeneracy
  convention used by `invz_chi0z`; and
- dominant versus remainder, defined by the electronuclear manifold retained
  in the \(1/z\) vertex construction.

They are not interchangeable. Resolve the narrow node-18 electronuclear peak
with an adaptive \(h_z\) scan, but do not treat it as evidence of a partition
artifact.

This work package decides whether the hybrid is a cause of the remaining outer
failures or only a quantitative approximation.

### Low-field representation checkpoint

The electronic/full-response mismatch is quantitatively large but is not, by
itself, the topological cause of the strict low-field mask:

- the framework's fixed-rank 16-state field-adapted manifold has a large
  16/17 gap, unit population mass, and continuous subspace overlap across the
  observed status boundary, but its connected-variance share falls to
  0.5--40% at the accepted high-\(h\) endpoints from 0.5 to 2.2 T;
- even 96 states remain below 90% connected-variance coverage at every
  non-predictor sample, so the existing 16-state dense-vertex construction is
  not a controlled comparator in this regime; and
- an internally closed electronic two-level model, evaluated over every
  source-grid node with accepted-state rollback, still has no admissible
  \(h=0\) predictor and exposes a shorter certified high-\(h\) component than
  the production hybrid.

The hybrid remains a quantitative systematic, particularly for the final
susceptibility, but replacing it with either cheap construction cannot restore
the required lower integration path. A definitive full-electronuclear vertex
would require a new feasible formulation; the existing dense full-space route
is budget-refused. The next equation-(45) work should therefore examine
anchoring/component topology rather than promote a low-rank bridge.

## Work package 4: lattice and cutoff convergence

Only after a physical node solver exists:

1. recompute complete coupled nodes on successively finer \(q\) meshes;
2. resolve the directional \(q\to0\) region analytically or with a controlled
   density-of-states quadrature;
3. compare \(\Phi(x)\), \(K_0\), \(\Sigma\), and pole margins, not merely frozen
   residual scans; and
4. add a high-frequency tail correction for the full multilevel \(G\), or
   demonstrate convergence with `Ecut`.

The \(q\)-mesh and Matsubara-cutoff errors must be reported separately.
Any density-of-states treatment must state the order of limits and whether it
computes an ordinary integral, a one-sided \(i0^\pm\) boundary value, or a
principal value. A principal-value prescription must not be used to admit an
in-band static pole as a physical ordered state.

## Work package 5: rebuild the \(H_{\mathrm{MF}}\) integral

Construct equation (45) only from certified node values. Replace the fixed
33-node trapezoidal rule with error-controlled quadrature once the node solver
is deterministic. The quadrature must adapt around narrow but continuous
single-ion features such as the node-18 electronuclear peak. An embedded ODE
integrator is a candidate only after branch selection and anchoring are fixed;
an endpoint coordinate such as \(h=h_c+t^2\) requires a measured square-root
edge and must not be assumed beforehand.

The \(h_z=0\) endpoint should initially remain a diagnostic boundary condition:

- if the admissible branch reaches it, integrate normally;
- if the branch terminates before it, report the field integral as undefined
  within the present approximation;
- use the boundary-linearized closure only where
  \(H_{\mathrm{MF}}\to0\) and verify its expected cubic error against the full
  integral where both exist.

Reconsider folds and the current “last crossing” rule only on the cleaned
profile. If more than one admissible ordered root survives, select it using a
validated free-energy or continuation criterion.

### Upper-anchor feasibility checkpoint

The exact upper form

\[
H_0(h)=h-\int_h^\infty[r(h')-1]\,dh'
\]

requires both \(\delta h(\infty)=0\) and an integrable \(r-1\) tail. The
production hybrid's certified high-\(h\) continuation fails the latter
prerequisite over the accessible asymptotic test window: factors 8--128 above
the ordinary profile ceiling give decay exponents 0.593 at 1 T and 0.811 at
4 T. The full electronuclear connected variance remains \(O(1)\), showing
that the retained single-ion space has not saturated. On the identical
window, a closed two-level response/vertex model gives exponents 1.733 and
1.568, so the construction is mathematically sound for Jensen's closed model
but is not transferable to the current hybrid.

No finite upper cutoff may be used as an implicit integration constant. A
production upper anchor now requires a controlled multilevel asymptotic
completion, not more tail nodes inside the truncated manifold.

### Temporary filtered-profile visual checkpoint

This experiment does not satisfy the WP5 certification prerequisite. It is
retained only for empirical morphology inspection.

A literal filtered rule that assigned \(H_0=0\) at the first finite accepted
positive-\(h\) node produced no root at 1, 2, 3, 3.5, 4, or 4.5 T: every
retained residual remained negative. That dead anchoring rule was removed.

The active visual rule instead preserves the two-endpoint experiment's
independent PM value \(r_0=1+\Sigma_0^{\rm PM}\), then inserts every finite
positive-\(h\) node from a geometric profile. In the original 65-node,
0.5 T-spaced census:

- 0.5 T remains `filtered_no_bracket`;
- 1.0--4.5 T return converged final ordered states;
- all positive-\(h\) nodes that entered those nine sampled quadratures were
  certified and formed one contiguous retained block; and
- for the successful 1--4.5 T roots, the only bridged interval was from the PM
  endpoint at \(h=0\) to the first certified ordered node, whose width was
  \(0.00607\)--\(0.01153\) meV.

At 4 T the root is stable over `nH=33,65,129`:
\(0.016026,0.016021,0.015978\) meV. At 1 T it moves materially:
\(0.028252,0.025673,0.024213\) meV, consistent with the already verified
rejected-state continuation-path defect. The 65-node mode is therefore a
useful visual interpolation, not resolution-converged low-field physics.

Human visual inspection of the 128-node production sweep with stronger
damping (`mix_outer=0.15`, `max_outer=3000`) found physically sensible
ordered output almost everywhere, except a masked sliver at 4.68 T; points
below about 1.8 T remained masked or physically implausible. This is an
observational morphology result, not a certification result.

For the next visual round, the production driver was advanced to `nH=256`
while retaining those damping and iteration settings. Exact probes give:

| \(B_x\) (T) | result | \(h_{\rm MF}\) (meV) | retained nodes | PM anchor |
|---:|---|---:|---:|---|
| 1.0 | converged | 0.041949 | 61 | unconverged |
| 1.5 | converged | 0.037035 | 74 | unconverged |
| 4.0 | converged | 0.015695 | 35 | converged |
| 4.68 | `filtered_no_bracket` | -- | 31 | converged |

No probe used an uncertified positive-\(h\) node. The 4.68 T failure therefore
survives the node doubling and is not explained by an isolated Inf/NaN sample.
The retained 128-node reference table used the earlier
`mix_outer=0.25`, `max_outer=2000` settings, so its numerical values are not a
controlled node-count comparator for this 256-node table. In particular, the
large 1 T change must not be assigned to resolution alone. The unconverged PM
anchors and finite lower bridges continue to preclude a low-field physical
claim.

For a direct output comparison, the production driver now temporarily
disables `filtered_profile_visual` without deleting its implementation. It
uses strict `full_profile`, the default 33-node profile, `mix_outer=0.15`, and
`max_outer=3000`. Direct probes give `node_failed` at 1, 4, and 4.68 T; the
strict \(h=0\) predictor fails at all three, and 10/33, 5/33, and 4/33 profile
nodes certify. This comparison separates the effect of bypassing the visual
filter from permanent code removal. It also shows that deleting only the
`filtered_no_bracket` status could not repair the underlying failure.
Exact census, endpoint provenance, residuals, and dispatcher/default checks
are retained in
`docs/diagnostics/invzp_outer_wp2/wp2_filtered_profile_visual_census.mat`.

## Deliverables and stopping rule

Each work package should leave:

- focused MATLAB tests;
- a compact machine-readable diagnostic table;
- the exact solver options and lattice provenance; and
- a short update to this document containing conclusions, not raw logs.

Update `invzp_convergence_dead_ends.md` whenever a proposed mechanism or fix is
falsified, including the rejecting observable and the condition under which it
would be legitimate to reconsider it.

The convergence problem is considered resolved only when ordered susceptibility
columns are either produced from fully certified states or masked with a
specific physical/numerical reason. A smooth spectrum alone is not an
acceptance criterion.
