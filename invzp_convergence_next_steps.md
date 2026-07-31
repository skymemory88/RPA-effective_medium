# Projected-spin ordered-state convergence: executed plan and evidence record

## Final status (31 July 2026)

The production objective is complete. The current driver uses
`missingAreaFactors = 1.0`, and the user-confirmed 101-point sweep is smooth
and finite from 0 to 9 T, including exact zero and 4.68 T. The separately
tested factors 0.75 and 1.5 remain incomplete at 4.68 T and are therefore not
enabled in that full sweep. Strict `full_profile` remains fail-closed.

The staged priorities below are retained as an executed investigation record;
their imperative wording is historical, not an open work list. The final
diagnosis is in `converg_diagnosis.md`.

## Original objective

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

An independent threshold audit confirms the limitation. The closed,
inelastic-only two-level identity would require
\(\Sigma_0>J_{\rm sup}M^2g_0-1\), but a certified production-hybrid node at
\(B_x=1\) T and \(h=0.0062416231096\) meV has
\(\Sigma_0=-0.29267261\) while that expression gives \(0.533395\); the exact
supremum mass is \(0.637269>0\). The hybrid uses a different full
electronuclear inelastic weight and a nonzero dressed elastic contribution.
Furthermore, the available triangle-inequality bound on \(|\Sigma_0|\), using
\(\max_q|J_q|\), is larger than the required closed-model predictor threshold
at every sampled 0.5--2.2 T field. No solver-free exclusion region has
therefore been certified, and this diagnostic must not displace the coupled
root search.

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

### Consolidated implementation order

The following order incorporates the retained conclusions of the convergence
diagnostics and subsequent independent review. It does not depend on the
deleted feedback documents.

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
   test:** the closed two-level control has fitted tail exponents above the
   integrability threshold on the sampled window, but the production hybrid
   does not on a certified 1/4 T window extending to 128 times the profile
   ceiling. Its fitted exponents are 0.593 and 0.811, while the full
   electronuclear connected variance is still far from saturation. Reaching
   the finite-Hilbert asymptote would leave the retained multiplet's physical
   validity. No upper-anchored production root is authorized.
4. **Controlled missing-area approximation.** Let \(h_e\) be the lower edge of
   a selected certified high-\(h\) component and isolate the unresolved interval
   as
   \[
   A(B_x,T)=\int_0^{h_e}r(s)\,ds.
   \]
   Solve
   \[
   F_A(h)=A+\int_{h_e}^{h}r(s)\,ds
          -J_{0,\mathrm{eff}}m(h)=0
   \]
   using only certified nodes above \(h_e\). Treat \(A\) as an explicit model
   quantity with an ensemble or interval of defensible completions; never hide
   it by deleting nonfinite nodes or assigning their contribution to zero.
   Propagate the \(A\) uncertainty into \(h_{\rm MF}\) and every reported
   susceptibility observable. This is the parallel practical-output path, not
   a claim that the missing constrained component has been solved.
   **Resumed and implemented as an opt-in approximation:** the census uses
   three continuous positive linear completions with
   \(r(0)/r(h_e)=0.5,1,2\), equivalently
   \(A/[h_er(h_e)]=0.75,1,1.5\). On the unchanged 33-node strict profiles,
   all three members bracket roots at every field with a terminal contiguous
   certified component in the sampled 0.2--1.2 T and 2.5--3 T ranges. The
   coarse 33-node census has no such component at 0.2 T, while the 129-node
   production profile and a 257-node check agree on one there. Conversely,
   0.5 T has no terminal component at 129 or 257 nodes and remains masked.
   The production map evaluates all members, exports their spread in
   \(h_{\rm MF}\), \(\Sigma_0\), moment,
   uniform mass, peak energy, and the full real-axis susceptibility, and labels
   the branch assumption as
   `picard_attracting_contiguous_high_h_component`. It neither uses a PM
   endpoint nor bridges a rejected node. The scalar mode is
   `missing_area_approx`; the map-level mode is `missing_area_ensemble`.
   The opt-in spectra driver uses a user-tunable 129 profile nodes; a 257-node
   spot ladder is retained as a separate resolution systematic rather than
   folded into the area band. A two-sided adjacent-field retry now recovers the
   0.36 and 0.45 T cold-pass masks at 129 nodes for every ensemble member. It
   freezes the original phase labels, seeds independently from the nearest
   accepted ordered fields on each side, and accepts only if both gated results
   agree. It therefore does not retry 0 T or bridge the boundary at 4.68 T.
   A separate ordered-boundary retry now handles only the final frozen
   `phase==0` ordered-side sample. It requires two untouched lower ordered
   sources, an accepted PM cold point above, an independently converged
   negative PM mass at the target, identical recovered components, and
   explicit \(D_{\rm uni}\) and \(F'\) margins. It never retries `phase==2` or
   seeds from a recovered state.
   `full_profile` remains the option-free strict default.
   **Near-boundary update:** after changing the visual grid to 111 points, the
   remaining mask occurs at 4.663636 T. Its cold factor-one support ends at
   0.991798 (0.99289 at 257 nodes), so factor one has no equation-(45) bracket;
   damping and iteration-budget increases do not change that fact. Factors
   0.85--0.95 converge cold. More importantly, independent descending sweeps
   seeded from accepted 4.5 and 4.581818 T factor-one states recover the same
   factor-one root at 4.663636 and 4.68 T. The initial direct-seed ladder then
   became seed-dependent at 4.70 T, which temporarily blocked wiring. The
   apparent derivative sign disagreement was a diagnostic substitution error:
   the exact identity uses \(\widetilde G_0\), not \(G_{\rm stat}\). Narrow
   fixed-\(h\) sections give \(F'=0.028518\) and 0.020014 at 4.663636 and
   4.68 T, agreeing with local secants within \(1.3\times10^{-5}\) relative.
   A direct long-jump 129/257-node ladder agrees through 4.69 T, loses
   two-seed agreement at 4.70 T, and fails from both seeds at 4.71 T. A new
   fine-step trace shows that this was a continuation-distance failure: two
   independent histories agree at both resolutions through 4.7188 T. At the
   current 101-grid 4.68 T target, direct untouched 4.50 and 4.59 T sources
   agree, the independent PM mass is negative, and the 129/257 peak shift is
   0.00110 GHz. The narrowly gated central-member boundary retry is now wired;
   noncentral area members remain independently incomplete.
   **Exact-zero update:** the separate \(B=0\) obstruction was an unhandled
   rejected predictor, not a non-unique final ordered state. The opt-in
   missing-area ensemble now rejects the unresolved \(h=0\) two-level node
   without committing it and integrates only the certified positive-\(h\)
   suffix. All three area members have finite, continuous exact-zero spectra;
   strict `full_profile` remains masked.
5. **Reduced simultaneous residual.** Determine whether each nonzero-frequency
   implicit \(K_n\)--\(\Sigma_n\) relation can be eliminated uniquely, leaving
   an exactly equivalent residual in
   \((\lambda_1,\lambda_2,\lambda_3,x)\). The reduction must retain all
   physical denominator gates and reproduce healthy 4 T fixtures before
   globalized Newton, multistart, or deflation is attempted. At fixed
   \((\lambda_1,\lambda_2,\lambda_3,K_0)\), each \(n>0\) equation is scalar
   because the ordered self-energy is affine in \(K_n\). This is a useful
   structural reduction only if existence, uniqueness, and singular cases of
   that scalar elimination are established over the full admissible search
   domain; it must not be assumed from the healthy branch alone.
   **Implemented and equivalence-verified:** for \(n>0\),
   \(\Sigma_n=c_n-d_nK_n\) reduces the medium equations to one bounded scalar
   equation. The diagnostic proves a derivative uniqueness bound when
   available and otherwise performs finite root enumeration without accepting
   multiple or unresolved results. The exact remaining residual is
   \[
   R=(\lambda-\lambda[K],\,x-\widetilde G_0),
   \]
   where the static form is pole-cancelled but algebraically equivalent to
   \(\Phi-G_{\rm stat}\). It reproduces the healthy 4 T dynamic state to
   \(4.94\times10^{-12}\), and an accepted 1 T node has reduced residual
   \(7.11\times10^{-15}\) and legacy outer-map residual
   \(2.02\times10^{-11}\).

   At 1 T and \(h=0\), forward/reverse continuation solves the three moment
   equations to at most \(5.59\times10^{-12}\), while
   \(x-\widetilde G_0\) remains positive from 1158.87 to
   2210.77 meV\(^{-1}\) across
   \(s=-J_{\rm sup}x=10^{-5}\ldots0.99999\). The apparent sign change of
   \(\Phi-G_{\rm stat}\) is a removable local \(G_{\rm stat}\) pole, not a
   root. An eight-start full solve reaches the search bounds with nonzero
   residual. At noncontractive node 22, two normalized starts reach nearly
   the same admissible nonzero-residual state. These results reject a nearby
   root on the sampled moment branch but do not yet exclude disconnected or
   off-profile roots.

   **Independent-method recheck:** a separate solver iterates the original
   740-component \(\Sigma\to K\to\lambda\to\Sigma\) map at fixed \(x\),
   without the scalar dynamic eliminator or least-squares moment solver. Its
   118-point \(h=0\) forward/reverse profiles agree in \(\Sigma\) to
   \(9.63\times10^{-14}\), agree with the reduced moments to
   \(2.57\times10^{-11}\), and retain positive dynamic lattice/medium masses
   of at least 0.46514/0.45182. Five disparate seeds at \(s=0.5\) agree to
   \(3.41\times10^{-14}\). The static residual stays in
   \([1158.79,2210.78]\ {\rm meV}^{-1}\).

   The same unreduced method supplies a positive control: it finds the known
   accepted-node crossing at \(s=0.362731098554\) with legacy outer residual
   \(1.65\times10^{-14}\). Direct EMT substitution gives dynamic closure
   residuals below \(3.55\times10^{-14}\), and the removable-pole identity
   agrees to \(9.09\times10^{-13}\).

   A 260-evaluation global surrogate/pattern search of the complete rigorous
   four-variable box again reaches the endpoint bounds with
   \(\|R_{\rm scaled}\|_\infty=1.083\). At node 22, global search also finds
   no root, while centrally differenced scaled Newton from the two prior
   candidates remains at \(4.053\times10^{-2}\). The local reduced Jacobians
   have reciprocal condition estimates \(2.92\times10^{-6}\) and
   \(5.86\times10^{-7}\). Rechecking every high-penalty global trial finds no
   multiple or unresolved dynamic scalar equation. These alternative methods
   materially strengthen the conclusion, but finite global budgets still do
   not certify root absence.
6. **Bordered pseudo-arclength continuation — completed locally at 1 T.**
   The certified root was followed through a regular saddle-node at
   \[
   (h,s)=(0.0054789314231\ {\rm meV},\,0.6565777265).
   \]
   A fixed-\(s\) refinement gives \(dh/ds=-8.0\times10^{-10}\),
   \(d^2h/ds^2=0.0220\), and a smallest fixed-\(h\) Jacobian singular value
   \(3.97\)--\(4.05\times10^{-8}\), while every physical mass remains
   positive. The return sheet reaches
   \((h,s)=(0.0097188640,0.999998661)\), where the supremum/uniform masses,
   rather than the dynamic masses, approach zero. The down-and-back
   normalized-state discrepancy is \(6.24\times10^{-9}\), and none of the 33
   accepted points requires dynamic-root enumeration.

   A separate 740-component full-\(\Sigma\) calculation independently finds
   two roots at \(h=0.006\) meV,
   \(s=0.412077554837\) and \(0.809166172289\), with legacy outer residuals
   \(9.25\times10^{-14}\) and \(4.61\times10^{-13}\). Thus the fold and
   return sheet are not artifacts of the reduction or corrector. Node 22 lies
   below this component's minimum \(h\); smaller fixed-\(h\) steps and damping
   cannot extend this component to \(h=0\).

   The original outer-map dominant eigenvalue is 0.4173036 on the low-\(s\)
   root, 2.6434068 on the return-sheet root, and 0.9999998 at the fold. The
   values agree across finite-difference steps \(10^{-6}\) and \(3\times
   10^{-6}\), with eigen-residual below \(1.3\times10^{-8}\). This completes
   the local iteration-stability classification: Picard loses contraction at
   the saddle-node. It does not supply a thermodynamic selector.

   Remaining obligations are bounded disconnected-component coverage below
   the fold, an independent thermodynamic consistency route, and a
   transverse-field repeat before treating the 1 T topology as universal.

   **Fixed-state coupling-path checkpoint.** At \(h=0.006\) meV, scale only
   the longitudinal fluctuation interaction by \(\rho\), holding
   \(B_x,h,T,J^{aa}_0\) and therefore the constrained single-ion state fixed.
   The interaction-energy normalization is
   \[
   \frac{\partial(\Phi/N)}{\partial\rho}
   =\frac{1}{2\beta\rho}\sum_n w_nK_n(\rho)G_n(\rho).
   \]
   The low-\(s\) root reaches \(\rho=0.001\), with anchor error
   \(8.66\times10^{-6}\), tail power 1.00002, and candidate
   \(\delta\Phi/N=-1.02656\times10^{-4}\) meV. Coarsening the coupling grid
   changes this integral by \(3.56\times10^{-7}\) meV.

   The return root instead approaches \(s=1\). Prescribing
   \(s=1-\epsilon\) and solving the exact four-variable residual for \(\rho\)
   gives
   \[
   \rho_*=0.7721812909\quad\text{(linear extrapolation)},\qquad
   0.7721812622\quad\text{(quadratic)},
   \]
   a spread of \(2.87\times10^{-8}\). Over the endpoint sequence, the
   uniform mass divided by \(\epsilon\) approaches 1.3274, the other mesh and
   dynamic masses remain at least 0.0103, and the endpoint Jacobian singular
   value remains 0.299. The sheet therefore ends on the uniform stability
   boundary rather than reaching the noninteracting anchor.

   This result supplies an adiabatic-connectivity discriminator, not a
   completed thermodynamic ranking. A direct common-anchor potential
   difference is undefined, and the hybrid closure has not been proved
   stationary or \(\Phi\)-derivable. The next thermodynamic task is an
   independent field/local-\(S_3\) derivation and a derivative check against
   equation (45). Do not assign the unavailable upper-sheet integration
   constant by truncation or by copying the lower-sheet anchor.

   **Fold-anchored field-potential checkpoint.** Jensen's differential
   identities imply, for two sheets meeting at one fold,
   \[
   \Phi_H(h)-\Phi_L(h)=
   \int_{h_f}^{h}[r_H(u)-r_L(u)][m(h)-m(u)]\,du .
   \]
   This cancels the unknown field-integration constant without using the
   paused \(h=0\) missing area. At 1 T and \(h=0.006\) meV the result is
   \(-1.38309\times10^{-6}\) meV per ion. A separately nested evaluation
   differs by \(1.72\times10^{-9}\) meV and alternate-node coarsening changes
   it by \(6.93\times10^{-9}\) meV. Extension to the return-sheet uniform
   boundary reaches \(-1.02210\times10^{-4}\) meV with no crossing.

   The sign conditionally favors the high-\(s\) sheet, opposite to the
   zero-coupling-connectivity discriminator. This is not a contradiction:
   the high sheet is interaction-born and lacks the coupling anchor. The
   blocking question is whether the production hybrid is a stationary
   \(\Phi\)-derivable approximation and whether both sheets share the same
   constrained field-potential constant at the fold. Until that is derived,
   neither diagnostic is a production selector.

   **Same-ion consistency disposition.** Jensen's ordered derivation obtains
   the moment, \(G_0\), elastic weights, and cumulant vertex from the same
   two-level ion. The production hybrid does not: its moment and \(G_0\) are
   136-state electronuclear objects, while its vertex and \(\xi\) are
   electronic two-level objects. The full bare response passes an independent
   derivative check to \(2.57\times10^{-7}\), but
   \(G_{0,\rm full}/G_{0,\rm twolevel}=0.147\)--0.889 and the corresponding
   reference moments differ by 0.095--1.326 over the sampled fields. At exact
   hybrid roots, using the Jensen-consistent two-level static weights changes
   \(G_{\rm stat}\) by 18.9% up to a factor 178.

   Consequently, equations (45)--(46) do not furnish an exact thermodynamic
   selector for the current hybrid. The fold integral remains a reproducible
   heuristic. Promoting it requires a derived multilevel connected four-point
   vertex or an explicit mixed functional that reproduces the existing
   stationary equations.

   **Transverse-field checkpoint.** Exact four-variable correction of the
   1.2 T high root, independent of inner Picard contraction, certifies two
   roots at 0.5, 0.8, 0.9, 1.1, and 1.2 T. The corresponding
   \(r_H-r_L\) values are -0.592, -0.716, -0.625, -0.352, and -0.177.
   The high branch reaches the \(s=1\) boundary near 1.21672 T. On an
   explicitly extrapolated diagnostic \(h(B_x)\) path, both roots persist to
   0.355 T but graph correction fails at 0.350 T; this is a bracketed
   merger/termination, not yet a certified regular fold.

   At 2.5, 2.7, and 2.9 T, by contrast, the exact moment branch extends to
   \(s=0.999999\) with static residual -33.01, -29.39, and -27.64. Fifteen
   constrained high-\(s\) full-residual starts find no additional root; their
   best residuals remain 1.45--1.55. This is strong finite one-root evidence
   on the sampled branch, not an interval uniqueness proof.

   Evidence:
   `docs/diagnostics/invzp_outer_wp2/wp2_reduced_pseudoarclength_1t.mat`,
   `wp2_reduced_fold_refinement.mat`, and
   `wp2_pseudoarclength_original_equations_audit.mat`, plus
   `wp2_fold_sheet_stability_audit.mat`,
   `wp2_coupling_constant_path_audit.mat`, and
   `wp2_coupling_endpoint_refinement.mat`,
   `wp2_local_field_potential_audit.mat`,
   `wp2_local_field_potential_extension.mat`,
   `wp2_hybrid_functional_consistency_audit.mat`,
   `wp2_multifield_crossfield_continuation.mat`,
   `wp2_lowfield_02_branch_extension.mat`, and
   `wp2_highfield_exact_section_audit.mat`.
7. **Representation and lattice audits.** Quantify how the
   full-electronuclear-\(G_0\)/two-level-vertex hybrid and controlled
   band-edge/DOS quadrature move the coupled component. The exact directional
   Gamma value is not a discrete member of the current \(\Phi\) average; its
   uniform susceptibility gate must not be removed as a mesh correction.
   Resolve the anisotropic \(q\to0\) density of states and measure its edge
   exponent rather than assuming a generic square-root cusp.

Items 1--3 have completed their diagnostic role: rejected-state history was
removed, small-\(M^2\) arithmetic and the two cheap representation replacements
were rejected as complete repairs, and the production upper anchor failed its
controlled-tail prerequisite. Item 4 is now an explicitly resumed, opt-in
approximation route. Its first census and end-to-end ensemble validation are
complete; they do not change the rigorous-route conclusions or acceptance
contract.
Items 5 and 7 remain active. Item 5 has passed derivation, healthy-root
equivalence, the local 1 T component-topology gate, and a fixed-state
coupling-path audit. It now also has a fold-anchored local-field potential and
cross-field exact-root continuation. The same-ion audit blocks transferring
Jensen's thermodynamic interpretation to the present hybrid; a rigorous
selector therefore requires a multilevel vertex or explicitly derived mixed
functional. Other remaining obligations are explicit disconnected-root
coverage below the 1 T fixed-\(h\) fold and matched topology repeats at
representative fields. Item 6 is complete for the known 1 T component and for the necessary
cross-field \(r_H<r_L\) sign at five fields; full field-specific fold/potential
repeats remain. Before any branch is used in equation (45), establish a
thermodynamic/constrained selector that does not assume a common coupling
anchor for the return sheet. Item 7 measures
representation dependence and supplies future approximation guidance. None is
authorized to relax the strict production acceptance contract. The resumed
missing-area mode is a separately labeled approximation whose absent-option
oracle is bitwise identical to the retained strict 1 T profile.

### Clean-slate decision

Do **not** reset the repository to `ckpt-00-start` (`a724862`) or discard the
post-baseline corrections before implementing item 5.

The current reduced-residual plan is independent of the previously rejected
algorithmic strategies: damping ladders, denser fixed profiles, filtered
integration, endpoint bridging, upper anchoring, and cheap representation
replacement are not inputs to its equations. It is intentionally **not**
independent of the following verified infrastructure:

- the reduced baseline `693c518`, which already removed the earlier
  experimental root/trace/functional machinery;
- the bounded physical-\(x\) static closure and deterministic outer map from
  `1e238a3`;
- the sign-root polishing correction from `51f80cc`;
- accepted-state-only profile commits from `bca2896`; and
- the healthy, boundary, and failed-node fixtures accumulated after those
  commits.

Returning to `a724862` would restore superseded solver machinery and
residual-only pseudo-roots, then require re-deriving and re-validating the
physical-domain contract before the simultaneous residual could be trusted.
Returning only to `693c518` would still discard the exact primitives and tests
that define the new residual's admissible domain. Neither rollback increases
mathematical independence.

Use a **derivational clean slate on the current codebase** instead:

1. derive the reduced residual directly from equations (16), (17), (23),
   (39)--(42), without copying the deleted historical Newton or arclength
   implementations;
2. expose it through a new diagnostic-only entry point with immutable node
   context and no production dispatcher changes;
3. prove numerical equivalence to `invz_ordered_outer_map` on healthy fixtures
   before searching elsewhere;
4. test every eliminated scalar equation for multiple roots, loss of
   admissibility, and singular Jacobians; and
5. keep the existing strict solver as the comparison oracle until the new
   formulation passes all gates.

Reconsider a historical reset only if an equivalence test shows that one of
the retained primitives changes the governing equations rather than enforcing
their physical domain.

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
  non-predictor sample on the tested lowest-energy rank ladder, so the existing
  16-state dense-vertex construction is not a controlled comparator in this
  regime; and
- an internally closed electronic two-level model, evaluated over every
  source-grid node with accepted-state rollback, still has no admissible
  \(h=0\) predictor and exposes a shorter certified high-\(h\) component than
  the production hybrid.

The hybrid remains a quantitative systematic, particularly for the final
susceptibility, but replacing it with either cheap construction cannot restore
the required lower integration path. A definitive full-electronuclear vertex
would require a new feasible formulation; the existing dense full-space route
is budget-refused. Connected variance is an early-warning coverage diagnostic,
not a four-point-vertex convergence norm, so this result does not exclude every
possible reduced factorization. The next work should therefore implement the
reduced simultaneous residual and examine component topology rather than
promote the tested low-rank bridge.

### Exact reduced representation checkpoint

The 1 T \(h=0\) reduced equations distinguish two failure mechanisms:

- In the production hybrid, every tested dynamic frequency is uniquely
  eliminable and the three moment equations define a reproducible branch, but
  that branch misses the static closure by at least
  \(1158.87\ {\rm meV}^{-1}\) on the sampled strict \(x\) domain.
- In the internally closed electronic two-level control, \(\lambda=0\) is
  outside the dynamic positive-mass domain. Along
  \(\lambda=a\lambda_{\rm hi}\), where \(\lambda_{\rm hi}\) is the rigorous
  componentwise moment bound, defined scalar roots first occur at
  \(a=0.7125\) for all five tested \(x\) values. The three moment residuals at
  and above that entry point point toward smaller, undefined \(\lambda\);
  no simultaneous root occurs on the ray. Thirty of 60 ray points are
  defined, no unresolved scalar frequency is accepted, and at most one
  frequency requires finite enumeration.

The ray census is not an off-ray root exclusion. An independent off-ray
surrogate/pattern search adds 320 evaluations over the rigorous closed-model
box, including six known feasible-ray points. Its best candidate is
nonadmissible with \(\|R_{\rm scaled}\|_\infty=0.399\), and the rejected
trials contain no multiple or unresolved scalar root. Together these controls
make the hybrid construction unsupported as a necessary cause of the lower
obstruction; representation changes its observed mechanism and scale. A matched
16-state reduced-residual control is still not authorized: the existing
16-state response split is not paired with a converged 16-state connected
four-point vertex, and the rank census already shows inadequate
connected-variance coverage. The finite searches still leave an unsampled
off-ray root logically possible.

### Relation of Hilbert-space truncation to \(r(h)\)

The present production projected solver is not simply a 16-state calculation.
Its bare \(G_0\), including the ordered static inelastic/elastic split, is
computed from the full retained 136-state electronuclear ion. Its self-energy
and ordered vertex are nevertheless electronic two-level objects. The
field-adapted lowest-16 electronuclear subspace was tested as a possible
replacement and was not promoted to production.

Hilbert-space/vertex truncation can still affect \(r(h)\) through two distinct
mechanisms:

1. High-energy states that have little thermal population or static
   susceptibility weight can remain important virtual intermediates in the
   connected four-point vertex. Dressing a full 136-state response with a
   two-level, or inadequately converged 16-state, vertex can therefore give the
   numerator and denominator of \(r=G_0/\widetilde G_0\) incompatible
   polarization and tail behavior.
2. The retained 136 states are themselves only the \(J=8\), \(I=7/2\)
   electronuclear manifold. At the very large longitudinal fields needed to
   force numerical saturation, omitted higher electronic multiplets become
   relevant. Eventual saturation of the finite 136-state matrix is therefore
   not automatically a controlled physical asymptote.

The evidence separates the lower-path and tail questions. A fully closed
electronic two-level control still loses the low-\(h\) component, so a
16/136 truncation is not required to produce that obstruction and cannot yet
be called its cause. Conversely, the closed control has much faster fitted
\(r-1\) decay than the production hybrid, while the full electronuclear
connected variance remains \(O(1)\). Representation closure is therefore a
credible contributor to the slow upper tail and to the quantitative value of
the inferred missing area.

The representation audit should consequently compare, at matched nodes and
fields:

- \(r(h)\), the component edge, and the inferred \(A(B_x,T)\);
- final approximate susceptibilities and their missing-area sensitivity;
- internally matched response/vertex pairs at each feasible rank, rather than
  changing only \(G_0\) or only the vertex; and
- actual connected four-point-vertex convergence wherever feasible, using
  connected variance only as an early-warning proxy.

This audit may explain or reduce the approximation uncertainty even if it does
not restore a certified lower component.

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

### Controlled missing-area approximation path

**Status: explicitly resumed as a separate opt-in approximation.** The
specification below remains its acceptance contract. It does not alter the
negative rigorous-route result or the option-free strict production default.

The incomplete-profile observations justify a separate approximation program,
provided the omitted integral is exposed rather than silently discarded. For a
selected certified component beginning at \(h_e\), write

\[
A(B_x,T)=\int_0^{h_e}r(s)\,ds,\qquad
F_A(h)=A+\int_{h_e}^{h}r(s)\,ds-J_{0,\rm eff}m(h).
\]

The certified upper-component quadrature and the local susceptibility at a
certified final root remain first-principles outputs of the present
approximation. The inferred \(A\), and therefore the selected \(h_{\rm MF}\),
carry the additional model uncertainty.

Construct an \(A\) ensemble from explicitly named prescriptions, initially:

1. a PM-endpoint bridge, retained as an empirical completion rather than a
   stable constrained path;
2. matched closed-model completions, used only after measuring their bias where
   the hybrid and closed calculations overlap;
3. field/temperature continuation from regions with a complete or independently
   calibrated path; and
4. bounded smooth completion families consistent with the available endpoint
   values, without inserting failed-node iterates.

Do not average mutually incompatible prescriptions into a falsely precise
central value. Report their spread and keep a prescription-specific result
when the ensemble is multimodal. If experimental data are used to determine
\(A\), label the result semi-empirical and validate it on data not used for the
fit.

For a simple root \(h_*\), the leading missing-area sensitivity is

\[
\Delta h_*\simeq-\frac{\Delta A}{F'(h_*)}
=-\frac{\Delta A}
 {r(h_*)[1+J_{0,\rm eff}\widetilde G_0(h_*)]}.
\]

This conditioning diagnostic is mandatory. Near a uniform-mass zero,
\(F'(h_*)\to0\), so apparent quadrature convergence can coexist with
uncontrolled physical root error.

An approximate susceptibility may be exported with status
`approximate_missing_area` only when:

- the quadrature over \([h_e,h_*]\) is converged using certified nodes;
- every \(A\) prescription being retained selects a root on the same certified
  component, away from a zero of \(F'\);
- the resulting \(h_{\rm MF}\), susceptibility, peak positions, and integrated
  weights are stable within declared tolerances over the \(A\) ensemble;
- the method passes available overlap/holdout tests where a complete path or a
  higher-fidelity reference exists; and
- the output records \(h_e\), every \(A\) prescription and range,
  \(F'(h_*)\), root spread, spectral spread, and all node/gate provenance.

Otherwise retain a diagnostic plot or a masked result with the specific
failure reason. Resolution convergence of the retained-node quadrature alone
is not an acceptance criterion.

### Executed remaining-mask plan

**Current visual status.** The driver currently uses the 101-point
`linspace(0,9,101)` grid, whose boundary-side sample is 4.68 T. The two-sided
adjacent-field retry removes the interior low-field slivers, the exact-zero
column is supplied directly by the opt-in missing-area machinery, and the new
ordered-boundary retry supplies the central factor-one 4.68 T column. The
driver currently requests only factor one. In the separate three-factor
diagnostic, the 0.75 and 1.5 members remain masked there, so no complete
sensitivity interval is claimed. The subsequent full production sweep
confirms a smooth, finite factor-one susceptibility across all 101 fields.
Option-free strict `full_profile` remains fail-closed.

#### Priority 1: ordered-side boundary near 4.663636--4.68 T

The next approximation work is restricted to the ordered side of the refined
PM mass zero, \(B_c=4.71897990927\) T. The 111-point driver samples
4.663636 T, 0.05534 T below this boundary. The mask has two separable causes:
the cold ascending solve reaches a shorter component whose fixed area factors
miss the root-support window, while ordered-side descending continuation
reaches a longer component and converges factor one.

Proceed in this order:

1. Re-evaluate a narrow fixed-\(B\), fixed-\(h\) section around each seeded root
   and compare direct \(dm/dh\), the bracket secant, \(G_{\rm stat}\), and
   the exact \(r(1+J_{0,\rm eff}\widetilde G_0)\). Resolve any sign disagreement
   before treating the continued samples as one differentiable component.
   **Completed:** direct \(dm/dh=-G_0^{\rm bare}\) holds to
   \(2.45\times10^{-8}\) relative, the two exact forms of \(F'\) agree to
   \(1.15\times10^{-16}\) absolute, and the narrow-section secants agree to
   \(1.24\times10^{-5}\) relative. The rejected
   \(r(1+J_{0,\rm eff}G_{\rm stat})\) substitution produced the false negative
   sign.
2. Repeat the continued factor-one branch on nested 129/257 node grids and a
   refined field ladder from 4.64 T toward \(B_c\). Record component topology,
   final residual, uniform and mesh margins, and the field at which independent
   lower seeds cease to agree. **Completed:** both resolutions and both direct
   lower-field seeds agree through 4.69 T. At 4.70 T only the 4.581818 T seed
   reaches the root; at 4.71 T both seeds fail refinement. Every accepted root
   has positive uniform, supremum, and mesh masses and an outer residual below
   \(7.1\times10^{-9}\).
   **Follow-up:** advancing each of two independent histories through smaller
   field steps removes the long-jump failure. Both paths converge to the same
   branch at 129 and 257 nodes through 4.7188 T, where
   \(D_{\rm uni}=1.12\times10^{-4}\), \(F'=8.97\times10^{-5}\), and the final
   residual remains below \(6.3\times10^{-9}\). Thus 4.70/4.71 did not locate
   a branch endpoint; they located the reach limit of the direct far seed.
3. Separate numerical reach from the area model. Map the root-support interval
   for every field and retain 0.85--0.95 only as diagnostic members unless an
   external or overlap calibration selects their missing area. Do not tune a
   factor merely to eliminate a plotted mask. **Completed for this boundary
   ladder:** factor-one support persists through 4.71 T on both grids even
   where refinement is seed-dependent or fails, so missing-area support does
   not explain the loss of numerical agreement above 4.69 T.
4. Only if steps 1--3 pass, prototype a distinct `ordered_boundary_retry` that
   targets frozen cold `phase==0` points with negative PM mass, uses at least
   two independently accepted lower-field cold states, forbids recovered
   states as seeds, and requires agreement plus explicit \(D_{\rm uni}\),
   \(F'\), component, and final-residual margins. It must never retry a cold
   paramagnetic (`phase==2`) point.
   **Completed:** the map-level retry uses two frozen lower cold states and
   never a recovered target, requires the next accepted cold state above to be
   PM, runs an independent PM solve with `mix_outer=0.25` and
   `max_outer=1000`, and admits only matching single-bracket components with
   \(D_{\rm uni},F'\ge10^{-3}\). The current 4.68 T factor-one candidate has
   \(D_{\rm uni}=0.023913\), \(F'=0.020014\), and PM mass
   \(-0.0178456\). Factors 0.75 and 1.5 remain incomplete.
5. Keep the present two-sided retry and production output unchanged until the
   boundary retry passes the focused field ladder and the full missing-area
   ensemble reports its incomplete members honestly.
   **Completed:** the focused map smoke recovers only the factor-one 4.68 T
   member, preserves 4.77 T as PM for every member, and exports the incomplete
   ensemble status. The old high-field smoke still masks 4.68 T when the new
   option is absent, and the 0.36/0.45 T adjacent-retry smoke is unchanged.

Stop and retain the mask if the derivative identity remains inconsistent, if
the accepted branch becomes seed-dependent, or if the stability margin tends
to zero faster than the requested field resolution can control.

**Disposition:** the derivative identity is consistent, and the fine-step
two-history 129/257 ladder supplies the new evidence that was previously
missing. The long-jump seed dependence is numerical reach, while the current
4.68 T production target is directly reproducible from two untouched lower
sources with comfortable margins. The distinct ordered-boundary retry is
enabled only in the opt-in missing-area driver and remains default-off in the
library. It does not turn the incomplete area ensemble into a certified
thermodynamic confidence band.

#### Priority 2: exact-zero-field limit

Treat 0 T as a separate endpoint problem, not as another missing-area sliver.
The present exact-zero call reaches an effectively zero splitting in
`invz_twolevel_ordered` (about \(1.78\times10^{-14}\) meV in the retained
audit), so the two-level representation becomes singular before an adjacent
field retry could be relevant. Do not copy the first positive-field spectrum
into the zero-field column or inject an arbitrary transverse field into a
result labeled \(B=0\).

After the high-field branch decision, proceed as follows:

1. Run a logarithmic and linear \(B\to0^+\) ladder, retaining the exact same
   interaction, missing-area, node, and acceptance settings. Record the
   two-level splitting, basis/projector diagnostics, component topology,
   self-energy, \(h_{\rm MF}\), uniform/mesh margins, and real-axis spectrum.
   **Completed:** a 19-point ladder from exact zero through 0.12 T is continuous
   in every basis-invariant single-ion observable and every accepted coupled
   quantity. The 129-node factor-one root stays near 0.038101 meV with final
   splitting about 0.41995 meV, positive masses, and residual below
   \(5.71\times10^{-9}\).
2. Compare the exact-zero construction with the positive-field limit using
   basis-invariant full-ion observables. Determine whether the discontinuity is
   only a degenerate-basis convention, a removable formula singularity, or a
   genuinely non-unique/soft limit. **Completed:** the electronic diagonal and
   off-diagonal \(J_z\) weights rotate strongly inside the degenerate doublet,
   but \(\tfrac12\operatorname{Tr}(PJ_zPJ_z)=30.372\), the combined static
   weight, the full-ion variance, and the full-ion path susceptibility are
   smooth. Exact-zero and \(10^{-6}\) T ensemble spectra agree within
   \(1.27\times10^{-12}\) absolute. The missing in-window peak is a censored
   zero-weight/window-edge limit, not a nonfinite spectrum.
3. If the limit is unique, derive and test explicit zero-splitting limiting
   expressions in the two-level response/vertex code. Validate them against
   successively smaller positive fields; do not use a numerical epsilon whose
   result depends materially on its chosen value. **Completed without an
   epsilon:** `domain_policy='return'` now marks the singular predictor
   `twolevel_domain_invalid` without mutating accepted state. No invalid node
   enters the missing-area quadrature. The exact \(M^2=0\) branch evaluates
   \(2m^2[\lambda_1-(1-n_{01}^2)K_0]/n_{01}^2\); positive-\(M^2\) arithmetic
   remains bitwise unchanged.
4. If no unique self-consistent limit is established, retain the exact-zero
   mask. A \(B\to0^+\) extrapolated or first-positive-field value may be offered
   only as a separately labeled plotting aid, never as a converged \(B=0\)
   production state.

The zero-field work is complete only when either an analytic/numerically stable
limit passes the ordinary final-state and susceptibility gates or the mask has
a precise symmetry-limit failure status. Visual smoothness alone is not an
acceptance criterion.

**Disposition:** the opt-in missing-area endpoint is complete. All three
members are ordered and finite at exact zero without adjacent retry, with
\(h_{\rm MF}=0.036749\)--0.038787 meV at 129 nodes. The 129/257 shifts are
\(5.51\times10^{-6}\)--\(1.86\times10^{-5}\) meV and full-spectrum differences
are below \(2.10\times10^{-11}\). Strict mode remains `node_failed` because its
complete-path predictor is explicitly outside the two-level domain. Thus the
approximation output is restored with honest provenance while the rigorous
acceptance contract is unchanged.

#### Execution order and production boundary

Both packets have now executed in that order. A follow-up high-field packet
resolved the earlier direct-seed stopping condition with fine-step two-history
continuation and wired the separately gated central-member retry. The
exact-zero packet passed its endpoint, ensemble, strict-mode, and 129/257
checks. Both exceptions are enabled only for the explicit missing-area
approximation; strict mode remains fail-closed.

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
1.568. This supports the tail-integrability prerequisite for Jensen's closed
model on the sampled window, but does not independently prove
\(\delta h(\infty)=0\); it is not transferable to the current hybrid.

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
