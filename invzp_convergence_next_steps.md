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

This work package is the immediate next implementation step.

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
single-ion features such as the node-18 electronuclear peak.

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
