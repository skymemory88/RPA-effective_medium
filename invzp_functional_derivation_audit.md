# Common-functional derivation audit

**Recorded:** 2026-07-28

**Scope:** existing projected resummed/Jensen equations and their retained diagnostic residual

**Verdict:** partial generating structures exist, but the repository does not contain a common
thermodynamic functional or a thermodynamically binding selector for the coupled ordered roots

## 1. Question and standard

The immediate question is whether the current equations already imply a functional that was overlooked.
Such a functional would need to generate, from one declared approximation:

- the order-parameter equation;
- the self-energy and medium closure;
- the static and dynamic susceptibilities;
- the free energy and internal energy;
- a comparison of multiple stationary roots.

It is not sufficient that one scalar equation can be integrated after every other variable has been
solved. In appropriate conjugate variables and Matsubara measure, the one- and two-point derivatives
must be mutually consistent, and the stationary Hessian must generate the susceptibility.

## 2. Functional content that is already present

### 2.1 Branch-conditioned Jensen/Landau path potential

J 2.31–2.33 and `invz_hmf_ordered` define

\[
h_0(h)=\int_0^h r(h')\,dh',\qquad
F(h)=h_0(h)-J_0m(h).
\]

Using `dm/dh = -G0bare` gives

\[
\frac{dF}{dh}=r+J_0G0bare=crit.
\]

After a single-valued internal root branch has already been selected, one may define

\[
\Phi_{\rm path}(m)=\int_0^m F(m')\,dm',
\]

so the spontaneous solution is stationary and `crit(hstar)>0` is its local-minimum condition.
`docs/INVZ-DESIGN-RATIONALE.md` §A records this derivation; production evaluates the last
negative-to-nonnegative `F` crossing and reports `crit_star`.

This construction additionally requires `h -> m` to be single-valued and locally invertible on the
chosen segment (`dm/dh=-G0bare>0` in the present convention) and every compared path to share the same
physical `m=0` reference. A folded or disconnected arclength curve does not automatically define
`Phi_path(m)`.

This is useful but not the required common functional:

- `r(h)` already depends on a selected `Sigma/K0` root at every `h`;
- disconnected internal roots need not share a justified `Phi_path` normalization;
- production does not compare `Phi_path` across internal branches;
- J 2.34 normalization is not implemented as a binding thermodynamic comparison.

The published J 2.34 relation is an independent field-route/free-energy consistency check. The
repository's `invz_deltaF_ordered` is explicitly a partial finite-domain hybrid diagnostic rather
than `delta F(m=0)`, and the closed-model comparison has a stable approximately 13.7% route
difference. This is evidence against claiming that an exact common functional is already implemented,
but it is not a formal proof of inconsistency because the two diagnostic routes use different
approximations.

### 2.2 Isolated scalar EMT leaf potential

At fixed local propagator `G`, define

\[
\mathcal V_{\rm EMT}(K;G)=KG+
\left\langle\log\!\left[1+(J-K)G\right]\right\rangle_J .
\]

Then, away from a denominator zero,

\[
\frac{\partial\mathcal V_{\rm EMT}}{\partial K}
=G-\left\langle\frac{G}{1+(J-K)G}\right\rangle_J
=G-\langle G_q\rangle_J .
\]

Thus stationarity in `K` at fixed `G` reproduces the scalar q-average closure.

There is also a separate reduced one-dimensional primitive after the Dyson substitution has eliminated
`G` from the closure equation. At fixed `D,G0`, define

\[
\mathcal V_{\rm red}(K;D)
=\log(D+KG0)-K\,\overline G(D),
\qquad
\overline G(D)=
\left\langle\frac{G0}{D+JG0}\right\rangle_J .
\]

Then `partial_K V_red = G0/(D+KG0)-Gbar(D)`. Wolfram Language independently verified this
nontrivial first-derivative identity. Equality of the mixed `D,K` derivatives is only Clairaut
regularity of the displayed smooth scalar function away from its denominator zeros; it is not an
independent physics or functional-consistency check.

`V_red` is not obtained by simply substituting the `K`-dependent `G` into `V_EMT`; that substitution
also differentiates the implicit `G(K)`. Relating the fixed-`G` and reduced primitives requires the
missing Legendre/local and double-counting construction.

This leaf does not establish a thermodynamic functional for the implemented theory:

- the required local/Dyson, self-energy, source, and double-counting terms are absent;
- the ordered `Gstat` depends on `K0`, `Sigma0`, `lambda`, and the resummed elastic factor `xi`;
- the dynamic and static leaves feed common Matsubara moments;
- the real logarithm changes branch when `1+(J-K)G` changes sign and is singular at the observed
  denominator events;
- no vacuum term has been derived for the repository's full-electronuclear/two-level-`xi` hybrid.

The EMT leaf is therefore a useful derivation clue, not a root-selection functional.

### 2.3 Paramagnetic local self-energy primitive—and where the ordered approximation breaks it

The paramagnetic two-level self-energy map has more functional structure than the coupled Newton
residual reveals. Write

\[
p=\frac{M^2}{n_{01}^2},\quad
c=1-n_{01}^2,\quad
A=\frac12[g(0)+\beta c],\quad
\lambda_q=\frac1\beta\sum_n w_nK_ng_n^q .
\]

Then

\[
\Sigma_n=p\{\lambda_2-A\lambda_1+[\lambda_1-cK_n]g_n\}
\]

is generated by the quadratic primitive

\[
\Phi_{\Sigma}^{\rm PM}[K]
=p\beta\left(\lambda_1\lambda_2-\frac A2\lambda_1^2\right)
-\frac{pc}{2}\sum_n w_nK_n^2g_n^2 ,
\]

because

\[
\frac{\partial\Phi_{\Sigma}^{\rm PM}}{\partial K_n}
=w_ng_n\Sigma_n .
\]

Wolfram Language verified this identity symbolically and verified that
`diag(w_n g_n) dSigma/dK` is symmetric. This is a promising local vacuum/self-energy pair for the
strict paramagnetic WP0 derivation, although it still lacks the lattice, source, and double-counting
pieces of a complete free energy.

The current ordered correction does not preserve that canonical mixed-derivative relation. For two
generic nonzero-frequency slots `i,j`, with `a=m^2/n01^2`, symbolic differentiation gives

\[
w_ig_i\frac{\partial\Delta\Sigma_i}{\partial K_j}
-w_jg_j\frac{\partial\Delta\Sigma_j}{\partial K_i}
=-\frac{a\,w_iw_jg_ig_j(g_j-g_i)
\left[4(g_i+g_j)-g(0)\right]}{\beta g(0)},
\]

which is nonzero generically. More strongly, define

\[
M_{ij}=\frac{\partial\Sigma_i^{\rm ord}}{\partial K_j}
=\frac{w_jg_j}{\beta}\{\phi(g_j)+(p-2a)g_i\},
\]

\[
\phi(x)=(p-a)x-\frac{4a}{g_0}x^2+ag_0-pA.
\]

Any finite nonzero diagonal symmetrizer requires the three-frequency cycle identity
`M_ij*M_jk*M_ki=M_ji*M_kj*M_ik`. Its difference is

\[
\frac{w_iw_jw_kg_ig_jg_k}{\beta^3}
\frac{a(2a-p)}{g_0^2}
(g_i-g_j)(g_i-g_k)(g_j-g_k)\mathcal P_{ijk},
\]

\[
\mathcal P_{ijk}
=3ag_0^2-4Ag_0p+g_0^2p-4ag_0(g_i+g_j+g_k)
+16a(g_ig_j+g_ig_k+g_jg_k).
\]

There is no additional standalone `(p-a)` factor. This does not exclude a non-diagonal change to
natural conjugate variables. It does show that the ordered formula cannot simply be appended to the
paramagnetic primitive with different Matsubara weights. The most likely derivation target is the
frequency dependence omitted when J 2.26–2.27 replaces factors in the internal sums by their
zero-frequency values, together with the matching ordered vacuum and one-point terms.

## 3. Direct integrability audit of the current Newton residual

The retained diagnostic uses

\[
u=[\Sigma_1,\ldots,\Sigma_{n_\omega},K_0],\qquad
R=[\Sigma_{\rm map}-\Sigma,\;(K_0-J_{\rm loc})/J_{\rm scale}].
\]

Its analytic Jacobian has an asymmetric structural zero pattern. The independently supplied `K0`
replaces dynamic slot 1, while dynamic `K_n` for `n>1` depends only on `Sigma_n`. Consequently,
`lambda_1`, `lambda_2`, and `lambda_3` couple the higher-frequency variables into every
`Sigma_map` component, but there is generically no reciprocal `Sigma1` coupling into the corresponding
higher-frequency equations.

At the saved `T=1.00 K`, `B=2.9 T` positive-control root:

- `J(1,2) = -1.66493599e-3`, while `J(2,1) = 0`;
- `J(1,3) = -1.13864097e-4`, while `J(3,1) = 0`;
- `norm(J-J.',Inf)/norm(J,Inf) = 22.44897`;
- `norm(J-J.','fro')/norm(J,'fro') = 1.32017`.

This rules out treating the present residual as the gradient of a scalar functional in its current
coordinates. The asymmetric zero pattern also cannot be repaired by a common scalar factor or finite
nonzero diagonal row/variable scaling: at a root, derivatives of state-dependent diagonal factors
multiply `R` and vanish, leaving the zero-pattern obstruction.

This result is deliberately limited. It does not prove that no linked-cluster, Legendre, or 2PI
functional exists in natural variables such as the source, moment, and propagator. A nonlinear,
non-diagonal change to conjugate variables can change the Jacobian test. It proves only that
functional derivability cannot be inferred from the current solver residual.

## 4. What the source derivation actually supplies

Jensen derives J 2.31–2.33 by equating the differential single-site response to the bulk static
susceptibility and integrating the resulting one-dimensional field relation. J 2.34 is then presented
as an independent comparison between a field-area free-energy change and the internal-energy route.
The paper does not supply a common stationary functional for the self-consistent resummed effective
medium used here.

The repository goes further from a single declared functional by combining:

- Dyson-resummed dynamics;
- the ordered moment-form self-energy;
- a separately resummed static elastic factor containing `tanh`;
- full electronuclear static weights with a two-level `xi`;
- an independently closed scalar effective medium;
- a first-order/partial free-energy diagnostic.

Several of these are reasonable approximations individually. Their mixture is precisely why
functional derivability must be demonstrated rather than presumed.

### 4.1 McKenzie--Stamp auxiliary-field functional: a valid starting point, not a selector for the current roots

The newly supplied paper
`References/Thermodynamics of a quantum Ising system coupled to a spin bath.pdf`
(McKenzie and Stamp, *Phys. Rev. B* **97**, 214430 (2018)) provides a useful independent route.
Its status must be separated into three levels.

1. Equations (19)--(36) rewrite the partition function with a Hubbard--Stratonovich field and an
   all-orders local-cumulant action. In the paper's notation,

   \[
   Z_\phi=Z/Z_{\rm MF}=\int\mathcal D\phi\,e^{-\beta H_{\rm eff}[\phi]},
   \qquad
   H_{\rm eff}[\phi]
   =\frac1\beta\sum_{n\geq1}\frac{u_n}{n!}\phi^n ,
   \]

   with the frequency and momentum labels and conservation deltas implicit. Subject to the paper's
   preceding low-energy Hamiltonian truncation, this transformation and the untruncated cumulant
   series are exact. Equations (23)--(26) derive the one- and two-point spin responses from the same
   source-dependent partition function. This is precisely the common-generator property required
   for thermodynamic consistency.

2. Equation (37) truncates that action at quartic order. Its quadratic term gives RPA, while its
   `u3` and `u4` vertices are fixed by the local connected three- and four-spin cumulants. The paper
   then uses these interactions perturbatively in a high-density expansion; equation (59) is the
   leading `O(1/z_c)` magnetization correction. The schematic double-well discussion in Fig. 6 is
   tied to that quartic/high-density framework. In particular, the reported `u2,u4>0` result is a
   zero-temperature result for the models evaluated there, not a general stability theorem for the
   present finite-temperature electronuclear closure.

3. The paper does **not** derive the repository's self-consistent resummed
   `Sigma <-> EMT <-> K0` equations by varying this action, does not construct a 2PI functional of a
   dressed propagator, and does not evaluate a common free energy on the multiple
   `[\Sigma,K0]` roots found here. `H_eff[phi]` is the auxiliary-field action inside a functional
   integral; away from a declared saddle approximation, comparing its raw value at two field
   configurations is not the same as comparing thermodynamic free energies. The physical one-point
   functional would be the source Legendre transform, and a stationary dressed-propagator theory
   requires the corresponding bilocal/2PI transform and its double-counting terms.

Consequently, the paper is a genuine derivational clue but not a drop-in branch selector. It supports
the existing WP0/WP4 direction:

- use one source-dependent local Hamiltonian to generate `C1--C4`;
- adopt the Hubbard--Stratonovich action as the normalization/sign starting point;
- derive either a strict vacuum/free-energy expansion without re-stationarizing it, or a complete
  1PI/2PI Legendre functional whose variations generate both the self-energy and one-point equation;
- compare stationary values only after the trace-log, skeleton, source, and double-counting pieces
  have been derived at one declared order.

The smallest new check suggested by the paper is therefore not another production solver. It is a
normalization cross-audit: translate equations (34)--(37) into the WP0 conventions and verify that
the already exact local `C2/C3/C4` oracles reproduce the paper's quadratic, cubic, and quartic
vertices, including Matsubara weights and the static degenerate limits. Passing that check would
anchor the starting action. It would not by itself repair the rejected varied-covariance WP4
candidate; the exact local bilocal Legendre kernel, same-order `C3-C3` sector, source derivatives,
and electronuclear Matsubara tail would remain necessary.

## 5. Existing branch-selection mechanics are not thermodynamic rules

The current code:

- accepts whichever resummed static root the Picard seed reaches;
- warm-starts the ascending `h` sweep and cold-retries failed nodes;
- selects the last sampled increasing crossing of the stitched `F(h)` profile.

These rules make a run deterministic when the local basin is simple. They do not compare the
thermodynamic weights of multiple `Sigma/K0` roots, and the new fold diagnostics show that basin
selection changes the Jensen path itself.

The smooth-`r(h)` backup is documented separately in `biased_convergence_solution.md`. It is an
explicit bias, not a latent rule recovered by this audit.

## 6. Minimal credible derivation

The smallest route that can establish a common functional is:

1. Freeze a scalar two-level Hamiltonian with a finite no-self-site coupling matrix, source `h`,
   centred operator `delta J`, the repository's `G=-chi` sign, and explicit high-density scaling.
2. Construct exact local `C1`–`C4` Matsubara cumulants, including the static/degenerate limits.
3. Enumerate the complete strict retained-order vacuum diagrams with signs and symmetry factors.
4. Differentiate the same vacuum object to obtain both:
   - the two-point self-energy;
   - the one-point fluctuation field vertex.
5. Verify the Maxwell pair, schematically

   \[
   \frac{\delta\Sigma_n}{\delta m}
   =
   \frac{\delta H_{\rm corr}}{\delta G_n},
   \]

   and weighted Hessian symmetry in the declared conjugate variables.
6. Either derive the EMT log/determinant contribution, its local term, and double counting so its
   variations give both Dyson and closure equations, or omit EMT initially and use the direct finite
   lattice propagator.
7. Evaluate the functional at every stationary solution and verify
   `m=-dF/dh`, `chi=-d2F/dh2`, and `U=d(beta F)/d beta`.
8. Compare the retained weak-coupling coefficients with a two- or four-site exact-diagonalization
   oracle.

Only after this strict construction closes should dressed lines or the static elastic resummation be
introduced through stationarity of the same functional.

## 7. Decision and next work package

The focused audit does **not** find an already implied physical branch selector. It supports the
existing staged recommendation:

1. retain the smooth-`r(h)` prescription only as a backup;
2. start WP0 on a closed scalar two-level model;
3. derive vacuum, one-point, and two-point terms before writing a new solver;
4. keep production Jensen code unchanged during that derivation.

The strict unicyclic WP0 scope is now closed in `invzp_functional_wp0_spec.md` and
`invzp_functional_wp0_ring_derivation.md`. The complete retained vacuum is the `C2` ring; source
derivatives supply the ordered `C3` and `C4+C3^2` inventory, and exact two-site gates pass. The
isolated `invz_functional/` prototype has therefore begun. The next theoretical scope,
`O(1/z^2)`, remains blocked from implementation until all `C3-C3` bicyclic classes, the one-`C4`
class, and their `C5/C6` source derivatives are derived together.
