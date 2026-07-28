# WP4 specification — stationary nonlocal-return skeleton

**Recorded:** 2026-07-28

**Status:** derivation/pre-registration; no implementation or production authorization

**Trigger:** the strict electronuclear ring functional has stable ordered states at
`Bx=1.5--4.6 T` and a stable symmetric state at `4.9 T`, but no admissible stationary state at
`4.7--4.85 T`. The missing interval terminates on the Gaussian determinant, not on a numerical
root-solver failure.

## 1. Scope and stop rules

WP4 may add a self-consistent skeleton only if it remains stationary under one declared functional.
It must not:

- add a floor, broadening, principal-value rule, or absolute value to a Matsubara determinant;
- dress only `omega=0`;
- retain the one-`C4` core while omitting the same-order `C3-C3` core at nonzero source;
- use the production `tanh/xi` static replacement;
- infer the `C3-C3` term from the coupling-eigenvalue multiset alone;
- dispatch from a production spectrum or phase solver before every gate in section 8 passes.

The first implementation remains scalar and translation invariant, but it needs q-resolved
sublattice coupling matrices and transforms, not only sorted eigenvalues.

## 2. Local generator and one-particle-irreducible vertices

At fixed total local source `h`, let `W0` be the exact dimensionless connected local generator and
let `C2`, `C3`, and `C4` be its connected, centred Matsubara derivatives. Compound indices below
include imaginary time/frequency and any retained local component. Repeated compound indices are
integrated/summed.

The local Legendre transform gives

\[
\gamma^{(2)}=C_2^{-1},
\]

\[
\gamma^{(3)}_{abc}
=-(C_2^{-1})_{aa'}(C_2^{-1})_{bb'}(C_2^{-1})_{cc'}C_{3,a'b'c'},
\]

\[
\begin{aligned}
\gamma^{(4)}_{abcd}
={}&-(C_2^{-1})_{aa'}(C_2^{-1})_{bb'}(C_2^{-1})_{cc'}
    (C_2^{-1})_{dd'}C_{4,a'b'c'd'}\\
&+\gamma^{(3)}_{abe}C_{2,ef}\gamma^{(3)}_{fcd}
 +\gamma^{(3)}_{ace}C_{2,ef}\gamma^{(3)}_{fbd}
 +\gamma^{(3)}_{ade}C_{2,ef}\gamma^{(3)}_{fbc}.
\end{aligned}
\]

For one scalar static variable these reduce to

\[
\gamma_3=-\frac{C_3}{C_2^3},\qquad
\gamma_4=\frac{3C_3^2-C_2C_4}{C_2^5}.
\]

Wolfram Language independently inverted the local series
`W(j)=C2*j^2/2+C3*j^3/6+C4*j^4/24` and reproduced both formulas. At a `Z2`-symmetric source,
`C3=gamma3=0`.

The multilevel service must build these vertices with exact KMS/Hermite limits and a rank ladder.
Finite source differences are acceptable as an oracle, not as the final dense `C4` engine.

## 3. Full propagator and nonlocal return

Let

\[
D^{-1}_{0}=C_2^{-1}-J
\]

in site-time compound space, and vary a positive Matsubara covariance `D` independently. Define the
nonlocal return

\[
\mathcal R_{ij}=D_{ij}-\delta_{ij}C_2 .
\]

At `J=0`, stationarity gives `D=C2` and hence `R=0`; this is not an off-shell identity because `D`
is an independent variational argument. In the stationary weak-coupling expansion, `R_ij=D_ij` for
`i != j` begins with one intersite bond, while `R_ii` begins with a two-bond return. This subtraction
prevents the exact local cumulants from being dressed again by purely local diagrams on shell.

The distinction is essential. A local scalar
`Rloc=<D_q>-C2` is enough for the one-`C4` double-return core, but not for two `C3` vertices joined by
three off-site lines. The latter depends on `R_ij^3`, or equivalently a momentum convolution with
sublattice eigenvectors. A `Jnu_flat` eigenvalue multiset contains insufficient information.

## 4. Stationary trace term

Per site, define

\[
\begin{aligned}
f_{\rm tr}[D;h]
=\frac{1}{2\beta N}\operatorname{Tr}\big[
&\log(C_2D^{-1})
+(C_2^{-1}-J)D-I+JC_2
\big].
\end{aligned}
\]

The trace covers site and Matsubara indices. The linear `JC2` subtraction is retained even when an
exact no-self-site lattice makes its trace vanish.

With no skeleton term, stationarity gives

\[
D=(C_2^{-1}-J)^{-1}
\]

and substitution gives exactly

\[
\frac{1}{2\beta N}\operatorname{Tr}
\{\log(I-JC_2)+JC_2\},
\]

the strict ring functional. At `J=0`, `D=C2`, `R=0`, and `f_tr=0`.

## 5. Smallest closed skeleton

The first nonlocal 2PI cores are

\[
\Phi_4[\mathcal R]
=\frac18\,
\gamma^{(4)}_{abcd}
\mathcal R_{ab}\mathcal R_{cd},
\]

\[
\Phi_{33}[\mathcal R]
=-\frac1{12}\,
\gamma^{(3)}_{abc}\gamma^{(3)}_{def}
\mathcal R_{ad}\mathcal R_{be}\mathcal R_{cf}.
\]

All local vertex site labels are contracted as required; the return lines carry the full site-time
structure. Fix the transform by
`X(tau)=beta^(-1) sum_n X_n exp(-i*omega_n*tau)`, so a frequency-space line is
`<X_n X_-n>=beta*D_n`, a cubic vertex carries `beta^-2`, and a quartic vertex carries `beta^-3`.
In the resulting energy-per-site signed-frequency convention, the scalar quartic term has the form

\[
\phi_4
=\frac{1}{8\beta^2}\sum_{n,m}
\gamma_4(n,-n,m,-m)\,
\mathcal R_{ii}(n)\mathcal R_{ii}(m),
\]

while the cubic term is

\[
\begin{aligned}
\phi_{33}
=-\frac{1}{12\beta^2N}
\sum_{i,j}\sum_{n,m}
&\gamma_3(n,m,-n-m)\gamma_3(-n,-m,n+m)\\
&\times
\mathcal R_{ij}(n)\mathcal R_{ij}(m)
\mathcal R_{ji}(n+m).
\end{aligned}
\]

The signed grid is mandatory for the convolution. Positive-frequency weights may be introduced only
after an exact reduction of the convolution and its Nyquist/cutoff boundary terms.

The factors follow directly from the interaction normalization
`gamma3*phi^3/3! + gamma4*phi^4/4!`: the quartic double bubble has `3/4!=1/8`; the two-cubic vacuum
has `(1/2)(1/3!)^2*3!=1/12`, and its sign in `-log Z` is negative. Wolfram Language independently
confirmed these combinatorial values.

The topology statement must be made only after re-expanding the 1PI vertices in connected local
cumulants. The `-C4` part of `gamma4` generates the one-`C4`/two-return topology, but away from a
`Z2`-symmetric source the three `gamma3*C2*gamma3` channels inside `gamma4` also carry `C3-C3`
content. The explicit cubic core supplies the leading triple-bond `C3-C3` topology, and dressing its
three return lines generates the same 2PI core with allowed ring decorations. Therefore neither
`Phi4` nor `Phi33` may be graded in isolation at nonzero source: their sum, after the 1PI-to-connected
conversion, must reproduce the complete declared `C4` plus `C3-C3` weak-coupling content.

## 6. Common functional and stationarity

The convention is fixed most compactly for the dimensionless extensive functional:

\[
\begin{aligned}
\beta N f_{\rm 2PI}(m,h,D;H)
={}&\beta N\left[f_0(h)+(h-H)m-\frac12J_0m^2\right]\\
&+\frac12\operatorname{Tr}\big[
\log(C_2D^{-1})+(C_2^{-1}-J)D-I+JC_2
\big]\\
&+\Phi_4[\mathcal R]+\Phi_{33}[\mathcal R].
\end{aligned}
\]

Thus the lowercase `phi4` and `phi33` displayed in section 5 are exactly
`Phi4/(beta*N)` and `Phi33/(beta*N)`. Their explicit `1/beta^2` factors follow from the stated
Matsubara transform and are part of the immutable normalization gate. No convention may
be adjusted after seeing the LiHoF4 field scan.

Stationarity in `m` gives

\[
h=H+J_0m.
\]

Stationarity in `h`, including the explicit source dependence of `C2`, `gamma3`, and `gamma4`, gives
the physical moment equation. It must be evaluated as a derivative of the same functional; dropping
vertex-source derivatives would omit the required `C5/C6` content.

Stationarity in `D` gives

\[
D^{-1}=C_2^{-1}-J+\Sigma_{\rm skel},\qquad
\Sigma_{\rm skel}=2\,\frac{\delta(\Phi_4+\Phi_{33})}{\delta D}
\]

in compound-index normalization. `Phi4` produces a local frequency-dependent mass. `Phi33`
generically produces momentum and sublattice dependence; replacing it by a local self-energy is a
new approximation and is outside this specification.

This is a skeleton resummation: its first weak-coupling terms are the complete declared
`O(1/z^2)` 2PI cores, while its stationary solution also resums higher orders. Results must therefore
be labelled `skeleton`, never strict `O(1/z^2)`.

## 7. Numerical representation

The first implementation should use:

1. a signed Matsubara grid with an explicit convolution-tail budget;
2. q-resolved Hermitian sublattice coupling matrices with frozen normalization and Gamma policy;
3. positive-domain variables for `D^{-1}`; leaving the domain is a reported outcome, not clipped;
4. a safeguarded Newton/trust-region solve of the complete stationarity residual;
5. independent cold starts, source continuation in both directions, and functional comparison of
   every stationary solution;
6. local-rank ladders for `C2`, `gamma3`, and `gamma4`, with separate discarded-weight estimates.

The existing production eigenvalue cache may still serve `Phi4` and the trace term, but `Phi33`
needs matrix/eigenvector provenance or an equivalent real-space transform.

## 8. Entry and exit gates

Implementation is authorized only after gates 1–4 have immutable fixtures:

1. **Local Legendre gate:** analytic two-level `gamma3/gamma4` agrees with numerical Legendre
   inversion, including `h=0`, nonzero source, and a near-degenerate fixture.
2. **Normalization gate:** `J=0` gives `D=C2`, `R=0`, the exact local free energy, and zero skeleton
   correction at every source.
3. **Ring-reduction gate:** setting `Phi4=Phi33=0` reproduces the existing ring value, gradient, and
   Hessian bit-for-bit.
4. **Weak-coupling topology gates:** the three-site mixed coefficient in
   `invzp_functional_wp0_spec.md` reproduces `+0.60239286138607826706` for its connected-`C4`
   residual. A nonzero-source cluster must first isolate the leading triple-bond `C3-C3` coefficient
   and then, at the next power where the `C3*C2*C3` part of `gamma4` contributes, verify the total
   `C3-C3` content of `Phi4+Phi33` after the 1PI-to-connected conversion. The latter total—not the
   explicit `Phi33` term alone—is the sign/routing/symmetry-factor oracle.
5. **Order gate:** exact-cluster minus skeleton has the first omitted power predicted by a fresh
   graph inventory; a fitted power alone cannot substitute for the inventory.
6. **Thermodynamic gate:** at every stationary root,
   `m=-dF/dH`, `chi=-d2F/dH2`, and `U=d(beta F)/d beta` agree with independent reoptimized
   differences.
7. **Domain/solver gate:** no denominator floor is used; all accepted roots have a declared positive
   covariance domain and converge from independent starts.
8. **Discretization gate:** local rank, Matsubara cutoff/convolution padding, q grid/offset, and
   dipolar backend errors are reported separately.
9. **LiHoF4 discriminator:** on the frozen scalar coupling input, stationary states form a continuous
   functional-minimum branch through `Bx=4.6--4.9 T`; otherwise WP4 does not fix the strict-ring
   obstruction.
10. **Production stop:** real-axis continuation and `invz_run_spectra` dispatch remain prohibited
    until gates 1–9 pass.

If the complete skeleton still has no stationary branch through the transition, the result is a
theory failure of this candidate. The backup smooth-`r(h)` prescription remains documented in
`biased_convergence_solution.md`; it must not be silently substituted into WP4.

## 9. Gate status at recording

- The scalar static Legendre formulas were checked by Wolfram series inversion.
- `invzf_local_1pi_static` agrees with an independently inverted nonzero-source two-level moment map
  to `2.7e-15`, `6.2e-11`, and `1.1e-8` for `gamma2`, `gamma3`, and `gamma4`.
- The near-degenerate symmetric gate initially failed because of cancellation in
  `invzf_twolevel_local`; its analytic Hermite series is now used and gives
  `gamma4=2/beta` to `5.6e-17`.
- Exact frequency-labelled two-level `C2/C3/C4` now use ordered-simplex matrix exponentials.
  They reproduce `C3(n,-n,0)=dC2(n)/dh` and `C4(n,-n,0,0)=d2C2(n)/dh2` below `3e-15`;
  the closed zero-source `C4(n,-n,m,-m)` formula, including its anomalous
  `beta*(delta_mn+delta_m,-n)` term, closes below `2e-15`.
- The ring-reduction implementation already exists, but has not yet been embedded in a varied-`D`
  skeleton object.
- The zero-source three-site `C4` oracle and the leading centred-pair nonzero-source `C3-C3` oracle
  exist.  The nonzero-source centred-chain mixed-`C4` cancellation oracle also exists, and the
  q-resolved coupling-matrix input is exposed for both dipolar backends.

The immutable two-level fixtures required by gates 1--4 are now present.  This authorizes an isolated
smallest-skeleton implementation, but not a LiHoF4 run: the full electronuclear frequency-labelled
vertex engine and the remaining exit gates are still absent.

## 10. Immutable centred-pair `C3-C3` oracle

The leading nonzero-source cubic core now has a contamination-free exact-cluster fixture.  Define

\[
\mathcal H_c(j)=
\sum_{i=1}^2[\Delta n_i-hX_i]
-j(X_1-m_0)(X_2-m_0),\qquad
m_0=\langle X\rangle_{j=0,h},
\]

and hold `m0` fixed throughout the signed `j` scan.  Since every interaction leg is the centred
operator `delta X=X-m0`, the linear coefficient vanishes.  The cubic coefficient is

\[
[j^3]F_c
=-\frac{1}{6\beta}
\int_0^\beta d\tau_1d\tau_2d\tau_3\,
C_3(\tau_1,\tau_2,\tau_3)^2 .
\]

It contains only the two-site, three-parallel-bond `C3-C3` core: a `C2` ring has even degree, a
`C4` vertex requires four legs, and all `C1` attachments vanish by centring.  In code the exact same
Hamiltonian is evaluated without adding a second diagonalizer:

```matlab
pair = invzf_two_site_exact(Delta,M,h-j*m0,j,beta,wn);
Fc = pair.F-j*m0^2;
```

For the frozen fixture

```text
Delta = 1.3, M = 1, beta = 1.7, h = 0.37,
m0 = 0.42257074773559524799819924078642109379,
C2(0) = 0.97505788500918161089657627947292331138,
C3(0,0,0) = dC2(0)/dh
           = -1.18442977486108158592225349632114137689,
```

120-digit matrix-exponential diagonalization and two levels of even-power Richardson cancellation
give

\[
[j^2]F_c=-0.2892028519514694014\ldots ,
\]

\[
\boxed{[j^3]F_c=-0.082768355766288255\ldots } .
\]

The cubic extraction uses

\[
d_3(s)=
\frac{F_c(2s)-2F_c(s)+2F_c(-s)-F_c(-2s)}{12s^3}
=[j^3]F_c+O(s^2),
\]

followed by Richardson factors `4/3` and `16/15`.  Results from
`s={1/20,1/40,1/80,1/160}` agree in the displayed digits.  This number grades the sign, routing, and
`1/12` two-cubic symmetry factor of the leading `Phi33` re-expansion.

This oracle does not grade the dressed same-order bicyclic continuations or the later
`C3*C2*C3` contribution carried by `gamma4`; the combined nonzero-source gate in section 8 remains
open until a higher-power cluster coefficient tests those terms together.

## 11. Frozen q-resolved coupling contract

The existing scalar production input flattened four sorted eigenvalues per q point and discarded
the matrices, which is insufficient for `Phi33`.  The fourth output of `invz_jq_modes` now captures
the exact Hermitian `Jcc(:,:,iq)` page immediately before that unchanged eigenvalue call.  It is
strictly `nargout`-gated:

- one- through three-output calls retain the existing cache lookup and allocate no matrix pages;
- a fourth-output request recomputes, because the frozen cache contains eigenvalues but no matrices;
- the brute-force and Ewald backends use their existing Lorentz/Gamma conventions unchanged;
- the opt-in ODD path fails explicitly for a matrix-page request rather than returning incomplete
  data.

The fourth output of `invz_bz_couplings` packages

```text
qvec, normalized q-row weights, Jcc pages,
Jnu_unflat [nq x 4], Jnu_flat, Juni, and info.
```

The flattening contract remains
`flat=(branch-1)*nq+q`.  Legacy endpoint-inclusive rows, including periodic duplicates, retain equal
literal-row weight after the historical coordinate-zero Gamma removal.  The phase-grid route uses
the weights and `invz_is_gamma_equiv` policy returned by `invz_phase1_qgrid`.  No eigenvectors are
stored: downstream diagonalization must record arbitrary phase and degenerate-subspace gauge, while
matrix-basis convolutions should avoid that gauge entirely when possible.

Fresh two-q gates found zero Hermiticity residual and reproduced the existing sorted eigenvalues
below `1e-14` for both brute-force and Ewald backends.  Repeated ordinary versus detailed calls gave
bitwise-identical `Jnu`, `info`, `Juni/Jaa0` outputs.  This closes the input-availability part of the
WP4 entry gate without enabling a skeleton or changing a production default.

## 12. Nonzero-source mixed-`C4` cancellation oracle

The higher-power nonzero-source gate uses a centred three-site open chain,

\[
\mathcal H_c(a,b)
=\sum_{i=1}^3[\Delta n_i-hX_i]
-a\,\delta X_1\delta X_2
-b\,\delta X_2\delta X_3,\qquad
\delta X=X-m_0 .
\]

The monomial `[a^2 b^2]` has local degrees `(2,4,2)`.  Its exact connected graph inventory is
therefore the fourth-order Gaussian ring plus one connected `C4` at site 2.  No `C1` attachment
survives the centring and no pair of degree-three vertices can match this monomial.  This makes the
fixture decisive away from symmetry: after converting connected vertices to 1PI vertices, the
`C3*C2*C3` channels inside `gamma4` and the compensating stationary `Phi33`/trace contributions must
cancel in the total re-expansion.  Grading `Phi4` alone would fail this oracle.

At the same frozen local point as section 10,

```text
Delta=1.3, M=1, beta=1.7, h=0.37,
```

120-digit `8x8` matrix-exponential evaluation and two Richardson levels give

\[
[a^2b^2]F_{\rm exact}
=-0.16480453047308016228870453332302872015\ldots .
\]

For a scalar local `C2`, `Tr J^4` has mixed coefficient `4a^2b^2`, so the ring contribution is

\[
[a^2b^2]F_{\rm ring}
=-\frac{1}{2\beta}\sum_{n=-\infty}^{\infty}C_2(i\omega_n)^4
=-0.26598229624005026010142443376005660462\ldots .
\]

The immutable connected-`C4` residual is consequently

\[
\boxed{
[a^2b^2](F_{\rm exact}-F_{\rm ring})
=+0.10117776576697009781271990043702788447\ldots } .
\]

An independent direct contraction of the exact dynamic vertex,

\[
-\frac{1}{4\beta^2}\sum_{n,m}
C_4(n,-n,m,-m)C_2(n)C_2(m),
\]

converges to the same residual.  A fresh double-precision centred-cluster scan gives the total
coefficient within `6.5e-10`; the signed Matsubara ring sum through `|n|=320` agrees with the displayed
ring value within `1e-15`.
