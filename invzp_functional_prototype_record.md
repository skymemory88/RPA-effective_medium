# Strict scalar functional prototype record

**Recorded:** 2026-07-28

**Theory:** `invzp_functional_wp0_ring_derivation.md`

**Scope:** isolated WP1/WP2/WP3 scalar fixture; no production dispatch

## 1. Implemented surface

`invz_functional/` now contains:

- `invzf_twolevel_local`: exact source-biased two-level `f0`, moment, static susceptibility,
  positive connected `C2(iwn)`, `dC2/dh`, and `d2C2/dh2`, including the analytic elastic zero mode;
- `invzf_ring_scalar`: finite-mode ring value and analytic first/second source derivatives, signed or
  nonnegative Matsubara measures, explicit `J^2` coefficients, fail-closed pole status, and a
  rigorous scalar free-energy tail bound;
- `invzf_cutoff_convergence`: successive-doubling/Richardson audits of `f`, `U`, the gradient, and
  Hessian;
- `invzf_scalar_functional`: the independently varied `(m,h)` functional, gradient, and Hessian;
- `invzf_stationary_scalar`: bounded root enumeration, reduced-curvature classification, and
  functional-value comparison without breaking degenerate minima;
- `invzf_cluster_exact`: dense exact one- through eight-site scalar clusters with arbitrary symmetric
  zero-diagonal coupling matrices, scalar or site-dependent sources, thermodynamics, and stable
  KMS/Hermite Lehmann response;
- `invzf_two_site_exact`: the pair-specialized wrapper used by the Jensen dynamic gate;
- `invzf_centered_pair_exact`: the fixed-reference centred-bond pair whose cubic coefficient isolates
  the leading nonzero-source `C3-C3` topology;
- `invzf_centered_cluster_exact`: arbitrary fixed-reference centred finite clusters, including the
  nonzero-source mixed-`C4` cancellation oracle;
- `invzf_projected_inputs`: a read-only bridge from the production transverse doublet, derived
  `Jcc0/Jaa0`, and scalar BZ spectrum into the isolated pilot;
- `invzf_stationary_convergence`: complete root re-enumeration across Matsubara cutoffs;
- `invzf_mode_grid_audit`: complete root re-enumeration across production BZ grids, with coupling
  provenance and a separate pole-domain gate;
- `invzf_electronuclear_local`: full source-biased 136-state local free energy and stable connected
  Lehmann response, with nested source/beta derivative stencils;
- `invzf_electronuclear_inputs`: safe separation of applied `Bz` from the total variational local
  source used by the full local oracle;
- `invzf_scalar_functional_local`: common-functional assembly from either analytic two-level or
  precomputed multilevel local data;
- `invzf_local_1pi_static`: the first static connected-to-1PI amputation gate for WP4;
- `invzf_twolevel_cumulant`: exact frequency-labelled two-level `C2/C3/C4` using analytic
  ordered-simplex matrix exponentials;
- `invzf_twolevel_1pi_vertex`: dynamic scalar amputation including all three
  `gamma3*C2*gamma3` quartic channels;
- `invzf_twolevel_vertex_table`: exact signed-grid freezing for future convolution or local-bilocal
  derivation gates;
- `invzf_local_bilocal_hessian_static`: exact `Q=0` local response Hessian from the connected
  zero-mode `C2/C3/C4` block, exposing both fixed-source and fixed-moment curvatures;
- `invzf_local_bilocal_hessian_modes`: the corresponding finite nonnegative-frequency response and
  bilocal Hessians, with explicit zero/nonzero-mode pairing multiplicities;
- `invzf_multilevel_cumulant`: `C2/C3/C4` exact for an explicitly retained multilevel local space,
  using block exponentials and a wide-spectrum dense-generator exponential-action fallback.

No production spectrum or phase solver dispatches into this directory.  The local exact/ring
fixtures are self-contained; the explicitly named production-input adapters call existing
`invz_common`/`invz_projected` services read-only.

## 2. Fresh verification

A temporary test driver was generated outside the worktree and is not retained. MATLAB R2025a
reported:

```text
local derivatives: d1 2.969e-13 d2 3.828e-08 db 9.268e-13
ring derivatives: d1 9.318e-16 d2 9.162e-12
common functional: grad 7.735e-11 Hessian 8.109e-11
internal-energy envelope: 6.032e-13
field thermodynamics: m 6.098e-13 chiF 1.018e-08 chim 2.436e-09
cutoff audit: actual 4.450e-11 bound 4.503e-11 order 2.963
stationary/domain gates: roots 3/3 signed 1.735e-18 pole -9.732e-01
edge stability: KMS 0.000e+00 weak 2.220e-16 nonfinite/symmetry PASS
two-site source oracle: series [ 8.67e-09 -3.89e-08 -9.67e-09], functional 1.03e-08
two-site internal energy: 1.748e-08
Jensen dynamic j^2 gate: 7.732e-11
cluster scaling H=0.00: p = 4.003, 4.001, 4.001
cluster scaling H=0.25: p = 2.844, 2.889, 2.922
three-site mixed C4: exact -0.0804438366493 ring -0.682836698356 C4 +0.602392861707
C4 beta-delta structure: 1.110e-16
INVZF_WP0_TESTS PASS
```

The checks cover:

1. `C2(iw0)=dm/dh`;
2. analytic `dC2/dh` and `d2C2/dh2` against five-point source differences;
3. ring first/second derivatives against independent differences of the ring value;
4. the common-functional gradient and Hessian against independent numerical derivatives;
5. `U=d(beta F)/d beta` by an independent derivative of reoptimized stationary free energies;
6. `m=-dF/dH` and `chi=-d2F/dH2=dm/dH` from independently reoptimized field differences;
7. the analytic free-energy tail bound against a `n=4096` reference and `N^-3` Richardson order;
8. equality of signed and positive-frequency Matsubara evaluations;
9. fail-closed behavior at a ring pole;
10. stable KMS/Hermite limits at high temperature and near degeneracy, preservation of the ring
   `J^2` term down to `J=1e-16`, nonfinite fail-closed states, and rejection of asymmetric signed
   frequency data;
11. three stationary roots below the scalar critical point: one unstable symmetric root and two
   exactly degenerate stable minima selected together;
12. unique lower functional value after a small symmetry-breaking field;
13. exact nonzero-source pair `F`, `m`, `chi`, and `U` and their `O(j^3)` residual after the retained
   `j^2` series;
14. the exact local Matsubara `j^2` coefficient against the expanded Jensen PM expression for
    `beta*Delta={0.5,1.7,5}` and `n=0,1,2`;
15. a four-site ring cluster whose exact-minus-functional free-energy residual scales as `epsilon^4`
    at the symmetric point and approaches `epsilon^3` at nonzero source, exactly matching the first
    omitted `C4` and `C3-C3` classes;
16. the three-site mixed-bond `C4` oracle and the anomalous `beta*delta_nm` component of the exact
    two-level connected `C4`.
17. the full electronuclear local `C2` against `invz_chi0z`, its local Maxwell identities, source-step
    halving, `Z2` source symmetry, and generic-functional parity with the analytic two-level route.
18. the exact-symmetry/near-degenerate static 1PI limit; this gate exposed and repaired a cancellation
    in the original two-level `d2C2/dh2` evaluation, after which `gamma4 -> 2/beta` to `5.6e-17`.
19. the centred nonzero-source pair coefficient
    `[j^3]F=-0.082768355766288255...`; a fresh double-precision signed-scan/Richardson extraction
    agrees within `1.6e-10`.
20. q-resolved `Jcc` pages reproduce the pre-existing sorted eigenvalues below `1e-14` for both
    brute-force and Ewald backends, are exactly Hermitian in the fixtures, and leave ordinary
    one- through three-output calls bitwise unchanged.
21. dynamic `C2` against the analytic local helper, `C3(n,-n,0)=dC2/dh`,
    `C4(n,-n,0,0)=d2C2/dh2`, permutation symmetry, the near-degenerate `gamma4=2/beta` limit, and the
    exact zero-source beta-delta `C4` formula; the largest displayed error is `5.4e-15`.
22. the centred nonzero-source chain coefficient
    `[a^2b^2]F=-0.164804530473080162...`, its ring part
    `-0.265982296240050260...`, and its connected-`C4` residual
    `+0.101177765766970098...`.
23. the exact static local bilocal Hessian.  For
    `Delta=1.3`, `M=1`, `beta=1.7`, and zero source, it gives
    `C2=1.23428900560563`, `C4=-4.46287005289486`, and reproduces the exact
    mixed-chain coefficient `-0.0944822130558562` with error `1.11e-16`.
    At source `h=0.37`, the connected response block remains positive
    (`lambda_min=0.1775`) and its fixed-moment Schur denominator is
    `0.621968104727707`.  The coordinate distinction is decisive there:
    the fixed-source partial-Legendre curvature reproduces the restricted static connected-`C4`
    coefficient exactly, while using the fixed-moment Schur curvature in the still-unreduced
    `(m,h,D)` functional adds the forbidden
    `C2*C3^2/(4*beta^2)=0.118329000835305`.
24. the finite-mode local bilocal response for labels `0:3` at the same nonzero source.  Subtracting
    the direct `A_n A_m` product from the exact full four-point function reproduces the assembled
    bilinear covariance within `2.22e-16`, including two crossed pairings at `n=m=0`, one at
    `n=m>0`, and none at unequal representatives.  An independent symmetric source difference
    reproduces its `C3(0,n,-n)` row within `1.92e-10`; the block-inverse and Schur-inverse fixed-m
    Hessians agree to `2.1e-17` relatively.  Re-stationarizing the leading three-site return using
    the fixed-source curvature reproduces the finite-cutoff dynamic connected-`C4` residual within
    `1.7e-15` at zero source and `3.4e-16` at `h=0.37`.  At the latter point the cutoff-3 ring and
    residual are `-0.265982292019735` and `+0.101148371627387`; the remaining difference from the
    frozen infinite-frequency residual is a cutoff tail, not a normalization error.
25. the multilevel block-cumulant engine.  Across rank-2/3/4 zero- and nonzero-frequency two-level
    fixtures, including a near-degenerate case, its connected values agree with the independently
    enumerated state-path oracle within `1.8e-15`.  For the full 136-state electronuclear fixture at
    `T=0.1 K`, `Bx=1.5 T`, and `h=0.0340929334913 meV`, its exact
    `C4(0,0,0,0)=26642.9548328` agrees with the stable large-step source-derivative plateau within
    `2.1e-5`.  For `C4(1,-1,2,-2)`, the rank-48 to rank-64 drift is `1.56e-4` relatively at
    `1.5 T`, `1.64e-4` at `4.6 T`, and `1.09e-4` at the symmetric `4.9 T` fixture.  These are
    empirical ladder errors, not certified discarded-state bounds.  The API rejects a rank cut
    through a roundoff-degenerate multiplet, records every dense/action contribution and similarity
    amplification, and leaves `functional_use_authorized=false`.

MATLAB `checkcode` returned no findings for all twenty-two isolated implementation files and the two
extended coupling APIs.

## 3. Demonstrated state selection

For the fixture

```text
Delta=1, M=1, beta=5, H=0, J0=0.7,
Jmodes=linspace(-0.12,0.12,41), n=0:500,
```

the prototype found:

| `h` | `m` | `f` | reduced curvature | classification |
|---:|---:|---:|---:|---|
| `-0.48622` | `-0.69461` | `-0.028999` | `+0.63151` | stable minimum |
| `0` | `0` | `-0.0026886` | `-0.18611` | unstable |
| `+0.48622` | `+0.69461` | `-0.028999` | `+0.63151` | stable minimum |

Both symmetry-related minima are returned with status `degenerate_minima`; no numerical seed or root
ordering chooses between them. At `H=0.02`, three roots remain but the positive-moment state has the
unique lowest functional value.

This is the capability missing from the present ordered Jensen fixed-point map: all stationary roots
are compared by one declared functional.

## 4. First production-input pilot

The read-only adapter was evaluated at `T=0.1 K`, pure transverse `B=1.5 T`, the actual legacy
`16^3`/brute-force scalar coupling multiset (`16384` modes), and its derived uniform couplings:

```text
Jmin = -0.00676310032 meV
Jmax = +0.00598513893 meV
J0   = +0.00642443566 meV
Jaa0 = +0.00351044621 meV
Delta = 0.0985077276 meV
M = 5.3330047
```

Complete root enumeration at `N_omega=128` finds two symmetry-related selected states. The positive
one is

```text
h = 0.0329775225243865 meV
m = 5.13313920359421
f = -0.0487592744634591 meV/site
reduced curvature = 0.0790074322453992 meV
minimum ring denominator = 0.932214180130684
analytic free-energy tail bound = 1.00049e-9 meV/site
```

The selected-root count remains two for `N_omega={32,64,128}`. Between `64` and `128`, the maximum
changes in `[|h|,|m|,|f|,curvature,min_den]` are approximately
`[2.7e-9,4.2e-7,6.7e-9,2.5e-7,1.5e-8]`. Thus the moderate-field state itself is easy to obtain from
the strict common functional; this supports the observation that the current Jensen failure at
`1.5 T` is not caused by absence of an ordered state.

On legacy `Nq={8,12,16}` grids at `N_omega=64`, every grid retains the same two minima. The final
`12 -> 16` changes in moment and reduced curvature are about `1e-4`. The extremal minimum
denominator changes by about `1.8e-3`; it is reported and required to stay positive, but is not
graded as a normalized BZ quadrature observable because the closest mesh point to an excluded Gamma
point changes with `Nq`. This diagnostic does not supersede the repository's existing conclusion
that the legacy endpoint-inclusive/brute-force coupling quadrature is not a rigorous production
limit.

The same strict pilot also gives a decisive limitation. On the actual `16^3` modes it has selected
ordered minima through about `3.1 T`, no admissible stationary root for `3.2--3.5 T`, and a stable
paramagnetic root from `3.6 T`; the intervening interval ends on the Gaussian ring pole
`min_q(1-J_q C_2)=0`. A production transition cannot pass through a no-state interval. Moreover,
at `1.5 T` the fixed electronic-doublet local and uniform-RPA energies are about `88.32 GHz` and
`85.05 GHz`, respectively, far above the requested `0--6 GHz` window. The hyperfine manifold is
therefore indispensable for the spectral problem.

The first full electronuclear local oracle was therefore added before any production wiring. It
deliberately uses `transverse_mf='none'`, because importing a self-consistent transverse molecular
field without also varying its conjugate moment would break the one-generator construction. At
`T=0.1 K`, `Bx=1.5 T`, and `h=0.033 meV`:

```text
max |C2 - invz_chi0z_cc| = 1.892e-13
|m + d f0/dh|            = 1.450e-10
|chi - d m/dh|           = 2.304e-10
|chi + d2 f0/dh2|        = 5.474e-6
```

The maximum source-step-halving changes in `dC2/dh` and `d2C2/dh2` were `2.48e-7` and `6.11e-3`;
the latter is about `2e-7` relative to the largest second derivative. Positive/negative source
symmetry closes to `1.4e-12` for `C2` and below `8e-8` relative for its second source derivative.

On the actual `16^3` coupling multiset, the positive electronuclear minimum at `1.5 T` is

```text
h = 0.0340929334913 meV
m = 5.30675931175
minimum ring denominator = 0.966617867
```

Narrow-bracket complete re-solves at `N_omega={128,256,512}` retain one positive minimum (the
negative partner follows from the verified `Z2` symmetry). The `256 -> 512` changes in `h`, `m`,
`f`, curvature, and minimum denominator are respectively about
`[1.5e-9,2.34e-7,1.67e-7,1.97e-7,2.8e-9]`; the analytic free-energy tail bound at `512` is
`5.39e-7 meV/site`.

Hyperfine structure narrows and shifts the electronic no-state interval: stable ordered roots were
found at `3.6,4.0,4.5,4.6 T`, with `h` decreasing continuously from `0.0233` to
`0.00912 meV`; the `4.6 T` moment is `1.4201`. No admissible root was found at
`4.7,4.8,4.85 T`, and the symmetric root is stable at `4.9 T`. Thus the full local spectrum puts the
transition on the expected field scale but does not remove the strict Gaussian pole/no-state gap.

## 5. Limits and next gates

This prototype does not yet justify a LiHoF4 production calculation. The next contained steps are:

1. do not revive the rejected Gaussian-local-trace plus 1PI-vertex skeleton;
2. freeze a multilevel local-rank/cutoff error budget, then derive the nonlinear fixed-source local
   partial-Legendre functional before attempting another stationary lattice functional; or activate
   the separately documented biased smooth-`r(h)` backup as an explicit experimental prescription;
3. do not wire a spectral backend until whichever route is chosen passes its thermodynamic, domain,
   discretization, and branch-identity gates.

The deferred `O(1/z^2)` non-Gaussian vacuum and the production ordered `tanh/xi` machinery remain out
of scope.
