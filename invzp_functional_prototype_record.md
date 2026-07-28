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
  zero-diagonal coupling matrices, thermodynamics, and stable KMS/Hermite Lehmann response;
- `invzf_two_site_exact`: the pair-specialized wrapper used by the Jensen dynamic gate.

Every function is disconnected from `invz_projected/`, `invz_common/`, and `invz_tensor/`.

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

MATLAB `checkcode` returned no findings for all seven implementation files.

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

## 4. Limits and next gates

This prototype does not yet justify a LiHoF4 production calculation. The next contained steps are:

1. add mode-grid comparison/error plumbing for lattice inputs;
2. add more coupling topologies before generalizing the exact local source derivatives beyond the
   scalar two-level model;
3. only afterward investigate a stationary EMT/2PI extension.

The deferred `O(1/z^2)` non-Gaussian vacuum and the production ordered `tanh/xi` machinery remain out
of scope.
