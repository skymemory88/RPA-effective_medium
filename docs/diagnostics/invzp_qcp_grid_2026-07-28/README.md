# QCP finite-grid and solver-edge diagnostics — 2026-07-28

## Question and frozen configuration

This packet tests the narrow hypothesis that the apparently intermittent ordered
Jensen coverage near the quantum critical point is controlled by the finite
Brillouin-zone coupling edge, rather than by physical critical slowing down.
It changes no production solver.

The frozen state gate uses:

- `T = 0.10 K`;
- legacy production grid construction (grid size present, all
  `gridConvention`/`gridOffset`/`gammaPolicy` fields absent);
- brute-force dipoles, `dpRng = 30`, cached couplings;
- grids `N = 12,16,20,24`;
- common ordered field mesh `4.000:0.025:4.700 T`;
- default resummed static medium, `legacy_x`, and unchanged solver tolerances;
- a converged-PM mass bisection on `[4.55,4.85] T`, stopped at a field
  bracket narrower than `2e-5 T`.

`W` below is explicitly `Jmax-Jmin`.  `2/W` is recorded as a moment-expansion
scale, not promoted to an acceptance gate.

## Coupling-only prediction

| N | Jmax (meV) | J0−Jmax (meV) | next-level gap (meV) | W (meV) | 2/W (meV⁻¹) | S_N(J0) (meV⁻¹) |
|---:|---:|---:|---:|---:|---:|---:|
| 12 | 0.00582319148013 | 0.000601244175573 | 3.90192520763e-5 | 0.0121180635243 | 165.042871412 | −196.730243280 |
| 16 | 0.00598513892920 | 0.000439296726499 | 2.27524792339e-5 | 0.0127482392457 | 156.884410580 | −198.006079765 |
| 20 | 0.00608577926539 | 0.000338656390309 | 2.82177980712e-5 | 0.0131568074433 | 152.012561453 | −198.572504603 |
| 24 | 0.00614300357268 | 0.000281432083025 | 1.11639199408e-5 | 0.0134043034205 | 149.205813779 | −199.384362995 |

The excluded-Γ gap fits approximately `N^-1.103`.  At fixed `delta=0`,
linear and quadratic fits in `1/N` give `S_inf(J0)=-201.8063` and
`-202.8384 meV^-1`; their `1.03 meV^-1` spread is retained as model
sensitivity, not called a continuum error bar.  Fixed positive deltas are in
`coupling_scan.mat`.

At the `16^3` PM mass root, `G(0)=-198.0070 meV^-1`, equal within the
mass-root bracket to `S_16(J0)=-198.0061 meV^-1`, while
`|G|/(2/W)=1.262`.  This is a useful identity/control diagnostic.  It is not
by itself a proof that the resummed theory must be rejected.

## State-only grid gate

| N | solver-grade Bc¹ᐟᶻ (T) | accepted total | stable total | QCP-contiguous count | QCP-contiguous width (T) | lowest field in suffix (T) |
|---:|---:|---:|---:|---:|---:|---:|
| 12 | 4.682284546 | 23 | 23 | 16 | 0.382284546 | 4.300 |
| 16 | 4.692758179 | 14 | 14 | 11 | 0.267758179 | 4.425 |
| 20 | 4.699093628 | 11 | 11 | 9 | 0.224093628 | 4.475 |
| 24 | 4.702957153 | 11 | 11 | 8 | 0.177957153 | 4.525 |

The QCP-contiguous suffix means the highest uninterrupted run of accepted
fields on the common mesh strictly below the independently measured PM mass
root.  Its widths times `N` are
`[4.5874,4.2841,4.4819,4.2710] T`, and a log fit gives
`width ~ N^-1.076`.  This strongly supports a finite-grid computability
sliver.

The suffix's lower edge is quantized by the declared `0.025 T` field mesh,
so the fitted exponent is descriptive rather than a precision estimate.
The `Bc` values themselves use the much narrower mass-root brackets stored
in `state_grid_gate.mat`.

The preregistered prediction that `Bc` would move by at most `0.01 T` from
`12^3` to `24^3` is falsified: the measured shift is `0.02067 T`.
Linear and quadratic `1/N` fits give `4.72386` and `4.72224 T`, respectively.
Those two close fits are encouraging but four grids do not supply a rigorous
continuum extrapolation.

Accepted islands below the contiguous suffix remain.  On the common mesh,
`A` means accepted and `x` means `node_failed`:

```text
N=12  AAAxxxAAxAAxAAAAAAAAAAAAAAAA.
N=16  xxAAAxxxxxxxxxxxxAAAAAAAAAAA.
N=20  xxxxxxxxAxxAxxxxxxxAAAAAAAAA.
N=24  xxxxxxxxxAAxxxAxxxxxxAAAAAAAA
```

The final dot is the field above the grid-specific `Bc` where no ordered
root is expected.  The islands demonstrate that the failure is not a
uniform iteration threshold.

## Phase-aligned susceptibility/peak grid gate

The user-priority observable was then evaluated at matched physical offsets
from each grid's solver-grade mass root:

```text
B-Bc = [-0.080,-0.040,-0.020,-0.010,-0.005,
         0.005, 0.010, 0.020, 0.040, 0.080] T.
```

The real-axis grid was `0:0.002:2.0 GHz` with the production
`eta=5e-5 meV`. All 40 susceptibility columns were finite. Every negative
offset was Jensen ordered and every positive offset was paramagnetic.

![Absolute-field shift and phase-aligned peak collapse](peak_grid_gate.png)

| N | −0.080 T | −0.020 T | −0.005 T | +0.005 T | +0.020 T | +0.080 T |
|---:|---:|---:|---:|---:|---:|---:|
| 12 | 0.938734 | 0.470125 | 0.235203 | 0.165653 | 0.330028 | 0.650630 |
| 16 | 0.936204 | 0.469271 | 0.234970 | 0.165070 | 0.329221 | 0.649491 |
| 20 | 0.935292 | 0.468637 | 0.234286 | 0.165112 | 0.328962 | 0.649001 |
| 24 | 0.934363 | 0.467917 | 0.234053 | 0.164782 | 0.328372 | 0.648158 |

Entries are the parabolically refined peak of the retained
`chi''_cc(omega)` column in GHz. Across all ten offsets the four-grid
relative spread is `0.38--0.53%`; the largest absolute spread is
`0.00437 GHz` at `-0.080 T`.

A `0.001 GHz` refinement at `N=12,24` and offsets
`[-0.020,-0.005,+0.005,+0.020] T` reproduces the coarse-grid extreme-grid
spreads as `[0.002212,0.001140,0.000855,0.001650] GHz`; individual peak
locations change by at most about `1e-5 GHz`. The peak-position comparison
is therefore not limited by the original `0.002 GHz` sampling.

This gives a useful separation:

- the phase-aligned soft-mode shape is already stable at the sub-percent
  level over `12^3--24^3`;
- the dominant grid effect is the horizontal `Bc(N)` shift, which is not
  yet converged.

`peak_grid_gate.mat` retains the full susceptibility columns for inspection.
Only peak position and phase provenance were graded here; spectral weight,
linewidth, and an independently solved analytic pole were not.

## Edge-pair trace and the monotone-inner-solve proposal

The representative `16^3` pair is 4.400 T (rejected) and 4.425 T
(accepted).  Every traced node and outer iterate at both fields has
`y_rank=16384` and `y>Jmax`: neither run changes discrete pole interval.
Thus a pole-index switch does not explain this edge pair.

At 4.400 T only the `h=0` predictor fails.  Its inner static closure remains
closed (`resid_static` about `5e-11`), but the outer Σ-map oscillates and
ends at `8.41e-6` after 200 iterations, with an asymptotic residual ratio
about `0.916`.  Its edge distance is only
`y-Jmax=1.4422e-5 meV` (`min|Dq|=0.00433`).  At 4.425 T the predictor has
`y-Jmax=3.6651e-5 meV` and reaches `6.37e-9` in 13 iterations.

A control with `max_outer=1000` accepts 4.400 T, but 4.300 T still returns
`node_failed` after 1000 iterations with Block-A residual `0.00833`.
Increasing the cap therefore moves the mask edge; it does not remove it.

Source inspection also limits the proposed P2 shortcut:
`invz_emt_static_ordered` recomputes `Gstat(K0)` on every inner K0 iterate,
and the failed edge already closes that inner residual.  Solving
`S_N(y)=constant` on the rightmost interval would instead freeze a quantity
that the implemented equation varies and would target the already-converged
block.  A branch-indexed scalar solver remains a possible diagnostic only
after deriving the coupled scalar equation and proving its monotonicity; it
is not a justified drop-in production fix.

The previously reported two-root state at
`h=0.01171732294388717 meV` cannot be replayed deterministically because its
full state vectors/seeds were retained only in temporary files.  This packet
does not infer their interval indices from rounded `r` values.

## Area-rule oracle

For a converged continuous branch,

```text
F_direct(h) = h0(h) - J0*m(h)
F_area(h)   = integral_0^h crit(h') dh'
```

are identical because `dm/dh=-G0bare` and
`crit=r+J0*G0bare`.  `invzp_area_rule_oracle.m` compares them on the same
stored profile; it does not fill a missing node.

At 4.05 T:

| h-grid | profile verdict | max \|F_direct−F_area\| (meV) | linear direct zero (meV) | linear area zero (meV) |
|---:|---|---:|---:|---:|
| 33 | ok | 2.43958e-5 | 0.0146834072 | 0.0146920938 |
| 65 | node_failed (64/65) | 6.13689e-6 | 0.0147132258 | 0.0147152452 |

The approximately fourfold route-difference reduction is consistent with
trapezoidal second-order convergence.  Near the 33-node crossing the direct
subtraction condition number reaches about 390, so the area form is a useful
cancellation meter.  The 65-node run nevertheless hits one newly sampled
failed state at `h=0.0033154657 meV`: Block A/C are about `0.0156/0.0150`
while the inner static residual is `4.24e-11`.  Refining the quadrature can
therefore expose, not solve, the same state-availability problem.
The two interpolated zeros on the failed 65-node profile are diagnostic
only and are not accepted Jensen endpoints.

## Decision

1. Keep the current accepted QCP spectrum as a finite-`16^3` preview. Its
   phase-aligned peak shape is sub-percent stable across the tested grids,
   but its absolute field alignment, ordered width, and `Bc` are not
   grid-converged.
2. Do not implement fixed-G rightmost-root bisection as P2: it does not match
   the current coupled static equation and the measured failing block is the
   outer Σ closure.
3. Do not use a larger iteration cap as the fix.
4. The next production-changing candidate must address the pole-conditioned
   outer map while preserving branch identity and must pass the same grid
   ladder.  A continuum-density treatment remains a theory change and needs
   a separately specified quadrature/edge prescription.

## Reproduction and retained data

- `invzp_qcp_coupling_scan.m` regenerates `coupling_scan.mat`.
- `invzp_qcp_state_grid_gate.m` regenerates `state_grid_gate.mat`; it performs
  no real-axis response calculation.
- `invzp_qcp_peak_grid_gate.m` regenerates `peak_grid_gate.mat` from the
  retained solver-grade `Bc(N)` values and stores the full response columns.
- `edge_pair_trace.mat` retains compact nodes/iterations with large lattice
  arrays removed.
- `area_rule_grid.mat` retains the two profiles and scalar `J0eff`;
  `invzp_area_rule_oracle.m` recomputes the comparison.

All MAT files use full-precision doubles.  Console values in this document
are rounded for readability.
