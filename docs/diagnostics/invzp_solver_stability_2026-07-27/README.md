# Ordered-node solver-stability pilot — 2026-07-27

This directory retains the reusable kernel from a solver-only diagnostic prompted by a first-hand
inspection of `invz_run_spectra`: convergence was non-uniform in the ordered field range, with usable
columns near the apparent QCP while moderate-field columns such as 1.5 T were masked.

The helper is deliberately **not on the production MATLAB path**. It solves the unchanged resummed
fixed-`h` equations with a defactored reciprocal-chart residual and analytic dense Jacobian, then
requires the independent production A–D audit. It does not choose a thermodynamic branch and must not
be inserted into Jensen quadrature without branch-tracked continuation.

## Frozen inputs

- MATLAB R2025a
- legacy 16³ production grid, brute-force dipoles, `dpRng = 30`
- exact coupling digest:
  `ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17`
- `static_medium = 'resummed'`, `Ecut = 40`
- Picard: `mix_outer = 0.7`, `max_outer = 200`, `tol_outer = 1e-8`
- diagnostic Newton events: `pole_margin > 1e-10`, `mean_margin > 1e-10`,
  Jacobian `rcond > 1e-12`; A–D audit remains binding

## Results

| T (K) | B (T) | experiment | result |
|---:|---:|---|---|
| 1.00 | 2.9 | positive-control reciprocal residual | residual `1.888049e-09`; analytic-vs-centred-FD Jacobian relative error `1.882e-11`; Newton-step relative error `5.632e-10` |
| 0.10 | 3.6 | production Picard | `node_failed`, 30/33 profile nodes accepted |
| 0.10 | 3.6 | reverse Picard traversal | the same nodes 22, 24, 27 fail; warm `h=0` succeeds, so traversal alone is not the cause |
| 0.10 | 3.6 | Picard plus residual corrector at only those nodes | 34/34 including `h=0`; failed nodes corrected in 4, 4, and 8 Newton iterations |
| 0.10 | 3.6 | temporary end-to-end experimental wiring | `status='ok'`, `hmf=0.01949623263696515` meV, 33/33 profile nodes, stable endpoint, 34.8 s |
| 0.10 | 1.5 | production Picard | `node_failed`, 11/33 profile nodes |
| 0.10 | 1.5 | same hybrid pilot | `h=0` and two additional nodes corrected, but only 13/33 profile nodes; the path remains incomplete |
| 0.31 | 3.6 | adaptive natural-`h` reciprocal continuation from a clean high-`h` node | reaches `h=0` in 12 accepted nodes; minimum pole margin `1.4827e-3` |
| 0.31 | 1.0 | same continuation | step collapse near `h=0.00793448728`; 258 accepted nodes, condition number rising to `9.84e7`, no pole/mean event |

The temporary HMF wiring used to obtain the end-to-end 3.6 T result was removed after review: a
fixed-`h` residual root is not automatically the same continuous branch required by the Jensen
integral. The result is evidence that Picard non-contraction is a numerical defect at that field, not
yet authorization to publish a corrected spectrum. The 1.5 T gap and the 0.31 K/1 T fold show why a
branch-identity check is still necessary.

One-off drivers and `.mat` files stayed under `/tmp` and were not retained. The surviving
`invz_ordered_node_newton.m` is the audited numerical kernel; a future retained continuation driver
must add forward/reverse branch tracking, cold checkpoints, grid/cutoff checks, and Jensen-integral
error control before promotion.
