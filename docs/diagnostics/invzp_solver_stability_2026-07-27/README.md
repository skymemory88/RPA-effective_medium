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
| 0.10 | 1.5 | first adaptive natural-`h`/secant experiment from the clean high-`h` branch | 300-node budget reached at `h=0.003414300661` without a pole/mean event; maximum residual `5.35e-12`. A later branch-controlled trace shows this experiment did not preserve branch identity through the fold |
| 0.10 | 1.5 | scaled tangent predictor with fixed-`h` correctors | stops at `h=0.005243548157786061` after 34 accepted nodes: `rcond(J)=1.027e-11` and scaled `abs(deta/ds)=1.771e-6`, with pole margin `6.640e-3`, mean margin `1.735`, and every A–D audit accepted |
| 0.10 | 1.5 | fold rank audit | `sigma_min(J)=9.672e-6` while `sigma_min([J R_eta])=0.4581` (`cond(J)=1.336e8`): the fixed-`h` Jacobian loses rank while the augmented tangent remains regular |
| 0.10 | 1.5 | scaled pseudo-arclength across that event | crosses in 24 accepted steps; tangent `deta/ds` changes from `-1.771e-6` to `+8.358e-4`, `min rcond(augmented)=1.953e-5`, maximum branch-correction ratio `2.663e-4`, and all A–D audits pass. A repeat at one-quarter the initial arc step brackets `deta/ds=0`, resolves `min(h)=0.005243548147122800`, and keeps successive oriented-tangent overlap above `0.9999999999` |
| 0.10 | 1.5 | retained full-state pseudo-arclength oracle | 149 accepted roots from the clean high-`h` endpoint through the fold; every A--D audit passes, `max residual=4.5353e-9`, `min bordered rcond=5.0733e-6`, and local Hermite interpolation gives `h_fold=0.0052435482911986821`, within `1.44e-10 meV` of the earlier one-off quarter-step result |
| 0.10 | 1.5 | independent low-`h` branch, three cold deterministic seeds | all three seeds converge to the same `h=0` root within scaled-state diameter `4.278e-12`; upward tangent continuation stops after 22 accepted nodes at `h=1.130400753628218e-5`, with `rcond(J)=2.041e-11` and healthy event margins |
| 0.10 | 1.5 | low-branch rank audit | `sigma_min(J)=1.558e-5` while `sigma_min([J R_eta])=0.5629` (`cond(J)=1.222e8`): this is a second regular fold, about 464 times below the high-branch fold |
| 0.10 | 1.5 | pseudo-arclength across the low-branch event | crosses in 11 accepted steps; `deta/ds` changes from `+2.677e-7` to `-5.243e-7`, directly resolving `max(h)=1.130402502867847e-5`; minimum oriented-tangent overlap is `0.9999999997`, `min rcond(augmented)=1.091e-6`, and all A–D audits pass |
| 0.31 | 3.6 | adaptive natural-`h` reciprocal continuation from a clean high-`h` node | reaches `h=0` in 12 accepted nodes; minimum pole margin `1.4827e-3` |
| 0.31 | 1.0 | same continuation | step collapse near `h=0.00793448728`; 258 accepted nodes, condition number rising to `9.84e7`, no pole/mean event |

The temporary HMF wiring used to obtain the end-to-end 3.6 T result was removed after review: a
fixed-`h` residual root is not automatically the same continuous branch required by the Jensen
integral. The result is evidence that Picard non-contraction is a numerical defect at that field, not
yet authorization to publish a corrected spectrum. At 1.5 T, the branch-controlled traces found two
regular folds rather than a pole or mean-cancellation event. The clean high-`h` component has a
minimum at `h=0.0052435481471` and turns toward increasing `h`; the independently cold-seeded low-`h`
component has a maximum at `h=1.1304025029e-5` and turns toward decreasing `h`. They therefore do not
overlap in the traced neighborhoods, and neither supplies a direct single-valued fixed-`h` path over
the Jensen interval. The earlier unconstrained secant experiment's roots between these folds are
necessarily on at least one other branch segment and cannot be counted as continuation of the clean
high-`h` branch. Whether still more folds join these segments globally is immaterial to the immediate
decision: the local folds already give multiple real roots at the same `h`. This is the
multiple-branch outcome of the Phase-1 fork: an explicit physical or preregistered numerical branch
prescription is required before production integration. The separate 0.31 K/1 T step collapse still
shows why branch/rank diagnostics are necessary.

One-off drivers and `.mat` files stayed under `/tmp`. The audited numerical kernel now consists of
`invz_ordered_node_context.m`, `invz_ordered_make_node.m`,
`invz_ordered_node_equations.m`, and `invz_ordered_node_newton.m`. The extraction leaves the
fixed-node equations unchanged while exposing the full residual/Jacobian and exact HMF node
construction to the separately isolated continuation oracle. A future retained enumerator must still
add forward/reverse branch tracking, cold checkpoints, grid/cutoff checks, event-crossing control,
and Jensen-integral error control before promotion.
