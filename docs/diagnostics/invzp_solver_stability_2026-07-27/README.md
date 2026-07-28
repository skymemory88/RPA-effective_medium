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
  selected Jacobian `rcond > 1e-12`; raw Jacobian by default, or the explicitly enabled
  row-equilibrated Jacobian; A–D audit remains binding
- `row_equilibrate=false` and `static_polish=false` by default. When explicitly enabled,
  equilibration affects only the linear solve/rank estimate. The physical-`K0` polish permits one
  capped proposal per eligible Newton iterate and returns immediately on acceptance; raw residual
  and A–D acceptance never change.

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
| 0.10 | 1.5 | retained 25-seed `h=0` root census | 16 A--D-accepted attempts cluster into seven distinct roots in the complete independent coordinates `[Sigma;K0]` (`K` and `lambda` are derived); nine attempts fail. Representative `r` ranges from `0.60295` to `1.20283`, so the earlier three-seed agreement sampled one basin rather than proving uniqueness |
| 0.10 | 1.5 | staged reconstruction of the low fold and returning leg | 21 fixed-`h` nodes plus 20 local pseudo-arclength steps recover `h_fold=1.1303917626651595e-5`; the returning leg exposes raw-closure amplification by `G*Gbar`, while the defactored equation remains accurate to about `2e-11`. The trace reaches `h=6.890625e-9` but does not establish a connection to any `h=0` root |
| 0.10 | 1.5 | all-root simple-first upward screen | the other zero-field roots fail fixed-`h` continuation nonuniformly rather than at a common physical event. A retained `8e-9` arclength tolerance was too close to the independent `1e-8` A/C Sigma gates; using `1e-9` in the defactored diagnostic coordinate removes the observed root-4/root-5 audit veto without changing a tolerance |
| 0.10 | 1.5 | root-4 defactored local trace | 97 accepted steps cross a fold between `5.5677706301911396e-6` and `5.568143040810256e-6 meV`; `min rcond(border)=4.35e-10`, with no audit or domain-event failure before the trace stops at a static-component line-search floor |
| 0.10 | 1.5 | root-6 strict trace and refined zero-field return | a fold occurs at `9.541141092330816e-6 meV`; the returning leg crosses `h=0` without audit failure. Refinement brackets zero within `6.76e-10 meV`, and the corrected endpoint matches census root 6 at frozen cluster distance `1.42e-14` after a `7.43e-12` predictor correction |
| 0.10 | 1.5 | root-6 opposite-leg fixed-`h` check | all 33 decreasing fixed-`h` nodes from the common handoff to zero pass A--D without an event failure; maximum frozen predictor distance is `0.0247702`, and the endpoint matches census root 6 at cluster distance `7.92e-14` |
| 0.10 | 1.5 | default-off physical-`K0` last-bit polish | at the positive-q `q=0.7625` stall, all Sigma equations pass and the static residual is `1.33548e-8`; one 656-physical-ULP proposal reduces the complete raw residual to `1.13e-11` and passes A–D |
| 0.10 | 1.5 | row-equilibrated, derivative-free final fixed-`h` descent | the static Jacobian row is `3.53e6` times larger than every Sigma row: raw/equilibrated `rcond` at the handoff are `8.27e-13`/`6.10e-9`. The retained fixed-`h` helper accepts 52 scheduled points down to `h=3.4852605192481022e-10`; a refined schedule reaches `3.4035747258282251e-10`, then fails line search at raw residual `3.83e-6`. A direct `h=0` solve from the last accepted state fails at `5.077e-2`; that state is at least `0.312` from every root in the seven-root census |
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

The subsequent 25-seed census corrects one limited inference in that account: agreement among the
three original cold seeds was only basin-local. Seven accepted `h=0` roots are now resolved under the
frozen clustering tolerances. The short staged trace still verifies the low regular fold, but its
returning leg appears to become a sharp, nearly rank-deficient boundary layer near `h=0`; no accepted
continuation from that leg to any enumerated endpoint has yet been demonstrated.

A later scaling audit shows that the apparent raw rank loss in the final boundary layer is mostly
row imbalance, not algebraic rank loss. The static row reaches millions of times the scale of the
self-energy rows, while row equilibration leaves `rcond≈6.1e-9`. The new options in
`invz_ordered_node_newton` are therefore default-off and deliberately narrow:
`row_equilibrate=true` divides each residual row and right-hand side by the same Jacobian infinity
norm only for the linear solve/rank test; `static_polish=true` proposes one rounded physical-`K0`
Newton update at each eligible Newton iterate only after all Sigma residuals pass. A failed proposal
does not prevent an ordinary full Newton step, while a successful proposal is fully gated and
returned immediately. Every accepted state still uses the raw residual, the existing domain events,
and the independent A–D audit.

Those mechanisms extend the positive branch by roughly a factor of five in `h`, but do not connect
it to zero. Refinement stops reproducibly near `3.4e-10 meV`, a direct zero solve is rejected, and
the last positive state is not within the frozen census tolerance of any zero-field root. This
strengthens the endpoint non-connection result; it does not authorize extrapolation across the
unresolved interval.

That unresolved boundary layer belongs to the previously staged root-7 component; it is not a
universal zero-field obstruction. A later strict trace of a different component (root 6) crossed a
regular fold, returned through `h=0`, and was locally refined to the root-6 census medoid within the
frozen same-root tolerance. A 33-node fixed-`h` check of the opposite root-6 leg also terminates at
the root-6 medoid with every A--D gate accepted. Root 4 supplies another distinct regular fold. This
is consistent with the user's visual observation that ordered-state nonconvergence is nonuniform
and numerical in part, while also
showing why a single fallback solve cannot define the Jensen path: several accepted branches coexist.
The simple numerical corrections are therefore retained only as branch diagnostics until a complete
single-valued candidate and its selection gates exist.

One-off drivers and `.mat` files stayed under `/tmp`. The audited numerical kernel now consists of
`invz_ordered_node_context.m`, `invz_ordered_make_node.m`,
`invz_ordered_node_equations.m`, `invz_ordered_node_jacobian_factors.m`, and
`invz_ordered_node_newton.m`. The extraction leaves the fixed-node equations unchanged while exposing
the full residual/Jacobian, its exact diagonal-plus-rank-two frequency block, and exact HMF node
construction to the separately isolated continuation oracle. The retained cold-seed enumerator now
records every attempt and clusters accepted `[Sigma;K0]` coordinates without order-dependent tie
breaking. It still needs a complete seed-coverage argument, forward/reverse branch tracking, cold
checkpoints, grid/cutoff checks, event-crossing control, and Jensen-integral error control before
promotion.
