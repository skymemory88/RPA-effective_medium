# Biased smooth-`r` diagnostic oracles

**Recorded:** 2026-07-28

**Status:** isolated diagnostic; not on the production MATLAB path and not wired into spectra

This directory freezes the decision-only part of the backup prescription in
`../../../biased_convergence_solution.md`. The function

```matlab
result = invzp_select_smooth_r_oracle(paths,spec)
```

accepts only already-enumerated candidate paths. It validates their certificates, computes the two
registered dimensionless smoothness metrics, and applies the registered lexicographic selector. It
does **not** find roots, connect branches, choose endpoints, differentiate data, interpolate,
resample, smooth, or repair missing nodes.

`invzp_trace_pseudo_arclength` is the separate, generic branch-tracing primitive. It follows one
regular scaled residual curve with a bordered pseudo-arclength corrector and retains every full
coordinate vector. It has no knowledge of the invZ equations and deliberately does not find initial
roots, join disconnected components, convert a folded curve into a Jensen path, or invoke the
selector. Its coordinate contract is `y=[scaled_state; scaled_parameter]`, with the continuation
parameter last. Every rejected adaptive attempt is retained. Persistent hard-event failures terminate
with `status='event'`; pointwise event callbacks alone cannot prove that a narrow boundary was not
crossed and re-entered between evaluations, so an invZ adapter must additionally report continuous
signed margins and bracket their changes.

`invzp_ordered_arclength_problem` is the first invZ-specific adapter. It uses the audited resummed
node equations and a Richardson-refined centred derivative in the fixed longitudinal field, with
explicit scales

```matlab
y = [Sigma(:)./sigma_scale; K0/K0_scale; h/h_scale].
```

Its accepted payload retains the full state, A--D audit, reciprocal-chart diagnostics, signed event
margins, and the field-Jacobian refinement drift normalized by
`max(1,norm(dR/dh,Inf))`. It still adapts one supplied initial root; it is not a root enumerator or
Jensen-path constructor.

`invz_ordered_node_jacobian_factors` and `invzp_solve_bordered_factors` expose and solve the exact
frequency-block form

```text
J_sigma,sigma = diag(d) + U*V',  rank(U*V') <= 2,
```

with `K0` and `h` left as a two-variable border. At the 742-variable 1.5 T fold, the factor solve
agrees with the dense solve to `1.66e-15` in the solution norm and `1.60e-14` backward residual.
It is an oracle, not currently dispatched by the tracer.

## Input boundary

Every path must contain:

| field | contract |
|---|---|
| `id` | unique nonempty character row or string scalar |
| `hstar` | finite positive endpoint found by that path's endpoint solver |
| `x` | the exact common mesh in `spec.x` |
| `r`, `drdx`, `d2rdx2` | finite real vectors on `x`; derivatives come from the tracked representation |
| `r_jump_free` | enumerator certificate that the path contains no finite jump in `r` |
| `complete` | exactly one accepted state exists throughout `[0,hstar]` |
| `admissible` | all residual, A--D, domain, event, and provenance gates passed |
| `endpoint_ok` | the unique increasing `F_path=0` and `crit>0` gates passed |
| `qcp_ok` | the external preregistered QCP-continuity gate passed |
| `max_state_curvature` | finite nonnegative tie metric, or `NaN` if unavailable |
| `qcp_distance` | finite nonnegative tie metric, or `NaN` if unavailable |
| `trace_agreement_available` | the forward/reverse/cold comparison was actually completed |
| `trace_agreement_ok` | completed forward/reverse/cold traces agree; must be false when unavailable |

The scalar `spec.version`, exact dimensionless `spec.x`, and absolute/relative tolerances
`spec.tie.{slope,curvature,state,qcp}_{abs,rel}` are mandatory. A mismatched mesh, duplicate ID,
nonfinite primary datum, nonlogical certificate, or invalid endpoint raises
`invzp:SmoothSelector:InvalidInput`; it is not converted into a physical branch status.

The primary metrics are exactly

```matlab
max_slope_change = max(abs(diff(path.drdx)));
integrated_r_curvature = trapz(spec.x,abs(path.d2rdx2).^2);
```

and lower values are retained within

```matlab
abs(a-b) <= abs_tol + rel_tol*max(abs(a),abs(b)).
```

The order is binding: admissibility and jump rejection, maximum resolved slope change, integrated
curvature, maximum state curvature, QCP distance, then trace agreement. The function returns
`selected`, `branch_ambiguous`, or `branch_unresolved`, together with all metrics, rejected reasons,
and the survivors after every stage. Missing tie information, including an unavailable trace audit,
causes ambiguity; it is never encoded as a failed comparison or assigned an artificial bad score.
Input order and path IDs never break a tie.

## Deliberate omissions

The generic pseudo-arclength primitive now exists, but the invZ-specific root enumerator,
bidirectional/cold-seed reconstruction, signed event-crossing audit, single-valued Jensen-section
builder, endpoint reconstruction, QCP field sequence, physical-`h` diagnostic metrics,
mesh-refinement drift, and cross-normalization winner audit still have to be built and frozen before
the global 1.5 T path-selection pilot. In particular, a unique result from this function means only
unique under the supplied dimensionless spec. It does not authorize production use unless the
omitted physical-`h`, normalization, and refinement audits agree under the registered grid and
Matsubara refinements.

Fresh in-memory MATLAB fixtures cover a known linear winner, lexicographic dominance, curvature and
each tie-break stage, finite-jump rejection, no accepted candidate, an exact unresolved tie, missing
tie data, and invalid-input cases. No persistent test script is kept in the worktree.

## Verification record

The generic tracer crossed both folds of

```matlab
R(u,h) = u^3-u+h
```

and recovered the exact fold fields `h=+/-2/(3*sqrt(3))` within the declared arc discretization.
Fresh guards also verified reverse orientation, retry retention, invalid-tangent rejection, and
distinct terminal statuses for callback and nonfinite-residual events.

On the frozen legacy 16³ LiHoF4 coupling at `T=0.10 K`, `Bx=1.5 T`, the invZ adapter reproduced the
previously observed high-branch fold with a retained full-state trace:

```text
initial h                         0.042371329941803114 meV
accepted unique roots            149
adaptive attempts                156
maximum residual_inf             4.535304520825179e-09
minimum bordered rcond           5.073311652228108e-06
minimum oriented tangent overlap 0.9501876552804539
minimum pole margin              0.006213739195997094
minimum mean margin              0.03106320582704512
maximum h-Jacobian scaled drift  7.825112291199848e-05
Hermite fold field               0.0052435482911986821 meV
```

Every corrected root passed the independent A--D audit. The fold value differs by
`1.44e-10 meV` from the earlier one-off quarter-step result
`0.005243548147122800 meV`. The dense/Richardson oracle took 1292.9 s over six bounded segments.
Subsequent profiling corrected the initial bottleneck attribution: at the fold, dense solve and
`rcond` take only `0.00745 s` and `0.00670 s`, while one Richardson equation/Jacobian evaluation takes
`1.771 s`. The exact factor solve takes `0.00290 s`, too small a saving to justify replacing the
dense oracle now. Exact last-point caching instead reduces a repeated evaluation from `1.748 s` to
`6.2e-5 s` and shortened a fresh three-step fold trace to `5.67 s/step` without changing arithmetic.
The next optimization target is redundant local node/field-derivative evaluation, not linear
algebra. Temporary `.mat` traces stayed under `/tmp` and are not retained.
