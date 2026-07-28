# Smooth-`r` selector implementation oracle

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

The path enumerator, pseudo-arclength continuation, endpoint reconstruction, QCP field sequence,
physical-`h` diagnostic metrics, mesh-refinement drift, and cross-normalization winner audit still
have to be built and frozen before the 1.5 T pilot. In particular, a unique result from this function
means only unique under the supplied dimensionless spec. It does not authorize production use unless
the omitted physical-`h`, normalization, and refinement audits agree under the registered grid and
Matsubara refinements.

Fresh in-memory MATLAB fixtures cover a known linear winner, lexicographic dominance, curvature and
each tie-break stage, finite-jump rejection, no accepted candidate, an exact unresolved tie, missing
tie data, and invalid-input cases. No persistent test script is kept in the worktree.
