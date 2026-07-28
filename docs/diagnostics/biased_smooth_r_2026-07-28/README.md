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
signed margins and bracket their changes. Each accepted-step search is also bounded deterministically
by `max_step_attempts` (default 8, matching the largest accepted retry count in the retained low
trace), and each corrector attempt by
`max_evaluations_per_attempt` (default 16). Exhausting the former returns
`status='budget_exhausted'` with the last numerical reason; it never promotes the last iterate. The
accepted corrector's already-gated `R`, Jacobian, and diagnostic record are carried directly into
the audit and tangent update, so there is no unmetered or unchecked post-audit equation call. A
line-search trial that already satisfies the residual, constraint, and event gates is promoted
directly with the same virtual iteration count used by the previous top-of-loop confirmation. Failed
attempt records retain their last residual, its MATLAB one-based component index and signed value,
constraint, bordered `rcond`, step norm, damping factor, and optional polish provenance. For the
ordered adapter, indices `1:nw` are the self-energy equations and `nw+1` is the defactored static
closure. These fields are diagnostics, not acceptance criteria. If the final permitted Newton update
itself lands on an accepted root, it is retained rather than discarded for lack of a redundant
confirmation iteration; the analytic
`R(u,h)=u+h^2` one-update oracle freezes that boundary semantics.

For the ordered adapter, the corrector tolerance is an accuracy control and must sit safely inside
the independent A--D gates. A retained `8e-9` setting was too close to the `1e-8` A/C Sigma gates:
predictors could satisfy the vector test before Newton corrected them enough for the independent
audit. The subsequent regular-coordinate pilot uses `1e-9`, while the stricter raw-coordinate
cross-check uses `1e-10`. These are tighter solve requirements, not changed A--D tolerances.

`invzp_ordered_arclength_problem` is the first invZ-specific adapter. It uses the audited resummed
node equations and a Richardson-refined centred derivative in the fixed longitudinal field, with
explicit scales

```matlab
y = [Sigma(:)./sigma_scale; K0/K0_scale; h/h_scale].
```

Its accepted payload retains the full state, A--D audit, reciprocal-chart diagnostics, signed event
margins, and the field-Jacobian refinement drift normalized by
`max(1,norm(dR/dh,Inf))`. It still adapts one supplied initial root; it is not a root enumerator or
Jensen-path constructor. The adapter now also exports the exact state gradient of
`w=z-K0`, a Richardson field derivative, and its refinement drift. These are used only by the
fixed-`w` conditioning probe.

`invzp_trace_fixed_w_sequence` is a deterministic diagnostic corrector at a caller-frozen monotone
sequence of `w=z-K0` targets. The augmented residual contains the unchanged ordered equations and
one scaled fixed-`w` row. Residual acceptance is subordinate to finite-system, bordered-rank, event,
and independent A--D gates; in particular, a small-residual point cannot bypass a failed rank gate.
The first Newton update may use one bounded full launch when its scaled RMS correction is already
inside the registered launch tube, after which the ordinary decreasing-residual line search applies.
This is a numerical reparameterization probe, not a continuous branch enclosure.

`invzp_scalar_event_interval_oracle` evaluates outward-rounded absolute pole and reciprocal-mean
bounds on a caller-supplied `w` interval. It can emit a passing graph event only when two independent
inputs are validated: a complete branch enclosure over the interval and a residual-derived
normalization enclosure bounding `|z|`. A bare caller number is not accepted as the latter.
`invzp_fixed_w_event_evidence` deliberately supplies neither proof: it matches the sampled fixed-`w`
and fixed-`h` records, exports conditional scalar bounds, and sets every
`event_bracket_ok=false`. Discrete targets cannot prove existence, uniqueness, nonsingularity, or
event exclusion between samples.

`invzp_enumerate_ordered_roots` is the retained fixed-`h` cold-seed driver. The caller supplies every
seed explicitly. The driver solves in lexical seed-ID order, retains submitted order and every
attempt, accepts only Newton plus independent A--D acceptance, and passes the accepted independent
state coordinates `[Sigma;K0]` to `invzp_cluster_ordered_full_states`. These are the complete
independent unknowns; `K` and `lambda` are derived. The cluster oracle uses a scaled RMS state
distance and an explicit three-way relation: same-root evidence, distinct-root evidence, or an
unresolved tolerance band. It exports actual-root medoids only when every connected tolerance
component is internally same-root evidence; an ambiguous component is reported, never split by input
order.

`invzp_trace_fixed_h_sequence` is a deliberately non-adaptive natural-parameter handoff primitive.
It follows a caller-frozen strictly monotone field sequence with initial, copy, and then secant
predictors. Every promoted state must pass Newton, the independent A--D audit, and a frozen
full-state predictor tube. It retains every scheduled point, including the first failure and all
trailing `not_run_after_stop` records. This is useful for cheap regular segments on either side of a
fold; it neither retries nor crosses a fixed-`h` rank loss.

`invzp_fixed_h_trace_graph_inputs` is the strict adapter from that retained trace schema to the
audited graph's input schema. It rechecks each emitted vertex against declared residual, A--D,
domain, pole, and mean gates and labels whether those thresholds came from the Newton record's
embedded normalized configuration or from an explicitly allowed legacy declaration. A passing edge
event certificate must be bound to a preregistered method, source trace, record pair, ordered margin
list, positive full-edge lower bounds, and nonempty proof payload; each lower bound is also checked
against both stored endpoint margins. Missing evidence sets `event_bracket_ok=false`. The adapter
never derives a continuous certificate from pointwise-positive nodes. The diagnostic Newton helper
now embeds its complete normalized option record so new traces can require
`require_embedded_config=true`.

`invzp_assemble_audited_graph` is the first evidence-only connectivity layer. It accepts only
caller-supplied, strictly accepted vertices and explicitly certified adjacent trace records. Vertex
certificates preserve the raw residual threshold, A--D payload, domain result, signed event margins,
source schema, and the continuity state `[Sigma;K0;m]`. An adjacency becomes an edge only when its
registered state-distance, predictor-tube, signed-event-bracket, and any required tangent-line gates
all pass. Tangent evidence is unconditionally required for equal-`h` adjacencies. The assembler
sorts IDs, exposes directed evidence and undirected connected components, marks nonmonotone fold
segments without splitting them, and retains failed adjacency/termination evidence as blockers.
Termination can be resolved only when its status belongs to the preregistered successful-status
set; `budget_exhausted`, events, rank loss, derivative drift, and unresolved `q=0` therefore remain
nonconnections. `status='ok'` means only that the supplied evidence assembled without blockers and
may still contain multiple disconnected components;
`complete_candidate_claim` is always false. The function never solves, clusters, merges nearby
roots, matches endpoints across traces, bridges a gap, constructs a Jensen section, or invokes the
selector.

`invzp_ordered_squared_field_problem` is a positive-field endpoint adapter with
`q=(h/h_reference)^2`. In the retained endpoint experiment,
`h_reference=6.890625e-9 meV`, the last accepted h-coordinate handoff, so the handoff is exactly
`q=1`. The reference is a coordinate scale only. For `q>0` the adapter is an exact
reparameterization of the original node equations;
it never averages positive and negative fields. Every positive point uses a centred
Richardson derivative with actual half-width `min(q_fd_step,q/2)`, a mandatory derivative-drift
margin, and a representability check on all stencil nodes. A fourth-order forward value is retained
only for rejected `q=0` diagnostics. A `q=0` point that otherwise evaluates normally is rejected as
`q_zero_endpoint_unresolved`: a finite forward stencil cannot prove that the residual contains no
`sqrt(q)` term. Thus the adapter can condition and trace a positive-q approach, but cannot certify
the endpoint. Its `q_domain` bounds central continuation points, like the h adapter's `h_domain`;
positive derivative-stencil nodes may lie outside a finite reporting bound and remain explicit
pointwise evaluations of the same equations.

The adapter also exposes a default-off, one-shot static-closure polish. It is eligible only when all
`Sigma` equations already pass, proposes only the analytic scalar Newton change in physical `K0`,
and rejects changes above `static_polish_max_ulps` (default 4096 physical-`K0` ULPs). The proposal
counts against the corrector budget and is then recomputed through the complete q adapter. It cannot
bypass the unchanged residual, arclength, event, A--D, rank, tangent, or correction-size gates.

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

The exact `spec.schema='invzp_smooth_r_selector/v2'`, scalar `spec.version`, exact dimensionless
`spec.x`, absolute/relative tolerances `spec.tie.{slope,curvature,state,qcp}_{abs,rel}`, and
nonnegative constant-path floors
`spec.degenerate.{derivative_energy_max,curvature_energy_max}` are mandatory. A mismatched mesh,
duplicate ID, nonfinite primary datum, nonlogical certificate, or invalid endpoint raises
`invzp:SmoothSelector:InvalidInput`; it is not converted into a physical branch status.

The primary metrics are exactly

```matlab
derivative_energy = trapz(spec.x,abs(path.drdx).^2);
normalized_max_slope_change = ...
    max(abs(diff(path.drdx)))/sqrt(derivative_energy);
normalized_r_curvature = ...
    trapz(spec.x,abs(path.d2rdx2).^2)/derivative_energy;
```

for a resolved nonconstant path. If derivative energy and raw curvature both lie below their
preregistered compatible floors, the path is treated as constant and both normalized scores are
zero. Low derivative energy with curvature above its floor is unresolved and cannot be selected.
The raw numerator metrics are exported for audit but never rank a path. The normalized scores are
invariant under `r -> a*r+b`, `a != 0`, and remove the candidate-endpoint/amplitude bias of the v1
raw metrics. Lower values are retained within

```matlab
abs(a-b) <= abs_tol + rel_tol*max(abs(a),abs(b)).
```

The order is binding: admissibility, jump rejection, and mandatory forward/reverse/cold-seed
agreement; normalized maximum slope change; normalized curvature; maximum state curvature; then QCP
distance. The function returns
`selected`, `branch_ambiguous`, or `branch_unresolved`, together with all metrics, rejected reasons,
and the survivors after every stage. Missing state/QCP tie information causes ambiguity; missing or
disagreeing trace reconstruction fails the mandatory reproducibility gate and cannot be rescued by
a tie breaker. Input order and path IDs never break a tie.

## Deliberate omissions

The cold-seed root enumerator, a one-direction fixed-`h` handoff primitive, a strict fixed-`h`
trace-input adapter, and an evidence-only graph assembler now exist. Complete bidirectional branch
reconstruction, continuous signed event-crossing production, arclength trace-to-graph conversion,
cross-trace endpoint matching, fold splitting, single-valued Jensen-section construction, endpoint
reconstruction, QCP field sequence, physical-`h` diagnostic metrics, mesh-refinement drift, and
cross-representation winner audit still have to be built and frozen before the global 1.5 T
path-selection pilot. In particular, a unique result from the selector means only unique under the
supplied normalized-shape spec. It does not authorize production use unless the omitted
physical-`h`, representation, and refinement audits agree under the registered grid and Matsubara
refinements.

Fresh in-memory MATLAB fixtures cover endpoint/amplitude invariance for identical normalized shapes,
a linear-over-curved winner, constant-path handling, unresolved low-energy/nonzero-curvature data,
mandatory missing/disagreeing trace rejection, the state-curvature tie break, exact ambiguity, and
v1/invalid-input rejection. `checkcode` is clean. No persistent test script is kept in the worktree.

A separate fresh graph fixture covers deterministic input ordering, a connected chain, an isolated
accepted root, a missing signed-event bracket, unresolved termination, a nonmonotone fold, a missing
required endpoint, equal-`h` tangent enforcement, strict residual equality rejection, and negative
`h` rejection. It also rejects a termination whose `resolved` flag contradicts the registered
successful-status set. No persistent graph test script is kept.

A fresh fixed-`h` adapter guard translates the retained 33-node root-6 trace into 33 vertices and,
with no signed-event certificates supplied, exactly 32 blocked adjacencies. Synthetic interval
certificate fixtures unlock exactly those 32 edges; an arbitrary nonempty payload, a source/pair
mismatch, missing embedded threshold provenance under strict mode, and an embedded `tol_outer`
mismatch all fail closed. The normalized Newton configuration is also checked on a fresh zero-field
correction. No persistent adapter test script is kept.

A fresh fixed-`w` fail-closed guard matches the same 33 fixed-`h` nodes to 32 positive fixed-`w`
targets and the separate zero-field anchor. All 32 scalar intervals have conditional absolute
clearance, but all 32 event records remain unvalidated and the graph unlocks exactly zero edges.
NaN corrector residual, NaN bordered `rcond`, NaN anchor data, a failed nested A--D audit, and an
initial-point provenance mismatch all reject the input. A real one-target replay with residual
already below tolerance and an intentionally raised rank floor stops at `rank_gate`, confirming that
small residual cannot bypass local regularity. `checkcode` is clean and no test script is retained.

## Verification record

The generic tracer crossed both folds of

```matlab
R(u,h) = u^3-u+h
```

and recovered the exact fold fields `h=+/-2/(3*sqrt(3))` within the declared arc discretization.
Fresh guards also verified reverse orientation, retry retention, invalid-tangent rejection, and
distinct terminal statuses for callback and nonfinite-residual events. A deliberately unsolvable
post-initial curve additionally freezes the retry envelope: two permitted attempts of three equation
evaluations each are both retained, then terminate as
`budget_exhausted/step_attempt_budget:evaluation_budget`. The cubic fold oracle remains unchanged
under the default envelope. Replaying the retained low-trace first step reproduced its state,
tangent, and record exactly: the eighth attempt was accepted and no attempt used more than four of
the permitted 16 equation evaluations.

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
fold-coordinate status           not certified
```

Every corrected root passed the independent A--D audit, and the oriented field tangent changes sign
with a regular bordered corrector, so the fold topology is established near
`h≈5.243548e-3 meV`. The two local interpolations are not a fold-coordinate oracle: the dense
Hermite calculation returned `0.0052435482911986821 meV`, which lies above an accepted state at
`0.005243548157786061 meV` and therefore cannot be a certified minimum. The earlier quarter-step
calculation returned `0.005243548147122800 meV`; neither estimate has a retained enclosure or
uncertainty. They are preserved here as rejected, method-dependent outputs, not as 17-digit fold
claims. A future coordinate result must retain the adjacent tangent-sign records, a fold bracket
consistent with every accepted node, and a separately labelled interpolation estimate with error
control. The dense/Richardson oracle took 1292.9 s over six bounded segments.
Subsequent profiling corrected the initial bottleneck attribution: at the fold, dense solve and
`rcond` take only `0.00745 s` and `0.00670 s`, while one Richardson equation/Jacobian evaluation takes
`1.771 s`. The exact factor solve takes `0.00290 s`, too small a saving to justify replacing the
dense oracle now. Exact last-point caching instead reduces a repeated evaluation from `1.748 s` to
`6.2e-5 s` and shortened a fresh three-step fold trace to `5.67 s/step` without changing arithmetic.
The next optimization target is redundant local node/field-derivative evaluation, not linear
algebra. Temporary `.mat` traces stayed under `/tmp` and are not retained.

### QCP-anchored simultaneous-audit pilot

A 4.05 T trace exposed an audit defect inherited by earlier continuation experiments. The
simultaneous unknowns are `[Sigma(:);K0]`, but nested Block A reran
`invz_emt_static_ordered` and could select another static root. Optional
`formulation='coupled'` now repeats the simultaneous Sigma and defactored-static equations
independently, and runs no nested static solver unless `debug_legacy_nested=true`.

At seven retained states, coupled Blocks A/B matched the corresponding components of
`invz_ordered_node_equations` with zero reported difference. One-ULP corruptions of derived
`lambda` and static `K(1)` were rejected. Profiling found no `invz_emt_static_ordered` call on the
ordinary coupled path, while the opt-in legacy diagnostic did call it. Implicit and explicit nested
audits were `isequaln`.

The QCP-connected trace retained 126 consecutive accepted states from
`h=0.00411428661370878` to `0.01172805209203838 meV`; its largest coupled residual was
`9.758e-10`. At `h=0.01171732294388717 meV`, continuation and natural-`h` Newton both produced
accepted roots, but with `r=0.768127507` and `0.822169537` despite scaled state distance
`0.00103702`. This is a direct branch-switch oracle: iteration count or a loose proximity tube
cannot decide the path. Full `.mat` traces remain temporary under `/tmp`.

A broader cold census at the frozen `h=0` node used 25 explicit product seeds:
five constant offsets of the seed `Sigma` crossed with five values of `K0/Jscale`. Sixteen attempts
were A--D accepted and clustered into **seven** distinct `[Sigma;K0]` roots; nine attempts failed.
Their representative `r` values span `0.60295` to `1.20283`. Thus the earlier agreement of three
cold seeds was basin-local evidence, not a uniqueness result.

The low-branch fold was then reconstructed with the staged policy intended for development economy:
21 cheap fixed-`h` nodes reached the fold neighborhood, and 20 local pseudo-arclength steps crossed
it with an oriented-tangent sign change and regular bordered corrector. This establishes a fold near
`h≈1.1304e-5 meV`, but not a certified coordinate. The short-trace interpolation returned
`1.1303917626651595e-5 meV`, below the accepted node at
`1.130400753628218e-5 meV`, so it cannot be a certified maximum. The earlier one-off interpolation
returned `1.130402502867847e-5 meV`; it too lacks a retained enclosure and uncertainty. Both raw
outputs are therefore method-dependent diagnostics, not bounds or precision claims. On the returning
leg, the raw closure audit became ill-conditioned while
the defactored equation remained accurate. At `h=1.72265625e-7 meV`,
`|Gbar-G|=1.14201498036e-7` narrowly failed its raw gate, whereas
`|K0-Jloc|/Jscale=1.93307590536e-11` and the direct vector residual was
`1.93299895588e-11`. This is the exact finite, nonzero-`G` identity
`Gbar-G = G*Gbar*(K0-Jloc)` exposing coordinate amplification by the large susceptibilities.

The optional `eso.audit_coordinate='defactored'` audit carried that returning segment to
`h=6.890625e-9 meV`, where the fixed-`h` Jacobian had `rcond=1.33e-11`. It did **not** establish a
connection to any of the seven `h=0` roots: fixed-`h` Newton then lost rank, and direct or polynomial
endpoint extrapolation failed the residual. The low endpoint topology therefore remains unresolved.
No production solver or Jensen path is authorized by this pilot. The bounded retry envelope is
intended to make any future endpoint trace fail promptly and reproducibly instead of consuming an
open-ended nested retry cost; it does not regularize the endpoint or change an acceptance gate.

At the last accepted positive-h point, the ordinary h-coordinate border has `rcond≈1.76e-11` while
its best tested Richardson column still drifts by `1.30e-3`; a bounded corrector stagnates at
residual `1.91e-5` and fails line search. Reparameterizing only the positive side by q first advanced
three A--D-audited states from `q=1` to `q=0.775` (`h=6.06609862186847e-9 meV`). The attempted next
point at `q=0.7625` then stalled at `1.3354847848219768e-8`, entirely in component 741, the final
static closure. Static-row scaling improved the bordered `rcond` but reached the identical floor and
was discarded. Direct, reversed, and compensated static sums agreed to about `1e-16`, excluding the
final Brillouin-zone summation order as the cause.

An ULP scan showed a quantized coupled-Newton direction: the static residual jumped from
`+1.33548e-8` to about `-1.96e-8`. Holding the already-converged `Sigma` and q fixed, however, a
656-physical-ULP `K0` proposal reduced the complete residual to `1.13e-11` and passed every A--D,
constraint, event, derivative-drift, rank, tangent, and correction-size gate. The default-off
one-shot polish freezes exactly this limited remedy. With it and the adaptive centred q stencil, the
same returning component reached `q=0.062004970571964718`
(`h=1.7158205633221428e-9 meV`) with residual `6.11e-10`, derivative drift `1.99e-3`, and bordered
`rcond=1.05e-12`. The next segment stopped on the derivative-drift event. At the last accepted state,
independently acceptable stencil caps `0.1` and `0.01` produced q columns differing by about 6%;
smaller caps drifted by up to 4.5%. This is a recorded double-precision derivative frontier, not an
endpoint or a connection to any enumerated h=0 root. At one enumerated h=0 root the forward stencil
is finite and has drift `1.13e-9`, but the independent `q_zero_endpoint_unresolved` event still
rejects it; small stencil drift is numerical evidence, not the missing smoothness proof.

A derivative-free fixed-`h` follow-up separates the remaining raw condition estimate from true rank.
At the positive-q handoff the static Jacobian row is `3.53e6` times larger than every Sigma row:
physical-coordinate raw `rcond=8.27e-13`, but row-equilibrated `rcond=6.10e-9`. With equilibration
used only for the linear solve/rank gate and the same raw residual/A--D acceptance, the retained
fixed-`h` sequence accepts 52 scheduled points down to
`h=3.4852605192481022e-10 meV`. Refining that final interval reaches
`h=3.4035747258282251e-10 meV` and then fails line search at raw residual `3.83e-6`. A direct `h=0`
solve from the last accepted state fails at `5.077e-2`, and the last state remains at frozen scaled
distance at least `0.312` from all seven enumerated `h=0` roots. Thus row equilibration removes a
false rank alarm but does not establish the missing endpoint connection.

A simple-first follow-up then screened every enumerated zero-field root upward before attempting
another global construction. The raw closure again proved to be the wrong numerical coordinate:
at a retained root-4 point, the accepted defactored vector residual was `2.78e-10`, while the raw
closure residual `2.39e-8` exceeded its `1.14e-8` gate. Switching the audit to the already-derived
exact defactored diagnostic coordinate removed every audit failure in the retained root-4 and
root-5 traces. This does not replace the frozen production/raw certification. Root 4 crossed a fold
between
`h=5.5677706301911396e-6` and `5.568143040810256e-6 meV`, with tangent-parameter sign change and
minimum bordered `rcond=4.35e-10`; it later stopped at a static-component line-search floor, not a
domain event. Root 5 accepted 100 local steps without an event or audit failure but advanced only to
`h=3.305411778307727e-6 meV`; that bounded trace does not locate its fold.

Root 6 gives the more decisive topology check. A strict raw-coordinate trace crossed a fold near
`h≈9.54114e-6 meV`; its raw turn estimate was `9.541141092330816e-6 meV`, but no retained enclosure
or uncertainty promotes that output to a certified coordinate. Its returning leg crossed zero field
with no audit failure and minimum bordered
`rcond=4.36e-7`. A local constant-step refinement bracketed zero by
`+5.270703260674207e-10` and `-1.492816168710522e-10 meV`. The interpolated seed needed a correction
of only `7.43e-12` in the frozen clustering metric, and the accepted zero-field solve matched census
root 6 within `1.42e-14`; both are inside the registered `1e-9` same-root threshold. Thus this
returning leg has a certified endpoint. A subsequent derivative-free check followed the opposite
leg on 33 decreasing fixed-`h` nodes from the common handoff to zero. All 33 nodes passed A--D with
no event failure, the largest frozen predictor distance was `0.0247702`, and the endpoint matched
root 6 at cluster distance `7.92e-14` (residual `7.52e-14`). Both observed root-6 legs therefore
have the same zero-field endpoint evidence. This fixed-`h` sequence is not by itself a continuous
signed-edge certificate, and it does not resolve the previously traced root-7 boundary layer or a
complete single-valued Jensen section.

The same opposite-leg data were replayed in the scalar coordinate `w=z-K0`. All 32 positive targets
were accepted in at most four Newton iterations:

```text
maximum residual_inf              9.4914e-11
minimum bordered rcond            1.6065e-9
maximum predictor correction RMS  3.4e-3
maximum fixed-h state mismatch    5.8624e-11
maximum h mismatch                1.6466e-14 meV
maximum w mismatch                2.7439e-15 meV
zero-field anchor residual        7.5218e-14
```

This is strong evidence that fixed `w` is a better-conditioned local solver coordinate on this leg.
It is not “32 certified edges.” The samples do not enclose the intervening root manifold, and the
normalized pole margin additionally needs a full-edge bound on the closure residual/`z`
normalization. The retained report therefore contains 32 conditional scalar clearances, zero
passing graph-event certificates, and zero unlocked adjacencies. A genuine proof would need a
validated interval-Newton/Krawczyk (or equivalent) cover over each target slab.

An h-adapter copy of the existing ULP-bounded q-coordinate `K0` polish was tested as an even simpler
remedy at the root-4 and root-5 frontiers. Both proposals correctly rejected themselves as outside
the 4096-ULP envelope; ordinary step shrink, not the polish, accepted the next points. The
unnecessary h-adapter extension was therefore removed. Future long diagnostics are one root per
process and save the trace before summaries; a two-root atomic batch consumed about 95 minutes and
nearly discarded both results without adding an acceptance guarantee.
