# Ordered \(1/z\) convergence: final diagnosis and production resolution

## Status

The convergence investigation is complete for the current production
objective. At \(T=0.1\) K, the factor-one missing-area member now produces a
smooth, finite susceptibility across the full 101-point
`linspace(0,9,101)` field grid. This was confirmed by the full production
sweep after the exact-zero and ordered-boundary repairs.

The result is deliberately labeled. It is the opt-in controlled missing-area
approximation, not a proof that the strict equation-(45) path exists at every
integration coordinate. Option-free `full_profile` remains fail-closed. The
current driver deliberately sets `missingAreaFactors = 1.0`; factors 0.75 and
1.5 remain separately tested sensitivity choices and can be incomplete where
factor one is accepted.

## Governing equation and the source of the difficulty

The ordered-field equation is

\[
H_0(h)=\int_0^h r(s)\,ds,
\qquad
r(h)=\frac{G_0(0;h)}{\widetilde G_0(0;h)},
\]

\[
\boxed{F(h)=H_0(h)-J_{0,\mathrm{eff}}m(h)=0},
\qquad h\equiv H_{\mathrm{MF}}.
\]

This is not one scalar root problem. It has three nested levels:

| Level | Numerical object | Acceptance requirement |
|---|---|---|
| Static closure | \(\widetilde G_0(0)\), \(K_0\) | bounded admissible root with positive physical masses |
| Coupled node | \(\Sigma,K,\lambda\) at fixed \(h\) | converged self-energy and static closure with no failed-state commit |
| Ordered field | the profile integral and \(F(h)=0\) | one declared component, justified lower area, one refined root |

The original implementation treated failures at these levels too much like a
single Picard-convergence problem. The investigation showed that they have
different causes and require different remedies.

## Verified root causes

### 1. Strict profiles do not provide a complete accepted path from \(h=0\)

The strict construction needs a valid endpoint and every positive-\(h\) node
joining it to the final root. In the problematic ordered regime, the low-\(h\)
nodes are outside the accepted coupled/static component while a contiguous
high-\(h\) suffix is well behaved. Consequently, `full_profile` cannot define
the complete integral and correctly returns `node_failed`.

This is the principal reason that increasing `nH`, changing damping, or
raising `max_outer` did not repair the blanket mask: those changes do not
create the missing lower component or its integration constant.

### 2. Failed profile nodes formerly contaminated later warm starts

The old profile sweep could return partially updated \(\Sigma\) and \(K_0\)
from a rejected node and use them to seed the next one. That made verdicts
depend on grid history. The solver now commits state transactionally: only a
fully accepted node becomes the next continuation carrier.

After this repair, shared nodes on nested grids agree. The change removed a
real numerical defect, but it did not by itself create the missing low-\(h\)
path.

### 3. The usable solution is component- and direction-dependent

The coupled equations contain folds and multiple admissible sheets. A cold
ascending sweep can reach a short terminal component, while a descending
sweep seeded from an accepted ordered field reaches a longer component. A
failed transfer therefore does not prove that the coupled root is absent, and
one successful warm start does not by itself select equilibrium.

This distinction caused both visible classes of interior holes:

- the 0.36 and 0.45 T masks were ordinary continuation-path misses between
  accepted ordered fields;
- the 4.68 T mask was the last ordered-side sample below the PM boundary, so
  it had no ordered source above and could not use the two-sided retry.

### 4. Exact zero field exposed a removable representation singularity

At \(B_x=0\), the lowest electronic doublet is degenerate and its numerical
basis is arbitrary. A basis-dependent transition weight can therefore be
exactly zero even though the invariant doublet response is finite. The old
arithmetic evaluated a structure equivalent to

\[
M^2\left(\frac{\text{finite}}{M^2}\right),
\]

which became `0*Inf = NaN`. This was not a non-unique final ordered state.

The exact \(M^2=0\) branch now evaluates the already-cancelled expression

\[
Q_0=\frac{2m^2}{n_{01}^2}
\left[\lambda_1-(1-n_{01}^2)K_0\right].
\]

Positive-\(M^2\) states preserve their historical arithmetic. The unresolved
exact endpoint is rejected without committing state, while the certified
positive-\(h\) component supplies the missing-area solution. The resulting
exact-zero spectrum is continuous to \(B_x=10^{-6}\) T.

### 5. The high-field sliver was a continuation-distance failure

The cold factor-one profile at 4.68 T follows a component whose algebraic
support does not reach factor one. Descending profiles seeded from untouched
4.50 and 4.59 T ordered states recover the same factor-one root:

\[
h_{\mathrm{MF}}=0.0036183\ \mathrm{meV},\quad
D_{\mathrm{uni}}=0.023913,\quad
F'=0.020014.
\]

The initial audit used large direct field jumps. One distant seed failed at
4.70 T and both failed at 4.71 T, which was initially treated as a stopping
condition. A fine-step continuation subsequently advanced two independent
histories at 129 and 257 profile nodes through 4.7188 T. Both histories reached
the same branch with positive margins and residuals below
\(6.3\times10^{-9}\). Thus the old 4.70/4.71 failures measured the reach of a
distant warm start, not termination of the ordered branch.

The derivative check was also corrected. The exact identity is

\[
F'(h)=r(h)\left[1+J_{0,\mathrm{eff}}\widetilde G_0(h)\right]
     =r(h)+J_{0,\mathrm{eff}}G_0^{\mathrm{bare}}(h),
\]

not \(r(1+J_{0,\mathrm{eff}}G_{\mathrm{stat}})\). Substituting
\(G_{\mathrm{stat}}\) produced the former false negative sign.

## Implemented production solution

### Controlled missing-area completion

Let \(h_e\) be the lower edge of the terminal contiguous certified component.
Only nodes at and above \(h_e\) enter quadrature. The unresolved contribution
is represented explicitly as

\[
A_f=f\,h_e r(h_e),
\qquad
H_0(h)=A_f+\int_{h_e}^{h}r(s)\,ds.
\]

For the linear completion used to interpret the ensemble,

\[
\frac{r(0)}{r(h_e)}=2f-1.
\]

The tested sensitivity factors are:

| Factor \(f\) | Implied linear endpoint ratio | Role |
|---:|---:|---|
| 0.75 | 0.5 | lower-area sensitivity member |
| 1.00 | 1.0 | current production member |
| 1.50 | 2.0 | upper-area sensitivity member |

These factors are not probabilities, fitted exchange constants, or confidence
limits. They vary only the missing lower integral. The full three-factor
ensemble is a diagnostic option; the successful 101-point run used only the
factor-one member because the noncentral high-field members do not satisfy the
boundary-retry contract. The declared numerical branch is
`picard_attracting_contiguous_high_h_component`.

### Frozen two-sided retry for interior holes

The map freezes the original cold-pass labels and states. A masked point
between two accepted ordered fields can be retried from each side. Admission
requires both direct results to pass all original gates and agree in
\(h_{\mathrm{MF}}\), self-energy, and component edge. Recovered points never
seed another recovery, so the output is traversal-order independent.

This policy recovers the 0.36 and 0.45 T masks for every area member at the
production 129-node resolution.

### Frozen ordered-boundary retry

The last ordered-side hole uses a distinct policy. A target is eligible only
when all of the following hold:

1. its frozen cold label is `phase==0`, never `phase==2`;
2. it lies above every accepted ordered cold point and below an accepted PM
   cold point;
3. two untouched lower ordered sources lie within the declared 0.20 T span;
4. an independent PM solve converges at the target with negative PM mass;
5. both descending target solves find one unbridged root on the same component;
6. the solutions agree and retain \(D_{\mathrm{uni}},F'\ge10^{-3}\), positive
   physical masses, and the original final-residual gate; and
7. the real-axis response is finite.

On the current grid, only factor one satisfies the complete contract at
4.68 T. Factor 1.5 lies outside the component's support, and the noncentral
members do not have the required source pair. The central spectrum is therefore
produced while the sensitivity interval remains honestly incomplete.

## Hypotheses that were ruled out

| Hypothesis | Verdict |
|---|---|
| “The plotted zeros are a real-axis or plotting problem.” | False. Masking occurred before `invz_chi_realaxis`; accepted states produce finite spectra. |
| “More profile nodes will repair the solver.” | False as a general fix. Nested grids diagnose edges but do not create a missing component. |
| “More damping or iterations is sufficient.” | False. Some states are noncontractive and some cold profiles lack algebraic root support. |
| “Every `no_admissible_static_root` proves no coupled root exists.” | False. It describes one trial state; continuation and the reduced residual find other coupled roots. |
| “Discarding failed/NaN nodes makes equation (45) valid.” | False. Filtering alone loses the lower area and can bridge components without a selector. |
| “Small \(n_{01}\) or small \(M^2\) explains the general field-range masks.” | False. \(n_{01}\) is near one at the measured failed nodes. Exact \(M^2=0\) matters only at the symmetry endpoint and has a removable algebraic limit. |
| “The electronic/electronuclear representation mismatch caused the masks.” | Not established. It limits exact thermodynamic interpretation but does not explain the observed continuation topology. |
| “The branch ends at 4.70–4.71 T.” | False on current evidence. Fine-step two-history continuation reaches 4.7188 T. |
| “Use \(G_{\mathrm{stat}}\) in the \(F'\) identity.” | False. The identity requires \(\widetilde G_0\) or equivalently \(G_0^{\mathrm{bare}}\). |
| “Tune the area factor until the plot is smooth.” | Rejected. Factors are a declared sensitivity ensemble; unsupported members remain incomplete. |

## Strict result versus production approximation

| Property | Strict `full_profile` | Opt-in production missing-area route |
|---|---|---|
| Lower path | requires every node from exact \(h=0\) | represents unresolved lower area explicitly |
| Failed nodes | any required failure masks the point | never integrated; only the certified terminal component is used |
| Branch selection | complete-path requirement | declared Picard-attracting high-\(h\) component |
| Field retries | none | frozen, independently seeded, fail-closed policies |
| Exact zero | remains `node_failed` | finite symmetry-limit solution |
| Current 101-grid central spectrum | not generally available | smooth and finite across the full field range |
| Thermodynamic claim | rigorous target, presently incomplete | controlled numerical approximation, not an equilibrium proof |

## Thermodynamic limitation that remains

The successful susceptibility does not retroactively prove a global Jensen
free-energy construction. The production hybrid combines the moment and bare
response of the full 136-state electronuclear ion with an electronic
two-level vertex. The fold-anchored sheet integrals and exact reduced residual
remain valuable branch diagnostics, but the same-ion stationary functional
required for an exact equilibrium ranking has not been derived for this hybrid.

This limitation is separate from numerical convergence. It does not invalidate
the labeled missing-area spectrum, and it must not be cited as the root cause
of the former masks.

## Production locations

- Driver configuration: `invz_projected/invz_run_spectra.m`.
- Map-level ensemble and retry policies: `invz_projected/invz_spectra_map.m`.
- Ordered profile and transactional node solver:
  `invz_projected/invz_solve_point_ordered.m`.
- Exact-zero Matsubara arithmetic: `invz_common/invz_sigma_ordered.m`.
- Exact-zero real-axis arithmetic: `invz_projected/invz_chi_realaxis.m`.
- Missing-area quadrature: `invz_common/invz_missing_area_integral.m`.

## Verification and evidence

The final checks include:

- the user-confirmed full 101-point production sweep;
- `invzp_ordered_boundary_retry_smoke`;
- `invzp_adjacent_retry_map_smoke` and
  `invzp_adjacent_retry_highfield_smoke`;
- `invzp_approximation_production_validation`, including
  `strict profile match 1`;
- `test_invzp_zero_field_failclosed`;
- `test_invzp_hmf_derivative_identity`;
- `test_invzp_missing_area_integral`;
- `test_invzp_reduced_residual`;
- all 13 `test_invzp_static_domain` gates; and
- MATLAB Code Analyzer and `git diff --check`.

The compact evidence index is
`docs/diagnostics/invzp_approximation_wp6/README.md`. The chronological record
is `docs/execution/invzp_convergence_journal.md`, especially checkpoints
26–31. Rejected explanations and their reconsideration conditions are retained
in `invzp_convergence_dead_ends.md`.
