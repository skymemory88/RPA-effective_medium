# Projected-spin QCP diagnosis and recommended correction

Date: 2026-07-21

> **Status banner added 2026-07-27 — the body below is unchanged and remains the authoritative
> record of the 2026-07-21 diagnosis. Three of its point-in-time observations have since been
> overtaken by work it recommended:**
>
> - *"The ordered elastic static sector is NOT implemented … J 2.28–2.29 exists nowhere."* Implemented
>   2026-07-22 (Stage 2, commits `1259ffc`, `1c8d64d`, `5fc7d6e`), together with the H_MF
>   self-consistency (J 2.31–2.33) and the free-energy route (J 2.34). Closure achieved at
>   `Bc_1z = 3.025 T` vs the PM `crit`-zero at `3.033 T`.
> - *"The fast suite is NOT currently green"* — the nominally-zero off-diagonal tolerance bug was
>   repaired in `3244f12` (AbsTol floor).
> - *Tensor-path handoff status* — since superseded by the committed stability-based dispatcher.
>
> The central analysis — that the non-closure follows from comparing two approximations with
> different static pole conditions, and must not be forced to close — was acted on and **stands**.
> Current project state: `docs/invzp_ordered_1z_state.md`. Execution history:
> `docs/INVZ-DEVELOPMENT-RECORD.md`.

Provenance: original diagnosis by Codex; revised the same day after code verification, then corrected again following the Codex plan review (`invzp_QCP_plan_review_Codex.md`), whose findings were independently re-verified: the free-energy section (framework §9.4) is restored to the citation record, the HEAD-vs-working-tree equation renumbering is documented, the `invz_sigma_ordered` implementation scope is corrected, and the tensor handoff status is downgraded from "production" to working-tree.

Status: the experimental projected-spin phase-handoff change has been reverted. The existing ordered-first behavior and its known finite-gap/discontinuity limitation are restored. The full-tensor logarithmic plotting change is intentionally retained.

## Verification record (2026-07-21, updated after the plan review)

- Denominator conventions match the analysis: `invz_chi_realaxis.m:8` implements `chi~0 = chi0/(1+Sigma)`, `chi = chi~0/(1 - J_nu chi~0)`, and the static 1/z mass `pt.crit = 1 + Sigma0 - J0eff*chi0cc0` is computed at `invz_solve_point.m:254` and `invz_solve_point_ordered.m:146`. The real-axis continuation seeds `K` with its static Matsubara value (`invz_chi_realaxis.m:49`), so quantitative real-axis peak positions carry a static-K approximation on top of the pole algebra.
- The shared dispatcher is as described: `invz_spectra_map.m` (`one_field`) selects one state via the ordered-first `invz_solve_auto`, and the RPA overlay is generated from that same state via `invz_zero_sigma_overlay`, sharing `chi0cc_w`.
- The bare-MF ordered-side limitation is self-documented: `invz_solve_point_ordered.m` defers the applied-field/H<sub>MF</sub> self-consistency, so the projected moment onsets at the bare boundary.
- The ordered elastic static sector is NOT implemented: `invz_common/invz_sigma_ordered.m` covers only the ordered self-energy and `alpha_m` (J 2.26--2.27); the elastic single-site function \(G(0)\) with the \(\xi\) resummation (J 2.28--2.29) exists nowhere -- `invz_solve_point_ordered.m:127` feeds the plain full electronuclear propagator into the generic EMT path. The simple ordered `pt.crit` is therefore not the complete Jensen ordered inverse response.
- The revert is clean: no residue of the linearized handoff exists in `invz_projected/`; the pending working-tree diffs there are driver knobs and README text only.
- Tensor-path status: the `a3d` ordered solve is committed (2026-07-19/20), but the PM-first stability dispatcher and the `hmf_J0z` linearized handoff exist ONLY as uncommitted working-tree edits (`hmf_J0z` has zero occurrences at `HEAD` in `invzt_solve_auto.m` / `invzt_solve_point_ordered.m`). The handoff is also `a1`-only (`invzt:hmfMode`); direct `a3d` calls retain the bare-MF state.
- The retained tensor plotting change is as described: `invz_tensor/invzt_run_spectra.m` computes a local 99.5th-percentile log10 clim inline, with no dependency on the projected-only `invz_robust_pct` helper.
- `jensen_1z_framework.html` carries an uncommitted renumbering of the Section 9 equations (+2 relative to `HEAD`: elastic \(G(0)\)/\(\xi\) now 41--42, the H<sub>MF</sub> chain 43--45, the free-energy test 46). Citations below therefore use Jensen's original equation numbers (J 2.x) and section anchors, which are stable across HTML revisions.
- The fast suite is NOT currently green: `test_invz_chiperp/test_anchors_and_symmetry` fails on a relative-tolerance comparison of nominally-zero off-diagonal elements (~1e-15 against diagonals ~17.6) -- numerical-test fragility, not a QCP regression. Repair it before using the suite as a regression gate (Stage-1 plan Task 0).

The reported ~1 GHz gap and the matched-field mode-count inspections are run-derived numbers: consistent with the verified structure, but not re-derived statically.

## Executive conclusion

The projected-spin 1/z zero-frequency pole failing to close at the field where the projected RPA mode closes is not primarily a peak-finding or plotting bug. The non-closure follows directly from comparing two approximations with different static pole conditions (the finite *numerical* real-axis gap quoted later additionally rests on the continuation and peak-tracking qualifications in the Interpretation section):

\[
\chi_{\mathrm{RPA}}(\omega)
=\frac{\chi_0(\omega)}{1-J\chi_0(\omega)},
\qquad
\chi_{1/z}(\omega)
=\frac{\chi_0(\omega)}{1+\Sigma(\omega)-J\chi_0(\omega)}.
\]

At the RPA/MF critical field,

\[
1-J\chi_0(0)=0,
\]

so the 1/z static denominator at the same state and field is

\[
D_{1/z}(0)=\Sigma(0),
\]

which is generally nonzero. Conversely, at the 1/z critical field,

\[
1+\Sigma(0)-J\chi_0(0)=0,
\]

the RPA denominator is generally \(-\Sigma(0)\), not zero. Thus RPA and 1/z should not be forced to close at the same physical field unless the self-energy vanishes there or the model parameters have been separately recalibrated.

The projected implementation compounds this conceptual difference by using one phase decision and one single-ion state for both plotted responses. In `invz_spectra_map.m`, the RPA curve is made by zeroing the self-energy of the already selected 1/z point. That is convenient away from the transition, but it cannot make both approximations thermodynamically consistent when their critical fields differ.

## What the reverted experiment changed

The attempted correction did two coupled things:

1. It selected PM or FM from the paramagnetic 1/z mass, so the shared phase switch occurred at the renormalized 1/z boundary rather than at the bare MF/RPA boundary.
2. On the ordered side it replaced the bare ordering coupling by the boundary-linearized value

   \[
   J_{\mathrm{MF}}=\frac{J(0)}{1+\Sigma_{\mathrm{PM}}(0)},
   \]

   and handed the rejected PM self-energy to the ordered solve as a seed.

The linearized relation does make the infinitesimal ordered moment appear at the PM 1/z instability. It does not, however, supply a complete finite-moment ordered 1/z theory.

## Why that experiment broke the spectra

### 1. One phase dispatcher was used for two theories

The 1/z stability test selected the single-ion state, after which the RPA overlay was generated from that same state with `Sigma = 0`. Between the 1/z and RPA critical fields, the chosen state is appropriate for one approximation but not the other. The RPA FM and PM pieces therefore do not meet at the RPA theory's own transition. This is the direct cause of the newly broken RPA soft-mode branch.

There is no generally valid shared `S.phase` for both overlays. RPA and 1/z require separate phase labels, states, and critical fields.

### 2. The modified field was only boundary-linearized

The factor `J/(1+Sigma_PM(0))` is a linear boundary relation. Freezing it into a finite-moment FM solve omits the nonlinear applied-field-to-molecular-field relation and the response-consistency conditions required by Jensen's ordered 1/z construction (Section 9.3 of `jensen_1z_framework.html`; see Stage 2 below for the precise equations). It therefore has no controlled guarantee away from an infinitesimal neighbourhood of the continuous transition.

Near the QCP the ordered moment, static self-energy, and pole denominator are all small and strongly coupled. Solving a bare-looking nonzero-moment root with a frozen PM correction makes the result sensitive to root tolerances, field spacing, self-energy iteration, and continuation direction. That explains why the attempted FM branch became very noisy even though its limiting boundary equation looked plausible.

### 3. A colormap maximum is a fragile mode estimator near criticality

Close to a soft pole, a finite frequency grid and Lorentzian broadening can move the brightest pixel between neighbouring bins or weak electronuclear satellites. This can add visual noise even to a smoothly varying denominator. A future correction should track the pole (or the minimum singular value/eigenvalue of the inverse response) continuously in field and use the intensity map only for display.

### 4. Scope: the same handoff runs on the tensor path (working tree) -- and does not exhibit this contradiction there

The pattern condemned above is not wrong in itself. The tensor sibling runs it in its current working tree: `invzt_solve_auto` dispatches PM-first by stability (its PM `crit > 0` is the tensor QPT criterion; review P0-1, 2026-07-19), and `invzt_solve_point_ordered` accepts the boundary-exact linearized modified field `opts.hmf_J0z` from the rejected PM leg precisely so the ordered moment vanishes at the renormalized instability instead of the higher bare-MF boundary. (Both are uncommitted edits: `hmf_J0z` does not exist at `HEAD`; the committed tensor ordered machinery is the bare-MF `a1` solve and the `a3d` hybrid.)

It works there because `invzt_run_spectra` plots no `Sigma = 0` RPA overlay: only one theory ever asks the dispatcher for a state, so the shared-dispatcher contradiction of point 1 never arises. The failure diagnosed here is therefore specific to the projected driver's shared-state RPA overlay, not to PM-first dispatch or the linearized field as such. The away-from-boundary limitation of point 2 does apply equally to the tensor `a1` handoff; a complete ordered treatment does not yet exist on either path (Stage 2 below).

## Interpretation of the original projected result

The restored projected ordered solver uses the bare MF order parameter. It therefore remains ordered up to the bare MF/RPA boundary. At that endpoint the RPA mode closes, while the projected 1/z denominator retains `Sigma(0)` and hence a finite gap (about 1 GHz in the reported run). Two qualifications on that number: the nonzero static mass proves only that the zero-frequency pole does not close -- the 1 GHz figure is a run-derived real-axis peak position that additionally carries the static-K continuation approximation (`invz_chi_realaxis.m:49`) and the display broadening, and should be reproduced by pole tracking under frequency-mesh/broadening refinement (Stage 3) before being quoted. And pole-closure statements here refer to the intrinsic response: with demag on, the strict-uniform measured observable saturates instead of diverging at the transition. This is internally understandable but is not a complete ordered-side 1/z prediction at its own renormalized QCP.

The paramagnetic projected 1/z mass

\[
\texttt{crit}=1+\Sigma(0)-J(0)\chi_0^{cc}(0)
\]

still provides a meaningful PM-side estimate of the renormalized instability. The unresolved part is constructing the corresponding thermodynamically consistent ordered state below that boundary.

## Recommended implementation

### Stage 1: separate the RPA and 1/z state paths

Refactor the spectra calculation so the two overlays no longer share one selected point:

- `phase_rpa`: use `Sigma = 0`, the bare-MF ordered state below `Bc_rpa`, and the PM state above `Bc_rpa`.
- `phase_1z`: use the 1/z stability condition and, once available, the consistent ordered 1/z state below `Bc_1z`.
- Return distinct diagnostics such as `S.phase_rpa`, `S.phase_1z`, `S.Bc_rpa`, and `S.Bc_1z`.
- Do not require the RPA and 1/z modes to close at the same field. Require each theory's FM and PM branches to meet at that theory's own critical field.

This separation fixes the architectural inconsistency without pretending that the present ordered 1/z state is complete.

Least-rework interim realisation (adopted by the Stage-1 plan; re-review finding 1): the auto-solve's ordered-first selection flips at the bare-MF boundary, so below it the `Sigma = 0` overlay it feeds *approximates* the RPA state -- exactly when the ordered 1/z EMT loop converges all the way up to that boundary. That is an approximation with known failure modes, not the required separation: the auto phase depends on 1/z convergence, so a failed ordered solve can mislabel a column PM or masked even though the bare-MF/RPA state exists. Stage 1 therefore ships the 1/z *dispatch and masking* (stability-gated PM state above `Bc_1z`; bare-MF ordered curve below it as a labelled diagnostic -- the complete ordered 1/z state remains deferred; columns with no accepted auto state or a suspect one masked out of BOTH overlays' spectra and peak curves), while reporting the overlay's boundary honestly as `Bc_auto` -- a dispatch diagnostic, NOT `Bc_rpa`. Implementation plan: `docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md`.

What Stage 1 does NOT discharge: regression 4 (an RPA dispatcher independent of 1/z convergence) stays OPEN. The independent bare-MF/RPA evaluator built directly on `1 - J chi0(0)` is a small, self-contained follow-up task scheduled immediately after Stage 1 -- deliberately NOT bundled into the Stage-2 ordered thermodynamics, which is unrelated and much larger work (re-review finding 1).

### Stage 2: implement the full ordered 1/z thermodynamics

STATUS: implemented on the projected path (docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md).

The J 2.34 two-route comparison in the closed 2x2 model shows a ~13.7% same-retained-order static-elastic approximation residual (scale-free in the coupling; Jensen's own full-machinery published check achieved 2-3% on the physical HoF3 lattice at low T, with his stated high-T elastic caveat). The implemented J 2.26-2.29 chain is not shown defective; thermodynamic cross-route closure is carried as a KNOWN APPROXIMATION LIMITATION, scientifically unvalidated pending the published-benchmark (HoF3) discriminator -- a follow-up task. Exact limits, convergence, and leading-order scaling ARE enforced in test_invz_deltaF.

For the 1/z FM side, two missing pieces must BOTH be implemented (citations use Jensen's original equation numbers -- the HTML numbers moved in an uncommitted revision, see the verification record):

- the ordered elastic static sector, J 2.28--2.29 (framework Section 9.2): the zero-frequency single-site function \(G(0)\) with the \(\xi\) resummation, inserted into the effective-medium form at \(\omega_n = 0\). `invz_sigma_ordered` implements only the ordered self-energy and \(\alpha_m\) (J 2.26--2.27); until J 2.28--2.29 exist, the simple ordered `pt.crit` must not be treated as the complete Jensen ordered inverse response;
- the nonlinear applied-field/H<sub>MF</sub> self-consistency, J 2.31--2.33 (Section 9.3): the moment integral, the differential response relation, and the boxed applied-field integral \(H_0=\int_0^{H_{\rm MF}} G_0(0;H')/\tilde G_0(0;H')\,dH'\), with the free-energy consistency test J 2.34 (Section 9.4) as the independent validation criterion.

The ordered moment and self-energy must be solved together rather than inserting a PM boundary value of `Sigma(0)` into a rescaled bare-MF equation.

Citation history (corrected twice; the lesson is the policy): the original diagnosis cited "the applied-field integral in equations 41--43" and "the free-energy consistency condition in equation 44" -- correct against the committed `HEAD` revision of the framework HTML. A first revision of this document, checking only the working-tree HTML (which renumbers Section 9 by +2), wrongly declared the free-energy condition nonexistent and "corrected" the numbers to 43--45. Both errors came from citing volatile HTML equation numbers: cite J-numbers and section anchors instead. Stale HTML-number comments remain in code (`invz_solve_point_ordered.m:27` "eqs 41-43"; `invz_common/invz_sigma_ordered.m:2` "eqs 37-38" -- both `HEAD` numbering) and should be switched to J-numbers when touched (Stage-1 plan Task 3).

Decision (2026-07-22, user): the PROJECTED path owns the ordered completion -- the observable being fixed is the projected spectra map, which the tensor route cannot supply (never stack the two routes). Implementation plan: docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md. The boundary-linearized cosmetic handoff is rejected. For context, the tensor a3d machinery remains NOT a complete vehicle by itself: The tensor branch's committed `a3d` solve improves the ordered self-energy treatment, but it is NOT by itself a complete Stage-2 vehicle: it has no H<sub>MF</sub> thermodynamic outer loop (the linearized `hmf_J0z` field is `a1`-only and rejected for `a3d` via `invzt:hmfMode`; direct `a3d` calls retain the bare-MF state), and the tensor PM-first handoff itself is presently an uncommitted working-tree edit. Whichever route is chosen -- projected or tensor -- must implement the elastic static sector (J 2.28--2.29) and the H<sub>MF</sub> relation (J 2.31--2.33) with the free-energy test (J 2.34) on top of its self-energy machinery. The projected README forbids stacking the two routes; the tensor branch keeps its own scope.

The solver should use field continuation from a well-converged ordered point, but convergence must be checked independently of the seed. Near the boundary it should verify all of the following:

- the ordered moment approaches zero continuously;
- the ordered inverse 1/z response approaches zero at `Bc_1z`;
- forward and reverse field sweeps converge to the same continuous branch;
- the result is stable against field step, Matsubara cutoff, q-grid, mixing, and root tolerance;
- the pole position is stable as the display broadening and frequency mesh are refined.

Until this stage is complete, the PM-side 1/z boundary is the controlled result; the ordered projected 1/z spectrum near that boundary should be labelled diagnostic/approximate rather than made artificially continuous.

### Stage 3: track modes by poles and continuity

For quantitative soft-mode curves, solve for poles or minima of the inverse response and connect them by field continuity and spectral overlap. Keep the full `chi''` map as the observable, but do not use only the maximum-intensity frequency bin as the mode identity near the QCP. This is particularly important for the electro-nuclear satellites and for the noisy FM side.

## Full-tensor versus projected mode count

The earlier apparent difference in the number of electro-nuclear modes was largely a display dynamic-range issue, not evidence that the full-tensor `a1` susceptibility mathematically contains only one transition. Matched inspections at representative fields found multiple low-frequency peaks in both responses, but the weaker tensor satellites were roughly orders of magnitude below the dominant line and were hidden by the former linear colour scale. The retained logarithmic colour scale in `invz_tensor/invzt_run_spectra.m` exposes those weak modes while using a local percentile calculation, so tensor plotting no longer depends on a helper that exists only on the projected path.

This display correction is independent of the projected QCP phase-consistency problem and is therefore intentionally not reverted.

## Required regression tests for a future fix

1. RPA FM and PM poles meet at `Bc_rpa` using RPA-consistent states, evaluated on the intrinsic (demag-free) response -- the demag-corrected strict-uniform observable saturates rather than diverging.
2. 1/z FM and PM poles meet at `Bc_1z` using the full nonlinear ordered 1/z state. (Delivered by the Stage-2 plan on the projected path: enforced by test_invz_qcp_closure, INVZ_SLOW-gated, via the static inverse response -- pole-based, per Stage 3's argmax warning.)
3. `Bc_rpa` and `Bc_1z` are permitted to differ and are reported separately.
4. Zeroing `Sigma` in a 1/z-selected state is not accepted as the RPA phase dispatcher near either boundary.
5. Pole trajectories are invariant, within tolerance, under forward/reverse field sweeps and frequency-grid refinement.
6. Weak electro-nuclear satellites remain visible and mode identities do not jump merely because their intensities cross.
7. At least one regression field sits inside the window between `Bc_1z` and the bare-MF boundary, asserting the auto/bare-MF leg ordered while `phase_1z` is PM, AND verifying the strict-PM state was actually used for the 1/z spectrum there -- the split must be exercised by a test, not merely typed. (Stage-1 plan Task 1b.)

## Practical recommendation

Retain the restored projected behavior as the reproducible baseline and explicitly annotate its bare-MF ordered-side limitation. Do not reintroduce the PM-first, boundary-linearized handoff into the projected driver **while it retains a shared-state RPA overlay** -- that combination, not the handoff mechanism, is what broke the spectra (Section 4 above; the tensor working tree uses the same handoff without an RPA overlay and is unaffected). Work order: first repair the failing fast-suite comparison so the regression gate is meaningful (verification record; plan Task 0); then ship the Stage-1 1/z-leg mitigation with honest `Bc_auto` naming and validity masks (including the window regression); then the small RPA-independent dispatcher follow-up that discharges regression 4; only then add and validate a complete ordered 1/z state -- J 2.28--2.29 plus J 2.31--2.34 on the explicitly chosen projected or tensor vehicle.
