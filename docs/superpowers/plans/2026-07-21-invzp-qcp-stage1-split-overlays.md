# Projected QCP Stage-1: Split Auto-Overlay / 1z State Paths — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Revision:** amended 2026-07-21 per `invzp_QCP_plan_review_Codex.md` round 1 (Task 0 baseline repair; crit gates; window regression; J-number citations), round 2 (findings re-verified, window probe reproduced on this machine: auto-dispatch naming `S.Bc_auto`, masking before peaks, probe-qualified anchors with a state-use assertion), and round 3 (the auto overlay is masked at phase 0 as well as suspect — the legacy Σ=0 fallback is no longer surfaced; suspect points are excluded from the `Bc_auto` predicate via the now-public `invz_boundary_field` helper with a direct unit check; residual RPA naming and the README contradiction removed).

**Goal:** Give the 1/z curve in `invz_spectra_map` its own stability-based phase leg (PM state above `Bc_1z`, bare-MF ordered curve kept as a labelled diagnostic below), with validity masks on untrusted columns. The `Sigma = 0` overlay keeps the auto-dispatch state and is *named* as such (`Bc_auto`): a truly RPA-independent dispatcher is a separate follow-up, so diagnosis regression 4 remains open after this plan.

**Architecture:** Least-rework *interim* realisation of Stage 1 of `invzp_QCP_diagnosis.md`. The restored ordered-first `invz_solve_auto` selection flips at the bare-MF boundary, so below it the `Sigma = 0` overlay approximates the RPA state — exactly when the ordered 1/z EMT loop converges up to that boundary. That approximation is kept unchanged and reported honestly as auto-dispatch. What actually changes is the 1/z curve: between `Bc_1z` and the bare boundary the strict-PM 1/z solve converges with `crit > 0` and is the consistent 1/z state. All changes are confined to `invz_spectra_map.m` (one extra PM solve attempt per moment-form field), plus display surfacing in `invz_run_spectra.m` and doc updates. No solver files change behavior; `invz_solve_auto`, `invz_solve_point`, `invz_solve_point_ordered` are reused as-is.

**Accepted interim coupling (re-review finding 1):** `S.phase` is the auto dispatch, which depends on ordered-1/z EMT convergence. Known failure modes: a spuriously converged below-`Bc` PM point (`crit <= 0` — flagged in `S.suspect` and masked); `phase = 0` although the bare-MF/RPA state exists (ordered EMT failure — a visibly masked column); and `phase = 2` with `crit > 0` below the bare boundary (ordered EMT failure inside the window — indistinguishable from a true window point without an RPA-independent check). This is why the outputs are named `S.phase`/`S.Bc_auto` rather than `phase_rpa`/`Bc_rpa`, and why regression 4 stays open. The independent bare-MF/RPA evaluator on `1 − Jχ0(0)` is a small self-contained follow-up task immediately after Stage 1 — not bundled into Stage-2 thermodynamics.

**Tech Stack:** MATLAB R2025a, matlab.unittest function-based tests (`functiontests(localfunctions)`), existing invz_projected solver stack.

## Global Constraints

- Run every MATLAB command from the repository root; the path contains spaces — keep the full binary path quoted exactly: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "..."`.
- Fast-suite gate after each task: `runtests('invz_projected/tests')` must report **0 failed**. The current tree is NOT green — `test_invz_chiperp/test_anchors_and_symmetry` fails on numerical-test fragility (verified 2026-07-21); Task 0 repairs it FIRST. Expected counts: after Task 0 **143 passed / 0 failed / 19 incomplete**; after Task 1 **144 passed**; after Task 1b **145 passed** (or 144 passed / 20 incomplete if the `INVZ_SLOW` fallback is taken).
- `S.phase` keeps its existing semantics (the auto-solve selection) — existing tests and plotters must not break.
- Error policy everywhere: absorb only `invz:*` identifiers; anything else rethrows (repo convention, visible throughout `invz_spectra_map.m`).
- The strict-PM solve must NEVER be attempted at `B(3) ~= 0` (its m = 0 gate raises `invz:orderedPhase` and would abort the parfor — review finding 5, already encoded in the current code).
- New `one_field` outputs must be parfor-sliced plain arrays (repo convention, header comment at `invz_spectra_map.m:94`).
- **Staging caution:** the working tree carries unrelated pending edits (`invz_projected/invz_run_spectra.m` driver knobs; `invz_projected/README.html` test-quickstart/caveats hunks; several `invz_tensor/` files). Stage ONLY the hunks this plan creates — never `git add -A`. Tasks 2 and 3 touch those same files; see the per-task staging notes.
- `fields` is assumed ascending for the boundary midpoints (both drivers use `linspace`); `invz_boundary_field` returns NaN when the sweep does not bracket a boundary, and its callers own the validity predicates (suspect points must never anchor `Bc_auto`).

---

### Task 0: Repair the fast-suite baseline (`test_invz_chiperp`)

**Files:**
- Modify: `invz_projected/tests/test_invz_chiperp.m:15`

**Interfaces:** none — test-only tolerance repair. The failure (verified by running the file): `test_anchors_and_symmetry` compares the computed 2×2 `Xp` against the pinned anchor `A.chiperp_1p53K_0T` with `'RelTol', 1e-9` only; the anchor's off-diagonals are nominally zero (~−1e-15 against diagonals ~17.64), and machine-precision noise on a nominal zero can never satisfy a relative tolerance. This is numerical-test fragility, not a QCP regression.

- [ ] **Step 1: Confirm the failure**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests/test_invz_chiperp.m'); disp(table(results)); assertSuccess(results)"
```
Expected: `test_anchors_and_symmetry` FAILS by verification (off-diagonal elements at ~1e-15), the other three pass.

- [ ] **Step 2: Add an absolute-tolerance floor**

Old (line 15):
```matlab
verifyEqual(testCase, Xp, A.chiperp_1p53K_0T, 'RelTol', 1e-9);
```
New:
```matlab
% AbsTol floor for the nominally-zero off-diagonals (machine noise on a nominal zero can
% never satisfy RelTol); the ~17.6 diagonals stay governed by RelTol (1e-9*17.6 >> 1e-12).
verifyEqual(testCase, Xp, A.chiperp_1p53K_0T, 'RelTol', 1e-9, 'AbsTol', 1e-12);
```
(MATLAB `verifyEqual` passes per element when EITHER tolerance is satisfied, so this weakens nothing on the diagonals.)

- [ ] **Step 3: Verify the file and the full fast suite pass**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests'); disp(results); assertSuccess(results)"
```
Expected: **143 passed / 0 failed / 19 incomplete** — this re-establishes the green baseline all later gates assume.

- [ ] **Step 4: Commit**

```bash
git add invz_projected/tests/test_invz_chiperp.m
git commit -m "test(invz): AbsTol floor for nominally-zero chiperp off-diagonals -- repairs the fast-suite baseline (RelTol can never pass on machine noise at a nominal zero)"
```

---

### Task 1: Two-leg phase dispatch in `invz_spectra_map`

**Files:**
- Modify: `invz_projected/invz_spectra_map.m` (docstring lines 9–19, sliced-array init lines 95–96, parfor body lines 106–113, pack lines 118–121, full `one_field` replacement lines 124–191)
- Create: `invz_projected/invz_boundary_field.m` (public boundary-midpoint helper — public so the suspect-skip predicate is directly unit-testable, re-review R3-2)
- Test: `invz_projected/tests/test_invz_spectra_map.m` (append one local function)

**Interfaces:**
- Consumes: `invz_solve_auto`, `invz_solve_point`, `invz_chi_realaxis`, `invz_zero_sigma_overlay`, `invz_twolevel`, `getf` — all unchanged, existing signatures.
- Produces (relied on by Tasks 1b/2): `S.phase_1z [1 x nB]` (0 = no valid 1/z label — masked column OR spurious `crit <= 0` PM point / 1 = moment-form diagnostic / 2 = stable PM, always `crit > 0` by construction), `S.crit_pm [1 x nB]` (PM 1/z mass, NaN where no PM solve returned), `S.suspect [1 x nB]` logical (auto `phase = 2` with finite `crit <= 0` — spurious below-`Bc` PM points, masked out of BOTH overlays), `S.Bc_auto`, `S.Bc_1z` (scalar midpoint estimates, NaN if unbracketed; `Bc_auto` is the auto-dispatch boundary — an RPA proxy only where the ordered EMT converged to the bare boundary, hence the name). Masking contract (re-review R3-1): `S.chiz` columns with `phase_1z = 0` are NaN; `S.chirpa` columns with NO accepted auto state (`phase = 0` — this retires the legacy Σ=0 fallback display, which one_field still computes but the pack no longer surfaces) OR `suspect` are NaN. Both applied BEFORE `Epeak`/`Epeak_rpa` extraction, so untrusted columns never reach plots or peak curves. `Bc_auto` uses only confirmed-stable PM anchors (`phase = 2` with finite `crit > 0`). `one_field` signature becomes `[chiz, chirpa, Sigma0, phase, phase_1z, crit_pm] = one_field(...)` (same inputs). `S.Sigma0` semantics change slightly: it is now always the Sigma0 of the state used for `S.chiz` (PM leg's value where `phase_1z = 2` under a moment-form `phase = 1`).

- [ ] **Step 1: Write the failing test**

Append to `invz_projected/tests/test_invz_spectra_map.m` (after `test_map_shape_and_phase_codes`, same synthetic-coupling pattern):

```matlab
function test_split_phase_diagnostics(testCase)
% Stage-1 QCP split (docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md):
% the auto overlay keeps the auto-solve (bare-MF) selection, the 1/z curve carries its own
% stability-based phase label, and the per-field PM 1/z mass plus midpoint boundary
% estimates are returned. Same cheap synthetic coupling as the structural test above.
ion = invz_ion();
T = 0.31;
fields = [2.0 5.5];                       % 2 T: ordered candidate;  5.5 T: paramagnet
w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
S = invz_spectra_map(ion, T, fields, w, ...
                     struct('Jnu', Jnu, 'info', info, 'verbose', false));

verifySize(testCase, S.phase_1z, [1 numel(fields)]);
verifySize(testCase, S.crit_pm,  [1 numel(fields)]);
verifyTrue(testCase, all(ismember(S.phase_1z, [0 1 2])));
% the phase_1z = 2 label is GATED on a stable PM mass -- unconditional by construction
lab2 = S.phase_1z == 2;
verifyTrue(testCase, all(isfinite(S.crit_pm(lab2))));
if any(lab2), verifyGreaterThan(testCase, S.crit_pm(lab2), 0); end
% a PM auto phase with stable mass carries the PM 1/z label; spurious PM (crit <= 0)
% must land in S.suspect with phase_1z = 0, never in phase_1z = 2
verifyTrue(testCase, isfield(S, 'suspect') && islogical(S.suspect));
verifyEqual(testCase, S.suspect, (S.phase == 2) & isfinite(S.crit_pm) & (S.crit_pm <= 0));
verifyTrue(testCase, all(S.phase_1z(S.phase == 2 & ~S.suspect & isfinite(S.crit_pm)) == 2));
% validity masks reach the returned spectra (re-review findings 2, R3-1): no valid 1/z
% label -> chiz column all-NaN; NO accepted auto state (phase 0) OR suspect -> chirpa
% column all-NaN (the legacy phase-0 Sigma-zero fallback is no longer surfaced)
verifyTrue(testCase, all(all(isnan(S.chiz(:,  S.phase_1z == 0)))));
verifyTrue(testCase, all(all(isnan(S.chirpa(:, S.phase == 0 | S.suspect)))));
% boundary estimates are always present as fields (NaN allowed: two fields cannot bracket
% a boundary that is not crossed between them)
verifyTrue(testCase, isfield(S, 'Bc_auto') && isfield(S, 'Bc_1z'));
% Bc_auto predicate (re-review R3-2): a suspect PM point between the last ordered and the
% first VALID PM point is excluded and WIDENS the bracket -- direct unit check on the
% public helper with crafted labels (index 3 plays the suspect point)
fieldsU = [1 2 3 4];
verifyEqual(testCase, invz_boundary_field(fieldsU, logical([1 1 0 0]), logical([0 0 0 1])), 3.0);
verifyEqual(testCase, invz_boundary_field(fieldsU, logical([1 1 0 0]), logical([0 0 1 1])), 2.5);
verifyTrue(testCase, isnan(invz_boundary_field(fieldsU, logical([1 1 0 0]), false(1, 4))));
% the 1/z map itself must still be a well-formed spectrum
verifyTrue(testCase, any(isfinite(S.chiz(:))));
end
```

- [ ] **Step 2: Run the test file to verify the new test fails**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests/test_invz_spectra_map.m'); disp(results); assertSuccess(results)"
```
Expected: `test_map_shape_and_phase_codes` PASSES, `test_split_phase_diagnostics` FAILS (unrecognized field `phase_1z` in struct `S`).

- [ ] **Step 3: Update the `invz_spectra_map.m` docstring**

Edit (old → new), keeping surrounding lines intact:

Old (line 10):
```matlab
%     S.chiz   [nw x nB]  1/z-renormalized chi''_cc  (moment-form below Bc, PM above)
```
New:
```matlab
%     S.chiz   [nw x nB]  1/z-renormalized chi''_cc at the 1/z theory's OWN phase: the
%                         strict-PM state wherever it is stable (crit > 0, i.e. above
%                         Bc_1z), the bare-MF moment-form state below Bc_1z (DIAGNOSTIC
%                         there -- the complete ordered 1/z state is Stage 2 of
%                         invzp_QCP_diagnosis.md and does not exist yet)
```

Old (lines 14–19):
```matlab
%   Per-field diagnostics: S.phase (1 = moment-form (spontaneous FM below Bc, or field-induced
%   under a longitudinal tilt -- a rounded crossover, no sharp Bc), 2 = strict paramagnet,
%   0 = masked), S.Sigma0, and S.Epeak/S.Epeak_rpa (censored, parabolic-refined peak energy;
```
New:
```matlab
%   Per-field diagnostics: S.phase (1 = moment-form (spontaneous FM below Bc, or field-induced
%   under a longitudinal tilt -- a rounded crossover, no sharp Bc), 2 = strict paramagnet,
%   0 = masked). RPA and 1/z are DIFFERENT theories with different critical fields
%   (invzp_QCP_diagnosis.md). S.phase stays the AUTO dispatch (ordered-first bare-MF); the
%   Sigma = 0 overlay built from it approximates the RPA state only where the ordered 1/z
%   EMT converged to the bare boundary -- an RPA-independent dispatcher is a scheduled
%   follow-up (diagnosis regression 4 stays OPEN), hence the boundary output is named
%   Bc_auto, never Bc_rpa. S.phase_1z is the 1/z leg's own label (2 = stable PM 1/z state
%   used for S.chiz, ALWAYS gated on crit > 0; 1 = bare-MF ordered diagnostic below
%   Bc_1z; 0 = no valid 1/z label: masked column OR a spurious converged PM point with
%   crit <= 0; mirrors S.phase under a longitudinal tilt -- no strict-PM phase there).
%   S.suspect flags the spurious PM columns (S.phase = 2 with finite crit <= 0). Masking
%   before peak extraction: chiz is NaN where phase_1z = 0; chirpa is NaN where there is
%   NO accepted auto state (S.phase = 0 -- the legacy Sigma-zero fallback is computed but
%   no longer surfaced) or S.suspect. S.crit_pm is the per-field PM 1/z mass
%   1 + Sigma0 - J0eff*chi0cc0 (NaN
%   where no PM solve returned); S.Bc_auto (anchored only on valid, non-suspect PM
%   points) / S.Bc_1z are sweep-midpoint boundary
%   estimates (NaN when the sweep does not bracket the flip; precision = half the gap
%   between the bracketing labelled fields, so masked or suspect columns between them
%   widen it -- refine with invz_critical for a solver-grade Bc). S.Sigma0 is the Sigma0
%   of the state used for
%   S.chiz. S.Epeak/S.Epeak_rpa (censored, parabolic-refined peak energy;
```

- [ ] **Step 4: Add the sliced arrays and extend the parfor body**

Old (lines 95–96):
```matlab
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);
```
New:
```matlab
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);
ph1z    = zeros(1, nB);  critPM  = nan(1, nB);
```

Old (lines 107–112):
```matlab
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k)] = ...
        one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol);
    if verbose
        ph = {'masked (no converged solve)', 'moment-form (FM or field-induced)', 'paramagnet'};
        fprintf('  |B| = %5.2f T : %-34s Sigma0 = %s\n', fields(k), ph{phaseC(k)+1}, num2str(Sig0(k)));
    end
```
New:
```matlab
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k), ph1z(k), critPM(k)] = ...
        one_field(ion, T, BvecM(k, :), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol);
    if verbose
        ph = {'masked', 'moment-form', 'paramagnet'};
        fprintf('  |B| = %5.2f T : auto-state %-11s 1z-state %-11s Sigma0 = %s\n', ...
                fields(k), ph{phaseC(k)+1}, ph{ph1z(k)+1}, num2str(Sig0(k)));
    end
```

Old (lines 118–121):
```matlab
S.chiz = chizM;  S.chirpa = chirpaM;
S.Sigma0 = Sig0;  S.phase = phaseC;
S.Epeak     = invz_peak_energy(chizM,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpaM, w, wmin);
```
New:
```matlab
S.Sigma0 = Sig0;  S.phase = phaseC;   % auto (ordered-first bare-MF) dispatch -- see header
S.phase_1z = ph1z;  S.crit_pm = critPM;
S.suspect = (phaseC == 2) & isfinite(critPM) & (critPM <= 0);   % spurious below-Bc auto-PM
% Validity masks BEFORE the spectra are packed and the peaks extracted (re-review
% findings 2, R3-1; R4-3 hardening): chiz is masked where there is no valid 1/z label;
% the auto overlay is masked wherever there is NO accepted auto state (phase 0 --
% including the legacy Sigma-zero fallbacks, still computed inside one_field but no
% longer surfaced) or the PM mass is not confirmed positive (suspect OR non-finite).
% The plotter then hides them (NaN -> transparent) and invz_peak_energy never sees them.
valid_auto_pm = (phaseC == 2) & isfinite(critPM) & (critPM > 0);
invalid_auto  = (phaseC == 0) | ((phaseC == 2) & ~valid_auto_pm);
chizM(:, ph1z == 0)      = NaN;
chirpaM(:, invalid_auto) = NaN;
S.chiz = chizM;  S.chirpa = chirpaM;
S.Bc_auto = invz_boundary_field(fields, phaseC == 1, valid_auto_pm);   % suspect never anchors
S.Bc_1z   = invz_boundary_field(fields, ph1z  == 1, ph1z  == 2);
S.Epeak     = invz_peak_energy(chizM,   w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpaM, w, wmin);
```

- [ ] **Step 5: Replace `one_field` wholesale and create `invz_boundary_field.m`**

Replace the entire `one_field` local function (current lines 125–191) with:

```matlab
function [chiz, chirpa, Sigma0, phase, phase_1z, crit_pm] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts, bztol)
%ONE_FIELD chi''_cc(omega) at one field -- TWO independently phased legs (QCP Stage 1,
% docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md):
%   auto/overlay leg (chirpa): the ordered-first invz_solve_auto selection (bare-MF
%     moment), UNCHANGED -- phase = 1 moment-form (FM or field-induced), 2 strict PM,
%     0 masked, historical semantics. The Sigma = 0 overlay built from it approximates
%     the RPA state only where the ordered EMT converged to the bare boundary; a
%     1/z-independent RPA dispatcher is a scheduled follow-up (diagnosis regression 4
%     stays open), hence the driver reports Bc_auto, never Bc_rpa.
%   1/z leg (chiz): the 1/z theory's own stability. Wherever the strict-PM solve converges
%     with crit = 1 + Sigma0 - J0eff*chi0cc0 > 0, the PM state is the consistent 1/z state
%     (phase_1z = 2) even though the bare-MF leg still orders; below Bc_1z the bare-MF
%     ordered 1/z curve is kept as a DIAGNOSTIC (phase_1z = 1) until the full ordered 1/z
%     state exists (diagnosis Stage 2). EVERY phase_1z = 2 label is gated on crit > 0;
%     phase_1z = 0 means no valid 1/z label (masked column OR a spurious converged PM
%     point with crit <= 0). Longitudinal tilts have no strict-PM phase
%     (rounded crossover): phase_1z mirrors phase there.
% crit_pm: PM 1/z mass at this field (NaN where no PM solve was attempted or returned).
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
phase_1z = 0;  crit_pm = NaN;
tmf = getf(sopts, 'transverse_mf', 'legacy_x');
% transverse_mf here so any fallback single-ion rebuild inside invz_chi_realaxis (when a
% branch below can't supply si) uses the SAME MF model as the solve, not the 'legacy_x'
% default (review finding I1; the C1 companion bug in invz_spectra_qpath.m).
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp, ...
               'transverse_mf', tmf);

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- moment-form branch (FM or induced) ---
    % 1/z leg: strict-PM stability probe, TRANSVERSE ONLY (the strict-PM solver's m = 0
    % gate would raise invz:orderedPhase under a longitudinal tilt -- review finding 5).
    ptp = [];
    if B(3) == 0
        try
            ptp = invz_solve_point(ion, T, B, Jnu, sopts);
        catch err
            if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        end
    end
    if ~isempty(ptp) && isfield(ptp, 'crit') && isfinite(ptp.crit), crit_pm = ptp.crit; end
    % crit_pm doubles as the field-guarded crit: pm_valid never touches ptp.crit directly
    % (an early-return ptp may lack the field).
    pm_valid = ~isempty(ptp) && ptp.converged && isfinite(ptp.Sigma0) && ...
               isfinite(crit_pm) && crit_pm > 0;

    if pm_valid                         % --- Bc_1z < |B| < bare boundary: 1/z theory is PM ---
        phase_1z = 2;  Sigma0 = ptp.Sigma0;
        copts1 = copts;  copts1.si = ptp.si;      % PM chi0 differs from the ordered leg's:
        o1 = invz_chi_realaxis(ion, T, B, ptp, w, copts1);   % no chi0cc_w sharing possible
        chiz = imag(o1.chi_cc_q(1, :)).';
        c0opts = copts;  c0opts.npass = 1;        % RPA overlay stays on the bare-MF state
        o0 = invz_chi_realaxis(ion, T, B, invz_zero_sigma_overlay(pt), w, c0opts);
        chirpa = imag(o0.chi_cc_q(1, :)).';
    else                                     % --- below Bc_1z: ordered 1/z, DIAGNOSTIC ---
        phase_1z = 1;  Sigma0 = pt.Sigma0;
        o  = invz_chi_realaxis(ion, T, B, pt, w, copts);   % reuses pt.si (moment-form eigenstates)
        chiz = imag(o.chi_cc_q(1, :)).';
        pt0 = invz_zero_sigma_overlay(pt);
        c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;   % share bare cc
        o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
        chirpa = imag(o0.chi_cc_q(1, :)).';
    end
    return;
end

if abs(B(3)) > bztol
    % --- longitudinal failure: NEVER the strict-paramagnet overlay (its m = 0 gate
    % would raise invz:orderedPhase and abort the parfor -- review finding 5). If the
    % failed moment-branch pt still carries valid si/tl, compute the RPA-only overlay
    % from the ordered-style pt0; otherwise leave the whole column masked.
    if ~isempty(pt) && ~isempty(pt.si) && isfield(pt, 'tl') && ~isempty(pt.tl)
        pt0 = invz_zero_sigma_overlay(pt);
        c0opts = copts;  c0opts.npass = 1;
        try
            o0 = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        catch err
            if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        end
    end
    if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
    return;
end

% --- transverse paramagnetic side: unchanged historical logic --------------------------
if phase == 2 && ~isempty(pt), tl0 = pt.tl;  si0 = pt.si;
else, tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0, 'transverse_mf', tmf));  si0 = []; end
chi0cc = [];
try
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    c0opts = copts;  c0opts.npass = 1;  c0opts.si = si0;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    chi0cc = o0.chi0cc_w;                          % share the bare cc with the 1/z call
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
end
if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
if ~isempty(pt) && isfield(pt, 'crit') && isfinite(pt.crit), crit_pm = pt.crit; end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    if isfinite(crit_pm) && crit_pm > 0
        phase_1z = 2;                             % the auto PM state IS the 1/z leg here
    end                                           % crit <= 0: spurious below-Bc PM point
                                                  % (documented tensor-side failure mode) --
                                                  % phase_1z stays 0, column flagged suspect
    copts1 = copts;  copts1.si = pt.si;           % reuse the solve's si (C1 companion: pt is
    if ~isempty(chi0cc), copts1.chi0cc_w = chi0cc; end   % not is_ordered, so without si this
    o = invz_chi_realaxis(ion, T, B, pt, w, copts1);     % would silently rebuild at 'legacy_x'
    chiz = imag(o.chi_cc_q(1, :)).';
end
end
```

Then create `invz_projected/invz_boundary_field.m` (public — the caller owns the validity predicates, and a public function makes the suspect-skip behavior directly unit-testable, re-review R3-2):

```matlab
function Bc = invz_boundary_field(fields, isBelow, isAbove)
%INVZ_BOUNDARY_FIELD midpoint estimate of a phase boundary from per-field labels.
% Bc between the last field labelled "below" and the first field labelled "above";
% fields assumed ascending (the drivers use linspace). NaN when the sweep does not
% bracket the flip. The CALLER owns the predicates: pass only VALID labels (e.g. exclude
% suspect auto-PM points, which then WIDEN the bracket instead of anchoring it -- QCP
% Stage 1). Precision is half the gap between the two bracketing labelled fields:
% masked/suspect/other-labelled columns between them widen it beyond half the field
% step. Overlapping/nonmonotone label sequences (reentrant or noisy near-boundary
% labels) return NaN whenever the first above-label precedes the last below-label --
% refine the sweep or inspect manually; use invz_critical for a solver-grade Bc.
Bc = NaN;
kb = find(isBelow, 1, 'last');  ka = find(isAbove, 1, 'first');
if ~isempty(kb) && ~isempty(ka) && ka > kb
    Bc = 0.5*(fields(kb) + fields(ka));
end
end
```

The diagnostic (else) branch reproduces the current moment-form code path verbatim; the PM-side section keeps the historical call structure but its label semantics change (the `phase_1z` gate). Do not restructure further — correctness is established by the behavioral tests, not by textual identity with the old code.

- [ ] **Step 6: Run the test file to verify both tests pass**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests/test_invz_spectra_map.m'); disp(results); assertSuccess(results)"
```
Expected: 2 passed, 0 failed. (The 2 T column may now attempt a non-convergent PM solve — bounded by the EMT iteration cap; the test stays fast.)

- [ ] **Step 7: Run the fast-suite regression gate**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests'); disp(results); assertSuccess(results)"
```
Expected: 144 passed / 0 failed / 19 incomplete (`INVZ_SLOW`-gated). Any failure in `test_invz_solve_auto`, `test_invz_ordered_phase`, `test_invz_field_angle_spectra`, or `test_invz_phi_spectra` means the shared solver call paths regressed behaviorally — fix before proceeding.

- [ ] **Step 8: Commit (only the three files this task touched)**

```bash
git add invz_projected/invz_spectra_map.m invz_projected/invz_boundary_field.m invz_projected/tests/test_invz_spectra_map.m
git commit -m "feat(invz): stage-1 QCP split -- 1/z leg gets its own stability-gated phase in invz_spectra_map; auto overlay honestly named (Bc_auto, valid-PM anchors only); phase-0/suspect columns masked from both overlays before peaks; S.phase_1z/S.crit_pm/S.suspect diagnostics (invzp_QCP_diagnosis.md)"
```

---

### Task 1b: Exercise the split window in the regression (re-review findings 3–4)

**Files:**
- Modify: `invz_projected/tests/test_invz_spectra_map.m` (add one local function)
- Scratch: window-probe script in the session scratchpad (throwaway; do not commit)

**Interfaces:**
- Consumes: Task 1's `S.phase`, `S.phase_1z`, `S.crit_pm`, `S.Bc_auto`, `S.Bc_1z`, plus `invz_solve_point`/`invz_chi_realaxis` directly for the state-use check.
- Produces: a test that FAILS if the two-leg behavior is absent — the structural test of Task 1 only proves the fields exist; this one proves a field where the legs disagree (`phase = 1`, `phase_1z = 2`) is produced AND that the strict-PM state was actually used for the 1/z spectrum there (a wiring bug that sets the label but keeps the ordered spectrum must fail).

- [ ] **Step 1: Probe and qualify THREE anchor fields**

Save as `probe_qcp.m` in the scratchpad (a multi-line `-batch` string does not survive the shell — use a script) and run with `-batch "run('<scratchpad>/probe_qcp.m')"`:

```matlab
root = '<absolute repo root>';
addpath(fullfile(root, 'invz_projected'));  addpath(fullfile(root, 'invz_common'));  addpath(root);
ion = invz_ion(); T = 0.31;
info = struct('Jcc0', 6.4e-3); Jnu = linspace(-2e-3, 6.0e-3, 24).';
sopts = struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'bz_tol', 1e-9);
for B = 2.5:0.05:5.5
    [pt, ph] = invz_solve_auto(ion, T, [B 0 0], Jnu, sopts);
    c = NaN; cv = false;
    try
        ptp = invz_solve_point(ion, T, [B 0 0], Jnu, sopts); c = ptp.crit; cv = ptp.converged;
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    end
    fprintf('B=%.2f  auto_phase=%d  pm_conv=%d  pm_crit=%.4g\n', B, ph, cv, c);
end
```

Select three independently qualified anchors WITH SIGN MARGINS (|crit| ≥ 0.03 — a crit near zero is the LEAST robust against solver/parameter drift, not the best-conditioned):

1. `Blow`: `auto_phase = 1`, PM converged, `pm_crit <= -0.03` → yields `phase_1z = 1`;
2. `Bwin`: `auto_phase = 1`, `pm_crit >= +0.03` → the disagreement point (`phase = 1`, `phase_1z = 2`);
3. `Bhigh`: `auto_phase = 2`, `pm_crit >= +0.03`.

Measured on this machine, 2026-07-21 (matches the re-review's independent run): `B = 2.00` → auto 0 (masked — NOT an anchor; a deep-ordered field does not guarantee `phase_1z = 1`); `2.85` → auto 1, crit −0.0579; `3.20` → auto 1, crit +0.0550; `3.30` → auto 1, crit +0.0867; `4.50`/`5.50` → auto 2, crit +0.354/+0.471. Anchors: **[2.85, 3.30, 5.50]** (3.30 preferred over 3.20 for the larger margin). Re-run the probe and update the hard-coded values if the synthetic model or solver defaults change.

- [ ] **Step 2: Add the window regression with the state-use check**

Append to `test_invz_spectra_map.m` (values from Step 1):

```matlab
function test_split_window_exercised(testCase)
% Re-review findings 3-4: the split must be EXERCISED at probe-qualified anchors
% (Blow 2.85 T: auto ordered, PM crit -0.058 -> 1/z diagnostic; Bwin 3.30 T: auto ordered,
% PM crit +0.087 -> 1/z PM = the disagreement; Bhigh 5.50 T: auto PM, crit +0.47), and the
% window column must BE the strict-PM spectrum, not merely carry the label. Probe log:
% docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md Task 1b.
ion = invz_ion();
T = 0.31;
fields = [2.85 3.30 5.50];               % Blow | Bwin | Bhigh
w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
S = invz_spectra_map(ion, T, fields, w, ...
                     struct('Jnu', Jnu, 'info', info, 'verbose', false));

verifyEqual(testCase, S.phase,    [1 1 2]);       % auto: ordered | ordered | PM
verifyEqual(testCase, S.phase_1z, [1 2 2]);       % 1/z: diagnostic | PM (window!) | PM
verifyLessThan(testCase,    S.crit_pm(1), -0.03); % sign margins, not near-zero values
verifyGreaterThan(testCase, S.crit_pm(2),  0.03);
verifyGreaterThan(testCase, S.crit_pm(3),  0.03);

% State-use check (re-review finding 4): reconstruct the strict-PM spectrum at Bwin
% through the same serial code path and require S.chiz to BE it.
sopts = struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'bz_tol', 1e-9);
ptp = invz_solve_point(ion, T, [3.30 0 0], Jnu, sopts);
verifyEqual(testCase, S.Sigma0(2), ptp.Sigma0, 'RelTol', 1e-10);   % PM leg's Sigma0 returned
copts = struct('Jsel', info.Jcc0, 'eta', 5e-3, 'Jxx0', ion.Jxx0, 'Jshape', 0, 'hyp', true, ...
               'transverse_mf', 'legacy_x', 'si', ptp.si);
o = invz_chi_realaxis(ion, T, [3.30 0 0], ptp, w, copts);
verifyEqual(testCase, S.chiz(:, 2), imag(o.chi_cc_q(1, :)).', 'RelTol', 1e-8, 'AbsTol', 1e-12);

% both boundaries bracketed by this 3-point sweep, renormalized below the auto/bare one:
% Bc_auto = (3.30+5.50)/2 = 4.4;  Bc_1z = (2.85+3.30)/2 = 3.075
verifyTrue(testCase, isfinite(S.Bc_auto) && isfinite(S.Bc_1z));
verifyLessThan(testCase, S.Bc_1z, S.Bc_auto);
end
```

(The `copts` values mirror `one_field`'s: `Jsel = info.Jcc0`, map-default `eta = 5e-3`, `Jxx0 = ion.Jxx0` and `Jshape = 0` because the synthetic `info` carries no `Jaa0`/`Jshape_cc`, default `transverse_mf`. Identical deterministic serial calls justify the tight tolerance.)

Fallback if a future model change closes the synthetic window (all `auto_phase = 1` probes have `crit < 0.03`): port the same three-anchor test to production couplings (`invz_bz_couplings`, `[16 16 16]` grid), gate it with the repo's `assumeTrue(... getenv('INVZ_SLOW') ...)` slow-suite convention, discover the anchors with the same probe on a production field sweep, and record them in the execution record. Boundary assertion becomes `S.Bc_1z < S.Bc_auto` over the discovered window, unchanged in form.

- [ ] **Step 3: Run the test file, then the fast-suite gate**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests/test_invz_spectra_map.m'); disp(results); assertSuccess(results)"
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests'); disp(results); assertSuccess(results)"
```
Expected: **145 passed / 0 failed / 19 incomplete** (fallback route: **144 passed / 0 failed / 20 incomplete**; run once with `INVZ_SLOW=1` to confirm the gated test actually passes before committing).

- [ ] **Step 4: Commit**

```bash
git add invz_projected/tests/test_invz_spectra_map.m
git commit -m "test(invz): stage-1 window regression -- probe-qualified anchors with sign margins, auto/1z disagreement exercised, strict-PM state-use asserted (re-review findings 3-4)"
```

---

### Task 2: Surface `Bc_auto` / `Bc_1z` in the field-sweep driver

**Files:**
- Modify: `invz_projected/invz_run_spectra.m:156-184` (field-sweep branch only)

**Interfaces:**
- Consumes: `S.Bc_auto`, `S.Bc_1z`, `S.suspect`, `S.fields` from Task 1 (scalar doubles, NaN allowed; `[1 x nB]` logical/double).
- Produces: display only — no struct changes.

**Staging caution:** this file already carries unrelated uncommitted knob edits (T, fields, w, eta, phi_ab, transverse_mf). Before committing, show the user `git diff invz_projected/invz_run_spectra.m` and either (a) they commit/stash their knob edits first, or (b) commit this task's hunks knowingly alongside them with their approval. Do not silently bundle.

- [ ] **Step 1: Add the boundary printout after the map solve**

After line 157 (the `S = invz_spectra_map(...)` call), insert:

```matlab
    fprintf('Bc_auto (bare-MF dispatch; RPA proxy only where the ordered EMT converged) ~ %s T | Bc_1z (renormalized) ~ %s T  (sweep midpoints; masked/suspect columns widen the bracket)\n', ...
            num2str(S.Bc_auto), num2str(S.Bc_1z));
    if any(S.suspect)
        fprintf('WARNING: %d suspect column(s) (auto-PM with crit <= 0 -- spurious below-Bc PM points; masked): |B| = %s T\n', ...
                nnz(S.suspect), mat2str(S.fields(S.suspect), 3));
    end
```

- [ ] **Step 2: Label the 1/z map panel**

Old (line 169):
```matlab
        ax1 = subplot(1, 2, 1);  invz_plot_spectra_map(ax1, Splot, Splot.chiz,   sprintf('1/z, T = %.2f K%s', T, tiltStr), eUnit);
```
New:
```matlab
        ax1 = subplot(1, 2, 1);  invz_plot_spectra_map(ax1, Splot, Splot.chiz,   sprintf('1/z (own phase; FM side diagnostic below B_c^{1/z}), T = %.2f K%s', T, tiltStr), eUnit);
```

- [ ] **Step 3: Mark both boundaries on the E_peak figure**

In the `showPeaks` block, after `legend show;` (line 181), insert:

```matlab
        if isfinite(S.Bc_auto), xline(S.Bc_auto, '--', 'B_c^{auto}', 'HandleVisibility', 'off'); end
        if isfinite(S.Bc_1z),   xline(S.Bc_1z,   ':',  'B_c^{1/z}', 'HandleVisibility', 'off'); end
```

- [ ] **Step 4: Verify the driver parses and the suite stays green**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "t = mtree('invz_projected/invz_run_spectra.m', '-file'); assert(isempty(mtfind(t, 'Kind', 'ERR')), 'parse error in invz_run_spectra.m'); results = runtests('invz_projected/tests'); disp(results); assertSuccess(results)"
```
Expected: no assertion (mtree finds no parse error — lint *warnings* are fine and not checked), then the post-Task-1b count: 145 passed / 0 failed / 19 incomplete (or 144 / 0 / 20 on the `INVZ_SLOW` fallback route). (No test executes the driver script; the mtree ERR check is the parse gate.)

- [ ] **Step 5: Commit (after resolving the staging caution above)**

```bash
git add invz_projected/invz_run_spectra.m
git commit -m "feat(invz): surface Bc_auto/Bc_1z in invz_run_spectra -- boundary printout, suspect-column warning, E_peak xlines, 1/z panel label"
```

---

### Task 3: Documentation — README callout, Jensen citation policy

**Files:**
- Modify: `invz_projected/README.html:378` (the "Both phases, automatically" callout)
- Modify: `invz_projected/invz_solve_point_ordered.m:27` (volatile HTML equation numbers)
- Modify: `invz_common/invz_sigma_ordered.m:2` (volatile HTML equation numbers)

**Interfaces:** none (prose only).

**Staging caution:** `README.html` already carries unrelated uncommitted hunks (test-quickstart rewrite, caveats callout). Same rule as Task 2: surface the pending diff to the user and agree on ordering before committing.

- [ ] **Step 1: Rewrite the README phase callout**

Old (line 378):
```html
<div class="callout note"><span class="tag">Both phases, automatically</span> At each field <code>invz_spectra_map</code> first tries the ferromagnetic (ordered) solve; if the point has a spontaneous moment it uses that, otherwise the paramagnet. So a sweep across \(B_c\) produces one continuous soft-mode map (the "V" dipping at the transition), FM below and PM above, labelled per field in <code>S.phase</code> (1 = FM, 2 = PM, 0 = no solution). Only genuinely intractable columns (the degenerate doublet at \(B_x\!\to\!0\) where the two-level treatment breaks down in <em>both</em> phases) come back <code>NaN</code> and masked.</div>
```
New:
```html
<div class="callout note"><span class="tag">Both phases, two theories</span> At each field <code>invz_spectra_map</code> first tries the ferromagnetic (ordered) solve; if the point has a spontaneous moment it uses that, otherwise the paramagnet (<code>S.phase</code>: 1 = FM, 2 = PM, 0 = no solution — the AUTO dispatch). RPA and 1/z are different approximations with <em>different critical fields</em> (<code>invzp_QCP_diagnosis.md</code>). The \(\Sigma=0\) overlay <code>S.chirpa</code> is built from the auto state, which approximates the RPA state only where the ordered 1/z EMT converged to the bare boundary — its boundary is therefore reported as <code>S.Bc_auto</code>, a dispatch diagnostic rather than a certified \(B_c^{\rm RPA}\) (the 1/z-independent RPA dispatcher, diagnosis regression 4, is a scheduled follow-up). <code>S.chiz</code> follows the 1/z theory's own stability (<code>S.phase_1z</code>: the strict-PM 1/z state wherever its mass <code>S.crit_pm</code> \(=1+\Sigma(0)-\Jc(0)\chi_0^{cc}(0)>0\); the bare-MF ordered curve below \(B_c^{1/z}\) as a flagged <em>diagnostic</em> — its ordered-side pole closure at \(B_c^{1/z}\) is a Stage-2 requirement, not a Stage-1 claim). The two boundaries are reported separately and are <em>not</em> forced to agree. Columns with no valid 1/z label are masked (<code>NaN</code>) in <code>S.chiz</code>; columns with no accepted auto state (phase&nbsp;0) or a spuriously converged below-\(B_c\) PM state (<code>S.suspect</code>, <code>crit</code> &le; 0) mask <em>both</em> overlays and never reach the peak curves. Genuinely intractable columns, including the degenerate doublet at \(B_x\!\to\!0\) where the two-level treatment breaks down in <em>both</em> phases, are likewise returned as <code>NaN</code> and masked.</div>
```

- [ ] **Step 2: Switch stale HTML equation numbers to Jensen J-numbers**

The framework HTML's equation numbers are volatile — the working tree renumbers Section 9 by +2 relative to `HEAD`, which is exactly how two citation errors already happened (see "Citation history" in `invzp_QCP_diagnosis.md` Stage 2). Cite Jensen's own numbers and section anchors instead; do NOT edit `jensen_1z_framework.html` itself (it carries the user's uncommitted revision).

First confirm the exact current texts:
```bash
sed -n '25,29p' invz_projected/invz_solve_point_ordered.m
sed -n '1,4p' invz_common/invz_sigma_ordered.m
```
Then, in order:
1. `invz_projected/invz_solve_point_ordered.m:27`: replace the fragment `the applied-field/H_MF self-consistency of HTML eqs 41-43 is deferred` with `the applied-field/H_MF self-consistency (framework SS9.3, J 2.31-2.33) is deferred` (keep the rest of the line byte-identical).
2. `invz_common/invz_sigma_ordered.m:2`: replace `HTML eqs 37-38 (J 2.26-2.27)` with `framework SS9.2, J 2.26-2.27` (the J-numbers already present are the stable part; drop the HTML numbers).
3. Sweep the REST of the files this task already touches (re-review finding 6 — a policy applied only to two comments is not a policy): `grep -n "HTML eq" invz_projected/invz_solve_point_ordered.m invz_common/invz_sigma_ordered.m` and `grep -nE "eqs? [0-9]+" invz_projected/README.html` — convert every remaining unpaired HTML equation reference IN THESE THREE FILES to a J-number (or add the J-number alongside). Known additional sites: `invz_solve_point_ordered.m:6` and the README scope / deferred-thermodynamics passages.
4. Record, without fixing, the out-of-scope hits: `grep -rn "HTML eq" invz_projected invz_common invz_tensor` — e.g. `invz_chi_realaxis.m` header equations (a file this plan does not otherwise touch) — list them in the execution record as follow-up work.

- [ ] **Step 3: Fast-suite gate (comment-only change, but cheap insurance)**

Run:
```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests'); disp(results); assertSuccess(results)"
```
Expected: the post-Task-1b count — 145 passed / 0 failed / 19 incomplete (or 144 / 0 / 20 on the `INVZ_SLOW` fallback route).

- [ ] **Step 4: Commit (after resolving the README staging caution; `invz_common/invz_sigma_ordered.m` has no pending edits and stages cleanly)**

```bash
git add invz_projected/README.html invz_projected/invz_solve_point_ordered.m invz_common/invz_sigma_ordered.m
git commit -m "docs(invz): two-theory phase semantics in README map callout; cite Jensen J-numbers instead of volatile HTML equation numbers (solve_point_ordered, sigma_ordered headers)"
```

---

## Out of scope (deliberately)

- **Stage 2** (full nonlinear ordered 1/z state — the elastic static sector J 2.28–2.29 plus the H_MF relation J 2.31–2.33 with the free-energy test J 2.34, on the explicitly chosen projected or tensor vehicle): a separate plan after the vehicle decision recorded in `invzp_QCP_diagnosis.md` Stage 2. Note the tensor `a3d` solve alone is NOT that vehicle (no H_MF outer loop; `hmf_J0z` is `a1`-only).
- **Independent bare-MF/RPA state evaluator** on `1 − Jχ0(0)`, decoupled from 1/z convergence: a small self-contained FOLLOW-UP task immediately after Stage 1 (re-review finding 1: not bundled into Stage-2 thermodynamics). Until it lands, diagnosis regression 4 remains open, `Bc_auto` is a dispatch diagnostic rather than `Bc_rpa`, and Stage 1's mitigation is the crit gates plus the `S.suspect` mask.
- **Stage 3** (pole tracking instead of colormap argmax): separate plan.
- **Slow physics tests** for diagnosis regression items 1 and 5 (RPA branch closure at `Bc_rpa` on a fine production intrinsic-response sweep; sweep-direction invariance): minutes-scale solves; add them `INVZ_SLOW`-gated when Stage 2 lands, since item 2 is waived until then. (Item 7, the window regression, IS in scope — Task 1b.)
- **Bc refinement** beyond sweep midpoints (`invz_critical` already exists for solver-grade values — YAGNI here).
- Any change to `invz_solve_auto`, `invz_solve_point`, `invz_solve_point_ordered` behavior, to `invz_tensor/`, or to `jensen_1z_framework.html` (carries the user's uncommitted revision).

## Cost note

Per moment-form transverse field the map now attempts one extra strict-PM solve. Between `Bc_1z` and the bare-MF boundary it converges (that is the point); below `Bc_1z` it burns one bounded non-convergent EMT iteration per field. For the production 301-point sweep expect roughly 1.3–1.8× the ordered-side per-field cost; the PM side above the bare boundary is unchanged.
