# Two-Regime Phase-Boundary Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Trace the LiHoF4 phase boundary in two regimes — the existing Bc(T) field bisection below a user-facing `Tsplit` knob, and a new Tc(B) temperature bisection (`invz_critical_T`) above it, where the boundary is nearly parallel to the field axis and the field bisection is ill-conditioned.

**Architecture:** One new function `invz_critical_T.m` mirrors the validated `invz_critical.m` line for line with the roles of T and B swapped (same `pt.crit` ordered/paramagnet classifier from `invz_solve_point`, same bracket orientation: low edge ordered, high edge paramagnetic). The driver `invz_run_phase_diagram.m` gains `Tsplit`/`Bs`/`Tmax` knobs and runs all boundary points (both regimes) as one flat `parfor` job list. Spec: `docs/superpowers/specs/2026-07-08-phase-diagram-two-regime-design.md`.

**Tech Stack:** MATLAB R2025a (`matlab -batch` CLI), matlab.unittest function-based tests, optional Parallel Computing Toolbox (`parfor` with serial fallback).

## Global Constraints

- The repo path contains spaces — ALWAYS double-quote paths in shell commands.
- MATLAB is NOT on PATH. Invoke it as `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "<code>"` from the repo root.
- Repo root: `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion` — every shell step below assumes `cd` there first.
- Slow tests are gated: `assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), ...)`. The fast suite (`runtests('invz/tests')` without the env var) must stay at **26 passed / 0 failed**, with the filtered count growing from 5 to 7 (the two new tests are SLOW-gated).
- Do NOT modify `invz/invz_critical.m`, `invz/invz_solve_point.m`, or anything below them in the solver layer. They are validated by the existing slow tests.
- Physics constants fixed by the approved spec: `Tsplit = 1.5` K, `Bs = [0.1 0.25 0.5 0.75 1.0 1.25 1.5]` T, `Tmax = 2.0` K, `invz_critical_T` defaults `window = [1.0 2.0]` K and `tol = 0.01` K, crossing-test tolerance 0.05 K, near-zero-field bounds (1.55, 1.76) K.
- Timing expectation with the warm `invz/cache/` on this machine: one EMT solve at T ∈ [1, 2] K takes seconds; a full bisection (~10 solves) takes on the order of a minute; each slow-test command below takes single-digit minutes. A cold cache adds ~15 min once (dipole grid recompute) — do not interpret that as a hang.
- Branch: `invz-1z-lihof4`. Commit after each task. End commit messages with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.

---

### Task 1: `invz_critical_T` mirror bisection + two SLOW-gated tests

**Files:**
- Create: `invz/invz_critical_T.m`
- Test: `invz/tests/test_invz_critical.m` (append two test functions; do not touch the existing ones)

**Interfaces:**
- Consumes: `pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts)` — returns `pt.crit` (scalar; finite and > 0 in the paramagnet, non-finite or ≤ 0 on the ordered side). `[Jf, J0] = lihof4_couplings()` — existing local helper in the test file (16³ q-grid couplings, cached).
- Produces: `tc = invz_critical_T(ion, Bx, Jnu_flat, opts)` — critical temperature (K) at fixed transverse field `Bx` (Tesla). `opts.window = [Tlo Thi]` (default `[1.0 2.0]`), `opts.tol` (default `0.01` K), all other fields passed through to `invz_solve_point` (e.g. `J0eff`). Throws error id `invz:bracket` when `[Tlo, Thi]` does not straddle the boundary. Task 2's driver calls exactly this signature.

- [ ] **Step 1: Write the two failing tests**

Append to the END of `invz/tests/test_invz_critical.m` (after `test_critical_field_at_310mK`):

```matlab
function test_tc_near_zero_field(testCase)
% Small-field limit of the fixed-B temperature bisection: Tc(0.2 T) must sit
% at or just below the closed-form Tc0 = 1.74 K. The undershoot direction is
% R 2007's own small-Bx caveat (the two-level Sigma overestimates the Tc
% suppression when the doublet splitting is below the hyperfine width). The
% window [1.5 2.0] is exactly the driver's high-T default [Tsplit Tmax], so
% this doubles as an integration check of the driver's bracket geometry. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
tc = invz_critical_T(ion, 0.2, Jf, struct('J0eff', J0, 'window', [1.5 2.0]));
verifyGreaterThan(testCase, tc, 1.55);
verifyLessThan(testCase, tc, 1.76);
end

function test_tc_at_fixed_field_crossing(testCase)
% Mirror consistency: at a mid-slope boundary point (T* = 1.4 K, where both
% cut directions are well-conditioned) the fixed-B temperature bisection must
% land back on the fixed-T field bisection's point. 0.05 K tolerance covers
% both bisection tolerances (0.02 T, 0.01 K) at the local boundary slope. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
Tstar = 1.4;
bstar = invz_critical(ion, Tstar, Jf, struct('J0eff', J0, 'window', [0.5 7]));
tc = invz_critical_T(ion, bstar, Jf, struct('J0eff', J0, 'window', [1.0 2.0]));
verifyEqual(testCase, tc, Tstar, 'AbsTol', 0.05);
end
```

Note: `invz_critical` gets an explicit `'window', [0.5 7]` because its default `[2 7]` low edge is NOT ordered at 1.4 K (Bc(1.4 K) ≈ 2 T) — same window the driver uses.

- [ ] **Step 2: Run the cheap test to verify it fails**

`test_tc_near_zero_field` calls `invz_critical_T` immediately, so it goes red fast; `test_tc_at_fixed_field_crossing` would spend minutes in `invz_critical` before failing on the identical missing symbol, so one red run suffices for both.

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "setenv('INVZ_SLOW','1'); results = runtests('invz/tests/test_invz_critical.m', 'ProcedureName', 'test_tc_near_zero_field'); assertSuccess(results)"
```

Expected: FAIL — the error output contains `Unrecognized function or variable 'invz_critical_T'`, and the command exits nonzero.

- [ ] **Step 3: Write the implementation**

Create `invz/invz_critical_T.m`:

```matlab
function tc = invz_critical_T(ion, Bx, Jnu_flat, opts)
%INVZ_CRITICAL_T Critical temperature at fixed transverse field Bx (paramagnetic side).
% Mirror of invz_critical with the roles of T and Bx swapped: bisection on
% pt.crit = 1 + Sigma(0) - J(0)*chi0_cc(0) over temperature at fixed field.
% crit is positive in the paramagnet (high T); inside the ordered phase (low T)
% the paramagnetic EMT fixed point does not exist and invz_solve_point returns
% non-finite crit; non-finite or <=0 is classified as the ordered side.
% Use where the boundary is nearly parallel to the field axis (T near Tc0 =
% 1.74 K): a fixed-B (horizontal) cut crosses it transversally, exactly where
% the fixed-T cut of invz_critical is ill-conditioned.
% opts.window = [Tlo Thi] (K, default [1.0 2.0]): Tlo must be on the ordered
% side (i.e. Bx < Bc(Tlo)) and Thi paramagnetic. opts.tol (K, default 0.01).
% Bx must split the doublet (>~ 0.1 T), or invz_twolevel raises
% invz:degenerateDoublet. Remaining opts pass through to invz_solve_point.
if nargin < 4, opts = struct(); end
win = [1.0 2.0]; if isfield(opts,'window'), win = opts.window; end
tol = 0.01;      if isfield(opts,'tol'),    tol = opts.tol;    end
is_ordered = @(c) ~isfinite(c) || c <= 0;
f = @(T) crit_of(ion, T, Bx, Jnu_flat, opts);
flo = f(win(1));  fhi = f(win(2));
assert(is_ordered(flo) && isfinite(fhi) && fhi > 0, 'invz:bracket', ...
    ['Boundary not bracketed: crit(%.3f K) = %.3g, crit(%.3f K) = %.3g. ' ...
     'Likely Bx = %.2f T exceeds Bc(%.3f K): lower the field or the window low edge.'], ...
    win(1), flo, win(2), fhi, Bx, win(1));
lo = win(1);  hi = win(2);
while hi - lo > tol
    mid = 0.5*(lo + hi);
    if is_ordered(f(mid)), lo = mid; else, hi = mid; end
end
tc = 0.5*(lo + hi);
end

function c = crit_of(ion, T, Bx, Jf, opts)
pt = invz_solve_point(ion, T, Bx, Jf, opts);
c = pt.crit;
end
```

- [ ] **Step 4: Run both new tests to verify they pass**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "setenv('INVZ_SLOW','1'); results = runtests('invz/tests/test_invz_critical.m', 'ProcedureName', 'test_tc_near_zero_field'); assertSuccess(results)"
```

Expected: PASS (~1–3 min warm cache).

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "setenv('INVZ_SLOW','1'); results = runtests('invz/tests/test_invz_critical.m', 'ProcedureName', 'test_tc_at_fixed_field_crossing'); assertSuccess(results)"
```

Expected: PASS (~2–6 min warm cache; it runs two full bisections).

- [ ] **Step 5: Run the fast suite to verify gating**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results); disp(results)"
```

Expected: PASS — 26 passed, 0 failed, 7 filtered (was 5; the two new tests are correctly SLOW-gated).

- [ ] **Step 6: Commit**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && git add invz/invz_critical_T.m invz/tests/test_invz_critical.m && git commit -m "feat(invz): invz_critical_T — Tc(B) bisection mirroring invz_critical

Fixed-B (horizontal) cut through the phase boundary for the regime near
Tc0 where the boundary is nearly parallel to the field axis and the
fixed-T field bisection is ill-conditioned. Same pt.crit classifier and
bracket orientation as invz_critical. Two SLOW-gated tests: crossing
consistency at 1.4 K against invz_critical, and the small-Bx limit
(expected slight undershoot of the closed-form Tc0, R2007 caveat).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Two-regime driver rework

**Files:**
- Modify: `invz/invz_run_phase_diagram.m` (full-file replacement below)

**Interfaces:**
- Consumes: `invz_critical(ion, T, Jf, struct('J0eff', info.Jcc0, 'window', [0.5 7]))` (unchanged) and Task 1's `invz_critical_T(ion, B, Jf, struct('J0eff', info.Jcc0, 'window', [Tsplit Tmax]))`.
- Produces: workspace variables `Bc` (fields at `Ts`, low-T branch), `TcB` (temperatures at `Bs`, high-T branch), `Tc0`, and `phase_boundary` (T-sorted `[T, B]` array of all finite boundary points); a figure with both branches and the closed-form `Tc0` marker. No other file depends on this script.

- [ ] **Step 1: Replace the driver**

Overwrite `invz/invz_run_phase_diagram.m` with exactly:

```matlab
%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
%
% Two-regime search, split by the Tsplit knob below:
%   low-T  (T <= Tsplit): bisect the critical field Bc(T) at each fixed T in
%                         Ts (invz_critical, vertical cuts);
%   high-T (boundary above Tsplit): bisect the critical temperature Tc(B) at
%                         each fixed B in Bs (invz_critical_T, horizontal cuts).
% Near the classical critical point (T -> Tc0 = 1.74 K, B -> 0) the boundary
% is nearly parallel to the field axis, so a vertical cut crosses it at a
% glancing angle: brackets fail and tiny T errors give huge Bc errors. A
% horizontal cut crosses it transversally and is well-conditioned there.
%
% Bracket geometry: at T = Tsplit a point is ordered exactly when
% B < Bc(Tsplit), so every Bs entry below Bc(Tsplit) brackets cleanly in
% [Tsplit, Tmax]; an entry above it fails its bracket assert to NaN without
% affecting other jobs.
%
% Bs starts at 0.1 T for two reasons: (1) invz_twolevel raises
% invz:degenerateDoublet when the field-induced doublet splitting is
% < 1e-4 meV (Bx ~ 0); (2) R2007's small-Bx caveat -- the two-level Sigma
% overestimates the Tc suppression when the doublet splitting is below the
% hyperfine width, so Tc(B->0) slightly undershoots the closed-form Tc0
% endpoint (expected physics, not a bug).
%
% Parallelism: all nT+nB boundary points are INDEPENDENT 1-D bisections, so a
% single flat parfor covers both regimes -- near-linear speedup up to the job
% count. The bisection *inside* one point cannot be parallelized (each bracket
% halving needs the sign of the previous EMT solve). parfor degrades to a
% serial loop automatically when the Parallel Computing Toolbox is absent
% (nWorkers = 0); with it, a local pool is auto-created on first use
% (~10-30 s once). The Jq lattice sum is computed ONCE up front, so workers
% do no disk I/O and never touch the invz/cache -- no contention. Progress
% lines interleave across workers when running in parallel; that is expected.

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Jf = Jnu(:);

% ---- knobs --------------------------------------------------------------
Tsplit = 1.5;   % regime split (K): boundary below Tsplit via Bc(T), above via Tc(B)
Ts = [0.05 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.64 1.68 1.72 1.74 1.76 1.78 1.8];
Bs = [0.1 0.25 0.5 0.75 1.0 1.25 1.5];   % fields (T) for the high-T Tc(B) regime
Tmax = 2.0;     % paramagnetic upper bracket edge for Tc(B) (above Tc0 = 1.74 K)
useParallel = true;                      % false -> force a serial run
% -------------------------------------------------------------------------

Ts = Ts(Ts <= Tsplit);              % master list; the knob trims it
nT = numel(Ts);  nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];           % one independent bisection per job
kind = [ones(1,nT) 2*ones(1,nB)];   % 1 = Bc(T) at fixed T, 2 = Tc(B) at fixed B

nWorkers = 0;                       % 0 = serial (also the no-toolbox fallback)
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

out = nan(1, nT + nB);
parfor (k = 1:nT+nB, nWorkers)
    t0 = tic;
    v = jobs(k);  val = NaN;
    if kind(k) == 1
        try
            % Wide bracket: the low-T points need the upper edge to reach
            % ~4-5 T, while points near Tsplit sit close to the lower edge.
            val = invz_critical(ion, v, Jf, struct('J0eff', info.Jcc0, 'window', [0.5 7]));
        catch err
            fprintf('  T = %.2f K: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    else
        try
            val = invz_critical_T(ion, v, Jf, struct('J0eff', info.Jcc0, 'window', [Tsplit Tmax]));
        catch err
            fprintf('  B = %.2f T: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] B = %.2f T  ->  Tc = %.3f K   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    end
    out(k) = val;
end
Bc  = out(1:nT);                    % low-T branch:  Bc at each Ts
TcB = out(nT+1:end);                % high-T branch: Tc at each Bs

Tc0 = invz_critical_T0field(ion, invz_sigma_crit(info.Jcc0, Jf), info.Jcc0);
fprintf('Zero-field Tc (1/z) = %.3f K  [target 1.74 K]\n', Tc0);

% Merged boundary, T-sorted, finite points only -- workspace export for
% downstream use ('boundary' would shadow the built-in of that name).
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :); %#ok<NASGU>

figure; hold on;
plot(Ts, Bc*10, 'o-');
plot(TcB, Bs*10, 's-');
plot(Tc0, 0, 'ks');
xlabel('T (K)'); ylabel('B_c (kOe)');
title('LiHoF_4 1/z phase boundary (paramagnetic side)');
legend({'B_c(T) bisection (T \leq Tsplit)', 'T_c(B) bisection', 'closed-form T_c(B=0)'}, 'Location', 'southwest');
```

Design notes baked into that code (do not "simplify" them away):
- `jobs`/`kind` arrays exist so every parfor slice `jobs(k)`, `kind(k)`, `out(k)` is indexed by the plain loop variable over its full range. Indexing `Ts(k)` / `Bs(k-nT)` inside branches instead can trip parfor's sliced-variable handling for out-of-range k.
- `Jf = Jnu(:)` is hoisted so workers broadcast the flat vector once.
- The high-T window is `[Tsplit Tmax]` — the same values validated by `test_tc_near_zero_field` in Task 1.

- [ ] **Step 2: Static-check the driver**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "iss = checkcode('invz/invz_run_phase_diagram.m'); arrayfun(@(s) fprintf('L%d: %s\n', s.line, s.message), iss); assert(isempty(iss), '%d checkcode message(s)', numel(iss))"
```

Expected: no `L<line>:` output, exit 0. Any message (especially parfor variable-classification warnings) is a defect — fix it before committing.

- [ ] **Step 3: Run the fast suite (regression guard)**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: PASS — 26 passed, 0 failed, 7 filtered.

(No driver smoke run is required here: the underlying bisections, including the driver's exact `[Tsplit Tmax]` window, are covered by Task 1's slow tests; a real driver run is the production step at the end of this plan.)

- [ ] **Step 4: Commit**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && git add invz/invz_run_phase_diagram.m && git commit -m "feat(invz): two-regime phase-diagram driver (Tsplit knob, flat parfor)

Below Tsplit: Bc(T) field bisection as before. Above: Tc(B) temperature
bisection at each field in the new Bs knob, window [Tsplit Tmax] --
well-conditioned where the boundary is nearly parallel to the field
axis. All points run as one flat parfor job list; merged T-sorted
phase_boundary array exported for downstream use.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Document `invz_critical_T` and the two-regime driver in the manual

**Files:**
- Modify: `invz/README.html` (5 surgical text edits; it is the sole user manual — `invz/README.md` was removed deliberately)

**Interfaces:**
- Consumes: names/signatures exactly as implemented in Tasks 1–2 (`invz_critical_T(ion,Bx,Jnu_flat,opts)`, knobs `Tsplit`, `Bs`, `Tmax`).
- Produces: documentation only.

- [ ] **Step 1: Update the driver section**

Find this paragraph (in the `<h3><code>invz_run_phase_diagram.m</code>` section, around line 306):

```html
<p>Bisects for the critical field at each of a list of temperatures and plots the R 2007 Fig. 1 paramagnetic-side boundary, plus the closed-form zero-field \(T_c\). Full single-core run ≈ 1–2 h; see the next section for the parallel speedup.</p>
```

Replace with:

```html
<p>Traces the R 2007 Fig. 1 paramagnetic-side boundary in <strong>two regimes</strong>, split by the <code>Tsplit</code> knob at the top of the script (default 1.5 K). Below <code>Tsplit</code> it bisects the critical <em>field</em> \(B_c(T)\) at each temperature in <code>Ts</code> (<code>invz_critical</code>, vertical cuts). Above — where the boundary runs nearly parallel to the field axis, so a fixed-\(T\) cut crosses it at a glancing angle, brackets fail, and tiny \(T\) errors give huge \(B_c\) errors — it bisects the critical <em>temperature</em> \(T_c(B)\) at each field in <code>Bs</code> (<code>invz_critical_T</code>, horizontal cuts, window <code>[Tsplit, Tmax]</code>). Keep <code>Bs</code> entries ≥ 0.1 T (the doublet is degenerate at \(B_x \to 0\), and R 2007's small-\(B_x\) caveat means \(T_c(B \to 0)\) slightly undershoots the closed-form endpoint) and below \(B_c(T_{\mathrm{split}})\) (bracket geometry; an entry above it fails cleanly to NaN). The closed-form zero-field \(T_c\) is plotted as before. Full single-core run ≈ 1–2 h; see the next section for the parallel speedup.</p>
```

- [ ] **Step 2: Update the parallelism table row**

Find (around line 313):

```html
<tr><td><code>invz_run_phase_diagram</code></td><td>temperature (<code>parfor</code> over <code>Ts</code>)</td><td><code>useParallel = true</code> (default)</td></tr>
```

Replace with:

```html
<tr><td><code>invz_run_phase_diagram</code></td><td>boundary points (one flat <code>parfor</code> over the <code>Ts</code> and <code>Bs</code> jobs)</td><td><code>useParallel = true</code> (default)</td></tr>
```

- [ ] **Step 3: Add the function-reference row**

Immediately AFTER this row in the section-6 reference table (around line 340):

```html
<tr><td><code>invz_critical</code></td><td><code>bx = invz_critical(ion,T,Jnu_flat,opts)</code></td><td>critical field \(B_c(T)\), bisection on <code>pt.crit</code></td></tr>
```

insert:

```html
<tr><td><code>invz_critical_T</code></td><td><code>tc = invz_critical_T(ion,Bx,Jnu_flat,opts)</code></td><td>critical temperature \(T_c(B)\), the same bisection over \(T\) at fixed field (near-vertical part of the boundary)</td></tr>
```

- [ ] **Step 4: Update the first SVG flowchart's phase-boundary box**

Two `<text>` lines in the tier-2 phase-boundary band (around lines 168–169). Find:

```html
  <text class="fc-fn" x="408" y="101">invz_critical</text>
  <text class="fc-note" x="408" y="119">bisect B&#8339; for the boundary at each T</text>
```

Replace with:

```html
  <text class="fc-fn" x="408" y="101">invz_critical (&middot;_T)</text>
  <text class="fc-note" x="408" y="119">bisect B&#8339; at fixed T &middot; T at fixed B&#8339;</text>
```

(The `(&middot;_T)` suffix matches the manual's existing `invz_solve_point (&middot;_ordered)` idiom. The second, band-7 flowchart is intentionally untouched: its `invz_critical` label is already width-constrained between neighbours, and the reference table is the authoritative catalog.)

- [ ] **Step 5: Verify the edits**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && grep -c "invz_critical_T<" invz/README.html && grep -c "Tsplit" invz/README.html
```

Expected: `2` then `1` (`grep -c` counts LINES: the driver paragraph and the new table row each match the first pattern once; only the driver paragraph mentions `Tsplit`). Then eyeball the render — flowchart text must not overflow its rounded box:

```bash
open "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion/invz/README.html"
```

- [ ] **Step 6: Commit**

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && git add invz/README.html && git commit -m "docs(invz): document invz_critical_T and the two-regime phase-diagram driver

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Final verification (after all tasks)

- [ ] Full slow suite, warm cache (~5–15 min):

```bash
cd "/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion" && "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "setenv('INVZ_SLOW','1'); results = runtests('invz/tests'); assertSuccess(results)"
```

Expected: all pass, 0 failed (33 passed / 0 filtered).

- [ ] Production run (user-initiated — hours-scale, NOT part of the implementation loop): run `invz_run_phase_diagram` in MATLAB and check the spec's acceptance criteria: finite `TcB` for every default `Bs` entry, monotonically decreasing with B; the two branches join smoothly near (Tsplit, Bc(Tsplit)); the high-T branch drops near-vertically to the `Tc0` marker at (≈1.74 K, 0), reproducing R 2007 Fig. 1.
