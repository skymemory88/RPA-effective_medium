# Tensor-branch drivers: `invzt_run_phase_diagram` / `invzt_run_spectra`

Date: 2026-07-18. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-18); scope and the peak-picker sharing decision were
each confirmed via an explicit multiple-choice pick — see Decisions below.
REVISED 2026-07-18 after a source-verification pass (subagent review) that
found a runtime-fatal defect in the first draft's `Tc0`-proxy handle plus
four smaller issues — see "Revision notes" at the end.

## Decisions (user-approved, govern over any conflicting prose below)

1. **Scope: tensor-only drivers, PM-side, as-is.** Ship two new scripts in
   `invz_tensor/`, honestly scoped to what `invzt_solve_point` /
   `invzt_critical` / `invzt_chi_realaxis` support today. No new tensor
   solver physics (no ordered-phase branch, no `invzt_critical_T`) as a
   prerequisite. `invz_projected`'s own drivers and module code stay
   unmodified except for the one `git mv` in item 2.
2. **`invz_peak_energy.m` moves to `invz_common/`** (`git mv`, zero logic
   change) so both `invz_projected` and the new `invzt_run_spectra.m` call
   the one shared copy, matching the precedent already set by the 16
   single-ion functions moved there for the same reason.

## Problem

`invz_projected` ships `invz_run_phase_diagram.m` and `invz_run_spectra.m` —
scripts a user just runs to get a phase diagram or a susceptibility
spectrum. `invz_tensor` has no equivalents: only the point-solver primitives
a driver would call (`invzt_solve_point`, `invzt_critical`,
`invzt_chi_realaxis`). A user on the tensor branch currently has to hand-roll
the sweep loop and plotting every time (see `invz_tensor/README.html` §2,
added as a stopgap in the prior session).

Porting the projected drivers verbatim is not possible: the tensor branch
has no ordered/FM solve at all (`invzt_solve_point` only ever returns a
paramagnetic fixed point — LOCKED convention 7) and no temperature-cut root
finder (`invzt_critical_T` does not exist, only the field-cut
`invzt_critical`). Both gaps are real physics/numerics work, not driver
plumbing, and are out of scope here (Decision 1).

## Solution overview

Two new files, `invz_tensor/invzt_run_phase_diagram.m` and
`invz_tensor/invzt_run_spectra.m`, structurally mirroring their projected
counterparts but reduced to the PM-side field-cut-only physics the tensor
branch actually has, plus a few tensor-only knobs (`mode`/`nlevels`/`dress`)
the projected driver has no equivalent of. One shared numeric helper
(`invz_peak_energy.m`) moves to `invz_common/` so `invzt_run_spectra.m` can
use it without depending on `invz_projected`.

Both drivers stay **self-contained by default**: neither requires
`invz_projected` on the path unless the user opts into the cross-model
comparison anchor (Component 1, `show_projected_anchor`).

## Verified facts this design rests on

Checked against source during the revision pass; no longer assumptions:

- **The `invz_peak_energy.m` move is safe.** All ten tests that invoke
  `invz_spectra_map` / `invz_spectra_qpath` / `invz_chi_tensor_ref` share an
  identical `setupOnce` that does
  `addpath(fullfile(here, '..', '..', 'invz_common'))`, as do the three
  driver callers (`invz_run_spectra.m`, `invz_run_ip_scan.m`,
  `invz_run_tensor_ref.m`) and the tensor interop test
  `tests/interop/test_invzt_rpa_parity.m` (which reaches it transitively via
  `invz_chi_tensor_ref`). `invz_common/` has no name collision.
- **`imag(out.chi_uniform(3,3,:))` is the correct, already-positive χ''**
  — no sign flip. The projected chain un-flips its internal `G`-language at
  `chit = (-G0)./(1+Sw)`, and `invz_spectra_map.m:143` builds the plotted
  `S.chiz` as `imag(o.chi_cc_q(1,:)).'` with no extra minus; the tensor
  chain (`invzt_chi0_split` → `invzt_chi_rpa` → `u'*X*u`) introduces no
  additional sign. `test_invzt_rpa_parity.m` compares both un-negated.
- **`invzt_tc_pm_extrap` does NOT guard its `critfun`.** It calls
  `[c(i), okv(i)] = critfun(Tg(i))` bare (line 46); any exception propagates
  and aborts the whole call. The projected analogue's effective critfun,
  `invz_crit_at.m`, supplies its own *selective* try/catch. The driver must
  do the same — see Component 1.
- **`invzt_solve_point` throws on a realistic sweep.** It contains no
  try/catch. Sweep-plausible identifiers: `invzt:a1ZeroField` (transverse
  field below 1e-6 T), and — propagated from `invz_common/invz_twolevel.m`
  — `invz:degenerateDoublet` (splitting < 1e-4 meV) and `invz:orderedPhase`
  (|m| > 1e-3). Note that a *failed PM fixed point* deep in the ordered
  phase does NOT throw: it returns `pt.converged = false` cleanly.
- **`invz_odd_zero_field` requires `ion.demag == 0`** for every `mode`
  (including `'off'`) — the guard lives in its callee `invz_odd_blocks.m`
  (`invz:oddDemag`) and fires before any mode branch. Its `mode` values are
  `'off' | 'full' | 'uniform_only' | 'qstruct_only'`; `mode='full'` solves
  **seven** variants per grid, so it is not a cheap anchor.

## Components

### 1. `invz_tensor/invzt_run_phase_diagram.m`

```matlab
%INVZT_RUN_PHASE_DIAGRAM  Full-tensor 1/z PM-side phase boundary Bc(T).
%
% SCOPE (read this before trusting the plot): FIELD-CUTS ONLY. This driver
% bisects Bc at each fixed T via invzt_critical. There is no temperature-cut
% finder on the tensor branch (no invzt_critical_T), so the near-vertical
% part of the boundary approaching Tc(0) is NOT resolved -- those Ts points
% simply fail their bracket and come back NaN. That region is left BLANK by
% design, the same way invz_run_phase_diagram's ODD-overlay branch leaves its
% own near-Tc0 region blank. There is also no ordered-phase tensor solve, so
% nothing below the boundary is computed at all.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs ------------------------------------------------------------------
Ts     = linspace(1.0, 2.2, 25);       % temperature grid (K)
Brange = [0.05 6.0];                   % [Blo Bhi] bisection bracket (T). Blo > 0:
                                       % mode 'a1' forbids exact zero transverse field.
gridN  = 16;  gridConv = 'halfopen';   % invzt_qgrid(gridN, gridConv)
dpRng  = 30;                           % invzt_jq_tensor coupling-sum range
useParallel = true;                    % false -> force serial
solve_opts  = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', 'dress', 'full');

show_proxy_anchor     = true;          % tensor-native small-Bx Tc(0) proxy
                                       % (invzt_tc_pm_extrap). Costs ~numel(Ts)
                                       % EXTRA serial point solves -- set false to skip.
Bproxy                = 0.05;          % transverse field (T) for that proxy
show_projected_anchor = false;         % OPT-IN cross-model comparator: the PROJECTED
                                       % closed-form Tc(0) (invz_odd_zero_field).
                                       % Off by default because it (a) puts
                                       % invz_projected on the path, (b) requires
                                       % ion.demag == 0, (c) solves 7 variants per
                                       % grid, and (d) is a DIFFERENT ODD treatment
                                       % (projected Tier-1 chi_perp-mediated dJ) from
                                       % this driver's structural tensor `odd` flag --
                                       % a comparator, NOT the same quantity.
% -----------------------------------------------------------------------------

if show_projected_anchor
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_projected'));
    assert(ion.demag == 0, 'invzt:anchorDemag', ...
        ['show_projected_anchor requires ion.demag == 0 (invz_odd_blocks is ' ...
         'intrinsic-only). Set ion.demag = 0 or disable the anchor.']);
end

g   = invzt_qgrid(gridN, gridConv);
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

nWorkers = 0;
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

nT = numel(Ts);
Bc = nan(1, nT);
parfor (k = 1:nT, nWorkers)
    t0 = tic;  val = NaN;
    % val-then-assign (NOT `Bc(k) = ...` inside the try): keeps the sliced
    % output variable unconditionally assigned on every iteration, matching
    % invz_run_phase_diagram's own parfor pattern.
    try
        val = invzt_critical(ion, Ts(k), lat, Brange, solve_opts);
    catch err
        % invzt:bracket is a genuine per-point outcome (the window did not
        % bracket a crossing at this T -- expected near Tc0, where the
        % boundary is near-vertical). Anything else (invzt:mode, a bad lat,
        % ...) is a MISCONFIGURATION that would otherwise silently produce a
        % whole sweep of NaNs, so it rethrows.
        if ~strcmp(err.identifier, 'invzt:bracket'), rethrow(err); end
        fprintf('  T = %.2f K: no bracket in [%.2f %.2f] T (left NaN)\n', ...
                Ts(k), Brange(1), Brange(2));
    end
    fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT, Ts(k), val, toc(t0));
    Bc(k) = val;
end

% ---- zero-field anchors ------------------------------------------------------
Tc0_proxy = NaN;  Tc0_closed = NaN;
if show_proxy_anchor
    % Pass the FULL Ts: invzt_tc_pm_extrap filters to converged-PM points
    % itself and extrapolates the two lowest, so no hand-picked cutoff.
    critfun = @(T) invzt_crit_at(ion, T, [Bproxy 0 0], lat, solve_opts);
    try
        Tc0_proxy = invzt_tc_pm_extrap(critfun, Ts);
    catch err
        if ~strcmp(err.identifier, 'invzt:tcNoWindow'), rethrow(err); end
        fprintf('  Tc(0) proxy: no PM extrapolation window on this Ts grid (skipped)\n');
    end
end
if show_projected_anchor
    % Track this driver's own odd setting -- comparing a tensor ODD-on curve
    % against a projected ODD-off anchor (or vice versa) would be misleading.
    anchor_mode = 'off';
    if ~isfield(solve_opts, 'odd') || ~isequal(solve_opts.odd, false), anchor_mode = 'full'; end
    Tc0_closed = invz_odd_zero_field(ion, struct('mode', anchor_mode, ...
        'grids', {{gridN}}, 'dpRng', dpRng, 'cache', true));   % same grid/dpRng as lat
end

figure; hold on;
plot(Bc, Ts, 'o-', 'DisplayName', 'tensor A1: field-cut B_c(T)');
if isfinite(Tc0_proxy)
    plot(0, Tc0_proxy, 'k^', 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('tensor small-B_x proxy T_c(0), %.2f T', Bproxy));
end
if isfinite(Tc0_closed)
    plot(0, Tc0_closed, 'ks', 'MarkerFaceColor', 'c', ...
         'DisplayName', 'projected closed-form T_c(0) (comparator)');
end
xlabel('B_c (T)'); ylabel('T (K)');
title('LiHoF_4 full-tensor 1/z phase boundary (PM side, field-cuts only)');
legend('Location', 'southwest');

phase_boundary = sortrows([Ts(:) Bc(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);

% =============================================================================
function [c, ok] = invzt_crit_at(ion, T, B, lat, opts)
% Tensor analogue of invz_projected's invz_crit_at: one criticality sample,
% with the SELECTIVE try/catch invzt_tc_pm_extrap does not provide itself.
% Physics signals (degenerate doublet, ordered toy doublet, zero-field guard)
% become ok = false; structural/config errors rethrow.
try
    pt = invzt_solve_point(ion, T, B, lat, opts);
    c  = pt.crit;
    % Same three-part sample-validity rule as invzt_critical's header.
    ok = pt.converged && isfinite(c) && pt.Sigma0 >= -0.5;
catch err
    switch err.identifier
        case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
            c = NaN;  ok = false;      % phase/physics signal, not an error
        otherwise
            rethrow(err);              % invzt:mode, invzt:a1OddGamma, ... = misconfiguration
    end
end
end
```

### 2. `invz_tensor/invzt_run_spectra.m`

```matlab
%INVZT_RUN_SPECTRA  Full-tensor 1/z chi''_cc spectra vs transverse field, or along a q-path.
%
% SCOPE: PARAMAGNETIC SIDE ONLY. invzt_solve_point has no ordered-phase branch,
% so unlike invz_run_spectra this driver cannot sweep ACROSS Bc to show the soft
% mode hardening on both sides. Fields that land on the ordered side (or fail the
% sample-validity rule) are masked to NaN with a console note.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs (mirroring invz_run_spectra's names where the concept carries over) --
T      = 1.6;                          % K
fields = linspace(0.3, 4.5, 40);       % T -- keep on the PM side (see SCOPE above)
w      = linspace(0, 0.6, 400).';      % meV, uniform spacing (invz_peak_energy assumes it)
eta    = 5e-3;                         % meV, real-axis Lorentzian HWHM
sliceMax  = 6;                         % field count at/below which line slices are drawn
peak_wmin = 0.02;                      % meV -- lower bound for the peak pick. NON-ZERO by
                                       % default to exclude the low-frequency hyperfine
                                       % line, mirroring invz_spectra_map's opts.peak_wmin;
                                       % set 0 to pick over the whole w grid.
theta_c = 0.0;  phi_ab = 0.0;          % field-direction tilt knobs (deg). Ported as-is:
                                       % invzt_solve_point already takes a full [Bx By Bz]
                                       % and forwards transverse_mf.
transverse_mf = 'legacy_x';            % 'legacy_x' | 'none' | 'vector_ab'
gridN = 16; gridConv = 'halfopen'; dpRng = 30;
useParallel = true;
solve_opts = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', 'dress', 'full', ...
                    'transverse_mf', transverse_mf);

% ---- q-path view: set qpath non-empty to switch views ------------------------
qpath = [];                            % [] = field sweep; [nq x 3] r.l.u. = q-path view
% qpath = [linspace(0, 0.5, 41).' zeros(41, 2)];   % Gamma -> (0.5,0,0)
Bq    = 2.0;                           % T -- fixed field for the q-path view
wq    = linspace(0, 0.6, 400).';       % meV -- own frequency grid for the q-path view
% -----------------------------------------------------------------------------

if phi_ab ~= 0 && strcmp(transverse_mf, 'legacy_x')
    error('invzt_run_spectra:transverseMF', ...
        ['phi_ab = %.3g deg needs the vector transverse mean field: set transverse_mf ' ...
         'to ''vector_ab'' (or ''none'' for a bare CF+Zeeman diagnostic). legacy_x is ' ...
         'x-only and C4-inconsistent for rotated fields.'], phi_ab);
end
dhat = [cosd(theta_c)*cosd(phi_ab), cosd(theta_c)*sind(phi_ab), sind(theta_c)];

g   = invzt_qgrid(gridN, gridConv);
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

if ~isempty(qpath)
    % ---------------- q-path dispersion at one fixed field --------------------
    pt = invzt_solve_point(ion, T, Bq*dhat, lat, solve_opts);
    assert(pt.converged && isfinite(pt.crit) && pt.crit > 0, 'invzt:qpathOrdered', ...
        'B = %.2f T at T = %.2f K is not a valid PM point (crit = %.4g).', Bq, T, pt.crit);
    out = invzt_chi_realaxis(ion, T, Bq*dhat, pt, wq, ...
            struct('qsel', qpath, 'dpRng', dpRng, 'eta', eta));
    chipp_q = imag(out.chi_cc_q);                 % [nq, nw], already positive chi''
    figure; imagesc(wq, 1:size(chipp_q, 1), chipp_q); set(gca, 'YDir', 'normal');
    xlabel('\omega (meV)'); ylabel('q index along path'); colorbar;
    title(sprintf('tensor 1/z \\chi''''_{cc}(q,\\omega), T = %.2f K, B = %.2f T', T, Bq));
    Epeak_q = invz_peak_energy(chipp_q.', wq, peak_wmin);   % columns must be per-q
    figure; plot(1:numel(Epeak_q), Epeak_q, 'o-');
    xlabel('q index along path'); ylabel('E_{peak} (meV)');
else
    % ---------------- field sweep at the uniform mode -------------------------
    nWorkers = 0;
    if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end
    nf = numel(fields);
    chipp = nan(numel(w), nf);
    parfor (k = 1:nf, nWorkers)
        col = nan(numel(w), 1);
        try
            pt = invzt_solve_point(ion, T, fields(k)*dhat, lat, solve_opts);
            % Same three-part sample-validity rule as invzt_critical's header
            % (converged, finite positive crit, Sigma0 floor) -- reused rather
            % than a looser ad hoc check, so the spurious negative-Sigma fixed
            % point invzt_critical warns about never reaches the plot.
            if pt.converged && isfinite(pt.crit) && pt.crit > 0 && pt.Sigma0 >= -0.5
                o = invzt_chi_realaxis(ion, T, fields(k)*dhat, pt, w, ...
                        struct('qsel', 'gamma_uniform', 'eta', eta));
                col = squeeze(imag(o.chi_uniform(3, 3, :)));   % already positive chi''
            else
                fprintf('  B = %.2f T: ordered / non-converged / invalid (masked)\n', fields(k));
            end
        catch err
            switch err.identifier
                case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
                    fprintf('  B = %.2f T: %s (masked)\n', fields(k), err.identifier);
                otherwise
                    rethrow(err);
            end
        end
        chipp(:, k) = col;
    end

    if nf <= sliceMax
        figure; hold on;  co = lines(nf);
        for k = 1:nf
            plot(w, chipp(:, k), '-', 'Color', co(k, :), ...
                 'DisplayName', sprintf('%.2f T', fields(k)));
        end
        xlabel('\omega (meV)'); ylabel('\chi''''_{cc}');
        title(sprintf('tensor 1/z, T = %.2f K', T)); legend show;
    else
        figure; imagesc(fields, w, chipp); set(gca, 'YDir', 'normal');
        xlabel('B (T)'); ylabel('\omega (meV)'); colorbar;
        title(sprintf('tensor 1/z \\chi''''_{cc}(B,\\omega), T = %.2f K', T));
    end

    Epeak = invz_peak_energy(chipp, w, peak_wmin);
    figure; plot(fields, Epeak, 'o-');
    xlabel('B (T)'); ylabel('E_{peak} (meV)');
    title(sprintf('\\chi''''_{cc} peak energy vs field, T = %.2f K', T));
    % Gaps are CENSORED peaks (boundary max / non-positive column) or masked
    % ordered points -- do not interpolate over them.
end
```

The q-path view is a real branch selected by `qpath`, with its own `Bq`/`wq`
knobs, rather than commented-out scaffolding — the first draft's commented
block referenced an undefined `Bq` and reused the field-sweep `w`, so it
would not actually have run if uncommented.

### 3. `invz_peak_energy.m` → `invz_common/`

`git mv invz_projected/invz_peak_energy.m invz_common/invz_peak_energy.m`,
no content changes. Verified safe (see "Verified facts"): every test and
driver that reaches it already has `invz_common` on the path, and there is no
name collision. No caller edits are expected; the PROJECTED regression gate
confirms it.

## Error handling

Summary of the conventions used above, all inherited from existing module
practice rather than newly invented:

- **Per-point boundary failure** (`invzt:bracket`): reported, left `NaN`,
  sweep continues. Every other identifier rethrows, so a misconfiguration
  fails loudly on the first point instead of producing a silent all-NaN plot.
- **Physics signals** (`invz:degenerateDoublet`, `invz:orderedPhase`,
  `invzt:a1ZeroField`): absorbed into `ok = false` / a masked column, never
  an error — mirroring `invz_crit_at.m`'s selective pattern.
- **Invalid-but-non-throwing samples** (`pt.converged == false`, non-finite
  or non-positive `crit`, `Sigma0` below the `-0.5` spurious-fixed-point
  floor): masked with a console note, reusing `invzt_critical`'s own
  three-part sample-validity rule verbatim.
- **`invzt:tcNoWindow`** from the `Tc0` proxy: caught, anchor skipped, rest
  of the plot still renders.
- **`phi_ab ~= 0` with `transverse_mf = 'legacy_x'`**: hard error, ported
  verbatim from the projected driver (pure input validation).

## Testing / verification

Neither new file gets a MATLAB unit test — matching precedent: the existing
`invz_run_phase_diagram.m` / `invz_run_spectra.m` have none either (drivers
are verified by running them, not by `runtests`). Verification for this task:

1. Run both scripts at reduced size (`gridN = 8`, `numel(Ts) ~ 4`,
   `numel(fields) ~ 5`, `useParallel = false`) via `matlab -batch` from the
   repo root; confirm they complete without error and produce non-empty
   figures. Exercise the q-path branch too (set `qpath` non-empty) — it is
   now a real code path, not a comment.
2. Exercise the opt-in anchor once (`show_projected_anchor = true`) to
   confirm the conditional `addpath` + `ion.demag` assert behave.
3. Re-run the CORE suite (`runtests('invz_tensor/tests')`, expect 47/0/1)
   and the PROJECTED regression gate (`runtests('invz_projected/tests')`,
   expect 143/0/19) after the `invz_peak_energy.m` move.
4. Update `invz_tensor/README.html` §2 ("Getting a phase diagram or a
   susceptibility spectrum") to point at the real scripts instead of the
   "no turnkey driver" framing, keeping the existing recipe prose as the
   "what these scripts do under the hood" explanation.

Full-resolution production sweeps (the real `gridN = 16`+ / dense-`Ts` runs)
are left to the user to kick off, same as the projected driver's own
"production sweep... LEFT TO THE USER" note.

## Success criteria

- `invzt_run_phase_diagram.m` and `invzt_run_spectra.m` exist in
  `invz_tensor/`, run end-to-end at reduced size with no runtime errors
  (both branches of the spectra driver), and produce the plots described.
- Both drivers run with `invz_projected` absent from the path under their
  default knobs.
- `invz_peak_energy.m` lives in `invz_common/`; `invz_projected/tests`
  (143/0/19) and `invz_tensor/tests` (47/0/1) both still pass unchanged.
- `invz_tensor/README.html` §2 reflects the new scripts.

## Out of scope

- `invzt_critical_T` (temperature-cut root finder) and any ordered-phase
  tensor solver — the actual blockers to full parity; explicitly deferred
  per Decision 1, not attempted as part of this task.
- Any change to `invz_projected`'s own drivers beyond the one `git mv`.
- Unifying the two branches' drivers behind a shared `invz_common` entry
  point / model toggle — considered and explicitly declined for now (the
  adapter work is only worth it once the tensor branch has FM + T-cut
  support; revisit then).
- Warm-starting `Sigma_seed` across sweep points (the projected drivers do
  not do this either; `invzt_critical` warm-starts only *within* one
  bisection).

## Revision notes (2026-07-18, post source-verification)

The first draft was reviewed against source. Changes:

1. **Fixed a runtime-fatal defect.** The draft's `Tc0` proxy used
   `@(T) crit_ok(invzt_solve_point(...))`, where `crit_ok` received an
   already-computed `pt` — so a throw inside `invzt_solve_point` escaped
   before any guard ran, and `invzt_tc_pm_extrap` has no try/catch of its
   own. Replaced with a proper `invzt_crit_at` local function carrying the
   selective try/catch, mirroring `invz_crit_at.m`.
2. **Inverted the error-catch rule.** The draft caught *all* `invzt:*` and
   NaN-ed, which would swallow config errors (`invzt:mode`) into a silent
   all-NaN sweep. Now only `invzt:bracket` is absorbed.
3. **Fixed the `parfor` slicing pattern** to the val-then-assign form the
   projected driver uses, so the sliced output is always assigned.
4. **Corrected the projected anchor**: made it opt-in (default off) with a
   conditional `addpath`, added the required `ion.demag == 0` assert, and
   tied its `mode` to the driver's own `odd` setting instead of hardcoding
   `'full'`. Also documented that it is a *different* ODD treatment — a
   comparator, not the same quantity.
5. **Smaller:** dropped the magic `Ts(Ts <= 1.6)` cutoff (pass the full
   grid; the function filters internally); made `peak_wmin` a knob with a
   non-zero default so the hyperfine line does not win the peak pick;
   promoted the q-path view from broken commented scaffolding to a real
   `qpath`-selected branch with its own `Bq`/`wq`; added `parfor` +
   `eta` to the spectra driver for parity with the projected one.
