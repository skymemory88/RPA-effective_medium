%INVZ_RUN_PHASE_DIAGRAM
% Two-regime search: Bc(T) at each fixed T in Ts (invz_critical, vertical field cuts), and
% Tc(B) at each fixed B in Bs (invz_critical_T, horizontal temperature cuts, self-adapting window).

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
% ion.demag = 1;  ion.alpha = 0.25;  % OPTIONAL sample-shape (demagnetization) knob; default off.
%   Ordering-channel criticality (info.Jcc0/Jnu) and Tc(B=0) are demag-INVARIANT (R2007); Bc(T)
%   vs APPLIED field can shift via the demag-aware info.Jaa0 (hoisted into Jxx0 below). See
%   invz_ion.m for the full explanation. demag = 0 (default) matches the R2007 benchmark.

ion.odd = 1;  % (off-diagonal dipolar) Tier-1 ODD knob

% ---- knobs ------------------------------------------------------------------
% Standard sweep (used when overlay_quick is false/unset; see below):
% Tc(B) search window is not set here: invz_critical_T self-adapts per field
% (anchors at Tc0+0.05 K, spans 0.5 K down; see its header for details).
% Ts = linspace(0.05, 1.95, 28);   % Temperature points to sample
% Bs = [];   % magnetic field points to sample
Ts = []; % Temperature points to sample
Bs = linspace(0.05,5,26); % magnetic field points to sample
useParallel = true;             % false -> force a serial run

% V4.1 quick ODD overlay (opt-in; see the header block below for what this
% does). Set `overlay_quick = true` in the base workspace before running this
% script to use these instead of the standard sweep above.
if ~exist('overlay_quick', 'var'), overlay_quick = false; end
Tq  = [0.8 1.2 1.4];              % B-cut temperatures (all < Tc_ODD(0) = 1.509 K)
win = [0.1 6];                    % base field window; top (6 T) is PM at every Tq
t12_budget_s = 15*60;            % Tier-2 wall-clock guard (drop remaining pts if over)
% -------------------------------------------------------------------------------

[Jf, info, Jxx0] = invz_bz_couplings(ion);   % shared BZ-grid coupling branches (Jaa0-aware)
J0 = info.Jcc0;   % scalar hoist: avoids broadcasting the whole info struct to workers
% Jxx0 is demag-aware; at demag = 0 it differs from the hardcoded ion.Jxx0 by <0.1% -- the live
% dipole sum supersedes the pasted constant. Tc0 below needs no Jxx0: at B = 0, <Jx> = 0.

% ===========================================================================
% V4.1 -- QUICK ODD OVERLAY (opt-in; the standard path below is byte-identical
% when overlay_quick is unset/false). Set `overlay_quick = true` in the base
% workspace before running this script to build the headline ODD phase-boundary
% overlay instead of the single-curve sweep. It draws FOUR boundary curves --
% 1/z baseline, 1/z + Tier-1 ODD, 1/z + Tier-1+2 ODD, and mean-field (RPA) --
% plus the closed-form zero-field endpoints (invz_odd_zero_field, Richardson
% 12/24) and the R 2007 Fig. 1 experimental anchors, with a second panel of the
% critical self-energy Sigma(0) along the boundary with/without ODD.
%
% Method (ODD-LOG V4.1): B-CUT points ONLY (invz_critical), at a few temperatures
% below Tc_ODD(0) = 1.509 K. invz_critical survives with ODD on through its
% `para_edge` fallback (the ordered side never re-converges with ODD on, ODD-LOG
% T2.2); invz_critical_T CANNOT bracket with ODD on (no metastable PM window
% below the boundary -> it lacks the para_edge fallback), so the near-Tc0 T-cut
% region of the ODD curves is LEFT BLANK here -- the missing `para_edge` analog
% for the T-cut finder is a documented V4-scope item (ODD-LOG T2), not fixed
% here. Runs in-session on warm caches (~5-15 min; the Tier-2 B-cuts are the
% slow pole -- if they exceed a ~15 min budget the remaining Tier-2 points are
% dropped and noted). The PRODUCTION sweep (dense T grid, full para-edge
% boundary, hours) is LEFT TO THE USER: keep overlay_quick = false, set
% ion.odd = 1 (and/or opts.odd_tier2), widen Ts, and run the standard parfor
% path once per config -- repo precedent for hours-long sweeps.
if overlay_quick
    assert(ion.demag == 0, 'invz:oddDemag', 'ODD overlay is intrinsic-only (requires demag = 0).');

    % --- ODD geometric blocks ONCE, pre-loop (P0.4), on the 16^3 driver grid ----
    Sodd  = build_driver_odd_blocks(ion);
    oOff  = struct('J0eff', J0, 'Jxx0', Jxx0);                        % 1/z baseline solver opts
    oT1   = struct('J0eff', J0, 'Jxx0', Jxx0, 'odd', true, 'odd_blocks', Sodd);
    oT12  = oT1;  oT12.odd_tier2 = true;

    nT = numel(Tq);
    Bc_off = nan(1, nT);  Bc_t1 = nan(1, nT);  Bc_t12 = nan(1, nT);  Bc_mf = nan(1, nT);
    Sg_off = nan(1, nT);  Sg_t1 = nan(1, nT);  Sg_t12 = nan(1, nT);
    t12_used = 0;  t12_drop = false;
    for it = 1:nT
        T = Tq(it);
        % 1/z baseline (real ordered/PM crossing).
        Bc_off(it) = try_bc(@() invz_critical(ion, T, Jf, setwin(oOff, win)));
        Sg_off(it) = sigma_at(ion, T, Bc_off(it), Jf, oOff);
        % 1/z + Tier-1 ODD (para_edge fallback). Narrow the window around the off
        % boundary to bound the coarse scan (T1 boundary < off boundary).
        wT1 = [max(win(1), Bc_off(it) - 1.4), Bc_off(it) + 0.3];
        Bc_t1(it) = try_bc(@() invz_critical(ion, T, [], setwin(oT1, wT1)));
        Sg_t1(it) = sigma_at(ion, T, Bc_t1(it), [], oT1);
        % Mean-field boundary: crit_MF = 1 - J0*chi0cc0 = pt.crit - pt.Sigma0.
        Bc_mf(it) = mf_boundary(ion, T, Jf, oOff, J0, win(2), Bc_off(it));
        % 1/z + Tier-1+2 ODD -- the slow pole. Tight window around the Tier-1
        % boundary (Tier-2 boundary sits at/just below Tier-1, ~97:3 split).
        if ~t12_drop && isfinite(Bc_t1(it))
            tt = tic;
            wT12 = [max(win(1), Bc_t1(it) - 0.4), Bc_t1(it) + 0.6];
            oT12w = setwin(oT12, wT12);  oT12w.fieldstep = 0.15;
            Bc_t12(it) = try_bc(@() invz_critical(ion, T, [], oT12w));
            Sg_t12(it) = sigma_at(ion, T, Bc_t12(it), [], oT12);
            t12_used = t12_used + toc(tt);
            if t12_used > t12_budget_s && it < nT
                t12_drop = true;
                warning('invz:oddOverlayT2Budget', ...
                    ['Tier-2 B-cuts used %.1f min after %d point(s) (> %.0f min budget); ' ...
                     'dropping the remaining Tier-2 points (ODD-LOG V4.1).'], t12_used/60, it, t12_budget_s/60);
            end
        end
        fprintf('  overlay T = %.2f K: Bc off/T1/T12/MF = %.3f / %.3f / %.3f / %.3f T\n', ...
                T, Bc_off(it), Bc_t1(it), Bc_t12(it), Bc_mf(it));
    end

    % --- closed-form zero-field endpoints (Richardson 12/24, ODD-LOG T1.5) ------
    [Tc0_off,  zOff]  = invz_odd_zero_field(ion, struct('mode', 'off'));
    [Tc0_full, zFull] = invz_odd_zero_field(ion, struct('mode', 'full'));
    Tc0_mf = invz_critical_T0field(ion, 0, J0);        % mean-field Tc(0) (Sc = 0) ~ 2.26 K
    Sc0_off = zOff.Sc_rich;  Sc0_full = zFull.Sc_rich; % Sigma_c(0): 0.298 (off) / 0.389 (ODD)

    % --- R 2007 Fig. 1 experimental anchors (Bitko 1996 / Babkevich 2016, as cited
    % in Rönnow 2007 Fig. 1): Tc(0) = 1.53 K; Hc(T -> 0) ~ 4.9-5.0 T. Two anchor
    % points only (the module ships no digitized experimental curve). --------------
    exp_pts = [1.53 0; 0.10 4.95];                     % [T(K) Hc(T)]

    fig = figure('Color', 'w', 'Position', [100 100 760 780]);
    % Panel 1: phase boundary B_c(T).
    ax1 = subplot(2, 1, 1);  hold(ax1, 'on');  box(ax1, 'on');
    plot(ax1, [Tq Tc0_mf],   [Bc_mf   0], 'k:',  'LineWidth', 1.2, 'DisplayName', 'mean field (\Sigma=0), T_c(0)=2.26 K');
    plot(ax1, [Tq Tc0_off],  [Bc_off  0], 'o-',  'Color', [0 0.45 0.74], 'LineWidth', 1.4, 'DisplayName', '1/z baseline (ODD off)');
    plot(ax1, [Tq Tc0_full], [Bc_t1   0], 's-',  'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, 'DisplayName', '1/z + Tier-1 ODD');
    plot(ax1, Tq,            Bc_t12,      '^--', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.4, 'DisplayName', '1/z + Tier-1+2 ODD (finite B only)');
    plot(ax1, exp_pts(:,1),  exp_pts(:,2), 'kp', 'MarkerFaceColor', 'y', 'MarkerSize', 11, 'DisplayName', 'exp. anchors (R2007 Fig.1)');
    xlabel(ax1, 'T (K)');  ylabel(ax1, 'B_c (T)');
    title(ax1, 'LiHoF_4 phase boundary: ODD overlay (quick, B-cuts + closed-form endpoints)');
    legend(ax1, 'Location', 'northeast');  xlim(ax1, [0 2.4]);  ylim(ax1, [0 6]);
    text(ax1, 0.05, 0.05, {'Near-T_c(0) T-cut region left blank for ODD curves', ...
        '(invz\_critical\_T cannot bracket with ODD on; ODD-LOG T2/V4.1).'}, ...
        'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
    % Panel 2: critical self-energy Sigma(0) along the boundary.
    ax2 = subplot(2, 1, 2);  hold(ax2, 'on');  box(ax2, 'on');
    plot(ax2, [Tq Tc0_off],  [Sg_off  Sc0_off],  'o-', 'Color', [0 0.45 0.74],  'LineWidth', 1.4, 'DisplayName', '\Sigma(0), ODD off');
    plot(ax2, [Tq Tc0_full], [Sg_t1   Sc0_full], 's-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, 'DisplayName', '\Sigma(0), Tier-1 ODD');
    plot(ax2, Tq,            Sg_t12,             '^--', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.4, 'DisplayName', '\Sigma(0), Tier-1+2 ODD');
    xlabel(ax2, 'T (K)');  ylabel(ax2, '\Sigma(0) at the boundary');
    title(ax2, 'Critical self-energy along the boundary (B=0 endpoints: Richardson 12/24)');
    legend(ax2, 'Location', 'northwest');  xlim(ax2, [0 2.4]);

    figpath = fullfile(fileparts(mfilename('fullpath')), '..', 'Data', 'Phase_ODD_overlay_quick.fig');
    savefig(fig, figpath);
    fprintf('\n==== V4.1 ODD overlay (quick) ====\n');
    fprintf('  Tc(0): MF %.3f K | off %.5f K | Tier-1 ODD %.5f K  (dTc = %.4f K)\n', ...
            Tc0_mf, Tc0_off, Tc0_full, Tc0_off - Tc0_full);
    fprintf('  Sigma_c(0): off %.5f -> Tier-1 ODD %.5f  (dSc = %+.4f)\n', Sc0_off, Sc0_full, Sc0_full - Sc0_off);
    for it = 1:nT
        fprintf('  T = %.2f K: Bc off %.4f | T1 %.4f | T1+2 %.4f | MF %.4f T | Sigma0 off/T1/T12 %.4f/%.4f/%.4f\n', ...
                Tq(it), Bc_off(it), Bc_t1(it), Bc_t12(it), Bc_mf(it), Sg_off(it), Sg_t1(it), Sg_t12(it));
    end
    if t12_drop, fprintf('  NOTE: Tier-2 points dropped past the 15 min budget (see warning above).\n'); end
    fprintf('  Figure saved: %s\n', figpath);
    return;
end
% ===========================================================================
% Zero-field Tc, computed ONCE up front (invz_sigma_crit warns once here rather
% than in every worker): it anchors the Tc(B) adaptive window and is the B=0
% endpoint on the plot.
Sodd = struct();                            % ODD blocks (empty when ion.odd is off; parfor-safe)
if ion.odd
    % ---- ODD (T1.4): geometric blocks ONCE, pre-parfor, on the SAME 16^3 Gamma-less
    % grid invz_bz_couplings just used (its defaults: grid [16 16 16], dpRng 30, cache).
    Sodd = build_driver_odd_blocks(ion);
    % Odd-aware zero-field anchor (T1.5): invz_odd_zero_field REPLACES the old inline-handle
    % SEAM -- identical mode-'full' governing algebra (Sc(J0-d, modes(Vcc+dJ)), J0(T)=Jcc0-d(T)),
    % here on the SINGLE driver grid (16^3) so the adaptive-window anchor matches the parfor
    % solves' mesh (not the {12,24} Richardson production value). invz_critical_T refuses to
    % anchor adaptively with opts.odd on (invz:oddTc0), so Tc0 MUST be odd-aware. The engine
    % suppresses the known 16^3 invz:sigmaCritExcluded warning internally (ODD-LOG T1.3), so no
    % wrapper is needed here.
    Tc0 = invz_odd_zero_field(ion, struct('mode', 'full', 'grids', {{16}}, 'dpRng', 30, 'cache', true));
else
Tc0 = invz_critical_T0field(ion, invz_sigma_crit(J0, Jf), J0);
end

nT = numel(Ts);  nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];           % one independent 1-D root find per job
kind = [ones(1,nT) 2*ones(1,nB)];   % 1 = Bc(T) at fixed T, 2 = Tc(B) at fixed B

nWorkers = 0;                       % 0 = serial (also the no-toolbox fallback)
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

out = nan(1, nT + nB);
parfor (k = 1:nT+nB, nWorkers)
    t0 = tic;
    v = jobs(k);  val = NaN;
    if kind(k) == 1
        try
            if ion.odd
                % ODD (T1.4): Jnu_flat = [] (modes are rebuilt in the solver from the
                % blocks + deltaJ(T,B)); J0 stays UNSHIFTED (solver applies the E5 -d).
                val = invz_critical(ion, v, [], struct('J0eff', J0, 'Jxx0', Jxx0, ...
                    'window', [0.1 6], 'odd', true, 'odd_blocks', Sodd));
            else
            % Field window [0.1 6]: the top (6 T) is paramagnetic for every Ts
            % (Bc(T=0) ~ 5 T); invz_critical scans DOWN from there to the
            % converged ordered/paramagnet crossing (see its header).
            val = invz_critical(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'window', [0.1 6]));
            end
        catch err
            fprintf('  T = %.2f K: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    else
        try
            if ion.odd
                % ODD (T1.4): Tc0 here is the odd-aware anchor computed above
                % (invz_critical_T errors invz:oddTc0 without it).
                val = invz_critical_T(ion, v, [], struct('J0eff', J0, 'Jxx0', Jxx0, ...
                    'Tc0', Tc0, 'odd', true, 'odd_blocks', Sodd));
            else
            val = invz_critical_T(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'Tc0', Tc0));
            end
        catch err
            fprintf('  B = %.2f T: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] B = %.2f T  ->  Tc = %.3f K   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    end
    out(k) = val;
end
Bc  = out(1:nT);                    % low-T branch:  Bc at each Ts
TcB = out(nT+1:end);                % high-T branch: Tc at each Bs

fprintf('Zero-field Tc (1/z) = %.3f K  [target 1.74 K]\n', Tc0);   % computed up front

% Merged boundary, T-sorted, finite points only -- workspace export for
% downstream use ('boundary' would shadow the built-in of that name).
% Columns [T(K) B(Tesla)]; the plot uses tesla directly. Near the regime join
% (~1.5-1.6 K) both branches contribute a point, so the curve is not
% strictly single-valued in T there.
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);

figure; hold on;
plot(Ts, Bc, 'o-');
plot(TcB, Bs, 's-');
plot(Tc0, 0, 'ks');
xlabel('T (K)'); ylabel('B_c (T)');
title('LiHoF_4 1/z phase boundary (paramagnetic side)');
legend({'B_c(T): fixed-T field cut', 'T_c(B): fixed-B temperature cut', 'closed-form T_c(B=0)'}, 'Location', 'southwest');

% ---------------------------------------------------------------------------
% (The T1.4 ODD anchor helpers odd_d_at / odd_Sc_at were removed at T1.5: the
% odd-aware zero-field anchor Tc0 above is now produced by invz_odd_zero_field
% (mode 'full', single 16^3 grid), which owns that mode-'full' governing algebra
% and its own invz:sigmaCritExcluded suppression.)

% ===========================================================================
% Shared ODD-blocks builder (the overlay branch AND the standard ion.odd branch
% build the identical struct; both call this). 16^3 Gamma-excluded driver mesh
% (same generator/range/filter as invz_bz_couplings' default grid) ->
% invz_odd_blocks -> the solver-opts blocks struct (P0.4: built ONCE, outside
% every solver loop).
function Sodd = build_driver_odd_blocks(ion)
[qodd, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
qodd = qodd(any(abs(qodd) > 1e-12, 2), :);
[VcaO, VcbO, VccO, infoO] = invz_odd_blocks(ion, qodd, struct('dpRng', 30, 'cache', true));
Sodd = struct('Vca', VcaO, 'Vcb', VcbO, 'Vcc', VccO, 'Jcc0', infoO.Jcc0);   % Jcc0 UNSHIFTED
end

% ===========================================================================
% V4.1 overlay helpers (used only when overlay_quick is true; local functions of
% this script). Each is defensive -- a single finder failure returns NaN and the
% overlay still renders every element that DID resolve.
function bx = try_bc(fn)
% Run a boundary finder, returning NaN on any error (incl. invz:bracket).
try
    bx = fn();
catch
    bx = NaN;
end
end

function o = setwin(o, w)
o.window = w;
end

function s = sigma_at(ion, T, B, Jf, opts)
% Critical self-energy Sigma(0) at (T, B): re-solve and read pt.Sigma0.
if ~isfinite(B), s = NaN; return; end
try
    pt = invz_solve_point(ion, T, B, Jf, opts);
    s  = pt.Sigma0;
catch
    s = NaN;
end
end

function bx = mf_boundary(ion, T, Jf, opts, J0, Bhi, Bc1z)
% Mean-field (RPA) boundary Bc(T): root of crit_MF = 1 - J0*chi0cc0 (Sigma = 0)
% over field, using the SAME converged single-ion cc propagator (pt.chi0cc0) as
% the 1/z solve, so crit_MF == pt.crit - pt.Sigma0 exactly. The MF boundary sits
% ABOVE the 1/z boundary (it needs chi0cc0 = 1/J0 < (1+Sigma_c)/J0), so bracket
% just above the 1/z crossing up to the PM window top.
if ~isfinite(Bc1z), bx = NaN; return; end
g = @(B) mf_crit(ion, T, B, Jf, opts, J0);
try
    bx = fzero(g, [Bc1z + 0.05, Bhi]);
catch
    bx = NaN;
end
end

function c = mf_crit(ion, T, B, Jf, opts, J0)
pt = invz_solve_point(ion, T, B, Jf, opts);
if pt.converged && isfinite(pt.chi0cc0)
    c = 1 - J0*pt.chi0cc0;            % == pt.crit - pt.Sigma0
else
    c = NaN;
end
end
% ===========================================================================
