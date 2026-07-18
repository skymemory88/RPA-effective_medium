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
%
% ERROR POLICY: the sweep absorbs ONLY per-point invzt:bracket (a genuine
% no-crossing outcome once Brange itself has been preflighted below); every
% other identifier that ESCAPES invzt_critical rethrows. invzt_critical's
% own internal sampler additionally classifies shared-engine physics signals
% (and other non-invzt errors) as ordered votes -- its committed T7 policy,
% documented in its header, not overridden here.

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
                                       % sigma_floor may be added here too; defaults to
                                       % invzt_critical's own -0.5 (single-sourced below)

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

% Brange preflight (invzt:bracket doubles as invzt_critical's arg-validation
% id; absorbing it per point without this check would turn a reversed/typo'd
% bracket into a silent all-NaN sweep).
assert(isnumeric(Brange) && isreal(Brange) && numel(Brange) == 2 && ...
    all(isfinite(Brange)) && Brange(2) > Brange(1) && Brange(1) > 0, ...
    'invzt_run_phase_diagram:Brange', ...
    'Brange must be finite [Blo Bhi] with 0 < Blo < Bhi (mode ''a1'' forbids B = 0); got %s.', ...
    mat2str(Brange));

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
        % invzt:bracket here is a genuine per-point outcome (Brange was
        % preflighted above, so this window simply lacks a crossing at this
        % T -- expected near Tc0, where the boundary is near-vertical).
        % Anything else is a MISCONFIGURATION that would otherwise silently
        % produce a whole sweep of NaNs, so it rethrows.
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
    % Same NOMINAL N and dipole cutoff as this driver's lat -- but NB the
    % projected engine always builds its own legacy endpoint-inclusive
    % [-0.5, 0.5] mesh (a different quadrature convention from 'halfopen'),
    % on a different model with a different ODD treatment: a cross-model
    % COMPARATOR, never the same quantity.
    Tc0_closed = invz_odd_zero_field(ion, struct('mode', anchor_mode, ...
        'grids', {{gridN}}, 'dpRng', dpRng, 'cache', true));
end

figure; hold on;
plot(Bc, Ts, 'o-', 'DisplayName', 'tensor A1: field-cut B_c(T)');
if isfinite(Tc0_proxy)
    plot(0, Tc0_proxy, 'k^', 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('tensor small-B_x proxy T_c(0), %.2f T', Bproxy));
end
if isfinite(Tc0_closed)
    plot(0, Tc0_closed, 'ks', 'MarkerFaceColor', 'c', ...
         'DisplayName', 'projected closed-form T_c(0) (comparator; legacy-inclusive mesh)');
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
    % Same VALIDITY floor as invzt_critical (single-sourced sigma_floor,
    % default -0.5). NB deliberately NO crit>0 term here, unlike the spectra
    % driver's otherwise-similar check: ok is validity-only by
    % invzt_tc_pm_extrap's contract -- IT applies the crit>0 PM filter, so
    % metastable ordered-side samples are FILTERED there, never asserted on.
    ok = pt.converged && isfinite(c) && pt.Sigma0 >= getf(opts, 'sigma_floor', -0.5);
catch err
    switch err.identifier
        case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
            c = NaN;  ok = false;      % phase/physics signal, not an error
        otherwise
            rethrow(err);              % invzt:mode, invzt:a1OddGamma, ... = misconfiguration
    end
end
end
