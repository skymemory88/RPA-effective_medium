%INVZ_RUN_SPECTRA chi''_cc spectra vs transverse field (uniform mode) or along a q-path.
%
% One driver, three views, selected by qpath/fields/Bq below:
%   qpath = []  (default) -- field sweep at q = 0 (invz_spectra_map; S.phase labels each
%     field 1=FM/2=PM/0=masked):
%       numel(fields) <= sliceMax  ->  1D line-slice overlay          cf. R 2007 Fig 2
%       numel(fields) >  sliceMax  ->  2D field-vs-frequency colormap cf. R 2007 Fig 2 / Kovacevic Fig 3d
%     showPeaks = true also plots S.Epeak/S.Epeak_rpa vs field (invz_peak_energy).
%   qpath = [nq x 3] r.l.u. -- FM-mode q-path view at fixed field(s) Bq (invz_spectra_qpath;
%     see its header for caveats), reproducing R 2007 Fig 3 TRENDS:
%       numel(Bq) == 1  ->  2D path-vs-frequency colormaps with censored peak overlay
%       numel(Bq) >  1  ->  E_peak(q) dispersion overlay, one colour per field
%
% Both result structs replot for free, e.g. invz_plot_spectra_map(gca, S, S.chiz, '1/z').
% Save with:  print(gcf, 'spectra.png', '-dpng', '-r150');

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));  addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();
% ion.demag = 1; ion.alpha = 0.5;    % OPTIONAL sample-shape (demagnetization) knob; default off.
%   demag = 0 (default): intrinsic response, the R 2007 benchmark. demag ~= 0 with spheroid
%   aspect ratio alpha (1 sphere, 0 c-needle, Inf disk) adds the sample shape:
%     - info.Jcc0, Jnu, and the ordering-channel criticality are demag-INVARIANT (R 2007:
%       the demagnetizing field cancels from the critical condition; ordering occurs at
%       q -> 0+, not strict q = 0); Tc(B = 0) is exactly demag-invariant too.
%     - Bc(T) vs APPLIED transverse field can shift through demag-aware info.Jaa0
%       (internal-vs-applied field relation).
%     - the strict-uniform (q = 0) chi''_cc is demag-corrected via info.Jshape_cc (saturates
%       instead of diverging); q-path spectra omit that transform (finite-q = intrinsic
%       response) but still see demag through info.Jaa0.
T = 0.1;                             % K
useParallel = true;
outerMix = 0.3;                     % smaller than 0.7: stronger damping;
outerMax = 2000;                     % exploratory budget; acceptance tolerances are unchanged
useMissingAreaApproximation = true; % OPT-IN production approximation; false preserves strict default
% missingAreaFactors = [0.75 1.0 1.5];   % continuous positive lower completions; factor 1 is central
missingAreaFactors = 1.0;            % factor 1: approximating constant r(x) ( =f(h_e) )
missingAreaNodes = 129;              % user-tunable profile resolution; 257 remains the refinement comparator
useAdjacentFieldRetry = true;        % retry cold-pass masks only when both ordered neighbours agree
adjacentRetryMaxSpan = 0.35;         % T; maximum separation of the two accepted seed fields
useOrderedBoundaryRetry = true;      % retry the final ordered-side mask from two untouched lower sources
orderedBoundaryRetryMaxSpan = 0.20;  % T; oldest of the two lower sources must lie within this span
orderedBoundaryMinDuni = 1e-3;       % explicit distance from the uniform static boundary
orderedBoundaryMinFprime = 1e-3;     % explicit H_MF-root slope margin
eUnit = 'GHz';                       % 'meV' or 'GHz' -- unit for the frequency INPUTS (w, wq) AND
                                     % the plotted axes. Computation always runs in meV; the driver
                                     % converts in/out (with 'meV' it is a no-op). eta is ALWAYS in
                                     % meV (below), independent of eUnit.

% Exploratory broad view containing the projected 1/z transition.  The retained
% finite-16^3 regression is the narrower linspace(4.60,4.90,61) subset; this
% 3--6 T view is expected to contain masked ordered columns and is not all-column certified.
fields = linspace(0, 9, 101);
% fields = 0.2:0.2:3.0;                % test points
% fields = [3.6 4.2 4.8 5.4 6.0];    % few -> line slices instead of a colormap
w = (0:0.005:6).';                    % eUnit -- field-sweep frequency grid
eta = 5e-5;                          % real-axis Lorentzian HWHM, ALWAYS in meV (1e-3 meV ~ 0.24 GHz),
                                     % independent of eUnit. Lower -> sharper peaks (resolves the
                                     % sub-6-GHz hyperfine lines); keep eta above the w/wq step
                                     % (converted to meV: step/eScale) or the peaks alias.
sliceMax = 6;                         % field count at/below which the line-slice view is used
showPeaks = true;                     % true -> ALSO line-plot chi''_cc peak energy vs field
                                     % (S.Epeak/S.Epeak_rpa, cf. the q-path E_peak(q) stream)

% ---- q-path view (R 2007 Fig 3 trends): set qpath non-empty to switch views -------------
qpath = [];                          % [] = field-sweep views; [nq x 3] r.l.u. = q-path view

% qh = linspace(0, 1, 101).';  
% qpath = [qh zeros(numel(qh), 2)];  % (h0,0,0)->(h1,0,0)
% qpath = [zeros(numel(qh), 1) qh zeros(numel(qh), 1)];  % (0,h0,0)->(0,h1,0)
% qpath = [zeros(numel(qh), 2) qh];  % (0,0,h0)->(0,0,h1) problematic

Bq = 4.95;                           % field(s), T, for the q-path view. One value -> colormaps;
% Bq = [4.75 4.85 4.95 5.05];
% Bq = [3.6 4.24 6.0];           % several -> E_peak(q) overlay (R 2007 Fig 3: [3.6 4.24 6.0])
wq = (0:0.02:6).';                 % eUnit -- q-path frequency grid (0-6 GHz ~ 0-0.025 meV). The
                                     % Fig-3 mode reaches ~0.75 meV (~180 GHz) at 60 kOe; keep the
                                     % top above the mode or the censoring peak picker NaNs it.
dispScale = 1;                       % dispersion display scale factor; R 2007 scales the
                                     % calculated energies by 1.15 to match experiment (Fig 3)

% ---- coupling ---------------------------------------------------------------
gridN = 16;                           % debug resolution; use a grid ladder for reported boundaries
dipoleBackend = 'ewald';             % 'ewald' (production) | 'bruteforce' (debug reference)
dpRng = 30;                          % brute-force cutoff; ignored by Ewald
ewaldOpts = struct('alpha', 0.3, 'r_cut', 16, 'g_cut', 3, ...
    'boundary', 'conducting_k0_omitted');

C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end

% Strict full_profile is the default.  The opt-in ensemble replaces only the
% unresolved lower area A=int_0^h_e r dh by declared positive completions and
% integrates the contiguous certified component above h_e.  It does not use a
% PM anchor, bridge rejected nodes, or claim thermodynamic branch selection.
solve_opts = struct('mix_outer', outerMix, 'max_outer', outerMax);
if useMissingAreaApproximation
    solve_opts.hmf_integral_mode = 'missing_area_ensemble';
    solve_opts.hmf_missing_area_factors = missingAreaFactors;
    solve_opts.nH = missingAreaNodes;
    solve_opts.hmf_approx_branch = ...
        'picard_attracting_contiguous_high_h_component';
    solve_opts.hmf_adjacent_retry = useAdjacentFieldRetry;
    solve_opts.hmf_adjacent_retry_max_span = adjacentRetryMaxSpan;
    solve_opts.hmf_ordered_boundary_retry = useOrderedBoundaryRetry;
    solve_opts.hmf_ordered_boundary_retry_max_span = ...
        orderedBoundaryRetryMaxSpan;
    solve_opts.hmf_ordered_boundary_retry_min_D_uni = ...
        orderedBoundaryMinDuni;
    solve_opts.hmf_ordered_boundary_retry_min_Fprime = ...
        orderedBoundaryMinFprime;
    warning('invz:missingAreaApproximation', ...
        ['Using the opt-in missing-area ensemble. Outputs are approximation ' ...
         'sensitivity bands, not certified thermodynamic branch selection.']);
end

% The quadrature convention is fixed; backend and resolution remain explicit.
coupling_opts = struct('grid', [gridN gridN gridN], 'dpRng', dpRng, ...
    'cache', false, 'dipole', dipoleBackend, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_drop');
if strcmp(dipoleBackend, 'ewald'), coupling_opts.ewald = ewaldOpts; end
qpOpts  = coupling_opts;  qpOpts.eta = eta;  qpOpts.solve_opts = solve_opts;
mapOpts = coupling_opts;  mapOpts.parallel = useParallel;  mapOpts.eta = eta;
mapOpts.solve_opts = solve_opts;

% w and wq are given in eUnit; solves run in meV, so convert on the way in (eScale is the
% meV->eUnit factor) and scale returned grids (S.w, Epeak) back up by eScale for plotting.
% eta is ALWAYS in meV, independent of eUnit, and passes through unconverted.
wMeV   = w   / eScale;
wqMeV  = wq  / eScale;

if ~isempty(qpath)
    if useMissingAreaApproximation
        error('invz:missingAreaQPath', ...
            ['The ensemble production mode is implemented for the field map. ' ...
             'Set useMissingAreaApproximation=false for q-path calculations.']);
    end
    % ---------------- FM-mode q-path view at fixed field(s) ----------------
    if isscalar(Bq) %#ok<UNRCH>
        S = invz_spectra_qpath(ion, T, Bq, qpath, wqMeV, qpOpts);
        Splot = S;   % display-only copy; the solve above always ran in meV
        Splot.w = S.w*eScale;  Splot.Epeak = S.Epeak*eScale;  Splot.Epeak_rpa = S.Epeak_rpa*eScale;
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);
        invz_plot_spectra_qpath(ax1, Splot, Splot.chiz, Splot.Epeak, ...
            sprintf('1/z FM-mode \\chi''''_{cc}, T = %.2f K, B = %.2f T', T, Bq), eUnit);
        ax2 = subplot(1, 2, 2);
        invz_plot_spectra_qpath(ax2, Splot, Splot.chirpa, Splot.Epeak_rpa, ...
            sprintf('RPA FM-mode \\chi''''_{cc}, T = %.2f K, B = %.2f T', T, Bq), eUnit);
    else
        figure; hold on;  co = lines(numel(Bq));
        for k = 1:numel(Bq)
            Sk = invz_spectra_qpath(ion, T, Bq(k), qpath, wqMeV, qpOpts);
            plot(Sk.x, Sk.Epeak*eScale*dispScale,     '-',  'Color', co(k, :), ...
                 'DisplayName', sprintf('1/z, %.2f T', Bq(k)));
            plot(Sk.x, Sk.Epeak_rpa*eScale*dispScale, '--', 'Color', co(k, :), ...
                 'DisplayName', sprintf('RPA, %.2f T', Bq(k)));
        end
        xlabel(Sk.xlab);
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('FM-mode dispersion, T = %.2f K, dispScale = %.2f', T, dispScale));
        legend show;
        % cf. R 2007 Fig 3 TRENDS: the x-axis shows the actual varying Miller component
        % (h = 1..2 for the (1,0,0)->(2,0,0) path); their theory lines are the calculated
        % energies scaled by 1.15 (set dispScale = 1.15). Gaps in the lines are CENSORED
        % peaks (mode outside the wq window) -- widen wq, do not interpolate over them.
    end
else
    % ---------------- field-sweep views at the uniform mode ----------------
    S = invz_spectra_map(ion, T, fields, wMeV, mapOpts);

    fprintf(['Bc_1z (sweep midpoint; use invz_critical for refinement) ~ %s T | ' ...
             'Bc_RPA (independent bare-MF/RPA dispatch) ~ %s T\n'], ...
            num2str(S.Bc_1z), num2str(S.Bc_rpa));

    if numel(fields) <= sliceMax
        figure; hold on;  co = lines(numel(fields));
        for k = 1:numel(fields)
            plot(w, S.chiz(:, k),   '-',  'Color', co(k, :), 'DisplayName', sprintf('1/z, %.2f T', fields(k)));
            plot(w, S.chirpa(:, k), '--', 'Color', co(k, :), 'DisplayName', sprintf('RPA, %.2f T', fields(k)));
        end
        xlabel(eLabel);  ylabel('\chi''''_{cc}');  title(sprintf('T = %.2f K', T));  legend show;
    else
        Splot = S;  Splot.w = S.w * eScale;    % display-only copy; solve above always ran in meV
        if useMissingAreaApproximation
            onezTitle = sprintf(['1/z MISSING-AREA ENSEMBLE central member ' ...
                '(bands in S.approximation_ensemble), T = %.2f K'],T);
        else
            onezTitle = sprintf('1/z STRICT 33-node Eq.-45 profile, T = %.2f K',T); %#ok<UNRCH>
        end
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);  invz_plot_spectra_map(ax1, Splot, Splot.chiz, onezTitle, eUnit);
        ax2 = subplot(1, 2, 2);  invz_plot_spectra_map(ax2, Splot, Splot.chirpa, sprintf('RPA, T = %.2f K', T), eUnit);
    end

    if showPeaks
        % ---- susceptibility peak energy vs field (toggle) --------------------------------
        figure; hold on; %#ok<UNRCH>
        plot(S.fields, S.Epeak*eScale,     '-',  'DisplayName', '1/z');
        plot(S.fields, S.Epeak_rpa*eScale, '--', 'DisplayName', 'RPA');
        xlabel('|B| (T)');
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('\\chi''''_{cc} peak energy vs field, T = %.2f K', T));
        legend show;
        if isfinite(S.Bc_1z), xline(S.Bc_1z, ':', 'B_c^{1/z}', 'HandleVisibility', 'off'); end
        if isfinite(S.Bc_rpa), xline(S.Bc_rpa, '--', 'B_c^{RPA}', 'HandleVisibility', 'off'); end
        % Gaps are CENSORED peaks (invz_peak_energy: boundary max or non-positive/non-finite
        % column) -- same convention as the q-path E_peak(q) stream, do not interpolate over them.
    end
end
