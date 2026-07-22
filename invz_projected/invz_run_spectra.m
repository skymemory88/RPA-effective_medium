%INVZ_RUN_SPECTRA chi''_cc spectra vs transverse field (uniform mode) or along a q-path.
%
% One driver, three views, selected by qpath/fields/Bq below:
%   qpath = []  (default) -- field sweep at q = 0 (invz_spectra_map; S.phase labels each
%     field 1=FM/2=PM/0=masked):
%       numel(fields) <= sliceMax  ->  1D line-slice overlay          cf. R 2007 Fig 2
%       numel(fields) >  sliceMax  ->  2D field-vs-frequency colormap cf. R 2007 Fig 2 / Kovacevic Fig 3d
%     showPeaks = true also plots S.Epeak/S.Epeak_rpa vs field (invz_peak_energy).
%     theta_c (deg) tilts the swept field toward c to model misalignment; see the knob comment for validity.
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
T = 0.001;                             % K
useParallel = true;                  % true -> parfor over fields (Parallel Computing Toolbox)
eUnit = 'GHz';                       % 'meV' or 'GHz' -- unit for the frequency INPUTS (w, wq) AND
                                     % the plotted axes. Computation always runs in meV; the driver
                                     % converts in/out (with 'meV' it is a no-op). eta is ALWAYS in
                                     % meV (below), independent of eUnit.

% fields = [3.6 4.2 4.8 5.4 6.0];       % few -> slices;  many -> colormap
fields = linspace(0,9,301);
w = (0:0.005:6).';                  % eUnit -- field-sweep frequency grid (0-108 GHz ~ 0-0.45 meV)
eta = 2e-5;                          % real-axis Lorentzian HWHM, ALWAYS in meV (1e-3 meV ~ 0.24 GHz),
                                     % independent of eUnit. Lower -> sharper peaks (resolves the
                                     % sub-6-GHz hyperfine lines); keep eta above the w/wq step
                                     % (converted to meV: step/eScale) or the peaks alias.
sliceMax = 6;                         % field count at/below which the line-slice view is used
theta_c = 0.0;                         % deg -- tilt of the field OUT of the transverse ab-plane
                                     % toward c (the Ising axis): models experimental field
                                     % misalignment. The direction is FIXED across the sweep;
                                     % x-axes stay the total magnitude |B| (the longitudinal
                                     % component is |B|*sind(theta_c)). theta_c = 0 reproduces
                                     % the pure transverse benchmark exactly. Convention matches
                                     % LiReF4_MF_Yikai theta at phi = 0 (spec 2026-07-16).
                                     % SCALAR-STAGE VALIDITY: Sigma dresses the cc channel only;
                                     % exact at theta_c = 0, uncontrolled O(theta^2) cross-
                                     % channel error beyond the tensor-referenced small-tilt
                                     % range (invz_run_tensor_ref); a longitudinal component
                                     % turns the sharp transition into a rounded crossover.
                                     % Full-tensor propagation: deferred (invz_tensor/). phi_ab: implemented below.
phi_ab = -11.0;                        % deg -- IN-PLANE rotation of the swept field, a -> b.
                                     % phi_ab = -11 deg reproduces the production experimental
                                     % geometry (external stack ion.cfRot(Ho) = -11 deg; SAME
                                     % sign, pinned by test_invz_cfrot_equiv:
                                     % cfRot = -11 deg <=> phi_ab = -11 deg).
                                     % Nonzero phi_ab REQUIRES transverse_mf = 'vector_ab'
                                     % below ('none' also passes, as a bare CF+Zeeman
                                     % diagnostic; the library errors only under
                                     % 'legacy_x', by design).
                                     % NOTE: vector_ab shifts even phi_ab = 0 results
                                     % slightly (~0.04 ueV at 4 T; grows at low field) --
                                     % never compare legacy_x and vector_ab runs as if only
                                     % the angle differed. Combined theta_c AND phi_ab is
                                     % NOT validated (tilt bound was measured under legacy_x).
transverse_mf = 'vector_ab';         % 'legacy_x' | 'none' | 'vector_ab'

showPeaks = true;                    % true -> ALSO line-plot chi''_cc peak energy vs field
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

C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end

if phi_ab ~= 0 && strcmp(transverse_mf, 'legacy_x')
    error('invz_run_spectra:transverseMF', ...
        ['phi_ab = %.3g deg needs the vector transverse mean field: set the transverse_mf ' ...
         'knob above to ''vector_ab'' (or ''none'' for a bare CF+Zeeman diagnostic). ' ...
         'legacy_x is x-only and C4-inconsistent for rotated fields.'], phi_ab);
end
dhat = [cosd(theta_c)*cosd(phi_ab), cosd(theta_c)*sind(phi_ab), sind(theta_c)];  % unit field direction
tiltStr = '';
if theta_c ~= 0, tiltStr = sprintf(', \\theta_c = %.2g\\circ', theta_c); end
if phi_ab  ~= 0, tiltStr = [tiltStr sprintf(', \\phi_{ab} = %.2g\\circ', phi_ab)]; end
if ~strcmp(transverse_mf, 'legacy_x'), tiltStr = [tiltStr sprintf(', %s', transverse_mf)]; end

solve_opts = struct('transverse_mf', transverse_mf);   % merged into every spectra call below

% w and wq are given in eUnit; solves run in meV, so convert on the way in (eScale is the
% meV->eUnit factor) and scale returned grids (S.w, Epeak) back up by eScale for plotting.
% eta is ALWAYS in meV, independent of eUnit, and passes through unconverted.
wMeV   = w   / eScale;
wqMeV  = wq  / eScale;

if ~isempty(qpath)
    % ---------------- FM-mode q-path view at fixed field(s) ----------------
    if isscalar(Bq)
        S = invz_spectra_qpath(ion, T, Bq*dhat, qpath, wqMeV, struct('eta', eta, 'solve_opts', solve_opts));
        Splot = S;   % display-only copy; the solve above always ran in meV
        Splot.w = S.w*eScale;  Splot.Epeak = S.Epeak*eScale;  Splot.Epeak_rpa = S.Epeak_rpa*eScale;
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);
        invz_plot_spectra_qpath(ax1, Splot, Splot.chiz, Splot.Epeak, ...
            sprintf('1/z FM-mode \\chi''''_{cc}, T = %.2f K, B = %.2f T%s', T, Bq, tiltStr), eUnit);
        ax2 = subplot(1, 2, 2);
        invz_plot_spectra_qpath(ax2, Splot, Splot.chirpa, Splot.Epeak_rpa, ...
            sprintf('RPA FM-mode \\chi''''_{cc}, T = %.2f K, B = %.2f T%s', T, Bq, tiltStr), eUnit);
    else
        figure; hold on;  co = lines(numel(Bq));
        for k = 1:numel(Bq)
            Sk = invz_spectra_qpath(ion, T, Bq(k)*dhat, qpath, wqMeV, struct('eta', eta, 'solve_opts', solve_opts));
            plot(Sk.x, Sk.Epeak*eScale*dispScale,     '-',  'Color', co(k, :), ...
                 'DisplayName', sprintf('1/z, %.2f T', Bq(k)));
            plot(Sk.x, Sk.Epeak_rpa*eScale*dispScale, '--', 'Color', co(k, :), ...
                 'DisplayName', sprintf('RPA, %.2f T', Bq(k)));
        end
        xlabel(Sk.xlab);
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('FM-mode dispersion, T = %.2f K, dispScale = %.2f%s', T, dispScale, tiltStr));
        legend show;
        % cf. R 2007 Fig 3 TRENDS: the x-axis shows the actual varying Miller component
        % (h = 1..2 for the (1,0,0)->(2,0,0) path); their theory lines are the calculated
        % energies scaled by 1.15 (set dispScale = 1.15). Gaps in the lines are CENSORED
        % peaks (mode outside the wq window) -- widen wq, do not interpolate over them.
    end
else
    % ---------------- field-sweep views at the uniform mode ----------------
    S = invz_spectra_map(ion, T, fields, wMeV, ...
            struct('parallel', useParallel, 'eta', eta, 'field_dir', dhat, 'solve_opts', solve_opts));

    fprintf('Bc_auto (bare-MF dispatch; RPA proxy only where the ordered EMT converged) ~ %s T | Bc_1z (renormalized) ~ %s T  (sweep midpoints; masked/suspect columns widen the bracket)\n', ...
            num2str(S.Bc_auto), num2str(S.Bc_1z));
    if any(S.suspect)
        fprintf('WARNING: %d suspect column(s) (auto-PM with crit <= 0 -- spurious below-Bc PM points; masked): |B| = %s T\n', ...
                nnz(S.suspect), mat2str(S.fields(S.suspect), 3));
    end

    if numel(fields) <= sliceMax
        figure; hold on;  co = lines(numel(fields));
        for k = 1:numel(fields)
            plot(w, S.chiz(:, k),   '-',  'Color', co(k, :), 'DisplayName', sprintf('1/z, %.2f T', fields(k)));
            plot(w, S.chirpa(:, k), '--', 'Color', co(k, :), 'DisplayName', sprintf('RPA, %.2f T', fields(k)));
        end
        xlabel(eLabel);  ylabel('\chi''''_{cc}');  title(sprintf('T = %.2f K%s', T, tiltStr));  legend show;
    else
        Splot = S;  Splot.w = S.w * eScale;    % display-only copy; solve above always ran in meV
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);  invz_plot_spectra_map(ax1, Splot, Splot.chiz,   sprintf('1/z (own phase; Stage-2 ordered below B_c^{1/z}), T = %.2f K%s', T, tiltStr), eUnit);
        ax2 = subplot(1, 2, 2);  invz_plot_spectra_map(ax2, Splot, Splot.chirpa, sprintf('RPA, T = %.2f K%s', T, tiltStr), eUnit);
    end

    if showPeaks
        % ---- susceptibility peak energy vs field (toggle) --------------------------------
        figure; hold on;
        plot(S.fields, S.Epeak*eScale,     '-',  'DisplayName', '1/z');
        plot(S.fields, S.Epeak_rpa*eScale, '--', 'DisplayName', 'RPA');
        xlabel('|B| (T)');
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('\\chi''''_{cc} peak energy vs field, T = %.2f K%s', T, tiltStr));
        legend show;
        if isfinite(S.Bc_auto), xline(S.Bc_auto, '--', 'B_c^{auto}', 'HandleVisibility', 'off'); end
        if isfinite(S.Bc_1z),   xline(S.Bc_1z,   ':',  'B_c^{1/z}', 'HandleVisibility', 'off'); end
        % Gaps are CENSORED peaks (invz_peak_energy: boundary max or non-positive/non-finite
        % column) -- same convention as the q-path E_peak(q) stream, do not interpolate over them.
    end
end
