%INVZ_RUN_SPECTRA chi''_cc spectra vs transverse field (uniform mode) or along a q-path.
%
% One driver, three views:
%   qpath = []  (default)  -- field sweep at the uniform mode (q = 0):
%     numel(fields) <= sliceMax  ->  1D line-slice overlay (1/z solid, RPA dashed, one
%                                    colour per field)                      cf. R 2007 Fig 2
%     numel(fields) >  sliceMax  ->  2D field-vs-frequency colormap, 1/z and RPA panels
%                                                            cf. R 2007 Fig 2 / Kovacevic Fig 3d
%   qpath = [nq x 3] r.l.u.  -- EXPLORATORY branch-resolved q-path view at fixed field(s)
%     Bq, for comparison with the TRENDS in R 2007 Fig 3 (branch susceptibility, not
%     neutron intensity; see invz_spectra_qpath header for the inherited caveats):
%     numel(Bq) == 1  ->  2D path-vs-frequency colormaps (1/z + RPA), censored peak overlay
%     numel(Bq) >  1  ->  E_peak(q) dispersion overlay, one colour per field
%
% The field sweep is done by invz_spectra_map (both phases: FM below Bc, PM above; S.phase
% labels each field 1=FM/2=PM/0=masked); the q-path by invz_spectra_qpath (one 1/z medium
% solve per field, then the pole formula along the guarded path couplings). Both result
% structs replot for free, e.g. invz_plot_spectra_map(gca, S, S.chiz, '1/z').
% Save with:  print(gcf, 'spectra.png', '-dpng', '-r150');

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

ion = invz_ion();
% ion.demag = 1; ion.alpha = 0.5;    % OPTIONAL sample-shape (demagnetization) knob; default off.
%   demag = 0 (default): intrinsic response -- the R 2007 benchmark. demag ~= 0 with spheroid
%   aspect ratio alpha (1 sphere, 0 c-needle, Inf disk) adds the sample shape as follows:
%     - info.Jcc0, Jnu, and the ordering-channel contribution to criticality are
%       demag-INVARIANT (R 2007: the demagnetizing field cancels from the critical
%       condition; ordering occurs at q -> 0+, not strict q = 0);
%     - Tc(B = 0) is EXACTLY demag-invariant (the transverse moment vanishes there);
%     - Bc(T) vs APPLIED transverse field CAN shift through the demag-aware transverse
%       coupling info.Jaa0 (internal-vs-applied field relation);
%     - the strict-uniform (q = 0) observable chi''_cc is demag-corrected via info.Jshape_cc
%       (chi_meas = chi/(1 + Jshape*chi)): the soft mode saturates instead of diverging;
%     - q-path spectra omit the Jshape_cc transform (finite-q probe = intrinsic longitudinal
%       response) but still see demag through info.Jaa0.
T = 0.31;                             % K
% fields = [3.6 4.2 4.8 5.4 6.0];       % few -> slices;  many -> colormap
fields = linspace(3,6.5,151);
w = (0:0.002:0.45).';               % meV -- field-sweep views
eta = 1e-3;                          % real-axis line width: Lorentzian HWHM in meV (5e-3 meV ~ 1.2 GHz).
                                     % Lower -> sharper peaks (resolves the sub-6-GHz hyperfine lines),
                                     % but keep eta >~ the w-spacing above or the peaks alias.
sliceMax = 6;                         % field count at/below which the line-slice view is used
useParallel = true;                  % true -> parfor over fields (Parallel Computing Toolbox)
eUnit = 'meV';                       % 'meV' or 'GHz' -- plotting only; computation always runs in meV

% ---- q-path view (R 2007 Fig 3 trends): set qpath non-empty to switch views -------------
% qpath = [];                          % [] = field-sweep views; [nq x 3] r.l.u. = q-path view

qh = linspace(1, 2, 51).';  
qpath = [qh zeros(numel(qh), 2)];  % (1,0,0)->(2,0,0)

Bq = 4.24;                           % field(s), T, for the q-path view. One value -> colormaps;
% Bq = [3.6 4.24 6.0];           % several -> E_peak(q) overlay (R 2007 Fig 3: [3.6 4.24 6.0])
wq = (0:0.004:0.85).';               % meV -- q-path grid. Fig 3 reaches ~0.75 meV near h = 1 at
                                     % 60 kOe (after their 1.15 scaling); 0.85 avoids clipping,
                                     % which the censoring peak picker would flag as NaN.
dispScale = 1;                       % dispersion display scale factor; R 2007 scales the
                                     % calculated energies by 1.15 to match experiment (Fig 3)

C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end

if ~isempty(qpath)
    % ---------------- exploratory q-path view at fixed field(s) ----------------
    if isscalar(Bq)
        S = invz_spectra_qpath(ion, T, Bq, qpath, wq, struct('eta', eta));
        Splot = S;   % display-only copy; the solve above always ran in meV
        Splot.w = S.w*eScale;  Splot.Epeak = S.Epeak*eScale;  Splot.Epeak_rpa = S.Epeak_rpa*eScale;
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);
        invz_plot_spectra_qpath(ax1, Splot, Splot.chiz, Splot.Epeak, ...
            sprintf('1/z branch \\chi''''_{cc} (exploratory), T = %.2f K, B = %.2f T', T, Bq), eUnit);
        ax2 = subplot(1, 2, 2);
        invz_plot_spectra_qpath(ax2, Splot, Splot.chirpa, Splot.Epeak_rpa, ...
            sprintf('RPA branch \\chi''''_{cc} (exploratory), T = %.2f K, B = %.2f T', T, Bq), eUnit);
    else
        figure; hold on;  co = lines(numel(Bq));
        for k = 1:numel(Bq)
            Sk = invz_spectra_qpath(ion, T, Bq(k), qpath, wq, struct('eta', eta));
            plot(Sk.x, Sk.Epeak*eScale*dispScale,     '-',  'Color', co(k, :), ...
                 'DisplayName', sprintf('1/z, %.2f T', Bq(k)));
            plot(Sk.x, Sk.Epeak_rpa*eScale*dispScale, '--', 'Color', co(k, :), ...
                 'DisplayName', sprintf('RPA, %.2f T', Bq(k)));
        end
        xlabel(Sk.xlab);
        ylabel(strrep(eLabel, '\omega', 'E_{peak}'));
        title(sprintf('branch dispersion (exploratory), T = %.2f K, dispScale = %.2f', T, dispScale));
        legend show;
        % cf. R 2007 Fig 3 TRENDS: the x-axis shows the actual varying Miller component
        % (h = 1..2 for the (1,0,0)->(2,0,0) path); their theory lines are the calculated
        % energies scaled by 1.15 (set dispScale = 1.15). Gaps in the lines are CENSORED
        % peaks (mode outside the wq window) -- widen wq, do not interpolate over them.
    end
else
    % ---------------- field-sweep views at the uniform mode ----------------
    S = invz_spectra_map(ion, T, fields, w, struct('parallel', useParallel, 'eta', eta));

    if numel(fields) <= sliceMax
        figure; hold on;  co = lines(numel(fields));
        for k = 1:numel(fields)
            plot(w*eScale, S.chiz(:, k),   '-',  'Color', co(k, :), 'DisplayName', sprintf('1/z, %.2f T', fields(k)));
            plot(w*eScale, S.chirpa(:, k), '--', 'Color', co(k, :), 'DisplayName', sprintf('RPA, %.2f T', fields(k)));
        end
        xlabel(eLabel);  ylabel('\chi''''_{cc}');  title(sprintf('T = %.2f K', T));  legend show;
    else
        Splot = S;  Splot.w = S.w * eScale;    % display-only copy; solve above always ran in meV
        figure('Position', [100 100 1150 460]);
        ax1 = subplot(1, 2, 1);  invz_plot_spectra_map(ax1, Splot, Splot.chiz,   sprintf('1/z, T = %.2f K', T), eUnit);
        ax2 = subplot(1, 2, 2);  invz_plot_spectra_map(ax2, Splot, Splot.chirpa, sprintf('RPA, T = %.2f K', T), eUnit);
    end
end
