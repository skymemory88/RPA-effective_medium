%INVZ_RUN_SPECTRA chi''_cc(omega) at the uniform mode vs transverse field at fixed T.
%
% One driver, two views, chosen automatically from the number of field samples:
%   numel(fields) <= sliceMax  ->  1D line-slice overlay (1/z solid, RPA dashed, one
%                                  colour per field)                      cf. R 2007 Fig 2
%   numel(fields) >  sliceMax  ->  2D field-vs-frequency colormap, 1/z and RPA panels
%                                  (B on x, omega on y, log10 chi'' in colour)
%                                                            cf. R 2007 Fig 2 / Kovacevic Fig 3d
%
% The heavy lifting is in invz_spectra_map, which covers BOTH phases: at each field it uses
% the ferromagnetic (ordered) solve below Bc and the paramagnetic solve above, so a sweep
% across Bc shows the soft mode dip at the quantum phase transition (S.phase labels each
% field 1=FM/2=PM/0=masked). The result struct S replots for free, e.g.
% invz_plot_spectra_map(gca, S, S.chiz, '1/z').
% Save with:  print(gcf, 'spectra.png', '-dpng', '-r150');

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

ion = invz_ion();
% ion.demag = 1; ion.alpha = 0.5;    % optional sample-shape (demag) correction to J(0), applied in
%                                    % invz_jq_modes and read consistently by MF/RPA/1z. demag=0
%                                    % (default) = intrinsic c-needle coupling / R2007 benchmark;
%                                    % alpha = a/c aspect ratio (1 sphere, 0 c-needle, Inf disk).
T = 0.31;                             % K
% fields = [3.6 4.2 4.8 5.4 6.0];       % few -> slices;  many -> colormap
fields = linspace(3,6.5,151);
w = (0:0.002:0.45).';               % meV
eta = 1e-3;                          % real-axis line width: Lorentzian HWHM in meV (5e-3 meV ~ 1.2 GHz).
                                     % Lower -> sharper peaks (resolves the sub-6-GHz hyperfine lines),
                                     % but keep eta >~ the w-spacing above or the peaks alias.
sliceMax = 6;                         % field count at/below which the line-slice view is used
useParallel = true;                  % true -> parfor over fields (Parallel Computing Toolbox)
eUnit = 'meV';                       % 'meV' or 'GHz' -- plotting only; computation always runs in meV

C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end

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
