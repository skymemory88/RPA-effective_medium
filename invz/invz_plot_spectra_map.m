function invz_plot_spectra_map(ax, S, chi, ttl, eUnit)
%INVZ_PLOT_SPECTRA_MAP Render one chi''(omega, B) colormap panel from invz_spectra_map output.
%   invz_plot_spectra_map(ax, S, chi, ttl) draws the map `chi` ([nw x nB], e.g. S.chiz or
%   S.chirpa) on axes `ax` against field S.fields (x) and frequency S.w (y).
%
%   invz_plot_spectra_map(ax, S, chi, ttl, eUnit) labels the y-axis with eUnit ('meV',
%   default, or 'GHz') instead; this only changes the label text, so the caller must
%   pre-scale S.w to match (invz_run_spectra does this via its eUnit knob).
%
%   Colour is log10(chi'') spanning three decades below the 99.5th-percentile peak, so the
%   dispersing mode stays visible at every field. Two greys are kept distinct on purpose:
%   transparent = NO solution (NaN, phase boundary reads off as the mask edge); darkest
%   colour (floored) = present but negligible/negative chi'' (the single-shot real-axis
%   continuation can dip slightly negative just above the mode -- a causality artifact,
%   not missing data).
if nargin < 4, ttl = ''; end
if nargin < 5, eUnit = 'meV'; end

finiteMask = isfinite(chi);            % transparent only where there is no solution at all
Z = log10(max(chi, realmin));          % negatives -> log10(realmin), clamped to the floor below
im = imagesc(ax, S.fields, S.w, Z);
set(im, 'AlphaData', double(finiteMask));
set(ax, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'Layer', 'top');

pos = chi(finiteMask & chi > 0);
if ~isempty(pos)
    hi = robust_pct(pos, 0.995);
    lo = hi / 1e3;
    clim(ax, [log10(lo) log10(hi)]);
end
colormap(ax, turbo);
xlabel(ax, 'B_x (T)');   ylabel(ax, sprintf('\\omega (%s)', eUnit));   title(ax, ttl);
cb = colorbar(ax);   cb.Label.String = 'log_{10} \chi''''_{cc}';
end

function v = robust_pct(x, p)
%ROBUST_PCT p-quantile of x (p in [0,1]) without the Statistics Toolbox.
x = sort(x(:));
v = x(max(1, min(numel(x), ceil(p * numel(x)))));
end
