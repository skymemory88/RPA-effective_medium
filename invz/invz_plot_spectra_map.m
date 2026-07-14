function invz_plot_spectra_map(ax, S, chi, ttl, eUnit)
%INVZ_PLOT_SPECTRA_MAP Render one chi''(omega, B) colormap panel from invz_spectra_map output.
%   invz_plot_spectra_map(ax, S, chi, ttl) draws the map `chi` ([nw x nB], e.g. S.chiz or
%   S.chirpa) on axes `ax` against field S.fields (x) and frequency S.w (y).
%
%   invz_plot_spectra_map(ax, S, chi, ttl, eUnit) labels the y-axis with eUnit ('meV',
%   default, or 'GHz') instead. This only changes the label text -- S.w is plotted as-is,
%   so the caller is responsible for pre-scaling S.w to match eUnit (invz_run_spectra does
%   this via its eUnit knob).
%
%   The colour is log10(chi'') because a soft mode near criticality spans several decades;
%   the scale spans three decades below the (99.5th-percentile) peak so the dispersing mode
%   stays visible at every field.
%
%   Two greys, kept distinct on purpose:
%     - transparent (grey background) = NO paramagnetic solution (NaN): the ordered /
%       degenerate-doublet columns, so the phase boundary reads off as the mask edge;
%     - the darkest colour (floored) = present but negligible/negative chi''. On the real
%       axis the single-shot 1/z continuation can dip slightly negative just above the
%       mode (a causality artifact of npass=1, not missing data), so it is floored into
%       the map rather than punched out as a hole.
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
