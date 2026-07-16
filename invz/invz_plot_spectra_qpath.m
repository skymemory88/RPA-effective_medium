function invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit)
%INVZ_PLOT_SPECTRA_QPATH Render one branch chi''(q, omega) colormap panel (exploratory).
%   invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit) draws `chi` ([nw x nq], e.g.
%   S.chiz or S.chirpa) against the plot coordinate S.x (varying Miller component for a
%   single-axis path, e.g. h = 1..2 on the R2007 Fig-3 path, else S.s) and frequency S.w
%   (y), overlaying the censored peak dispersion `Epeak` (S.Epeak or S.Epeak_rpa) in white;
%   NaN peak points leave gaps. This is a BRANCH susceptibility, not neutron intensity.
%   Colour conventions match invz_plot_spectra_map (log10 scale, NaN transparent, negative
%   chi'' floored). The caller pre-scales S.w and Epeak when plotting in GHz.
if nargin < 5, ttl = ''; end
if nargin < 6, eUnit = 'meV'; end

if isfield(S, 'x')                       % S structs saved before the S.x field existed
    x = S.x;  xlab = S.xlab;
else
    x = S.s;  xlab = sprintf('s along path from Q = [%g %g %g] (index r.l.u.)', S.qpath(1,:));
end

finiteMask = isfinite(chi);
Z = log10(max(chi, realmin));
im = imagesc(ax, x, S.w, Z);
set(im, 'AlphaData', double(finiteMask));
set(ax, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'Layer', 'top');

pos = chi(finiteMask & chi > 0);
if ~isempty(pos)
    pos = sort(pos(:));
    hi = pos(max(1, min(numel(pos), ceil(0.995*numel(pos)))));
    clim(ax, [log10(hi/1e3) log10(hi)]);
end
colormap(ax, turbo);
hold(ax, 'on');  plot(ax, x, Epeak, 'w.-', 'MarkerSize', 8);  hold(ax, 'off');
xlabel(ax, xlab);
ylabel(ax, sprintf('\\omega (%s)', eUnit));
title(ax, ttl);
cb = colorbar(ax);  cb.Label.String = 'log_{10} \chi''''_{cc}';
end
