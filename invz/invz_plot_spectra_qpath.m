function invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit)
%INVZ_PLOT_SPECTRA_QPATH Render one branch chi''(q, omega) colormap panel (exploratory).
%   invz_plot_spectra_qpath(ax, S, chi, Epeak, ttl, eUnit) draws `chi` ([nw x nq], e.g.
%   S.chiz or S.chirpa) against the INDEX-coordinate path distance S.s (x; r.l.u. Miller
%   coordinates, not Cartesian reciprocal distance -- S.s_cart carries that) and frequency
%   S.w (y), and overlays the censored peak dispersion `Epeak` (pass S.Epeak or
%   S.Epeak_rpa) in white; censored (NaN) peak points simply leave gaps. This is a BRANCH
%   susceptibility, not neutron intensity (no structure-factor/form-factor weights).
%   Colour conventions match invz_plot_spectra_map: log10 scale spanning three decades
%   below the robust (99.5th-percentile) peak; NaN transparent on grey; present-but-
%   negative chi'' floored into the darkest colour. As there, the caller pre-scales S.w
%   and Epeak when plotting in GHz (see invz_run_spectra's eUnit knob).
if nargin < 5, ttl = ''; end
if nargin < 6, eUnit = 'meV'; end

finiteMask = isfinite(chi);
Z = log10(max(chi, realmin));
im = imagesc(ax, S.s, S.w, Z);
set(im, 'AlphaData', double(finiteMask));
set(ax, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'Layer', 'top');

pos = chi(finiteMask & chi > 0);
if ~isempty(pos)
    x = sort(pos(:));
    hi = x(max(1, min(numel(x), ceil(0.995*numel(x)))));
    clim(ax, [log10(hi/1e3) log10(hi)]);
end
colormap(ax, turbo);
hold(ax, 'on');  plot(ax, S.s, Epeak, 'w.-', 'MarkerSize', 8);  hold(ax, 'off');
xlabel(ax, sprintf('s along path from Q = [%g %g %g] (index r.l.u.)', S.qpath(1,:)));
ylabel(ax, sprintf('\\omega (%s)', eUnit));
title(ax, ttl);
cb = colorbar(ax);  cb.Label.String = 'log_{10} \chi''''_{cc}';
end
