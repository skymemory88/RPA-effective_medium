function v = invz_robust_pct(x, p)
%INVZ_ROBUST_PCT p-quantile of x (p in [0,1]) without the Statistics Toolbox.
% Shared by invz_plot_spectra_map / invz_plot_spectra_qpath for the colour-limit clip.
x = sort(x(:));
v = x(max(1, min(numel(x), ceil(p * numel(x)))));
end
