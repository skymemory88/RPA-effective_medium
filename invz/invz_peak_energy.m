function E = invz_peak_energy(chi, w, wmin)
%INVZ_PEAK_ENERGY Per-column single peak of chi''(w) at w >= wmin, parabolic-refined.
%   E = invz_peak_energy(chi, w, wmin) returns the CENSORED (NaN), parabolic-refined peak
%   frequency of each column of chi ([nw x nB], e.g. S.chiz/S.chirpa from invz_spectra_map
%   or invz_spectra_qpath), restricted to w >= wmin (e.g. to exclude a known low-frequency
%   line such as the hyperfine pole). A column is censored to NaN when its maximum is
%   non-positive/non-finite or sits in the first/last usable bin -- a boundary maximum means
%   the true peak lies outside the sampled window, so reporting the grid edge would
%   fabricate a flat dispersion. Assumes uniform w spacing.
%
%   Shared by invz_spectra_qpath (per q-point) and invz_spectra_map (per field).
E = nan(1, size(chi, 2));
mask = w >= wmin;
wm = w(mask);
if numel(wm) < 3
    warning('invz:peakWindowEmpty', ['invz_peak_energy: wmin = %.4g excludes (nearly) the ' ...
        'entire w grid (max(w) = %.4g) -- every column censored to NaN. Lower wmin or widen w.'], ...
        wmin, max(w));
    return;
end
dw = wm(2) - wm(1);
for k = 1:size(chi, 2)
    c = chi(mask, k);
    [cmax, i] = max(c);
    if ~isfinite(cmax) || cmax <= 0 || i == 1 || i == numel(c), continue; end
    d = (c(i-1) - c(i+1)) / (2*(c(i-1) - 2*c(i) + c(i+1)));   % parabolic vertex offset
    if ~isfinite(d) || abs(d) > 1, d = 0; end
    E(k) = wm(i) + d*dw;
end
end
