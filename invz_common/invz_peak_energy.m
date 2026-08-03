function E = invz_peak_energy(chi, w, wmin, mode)
%INVZ_PEAK_ENERGY Per-column single peak of chi''(w) at w >= wmin, parabolic-refined.
%   E = invz_peak_energy(chi, w, wmin) returns the CENSORED (NaN), parabolic-refined peak
%   frequency of each column of chi ([nw x nB], e.g. S.chiz/S.chirpa from invz_spectra_map
%   or invz_spectra_qpath), restricted to w >= wmin (e.g. to exclude a known low-frequency
%   line such as the hyperfine pole). A column is censored to NaN when its maximum is
%   non-positive/non-finite or sits in the first/last usable bin -- a boundary maximum means
%   the true peak lies outside the sampled window, so reporting the grid edge would
%   fabricate a flat dispersion. Assumes uniform w spacing.
%
%   E = invz_peak_energy(...,mode) uses mode='strongest' (default, preserving
%   the established electro-nuclear route) or mode='lowest_local', which
%   follows the lowest positive interior local maximum for an electronic
%   soft-mode plot when hyperfine is disabled.
%
%   Shared by invz_spectra_qpath (per q-point) and invz_spectra_map (per field).
if nargin < 4, mode = 'strongest'; end
mode = char(mode);
if ~ismember(mode, {'strongest', 'lowest_local'})
    error('invz:peakMode', ...
        'Peak mode must be ''strongest'' or ''lowest_local''.');
end
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
    if strcmp(mode, 'strongest')
        [cmax, i] = max(c);
    else
        candidates = find(isfinite(c(2:end-1)) & c(2:end-1) > 0 & ...
            c(2:end-1) >= c(1:end-2) & c(2:end-1) > c(3:end)) + 1;
        if isempty(candidates), continue; end
        i = candidates(1);
        cmax = c(i);
    end
    if ~isfinite(cmax) || cmax <= 0 || i == 1 || i == numel(c), continue; end
    d = (c(i-1) - c(i+1)) / (2*(c(i-1) - 2*c(i) + c(i+1)));   % parabolic vertex offset
    if ~isfinite(d) || abs(d) > 1, d = 0; end
    E(k) = wm(i) + d*dw;
end
end
