function pt0 = invz_zero_sigma_overlay(pt)
%INVZ_ZERO_SIGMA_OVERLAY Sigma=0 ordered-style point for the bare-RPA overlay (same tl/si as pt).
% Shared by invz_spectra_map / invz_spectra_qpath: the ordered-form invz_chi_realaxis call
% with alpha = alpha_m = 0, lambda = 0, K = [] reduces to the bare (Sigma = 0) RPA response,
% evaluated on the SAME single-ion state (pt.si) and two-level params (pt.tl) as the 1/z pass.
pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
             'K', [], 'is_ordered', true, 'si', pt.si);
end
