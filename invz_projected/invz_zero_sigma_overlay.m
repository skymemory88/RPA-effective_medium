function pt0 = invz_zero_sigma_overlay(pt)
%INVZ_ZERO_SIGMA_OVERLAY Sigma=0 ordered-style point for the bare-RPA overlay (same tl/si as pt).
% This is retained for explicit same-state comparisons. Production projected spectra use
% invz_chi_realaxis opts.bare_rpa instead, because a physical RPA curve must select its own
% bare-MF state rather than inherit the 1/z state.
pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
             'K', [], 'is_ordered', true, 'si', pt.si);
end
