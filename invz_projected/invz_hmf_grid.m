function [hgrid, ratio] = invz_hmf_grid(hmax, nH, hmin_frac)
%INVZ_HMF_GRID The initial geometric H_MF profile grid, clustered at 0 (P1-4).
% Extracted VERBATIM from invz_hmf_ordered.m (no line number: the reference went stale twice
% in two consecutive commits as the caller grew; the call site is unambiguous by name) so the
% prospective Gate-0 scanner
% (invz_static_domain_scan) and the solver build the SAME initial grid from one definition
% rather than two implementations that happen to agree (spec SS7.2).
%
% The scanner must NOT reproduce invz_hmf_ordered's adaptive extension or redensification: this
% helper covers the INITIAL grid only, and solved-path margins are read off prof.hgrid AFTER
% adaptation.
%   hgrid = hmax * ratio.^((nH-1):-1:0),  ratio = hmin_frac^(1/(nH-1))
% ascending, hgrid(end) = hmax, hgrid(1) = hmax*hmin_frac.
if ~(isscalar(hmax) && isfinite(hmax) && hmax > 0)
    error('invz:hmfGrid', 'hmax must be a positive finite scalar; got %s', mat2str(hmax));
end
if ~(isscalar(nH) && nH == round(nH) && nH >= 2)
    error('invz:hmfGrid', 'nH must be an integer >= 2; got %s', mat2str(nH));
end
if ~(isscalar(hmin_frac) && isfinite(hmin_frac) && hmin_frac > 0 && hmin_frac < 1)
    error('invz:hmfGrid', 'hmin_frac must be in (0,1); got %s', mat2str(hmin_frac));
end
ratio = hmin_frac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);
end
