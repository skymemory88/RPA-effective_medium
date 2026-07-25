function [K0, info] = invz_medium_moment_closure(Gref, mom, scheme)
%INVZ_MEDIUM_MOMENT_CLOSURE One-shot strict-O(1/z) static effective medium (spec SS4.1, SSB).
% G = -chi (meV^-1), ferromagnetic positive J.
%
%   K0 = Jbar - mu2*Gref                                          (ONE-SHOT, no feedback)
%
% Both legs call THIS function with an already-constructed Gref (invz_static_medium_reference),
% which is what makes same-order truncation structural rather than a coincidence of two edits:
% the exact closures of invz_emt_scalar and invz_emt_static_ordered are the SAME function of
% (local G; moments) --
%   K = Jbar - mu2*G + mu3*G^2 + (2*mu2^2 - mu4)*G^3 + (mu5 - 5*mu2*mu3)*G^4 + O(G^5)
% -- verified term-by-term to O(G^4). Truncating both at mu2 with the same Gref makes the
% m -> 0 cross-phase identity exact within the scheme.
%
% ORDER ACCOUNTING: under the high-density counting mu2 ~ 1/z, mu2*Gref is the O(1/z) medium
% correction. Solving K0 = Jbar - mu2*G(K0) self-consistently would re-admit denominator
% feedback of the same class that exceeds retained order, so it is NOT done here (the
% self-consistent quadratic is a separately named diagnostic comparator, not a scheme).
%
% info.omit_mu3   = |mu3*Gref^2|            / |mu2*Gref|   <-- the FIRST omitted term
% info.omit_cubic = |(2*mu2^2-mu4)*Gref^3|  / |mu2*Gref|
% info.omit_max   = max of the two, but NaN whenever EITHER ratio is NaN (fails closed: MATLAB's
%                   max ignores a NaN operand, which would otherwise mask a corrupted moment
%                   behind a finite omit_max on an otherwise-'ok' node). Inf still propagates
%                   through omit_max unchanged -- isnan(Inf) is false, so the zero-denominator
%                   convention above is preserved.
% Both are ALWAYS reported. mu3 is near zero on the production multiset, but that is a measured
% property of one grid/cutoff/backend -- generalising it is the same inference error that
% produced the synthetic-Jnu defect this work repairs. Zero convention: if mu2*Gref == 0, a
% ratio is 0 when its own numerator is also 0, else Inf.
%
% This leaf NEVER rejects a node on a large ratio: the truncated polynomial stays defined, and
% the frozen omit_report/omit_promote thresholds (docs/invzp_strict_medium_prereg.md SS4) are the
% CALLER's promotion gate. A large ratio must never trigger a scheme switch.
% An unknown or 'resummed' scheme is a wiring error (invz:staticMedium); the resummed path
% bypasses this primitive entirely.
if ~any(strcmp(scheme, {'strict_1z_dyson_ref', 'strict_1z_bare_ref'}))
    error('invz:staticMedium', ['invz_medium_moment_closure: scheme must be ' ...
        '''strict_1z_dyson_ref'' or ''strict_1z_bare_ref''; got ''%s''.'], scheme);
end
req = {'Jbar', 'mu2', 'mu3', 'mu4'};
for k = 1:numel(req)
    if ~isfield(mom, req{k}) || ~isscalar(mom.(req{k}))
        error('invz:staticMedium', ['invz_medium_moment_closure: mom.%s must be a SCALAR. ' ...
            'For a [nJ,nw] retardation moment set, pass the static column (index 1) -- ' ...
            'never the whole set (spec SS4.3).'], req{k});
    end
end
if ~isnumeric(Gref) || ~isscalar(Gref)
    error('invz:staticMedium', ['invz_medium_moment_closure: Gref must be a SCALAR. ' ...
        'For a [nJ,nw] retardation reference set, pass the static column (index 1) -- ' ...
        'never the whole set (spec SS4.3).']);
end

corr = mom.mu2*Gref;                                   % the retained O(1/z) correction
num3 = abs(mom.mu3*Gref^2);
num4 = abs((2*mom.mu2^2 - mom.mu4)*Gref^3);
den  = abs(corr);
info = struct('scheme', scheme, 'retained', 'mu2', 'Kstrict', NaN, ...
              'omit_mu3', ratio(num3, den), 'omit_cubic', ratio(num4, den), ...
              'omit_max', NaN, 'status', 'ok');
if any(isnan([info.omit_mu3, info.omit_cubic]))
    info.omit_max = NaN;
else
    info.omit_max = max(info.omit_mu3, info.omit_cubic);
end

if ~isfinite(Gref) || ~isfinite(corr)
    info.status = 'nonfinite';
    K0 = NaN;
else
    K0 = mom.Jbar - corr;
    if ~isfinite(K0), info.status = 'nonfinite'; end
end
info.Kstrict = K0;      % one-shot: the checked value IS the returned value, by construction
end

% ---------------------------------------------------------------------------------------------
function r = ratio(num, den)
%RATIO Omitted-term ratio with the explicit zero convention (spec SS4.1): 0/0 is 0, x/0 is Inf.
if den == 0
    if num == 0, r = 0; else, r = Inf; end
else
    r = num/den;
end
end
