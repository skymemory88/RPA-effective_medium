function [Gref, ref] = invz_static_medium_reference(G0bare, Sigma0, scheme, opts)
%INVZ_STATIC_MEDIUM_REFERENCE Reference propagator for the strict-O(1/z) static medium
% (spec SS4.1). G = -chi (meV^-1), ferromagnetic positive J.
%   'strict_1z_dyson_ref' : denom = 1 + Sigma0;  Gref = G0bare/denom   (SELECTED convention)
%   'strict_1z_bare_ref'  : denom = 1;           Gref = G0bare         (systematic comparator)
%
% WHY THIS IS A SEPARATE PRIMITIVE: invz_medium_moment_closure receives only the QUOTIENT
% Gref, so it cannot reconstruct the denominator or report its margin. Reference construction
% therefore cannot live inside the closure.
%
% The Dyson convention is O(1/z)-equivalent to the bare one (Sigma0 is itself O(1/z)), so it
% is a Dyson-IMPROVED scheme choice, not uniquely the literal first-order expansion. It is
% selected because it makes the PM-leg expression exactly G0/D with invz_emt_scalar's own
% D = 1 + Sigma, so both legs' expressions are textually the same object. Gref carries NO K0,
% lambda or xi dependence -- that is what makes the closure one-shot.
%
% Domain events RETURN a status and Gref = NaN; they never throw (spec SS5.2). A caller must
% not evaluate the closure on a non-'ok' reference. An unknown or 'resummed' scheme IS a
% wiring error and throws invz:staticMedium -- the resummed path must bypass this primitive
% entirely.
% opts.ref_margin (default 1e-6; named ref_floor in the preregistration): denom at or below
% this is 'ref_denom_small'. 1 + Sigma0 is O(1) because Sigma0 is O(1/z), so a
% denominator this small means the reference is degenerate.
if nargin < 4 || isempty(opts), opts = struct(); end
margin = getf(opts, 'ref_margin', 1e-6);
if ~(isscalar(margin) && isfinite(margin) && margin > 0)
    error('invz:staticMedium', 'opts.ref_margin must be a positive finite scalar.');
end
switch scheme
    case 'strict_1z_dyson_ref'
        denom = 1 + Sigma0;
    case 'strict_1z_bare_ref'
        denom = 1;
    otherwise
        error('invz:staticMedium', ['invz_static_medium_reference: scheme must be ' ...
            '''strict_1z_dyson_ref'' or ''strict_1z_bare_ref''; got ''%s''. The ' ...
            '''resummed'' path must not call this primitive.'], scheme);
end
ref = struct('denom', denom, 'floor', margin, 'margin', denom - margin, ...
             'status', 'ok', 'scheme', scheme);
if ~isfinite(G0bare) || ~isfinite(denom)
    ref.status = 'nonfinite';
elseif denom <= 0
    ref.status = 'ref_denom_nonpositive';
elseif denom <= margin
    ref.status = 'ref_denom_small';
end
if strcmp(ref.status, 'ok')
    Gref = G0bare/denom;
else
    Gref = NaN;
end
end
