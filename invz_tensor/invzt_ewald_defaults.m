function eopts = invzt_ewald_defaults(ion)
%INVZT_EWALD_DEFAULTS  Certified full-tensor Ewald production controls.
%
%   eopts = invzt_ewald_defaults(ion)
%
% Returns the exact four-field control struct consumed by
% INVZ_DIPOLE_EWALD/INVZT_JQ_TENSOR.  The dimensionless cutoffs are the
% frozen primitive-calibration defaults Cr=5.5 and Cg=11; alpha is derived
% deterministically from the crystallographic cell volume.

if nargin ~= 1 || ~isstruct(ion) || ~isscalar(ion) ...
        || ~isfield(ion, 'Vc') || ~isnumeric(ion.Vc) || ~isreal(ion.Vc) ...
        || ~isscalar(ion.Vc) || ~isfinite(ion.Vc) || ion.Vc <= 0
    error('invzt:ewaldDefaultsIon', ...
        'ion must be a scalar struct with a finite positive scalar Vc.');
end

alpha0 = sqrt(pi) / ion.Vc^(1/3);
eopts = struct('alpha', alpha0, 'r_cut', 5.5/alpha0, ...
    'g_cut', 11*alpha0, 'boundary', 'conducting_k0_omitted');
end
