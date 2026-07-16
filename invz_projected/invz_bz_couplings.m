function [Jnu, info, Jaa0] = invz_bz_couplings(ion, opts)
%INVZ_BZ_COUPLINGS Shared BZ-grid coupling branches + demag-aware transverse J(0).
% Builds the standard BZ q-grid (Gamma dropped), evaluates invz_jq_modes with the
% production dpRng/cache defaults, and hoists Jaa0 = info.Jaa0 (falling back to
% ion.Jxx0 when invz_jq_modes did not supply it). Shared by the sweep/diagnostic
% drivers so the grid and coupling evaluation are defined in exactly one place.
%   opts.grid  ([16 16 16])  BZ q-grid
%   opts.dpRng (30)          real-space dipole cutoff
%   opts.cache (true)        invz_jq_modes file cache
if nargin < 2, opts = struct(); end
grid  = getf(opts, 'grid', [16 16 16]);
dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', true);
[qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
qc = qc(any(abs(qc) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', cache));
Jnu = Jnu(:);
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end
end
