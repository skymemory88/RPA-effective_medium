function [Jnu, info, Jaa0, detail] = invz_bz_couplings(ion, opts)
%INVZ_BZ_COUPLINGS Shared BZ-grid coupling branches + demag-aware transverse J(0).
% Builds the standard BZ q-grid (Gamma dropped), evaluates invz_jq_modes with the
% production dpRng/cache defaults, and hoists Jaa0 = info.Jaa0 (falling back to
% ion.Jxx0 when invz_jq_modes did not supply it). Shared by the sweep/diagnostic
% drivers so the grid and coupling evaluation are defined in exactly one place.
%   opts.grid  ([16 16 16])  BZ q-grid
%   opts.dpRng (30)          real-space dipole cutoff
%   opts.cache (true)        invz_jq_modes file cache
%
% Optional fourth output DETAIL exposes the q-resolved data required by
% nonlocal skeleton diagnostics: qvec, normalized row weights, Hermitian
% Jcc pages, unflattened/flattened sorted eigenvalues, uniform projections,
% and the returned info provenance.  It is nargout-gated; existing one-
% through three-output callers retain the previous invz_jq_modes call.
%
% Ewald Step-5 Task 4 (2026-07-24; operative contract retained here):
%   opts.dipole/opts.ewald   forwarded into invz_jq_modes BY PRESENCE only (absent here means
%                            absent there too -- see invz_jq_modes.m for their exact contract).
%   opts.gridConvention      'legacy_inclusive' | 'halfopen'
%   opts.gridOffset          [1x3] logical/0-1 {0,1/2}^3 phase offset
%   opts.gammaPolicy         'P_complete' | 'P_drop'
% ANY of the three grid-policy fields being PRESENT (regardless of value, including a value
% matching a below-listed default) switches q-grid construction from the legacy
% qVec_generator+Gamma-drop call to invz_phase1_qgrid(ion,grid(1),gridOffset,gridConvention,
% gammaPolicy) -- a cubic grid whose convention/offset/policy VALUES invz_phase1_qgrid alone
% validates; omitted grid-policy fields default to gridConvention='legacy_inclusive',
% gridOffset=[0 0 0], gammaPolicy='P_drop'. opts.dipole/opts.ewald never trigger this switch
% by themselves (grid convention and dipolar backend are orthogonal knobs, freely combinable).
% With NO grid-policy field present, q-grid construction is UNCHANGED and BIT-IDENTICAL to the
% pre-Task-4 behavior, and info.grid remains ABSENT. On the new route, info.grid carries
% complete provenance (schema/convention/offset/gammaPolicy/requested/nominal/retained/
% n_gamma/qhash); explicit gridConvention='legacy_inclusive' still takes this new route and is
% NOT bit-identical to the absent-field route (invz_phase1_qgrid wraps every point via
% mod(q+0.5,1)-0.5, folding the endpoint-inclusive convention's +0.5 face onto -0.5). A
% non-cubic opts.grid on the new route errors invz:bzCouplingsAnisotropicHalfopen
% (invz_phase1_qgrid takes a single scalar N); the absent-grid route remains anisotropic-safe.
%
% A private cacheContext is always forwarded to invz_jq_modes (never exported): kind=
% 'legacy_bz' (absent-grid route; carries the exact opts.grid dimensions, an exact-byte q
% digest, and explicit absent-policy sentinels for gridConvention/gridOffset/gammaPolicy) or
% kind='phase1_qgrid' (new route; carries the fully-resolved grid provenance, identical in
% content to info.grid). This keeps both BZ routes -- and a bare direct invz_jq_modes call,
% which gets invz_jq_modes' own canonical 'direct_call' sentinel -- cache-identity-distinct
% even when they happen to submit byte-identical qvec/physical-parameter payloads.
if nargin < 2, opts = struct(); end
wantDetail = nargout >= 4;
detail = [];
grid  = getf(opts, 'grid', [16 16 16]);
dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', true);

useNewGrid = isfield(opts, 'gridConvention') || isfield(opts, 'gridOffset') || isfield(opts, 'gammaPolicy');
gridInfo = [];
if ~useNewGrid
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    if wantDetail, qweights = ones(size(qc,1),1)/size(qc,1); end
    cacheContext = struct('kind', 'legacy_bz', 'grid', grid(:).', 'qhash', local_qhash(qc), ...
        'gridConvention', '<absent>', 'gridOffset', [], 'gammaPolicy', '<absent>');
else
    if any(grid ~= grid(1))
        error('invz:bzCouplingsAnisotropicHalfopen', ['invz_bz_couplings: gridConvention/gridOffset/' ...
            'gammaPolicy require a cubic opts.grid (invz_phase1_qgrid takes a single scalar N); got ' ...
            'opts.grid = %s.'], mat2str(grid));
    end
    convention  = getf(opts, 'gridConvention', 'legacy_inclusive');
    gridOffset  = getf(opts, 'gridOffset',     [0 0 0]);
    gammaPolicy = getf(opts, 'gammaPolicy',    'P_drop');
    g  = invz_phase1_qgrid(ion, grid(1), gridOffset, convention, gammaPolicy);
    qc = g.qvec;
    if wantDetail, qweights = g.w; end
    gridInfo = struct('schema', 'invz_bz_couplings.grid/v1', 'convention', g.convention, ...
        'offset', g.offsetFlags, 'gammaPolicy', g.gammaPolicy, 'requested', [g.N g.N g.N], ...
        'nominal', g.nominal, 'retained', size(qc, 1), 'n_gamma', g.n_gamma, 'qhash', local_qhash(qc));
    cacheContext = gridInfo;
    cacheContext.kind = 'phase1_qgrid';
end

jqOpts = struct('dpRng', dpRng, 'cache', cache, 'cacheContext', cacheContext);
if isfield(opts, 'dipole'), jqOpts.dipole = opts.dipole; end
if isfield(opts, 'ewald'),  jqOpts.ewald  = opts.ewald;  end
if wantDetail
    [Jnu, info, Juni, Jcc] = invz_jq_modes(ion, qc, jqOpts);
    Jnu_unflat = Jnu;
else
    [Jnu, info] = invz_jq_modes(ion, qc, jqOpts);
end
if ~isempty(gridInfo), info.grid = gridInfo; end
Jnu = Jnu(:);
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end
if wantDetail
    detail = struct('schema','invz_bz_couplings.detail/v1', ...
        'qvec',qc,'weights',qweights,'Jcc',Jcc, ...
        'Jnu_unflat',Jnu_unflat,'Jnu_flat',Jnu, ...
        'Juni',Juni,'info',info, ...
        'flattening','column-major: flat=(branch-1)*nq+q', ...
        'eigenvectors_included',false);
end
end

function h = local_qhash(q)
%LOCAL_QHASH Exact-byte SHA-256 (java.security.MessageDigest, toolbox-free) digest of a
% q-array's class, shape, and exact IEEE-754 contents (Step-5 Task 4 "Cache and
% provenance"): info.grid.qhash must not reuse the weak single-precision-sum hash_vec used
% elsewhere for cache filenames. PROVENANCE ONLY -- invz_jq_modes' own cache-hit validation
% rests on an exact isequaln comparison of the FULL stored qvec, never on this digest's
% collision resistance.
% NOTE: this runs on BOTH BZ routes (including the legacy/absent-grid default), so it requires
% the JVM (as does this plotting-oriented codebase generally); a `matlab -nojvm` batch run is
% unsupported here. It only feeds the private cacheContext, never Jnu/info/Juni.
clsBytes   = uint8(class(q));
shapeBytes = typecast(size(q), 'uint8');   % size(q) is always double
dataBytes  = typecast(reshape(q, 1, []), 'uint8');
allBytes   = typecast([clsBytes(:); shapeBytes(:); dataBytes(:)], 'int8');   % Java byte[] is signed
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(allBytes);
h = sprintf('%02x', typecast(md.digest(), 'uint8'));
end
