function ctx = invzt_a1_ordered_context(si, tl, T, lat, opts)
%INVZT_A1_ORDERED_CONTEXT Freeze inputs for a side-effect-free ordered A1 map.
%   CTX = INVZT_A1_ORDERED_CONTEXT(SI, TL, T, LAT, OPTS) constructs one
%   immutable numerical context for INVZT_A1_ORDERED_MAP. This is a diagnostic
%   and residual-building primitive; it does not solve a fixed point or alter a
%   production representation.
%
%   opts.dominant_count is REQUIRED so the selector cannot drift or fall back
%   to an implicit default. The fixed-rank split candidate uses 16 for the full
%   electronuclear state. Setting it to numel(si.E) provides the full-rank
%   reference that algebraically reproduces the current whole-response map;
%   there is deliberately no runtime closure switch. opts.selector_source is
%   also required provenance (normally 'same_point_pm_fixed_rank'). opts.Ecut
%   (40), odd (true), and rank_tol (1e-12) follow the solver conventions.
if nargin < 5, opts = struct(); end
if ~isfield(opts, 'dominant_count') || isempty(opts.dominant_count)
    error('invzt:orderedContextSelector', ...
        'opts.dominant_count is required and must be explicit.');
end
if ~isfield(opts, 'selector_source') || isempty(opts.selector_source)
    error('invzt:orderedContextSelectorSource', ...
        'opts.selector_source is required and must identify the frozen selector producer.');
end
selector_source = opts.selector_source;
if isstring(selector_source) && isscalar(selector_source)
    selector_source = char(selector_source);
end
if ~(ischar(selector_source) && isrow(selector_source) && ~isempty(selector_source))
    error('invzt:orderedContextSelectorSource', ...
        'opts.selector_source must be a nonempty character row or string scalar.');
end
Ecut = getf(opts, 'Ecut', 40);
odd = ~isfield(opts, 'odd') || isempty(opts.odd) || ~isequal(opts.odd, false);
rank_tol = getf(opts, 'rank_tol', 1e-12);
if ~(isstruct(si) && isscalar(si) && isfield(si,'E'))
    error('invzt:orderedContextSI', 'si must be a scalar single-ion state struct.');
end
if ~(isstruct(tl) && isscalar(tl))
    error('invzt:orderedContextTL', 'tl must be a scalar ordered two-level struct.');
end
if ~(isscalar(T) && isreal(T) && isfinite(T) && T > 0)
    error('invzt:orderedContextT', 'T must be a finite positive scalar.');
end
if ~(isstruct(lat) && isscalar(lat) && isfield(lat,'Jt') ...
        && isfield(lat,'JtGamma') && isfield(lat,'w') && isfield(lat,'info') ...
        && isfield(lat,'conv'))
    error('invzt:orderedContextLat', 'lat must be a full INVZT_JQ_TENSOR result.');
end

[wn, wts, beta] = invz_matsubara(T, Ecut);
cfull = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
split_opts = struct('elastic', true, 'dominant_count', opts.dominant_count);
[cdom, crest, mspec] = invzt_chi0_split(si, T, 1i*wn, split_opts);
if mspec.ndom == numel(si.E)
    representation = 'full_rank_reference';
else
    representation = 'fixed_rank_split_candidate';
end
if max(abs(cfull(:) - cdom(:) - crest(:))) > 1e-12*max(1,max(abs(cfull(:))))
    error('invzt:orderedContextDecomposition', ...
        'The selected response decomposition does not reconstruct chi_full.');
end

lat_eff = lat;
if ~odd, lat_eff.Jt = invzt_odd_mask(lat.Jt); end
JG = lat.JtGamma;
odd_ca = JG(3:3:12,1:3:12); odd_cb = JG(3:3:12,2:3:12);
oddmag = max(abs([odd_ca(:);odd_cb(:)]));
if ~(oddmag < 1e-10*abs(lat.info.Jcc0))
    error('invzt:orderedContextOddGamma', ...
        'lat.JtGamma violates the ordered A1 C2 Gamma gate.');
end

ctx = struct();
ctx.schema = 'invzt_a1_ordered_context/v1';
ctx.representation = representation;
ctx.dominant_count = mspec.ndom;
ctx.selector_source = selector_source;
ctx.si = si;
ctx.tl = tl;
ctx.T = T;
ctx.Ecut = Ecut;
ctx.wn = wn;
ctx.wts = wts;
ctx.beta = beta;
ctx.g = real(invz_g(tl, 1i*wn));
ctx.cfull = cfull;
ctx.cdom = cdom;
ctx.crest = crest;
ctx.cdom_cc = real(squeeze(cdom(3,3,:)));
ctx.crest_cc = real(squeeze(crest(3,3,:)));
ctx.mspec = mspec;
ctx.lat = lat_eff;
ctx.odd = odd;
ctx.rank_tol = rank_tol;
ctx.provenance = struct('T',T,'Ecut',Ecut,'odd',odd,'rank_tol',rank_tol, ...
    'dominant_count',mspec.ndom,'selector_source',selector_source, ...
    'lattice_conv',lat.conv,'lattice_info',lat.info);
end

function v = getf(s, name, default)
if isfield(s,name) && ~isempty(s.(name)), v=s.(name); else, v=default; end
end
