function out = invzt_a1_ordered_realaxis_diagnostic(ion, ctx, map, w, opts)
%INVZT_A1_ORDERED_REALAXIS_DIAGNOSTIC Same-representation A1 continuation.
%   OUT = INVZT_A1_ORDERED_REALAXIS_DIAGNOSTIC(ION, CTX, MAP, W, OPTS)
%   continues a converged frozen INVZT_A1_ORDERED_MAP co-state to
%   z=W+i*eta using the same explicit fixed-rank split carried by CTX:
%       chi_tilde(z) = chi_dom(z)/(1+Sigma_R(z)) + chi_rest(z).
%   Full rank is the current whole-response reference. The function is a WP1
%   proof object: it cannot select a phase, solve a field/state, or populate
%   production spectra.
%
%   opts: eta (5e-3), qsel ('gamma_uniform'|'gamma'|[nq,3]), cache
%   (false), residual_tol (1e-8), decomposition_tol (1e-12), denom_tol
%   (1e-12), and force_sigma0 (false, exact bare-response test hook).
if nargin < 5, opts = struct(); end
if ~(isstruct(ctx) && isscalar(ctx) && isfield(ctx,'schema') ...
        && strcmp(ctx.schema,'invzt_a1_ordered_context/v1'))
    error('invzt:orderedRealaxisContext', ...
        'ctx must be an invzt_a1_ordered_context/v1 struct.');
end
if ~(isstruct(map) && isscalar(map) && isfield(map,'schema') ...
        && strcmp(map.schema,'invzt_a1_ordered_map/v1') && map.valid)
    error('invzt:orderedRealaxisMap', ...
        'map must be a valid invzt_a1_ordered_map/v1 result.');
end
if map.dominant_count~=ctx.dominant_count ...
        || ~strcmp(map.representation,ctx.representation)
    error('invzt:orderedRealaxisRepresentation', ...
        'map representation and dominant rank must match ctx exactly.');
end
residual_tol = getf(opts,'residual_tol',1e-8);
if ~(isfinite(map.residual_inf) && map.residual_inf<=residual_tol)
    error('invzt:orderedRealaxisResidual', ...
        'The frozen map residual %.3e exceeds residual_tol %.3e.', ...
        map.residual_inf,residual_tol);
end
if ~(isnumeric(w) && isreal(w) && all(isfinite(w(:))))
    error('invzt:orderedRealaxisFrequency','w must contain finite real frequencies.');
end
eta = getf(opts,'eta',5e-3);
decomp_tol = getf(opts,'decomposition_tol',1e-12);
denom_tol = getf(opts,'denom_tol',1e-12);
cacheq = getf(opts,'cache',false);
qsel = getf(opts,'qsel','gamma_uniform');
force_sigma0 = isfield(opts,'force_sigma0') && ~isempty(opts.force_sigma0) ...
    && ~isequal(opts.force_sigma0,false);
if ~(isscalar(eta) && isreal(eta) && isfinite(eta) && eta>0)
    error('invzt:orderedRealaxisEta','opts.eta must be a finite positive scalar.');
end
if ~(isscalar(decomp_tol) && isreal(decomp_tol) && isfinite(decomp_tol) ...
        && decomp_tol>0 && isscalar(denom_tol) && isreal(denom_tol) ...
        && isfinite(denom_tol) && denom_tol>0)
    error('invzt:orderedRealaxisTolerance', ...
        'decomposition_tol and denom_tol must be finite positive scalars.');
end
if isfield(opts,'odd') && ~isequal(opts.odd,ctx.odd)
    error('invzt:orderedRealaxisOdd', ...
        'opts.odd must equal the frozen context convention.');
end

w = w(:); z = w+1i*eta; nw = numel(w);
tl = ctx.tl; g = invz_g(tl,z);
if force_sigma0
    Sw = zeros(nw,1);
else
    pref = tl.M2/tl.n01^2;
    K0 = map.K(1);
    gamma_w = pref*(map.lambda(1)-(1-tl.n01^2)*K0);
    gamma0 = gamma_w;
    if tl.M2==0
        Q0 = (2*tl.m^2/tl.n01^2) ...
            *(map.lambda(1)-(1-tl.n01^2)*K0);
    else
        Q0 = (2*tl.m^2/tl.M2)*gamma0;
    end
    Sw = (map.alpha-map.alpha_m)+(gamma_w-Q0).*g;
end
den = 1+Sw;
if any(~isfinite(den)) || any(abs(den)<=denom_tol)
    error('invzt:orderedRealaxisDomain', ...
        'The real-axis 1+Sigma denominator is nonfinite or below denom_tol.');
end

split_opts = struct('elastic',false,'dominant_count',ctx.dominant_count);
[cdom,crest,mspec] = invzt_chi0_split(ctx.si,ctx.T,z,split_opts);
cfull = invz_chi0z(ctx.si,ctx.T,z,struct('elastic',false));
decomp_error = max(abs(cfull(:)-cdom(:)-crest(:)));
decomp_scale = max(1,max(abs(cfull(:))));
if decomp_error>decomp_tol*decomp_scale
    error('invzt:orderedRealaxisDecomposition', ...
        'chi_full=chi_dom+chi_rest failed (error %.3e).',decomp_error);
end
ctil = cdom./reshape(den,1,1,nw)+crest;

Jaa0 = ctx.lat.info.Jaa0; Jcc0 = ctx.lat.info.Jcc0;
u = kron(ones(4,1)/2,eye(3));
Jd = kron(ones(4)/4,diag([Jaa0,Jaa0,Jcc0]));
Jpage = Jd; explicitq = false; latq = struct(); nq = 0;
if ischar(qsel) || (isstring(qsel) && isscalar(qsel))
    qlabel = char(qsel);
    switch qlabel
        case 'gamma_uniform'
            Jpage = Jd;
        case 'gamma'
            Jpage = ctx.lat.JtGamma;
            if ~ctx.odd, Jpage=invzt_odd_mask(Jpage); end
        otherwise
            error('invzt:orderedRealaxisQsel', ...
                'qsel must be ''gamma_uniform'', ''gamma'', or an [nq,3] array.');
    end
elseif isnumeric(qsel) && ismatrix(qsel) && size(qsel,2)==3 ...
        && all(isfinite(qsel(:)))
    qlabel = 'explicit'; explicitq = true;
    nq = size(qsel,1);
    jqopts = lattice_options(ctx.lat,opts,cacheq);
    latq = invzt_jq_tensor(ion,qsel,jqopts);
    if ~ctx.odd, latq.Jt=invzt_odd_mask(latq.Jt); end
else
    error('invzt:orderedRealaxisQsel', ...
        'qsel must be ''gamma_uniform'', ''gamma'', or a finite [nq,3] array.');
end

chi_uniform = complex(zeros(3,3,nw));
chi_cc_q = complex(zeros(nq,nw));
for k=1:nw
    X = invzt_chi_rpa(ctil(:,:,k),Jpage);
    chi_uniform(:,:,k)=u'*X*u;
    if explicitq
        Xq=invzt_chi_rpa(ctil(:,:,k),latq.Jt);
        for iq=1:nq
            acc=0;
            for s=1:4
                acc=acc+Xq(3*(s-1)+3,3*(s-1)+3,iq);
            end
            chi_cc_q(iq,k)=acc/4;
        end
    end
end

out = struct();
out.schema = 'invzt_a1_ordered_realaxis_diagnostic/v1';
out.classification = 'same_rep_candidate';
out.representation = ctx.representation;
out.selector_source = ctx.selector_source;
out.dominant_count = ctx.dominant_count;
out.mspec = mspec;
out.map_residual_inf = map.residual_inf;
out.w = w; out.eta = eta; out.qsel = qlabel;
out.Sigma_w = Sw;
out.min_abs_one_plus_Sigma_w = min(abs(den));
out.cdom = cdom; out.crest = crest; out.chi_tilde = ctil;
out.decomposition_error = decomp_error;
out.chi_uniform = chi_uniform; out.chi_cc_q = chi_cc_q;
out.lattice_provenance = ctx.provenance;
out.explicit_q_dipole = [];
if explicitq, out.explicit_q_dipole=latq.info.dipole; end
out.digests = struct('Sigma_w',complex_digest(Sw), ...
    'cdom',complex_digest(cdom),'crest',complex_digest(crest), ...
    'chi_tilde',complex_digest(ctil));
end

function jqopts = lattice_options(lat,opts,cacheq)
jqopts = struct('cache',cacheq);
if isfield(lat.info,'dipole')
    p=lat.info.dipole; jqopts.dipole=p.backend;
    if strcmp(p.backend,'ewald')
        jqopts.ewald=p.ewald;
    else
        if isfield(opts,'dpRng') && ~isempty(opts.dpRng) ...
                && opts.dpRng~=lat.info.dpRng
            error('invzt:orderedRealaxisDipoleMismatch', ...
                'opts.dpRng does not match the frozen lattice.');
        end
        jqopts.dpRng=lat.info.dpRng;
    end
else
    error('invzt:orderedRealaxisDipoleProvenance', ...
        'Explicit q requires recorded dipole backend provenance.');
end
end

function d = complex_digest(x)
d=struct('size',size(x),'real',invz_exact_numeric_digest(real(x)), ...
    'imag',invz_exact_numeric_digest(imag(x)));
end

function v = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name)), v=s.(name); else, v=default; end
end
