function out = invzt_a1_ordered_map(ctx, Sigma)
%INVZT_A1_ORDERED_MAP Evaluate one ordered scalar-A1 map transactionally.
%   OUT = INVZT_A1_ORDERED_MAP(CTX, SIGMA) evaluates, without mixing or state
%   mutation, the exact ordered moment-form map declared by CTX:
%       chi_tilde = chi_dom/(1+Sigma) + chi_rest
%       K = 1/G_loc - 1/G0_tilde
%       Sigma_next = INVZ_SIGMA_ORDERED(...)
%       residual = Sigma_next - Sigma.
%   CTX is produced by INVZT_A1_ORDERED_CONTEXT and carries an explicit fixed
%   rank; full rank is the current whole-response reference. Domain failures
%   return a classified OUT rather than a finite pseudo-root. Structural input
%   errors still throw.
if ~(isstruct(ctx) && isscalar(ctx) && isfield(ctx,'schema') ...
        && strcmp(ctx.schema,'invzt_a1_ordered_context/v1'))
    error('invzt:orderedMapContext', ...
        'ctx must be an invzt_a1_ordered_context/v1 struct.');
end
nwn = numel(ctx.wn);
if ~(isnumeric(Sigma) && isreal(Sigma) && numel(Sigma)==nwn)
    error('invzt:orderedMapSigma', ...
        'Sigma must contain %d real entries.', nwn);
end
Sigma = Sigma(:);
out = empty_out(ctx, Sigma);
if any(~isfinite(Sigma))
    out.status = 'sigma_nonfinite'; out.failure = out.status; return;
end
out.min_one_plus_Sigma = min(1 + Sigma);
if out.min_one_plus_Sigma <= 0
    out.status = 'sigma_domain'; out.failure = out.status; return;
end

denom = reshape(1 + Sigma, 1, 1, nwn);
ctil = ctx.cdom ./ denom + ctx.crest;
if any(~isfinite(ctil(:)))
    out.status = 'response_nonfinite'; out.failure = out.status; return;
end
[Gcc, diag4] = invzt_gcc_lattice(ctil, ctx.lat);
G = -Gcc(:);
G0til = -(ctx.cdom_cc ./ (1 + Sigma) + ctx.crest_cc);
if any(~isfinite(G)) || any(~isfinite(G0til)) || any(G==0) || any(G0til==0)
    out.status = 'medium_nonfinite'; out.failure = out.status; return;
end
K = 1 ./ G - 1 ./ G0til;
if any(~isfinite(K))
    out.status = 'medium_nonfinite'; out.failure = out.status; return;
end
lambda = invz_lambdas(K, ctx.g, ctx.wts, ctx.beta, [1 2 3]);
sg = invz_sigma_ordered(ctx.tl, lambda, K, ctx.g, ctx.beta);
residual = sg.Sigma(:) - Sigma;
if any(~isfinite(residual))
    out.status = 'map_nonfinite'; out.failure = out.status; return;
end

ctil0 = real((ctil(:,:,1)+ctil(:,:,1)')/2);
[crit, cmass, arank] = invzt_crit_static( ...
    ctil0, ctx.lat.JtGamma, ctx.rank_tol);
out.status = 'evaluated';
out.valid = true;
out.failure = '';
out.cdom = ctx.cdom;
out.crest = ctx.crest;
out.chi_tilde = ctil;
out.G = G;
out.G0_tilde = G0til;
out.K = K;
out.lambda = lambda;
out.alpha = sg.alpha;
out.alpha_m = sg.alpha_m;
out.gamma = sg.gamma;
out.Sigma_next = sg.Sigma(:);
out.residual = residual;
out.residual_inf = max(abs(residual));
out.residual_scaled_inf = out.residual_inf / ...
    max([1; abs(Sigma); abs(out.Sigma_next)]);
out.diag4 = diag4;
out.diag4_spread = max(max(diag4,[],1)-min(diag4,[],1));
out.crit = crit;
out.crit_clipped_mass = cmass;
out.crit_active_rank = arank;
out.local_static_eigenvalues = sort(real(eig(ctil0)));
out.dominant_count = ctx.dominant_count;
out.mspec = ctx.mspec;
out.provenance = ctx.provenance;
end

function out = empty_out(ctx, Sigma)
out = struct('schema','invzt_a1_ordered_map/v1', ...
    'representation',ctx.representation,'valid',false,'failure','', ...
    'status','not_evaluated','Sigma',Sigma,'min_one_plus_Sigma',NaN, ...
    'cdom',[],'crest',[],'chi_tilde',[],'G',[],'G0_tilde',[],'K',[], ...
    'lambda',[],'alpha',NaN,'alpha_m',NaN,'gamma',[],'Sigma_next',[], ...
    'residual',[],'residual_inf',NaN,'residual_scaled_inf',NaN, ...
    'diag4',[],'diag4_spread',NaN,'crit',NaN,'crit_clipped_mass',NaN, ...
    'crit_active_rank',NaN,'local_static_eigenvalues',[], ...
    'dominant_count',ctx.dominant_count,'mspec',ctx.mspec, ...
    'provenance',ctx.provenance);
end
