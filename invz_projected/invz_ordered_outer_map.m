function out = invz_ordered_outer_map(Sigma, ctx, opts)
%INVZ_ORDERED_OUTER_MAP Deterministic ordered self-energy map at one fixed node.
%   OUT = INVZ_ORDERED_OUTER_MAP(SIGMA,CTX,OPTS) evaluates F[Sigma] without
%   hidden lambda/K0 seeds:
%     1. dynamic K(n>0) is algebraic through invz_emt_scalar;
%     2. lambda_p(K0)=lambda_dyn_p+w0*g0^p*K0/beta is constructed exactly;
%     3. the bounded physical-x static solver evaluates that affine lambda at
%        every K0(x) candidate; and
%     4. a unique admissible root yields Sigma_map=invz_sigma_ordered(...).
%
% If the static map has zero or multiple admissible roots, F is respectively
% undefined or multivalued. OUT.status carries the static status and no Sigma
% map/residual is exported.
%
% Required CTX fields:
%   G0, g, tl, wts, beta, Jnu_flat, J0eff, G0inel0, G0el0
% OPTS:
%   emt                 options for invz_emt_scalar
%   emt_static          bounded-static options (lambda_affine is caller-owned)
%   dynamic_diagnostics compute all n>0 mesh denominator minima (false)
if nargin < 3, opts = struct(); end
required = {'G0','g','tl','wts','beta','Jnu_flat','J0eff','G0inel0','G0el0'};
if ~(isstruct(ctx) && isscalar(ctx) && all(isfield(ctx,required)))
    error('invz:outerMapContext', ...
        'ctx must be a scalar struct with fields: %s.',strjoin(required,', '));
end
Sigma = Sigma(:);
G0 = ctx.G0(:);
g = ctx.g(:);
wts = ctx.wts(:);
nw = numel(Sigma);
if numel(G0) ~= nw || numel(g) ~= nw || numel(wts) ~= nw
    error('invz:outerMapShape', ...
        'Sigma, ctx.G0, ctx.g, and ctx.wts must have equal lengths.');
end
if any(~isfinite(Sigma)) || ~isreal(Sigma)
    error('invz:outerMapShape', 'Sigma must be a finite real vector.');
end

eopts = getf(opts,'emt',struct());
med = invz_emt_scalar(G0,Sigma,ctx.Jnu_flat,eopts);
out = struct('status','dynamic_medium_failed','defined',false, ...
    'Sigma_map',nan(size(Sigma)),'residual',nan(size(Sigma)), ...
    'residual_norm',NaN,'K',nan(size(Sigma)),'lambda',nan(3,1), ...
    'G',nan(size(Sigma)),'static',[],'medium',med, ...
    'lambda_affine',[],'lambda_consistency',NaN, ...
    'dynamic_min_abs',NaN,'dynamic_nonpositive_count',NaN);
if ~med.converged
    return;
end

Kbase = med.K(:);
Kbase(1) = 0;
plist = [1 2 3];
lambda_base = invz_lambdas(Kbase,g,wts,ctx.beta,plist);
lambda_slope = (wts(1)/ctx.beta) * g(1).^plist(:);
la = struct('base',lambda_base,'slope',lambda_slope);

sopts = getf(opts,'emt_static',struct());
if isfield(sopts,'lambda_affine')
    error('invz:outerMapOptions', ...
        'opts.emt_static.lambda_affine is owned by invz_ordered_outer_map.');
end
sopts.lambda_affine = la;
[K0,Gstat,so] = invz_emt_static_ordered(ctx.tl,lambda_base(1:2),Sigma(1), ...
    ctx.Jnu_flat,0,ctx.beta,ctx.J0eff,ctx.G0inel0,ctx.G0el0,sopts);
out.static = so;
out.lambda_affine = la;
out.status = so.status;
if ~so.converged
    return;
end

K = med.K(:);
K(1) = K0;
lambda = lambda_base + lambda_slope*K0;
lambda_check = invz_lambdas(K,g,wts,ctx.beta,plist);
lambda_consistency = max(abs(lambda-lambda_check));
if ~(isfinite(lambda_consistency) && lambda_consistency <= ...
        64*eps(max(1,max(abs(lambda_check)))))
    out.status = 'lambda_affine_inconsistent';
    out.lambda_consistency = lambda_consistency;
    return;
end
sg = invz_sigma_ordered(ctx.tl,lambda,K,g,ctx.beta);
Sigma_map = sg.Sigma(:);
if any(~isfinite(Sigma_map))
    out.status = 'outer_map_nonfinite';
    return;
end

G = med.G(:);
G(1) = Gstat;
out.status = 'ok';
out.defined = true;
out.Sigma_map = Sigma_map;
out.residual = Sigma-Sigma_map;
out.residual_norm = max(abs(out.residual));
out.K = K;
out.lambda = lambda;
out.G = G;
out.sigma_diagnostics = sg;
out.lambda_consistency = lambda_consistency;

if getf(opts,'dynamic_diagnostics',false)
    Jf = ctx.Jnu_flat(:);
    minabs = Inf;
    nnonpos = 0;
    block = getf(opts,'dynamic_block',256);
    for i0 = 2:block:nw
        idx = i0:min(i0+block-1,nw);
        D = 1 + (Jf-K(idx).').*G(idx).';
        minabs = min(minabs,min(abs(D),[],'all'));
        nnonpos = nnonpos + nnz(D <= 0);
    end
    out.dynamic_min_abs = minabs;
    out.dynamic_nonpositive_count = nnonpos;
end
end
