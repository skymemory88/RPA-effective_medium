function out = invz_ordered_reduced_residual(y, ctx, opts)
%INVZ_ORDERED_REDUCED_RESIDUAL Four-variable simultaneous ordered residual.
% Diagnostic-only formulation with y=[lambda_1;lambda_2;lambda_3;x].
%
% The static coordinate x fixes
%   Phi(x)=mean_q x/(1+J_q x),
%   K0(x)=mean_q J_q/(1+J_q x) / mean_q 1/(1+J_q x).
% At fixed (lambda,K0), INVZ_ORDERED_DYNAMIC_ELIMINATE solves every n>0
% medium/self-energy relation. The remaining exact residual is
%   R_lambda = lambda - beta^-1 sum_n w_n K_n g_n^p, p=1,2,3,
%   R_static = x-Gtilde0(lambda,K0,Sigma0).
% The latter is exactly equivalent to Phi(x)-Gstat=0 on the physical
% domain, but remains finite across the removable Gstat local pole.
%
% Required CTX fields:
%   G0,g,tl,wts,beta,Jnu_flat,J0eff,G0inel0,G0el0.
% OPTS:
%   Jsup       verified uniform/lattice supremum (default max(J0eff,max Jstatic))
%   mass_tol   strict mass floor (default 0)
%   dynamic    options for INVZ_ORDERED_DYNAMIC_ELIMINATE
if nargin < 3, opts = struct(); end
required = {'G0','g','tl','wts','beta','Jnu_flat','J0eff','G0inel0','G0el0'};
if ~(isstruct(ctx) && isscalar(ctx) && all(isfield(ctx,required)))
    error('invz:reducedResidualContext', ...
        'ctx must be a scalar struct with fields: %s.',strjoin(required,', '));
end
y = y(:);
if numel(y) ~= 4 || any(~isfinite(y)) || ~isreal(y)
    error('invz:reducedResidualShape', 'y must be four finite real values [lambda(1:3);x].');
end
lambda = y(1:3);
x = y(4);
mass_tol = getf(opts,'mass_tol',0);

retarded = ~isvector(ctx.Jnu_flat);
if retarded
    if size(ctx.Jnu_flat,2) ~= numel(ctx.G0)
        error('invz:reducedResidualJ', ...
            'Matrix Jnu_flat must have one column per Matsubara frequency.');
    end
    Jstatic = ctx.Jnu_flat(:,1);
else
    Jstatic = ctx.Jnu_flat(:);
end
Jsup = getf(opts,'Jsup',max([ctx.J0eff;Jstatic]));

out = empty_output();
out.y = y;
out.lambda = lambda;
out.x = x;
out.Jsup = Jsup;
if ~(isfinite(x) && x < 0 && isfinite(1+Jsup*x) && 1+Jsup*x > mass_tol)
    out.status = 'static_x_outside_domain';
    return;
end
Dx = 1+Jstatic*x;
if any(~isfinite(Dx)) || min(Dx) <= mass_tol
    out.status = 'static_lattice_x_mass';
    return;
end

invDx = 1./Dx;
M = mean(invDx);
Phi = x*M;
K0 = mean(Jstatic.*invDx)/M;
out.Phi = Phi;
out.K0 = K0;
out.static_x_min_mass = min(Dx);
out.supremum_mass = 1+Jsup*x;

dopts = getf(opts,'dynamic',struct());
if ~isfield(dopts,'mass_tol'), dopts.mass_tol = mass_tol; end
dyn = invz_ordered_dynamic_eliminate(lambda,K0,ctx,dopts);
out.dynamic = dyn;
out.status = dyn.status;
if ~dyn.defined
    return;
end

Sigma = dyn.Sigma(:);
K = dyn.K(:);
[Gstat,go] = invz_gstat_ordered(ctx.tl,lambda(1:2),K0,Sigma(1), ...
    ctx.beta,ctx.G0inel0,ctx.G0el0,struct('stable_form',true));
Duni = 1+(ctx.J0eff-K0)*Gstat;
Dmesh = 1+(Jstatic-K0)*Gstat;
lambda_residual = lambda-dyn.lambda_check;
raw_static_residual = Phi-Gstat;
static_residual = x-go.Gtil0;
residual = [lambda_residual;static_residual];

xi_physical = isfinite(go.xi) && go.xi >= 0 && ...
    isfinite(go.xi_denom) && go.xi_denom > 0;
finite_reduced_static = all(isfinite([go.Gtil0,go.r,K0,Phi,static_residual]));
finite_physical_static = finite_reduced_static && ...
    all(isfinite([Gstat,go.gstat_local_denom]));
trial_admissible = finite_physical_static && xi_physical && Gstat < 0 && ...
    Duni > mass_tol && min(Dmesh) > mass_tol && ...
    dyn.dynamic_min_lattice_mass > mass_tol && ...
    dyn.dynamic_min_medium_mass > mass_tol;

G = dyn.G(:);
G(1) = Gstat;
if finite_reduced_static
    out.status = 'ok';
    out.defined = true;
else
    out.status = 'static_reduced_undefined';
    out.defined = false;
end
out.residual = residual;
out.residual_norm = max(abs(residual));
out.lambda_residual = lambda_residual;
out.static_residual = static_residual;
out.raw_static_residual = raw_static_residual;
out.Sigma = Sigma;
out.K = K;
out.G = G;
out.Gstat = Gstat;
out.Gtil0 = go.Gtil0;
out.r = go.r;
out.static = go;
out.D_uni = Duni;
out.static_mesh_min_mass = min(Dmesh);
out.x_consistency = static_residual;
out.trial_admissible = trial_admissible;
out.gates = struct('finite_reduced_static',finite_reduced_static, ...
    'finite_physical_static',finite_physical_static,'xi_physical',xi_physical, ...
    'negative_Gstat',isfinite(Gstat) && Gstat < 0, ...
    'positive_uniform_mass',isfinite(Duni) && Duni > mass_tol, ...
    'positive_static_mesh_mass',all(isfinite(Dmesh)) && min(Dmesh) > mass_tol, ...
    'positive_dynamic_lattice_mass', ...
        dyn.dynamic_min_lattice_mass > mass_tol, ...
    'positive_dynamic_medium_mass', ...
        dyn.dynamic_min_medium_mass > mass_tol);

    function s = empty_output()
        s = struct('status','not_evaluated','defined',false,'y',nan(4,1), ...
            'lambda',nan(3,1),'x',NaN,'Jsup',NaN,'Phi',NaN,'K0',NaN, ...
            'static_x_min_mass',NaN,'supremum_mass',NaN,'dynamic',[], ...
            'residual',nan(4,1),'residual_norm',NaN, ...
            'lambda_residual',nan(3,1),'static_residual',NaN, ...
            'raw_static_residual',NaN, ...
            'Sigma',nan(size(ctx.G0(:))),'K',nan(size(ctx.G0(:))), ...
            'G',nan(size(ctx.G0(:))),'Gstat',NaN,'Gtil0',NaN,'r',NaN, ...
            'static',[],'D_uni',NaN,'static_mesh_min_mass',NaN, ...
            'x_consistency',NaN,'trial_admissible',false,'gates',[]);
    end
end
