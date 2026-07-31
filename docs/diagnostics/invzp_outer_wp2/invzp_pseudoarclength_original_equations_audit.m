function result = invzp_pseudoarclength_original_equations_audit()
%INVZP_PSEUDOARCLENGTH_ORIGINAL_EQUATIONS_AUDIT Unreduced fold/branch check.
% Diagnostic only.  Construct two roots at one fixed h by iterating the full
% Sigma vector with the original EMT and self-energy equations at fixed
% s=-J_sup*x, then bisect only the pole-cancelled static residual.  The
% reduced scalar dynamic eliminator and its least-squares moment solve do not
% construct these roots.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
A = load(fullfile(here,'wp2_reduced_pseudoarclength_1t.mat'));
R = load(fullfile(here,'wp2_reduced_fold_refinement.mat'));
T = A.result.T;
Bx = A.result.Bx;
Jsup = A.result.Jsup;
opts = struct('mix',0.5,'tol',1e-12,'maxit',250, ...
    'direct_dynamic_checks',false);
checkopts = opts;
checkopts.direct_dynamic_checks = true;

% At h=0.006 the traced component crosses the fixed-h section on both sides
% of the fold.  Broad disjoint brackets prevent the two roots from being
% seeded by the same local continuation point.
hsection = 0.006;
ctx = make_context(Bx,hsection,T,F);
brackets = [0.35 0.50;0.75 0.83];
roots = repmat(empty_record(numel(ctx.G0)),size(brackets,1),1);
bracket_records = repmat(struct('left',[],'right',[]),size(brackets,1),1);
Sigma = zeros(size(ctx.G0));
for k = 1:size(brackets,1)
    left = direct_solve(ctx,brackets(k,1),Jsup,Sigma,opts);
    right = direct_solve(ctx,brackets(k,2),Jsup,left.Sigma,opts);
    fprintf('bracket %d: converged %d/%d, Rstatic %.9g / %.9g\n', ...
        k,left.converged,right.converged,left.static_residual, ...
        right.static_residual);
    assert(left.converged && right.converged && ...
        sign(left.static_residual) ~= sign(right.static_residual));
    bracket_records(k) = struct('left',left,'right',right);
    roots(k) = bisect_static(ctx,Jsup,brackets(k,1),brackets(k,2), ...
        left.Sigma,opts);
    roots(k) = direct_solve(ctx,-Jsup*roots(k).x,Jsup, ...
        roots(k).Sigma,checkopts);
    Sigma = roots(k).Sigma;
end

comparison = repmat(struct('reduced_status',"", ...
    'reduced_residual_norm',NaN,'Sigma_delta',NaN,'K_delta',NaN, ...
    'outer_status',"",'outer_residual_norm',NaN, ...
    'supremum_mass',NaN,'uniform_mass',NaN),numel(roots),1);
for k = 1:numel(roots)
    q = invz_ordered_reduced_residual([roots(k).lambda;roots(k).x], ...
        ctx,struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12)));
    outer = invz_ordered_outer_map(roots(k).Sigma,ctx, ...
        struct('emt_static',struct('Jsup',Jsup,'warn',false)));
    comparison(k) = struct('reduced_status',string(q.status), ...
        'reduced_residual_norm',q.residual_norm, ...
        'Sigma_delta',max(abs(q.Sigma-roots(k).Sigma)), ...
        'K_delta',max(abs(q.K-roots(k).K)), ...
        'outer_status',string(outer.status), ...
        'outer_residual_norm',outer.residual_norm, ...
        'supremum_mass',q.supremum_mass,'uniform_mass',q.D_uni);
end

% At the independently refined turning point, verify that full-Sigma Picard
% from both the reduced state and a zero state reaches the same original-
% equation solution.
ctxf = make_context(Bx,R.result.h,T,F);
qf = invz_ordered_reduced_residual(R.result.fold.y,ctxf, ...
    struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12)));
fold_from_reduced = direct_solve(ctxf,R.result.s,Jsup,qf.Sigma,checkopts);
fold_from_zero = direct_solve(ctxf,R.result.s,Jsup, ...
    zeros(size(ctxf.G0)),checkopts);
outerf = invz_ordered_outer_map(qf.Sigma,ctxf, ...
    struct('emt_static',struct('Jsup',Jsup,'warn',false)));
fold_check = struct('s',R.result.s,'h',R.result.h, ...
    'reduced_residual_norm',qf.residual_norm, ...
    'from_reduced',fold_from_reduced,'from_zero',fold_from_zero, ...
    'seed_sigma_delta',max(abs( ...
        fold_from_reduced.Sigma-fold_from_zero.Sigma)), ...
    'reduced_sigma_delta',max(abs(qf.Sigma-fold_from_reduced.Sigma)), ...
    'outer_status',string(outerf.status), ...
    'outer_residual_norm',outerf.residual_norm);

result = struct('T',T,'Bx',Bx,'Jsup',Jsup,'hsection',hsection, ...
    'brackets',brackets,'bracket_records',bracket_records,'roots',roots, ...
    'comparison',comparison,'fold_check',fold_check, ...
    'sources',struct('continuation',fullfile(here, ...
        'wp2_reduced_pseudoarclength_1t.mat'),'fold',fullfile(here, ...
        'wp2_reduced_fold_refinement.mat')), ...
    'note',['The two fixed-h roots and the fold state are checked through ' ...
    'full-Sigma Picard and direct original EMT substitution. The two brackets ' ...
    'prove at least two roots in the sampled section, not global uniqueness.']);
save(fullfile(here, ...
    'wp2_pseudoarclength_original_equations_audit.mat'),'result','-v7');

fprintf(['original-equation branch audit: h %.6g roots s %.12g / %.12g, ' ...
    'Rstatic %.3g / %.3g, outer %.3g / %.3g; fold seed delta %.3g, ' ...
    'fold outer %.3g\n'],hsection,-Jsup*roots(1).x,-Jsup*roots(2).x, ...
    roots(1).static_residual,roots(2).static_residual, ...
    comparison(1).outer_residual_norm,comparison(2).outer_residual_norm, ...
    fold_check.seed_sigma_delta,fold_check.outer_residual_norm);
end

function rec = bisect_static(ctx,Jsup,a,b,Sigma,opts)
ra = direct_solve(ctx,a,Jsup,Sigma,opts);
rb = direct_solve(ctx,b,Jsup,ra.Sigma,opts);
assert(ra.converged && rb.converged && ...
    sign(ra.static_residual) ~= sign(rb.static_residual));
rec = ra;
for it = 1:100
    m = a+(b-a)/2;
    rm = direct_solve(ctx,m,Jsup,rec.Sigma,opts);
    assert(rm.converged);
    if abs(rm.static_residual) < abs(rec.static_residual), rec = rm; end
    if sign(ra.static_residual) ~= sign(rm.static_residual)
        b = m;
        rb = rm; %#ok<NASGU>
    else
        a = m;
        ra = rm;
    end
    if abs(rec.static_residual) <= 1e-11 || b-a <= 1e-12, break; end
end
end

function rec = direct_solve(ctx,s,Jsup,Sigma,opts)
J = static_J(ctx);
x = -s/Jsup;
dx = 1+J*x;
M = mean(1./dx);
Phi = x*M;
K0 = mean(J./dx)/M;
map_residual = Inf;
for it = 1:opts.maxit
    med = invz_emt_scalar(ctx.G0,Sigma,ctx.Jnu_flat,struct());
    if ~med.converged, break; end
    K = med.K(:);
    K(1) = K0;
    lambda = invz_lambdas(K,ctx.g,ctx.wts,ctx.beta,[1 2 3]);
    sig = invz_sigma_ordered(ctx.tl,lambda,K,ctx.g,ctx.beta);
    map_residual = max(abs(sig.Sigma-Sigma));
    if map_residual <= opts.tol, break; end
    Sigma = (1-opts.mix)*Sigma+opts.mix*sig.Sigma;
end
converged = isfinite(map_residual) && map_residual <= opts.tol;
if ~converged
    rec = empty_record(numel(ctx.G0));
    rec.s = s;
    rec.x = x;
    rec.iterations = it;
    rec.map_residual = map_residual;
    return;
end

med = invz_emt_scalar(ctx.G0,Sigma,ctx.Jnu_flat,struct());
K = med.K(:);
K(1) = K0;
lambda = invz_lambdas(K,ctx.g,ctx.wts,ctx.beta,[1 2 3]);
sig = invz_sigma_ordered(ctx.tl,lambda,K,ctx.g,ctx.beta);
map_residual = max(abs(sig.Sigma-Sigma));
Sigma = sig.Sigma(:);
med = invz_emt_scalar(ctx.G0,Sigma,ctx.Jnu_flat,struct());
K = med.K(:);
K(1) = K0;
lambda = invz_lambdas(K,ctx.g,ctx.wts,ctx.beta,[1 2 3]);
sig = invz_sigma_ordered(ctx.tl,lambda,K,ctx.g,ctx.beta);
self_energy_residual = max(abs(sig.Sigma-Sigma));

[Gstat,go] = invz_gstat_ordered(ctx.tl,lambda(1:2),K0,Sigma(1), ...
    ctx.beta,ctx.G0inel0,ctx.G0el0,struct('stable_form',true));
static_residual = x-go.Gtil0;
raw_static_residual = Phi-Gstat;
scaled_raw = (1+K0*x)*raw_static_residual/(1-K0*Gstat);
Gq0 = Gstat./(1+(J-K0)*Gstat);
static_closure = mean(Gq0)-Gstat;

if opts.direct_dynamic_checks
    [dynamic_closure,dynamic_lattice_mass,dynamic_medium_mass] = ...
        dynamic_direct_checks(ctx,Sigma,K,med.G);
else
    dynamic_closure = NaN;
    dynamic_lattice_mass = NaN;
    dynamic_medium_mass = NaN;
end
rec = struct('s',s,'x',x,'converged',true,'iterations',it, ...
    'map_residual',map_residual,'self_energy_residual',self_energy_residual, ...
    'Sigma',Sigma,'K',K,'lambda',lambda,'G',med.G(:), ...
    'K0',K0,'Phi',Phi,'Gstat',Gstat,'Gtil0',go.Gtil0, ...
    'static_residual',static_residual, ...
    'raw_static_residual',raw_static_residual, ...
    'pole_cancel_identity_error',abs(static_residual-scaled_raw), ...
    'static_closure_residual',static_closure, ...
    'dynamic_closure_residual',dynamic_closure, ...
    'dynamic_lattice_mass',dynamic_lattice_mass, ...
    'dynamic_medium_mass',dynamic_medium_mass);
end

function [closure,minL,minD] = dynamic_direct_checks(ctx,Sigma,K,G)
nw = numel(G);
closure = 0;
minL = Inf;
minD = Inf;
retarded = ~isvector(ctx.Jnu_flat);
block = 64;
for i0 = 2:block:nw
    idx = i0:min(i0+block-1,nw);
    if retarded, J = ctx.Jnu_flat(:,idx); else, J = ctx.Jnu_flat(:); end
    D = 1+(J-K(idx).').*G(idx).';
    Gq = G(idx).'./D;
    closure = max(closure,max(abs(mean(Gq,1).'-G(idx))));
    L = 1+Sigma(idx).'+J.*ctx.G0(idx).';
    minL = min(minL,min(L,[],'all'));
    minD = min(minD,min(D,[],'all'));
end
end

function J = static_J(ctx)
if isvector(ctx.Jnu_flat), J = ctx.Jnu_flat(:);
else, J = ctx.Jnu_flat(:,1);
end
end

function rec = empty_record(nw)
rec = struct('s',NaN,'x',NaN,'converged',false,'iterations',NaN, ...
    'map_residual',NaN,'self_energy_residual',NaN, ...
    'Sigma',nan(nw,1),'K',nan(nw,1),'lambda',nan(3,1), ...
    'G',nan(nw,1),'K0',NaN,'Phi',NaN,'Gstat',NaN,'Gtil0',NaN, ...
    'static_residual',NaN,'raw_static_residual',NaN, ...
    'pole_cancel_identity_error',NaN,'static_closure_residual',NaN, ...
    'dynamic_closure_residual',NaN,'dynamic_lattice_mass',NaN, ...
    'dynamic_medium_mass',NaN);
end

function ctx = make_context(Bx,h,T,F)
ion = invz_ion();
[wn,wts,beta] = invz_matsubara(T,40);
si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,T,1i*wn(1),struct('elastic',false));
G0i = -real(c0i(3,3,1));
X = real(c0(:,:,1));
feedback = X(3,1)*(F.info.Jaa0/(1-F.info.Jaa0*X(1,1)))*X(1,3);
G0e = -(X(3,3)+feedback)-G0i;
g = real(invz_g(tl,1i*wn));
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',F.J,'J0eff',F.info.Jcc0, ...
    'G0inel0',G0i,'G0el0',G0e);
end
