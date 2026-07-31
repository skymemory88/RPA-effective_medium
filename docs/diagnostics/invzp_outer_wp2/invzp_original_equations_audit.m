function result = invzp_original_equations_audit()
%INVZP_ORIGINAL_EQUATIONS_AUDIT Independent check of the reduced conclusions.
% This diagnostic does not use INVZ_ORDERED_DYNAMIC_ELIMINATE to construct its
% branches. At fixed static x it iterates the original full Sigma vector:
%   Sigma -> invz_emt_scalar -> K -> lambda -> invz_sigma_ordered.
% It then evaluates the original EMT closures by direct substitution.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
H = load(fullfile(here,'wp2_reduced_endpoint_profile.mat'));
T = H.result.T;
Bx = 1;
Jsup = F.info.Jcc0;
s_dense = sort(unique([H.result.s(:);linspace(1e-6,1-1e-6,101).'; ...
    1-10.^(-(2:7)).']));
opts = struct('mix',0.5,'tol',1e-12,'maxit',250, ...
    'direct_dynamic_checks',false);

ctx0 = make_context(Bx,0,T,F);
forward0 = direct_profile(ctx0,s_dense,Jsup,zeros(size(ctx0.G0)),opts);
reverse0 = flip_profile(direct_profile(ctx0,flipud(s_dense),Jsup, ...
    zeros(size(ctx0.G0)),opts));
fprintf('completed original-equation h0 forward/reverse profiles\n');

assert(all(forward0.converged) && all(reverse0.converged));
forward_reverse_sigma = max(abs(forward0.Sigma-reverse0.Sigma),[],'all');
forward_reverse_lambda = max(abs(forward0.lambda-reverse0.lambda),[],'all');
forward_reverse_static = max(abs( ...
    forward0.static_residual-reverse0.static_residual));

[tf,idx] = ismember(H.result.s,s_dense);
assert(all(tf));
saved_lambda = reshape([H.result.forward.lambda],3,[]);
saved_static = [H.result.forward.static_residual];
reduced_profile_lambda_delta = max(abs( ...
    forward0.lambda(:,idx)-saved_lambda),[],'all');
reduced_profile_static_delta = max(abs( ...
    forward0.static_residual(idx).'-saved_static),[],'all');

imid = find(abs(s_dense-0.5) == min(abs(s_dense-0.5)),1);
seed_matrix = [zeros(numel(ctx0.G0),1),forward0.Sigma(:,imid), ...
    -forward0.Sigma(:,imid),2*forward0.Sigma(:,imid), ...
    0.25*ones(numel(ctx0.G0),1)];
multiseed0 = repmat(empty_record(numel(ctx0.G0)),size(seed_matrix,2),1);
for k = 1:size(seed_matrix,2)
    multiseed0(k) = direct_solve(ctx0,s_dense(imid),Jsup, ...
        seed_matrix(:,k),opts);
end
multiseed_sigma_delta = max(abs( ...
    reshape([multiseed0.Sigma],numel(ctx0.G0),[])- ...
    forward0.Sigma(:,imid)),[],'all');
checkopts = opts;
checkopts.direct_dynamic_checks = true;
audit_indices = unique([1 imid numel(s_dense)]);
direct_substitution0 = repmat(empty_record(numel(ctx0.G0)), ...
    numel(audit_indices),1);
for k = 1:numel(audit_indices)
    j = audit_indices(k);
    direct_substitution0(k) = direct_solve(ctx0,s_dense(j),Jsup, ...
        forward0.Sigma(:,j),checkopts);
end

h_accepted = 0.00624162310965;
ctxa = make_context(Bx,h_accepted,T,F);
s_accept = linspace(0.2,0.5,81).';
profilea = direct_profile(ctxa,s_accept,Jsup,zeros(size(ctxa.G0)),opts);
assert(all(profilea.converged));
brackets = find(profilea.static_residual(1:end-1).* ...
    profilea.static_residual(2:end) < 0);
roots = repmat(empty_record(numel(ctxa.G0)),numel(brackets),1);
for k = 1:numel(brackets)
    roots(k) = bisect_static(ctxa,Jsup,s_accept(brackets(k)), ...
        s_accept(brackets(k)+1),profilea.Sigma(:,brackets(k)),opts);
end
direct_substitution_accepted = repmat(empty_record(numel(ctxa.G0)), ...
    numel(roots),1);
for k = 1:numel(roots)
    direct_substitution_accepted(k) = direct_solve(ctxa, ...
        -Jsup*roots(k).x,Jsup,roots(k).Sigma,checkopts);
end

root_comparison = repmat(struct('reduced_status',"", ...
    'reduced_residual_norm',NaN,'K_delta',NaN,'Sigma_delta',NaN, ...
    'outer_status',"",'outer_residual_norm',NaN),numel(roots),1);
for k = 1:numel(roots)
    q = invz_ordered_reduced_residual( ...
        [roots(k).lambda;roots(k).x],ctxa,struct('Jsup',Jsup));
    outer = invz_ordered_outer_map(roots(k).Sigma,ctxa, ...
        struct('emt_static',struct('Jsup',Jsup,'warn',false)));
    root_comparison(k) = struct( ...
        'reduced_status',string(q.status), ...
        'reduced_residual_norm',q.residual_norm, ...
        'K_delta',max(abs(q.K-roots(k).K)), ...
        'Sigma_delta',max(abs(q.Sigma-roots(k).Sigma)), ...
        'outer_status',string(outer.status), ...
        'outer_residual_norm',outer.residual_norm);
end

result = struct('T',T,'Bx',Bx,'Jsup',Jsup, ...
    'h0',struct('s',s_dense,'forward',forward0,'reverse',reverse0, ...
        'forward_reverse_sigma_delta',forward_reverse_sigma, ...
        'forward_reverse_lambda_delta',forward_reverse_lambda, ...
        'forward_reverse_static_delta',forward_reverse_static, ...
        'reduced_profile_lambda_delta',reduced_profile_lambda_delta, ...
        'reduced_profile_static_delta',reduced_profile_static_delta, ...
        'multiseed',multiseed0, ...
        'multiseed_sigma_delta',multiseed_sigma_delta, ...
        'direct_substitution',direct_substitution0), ...
    'accepted',struct('h',h_accepted,'profile',profilea, ...
        'brackets',brackets,'roots',roots, ...
        'root_comparison',root_comparison, ...
        'direct_substitution',direct_substitution_accepted), ...
    'note',['Branches are constructed by full-Sigma Picard using the original ' ...
    'direct EMT medium, independently of the scalar dynamic eliminator and ' ...
    'least-squares moment solver. Dense scans remain numerical evidence, ' ...
    'not interval proofs.']);
save(fullfile(here,'wp2_original_equations_audit.mat'),'result','-v7');

fprintf(['original-equation audit: h0 Rstatic [%.9g, %.9g], ' ...
    'forward/reverse dSigma %.3g, saved-profile dLambda %.3g; ' ...
    'accepted brackets %d\n'],min(forward0.static_residual), ...
    max(forward0.static_residual),forward_reverse_sigma, ...
    reduced_profile_lambda_delta,numel(brackets));
end

function profile = direct_profile(ctx,svals,Jsup,Sigma,opts)
nw = numel(ctx.G0);
ns = numel(svals);
records = repmat(empty_record(nw),ns,1);
for k = 1:ns
    records(k) = direct_solve(ctx,svals(k),Jsup,Sigma,opts);
    if records(k).converged, Sigma = records(k).Sigma; end
end
profile = records_to_profile(records,svals,nw);
[profile.dynamic_lattice_mass,profile.dynamic_medium_mass] = ...
    fast_dynamic_masses(ctx,profile);
end

function profile = flip_profile(profile)
fields = fieldnames(profile);
for k = 1:numel(fields)
    v = profile.(fields{k});
    if isvector(v)
        profile.(fields{k}) = flip(v);
    else
        profile.(fields{k}) = fliplr(v);
    end
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

% Re-evaluate every original equation at the returned Sigma.
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

function root = bisect_static(ctx,Jsup,a,b,Sigma,opts)
ra = direct_solve(ctx,a,Jsup,Sigma,opts);
rb = direct_solve(ctx,b,Jsup,ra.Sigma,opts);
assert(ra.converged && rb.converged && ...
    sign(ra.static_residual) ~= sign(rb.static_residual));
root = ra;
for it = 1:100
    m = a+(b-a)/2;
    rm = direct_solve(ctx,m,Jsup,root.Sigma,opts);
    assert(rm.converged);
    if abs(rm.static_residual) < abs(root.static_residual), root = rm; end
    if sign(ra.static_residual) ~= sign(rm.static_residual)
        b = m;
        rb = rm; %#ok<NASGU>
    else
        a = m;
        ra = rm;
    end
    if abs(root.static_residual) <= 1e-11 || b-a <= 1e-13, break; end
end
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

function profile = records_to_profile(records,svals,nw)
profile = struct('s',svals(:),'converged',[records.converged].', ...
    'iterations',[records.iterations].','map_residual',[records.map_residual].', ...
    'self_energy_residual',[records.self_energy_residual].', ...
    'Sigma',reshape([records.Sigma],nw,[]), ...
    'K',reshape([records.K],nw,[]), ...
    'G',reshape([records.G],nw,[]), ...
    'lambda',reshape([records.lambda],3,[]), ...
    'static_residual',[records.static_residual].', ...
    'raw_static_residual',[records.raw_static_residual].', ...
    'pole_cancel_identity_error',[records.pole_cancel_identity_error].', ...
    'static_closure_residual',[records.static_closure_residual].', ...
    'dynamic_closure_residual',[records.dynamic_closure_residual].', ...
    'dynamic_lattice_mass',[records.dynamic_lattice_mass].', ...
    'dynamic_medium_mass',[records.dynamic_medium_mass].');
end

function [minL,minD] = fast_dynamic_masses(ctx,profile)
nw = numel(ctx.G0);
if isvector(ctx.Jnu_flat)
    jlo = repmat(min(ctx.Jnu_flat),nw,1);
    jhi = repmat(max(ctx.Jnu_flat),nw,1);
else
    jlo = min(ctx.Jnu_flat,[],1).';
    jhi = max(ctx.Jnu_flat,[],1).';
end
jL = jlo;
jL(ctx.G0 < 0) = jhi(ctx.G0 < 0);
L = 1+profile.Sigma+jL.*ctx.G0;
jD = repmat(jlo,1,size(profile.G,2));
jhi_matrix = repmat(jhi,1,size(profile.G,2));
jD(profile.G < 0) = jhi_matrix(profile.G < 0);
D = 1+(jD-profile.K).*profile.G;
minL = min(L(2:end,:),[],1).';
minD = min(D(2:end,:),[],1).';
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
