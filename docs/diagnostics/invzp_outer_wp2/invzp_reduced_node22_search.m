function result = invzp_reduced_node22_search()
%INVZP_REDUCED_NODE22_SEARCH Bounded exact-residual search at noncontractive node 22.
% Diagnostic only. Uses the last admissible Picard state plus independent
% starts in normalized [lambda(1:3),s] coordinates.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
N = load(fullfile(here,'wp2_node22_step_domain.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'),'T');
T = D.T;
Bx = 1;
h = N.result.h;
ctx = make_context(Bx,h,T,F);
Jsup = F.info.Jcc0;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));
[lambda_lo,lambda_hi] = lambda_bounds(ctx);
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;

M = N.result.probe_mix05.last_admissible_map;
y_last = [M.lambda(:);M.static.Gtil0];
q_last = invz_ordered_reduced_residual(y_last,ctx,ropts);
lambda_moment = min(max(q_last.dynamic.lambda_check,lambda_lo),lambda_hi);
z_seed = [encode(y_last).'; ...
    encode([lambda_moment;y_last(4)]).'];

lb = [-ones(3,1)+1e-9;1e-7];
ub = [ ones(3,1)-1e-9;1-1e-7];
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-22,'StepTolerance',1e-13, ...
    'OptimalityTolerance',1e-12,'MaxIterations',30, ...
    'MaxFunctionEvaluations',180,'FiniteDifferenceType','forward');

nstart = size(z_seed,1);
records = repmat(empty_record(),nstart,1);
for k = 1:nstart
    tic;
    [z,resnorm,rscaled,exitflag,output] = lsqnonlin( ...
        @scaled_vector,z_seed(k,:).',lb,ub,lsopts);
    elapsed = toc;
    y = decode(z);
    q = invz_ordered_reduced_residual(y,ctx,ropts);
    records(k) = struct('start',z_seed(k,:).','solution',z,'y',y, ...
        'resnorm',resnorm,'scaled_residual',rscaled(:), ...
        'scaled_residual_norm',max(abs(rscaled)), ...
        'physical_residual',q.residual, ...
        'physical_residual_norm',q.residual_norm, ...
        'status',string(q.status),'defined',q.defined, ...
        'admissible',q.trial_admissible,'exitflag',exitflag, ...
        'output',output,'elapsed_s',elapsed, ...
        'supremum_mass',q.supremum_mass,'uniform_mass',q.D_uni, ...
        'static_mesh_mass',q.static_mesh_min_mass, ...
        'dynamic_lattice_mass',q.dynamic.dynamic_min_lattice_mass, ...
        'dynamic_medium_mass',q.dynamic.dynamic_min_medium_mass, ...
        'dynamic_max_derivative',max(q.dynamic.derivative(2:end)), ...
        'dynamic_max_bound',max(q.dynamic.derivative_bound(2:end)), ...
        'dynamic_unresolved_frequencies', ...
            sum(q.dynamic.fallback_unresolved(2:end)));
    fprintf(['node22 start %d/%d: exit %d, ||Rscaled||inf %.9g, ' ...
        'admissible %d, s %.9g, %.3fs\n'],k,nstart,exitflag, ...
        records(k).scaled_residual_norm,q.trial_admissible,z(4),elapsed);
end

result = struct('T',T,'Bx',Bx,'h',h,'Jsup',Jsup, ...
    'lambda_lo',lambda_lo,'lambda_hi',lambda_hi, ...
    'last_admissible_y',y_last,'last_admissible_reduced',q_last, ...
    'z_seed',z_seed,'records',records, ...
    'source',fullfile(here,'wp2_node22_step_domain.mat'), ...
    'note',['Bounded multi-start least-squares census of the exact reduced ' ...
    'residual at the previously noncontractive Picard fixture; not a ' ...
    'completeness proof.']);
save(fullfile(here,'wp2_reduced_node22_search.mat'),'result','-v7');

    function z = encode(y)
        z = [(y(1:3)-lambda_mid)./lambda_half;-Jsup*y(4)];
    end

    function y = decode(z)
        y = [lambda_mid+lambda_half.*z(1:3);-z(4)/Jsup];
    end

    function r = scaled_vector(z)
        q = invz_ordered_reduced_residual(decode(z),ctx,ropts);
        if q.defined && all(isfinite(q.residual))
            r = [q.lambda_residual./lambda_half; ...
                q.static_residual/x_scale];
        else
            r = 1e4*(1+abs(z));
        end
    end
end

function [lo,hi] = lambda_bounds(ctx)
nw = numel(ctx.g);
if isvector(ctx.Jnu_flat)
    jlo = repmat(min(ctx.Jnu_flat),nw,1);
    jhi = repmat(max(ctx.Jnu_flat),nw,1);
else
    jlo = min(ctx.Jnu_flat,[],1).';
    jhi = max(ctx.Jnu_flat,[],1).';
end
lo = nan(3,1);
hi = nan(3,1);
for p = 1:3
    a = ctx.wts(:).*ctx.g(:).^p/ctx.beta;
    lo(p) = sum(min(a.*jlo,a.*jhi));
    hi(p) = sum(max(a.*jlo,a.*jhi));
end
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

function rec = empty_record()
rec = struct('start',nan(4,1),'solution',nan(4,1),'y',nan(4,1), ...
    'resnorm',NaN,'scaled_residual',nan(4,1), ...
    'scaled_residual_norm',NaN,'physical_residual',nan(4,1), ...
    'physical_residual_norm',NaN,'status',"",'defined',false, ...
    'admissible',false,'exitflag',NaN,'output',struct(),'elapsed_s',NaN, ...
    'supremum_mass',NaN,'uniform_mass',NaN,'static_mesh_mass',NaN, ...
    'dynamic_lattice_mass',NaN,'dynamic_medium_mass',NaN, ...
    'dynamic_max_derivative',NaN,'dynamic_max_bound',NaN, ...
    'dynamic_unresolved_frequencies',NaN);
end
