function result = invzp_reduced_endpoint_roots()
%INVZP_REDUCED_ENDPOINT_ROOTS Diagnostic multi-start solve of the exact 4D residual.
% This is a root census, not a completeness proof. Dynamic K_n variables are
% eliminated exactly by INVZ_ORDERED_DYNAMIC_ELIMINATE. The four retained
% variables are lambda(1:3) and s=-J_sup*x in the strict interval (0,1).
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'),'T');
T = D.T;
Bx = 1;
h = 0;
ctx = make_context(Bx,h,T,F);
Jsup = F.info.Jcc0;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));

[lambda_lo,lambda_hi] = lambda_bounds(ctx);
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;

% Each s gets a zero-lambda start and a one-step moment-consistent start.
s_seed = [0.1 0.35 0.65 0.9];
z_seed = nan(2*numel(s_seed),4);
for j = 1:numel(s_seed)
    q0 = invz_ordered_reduced_residual( ...
        [zeros(3,1);-s_seed(j)/Jsup],ctx,ropts);
    lambda_seed = min(max(q0.dynamic.lambda_check,lambda_lo),lambda_hi);
    z_seed(2*j-1,:) = [normalize_lambda(zeros(3,1)).' s_seed(j)];
    z_seed(2*j,:) = [normalize_lambda(lambda_seed).' s_seed(j)];
end

lb = [-ones(3,1)+1e-9;1e-7];
ub = [ ones(3,1)-1e-9;1-1e-7];
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-20,'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-11,'MaxIterations',120, ...
    'MaxFunctionEvaluations',1200,'FiniteDifferenceType','forward');

nstart = size(z_seed,1);
records = repmat(empty_record(),nstart,1);
for k = 1:nstart
    z0 = z_seed(k,:).';
    tic;
    [z,resnorm,scaled_residual,exitflag,output] = lsqnonlin( ...
        @scaled_vector,z0,lb,ub,lsopts);
    elapsed = toc;
    y = decode(z);
    q = invz_ordered_reduced_residual(y,ctx,ropts);
    if isempty(q.static)
        local_denom = NaN;
    else
        local_denom = q.static.gstat_local_denom;
    end
    if isempty(q.dynamic)
        dynamic_max_derivative = NaN;
        dynamic_max_bound = NaN;
        dynamic_max_root_count = NaN;
    else
        dynamic_max_derivative = max(q.dynamic.derivative(2:end));
        dynamic_max_bound = max(q.dynamic.derivative_bound(2:end));
        dynamic_max_root_count = max(q.dynamic.root_count(2:end));
    end
    records(k) = struct('start',z0,'solution',z,'y',y, ...
        'resnorm',resnorm,'scaled_residual',scaled_residual(:), ...
        'scaled_residual_norm',max(abs(scaled_residual)), ...
        'exitflag',exitflag,'output',output,'elapsed_s',elapsed, ...
        'status',string(q.status),'defined',q.defined, ...
        'admissible',q.trial_admissible,'physical_residual',q.residual, ...
        'physical_residual_norm',q.residual_norm, ...
        'raw_static_residual',q.raw_static_residual, ...
        'local_denom',local_denom, ...
        'dynamic_max_derivative',dynamic_max_derivative, ...
        'dynamic_max_bound',dynamic_max_bound, ...
        'dynamic_max_root_count',dynamic_max_root_count);
    fprintf(['start %d/%d: exit %d, ||Rscaled||inf %.6g, ' ...
        'admissible %d, s %.9g, %.3fs\n'],k,nstart,exitflag, ...
        records(k).scaled_residual_norm,q.trial_admissible,z(4),elapsed);
end

Z = reshape([records.solution],4,[]).';
cluster = zeros(nstart,1);
ncluster = 0;
for k = 1:nstart
    prior = find(vecnorm(Z(1:k-1,:)-Z(k,:),Inf,2) < 2e-6,1);
    if isempty(prior)
        ncluster = ncluster+1;
        cluster(k) = ncluster;
    else
        cluster(k) = cluster(prior);
    end
end

result = struct('T',T,'Bx',Bx,'h',h,'Jsup',Jsup, ...
    'lambda_lo',lambda_lo,'lambda_hi',lambda_hi, ...
    'z_seed',z_seed,'records',records,'cluster',cluster, ...
    'ncluster',ncluster, ...
    'note',['Bounded multi-start least-squares census of the exact reduced ' ...
    'residual; absence of a root is not a completeness proof.']);
save(fullfile(here,'wp2_reduced_endpoint_roots.mat'),'result','-v7');

    function zlam = normalize_lambda(lambda)
        zlam = (lambda-lambda_mid)./lambda_half;
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
    'scaled_residual_norm',NaN,'exitflag',NaN,'output',struct(), ...
    'elapsed_s',NaN,'status',"",'defined',false,'admissible',false, ...
    'physical_residual',nan(4,1),'physical_residual_norm',NaN, ...
    'raw_static_residual',NaN,'local_denom',NaN, ...
    'dynamic_max_derivative',NaN,'dynamic_max_bound',NaN, ...
    'dynamic_max_root_count',NaN);
end
