function result = invzp_reduced_endpoint_profile()
%INVZP_REDUCED_ENDPOINT_PROFILE Solve moment closure along the h=0 static domain.
% For each fixed s=-J_sup*x, solve the three exact lambda moment equations
% after dynamic K_n elimination, then evaluate the pole-cancelled static
% residual x-Gtilde0. Forward and reverse continuation are compared.
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

s = [1e-5 0.02 0.05 0.1 0.2 0.35 0.5 0.65 0.8 0.9 0.95 0.98 ...
    0.999 0.99999].';
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-22,'StepTolerance',1e-13, ...
    'OptimalityTolerance',1e-12,'MaxIterations',80, ...
    'MaxFunctionEvaluations',600,'FiniteDifferenceType','forward');
lb = -ones(3,1)+1e-10;
ub = ones(3,1)-1e-10;

forward = run_pass(s,normalize_lambda(zeros(3,1)));
reverse = flipud(run_pass(flipud(s),forward(end).z));

lambda_delta = nan(size(s));
static_delta = nan(size(s));
for k = 1:numel(s)
    lambda_delta(k) = max(abs(forward(k).lambda-reverse(k).lambda));
    static_delta(k) = abs(forward(k).static_residual- ...
        reverse(k).static_residual);
end

result = struct('T',T,'Bx',Bx,'h',h,'Jsup',Jsup,'s',s, ...
    'lambda_lo',lambda_lo,'lambda_hi',lambda_hi, ...
    'forward',forward,'reverse',reverse, ...
    'forward_reverse_lambda_delta',lambda_delta, ...
    'forward_reverse_static_delta',static_delta, ...
    'note',['Moment closure is solved at fixed s in both directions. ' ...
    'This is a profile census, not an interval-complete root proof.']);
save(fullfile(here,'wp2_reduced_endpoint_profile.mat'),'result','-v7');

fprintf(['profile: static residual range [%.9g, %.9g], max lambda ' ...
    'closure %.3g, max forward/reverse delta %.3g\n'], ...
    min([forward.static_residual]),max([forward.static_residual]), ...
    max([forward.lambda_residual_norm]),max(lambda_delta));

    function records = run_pass(svals,z0)
        records = repmat(empty_record(),numel(svals),1);
        for j = 1:numel(svals)
            sj = svals(j);
            tic;
            [z,resnorm,rscaled,exitflag,output] = lsqnonlin( ...
                @(zz) lambda_vector(zz,sj),z0,lb,ub,lsopts);
            elapsed = toc;
            lambda = decode_lambda(z);
            q = invz_ordered_reduced_residual( ...
                [lambda;-sj/Jsup],ctx,ropts);
            records(j) = struct('s',sj,'z',z,'lambda',lambda, ...
                'resnorm',resnorm,'scaled_lambda_residual',rscaled(:), ...
                'lambda_residual_norm',max(abs(q.lambda_residual)), ...
                'static_residual',q.static_residual, ...
                'raw_static_residual',q.raw_static_residual, ...
                'Gstat',q.Gstat,'Gtil0',q.Gtil0, ...
                'local_denom',q.static.gstat_local_denom, ...
                'admissible',q.trial_admissible,'exitflag',exitflag, ...
                'output',output,'elapsed_s',elapsed, ...
                'dynamic_max_derivative',max(q.dynamic.derivative(2:end)), ...
                'dynamic_max_bound',max(q.dynamic.derivative_bound(2:end)));
            fprintf(['s %.6g: exit %d, ||Rlam||inf %.3g, ' ...
                'Rstatic %.9g, admissible %d, %.3fs\n'], ...
                sj,exitflag,records(j).lambda_residual_norm, ...
                q.static_residual,q.trial_admissible,elapsed);
            z0 = z;
        end
    end

    function r = lambda_vector(z,sj)
        q = invz_ordered_reduced_residual( ...
            [decode_lambda(z);-sj/Jsup],ctx,ropts);
        if q.defined && all(isfinite(q.lambda_residual))
            r = q.lambda_residual./lambda_half;
        else
            r = 1e4*(1+abs(z));
        end
    end

    function z = normalize_lambda(lambda)
        z = (lambda-lambda_mid)./lambda_half;
    end

    function lambda = decode_lambda(z)
        lambda = lambda_mid+lambda_half.*z;
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
rec = struct('s',NaN,'z',nan(3,1),'lambda',nan(3,1), ...
    'resnorm',NaN,'scaled_lambda_residual',nan(3,1), ...
    'lambda_residual_norm',NaN,'static_residual',NaN, ...
    'raw_static_residual',NaN,'Gstat',NaN,'Gtil0',NaN, ...
    'local_denom',NaN,'admissible',false,'exitflag',NaN, ...
    'output',struct(),'elapsed_s',NaN,'dynamic_max_derivative',NaN, ...
    'dynamic_max_bound',NaN);
end
