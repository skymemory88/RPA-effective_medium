function result = invzp_reduced_representation_profile()
%INVZP_REDUCED_REPRESENTATION_PROFILE Exact reduced h=0 hybrid/closed control.
% Reuse the saved moment-closed hybrid profile and probe an internally closed
% electronic two-level response/vertex pair along a bounded positive-lambda
% feasibility ray. The ray census exposes the dynamic-domain boundary and
% residual direction; it is not a global root-exclusion proof.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

hybrid_path = fullfile(repo,'docs','diagnostics','invzp_outer_wp2', ...
    'wp2_reduced_endpoint_profile.mat');
fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
H = load(hybrid_path);
F = load(fixture_path);
T = H.result.T;
Bx = H.result.Bx;
h = H.result.h;
Jsup = H.result.Jsup;
ctx = make_closed_context(Bx,h,T,F);
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));
[lambda_lo,lambda_hi] = lambda_bounds(ctx);

s_values = [1e-5 0.25 0.5 0.75 0.99999].';
ray_values = [0 0.5 0.6 0.65 0.675 0.7 0.7125 0.725 0.75 0.8 0.9 1].';
nrow = numel(s_values)*numel(ray_values);
s = nan(nrow,1);
ray = nan(nrow,1);
status = strings(nrow,1);
defined = false(nrow,1);
admissible = false(nrow,1);
lambda1_residual = nan(nrow,1);
lambda2_residual = nan(nrow,1);
lambda3_residual = nan(nrow,1);
static_residual = nan(nrow,1);
dynamic_max_derivative = nan(nrow,1);
dynamic_max_bound = nan(nrow,1);
dynamic_enumerated_frequencies = nan(nrow,1);
dynamic_unresolved_frequencies = nan(nrow,1);
first_failed_frequency = nan(nrow,1);

row = 0;
for is = 1:numel(s_values)
    for ia = 1:numel(ray_values)
        row = row+1;
        s(row) = s_values(is);
        ray(row) = ray_values(ia);
        lambda = ray(row)*lambda_hi;
        q = invz_ordered_reduced_residual( ...
            [lambda;-s(row)/Jsup],ctx,ropts);
        status(row) = string(q.status);
        defined(row) = q.defined;
        admissible(row) = q.trial_admissible;
        if ~isempty(q.dynamic)
            dynamic_enumerated_frequencies(row) = ...
                sum(q.dynamic.proof(2:end) == "finite_enumeration");
            dynamic_unresolved_frequencies(row) = ...
                sum(q.dynamic.fallback_unresolved(2:end));
            jf = find(q.dynamic.root_count(2:end) ~= 1,1);
            if ~isempty(jf), first_failed_frequency(row) = jf+1; end
            if q.dynamic.defined
                dynamic_max_derivative(row) = ...
                    max(q.dynamic.derivative(2:end));
                dynamic_max_bound(row) = ...
                    max(q.dynamic.derivative_bound(2:end));
            end
        end
        if q.defined
            lambda1_residual(row) = q.lambda_residual(1);
            lambda2_residual(row) = q.lambda_residual(2);
            lambda3_residual(row) = q.lambda_residual(3);
            static_residual(row) = q.static_residual;
        end
    end
    rows = find(s == s_values(is));
    first_defined = rows(find(defined(rows),1));
    if isempty(first_defined)
        fprintf('closed s %.6g: no defined point on feasibility ray\n', ...
            s_values(is));
    else
        fprintf(['closed s %.6g: first defined ray %.6g, ' ...
            'Rlambda=[%.6g %.6g %.6g], Rstatic %.6g\n'], ...
            s_values(is),ray(first_defined), ...
            lambda1_residual(first_defined), ...
            lambda2_residual(first_defined), ...
            lambda3_residual(first_defined), ...
            static_residual(first_defined));
    end
end

closed_ray = table(s,ray,status,defined,admissible, ...
    lambda1_residual,lambda2_residual,lambda3_residual,static_residual, ...
    dynamic_max_derivative,dynamic_max_bound, ...
    dynamic_enumerated_frequencies,dynamic_unresolved_frequencies, ...
    first_failed_frequency);
hybrid = table(H.result.s,[H.result.forward.static_residual].', ...
    [H.result.forward.lambda_residual_norm].', ...
    [H.result.forward.dynamic_max_bound].', ...
    'VariableNames',{'s','static_residual','lambda_residual_norm', ...
    'dynamic_max_bound'});
result = struct('T',T,'Bx',Bx,'h',h,'Jsup',Jsup, ...
    'lambda_lo',lambda_lo,'lambda_hi',lambda_hi, ...
    'hybrid_moment_closed_profile',hybrid,'closed_feasibility_ray',closed_ray, ...
    'sources',struct('hybrid_profile',hybrid_path,'fixture',fixture_path), ...
    'note',['Closed ray uses lambda=a*lambda_hi, where lambda_hi is the ' ...
    'rigorous componentwise moment bound and 0<=a<=1. It diagnoses the ' ...
    'dynamic-domain boundary and local residual direction, but does not ' ...
    'exclude off-ray roots.']);
save(fullfile(here,'wp3_reduced_representation_profile.mat'),'result','-v7');
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

function ctx = make_closed_context(Bx,h,T,F)
ion = invz_ion();
[wn,wts,beta] = invz_matsubara(T,40);
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
g = real(invz_g(tl,1i*wn));
G0 = -tl.M2*g;
G0(1) = G0(1)-tl.m^2*tl.h0;
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',F.J,'J0eff',F.info.Jcc0, ...
    'G0inel0',-tl.M2*tl.g0,'G0el0',-tl.m^2*tl.h0);
end
