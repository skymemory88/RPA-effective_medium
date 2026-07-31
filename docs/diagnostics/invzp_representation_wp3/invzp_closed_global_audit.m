function result = invzp_closed_global_audit()
%INVZP_CLOSED_GLOBAL_AUDIT Off-ray global root search in the closed 1 T h=0 model.
% Uses the exact reduced residual with an internally matched electronic
% two-level response/vertex pair. Surrogate and pattern searches cover the
% rigorous lambda box and strict static domain. This remains budgeted evidence.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
T = 0.1;
Bx = 1;
h = 0;
Jsup = F.info.Jcc0;
ctx = make_closed_context(Bx,h,T,F);
[lambda_lo,lambda_hi] = lambda_bounds(ctx);
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));
lb = [-ones(3,1)+1e-8;1e-6];
ub = [ ones(3,1)-1e-8;1-1e-6];

initial = zeros(6,4);
row = 0;
for a = [0.75 1]
    for s = [0.25 0.5 0.999]
        row = row+1;
        initial(row,:) = [((a*lambda_hi-lambda_mid)./lambda_half).' s];
    end
end
sopts = optimoptions('surrogateopt','Display','iter','PlotFcn',[], ...
    'MaxFunctionEvaluations',120,'MinSurrogatePoints',20, ...
    'BatchUpdateInterval',1,'InitialPoints',initial);
tic;
[zglobal,fglobal,exitglobal,outglobal,trials] = surrogateopt( ...
    @objective,lb,ub,sopts);
global_elapsed = toc;

popts = optimoptions('patternsearch','Display','final', ...
    'MaxIterations',100,'MaxFunctionEvaluations',200, ...
    'FunctionTolerance',1e-12,'MeshTolerance',1e-8, ...
    'UseCompletePoll',true,'UseCompleteSearch',false);
tic;
[zpattern,fpattern,exitpattern,outpattern] = patternsearch( ...
    @objective,zglobal,[],[],[],[],lb,ub,[],popts);
pattern_elapsed = toc;

qglobal = invz_ordered_reduced_residual(decode(zglobal),ctx,ropts);
qpattern = invz_ordered_reduced_residual(decode(zpattern),ctx,ropts);
result = struct('T',T,'Bx',Bx,'h',h,'Jsup',Jsup, ...
    'lambda_lo',lambda_lo,'lambda_hi',lambda_hi,'initial_points',initial, ...
    'global_solution',zglobal,'global_y',decode(zglobal), ...
    'global_objective',fglobal,'global_exitflag',exitglobal, ...
    'global_output',outglobal,'global_trials',trials, ...
    'global_elapsed_s',global_elapsed,'global_status',string(qglobal.status), ...
    'global_admissible',qglobal.trial_admissible, ...
    'global_residual',qglobal.residual, ...
    'pattern_solution',zpattern,'pattern_y',decode(zpattern), ...
    'pattern_objective',fpattern,'pattern_exitflag',exitpattern, ...
    'pattern_output',outpattern,'pattern_elapsed_s',pattern_elapsed, ...
    'pattern_status',string(qpattern.status), ...
    'pattern_admissible',qpattern.trial_admissible, ...
    'pattern_residual',qpattern.residual, ...
    'pattern_scaled_residual',scaled_residual(qpattern), ...
    'note',['Off-ray global search with six feasible-ray initial points. ' ...
    'Finite budgets do not prove root absence.']);
save(fullfile(here,'wp3_closed_global_audit.mat'),'result','-v7');
fprintf(['closed global audit: surrogate f %.9g, pattern f %.9g, ' ...
    '||Rscaled||inf %.9g, admissible %d\n'],fglobal,fpattern, ...
    max(abs(result.pattern_scaled_residual)),qpattern.trial_admissible);

    function y = decode(z)
        z = z(:);
        y = [lambda_mid+lambda_half.*z(1:3);-z(4)/Jsup];
    end

    function r = scaled_residual(q)
        r = [q.lambda_residual./lambda_half; ...
            q.static_residual/x_scale];
    end

    function f = objective(z)
        q = invz_ordered_reduced_residual(decode(z),ctx,ropts);
        if q.defined && all(isfinite(q.residual))
            r = scaled_residual(q);
            f = sum(r.^2);
            if ~q.trial_admissible, f = f+1e-6; end
        else
            f = 1e3+sum(z(:).^2);
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
