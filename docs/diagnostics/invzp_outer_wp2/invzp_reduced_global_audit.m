function result = invzp_reduced_global_audit()
%INVZP_REDUCED_GLOBAL_AUDIT Independent global search of the exact 4D residual.
% Uses surrogateopt followed by patternsearch in normalized bounded variables.
% This is algorithmically independent of the prior least-squares multistarts;
% a finite evaluation budget remains evidence, not a root-completeness proof.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'),'T');
N = load(fullfile(here,'wp2_node22_step_domain.mat'));
T = D.T;
Bx = 1;
cases = struct('name',{"endpoint_h0","node22"}, ...
    'h',{0,N.result.h});
records = cell(numel(cases),1);

for ic = 1:numel(cases)
    ctx = make_context(Bx,cases(ic).h,T,F);
    records{ic} = search_case(cases(ic).name,ctx,F.info.Jcc0);
end

result = struct('T',T,'Bx',Bx,'cases',cases,'records',{records}, ...
    'note',['Global surrogate search plus derivative-free pattern refinement ' ...
    'in the rigorous lambda box and strict physical x interval. Finite ' ...
    'budgets do not prove root absence.']);
save(fullfile(here,'wp2_reduced_global_audit.mat'),'result','-v7');
end

function rec = search_case(name,ctx,Jsup)
[lambda_lo,lambda_hi] = lambda_bounds(ctx);
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));
lb = [-ones(3,1)+1e-8;1e-6];
ub = [ ones(3,1)-1e-8;1-1e-6];

sopts = optimoptions('surrogateopt','Display','iter', ...
    'PlotFcn',[],'MaxFunctionEvaluations',100, ...
    'MinSurrogatePoints',20,'BatchUpdateInterval',1);
tic;
[zglobal,fglobal,exitglobal,outglobal,trials] = surrogateopt( ...
    @objective,lb,ub,sopts);
global_elapsed = toc;

popts = optimoptions('patternsearch','Display','final', ...
    'MaxIterations',80,'MaxFunctionEvaluations',160, ...
    'FunctionTolerance',1e-12,'MeshTolerance',1e-8, ...
    'UseCompletePoll',true,'UseCompleteSearch',false);
tic;
[zpattern,fpattern,exitpattern,outpattern] = patternsearch( ...
    @objective,zglobal,[],[],[],[],lb,ub,[],popts);
pattern_elapsed = toc;

qglobal = invz_ordered_reduced_residual(decode(zglobal),ctx,ropts);
qpattern = invz_ordered_reduced_residual(decode(zpattern),ctx,ropts);
rec = struct('name',name,'lambda_lo',lambda_lo,'lambda_hi',lambda_hi, ...
    'global_solution',zglobal,'global_y',decode(zglobal), ...
    'global_objective',fglobal,'global_exitflag',exitglobal, ...
    'global_output',outglobal,'global_trials',trials, ...
    'global_elapsed_s',global_elapsed, ...
    'global_status',string(qglobal.status), ...
    'global_admissible',qglobal.trial_admissible, ...
    'global_residual',qglobal.residual, ...
    'pattern_solution',zpattern,'pattern_y',decode(zpattern), ...
    'pattern_objective',fpattern,'pattern_exitflag',exitpattern, ...
    'pattern_output',outpattern,'pattern_elapsed_s',pattern_elapsed, ...
    'pattern_status',string(qpattern.status), ...
    'pattern_admissible',qpattern.trial_admissible, ...
    'pattern_residual',qpattern.residual, ...
    'pattern_scaled_residual',scaled_residual(qpattern));
fprintf(['%s global audit: surrogate f %.9g, pattern f %.9g, ' ...
    '||Rscaled||inf %.9g, admissible %d\n'],name,fglobal,fpattern, ...
    max(abs(rec.pattern_scaled_residual)),qpattern.trial_admissible);

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
            if ~q.trial_admissible
                % Keep the exact residual visible while weakly preferring the
                % physical trial domain. Roots are checked without this term.
                f = f+1e-6;
            end
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
