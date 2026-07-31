function result = invzp_reduced_fold_refinement()
%INVZP_REDUCED_FOLD_REFINEMENT Refine the 1 T continuation turning point.
% Diagnostic only.  Parameterize the certified component by
% s=-J_sup*x, solve the other four variables at fixed s, and refine the
% minimum of h(s) with a bounded local quadratic stencil.
% This is deliberately different from the bordered corrector used to trace
% the component.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
A = load(fullfile(here,'wp2_reduced_pseudoarclength_1t.mat'));
p = A.result.points;
T = A.result.T;
Bx = A.result.Bx;
Jsup = A.result.Jsup;
h_scale = A.result.h_scale;
lambda_lo = A.result.lambda_lo;
lambda_hi = A.result.lambda_hi;
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));

sdata = [p.s].';
udata = reshape([p.u],5,[]).';
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-18,'StepTolerance',1e-11, ...
    'OptimalityTolerance',1e-10,'MaxIterations',35, ...
    'MaxFunctionEvaluations',250,'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',5e-6);

% A three-point local stencil gives a bounded alternative to a nested scalar
% minimizer.  Refine once at the fitted vertex, then test stationarity there.
sgrid = [0.6515 0.6555 0.6595].';
grid = repmat(struct('s',NaN,'solution',[]),numel(sgrid),1);
hgrid = nan(size(sgrid));
for k = 1:numel(sgrid)
    [grid(k).solution,ok] = solve_fixed_s(sgrid(k));
    assert(ok);
    grid(k).s = sgrid(k);
    hgrid(k) = grid(k).solution.h;
end
quadratic = polyfit(sgrid,hgrid,2);
sfit = -quadratic(2)/(2*quadratic(1));
assert(sfit > min(sgrid) && sfit < max(sgrid));
[fit_solution,ok] = solve_fixed_s(sfit);
assert(ok);

d0 = 2e-4;
[fit_minus,okm] = solve_fixed_s(sfit-d0);
[fit_plus,okp] = solve_fixed_s(sfit+d0);
assert(okm && okp);
fit_d1 = (fit_plus.h-fit_minus.h)/(2*d0);
fit_d2 = (fit_plus.h-2*fit_solution.h+fit_minus.h)/(d0*d0);
sstar = sfit-fit_d1/fit_d2;
assert(sstar > min(sgrid) && sstar < max(sgrid) && fit_d2 > 0);
[fold,ok] = solve_fixed_s(sstar);
assert(ok);

deltas = 1e-4;
stencil = repmat(struct('delta',NaN,'hminus',NaN,'hplus',NaN, ...
    'first_derivative',NaN,'second_derivative',NaN),numel(deltas),1);
for k = 1:numel(deltas)
    d = deltas(k);
    [qm,okm] = solve_fixed_s(sstar-d);
    [qp,okp] = solve_fixed_s(sstar+d);
    assert(okm && okp);
    stencil(k) = struct('delta',d,'hminus',qm.h,'hplus',qp.h, ...
        'first_derivative',(qp.h-qm.h)/(2*d), ...
        'second_derivative',(qp.h-2*fold.h+qm.h)/(d*d));
end

fd_steps = [1e-5 5e-6];
jcheck = repmat(struct('step',NaN,'singular_values',nan(4,1), ...
    'svmin',NaN,'rcond',NaN,'determinant',NaN),numel(fd_steps),1);
for k = 1:numel(fd_steps)
    J = fixed_h_jacobian(fold.u,fd_steps(k));
    sv = svd(J);
    jcheck(k) = struct('step',fd_steps(k),'singular_values',sv, ...
        'svmin',min(sv),'rcond',rcond(J),'determinant',det(J));
end

result = struct('T',T,'Bx',Bx,'Jsup',Jsup, ...
    'source',fullfile(here,'wp2_reduced_pseudoarclength_1t.mat'), ...
    'search_interval',[min(sgrid) max(sgrid)],'grid',grid, ...
    'quadratic_coefficients',quadratic, ...
    'fitted_vertex',struct('s',sfit,'solution',fit_solution, ...
        'delta',d0,'hminus',fit_minus.h,'hplus',fit_plus.h, ...
        'first_derivative',fit_d1,'second_derivative',fit_d2), ...
    's',sstar,'h',fold.h, ...
    'fold',fold,'stencil',stencil, ...
    'fixed_h_jacobian_checks',jcheck, ...
    'note',['The fold is refined by fixed-s solves plus a local quadratic ' ...
    'fit, not by the bordered pseudo-arclength corrector. Central finite-difference ' ...
    'stencils test stationarity and singularity independently at three steps.']);
save(fullfile(here,'wp2_reduced_fold_refinement.mat'),'result','-v7');
fprintf(['fold refinement: s %.12g, h %.12g, ||R|| %.3g, ' ...
    'dh/ds %.3g, d2h/ds2 %.3g, svmin %.3g\n'], ...
    sstar,fold.h,fold.q.residual_norm,stencil(end).first_derivative, ...
    stencil(end).second_derivative,jcheck(end).svmin);

    function [sol,valid] = solve_fixed_s(s)
        uinit = nan(5,1);
        for jj = 1:5
            uinit(jj) = interp1(sdata,udata(:,jj),s,'pchip');
        end
        z0 = [uinit(1:3);uinit(5)];
        fun = @(z) scaled_fixed_s(z,s);
        [z,~,f,ef,out] = lsqnonlin(fun,z0, ...
            [-2*ones(3,1);0],[2*ones(3,1);2],lsopts);
        u = [z(1:3);s;z(4)];
        q = evaluate_q(u);
        valid = ef > 0 && q.defined && q.trial_admissible && ...
            max(abs(f)) <= 2e-8;
        sol = struct('u',u,'y',decode_y(u),'h',h_scale*u(5), ...
            'q',q,'scaled_residual',f(:),'exitflag',ef,'output',out);
    end

    function f = scaled_fixed_s(z,s)
        u = [z(1:3);s;z(4)];
        q = evaluate_q(u);
        if q.defined && all(isfinite(q.residual))
            f = [q.lambda_residual./lambda_half; ...
                q.static_residual/x_scale];
        else
            f = 100*(1+abs(z));
        end
    end

    function J = fixed_h_jacobian(u,step)
        f0 = scaled_full(u);
        J = nan(4,4);
        for jj = 1:4
            du = zeros(5,1);
            du(jj) = step*max(1,abs(u(jj)));
            fp = scaled_full(u+du);
            fm = scaled_full(u-du);
            J(:,jj) = (fp-fm)/(2*du(jj));
        end
        assert(all(isfinite(J),'all') && all(isfinite(f0)));
    end

    function f = scaled_full(u)
        q = evaluate_q(u);
        if q.defined && all(isfinite(q.residual))
            f = [q.lambda_residual./lambda_half; ...
                q.static_residual/x_scale];
        else
            f = nan(4,1);
        end
    end

    function q = evaluate_q(u)
        h = h_scale*u(5);
        if ~(isfinite(h) && h >= 0)
            q = struct('defined',false,'trial_admissible',false, ...
                'residual',nan(4,1));
            return;
        end
        try
            ctx = make_context(Bx,h,T,F);
            q = invz_ordered_reduced_residual(decode_y(u),ctx,ropts);
        catch
            q = struct('defined',false,'trial_admissible',false, ...
                'residual',nan(4,1));
        end
    end

    function y = decode_y(u)
        y = [lambda_mid+lambda_half.*u(1:3);-u(4)/Jsup];
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
