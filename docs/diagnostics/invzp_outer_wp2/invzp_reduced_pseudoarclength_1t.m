function result = invzp_reduced_pseudoarclength_1t()
%INVZP_REDUCED_PSEUDOARCLENGTH_1T Bordered continuation of the exact 1 T root.
% Diagnostic only. Continues the verified four-variable reduced root in
% u=[normalized lambda(1:3);s=-J_sup*x;eta=h/h_scale], using a bordered
% pseudo-arclength corrector. Production dispatch and acceptance rules are
% unchanged.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
A = load(fullfile(here,'wp2_original_equations_audit.mat'));
T = A.result.T;
Bx = A.result.Bx;
Jsup = A.result.Jsup;
h0 = A.result.accepted.h;
root0 = A.result.accepted.roots(1);
ctx0 = make_context(Bx,h0,T,F);
[lambda_lo,lambda_hi] = lambda_bounds(ctx0);
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;
h_scale = 0.01;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));

u0 = encode([root0.lambda;root0.x],h0);
q0 = evaluate_q(u0);
assert(q0.defined && q0.trial_admissible && q0.residual_norm < 1e-9);
[J0,jstatus] = finite_jacobian(u0);
assert(jstatus == "ok");
[t0,sv0] = tangent_from_jacobian(J0,[]);
if t0(5) > 0, t0 = -t0; end

target_dh = 5e-5;
ds0 = min(0.03,max(0.001,(target_dh/h_scale)/max(abs(t0(5)),0.1)));
lb = [-2*ones(3,1);1e-8;0];
ub = [ 2*ones(3,1);1-1e-8;2];
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-20,'StepTolerance',1e-11, ...
    'OptimalityTolerance',1e-10,'MaxIterations',35, ...
    'MaxFunctionEvaluations',300,'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',1e-5);

% Reversibility gate: one small downward correction, then reverse orientation
% and return to the starting section.
[down,down_ok] = correct_step(u0,t0,ds0);
assert(down_ok);
[tdown,svdown] = tangent_from_jacobian(down.J,t0);
[back,back_ok] = correct_step(down.u,-tdown,ds0);
assert(back_ok);
reversibility_delta = max(abs(back.u-u0));

max_points = 50;
points = repmat(empty_point(),max_points,1);
points(1) = point_record(u0,q0,J0,t0,sv0,0,0);
qdown = evaluate_q(down.u);
points(2) = point_record(down.u,qdown,down.J,tdown,svdown, ...
    down.output.iterations,ds0);
npoint = 2;
step = ds0;
attempts = repmat(empty_attempt(),0,1);
stop_reason = "max_points";

for ip = 3:max_points
    accepted = false;
    for retry = 1:8
        [candidate,ok] = correct_step(points(ip-1).u, ...
            points(ip-1).tangent,step);
        attempts(end+1) = attempt_record(ip,retry,step,candidate,ok); %#ok<AGROW>
        if ok
            accepted = true;
            break;
        end
        step = step/2;
        if step < 2e-5, break; end
    end
    if ~accepted
        stop_reason = "corrector_failed";
        break;
    end

    q = evaluate_q(candidate.u);
    [tnew,sv] = tangent_from_jacobian(candidate.J, ...
        points(ip-1).tangent);
    npoint = ip;
    points(ip) = point_record(candidate.u,q,candidate.J,tnew,sv, ...
        candidate.output.iterations,step);
    fprintf(['arc %d: h %.12g, s %.9g, ||R|| %.3g, ' ...
        'svmin %.3g, masses %.3g/%.3g, step %.3g\n'], ...
        ip,points(ip).h,points(ip).s,points(ip).residual_norm, ...
        points(ip).unbordered_svmin,points(ip).supremum_mass, ...
        points(ip).uniform_mass,step);

    if points(ip).h <= 1e-10
        stop_reason = "reached_h0";
        break;
    end
    if min([points(ip).supremum_mass,points(ip).uniform_mass, ...
            points(ip).static_mesh_mass,points(ip).dynamic_lattice_mass, ...
            points(ip).dynamic_medium_mass]) <= 1e-7
        stop_reason = "physical_mass_boundary";
        break;
    end
    if candidate.output.iterations <= 5
        step = min(0.04,1.25*step);
    elseif candidate.output.iterations >= 12
        step = max(2e-5,0.7*step);
    end
end
points = points(1:npoint);

if stop_reason == "corrector_failed"
    last = points(end);
    if min([last.supremum_mass,last.uniform_mass,last.static_mesh_mass, ...
            last.dynamic_lattice_mass,last.dynamic_medium_mass]) < 1e-4
        endpoint_class = "near_physical_mass_boundary";
    elseif last.unbordered_svmin < 1e-4
        endpoint_class = "near_unbordered_singularity";
    else
        endpoint_class = "unresolved_corrector_domain_boundary";
    end
else
    endpoint_class = stop_reason;
end

result = struct('T',T,'Bx',Bx,'Jsup',Jsup,'h_scale',h_scale, ...
    'lambda_lo',lambda_lo,'lambda_hi',lambda_hi, ...
    'start_source',fullfile(here,'wp2_original_equations_audit.mat'), ...
    'reversibility_delta',reversibility_delta,'initial_ds',ds0, ...
    'points',points,'attempts',attempts,'stop_reason',stop_reason, ...
    'endpoint_class',endpoint_class, ...
    'note',['Bordered pseudo-arclength in fixed normalized coordinates. ' ...
    'Every accepted point passes the exact reduced residual and all current ' ...
    'physical gates; failed correctors are not roots or endpoint proofs.']);
save(fullfile(here,'wp2_reduced_pseudoarclength_1t.mat'),'result','-v7');
fprintf(['pseudo-arclength stopped: %s / %s after %d points; ' ...
    'reversibility delta %.3g\n'],stop_reason,endpoint_class,npoint, ...
    reversibility_delta);

    function u = encode(y,h)
        u = [(y(1:3)-lambda_mid)./lambda_half;-Jsup*y(4);h/h_scale];
    end

    function [y,h] = decode(u)
        y = [lambda_mid+lambda_half.*u(1:3);-u(4)/Jsup];
        h = h_scale*u(5);
    end

    function q = evaluate_q(u)
        [y,h] = decode(u);
        if ~(isfinite(h) && h >= 0)
            q = struct('defined',false,'trial_admissible',false, ...
                'residual',nan(4,1),'residual_norm',NaN,'status','h_domain');
            return;
        end
        try
            ctx = make_context(Bx,h,T,F);
            q = invz_ordered_reduced_residual(y,ctx,ropts);
        catch ME
            q = struct('defined',false,'trial_admissible',false, ...
                'residual',nan(4,1),'residual_norm',NaN, ...
                'status',['error:' ME.identifier]);
        end
    end

    function f = scaled_residual(u)
        q = evaluate_q(u);
        if q.defined && all(isfinite(q.residual))
            f = [q.lambda_residual./lambda_half; ...
                q.static_residual/x_scale];
        else
            f = 100*(1+abs(u(1:4)));
        end
    end

    function [J,status] = finite_jacobian(u)
        f0 = scaled_residual(u);
        J = nan(4,5);
        status = "ok";
        for j = 1:5
            du = zeros(5,1);
            du(j) = 1e-5*max(1,abs(u(j)));
            fp = scaled_residual(u+du);
            if max(abs(fp)) < 50
                J(:,j) = (fp-f0)/du(j);
            else
                fm = scaled_residual(u-du);
                if max(abs(fm)) < 50
                    J(:,j) = (f0-fm)/du(j);
                else
                    status = "undefined_perturbation";
                    return;
                end
            end
        end
    end

    function [corr,ok] = correct_step(u,tangent,ds)
        predictor = u+ds*tangent;
        corr = struct('u',predictor,'J',nan(4,5),'residual',nan(5,1), ...
            'exitflag',NaN,'output',struct(),'predictor',predictor, ...
            'status',"not_run");
        ok = false;
        if any(predictor <= lb) || any(predictor >= ub), return; end
        bordered = @(v) [scaled_residual(v);tangent.'*(v-predictor)];
        try
            [uc,~,fc,exitflag,output,~,JB] = lsqnonlin( ...
                bordered,predictor,lb,ub,lsopts);
        catch ME
            corr.status = "error:"+string(ME.identifier);
            return;
        end
        corr.u = uc;
        corr.J = full(JB(1:4,:));
        corr.residual = fc(:);
        corr.exitflag = exitflag;
        corr.output = output;
        corr.status = "solver_returned";
        q = evaluate_q(uc);
        ok = q.defined && q.trial_admissible && ...
            max(abs(fc(1:4))) <= 1e-8 && abs(fc(5)) <= 1e-8 && ...
            all(uc > lb) && all(uc < ub);
        if ok, corr.status = "accepted"; else, corr.status = "rejected"; end
    end

    function [t,sv] = tangent_from_jacobian(J,orientation)
        J = full(J);
        [~,S,V] = svd(J);
        t = V(:,end);
        t = t/norm(t);
        if ~isempty(orientation) && dot(t,orientation) < 0, t = -t; end
        sv = diag(S);
    end

    function p = point_record(u,q,J,tangent,sv,iterations,ds)
        J = full(J);
        [y,h] = decode(u);
        sy = svd(J(:,1:4));
        p = struct('u',u,'y',y,'h',h,'s',u(4), ...
            'residual_norm',q.residual_norm,'status',string(q.status), ...
            'admissible',q.trial_admissible,'tangent',tangent, ...
            'full_jacobian_singular_values',sv, ...
            'unbordered_svmin',min(sy),'unbordered_rcond',rcond(J(:,1:4)), ...
            'supremum_mass',q.supremum_mass,'uniform_mass',q.D_uni, ...
            'static_mesh_mass',q.static_mesh_min_mass, ...
            'dynamic_lattice_mass',q.dynamic.dynamic_min_lattice_mass, ...
            'dynamic_medium_mass',q.dynamic.dynamic_min_medium_mass, ...
            'xi_denom',q.static.xi_denom, ...
            'gstat_local_denom',q.static.gstat_local_denom, ...
            'dynamic_max_derivative',max(q.dynamic.derivative(2:end)), ...
            'dynamic_max_bound',max(q.dynamic.derivative_bound(2:end)), ...
            'dynamic_enumerated_frequencies', ...
                sum(q.dynamic.proof(2:end) == "finite_enumeration"), ...
            'dynamic_unresolved_frequencies', ...
                sum(q.dynamic.fallback_unresolved(2:end)), ...
            'corrector_iterations',iterations,'ds',ds);
    end
end

function a = attempt_record(index,retry,ds,candidate,ok)
a = struct('index',index,'retry',retry,'ds',ds,'accepted',ok, ...
    'status',candidate.status,'predictor',candidate.predictor, ...
    'candidate',candidate.u,'residual',candidate.residual, ...
    'exitflag',candidate.exitflag);
end

function p = empty_point()
p = struct('u',nan(5,1),'y',nan(4,1),'h',NaN,'s',NaN, ...
    'residual_norm',NaN,'status',"",'admissible',false, ...
    'tangent',nan(5,1),'full_jacobian_singular_values',nan(4,1), ...
    'unbordered_svmin',NaN,'unbordered_rcond',NaN, ...
    'supremum_mass',NaN,'uniform_mass',NaN,'static_mesh_mass',NaN, ...
    'dynamic_lattice_mass',NaN,'dynamic_medium_mass',NaN, ...
    'xi_denom',NaN,'gstat_local_denom',NaN, ...
    'dynamic_max_derivative',NaN,'dynamic_max_bound',NaN, ...
    'dynamic_enumerated_frequencies',NaN, ...
    'dynamic_unresolved_frequencies',NaN, ...
    'corrector_iterations',NaN,'ds',NaN);
end

function a = empty_attempt()
a = struct('index',NaN,'retry',NaN,'ds',NaN,'accepted',false, ...
    'status',"",'predictor',nan(5,1),'candidate',nan(5,1), ...
    'residual',nan(5,1),'exitflag',NaN);
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
