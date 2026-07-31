function result = invzp_coupling_endpoint_refinement()
%INVZP_COUPLING_ENDPOINT_REFINEMENT Refine the upper sheet at s -> 1.
% Diagnostic only.  This is independent of fixed-rho continuation: prescribe
% the positive supremum mass epsilon=1-s and solve the four exact reduced
% equations for lambda(1:3) and the coupling scale rho.  A finite limiting
% rho with uniform mass -> 0, while the mesh/dynamic masses stay positive,
% identifies a physical stability-boundary endpoint rather than a fold in
% the interior of the admissible domain.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
C = load(fullfile(here,'wp2_coupling_constant_path_audit.mat'));
T = C.result.T;
Bx = C.result.Bx;
h = C.result.h;
Jbase = F.J(:);
J0base = F.info.Jcc0;
Jsupbase = C.result.Jsupbase;
base = make_context(Bx,h,T,F);
[lambda_lo,lambda_hi] = lambda_bounds(base);
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;

seed = C.result.branches(2).failure.last_accepted;
assert(seed.accepted && isfinite(seed.rho));
z = [(seed.y(1:3)/seed.rho-lambda_mid)./lambda_half;seed.rho];
epsilon_grid = [seed.supremum_mass,5e-5,2e-5,1e-5, ...
    5e-6,2e-6,1e-6].';
lb = [-ones(3,1)+1e-9;0.70];
ub = [ ones(3,1)-1e-9;0.85];
lsopts = optimoptions('lsqnonlin','Display','off', ...
    'FunctionTolerance',1e-16,'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-10,'MaxIterations',20, ...
    'MaxFunctionEvaluations',120,'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',2e-6);
records = repmat(empty_record(),numel(epsilon_grid),1);

for k = 1:numel(epsilon_grid)
    epsilon = epsilon_grid(k);
    zstart = z;
    rstart = scaled_vector(zstart,epsilon);
    fprintf('  endpoint start eps %.3g: scaled R %.3g\n', ...
        epsilon,max(abs(rstart)));
    if max(abs(rstart)) <= 2e-8
        z = zstart;
        rscaled = rstart;
        resnorm = sum(rscaled.^2);
        exitflag = 1;
        output = struct('message','Accepted supplied exact root.');
    else
        [z,resnorm,rscaled,exitflag,output] = lsqnonlin( ...
            @(zz) scaled_vector(zz,epsilon),zstart,lb,ub,lsopts);
    end
    [qroot,yroot] = evaluate(z,epsilon);
    J = finite_jacobian(@(zz) scaled_vector(zz,epsilon),z);
    sv = svd(J);
    accepted = exitflag > 0 && qroot.defined && qroot.trial_admissible && ...
        qroot.residual_norm <= 2e-8 && max(abs(rscaled)) <= 2e-8 && ...
        sum(qroot.dynamic.fallback_unresolved(2:end)) == 0;
    records(k) = struct('epsilon',epsilon,'s',1-epsilon, ...
        'start',zstart,'z',z,'rho',z(4),'y',yroot,'q',qroot, ...
        'resnorm',resnorm,'scaled_residual',rscaled(:), ...
        'scaled_residual_norm',max(abs(rscaled)), ...
        'physical_residual_norm',qroot.residual_norm,'exitflag',exitflag, ...
        'output',output,'accepted',accepted, ...
        'jacobian_singular_values',sv,'jacobian_svmin',min(sv), ...
        'supremum_mass',qroot.supremum_mass,'uniform_mass',qroot.D_uni, ...
        'static_mesh_mass',qroot.static_mesh_min_mass, ...
        'dynamic_lattice_mass',qroot.dynamic.dynamic_min_lattice_mass, ...
        'dynamic_medium_mass',qroot.dynamic.dynamic_min_medium_mass);
    fprintf(['coupling endpoint eps %.3g: rho %.12g, R %.3g, ' ...
        'masses [%.3g %.3g %.3g %.3g %.3g], svmin %.3g, ok %d\n'], ...
        epsilon,z(4),qroot.residual_norm,qroot.supremum_mass,qroot.D_uni, ...
        qroot.static_mesh_min_mass, ...
        qroot.dynamic.dynamic_min_lattice_mass, ...
        qroot.dynamic.dynamic_min_medium_mass,min(sv),accepted);
    assert(accepted,'Endpoint solve failed at epsilon %.9g.',epsilon);
end

eps_tail = epsilon_grid(end-4:end);
rho_tail = [records(end-4:end).rho].';
p1 = polyfit(eps_tail,rho_tail,1);
p2 = polyfit(eps_tail,rho_tail,2);
rho_endpoint_linear = polyval(p1,0);
rho_endpoint_quadratic = polyval(p2,0);
extrapolation_spread = abs(rho_endpoint_linear-rho_endpoint_quadratic);
uniform_ratio = [records(end-4:end).uniform_mass].'./eps_tail;
interior_other_masses = [[records.static_mesh_mass].', ...
    [records.dynamic_lattice_mass].',[records.dynamic_medium_mass].'];

result = struct('T',T,'Bx',Bx,'h',h,'J0base',J0base, ...
    'Jsupbase',Jsupbase,'epsilon_grid',epsilon_grid, ...
    'records',records,'rho_endpoint_linear',rho_endpoint_linear, ...
    'rho_endpoint_quadratic',rho_endpoint_quadratic, ...
    'extrapolation_spread',extrapolation_spread, ...
    'tail_uniform_mass_over_epsilon',uniform_ratio, ...
    'minimum_other_mass',min(interior_other_masses,[],'all'), ...
    'minimum_endpoint_jacobian_sv',min([records.jacobian_svmin]), ...
    'source',fullfile(here,'wp2_coupling_constant_path_audit.mat'), ...
    'interpretation',['The prescribed-epsilon solve distinguishes an ' ...
        's=1 stability-boundary endpoint from an interior fixed-rho fold. ' ...
        'A finite rho limit is not a zero-coupling thermodynamic anchor.']);
save(fullfile(here,'wp2_coupling_endpoint_refinement.mat'), ...
    'result','-v7');
fprintf(['coupling endpoint: rho*=%.12g (linear), %.12g ' ...
    '(quadratic), spread %.3g; min other mass %.3g, min sv %.3g\n'], ...
    rho_endpoint_linear,rho_endpoint_quadratic,extrapolation_spread, ...
    result.minimum_other_mass,result.minimum_endpoint_jacobian_sv);

    function r = scaled_vector(zz,epsilon)
        [qscale,~] = evaluate(zz,epsilon);
        rho = zz(4);
        if qscale.defined && all(isfinite(qscale.residual))
            r = [qscale.lambda_residual./(rho*lambda_half); ...
                J0base*qscale.static_residual];
        else
            r = 1e3*(1+abs(zz));
        end
    end

    function [qout,yout] = evaluate(zz,epsilon)
        rho = zz(4);
        lambda = rho*(lambda_mid+lambda_half.*zz(1:3));
        x = -(1-epsilon)/(rho*Jsupbase);
        yout = [lambda;x];
        ctx = base;
        ctx.Jnu_flat = rho*Jbase;
        ctx.J0eff = rho*J0base;
        qout = invz_ordered_reduced_residual(yout,ctx, ...
            struct('Jsup',rho*Jsupbase, ...
            'dynamic',struct('resid_tol',1e-12)));
    end
end

function J = finite_jacobian(fun,z)
f0 = fun(z);
J = nan(numel(f0),numel(z));
for k = 1:numel(z)
    dz = 2e-5*max(1,abs(z(k)));
    zp = z;
    zp(k) = zp(k)+dz;
    J(:,k) = (fun(zp)-f0)/dz;
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

function r = empty_record()
r = struct('epsilon',NaN,'s',NaN,'start',nan(4,1), ...
    'z',nan(4,1),'rho',NaN,'y',nan(4,1),'q',struct(), ...
    'resnorm',NaN,'scaled_residual',nan(4,1), ...
    'scaled_residual_norm',NaN,'physical_residual_norm',NaN, ...
    'exitflag',NaN,'output',struct(),'accepted',false, ...
    'jacobian_singular_values',nan(4,1),'jacobian_svmin',NaN, ...
    'supremum_mass',NaN,'uniform_mass',NaN,'static_mesh_mass',NaN, ...
    'dynamic_lattice_mass',NaN,'dynamic_medium_mass',NaN);
end
