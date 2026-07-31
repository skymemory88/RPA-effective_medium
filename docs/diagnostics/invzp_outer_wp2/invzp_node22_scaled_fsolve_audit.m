function result = invzp_node22_scaled_fsolve_audit()
%INVZP_NODE22_SCALED_FSOLVE_AUDIT Newton audit from bounded-search candidates.
% Uses normalized variables/residuals but no trust-region bounds. Final states
% are independently checked against the rigorous box and all physical gates.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'),'T');
N = load(fullfile(here,'wp2_reduced_node22_search.mat'));
T = D.T;
Bx = N.result.Bx;
h = N.result.h;
Jsup = N.result.Jsup;
ctx = make_context(Bx,h,T,F);
lambda_lo = N.result.lambda_lo;
lambda_hi = N.result.lambda_hi;
lambda_mid = (lambda_lo+lambda_hi)/2;
lambda_half = (lambda_hi-lambda_lo)/2;
x_scale = 1/Jsup;
ropts = struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12));

fsopts = optimoptions('fsolve','Display','off', ...
    'FunctionTolerance',1e-20,'StepTolerance',1e-13, ...
    'OptimalityTolerance',1e-12,'MaxIterations',20, ...
    'MaxFunctionEvaluations',200,'FiniteDifferenceType','central');
nstart = numel(N.result.records);
records = repmat(empty_record(),nstart,1);
for k = 1:nstart
    z0 = N.result.records(k).solution;
    tic;
    [z,fval,exitflag,output,jacobian] = fsolve(@scaled_vector,z0,fsopts);
    elapsed = toc;
    q = invz_ordered_reduced_residual(decode(z),ctx,ropts);
    in_box = all(z > [-ones(3,1);0]) && all(z < ones(4,1));
    records(k) = struct('start',z0,'solution',z,'y',decode(z), ...
        'scaled_residual',fval(:),'scaled_residual_norm',max(abs(fval)), ...
        'exitflag',exitflag,'output',output,'jacobian',jacobian, ...
        'jacobian_rcond',rcond(jacobian),'elapsed_s',elapsed, ...
        'in_rigorous_box',in_box,'status',string(q.status), ...
        'admissible',q.trial_admissible,'physical_residual',q.residual, ...
        'supremum_mass',q.supremum_mass,'uniform_mass',q.D_uni, ...
        'static_mesh_mass',q.static_mesh_min_mass, ...
        'dynamic_lattice_mass',q.dynamic.dynamic_min_lattice_mass, ...
        'dynamic_medium_mass',q.dynamic.dynamic_min_medium_mass);
    fprintf(['scaled fsolve start %d: exit %d, ||Rscaled||inf %.9g, ' ...
        'box %d, admissible %d, rcond(J) %.3g, %.3fs\n'], ...
        k,exitflag,records(k).scaled_residual_norm,in_box, ...
        q.trial_admissible,records(k).jacobian_rcond,elapsed);
end

result = struct('T',T,'Bx',Bx,'h',h,'records',records, ...
    'source',fullfile(here,'wp2_reduced_node22_search.mat'), ...
    'note',['Unbounded normalized Newton audit; final box and physical gates ' ...
    'are checked rather than assumed.']);
save(fullfile(here,'wp2_node22_scaled_fsolve_audit.mat'),'result','-v7');

    function y = decode(z)
        y = [lambda_mid+lambda_half.*z(1:3);-z(4)/Jsup];
    end

    function r = scaled_vector(z)
        q = invz_ordered_reduced_residual(decode(z),ctx,ropts);
        if q.defined && all(isfinite(q.residual))
            r = [q.lambda_residual./lambda_half; ...
                q.static_residual/x_scale];
        else
            r = nan(4,1);
        end
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
    'scaled_residual',nan(4,1),'scaled_residual_norm',NaN, ...
    'exitflag',NaN,'output',struct(),'jacobian',nan(4), ...
    'jacobian_rcond',NaN,'elapsed_s',NaN,'in_rigorous_box',false, ...
    'status',"",'admissible',false,'physical_residual',nan(4,1), ...
    'supremum_mass',NaN,'uniform_mass',NaN,'static_mesh_mass',NaN, ...
    'dynamic_lattice_mass',NaN,'dynamic_medium_mass',NaN);
end
