function result = invzp_outer_contracting_probes()
%INVZP_OUTER_CONTRACTING_PROBES Undamped probes only where the measured map contracts.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'));
ion = invz_ion();
T = D.T;
[wn,wts,beta] = invz_matsubara(T,40);
mapopts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false));
fixtures = [1 22; 3 28; 3 33];

n = size(fixtures,1);
Bx = fixtures(:,1);
node = fixtures(:,2);
h = nan(n,1);
probe_status = strings(n,1);
iterations = nan(n,1);
final_residual = nan(n,1);
lambda_dom = nan(n,1);
eigen_residual = nan(n,1);
D_uni = nan(n,1);
dynamic_min_abs = nan(n,1);
half_start_delta = nan(n,1);
mixed_delta = nan(n,1);
half_start_iterations = nan(n,1);
mixed_iterations = nan(n,1);
details = cell(n,1);
for k = 1:n
    ib = find([D.out.Bx] == Bx(k),1);
    h(k) = D.out(ib).nodes{node(k)+1}.h;
    si = invz_single_ion(ion,T,[Bx(k) 0 0], ...
        struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h(k)));
    tl = invz_twolevel_ordered(ion,T,Bx(k),h(k),struct('Jxx0',F.info.Jaa0));
    ctx = make_context(si,tl,T,wn,wts,beta,F.J,F.info);
    mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
    probe = invz_outer_picard_diagnostic(mapfun,zeros(size(wn)), ...
        struct('mix',1,'tol',1e-8,'maxit',200));
    probe_status(k) = string(probe.status);
    iterations(k) = probe.iterations;
    if probe.converged
        final = invz_ordered_outer_map(probe.Sigma,ctx, ...
            struct('emt_static',mapopts.emt_static,'dynamic_diagnostics',true));
        jac = invz_outer_dominant_eigen(mapfun,probe.Sigma, ...
            struct('fd_step',3e-6,'tol',1e-5,'maxit',40));
        final_residual(k) = final.residual_norm;
        lambda_dom(k) = jac.lambda;
        eigen_residual(k) = jac.eigen_residual;
        D_uni(k) = final.static.D_uni;
        dynamic_min_abs(k) = final.dynamic_min_abs;
        half_probe = invz_outer_picard_diagnostic(mapfun,0.5*probe.Sigma, ...
            struct('mix',1,'tol',1e-8,'maxit',200));
        mixed_probe = invz_outer_picard_diagnostic(mapfun,zeros(size(wn)), ...
            struct('mix',0.5,'tol',1e-8,'maxit',200));
        if half_probe.converged
            half_start_delta(k) = max(abs(half_probe.Sigma-probe.Sigma));
        end
        if mixed_probe.converged
            mixed_delta(k) = max(abs(mixed_probe.Sigma-probe.Sigma));
        end
        half_start_iterations(k) = half_probe.iterations;
        mixed_iterations(k) = mixed_probe.iterations;
    else
        final = [];
        jac = [];
        half_probe = [];
        mixed_probe = [];
    end
    details{k} = struct('probe',probe,'final',final,'jacobian',jac, ...
        'half_start_probe',half_probe,'mixed_probe',mixed_probe);
end
tab = table(Bx,node,h,probe_status,iterations,final_residual,lambda_dom, ...
    eigen_residual,D_uni,dynamic_min_abs,half_start_delta,mixed_delta, ...
    half_start_iterations,mixed_iterations);
result = struct('table',tab,'details',{details}, ...
    'note',['Existence probes only at fixtures preselected by a converged ' ...
            'contractive Sigma=0 Jacobian; not a universal solver test.']);
save(fullfile(here,'wp2_outer_contracting_probes.mat'),'result','-v7');
disp(tab);
end

function ctx = make_context(si,tl,T,wn,wts,beta,J,info)
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,T,1i*wn(1),struct('elastic',false));
G0i = -real(c0i(3,3,1));
X = real(c0(:,:,1));
feedback = X(3,1)*(info.Jaa0/(1-info.Jaa0*X(1,1)))*X(1,3);
G0e = -(X(3,3)+feedback)-G0i;
g = real(invz_g(tl,1i*wn));
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',J,'J0eff',info.Jcc0,'G0inel0',G0i,'G0el0',G0e);
end
