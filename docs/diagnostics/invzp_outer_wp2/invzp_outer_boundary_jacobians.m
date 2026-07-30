function result = invzp_outer_boundary_jacobians()
%INVZP_OUTER_BOUNDARY_JACOBIANS Local map diagnostics at selected Sigma=0 nodes.
% Measurement only; no outer iteration or root solver.
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
fixtures = [1 20; 1 22; 3 28; 3 33];

n = size(fixtures,1);
Bx = fixtures(:,1);
node = fixtures(:,2);
h = nan(n,1);
map_status = strings(n,1);
residual_norm = nan(n,1);
jac_status = strings(n,1);
lambda_dom = nan(n,1);
eigen_residual = nan(n,1);
iterations = nan(n,1);
details = cell(n,1);
for k = 1:n
    ib = find([D.out.Bx] == Bx(k),1);
    h(k) = D.out(ib).nodes{node(k)+1}.h;
    si = invz_single_ion(ion,T,[Bx(k) 0 0], ...
        struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h(k)));
    tl = invz_twolevel_ordered(ion,T,Bx(k),h(k),struct('Jxx0',F.info.Jaa0));
    ctx = make_context(si,tl,T,wn,wts,beta,F.J,F.info);
    mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
    base = mapfun(zeros(size(wn)));
    jac = invz_outer_dominant_eigen(mapfun,zeros(size(wn)), ...
        struct('fd_step',3e-6,'tol',1e-5,'maxit',40));
    map_status(k) = string(base.status);
    residual_norm(k) = base.residual_norm;
    jac_status(k) = string(jac.status);
    lambda_dom(k) = jac.lambda;
    eigen_residual(k) = jac.eigen_residual;
    iterations(k) = jac.iterations;
    details{k} = struct('base',base,'jacobian',jac);
end
tab = table(Bx,node,h,map_status,residual_norm,jac_status,lambda_dom, ...
    eigen_residual,iterations);
result = struct('table',tab,'details',{details}, ...
    'note','Local Sigma=0 measurements only; not coupled-root verdicts.');
save(fullfile(here,'wp2_outer_boundary_jacobians.mat'),'result','-v7');
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
