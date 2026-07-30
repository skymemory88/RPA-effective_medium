function result = invzp_node22_step_domain()
%INVZP_NODE22_STEP_DOMAIN Resolve the admissible boundary on failed Picard step 2.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));
F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'));
ion = invz_ion();
T = D.T; Bx = 1; node = 22;
[wn,wts,beta] = invz_matsubara(T,40);
h = D.out(1).nodes{node+1}.h;
si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
ctx = make_context(si,tl,T,wn,wts,beta,F.J,F.info);
opts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false));
mapfun = @(s) invz_ordered_outer_map(s,ctx,opts);

S0 = zeros(size(wn));
M0 = mapfun(S0);
S1 = M0.Sigma_map;
M1 = mapfun(S1);
direction = M1.Sigma_map-S1;
fraction = linspace(0,1,41).';
defined = false(size(fraction));
status = strings(size(fraction));
residual_norm = nan(size(fraction));
D_uni = nan(size(fraction));
for k = 1:numel(fraction)
    q = mapfun(S1+fraction(k)*direction);
    defined(k) = q.defined;
    status(k) = string(q.status);
    residual_norm(k) = q.residual_norm;
    if q.defined, D_uni(k) = q.static.D_uni; end
end
tab = table(fraction,defined,status,residual_norm,D_uni);
probe05 = invz_outer_picard_diagnostic(mapfun,S0, ...
    struct('mix',0.5,'tol',1e-8,'maxit',200));
if probe05.converged
    final05 = invz_ordered_outer_map(probe05.Sigma,ctx, ...
        struct('emt_static',opts.emt_static,'dynamic_diagnostics',true));
    jac05 = invz_outer_dominant_eigen(mapfun,probe05.Sigma, ...
        struct('fd_step',3e-6,'tol',1e-5,'maxit',40));
else
    final05 = [];
    jac05 = [];
end
if ~isempty(probe05.last_admissible_map)
    jac_last = invz_outer_dominant_eigen(mapfun,probe05.last_admissible_Sigma, ...
        struct('fd_step',3e-6,'tol',1e-5,'maxit',40));
else
    jac_last = [];
end
result = struct('Bx',Bx,'node',node,'h',h,'table',tab, ...
    'first_residual',M0.residual_norm,'second_residual',M1.residual_norm, ...
    'probe_mix05',probe05,'final_mix05',final05,'jacobian_mix05',jac05, ...
    'jacobian_last_admissible',jac_last, ...
    'note','Line scan from first Picard state along its full second update.');
save(fullfile(here,'wp2_node22_step_domain.mat'),'result','-v7');
fprintf('node 22 admissible fractions:');
fprintf(' %.3g',fraction(defined));
fprintf('\nfirst undefined fraction %.3g\n',fraction(find(~defined,1)));
fprintf('mix=0.5 probe: %s after %d iterations\n',probe05.status,probe05.iterations);
if ~isempty(jac_last)
    fprintf('last-admissible Jacobian: %s, lambda %.9g, eigen-residual %.3g\n', ...
        jac_last.status,jac_last.lambda,jac_last.eigen_residual);
end
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
