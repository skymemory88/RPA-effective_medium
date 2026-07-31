function result = invzp_fold_sheet_stability_audit()
%INVZP_FOLD_SHEET_STABILITY_AUDIT Outer-map stability on both 1 T sheets.
% Diagnostic only.  Evaluate the original 740-component outer map and its
% matrix-free dominant Jacobian eigenvalue at the two independently
% constructed h=0.006 roots and at the refined saddle-node.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
A = load(fullfile(here, ...
    'wp2_pseudoarclength_original_equations_audit.mat'));
R = load(fullfile(here,'wp2_reduced_fold_refinement.mat'));
T = A.result.T;
Bx = A.result.Bx;
Jsup = A.result.Jsup;
mapopts = struct('emt_static',struct('Jsup',Jsup,'warn',false));
fd_steps = [1e-6 3e-6];

ctx = make_context(Bx,A.result.hsection,T,F);
mapfun = @(Sigma) invz_ordered_outer_map(Sigma,ctx,mapopts);
sheet = repmat(struct('s',NaN,'base',[],'eigen',[]),2,1);
for k = 1:2
    Sigma = A.result.roots(k).Sigma;
    base = mapfun(Sigma);
    e1 = invz_outer_dominant_eigen(mapfun,Sigma, ...
        struct('fd_step',fd_steps(1),'tol',2e-6,'maxit',80));
    eigen = repmat(e1,numel(fd_steps),1);
    for j = 2:numel(fd_steps)
        eigen(j) = invz_outer_dominant_eigen(mapfun,Sigma, ...
            struct('fd_step',fd_steps(j),'tol',2e-6,'maxit',80));
    end
    sheet(k) = struct('s',-Jsup*A.result.roots(k).x, ...
        'base',base,'eigen',eigen);
end

ctxf = make_context(Bx,R.result.h,T,F);
qf = invz_ordered_reduced_residual(R.result.fold.y,ctxf, ...
    struct('Jsup',Jsup,'dynamic',struct('resid_tol',1e-12)));
mapfunf = @(Sigma) invz_ordered_outer_map(Sigma,ctxf,mapopts);
basef = mapfunf(qf.Sigma);
e1 = invz_outer_dominant_eigen(mapfunf,qf.Sigma, ...
    struct('fd_step',fd_steps(1),'tol',5e-6,'maxit',100));
eigenf = repmat(e1,numel(fd_steps),1);
for j = 2:numel(fd_steps)
    eigenf(j) = invz_outer_dominant_eigen(mapfunf,qf.Sigma, ...
        struct('fd_step',fd_steps(j),'tol',5e-6,'maxit',100));
end
fold = struct('s',R.result.s,'h',R.result.h,'base',basef, ...
    'eigen',eigenf);

result = struct('T',T,'Bx',Bx,'Jsup',Jsup,'hsection',A.result.hsection, ...
    'fd_steps',fd_steps,'sheet',sheet,'fold',fold, ...
    'sources',struct('roots',fullfile(here, ...
        'wp2_pseudoarclength_original_equations_audit.mat'), ...
        'fold',fullfile(here,'wp2_reduced_fold_refinement.mat')), ...
    'note',['Matrix-free central-difference eigenvalues diagnose local ' ...
    'iteration stability only. They neither prove global basin geometry nor ' ...
    'select thermodynamic equilibrium.']);
save(fullfile(here,'wp2_fold_sheet_stability_audit.mat'),'result','-v7');

fprintf(['sheet stability at h %.6g: s %.12g lambda %.9g / %.9g; ' ...
    's %.12g lambda %.9g / %.9g; fold lambda %.9g / %.9g\n'], ...
    A.result.hsection,sheet(1).s,sheet(1).eigen(1).lambda, ...
    sheet(1).eigen(2).lambda,sheet(2).s,sheet(2).eigen(1).lambda, ...
    sheet(2).eigen(2).lambda,eigenf(1).lambda,eigenf(2).lambda);
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
