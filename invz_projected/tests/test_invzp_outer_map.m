function results = test_invzp_outer_map()
%TEST_INVZP_OUTER_MAP Focused deterministic outer-residual gates.
here = fileparts(mfilename('fullpath'));
projected = fileparts(here);
repo = fileparts(projected);
addpath(projected,fullfile(repo,'invz_common'),repo);

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
ion = invz_ion();
T = F.provenance.T;
Bx = F.provenance.Bx;
[wn,wts,beta] = invz_matsubara(T,F.provenance.solve_opts.Ecut);
ctx4 = make_context(F.pt.si,F.pt.tl,T,wn,wts,beta,F.J,F.info);
opts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false), ...
    'dynamic_diagnostics',true);

healthy = invz_ordered_outer_map(F.pt.Sigma,ctx4,opts);
if ~(healthy.defined && strcmp(healthy.status,'ok') && ...
        healthy.static.n_admissible_roots == 1)
    error('invz:testOuterMap','healthy 4 T state did not define a unique outer map.');
end
resid_delta = abs(healthy.residual_norm-F.pt.final_resid);
K0_delta = abs(healthy.K(1)-F.pt.K(1));
lambda_delta = max(abs(healthy.lambda-F.pt.lambda));
if resid_delta > 1e-10 || K0_delta > 1e-10 || lambda_delta > 1e-9 || ...
        healthy.lambda_consistency > 1e-14 || ...
        healthy.dynamic_nonpositive_count ~= 0
    error('invz:testOuterMap', ...
        'healthy outer-map regression failed (dR %.3g, dK %.3g, dlam %.3g).', ...
        resid_delta,K0_delta,lambda_delta);
end

% The map is a function of Sigma, not of a hidden prior lambda/K0 seed.
zero_at_healthy_node = invz_ordered_outer_map(zeros(size(wn)),ctx4,opts);
if ~(zero_at_healthy_node.defined && ...
        zero_at_healthy_node.static.n_admissible_roots == 1 && ...
        zero_at_healthy_node.lambda_consistency <= 1e-14)
    error('invz:testOuterMap','healthy-node Sigma=0 map was not deterministic/defined.');
end

% The h=0 predictor at the same Bx is outside the map domain at Sigma=0.
si0 = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',0));
tl0 = invz_twolevel_ordered(ion,T,Bx,0,struct('Jxx0',F.info.Jaa0));
ctx0 = make_context(si0,tl0,T,wn,wts,beta,F.J,F.info);
predictor = invz_ordered_outer_map(zeros(size(wn)),ctx0,opts);
if predictor.defined || ~strcmp(predictor.status,'no_admissible_static_root')
    error('invz:testOuterMap','4 T h=0 predictor domain classification changed.');
end

% Matrix-free Jacobian diagnostic: verify the kernel on a known linear map,
% then measure (without solving) the healthy production map.
A = diag([2 -.5 .1]);
linear_map = @(s) struct('status','ok','defined',true, ...
    'Sigma_map',A*s(:)+[1;2;3]);
linear_jacobian = invz_outer_dominant_eigen(linear_map,zeros(3,1), ...
    struct('fd_step',1e-6,'tol',1e-9,'maxit',50));
if ~(linear_jacobian.converged && abs(linear_jacobian.lambda-2) < 1e-8)
    error('invz:testOuterMap','matrix-free Jacobian failed the known linear map.');
end

% Power iteration need not settle on a leading complex pair. Verify that the
% Arnoldi diagnostic recovers the radius of a known real rotation map.
rho = 0.8;
theta = 0.37;
R = rho*[cos(theta) -sin(theta); sin(theta) cos(theta)];
Arot = blkdiag(R,zeros(4));
rotation_map = @(s) struct('status','ok','defined',true, ...
    'Sigma_map',Arot*s(:)+(1:6).');
rotation_arnoldi = invz_outer_arnoldi_diagnostic( ...
    rotation_map,zeros(6,1), ...
    struct('n_eigs',4,'fd_step',1e-6,'tol',1e-8,'maxit',100));
if ~(rotation_arnoldi.converged && ...
        abs(rotation_arnoldi.spectral_radius-rho) < 1e-8 && ...
        nnz(rotation_arnoldi.active_mode_mask) == 2 && ...
        rotation_arnoldi.max_active_residual < 1e-8)
    error('invz:testOuterMap', ...
        'Arnoldi diagnostic failed the known complex-pair rotation map.');
end
boundary_arnoldi = invz_outer_arnoldi_diagnostic( ...
    @bounded_test_map,zeros(6,1), ...
    struct('n_eigs',2,'fd_step',1e-6,'tol',1e-8,'maxit',20));
if ~strcmp(boundary_arnoldi.status,'domain_boundary') || ...
        any(~isfinite(boundary_arnoldi.boundary_direction)) || ...
        ~isfinite(boundary_arnoldi.boundary_delta)
    error('invz:testOuterMap', ...
        'Arnoldi diagnostic did not stop at a synthetic map-domain boundary.');
end
jac_opts = struct('emt_static',opts.emt_static);
mapfun = @(s) invz_ordered_outer_map(s,ctx4,jac_opts);
healthy_jacobian = invz_outer_dominant_eigen(mapfun,F.pt.Sigma, ...
    struct('fd_step',3e-6,'tol',1e-5,'maxit',40));
if ~(healthy_jacobian.converged && healthy_jacobian.eigen_residual < 1e-4 && ...
        healthy_jacobian.lambda > -.02 && healthy_jacobian.lambda < 0)
    error('invz:testOuterMap','healthy 4 T dominant Jacobian diagnostic changed.');
end

results = struct();
results.healthy = healthy;
results.zero_at_healthy_node = zero_at_healthy_node;
results.predictor = predictor;
results.linear_jacobian = linear_jacobian;
results.rotation_arnoldi = rotation_arnoldi;
results.boundary_arnoldi = boundary_arnoldi;
results.healthy_jacobian = healthy_jacobian;
results.regression = struct('residual_delta',resid_delta,'K0_delta',K0_delta, ...
    'lambda_delta',lambda_delta);
results.provenance = struct('fixture', ...
    fullfile(repo,'docs','diagnostics','invzp_static_wp1','legacy_4T_fixture.mat'), ...
    'test',mfilename);
diag_path = fullfile(repo,'docs','diagnostics','invzp_outer_wp2', ...
    'wp2_outer_map_gate.mat');
save(diag_path,'results','-v7');
fprintf(['test_invzp_outer_map: healthy residual %.9g, dominant lambda %.9g; ' ...
    'h=0 predictor=%s\n'],healthy.residual_norm,healthy_jacobian.lambda, ...
    predictor.status);
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

function out = bounded_test_map(Sigma)
defined = norm(Sigma,inf) <= 1e-9;
if defined
    out = struct('status','ok','defined',true,'Sigma_map',0.2*Sigma(:));
else
    out = struct('status','synthetic_boundary','defined',false, ...
        'Sigma_map',nan(size(Sigma(:))));
end
end
