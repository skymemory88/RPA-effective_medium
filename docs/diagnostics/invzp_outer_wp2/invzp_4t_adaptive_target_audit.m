function result = invzp_4t_adaptive_target_audit(source_name,output_name)
%INVZP_4T_ADAPTIVE_TARGET_AUDIT Resolution and branch gates at a 4 T target.
% Re-evaluates the adaptively continued target root over an independent
% scan-density/endpoint grid, a finite-difference-step Jacobian ladder, and
% one halfway-to-root undamped start. No production state is modified.
if nargin < 1
    source_name = "wp2_4t_adaptive_boundary_continuation.mat";
end
if nargin < 2
    output_name = "wp2_4t_adaptive_target_audit.mat";
end
source_name = string(source_name);
output_name = string(output_name);
if ~(isscalar(source_name) && isscalar(output_name))
    error('invz:adaptiveTargetAudit', ...
        'source_name and output_name must be scalar strings.');
end
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
adaptive_path = fullfile(here,source_name);
F = load(fixture_path);
A = load(adaptive_path);
if ~A.result.summary.target_reached
    error('invz:adaptiveTargetAudit', ...
        'adaptive-continuation artifact did not reach its target.');
end

T = A.result.config.T;
Bx = A.result.config.Bx;
h = A.result.summary.target_h;
Sigma = A.result.accepted_Sigma(:,end);
Sigma_previous = A.result.accepted_Sigma(:,end-1);
[wn,wts,beta] = invz_matsubara(T,A.result.config.Ecut);
ion = invz_ion();
ctx = make_context(ion,T,Bx,h,wn,wts,beta,F.J,F.info);

[scan_grid,margin_grid] = ndgrid([2049 4097 8193],[1e-8 1e-10 1e-12]);
scan_points = scan_grid(:);
endpoint_margin = margin_grid(:);
n = numel(scan_points);
status = strings(n,1);
defined = false(n,1);
n_total_roots = nan(n,1);
n_admissible_roots = nan(n,1);
outer_residual = nan(n,1);
root_residual = nan(n,1);
closure_residual = nan(n,1);
supremum_mass = nan(n,1);
D_uni = nan(n,1);
mesh_x_margin = nan(n,1);
mesh_medium_margin = nan(n,1);
dynamic_min_abs = nan(n,1);
dynamic_nonpositive_count = nan(n,1);
for k = 1:n
    opts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false, ...
        'scan_points',scan_points(k),'endpoint_margin',endpoint_margin(k), ...
        'resid_tol',1e-10),'dynamic_diagnostics',true);
    q = invz_ordered_outer_map(Sigma,ctx,opts);
    status(k) = string(q.status);
    defined(k) = q.defined;
    n_total_roots(k) = q.static.n_roots;
    n_admissible_roots(k) = q.static.n_admissible_roots;
    outer_residual(k) = q.residual_norm;
    if q.defined
        row = q.static.selected_index;
        rt = q.static.root_table;
        root_residual(k) = rt.root_resid(row);
        closure_residual(k) = rt.closure_resid(row);
        supremum_mass(k) = rt.supremum_mass(row);
        D_uni(k) = q.static.D_uni;
        mesh_x_margin(k) = rt.min_mesh_x_signed(row);
        mesh_medium_margin(k) = rt.min_mesh_medium_signed(row);
        dynamic_min_abs(k) = q.dynamic_min_abs;
        dynamic_nonpositive_count(k) = q.dynamic_nonpositive_count;
    end
end
resolution = table(scan_points,endpoint_margin,status,defined, ...
    n_total_roots,n_admissible_roots,outer_residual,root_residual, ...
    closure_residual,supremum_mass,D_uni,mesh_x_margin, ...
    mesh_medium_margin,dynamic_min_abs,dynamic_nonpositive_count);
resolution_stable = all(defined) && all(n_total_roots == 1) && ...
    all(n_admissible_roots == 1) && ...
    max(outer_residual) <= A.result.config.probeopts.tol && ...
    all(supremum_mass > 0) && all(D_uni > 0) && ...
    all(mesh_x_margin > 0) && all(mesh_medium_margin > 0) && ...
    all(dynamic_nonpositive_count == 0);

mapopts = A.result.config.mapopts;
mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
Sigma_half = 0.5*(Sigma_previous+Sigma);
half_probe = invz_outer_picard_diagnostic(mapfun,Sigma_half, ...
    A.result.config.probeopts);
half_delta = NaN;
if half_probe.converged
    half_delta = max(abs(half_probe.Sigma-Sigma));
end

fd_step = [1e-5;3e-6;1e-6];
jac_status = strings(size(fd_step));
jac_converged = false(size(fd_step));
jac_lambda = nan(size(fd_step));
jac_eigen_residual = nan(size(fd_step));
jac_details = cell(size(fd_step));
for k = 1:numel(fd_step)
    q = invz_outer_dominant_eigen(mapfun,Sigma, ...
        struct('fd_step',fd_step(k),'tol',1e-5,'maxit',40));
    jac_status(k) = string(q.status);
    jac_converged(k) = q.converged;
    jac_lambda(k) = q.lambda;
    jac_eigen_residual(k) = q.eigen_residual;
    jac_details{k} = rmfield(q,'vector');
end
jacobian_ladder = table(fd_step,jac_status,jac_converged,jac_lambda, ...
    jac_eigen_residual);

result = struct('h',h,'Sigma',Sigma,'resolution',resolution, ...
    'resolution_stable',resolution_stable,'half_probe', ...
    compact_probe(half_probe),'half_delta',half_delta, ...
    'jacobian_ladder',jacobian_ladder,'jacobian_details',{jac_details}, ...
    'provenance',struct('fixture',fixture_path, ...
                        'adaptive_continuation',adaptive_path), ...
    'note',['Target-root numerical and branch-consistency audit only; it ' ...
            'does not evaluate the H_MF integral or select equilibrium.']);
save(fullfile(here,output_name),'result','-v7');
disp(jacobian_ladder);
fprintf(['target resolution stable=%d, half-start=%s, delta %.3g, ' ...
    'min masses [sup %.6g, Duni %.6g, mesh %.6g]\n'], ...
    resolution_stable,half_probe.status,half_delta,min(supremum_mass), ...
    min(D_uni),min(mesh_medium_margin));
end

function ctx = make_context(ion,T,Bx,h,wn,wts,beta,J,info)
si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',info.Jaa0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',info.Jaa0));
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

function out = compact_probe(q)
out = struct('status',q.status,'converged',q.converged, ...
    'Sigma',q.Sigma,'iterations',q.iterations, ...
    'residual_history',q.residual_history,'map_status',q.map_status, ...
    'mix',q.mix,'tol',q.tol);
end
