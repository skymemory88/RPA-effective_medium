function result = invzp_4t_adaptive_component_edge_audit()
%INVZP_4T_ADAPTIVE_COMPONENT_EDGE_AUDIT Audit the node-28 branch endpoint.
% Measures the last certified root and first retained below-edge proposal over
% the scan/margin grid, estimates the common zero of the uniform/supremum
% masses, and checks local branch/Jacobian consistency. The zero estimate is
% an empirical local extrapolation, not a thermodynamic selector.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
adaptive_path = fullfile(here,'wp2_4t_adaptive_node28_to27.mat');
F = load(fixture_path);
A = load(adaptive_path);
if A.result.summary.target_reached || numel(A.result.accepted_h) < 4
    error('invz:componentEdgeAudit', ...
        'expected an unresolved gap with at least four certified roots.');
end

T = A.result.config.T;
Bx = A.result.config.Bx;
Ecut = A.result.config.Ecut;
[wn,wts,beta] = invz_matsubara(T,Ecut);
ion = invz_ion();
h = A.result.accepted_h(end);
Sigma = A.result.accepted_Sigma(:,end);
Sigma_previous = A.result.accepted_Sigma(:,end-1);
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
    end
end
edge_root_resolution = table(scan_points,endpoint_margin,status,defined, ...
    n_total_roots,n_admissible_roots,outer_residual,root_residual, ...
    closure_residual,supremum_mass,D_uni,mesh_x_margin, ...
    mesh_medium_margin,dynamic_min_abs);
edge_root_stable = all(defined) && all(n_total_roots == 1) && ...
    all(n_admissible_roots == 1) && all(supremum_mass > 0) && ...
    all(D_uni > 0) && all(mesh_x_margin > 0) && ...
    all(mesh_medium_margin > 0);

accepted_table = A.result.table(A.result.table.accepted,:);
fit_window = [3;4;5];
h_zero_D_uni = nan(size(fit_window));
h_zero_supremum = nan(size(fit_window));
slope_D_uni = nan(size(fit_window));
slope_supremum = nan(size(fit_window));
R2_D_uni = nan(size(fit_window));
R2_supremum = nan(size(fit_window));
for k = 1:numel(fit_window)
    idx = height(accepted_table)-fit_window(k)+1:height(accepted_table);
    hv = accepted_table.candidate_h(idx);
    Dv = accepted_table.final_D_uni(idx);
    Sv = accepted_table.final_supremum_mass(idx);
    [h_zero_D_uni(k),slope_D_uni(k),R2_D_uni(k)] = linear_zero(hv,Dv);
    [h_zero_supremum(k),slope_supremum(k),R2_supremum(k)] = ...
        linear_zero(hv,Sv);
end
mass_zero_fit = table(fit_window,h_zero_D_uni,h_zero_supremum, ...
    slope_D_uni,slope_supremum,R2_D_uni,R2_supremum);

rejected_rows = find(~A.result.table.accepted);
last_rejected = rejected_rows(end);
failed_h = A.result.table.candidate_h(last_rejected);
if abs(A.result.table.from_h(last_rejected)-h) > 10*eps(h)
    error('invz:componentEdgeAudit', ...
        'last retained rejection was not seeded from the last certified root.');
end
failed_ctx = make_context(ion,T,Bx,failed_h,wn,wts,beta,F.J,F.info);
failed_status = strings(n,1);
failed_total_roots = nan(n,1);
failed_admissible_roots = nan(n,1);
failed_min_sample_abs_residual = nan(n,1);
for k = 1:n
    opts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false, ...
        'scan_points',scan_points(k),'endpoint_margin',endpoint_margin(k), ...
        'resid_tol',1e-10));
    q = invz_ordered_outer_map(Sigma,failed_ctx,opts);
    failed_status(k) = string(q.status);
    failed_total_roots(k) = q.static.n_roots;
    failed_admissible_roots(k) = q.static.n_admissible_roots;
    fg = abs(q.static.search.fgrid(isfinite(q.static.search.fgrid)));
    if ~isempty(fg), failed_min_sample_abs_residual(k) = min(fg); end
end
below_edge_resolution = table(scan_points,endpoint_margin,failed_status, ...
    failed_total_roots,failed_admissible_roots, ...
    failed_min_sample_abs_residual);

mapfun = @(s) invz_ordered_outer_map(s,ctx,A.result.config.mapopts);
half_probe = invz_outer_picard_diagnostic(mapfun, ...
    0.5*(Sigma_previous+Sigma),A.result.config.probeopts);
half_delta = NaN;
if half_probe.converged
    half_delta = max(abs(half_probe.Sigma-Sigma));
end
fd_step = [1e-5;3e-6;1e-6];
jac_status = strings(size(fd_step));
jac_converged = false(size(fd_step));
jac_lambda = nan(size(fd_step));
jac_eigen_residual = nan(size(fd_step));
for k = 1:numel(fd_step)
    q = invz_outer_dominant_eigen(mapfun,Sigma, ...
        struct('fd_step',fd_step(k),'tol',1e-5,'maxit',40));
    jac_status(k) = string(q.status);
    jac_converged(k) = q.converged;
    jac_lambda(k) = q.lambda;
    jac_eigen_residual(k) = q.eigen_residual;
end
jacobian_ladder = table(fd_step,jac_status,jac_converged,jac_lambda, ...
    jac_eigen_residual);

result = struct('h',h,'failed_h',failed_h,'Sigma',Sigma, ...
    'edge_root_resolution',edge_root_resolution, ...
    'edge_root_stable',edge_root_stable, ...
    'below_edge_resolution',below_edge_resolution, ...
    'mass_zero_fit',mass_zero_fit,'half_probe_status',half_probe.status, ...
    'half_probe_converged',half_probe.converged,'half_delta',half_delta, ...
    'jacobian_ladder',jacobian_ladder, ...
    'provenance',struct('fixture',fixture_path, ...
                        'adaptive_continuation',adaptive_path), ...
    'note',['Local component-edge evidence only. Linear mass-zero fits do ' ...
            'not prove a global fold or choose thermodynamic equilibrium.']);
save(fullfile(here,'wp2_4t_adaptive_component_edge_audit.mat'), ...
    'result','-v7');
disp(mass_zero_fit);
fprintf(['edge root stable=%d, last h %.12g, failed h %.12g, ' ...
    'min masses [sup %.3g, Duni %.3g], half delta %.3g\n'], ...
    edge_root_stable,h,failed_h,min(supremum_mass),min(D_uni),half_delta);
end

function [hzero,slope,R2] = linear_zero(h,y)
p = polyfit(h,y,1);
slope = p(1);
hzero = -p(2)/p(1);
yfit = polyval(p,h);
den = sum((y-mean(y)).^2);
R2 = 1-sum((y-yfit).^2)/den;
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
