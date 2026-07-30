function result = invzp_4t_reverse_continuation()
%INVZP_4T_REVERSE_CONTINUATION Diagnostic high-to-low H_MF branch probe.
% This does not evaluate the ordered H_MF integral or modify the production
% solver. At each 4 T profile node it:
%   1. seeds Sigma from the immediately preceding certified higher-h node;
%   2. requires a uniquely defined bounded outer map at that seed;
%   3. requires a resolved leading-mode measurement with spectral radius < 1;
%   4. runs undamped, domain-gated Picard; and
%   5. certifies the final static, dynamic, finite-value, and residual gates.
% Power iteration is used first. If it cannot resolve a rotating/clustered
% subspace, a three-step finite-difference Arnoldi ladder is the only fallback.
% The sweep stops at the first failed requirement.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
failed_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'new_4T_fullsolve_failed.mat');
F = load(fixture_path);
N = load(failed_path);

ion = invz_ion();
T = F.provenance.T;
Bx = F.provenance.Bx;
Ecut = F.provenance.solve_opts.Ecut;
[wn,wts,beta] = invz_matsubara(T,Ecut);
h = N.pt_new.hmf_prof.hgrid(:);
legacy_certified = N.pt_new.hmf_prof.node_conv(:);
legacy_static_status = N.pt_new.hmf_prof.static_status(:);
nH = numel(h);

mapopts = struct('emt_static', ...
    struct('Jsup',F.info.Jcc0,'warn',false,'scan_points',4097, ...
           'endpoint_margin',1e-10,'resid_tol',1e-10));
jacopts = struct('fd_step',3e-6,'tol',1e-5,'maxit',40);
probeopts = struct('mix',1,'tol',1e-8,'maxit',200);

attempted = false(nH,1);
certified = false(nH,1);
seed_source = strings(nH,1);
seed_status = repmat("not_attempted",nH,1);
seed_residual = nan(nH,1);
seed_D_uni = nan(nH,1);
seed_static_roots = nan(nH,1);
jac_status = repmat("not_attempted",nH,1);
jac_lambda = nan(nH,1);
jac_eigen_residual = nan(nH,1);
jac_iterations = nan(nH,1);
linear_gate_method = repmat("not_attempted",nH,1);
linear_gate_status = repmat("not_attempted",nH,1);
linear_spectral_radius = nan(nH,1);
arnoldi_radius_spread = nan(nH,1);
probe_status = repmat("not_attempted",nH,1);
probe_iterations = nan(nH,1);
final_residual = nan(nH,1);
final_D_uni = nan(nH,1);
final_supremum_mass = nan(nH,1);
final_mesh_x_margin = nan(nH,1);
final_mesh_medium_margin = nan(nH,1);
final_dynamic_min_abs = nan(nH,1);
final_dynamic_nonpositive_count = nan(nH,1);
stop_reason = strings(nH,1);
Sigma_solution = nan(numel(wn),nH);
arnoldi_details = cell(nH,1);
stop_detail = struct([]);

for inode = nH:-1:1
    attempted(inode) = true;
    if inode == nH
        Sigma_seed = zeros(size(wn));
        seed_source(inode) = "zero";
    else
        Sigma_seed = Sigma_solution(:,inode+1);
        seed_source(inode) = "certified_node_" + string(inode+1);
    end

    si = invz_single_ion(ion,T,[Bx 0 0], ...
        struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h(inode)));
    tl = invz_twolevel_ordered(ion,T,Bx,h(inode), ...
        struct('Jxx0',F.info.Jaa0));
    ctx = make_context(si,tl,T,wn,wts,beta,F.J,F.info);
    mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);

    base = mapfun(Sigma_seed);
    seed_status(inode) = string(base.status);
    seed_residual(inode) = base.residual_norm;
    if ~isempty(base.static)
        seed_static_roots(inode) = base.static.n_admissible_roots;
        seed_D_uni(inode) = base.static.D_uni;
    end
    if ~base.defined
        stop_reason(inode) = "seed_map_" + string(base.status);
        stop_detail = struct('node',inode,'h',h(inode),'stage',"seed_map", ...
            'base',base,'jacobian',[],'arnoldi',[],'probe',[],'final',[]);
        fprintf('4 T reverse node %d h=%.9g: STOP %s\n', ...
            inode,h(inode),stop_reason(inode));
        break;
    end

    jac = invz_outer_dominant_eigen(mapfun,Sigma_seed,jacopts);
    jac_status(inode) = string(jac.status);
    jac_lambda(inode) = jac.lambda;
    jac_eigen_residual(inode) = jac.eigen_residual;
    jac_iterations(inode) = jac.iterations;

    if jac.converged && isfinite(jac.lambda)
        linear_gate_method(inode) = "power";
        linear_gate_status(inode) = "ok";
        linear_spectral_radius(inode) = abs(jac.lambda);
    else
        fd_steps = [1e-5;3e-6;1e-6];
        aruns = cell(numel(fd_steps),1);
        aradius = nan(size(fd_steps));
        aconverged = false(size(fd_steps));
        for ia = 1:numel(fd_steps)
            aopts = struct('n_eigs',6,'fd_step',fd_steps(ia), ...
                'tol',1e-5,'maxit',100,'subspace_dimension',20);
            aruns{ia} = invz_outer_arnoldi_diagnostic( ...
                mapfun,Sigma_seed,aopts);
            aconverged(ia) = aruns{ia}.converged;
            aradius(ia) = aruns{ia}.spectral_radius;
        end
        radius_spread = max(aradius)-min(aradius);
        radius_stable = all(isfinite(aradius)) && ...
            radius_spread <= max(1e-4,0.05*median(aradius));
        arnoldi_details{inode} = struct('fd_steps',fd_steps, ...
            'runs',{aruns},'spectral_radius',aradius, ...
            'radius_spread',radius_spread,'radius_stable',radius_stable);
        linear_gate_method(inode) = "arnoldi_fd_ladder";
        arnoldi_radius_spread(inode) = radius_spread;
        if all(aconverged) && radius_stable
            linear_gate_status(inode) = "ok";
            linear_spectral_radius(inode) = max(aradius);
        elseif ~all(aconverged)
            linear_gate_status(inode) = "arnoldi_unresolved";
        else
            linear_gate_status(inode) = "arnoldi_fd_unstable";
        end
    end

    if ~strcmp(linear_gate_status(inode),"ok") || ...
            ~(isfinite(linear_spectral_radius(inode)) && ...
              linear_spectral_radius(inode) < 1)
        if strcmp(linear_gate_status(inode),"ok")
            stop_reason(inode) = "seed_leading_modes_noncontractive";
        else
            stop_reason(inode) = "seed_linear_gate_" + ...
                linear_gate_status(inode);
        end
        stop_detail = struct('node',inode,'h',h(inode),'stage',"jacobian", ...
            'base',base,'jacobian',jac,'arnoldi',arnoldi_details{inode}, ...
            'probe',[],'final',[]);
        fprintf('4 T reverse node %d h=%.9g: STOP %s (radius %.6g)\n', ...
            inode,h(inode),stop_reason(inode),linear_spectral_radius(inode));
        break;
    end

    probe = invz_outer_picard_diagnostic(mapfun,Sigma_seed,probeopts);
    probe_status(inode) = string(probe.status);
    probe_iterations(inode) = probe.iterations;
    if ~probe.converged
        stop_reason(inode) = "picard_" + string(probe.status);
        stop_detail = struct('node',inode,'h',h(inode),'stage',"picard", ...
            'base',base,'jacobian',jac,'arnoldi',arnoldi_details{inode}, ...
            'probe',probe,'final',[]);
        fprintf('4 T reverse node %d h=%.9g: STOP %s after %d iterations\n', ...
            inode,h(inode),stop_reason(inode),probe.iterations);
        break;
    end

    final = invz_ordered_outer_map(probe.Sigma,ctx, ...
        struct('emt_static',mapopts.emt_static,'dynamic_diagnostics',true));
    final_residual(inode) = final.residual_norm;
    if final.defined
        row = final.static.selected_index;
        rt = final.static.root_table;
        final_D_uni(inode) = final.static.D_uni;
        final_supremum_mass(inode) = rt.supremum_mass(row);
        final_mesh_x_margin(inode) = rt.min_mesh_x_signed(row);
        final_mesh_medium_margin(inode) = rt.min_mesh_medium_signed(row);
        final_dynamic_min_abs(inode) = final.dynamic_min_abs;
        final_dynamic_nonpositive_count(inode) = ...
            final.dynamic_nonpositive_count;
    end
    gate_ok = final.defined && final.static.n_admissible_roots == 1 && ...
        isfinite(final.residual_norm) && final.residual_norm <= probeopts.tol && ...
        all(isfinite(probe.Sigma)) && final.static.D_uni > 0 && ...
        final_supremum_mass(inode) > 0 && ...
        final_mesh_x_margin(inode) > 0 && ...
        final_mesh_medium_margin(inode) > 0 && ...
        final.dynamic_nonpositive_count == 0 && ...
        isfinite(final.dynamic_min_abs) && final.dynamic_min_abs > 0;
    if ~gate_ok
        stop_reason(inode) = "final_acceptance_gate";
        stop_detail = struct('node',inode,'h',h(inode),'stage',"final_gate", ...
            'base',base,'jacobian',jac,'arnoldi',arnoldi_details{inode}, ...
            'probe',probe,'final',final);
        fprintf('4 T reverse node %d h=%.9g: STOP final acceptance gate\n', ...
            inode,h(inode));
        break;
    end

    certified(inode) = true;
    Sigma_solution(:,inode) = probe.Sigma;
    fprintf(['4 T reverse node %d h=%.9g: certified, seed radius %.6g ' ...
        '(%s), Picard %d, residual %.3g\n'],inode,h(inode), ...
        linear_spectral_radius(inode),linear_gate_method(inode), ...
        probe.iterations,final.residual_norm);
end

if all(certified)
    overall_status = "completed_all_nodes";
    failed_node = NaN;
    failed_h = NaN;
else
    overall_status = "stopped_at_first_failure";
    failed_node = find(attempted & ~certified,1,'first');
    failed_h = h(failed_node);
end

node = (1:nH).';
tab = table(node,h,legacy_certified,legacy_static_status,attempted,certified, ...
    seed_source,seed_status,seed_residual,seed_static_roots,seed_D_uni, ...
    jac_status,jac_lambda,jac_eigen_residual,jac_iterations,probe_status, ...
    linear_gate_method,linear_gate_status,linear_spectral_radius, ...
    arnoldi_radius_spread,probe_iterations,final_residual,final_D_uni, ...
    final_supremum_mass, ...
    final_mesh_x_margin,final_mesh_medium_margin,final_dynamic_min_abs, ...
    final_dynamic_nonpositive_count,stop_reason);
summary = table(nnz(legacy_certified),nnz(certified), ...
    nnz(certified & ~legacy_certified),failed_node,failed_h,overall_status, ...
    'VariableNames',{'legacy_certified_nodes','reverse_certified_nodes', ...
    'newly_certified_nodes','failed_node','failed_h','status'});
result = struct('table',tab,'summary',summary, ...
    'Sigma_solution',Sigma_solution,'arnoldi_details',{arnoldi_details}, ...
    'stop_detail',stop_detail, ...
    'config',struct('T',T,'Bx',Bx,'Ecut',Ecut,'mapopts',mapopts, ...
                    'jacopts',jacopts,'probeopts',probeopts), ...
    'provenance',struct('legacy_fixture',fixture_path, ...
                        'failed_profile_fixture',failed_path, ...
                        'coupling_opts',F.provenance.coupling_opts), ...
    'note',['Diagnostic branch/component evidence only. It does not evaluate ' ...
            'the H_MF integral, prove root uniqueness outside the bounded ' ...
            'static scan, or select thermodynamic equilibrium.']);
save(fullfile(here,'wp2_4t_reverse_continuation.mat'),'result','-v7');
disp(summary);
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
