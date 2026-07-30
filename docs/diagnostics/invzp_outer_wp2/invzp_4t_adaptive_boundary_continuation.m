function result = invzp_4t_adaptive_boundary_continuation()
%INVZP_4T_ADAPTIVE_BOUNDARY_CONTINUATION Resolve the node-29 -> node-28 gap.
% Diagnostic only. The controller keeps the deterministic bounded outer map
% and undamped Picard fixed while adapting only the downward H_MF step:
% rejected proposals halve the step; accepted proposals grow it by 1.5 up to
% the remaining target distance. Every accepted root passes the same static,
% leading-mode, residual, and dynamic-denominator gates as the coarse reverse
% continuation. No H_MF integral or production output is evaluated.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
reverse_path = fullfile(here,'wp2_4t_reverse_continuation.mat');
F = load(fixture_path);
R = load(reverse_path);

start_node = 29;
target_node = 28;
h_start = R.result.table.h(start_node);
h_target = R.result.table.h(target_node);
Sigma_current = R.result.Sigma_solution(:,start_node);
if ~(R.result.table.certified(start_node) && ...
        all(isfinite(Sigma_current)) && h_start > h_target)
    error('invz:adaptiveContinuation', ...
        'the retained reverse-continuation artifact has no certified start.');
end

ion = invz_ion();
T = R.result.config.T;
Bx = R.result.config.Bx;
Ecut = R.result.config.Ecut;
[wn,wts,beta] = invz_matsubara(T,Ecut);
mapopts = struct('emt_static', ...
    struct('Jsup',F.info.Jcc0,'warn',false,'scan_points',4097, ...
           'endpoint_margin',1e-10,'resid_tol',1e-10));
probeopts = struct('mix',1,'tol',1e-8,'maxit',200);
controller = struct('initial_step',h_target-h_start,'step_growth',1.5, ...
    'min_step',(h_start-h_target)/2^14,'h_tolerance',1e-12, ...
    'max_attempts',60);

start_ctx = make_context(ion,T,Bx,h_start,wn,wts,beta,F.J,F.info);
start_map = invz_ordered_outer_map(Sigma_current,start_ctx,mapopts);
if ~(start_map.defined && start_map.residual_norm <= probeopts.tol)
    error('invz:adaptiveContinuation', ...
        'the certified starting Sigma no longer passes its outer-map gate.');
end

nmax = controller.max_attempts;
attempt = (1:nmax).';
from_h = nan(nmax,1);
candidate_h = nan(nmax,1);
step_size = nan(nmax,1);
base_status = repmat("not_attempted",nmax,1);
base_residual = nan(nmax,1);
base_static_roots = nan(nmax,1);
base_min_sample_abs_residual = nan(nmax,1);
linear_method = repmat("not_attempted",nmax,1);
linear_status = repmat("not_attempted",nmax,1);
linear_spectral_radius = nan(nmax,1);
linear_radius_spread = nan(nmax,1);
probe_status = repmat("not_attempted",nmax,1);
probe_iterations = nan(nmax,1);
accepted = false(nmax,1);
final_residual = nan(nmax,1);
final_D_uni = nan(nmax,1);
final_supremum_mass = nan(nmax,1);
final_mesh_x_margin = nan(nmax,1);
final_mesh_medium_margin = nan(nmax,1);
final_dynamic_min_abs = nan(nmax,1);
final_dynamic_nonpositive_count = nan(nmax,1);
outcome = repmat("not_attempted",nmax,1);
details = cell(nmax,1);

accepted_h = nan(nmax+1,1);
accepted_Sigma = nan(numel(wn),nmax+1);
accepted_h(1) = h_start;
accepted_Sigma(:,1) = Sigma_current;
naccepted = 1;
h_current = h_start;
step = controller.initial_step;
overall_status = "max_attempts";
nattempt = 0;

for k = 1:nmax
    nattempt = k;
    remaining = h_current-h_target;
    if remaining <= controller.h_tolerance
        overall_status = "reached_target";
        nattempt = k-1;
        break;
    end
    step = -min(abs(step),remaining);
    h_candidate = h_current+step;
    from_h(k) = h_current;
    candidate_h(k) = h_candidate;
    step_size(k) = step;

    ctx = make_context(ion,T,Bx,h_candidate,wn,wts,beta,F.J,F.info);
    mapfun = @(s) invz_ordered_outer_map(s,ctx,mapopts);
    base = mapfun(Sigma_current);
    base_status(k) = string(base.status);
    base_residual(k) = base.residual_norm;
    if ~isempty(base.static)
        base_static_roots(k) = base.static.n_admissible_roots;
        fg = abs(base.static.search.fgrid( ...
            isfinite(base.static.search.fgrid)));
        if ~isempty(fg), base_min_sample_abs_residual(k) = min(fg); end
    end

    if ~base.defined
        outcome(k) = "seed_map_" + string(base.status);
        details{k} = struct('stage',"seed_map",'linear',[], ...
            'probe',[],'final',[]);
        [step,overall_status] = reject_step(step,controller);
        fprintf(['adaptive attempt %d: %.9g -> %.9g rejected at seed map ' ...
            '(%s), next |dh| %.3g\n'],k,h_current,h_candidate, ...
            base.status,abs(step));
        if overall_status ~= "continue", break; end
        continue;
    end

    linear = measure_linear_gate(mapfun,Sigma_current);
    details{k} = struct('stage',"linear_gate",'linear',linear, ...
        'probe',[],'final',[]);
    linear_method(k) = linear.method;
    linear_status(k) = linear.status;
    linear_spectral_radius(k) = linear.spectral_radius;
    linear_radius_spread(k) = linear.radius_spread;
    if ~(linear.status == "ok" && isfinite(linear.spectral_radius) && ...
            linear.spectral_radius < 1)
        outcome(k) = "linear_gate_" + linear.status;
        [step,overall_status] = reject_step(step,controller);
        fprintf(['adaptive attempt %d: %.9g -> %.9g rejected by %s ' ...
            '(radius %.6g), next |dh| %.3g\n'],k,h_current,h_candidate, ...
            linear.status,linear.spectral_radius,abs(step));
        if overall_status ~= "continue", break; end
        continue;
    end

    probe = invz_outer_picard_diagnostic(mapfun,Sigma_current,probeopts);
    probe_status(k) = string(probe.status);
    probe_iterations(k) = probe.iterations;
    details{k}.stage = "picard";
    details{k}.probe = compact_probe(probe);
    if ~probe.converged
        outcome(k) = "picard_" + string(probe.status);
        [step,overall_status] = reject_step(step,controller);
        fprintf(['adaptive attempt %d: %.9g -> %.9g rejected by Picard ' ...
            '(%s), next |dh| %.3g\n'],k,h_current,h_candidate, ...
            probe.status,abs(step));
        if overall_status ~= "continue", break; end
        continue;
    end

    final = invz_ordered_outer_map(probe.Sigma,ctx, ...
        struct('emt_static',mapopts.emt_static,'dynamic_diagnostics',true));
    final_residual(k) = final.residual_norm;
    if final.defined
        row = final.static.selected_index;
        rt = final.static.root_table;
        final_D_uni(k) = final.static.D_uni;
        final_supremum_mass(k) = rt.supremum_mass(row);
        final_mesh_x_margin(k) = rt.min_mesh_x_signed(row);
        final_mesh_medium_margin(k) = rt.min_mesh_medium_signed(row);
        final_dynamic_min_abs(k) = final.dynamic_min_abs;
        final_dynamic_nonpositive_count(k) = ...
            final.dynamic_nonpositive_count;
    end
    gate_ok = final.defined && final.static.n_admissible_roots == 1 && ...
        isfinite(final.residual_norm) && final.residual_norm <= probeopts.tol && ...
        all(isfinite(probe.Sigma)) && final.static.D_uni > 0 && ...
        final_supremum_mass(k) > 0 && final_mesh_x_margin(k) > 0 && ...
        final_mesh_medium_margin(k) > 0 && ...
        final.dynamic_nonpositive_count == 0 && ...
        isfinite(final.dynamic_min_abs) && final.dynamic_min_abs > 0;
    details{k}.stage = "final_gate";
    details{k}.final = compact_final(final);
    if ~gate_ok
        outcome(k) = "final_acceptance_gate";
        [step,overall_status] = reject_step(step,controller);
        fprintf(['adaptive attempt %d: %.9g -> %.9g rejected at final gate, ' ...
            'next |dh| %.3g\n'],k,h_current,h_candidate,abs(step));
        if overall_status ~= "continue", break; end
        continue;
    end

    accepted(k) = true;
    outcome(k) = "accepted";
    h_current = h_candidate;
    Sigma_current = probe.Sigma;
    naccepted = naccepted+1;
    accepted_h(naccepted) = h_current;
    accepted_Sigma(:,naccepted) = Sigma_current;
    fprintf(['adaptive attempt %d: accepted h=%.9g, |dh| %.3g, ' ...
        'radius %.6g, Picard %d, residual %.3g\n'],k,h_current, ...
        abs(step),linear.spectral_radius,probe.iterations,final.residual_norm);

    remaining = h_current-h_target;
    if remaining <= controller.h_tolerance
        overall_status = "reached_target";
        break;
    end
    step = -min(controller.step_growth*abs(step),remaining);
end

if overall_status == "continue"
    overall_status = "max_attempts";
end
attempt = attempt(1:nattempt);
from_h = from_h(1:nattempt);
candidate_h = candidate_h(1:nattempt);
step_size = step_size(1:nattempt);
base_status = base_status(1:nattempt);
base_residual = base_residual(1:nattempt);
base_static_roots = base_static_roots(1:nattempt);
base_min_sample_abs_residual = base_min_sample_abs_residual(1:nattempt);
linear_method = linear_method(1:nattempt);
linear_status = linear_status(1:nattempt);
linear_spectral_radius = linear_spectral_radius(1:nattempt);
linear_radius_spread = linear_radius_spread(1:nattempt);
probe_status = probe_status(1:nattempt);
probe_iterations = probe_iterations(1:nattempt);
accepted = accepted(1:nattempt);
final_residual = final_residual(1:nattempt);
final_D_uni = final_D_uni(1:nattempt);
final_supremum_mass = final_supremum_mass(1:nattempt);
final_mesh_x_margin = final_mesh_x_margin(1:nattempt);
final_mesh_medium_margin = final_mesh_medium_margin(1:nattempt);
final_dynamic_min_abs = final_dynamic_min_abs(1:nattempt);
final_dynamic_nonpositive_count = ...
    final_dynamic_nonpositive_count(1:nattempt);
outcome = outcome(1:nattempt);
details = details(1:nattempt);
tab = table(attempt,from_h,candidate_h,step_size,base_status, ...
    base_residual,base_static_roots,base_min_sample_abs_residual, ...
    linear_method,linear_status,linear_spectral_radius,linear_radius_spread, ...
    probe_status,probe_iterations,accepted,final_residual,final_D_uni, ...
    final_supremum_mass,final_mesh_x_margin,final_mesh_medium_margin, ...
    final_dynamic_min_abs,final_dynamic_nonpositive_count,outcome);

accepted_h = accepted_h(1:naccepted);
accepted_Sigma = accepted_Sigma(:,1:naccepted);
target_reached = h_current-h_target <= controller.h_tolerance;
summary = table(overall_status,target_reached,nattempt,nnz(accepted), ...
    h_start,h_target,h_current,min(abs(step_size)), ...
    'VariableNames',{'status','target_reached','attempts', ...
    'accepted_steps','start_h','target_h','final_h','minimum_attempted_step'});
result = struct('table',tab,'summary',summary,'accepted_h',accepted_h, ...
    'accepted_Sigma',accepted_Sigma,'details',{details}, ...
    'start_map',compact_final(start_map),'config', ...
    struct('T',T,'Bx',Bx,'Ecut',Ecut,'start_node',start_node, ...
           'target_node',target_node,'mapopts',mapopts, ...
           'probeopts',probeopts,'controller',controller), ...
    'provenance',struct('fixture',fixture_path, ...
                        'reverse_continuation',reverse_path), ...
    'note',['Adaptive parameter continuation on one coupled component only. ' ...
            'Rejected proposals are step-size observations within this one ' ...
            'method; success does not select thermodynamic equilibrium.']);
save(fullfile(here,'wp2_4t_adaptive_boundary_continuation.mat'), ...
    'result','-v7');
disp(summary);
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

function out = measure_linear_gate(mapfun,Sigma)
power = invz_outer_dominant_eigen(mapfun,Sigma, ...
    struct('fd_step',3e-6,'tol',1e-5,'maxit',40));
out = struct('method',"power",'status',"unresolved", ...
    'spectral_radius',NaN,'radius_spread',NaN, ...
    'power',compact_power(power),'arnoldi',[]);
if power.converged && isfinite(power.lambda)
    out.status = "ok";
    out.spectral_radius = abs(power.lambda);
    return;
end

fd_steps = [1e-5;3e-6;1e-6];
runs = cell(numel(fd_steps),1);
radius = nan(size(fd_steps));
converged = false(size(fd_steps));
for k = 1:numel(fd_steps)
    q = invz_outer_arnoldi_diagnostic(mapfun,Sigma, ...
        struct('n_eigs',6,'fd_step',fd_steps(k),'tol',1e-5, ...
               'maxit',100,'subspace_dimension',20));
    runs{k} = compact_arnoldi(q);
    radius(k) = q.spectral_radius;
    converged(k) = q.converged;
end
spread = max(radius)-min(radius);
stable = all(isfinite(radius)) && ...
    spread <= max(1e-4,0.05*median(radius));
out.method = "arnoldi_fd_ladder";
out.radius_spread = spread;
out.arnoldi = struct('fd_steps',fd_steps,'runs',{runs}, ...
    'spectral_radius',radius,'radius_stable',stable);
if all(converged) && stable
    out.status = "ok";
    out.spectral_radius = max(radius);
elseif ~all(converged)
    out.status = "arnoldi_unresolved";
else
    out.status = "arnoldi_fd_unstable";
end
end

function [step,status] = reject_step(step,controller)
step = step/2;
if abs(step) < controller.min_step
    status = "minimum_step_reached";
else
    status = "continue";
end
end

function out = compact_power(q)
out = rmfield(q,'vector');
end

function out = compact_arnoldi(q)
out = rmfield(q,'vectors');
end

function out = compact_probe(q)
out = struct('status',q.status,'converged',q.converged, ...
    'iterations',q.iterations,'residual_history',q.residual_history, ...
    'map_status',q.map_status,'mix',q.mix,'tol',q.tol);
end

function out = compact_final(q)
out = struct('status',q.status,'defined',q.defined, ...
    'residual_norm',q.residual_norm,'Sigma_map',q.Sigma_map, ...
    'K',q.K,'lambda',q.lambda,'G',q.G, ...
    'lambda_consistency',q.lambda_consistency, ...
    'dynamic_min_abs',q.dynamic_min_abs, ...
    'dynamic_nonpositive_count',q.dynamic_nonpositive_count);
if ~isempty(q.static)
    out.static_status = q.static.status;
    out.static_D_uni = q.static.D_uni;
    out.static_n_admissible_roots = q.static.n_admissible_roots;
else
    out.static_status = "not_evaluated";
    out.static_D_uni = NaN;
    out.static_n_admissible_roots = NaN;
end
end
