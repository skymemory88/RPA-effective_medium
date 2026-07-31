function result = invzp_highfield_derivative_section_audit()
%INVZP_HIGHFIELD_DERIVATIVE_SECTION_AUDIT Resolve the near-Bc F' discrepancy.
% Diagnostic only. Reconstruct a narrow fixed-(B,h) section around the
% continued factor-one roots and compare
%   dm/dh = -G0bare,
%   F' = r + J0eff*G0bare = r*(1+J0eff*Gtilde0),
% with local and retained coarse bracket secants. Gstat is recorded to show
% why substituting it for Gtilde0 is not the derivative identity.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));
F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
ion = invz_ion();
T = 0.1;
J = F.J;
J0eff = F.info.Jcc0;
Jxx0 = F.info.Jaa0;
driver_fields = linspace(0,9,111);
[~,target_index] = min(abs(driver_fields-4.68));
source_B = driver_fields(target_index-1);
target_B = [driver_fields(target_index) 4.68];
opts = make_opts(J0eff,Jxx0);

source = invz_solve_point_ordered(ion,T,source_B,J,opts);
assert(source.converged && source.is_ordered, ...
    'The independently accepted ordered-side source must converge.');
source_seed = struct('Sigma',source.Sigma,'lambda',source.lambda, ...
    'K0',source.K(1));

points = cell(size(target_B));
sections = cell(size(target_B));
summary_records = repmat(empty_summary(),numel(target_B),1);
for ib = 1:numel(target_B)
    seeded_opts = opts;
    seeded_opts.hmf_profile_state_seed = source_seed;
    seeded_opts.hmf_profile_sweep_direction = 'descending';
    point = invz_solve_point_ordered(ion,T,target_B(ib),J,seeded_opts);
    assert(point.converged && point.is_ordered, ...
        'Continued factor-one target %.12g T did not converge.',target_B(ib));
    points{ib} = point;

    p = point.hmf_prof;
    ia = p.root_bracket_indices(1);
    ic = p.root_bracket_indices(2);
    assert(all(isfinite([ia ic])) && ic == ia+1, ...
        'Expected one adjacent retained root bracket.');
    coarse_width = p.hgrid(ic)-p.hgrid(ia);
    section_step = min(2e-5,0.1*coarse_width);
    hsection = point.hmf + (-2:2)*section_step;
    fixed_seed = struct('Sigma',point.Sigma,'lambda',point.lambda, ...
        'K0',point.K(1));
    node_records = repmat(empty_node(),numel(hsection),1);
    for ih = 1:numel(hsection)
        q = eval_fixed_node(ion,T,target_B(ib),hsection(ih),J, ...
            seeded_opts,fixed_seed);
        fd_step = min(1e-7,section_step/20);
        mp = fixed_moment(ion,T,target_B(ib),hsection(ih)+fd_step,Jxx0);
        mm = fixed_moment(ion,T,target_B(ib),hsection(ih)-fd_step,Jxx0);
        q.dm_dh_fd = (mp-mm)/(2*fd_step);
        q.dm_response_error = q.dm_dh_fd+q.G0bare;
        q.B = target_B(ib);
        q.h_offset = hsection(ih)-point.hmf;
        node_records(ih) = q;
    end
    section = struct2table(node_records);
    assert(all(section.converged), ...
        'A narrow fixed-h section left the certified component.');

    section_integral = cumtrapz(section.h,section.r);
    center = 3;
    section.F_local = section_integral-section_integral(center) - ...
        J0eff*(section.m-section.m(center));
    section_secant = (section.F_local(center+1)- ...
        section.F_local(center-1))/(2*section_step);
    coarse_secant = (p.F(ic)-p.F(ia))/coarse_width;
    exact = section.Fprime_tilde(center);
    wrong = section.Fprime_Gstat(center);
    dm_scale = max(1,abs(section.G0bare));
    response_rel_error = max(abs(section.dm_response_error)./dm_scale);
    identity_error = max(abs(section.Fprime_direct-section.Fprime_tilde));
    section_rel_error = abs(section_secant-exact)/max(abs(exact),eps);
    coarse_rel_error = abs(coarse_secant-exact)/max(abs(exact),eps);

    summary_records(ib) = struct('B',target_B(ib), ...
        'source_B',source_B,'hmf',point.hmf,'section_step',section_step, ...
        'coarse_bracket_width',coarse_width,'D_uni',point.D_uni, ...
        'Fprime_exact',exact,'Fprime_Gstat_substitution',wrong, ...
        'section_secant',section_secant,'coarse_bracket_secant',coarse_secant, ...
        'response_rel_error',response_rel_error, ...
        'identity_abs_error',identity_error, ...
        'section_rel_error',section_rel_error, ...
        'coarse_rel_error',coarse_rel_error, ...
        'min_supremum_mass',min(section.supremum_mass), ...
        'min_mesh_x_mass',min(section.min_mesh_x_mass), ...
        'min_mesh_medium_mass',min(section.min_mesh_medium_mass), ...
        'max_outer_residual',max(section.outer_residual));
    sections{ib} = section;
end

summary = struct2table(summary_records);
assert(all(summary.Fprime_exact > 0) && ...
    all(summary.section_secant > 0) && ...
    all(summary.coarse_bracket_secant > 0), ...
    'Exact and sampled F derivatives must have the same positive sign.');
assert(all(summary.Fprime_Gstat_substitution < 0), ...
    'The retained Gstat substitution should reproduce the rejected sign.');
assert(all(summary.response_rel_error < 2e-6), ...
    'Finite-difference dm/dh does not agree with -G0bare.');
assert(all(summary.identity_abs_error < 1e-12), ...
    'The direct and Gtilde forms of F'' are not algebraically consistent.');
assert(all(summary.section_rel_error < 5e-3), ...
    'The narrow-section secant does not agree with the exact F'' identity.');
assert(all(summary.coarse_rel_error < 1.5e-1), ...
    'The retained root-bracket secant is inconsistent with the exact F'' identity.');

result = struct('T',T,'nH',opts.nH,'area_factor',1, ...
    'source_B',source_B,'target_B',target_B,'summary',summary, ...
    'sections',{sections}, ...
    'identity',['Fprime = r + J0eff*G0bare = ' ...
    'r*(1+J0eff*Gtilde0)'], ...
    'resolved_cause',['The prior diagnostic substituted Gstat for Gtilde0. ' ...
    'That substitution is not algebraically valid and reverses the sign ' ...
    'near the soft uniform boundary.'], ...
    'production_changed',false);
save(fullfile(here,'wp6_highfield_derivative_section_audit.mat'), ...
    'result','-v7');
disp(summary);
end

function opts = make_opts(J0eff,Jxx0)
opts = struct('J0eff',J0eff,'Jxx0',Jxx0,'mix_outer',0.25, ...
    'max_outer',1000,'nH',129, ...
    'hmf_integral_mode','missing_area_approx', ...
    'hmf_missing_area_factor',1, ...
    'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
end

function out = eval_fixed_node(ion,T,Bx,h,Jnu_flat,opts,seed)
J0eff = opts.J0eff;
Jxx0 = opts.Jxx0;
mixo = get_opt(opts,'mix_outer',0.7);
tolo = get_opt(opts,'tol_outer',1e-8);
maxo = get_opt(opts,'max_outer',200);
Ecut = get_opt(opts,'Ecut',40);
hyp = get_opt(opts,'hyp',true);
eopts = get_opt(opts,'emt',struct());
eso = get_opt(opts,'emt_static',struct());
eso.warn = false;
[wn,wts,beta] = invz_matsubara(T,Ecut);

si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',hyp,'Jxx0',Jxx0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',Jxx0));
m = si.Jexp(3);
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,T,1i*wn(1),struct('elastic',false));
G0inel0 = -real(c0i(3,3,1));
X = real(c0(:,:,1));
feedback = X(3,1)*(Jxx0/(1-Jxx0*X(1,1)))*X(1,3);
G0bare = -(X(3,3)+feedback);
G0el0 = G0bare-G0inel0;
g = real(invz_g(tl,1i*wn));

Sigma = seed.Sigma(:);
K0 = seed.K0;
lambda = seed.lambda(:);
K = zeros(size(wn));
outer_converged = false;
med = struct('converged',false);
for outer = 1:maxo
    eopts.K0 = K;
    med = invz_emt_scalar(G0,Sigma,Jnu_flat,eopts);
    K = med.K;
    [K0,~,so] = invz_emt_static_ordered(tl,lambda(1:2),Sigma(1), ...
        Jnu_flat,K0,beta,J0eff,G0inel0,G0el0,eso);
    if ~so.converged
        break;
    end
    K(1) = K0;
    lambda = invz_lambdas(K,g,wts,beta,[1 2 3]);
    sg = invz_sigma_ordered(tl,lambda,K,g,beta);
    dSigma = max(abs(sg.Sigma-Sigma));
    Sigma = Sigma+mixo*(sg.Sigma-Sigma);
    if dSigma < tolo
        outer_converged = true;
        break;
    end
end
[K0,Gstat,so] = invz_emt_static_ordered(tl,lambda(1:2),Sigma(1), ...
    Jnu_flat,K0,beta,J0eff,G0inel0,G0el0,eso);
K(1) = K0;
lambda_check = invz_lambdas(K,g,wts,beta,[1 2 3]);
sigma_check = invz_sigma_ordered(tl,lambda_check,K,g,beta);
outer_residual = max(abs(sigma_check.Sigma-Sigma));
ctol = get_opt(eso,'resid_tol',1e-10);
converged = outer_converged && med.converged && so.converged && ...
    isfinite(so.resid) && so.resid < ctol && outer_residual < tolo;

min_mesh_x_mass = NaN;
min_mesh_medium_mass = NaN;
supremum_mass = NaN;
if so.converged
    row = so.root_table(so.selected_index,:);
    min_mesh_x_mass = row.min_mesh_x_signed;
    min_mesh_medium_mass = row.min_mesh_medium_signed;
    supremum_mass = row.supremum_mass;
end
Fprime_direct = so.r+J0eff*G0bare;
Fprime_tilde = so.r*(1+J0eff*so.Gtil0);
Fprime_Gstat = so.r*(1+J0eff*Gstat);
out = struct('B',NaN,'h',h,'h_offset',NaN,'converged',converged, ...
    'outer_iterations',outer,'outer_residual',outer_residual, ...
    'm',m,'dm_dh_fd',NaN,'G0bare',G0bare,'dm_response_error',NaN, ...
    'Gstat',Gstat,'K0',K0,'Gtilde0',so.Gtil0,'r',so.r, ...
    'D_uni',so.D_uni,'supremum_mass',supremum_mass, ...
    'min_mesh_x_mass',min_mesh_x_mass, ...
    'min_mesh_medium_mass',min_mesh_medium_mass, ...
    'Fprime_direct',Fprime_direct,'Fprime_tilde',Fprime_tilde, ...
    'Fprime_Gstat',Fprime_Gstat);
end

function m = fixed_moment(ion,T,Bx,h,Jxx0)
si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',Jxx0,'hz_fixed',h));
m = si.Jexp(3);
end

function out = empty_node()
out = struct('B',NaN,'h',NaN,'h_offset',NaN,'converged',false, ...
    'outer_iterations',NaN,'outer_residual',NaN,'m',NaN, ...
    'dm_dh_fd',NaN,'G0bare',NaN,'dm_response_error',NaN, ...
    'Gstat',NaN,'K0',NaN,'Gtilde0',NaN,'r',NaN,'D_uni',NaN, ...
    'supremum_mass',NaN,'min_mesh_x_mass',NaN, ...
    'min_mesh_medium_mass',NaN,'Fprime_direct',NaN, ...
    'Fprime_tilde',NaN,'Fprime_Gstat',NaN);
end

function out = empty_summary()
out = struct('B',NaN,'source_B',NaN,'hmf',NaN,'section_step',NaN, ...
    'coarse_bracket_width',NaN,'D_uni',NaN,'Fprime_exact',NaN, ...
    'Fprime_Gstat_substitution',NaN,'section_secant',NaN, ...
    'coarse_bracket_secant',NaN,'response_rel_error',NaN, ...
    'identity_abs_error',NaN,'section_rel_error',NaN, ...
    'coarse_rel_error',NaN,'min_supremum_mass',NaN, ...
    'min_mesh_x_mass',NaN,'min_mesh_medium_mass',NaN, ...
    'max_outer_residual',NaN);
end

function v = get_opt(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
