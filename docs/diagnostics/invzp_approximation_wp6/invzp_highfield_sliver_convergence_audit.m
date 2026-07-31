function result = invzp_highfield_sliver_convergence_audit()
%INVZP_HIGHFIELD_SLIVER_CONVERGENCE_AUDIT Classify the moving Bc-side mask.
% Compare fixed-area production members, supported intermediate areas,
% resolution/damping changes, and one-sided ordered-state continuation.  No
% diagnostic result changes the production area prescription or its gates.

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
target_B = driver_fields(target_index);
left_B = driver_fields(target_index-2:target_index-1);
right_B = driver_fields(target_index+1);

name = ["grid_4p5_factor1";"grid_left_factor1"; ...
    "grid_target_factor1";"fixed_4p68_factor1"; ...
    "fixed_4p70_factor1";"grid_right_factor1"; ...
    "target_factor1_n257";"target_factor1_damped"; ...
    "target_factor085_n129";"target_factor090_n129"; ...
    "target_factor095_n129";"target_factor090_n257"];
B = [left_B(1);left_B(2);target_B;4.68;4.70;right_B; ...
    target_B;target_B;target_B;target_B;target_B;target_B];
nH = [129;129;129;129;129;129;257;129;129;129;129;257];
mix_outer = [0.25;0.25;0.25;0.25;0.25;0.25;0.25;0.15; ...
    0.25;0.25;0.25;0.25];
max_outer = [1000;1000;1000;1000;1000;1000;1000;3000; ...
    1000;1000;1000;1000];
area_factor = [1;1;1;1;1;1;1;1;0.85;0.90;0.95;0.90];

points = cell(size(B));
parfor k = 1:numel(B)
    opts = make_opts(J0eff,Jxx0,nH(k),mix_outer(k), ...
        max_outer(k),area_factor(k));
    points{k} = invz_solve_point_ordered(ion,T,B(k),J,opts);
end
records = repmat(empty_record(),numel(B),1);
for k = 1:numel(B)
    records(k) = summarize_point(name(k),B(k),nH(k),mix_outer(k), ...
        max_outer(k),area_factor(k),NaN,points{k},J0eff);
end

% The last two accepted cold ordered grid points are independent same-side
% sources. This is diagnostic only: unlike the production retry, no accepted
% ordered source exists above the target field.
seed_source_index = [1 2];
seeded_points = cell(1,numel(seed_source_index));
for k = 1:numel(seed_source_index)
    source = points{seed_source_index(k)};
    if ~(source.converged && source.is_ordered)
        continue;
    end
    opts = make_opts(J0eff,Jxx0,129,0.25,1000,1);
    opts.hmf_profile_state_seed = struct('Sigma',source.Sigma, ...
        'lambda',source.lambda,'K0',source.K(1));
    opts.hmf_profile_sweep_direction = 'descending';
    seeded_points{k} = invz_solve_point_ordered(ion,T,target_B,J,opts);
end
seeded_records = repmat(empty_record(),numel(seed_source_index),1);
for k = 1:numel(seed_source_index)
    if isempty(seeded_points{k})
        seeded_records(k).name = "source_not_accepted";
        seeded_records(k).B = target_B;
        seeded_records(k).seed_B = B(seed_source_index(k));
    else
        seeded_records(k) = summarize_point( ...
            "target_factor1_seeded",target_B,129,0.25,1000,1, ...
            B(seed_source_index(k)),seeded_points{k},J0eff);
    end
end

pm_fields = [left_B target_B 4.68 4.70 right_B];
pm_crit = nan(size(pm_fields));
pm_converged = false(size(pm_fields));
parfor k = 1:numel(pm_fields)
    popts = struct('J0eff',J0eff,'Jxx0',Jxx0, ...
        'mix_outer',0.25,'max_outer',1000);
    pm = invz_solve_point(ion,T,pm_fields(k),J,popts);
    pm_crit(k) = pm.crit;
    pm_converged(k) = pm.converged;
end
pm_table = table(pm_fields(:),pm_converged(:),pm_crit(:), ...
    'VariableNames',{'B_T','converged','pm_mass'});

summary = struct2table(records);
seeded_summary = struct2table(seeded_records);
result = struct('T',T,'driver_field_count',numel(driver_fields), ...
    'driver_field_step',driver_fields(2)-driver_fields(1), ...
    'target_index',target_index,'target_B',target_B, ...
    'left_B',left_B,'right_B',right_B,'summary',summary, ...
    'seeded_summary',seeded_summary,'pm_table',pm_table, ...
    'known_refined_pm_mass_zero_T',4.71897990927, ...
    'interpretation',['The moving high-field sliver is classified by ' ...
    'separating nonlinear convergence from missing-area root support. ' ...
    'Intermediate area factors are diagnostic prescriptions, not a ' ...
    'thermodynamic calibration, and same-side seeds are not production-' ...
    'accepted without an additional boundary-continuation policy.']);
save(fullfile(here,'wp6_highfield_sliver_convergence_audit.mat'), ...
    'result','-v7');

disp(summary(:,{'name','B','nH','mix_outer','area_factor','seed_B', ...
    'status','converged','hmf','D_uni','final_resid','component_nodes', ...
    'support_lo','support_hi'}));
disp(seeded_summary(:,{'name','B','area_factor','seed_B','status', ...
    'converged','hmf','component_nodes','support_lo','support_hi'}));
disp(pm_table);
end

function opts = make_opts(J0eff,Jxx0,nH,mix_outer,max_outer,area_factor)
opts = struct('J0eff',J0eff,'Jxx0',Jxx0,'mix_outer',mix_outer, ...
    'max_outer',max_outer,'nH',nH, ...
    'hmf_integral_mode','missing_area_approx', ...
    'hmf_missing_area_factor',area_factor, ...
    'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
end

function out = summarize_point(name,B,nH,mix_outer,max_outer, ...
        area_factor,seed_B,point,J0eff)
out = empty_record();
out.name = name;
out.B = B;
out.nH = nH;
out.mix_outer = mix_outer;
out.max_outer = max_outer;
out.area_factor = area_factor;
out.seed_B = seed_B;
out.status = string(point.hmf_status);
out.converged = point.converged;
out.hmf = point.hmf;
if isfield(point,'D_uni'), out.D_uni = point.D_uni; end
if isfield(point,'final_resid'), out.final_resid = point.final_resid; end
p = point.hmf_prof;
if isfield(p,'r_star') && isfinite(p.r_star) && ...
        isfield(p,'Gstat_star') && isfinite(p.Gstat_star) && ...
        isfield(point,'K') && isfinite(point.K(1))
    Gtilde_star = p.Gstat_star/(1-point.K(1)*p.Gstat_star);
    out.Fprime = p.r_star*(1+J0eff*Gtilde_star);
end
if ~(isfield(p,'missing_area_integral') && ...
        isfield(p.missing_area_integral,'node_count'))
    return;
end
m = p.missing_area_integral;
out.component_nodes = m.node_count;
out.component_edge = m.component_edge;
if m.node_count < 2 || ~isfinite(m.component_edge) || ...
        ~isfinite(m.edge_r) || m.edge_r <= 0
    return;
end
ix = m.selected_indices(:).';
upper_integral = cumtrapz(p.hgrid(ix),p.r(ix));
required_factor = (J0eff*p.m(ix)-upper_integral)/ ...
    (m.component_edge*m.edge_r);
out.support_lo = min(required_factor);
out.support_hi = max(required_factor);
end

function out = empty_record()
out = struct('name',"",'B',NaN,'nH',NaN,'mix_outer',NaN, ...
    'max_outer',NaN,'area_factor',NaN,'seed_B',NaN,'status',"", ...
    'converged',false,'hmf',NaN,'D_uni',NaN,'Fprime',NaN, ...
    'final_resid',NaN,'component_nodes',0,'component_edge',NaN, ...
    'support_lo',NaN,'support_hi',NaN);
end
