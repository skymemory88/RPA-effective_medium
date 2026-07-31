function result = invzp_highfield_boundary_continuation_ladder()
%INVZP_HIGHFIELD_BOUNDARY_CONTINUATION_LADDER Nested ordered-side B ladder.
% Diagnostic only. At each nH, two independently accepted lower-field states
% seed every target directly. Recovered targets are never reused as seeds.
% This separates numerical reach from the missing-area prescription and does
% not authorize an ordered-boundary production retry.

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
source_B = driver_fields(target_index-2:target_index-1);
target_B = [4.64 driver_fields(target_index) 4.68 4.69 4.70 4.71];
nH_values = [129 257];

[source_resolution_grid,source_field_grid] = ndgrid(nH_values,source_B);
source_points = cell(size(source_resolution_grid));
parfor k = 1:numel(source_points)
    opts = make_opts(J0eff,Jxx0,source_resolution_grid(k));
    source_points{k} = invz_solve_point_ordered( ...
        ion,T,source_field_grid(k),J,opts);
end
assert(all(cellfun(@(p)p.converged && p.is_ordered,source_points),'all'), ...
    'Every lower-field source must be independently accepted at both nH.');

[target_grid,seed_grid,resolution_grid] = ndgrid(target_B,source_B,nH_values);
points = cell(size(target_grid));
seed_states = cell(size(target_grid));
for k = 1:numel(target_grid)
    ir = find(nH_values == resolution_grid(k),1);
    is = find(source_B == seed_grid(k),1);
    source = source_points{ir,is};
    seed_states{k} = struct('Sigma',source.Sigma, ...
        'lambda',source.lambda,'K0',source.K(1));
end
parfor k = 1:numel(target_grid)
    opts = make_opts(J0eff,Jxx0,resolution_grid(k));
    opts.hmf_profile_state_seed = seed_states{k};
    opts.hmf_profile_sweep_direction = 'descending';
    points{k} = invz_solve_point_ordered( ...
        ion,T,target_grid(k),J,opts);
end

records = repmat(empty_record(),numel(points),1);
for k = 1:numel(points)
    records(k) = summarize(target_grid(k),seed_grid(k), ...
        resolution_grid(k),points{k},J,J0eff);
end
summary = struct2table(records);

agreement_records = repmat(empty_agreement(), ...
    numel(target_B)*numel(nH_values),1);
ia = 0;
for ir = 1:numel(nH_values)
    for ib = 1:numel(target_B)
        ia = ia+1;
        p1 = points{ib,1,ir};
        p2 = points{ib,2,ir};
        r1 = records(sub2ind(size(points),ib,1,ir));
        r2 = records(sub2ind(size(points),ib,2,ir));
        both = p1.converged && p2.converged;
        hmf_delta = NaN;
        Sigma_delta = NaN;
        D_uni_delta = NaN;
        if both
            hmf_delta = abs(p1.hmf-p2.hmf);
            Sigma_delta = max(abs(p1.Sigma-p2.Sigma));
            D_uni_delta = abs(p1.D_uni-p2.D_uni);
        end
        agrees = both && hmf_delta <= 1e-6 && ...
            Sigma_delta <= 1e-7 && D_uni_delta <= 1e-6 && ...
            r1.root_bracket_count == r2.root_bracket_count && ...
            r1.component_nodes == r2.component_nodes;
        agreement_records(ia) = struct('B',target_B(ib), ...
            'nH',nH_values(ir),'both_converged',both,'agrees',agrees, ...
            'hmf_delta',hmf_delta,'Sigma_delta',Sigma_delta, ...
            'D_uni_delta',D_uni_delta, ...
            'component_nodes_delta',abs(r1.component_nodes-r2.component_nodes), ...
            'support_lo_delta',abs(r1.support_lo-r2.support_lo), ...
            'support_hi_delta',abs(r1.support_hi-r2.support_hi));
    end
end
agreement = struct2table(agreement_records);

cease_field = nan(size(nH_values));
for ir = 1:numel(nH_values)
    rows = agreement.nH == nH_values(ir);
    first = find(~agreement.agrees(rows),1);
    fields = agreement.B(rows);
    if ~isempty(first), cease_field(ir) = fields(first); end
end
cease_table = table(nH_values(:),cease_field(:), ...
    'VariableNames',{'nH','first_nonagreement_B'});

pm_converged = false(size(target_B));
pm_mass = nan(size(target_B));
parfor k = 1:numel(target_B)
    pm = invz_solve_point(ion,T,target_B(k),J, ...
        struct('J0eff',J0eff,'Jxx0',Jxx0, ...
        'mix_outer',0.25,'max_outer',1000));
    pm_converged(k) = pm.converged;
    pm_mass(k) = pm.crit;
end
pm_table = table(target_B(:),pm_converged(:),pm_mass(:), ...
    'VariableNames',{'B_T','converged','pm_mass'});

accepted = summary.converged;
assert(all(summary.D_uni(accepted) > 0) && ...
    all(summary.supremum_mass(accepted) > 0) && ...
    all(summary.min_mesh_x_mass(accepted) > 0) && ...
    all(summary.min_mesh_medium_mass(accepted) > 0), ...
    'Every accepted continued root must retain positive static masses.');
assert(all(summary.Fprime(accepted) > 0), ...
    'Every accepted continued root must remain away from a zero of F''.');
assert(all(summary.final_resid(accepted) < 1e-8), ...
    'Every accepted continued root must pass the outer residual gate.');

result = struct('T',T,'area_factor',1,'nH_values',nH_values, ...
    'source_B',source_B,'target_B',target_B,'summary',summary, ...
    'agreement',agreement,'cease_table',cease_table,'pm_table',pm_table, ...
    'known_refined_pm_mass_zero_T',4.71897990927, ...
    'seed_policy',['Every target uses one of two independently accepted ' ...
    'lower-field cold states directly; recovered targets are not seeds.'], ...
    'production_changed',false);
save(fullfile(here,'wp6_highfield_boundary_continuation_ladder.mat'), ...
    'result','-v7');
disp(summary(:,{'B','seed_B','nH','status','converged','hmf','D_uni', ...
    'Fprime','final_resid','component_nodes','support_lo','support_hi'}));
disp(agreement);
disp(cease_table);
disp(pm_table);
end

function opts = make_opts(J0eff,Jxx0,nH)
opts = struct('J0eff',J0eff,'Jxx0',Jxx0,'mix_outer',0.25, ...
    'max_outer',1000,'nH',nH, ...
    'hmf_integral_mode','missing_area_approx', ...
    'hmf_missing_area_factor',1, ...
    'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
end

function out = summarize(B,seed_B,nH,point,J,J0eff)
out = empty_record();
out.B = B;
out.seed_B = seed_B;
out.nH = nH;
out.status = string(point.hmf_status);
out.converged = point.converged;
out.hmf = point.hmf;
out.D_uni = point.D_uni;
out.final_resid = point.final_resid;
p = point.hmf_prof;
out.component_nodes = p.missing_area_integral.node_count;
out.component_edge = p.missing_area_integral.component_edge;
out.root_bracket_count = p.root_bracket_count;
if out.component_nodes >= 2 && isfinite(out.component_edge) && ...
        p.missing_area_integral.edge_r > 0
    ix = p.missing_area_integral.selected_indices(:).';
    upper_integral = cumtrapz(p.hgrid(ix),p.r(ix));
    required_factor = (J0eff*p.m(ix)-upper_integral)/ ...
        (out.component_edge*p.missing_area_integral.edge_r);
    out.support_lo = min(required_factor);
    out.support_hi = max(required_factor);
end
if ~(point.converged && isfinite(p.r_star) && ...
        isfinite(p.Gstat_star) && isfinite(point.K(1)))
    return;
end
Gstat = p.Gstat_star;
K0 = point.K(1);
Gtilde = Gstat/(1-K0*Gstat);
out.Gtilde = Gtilde;
out.Fprime = p.r_star*(1+J0eff*Gtilde);
out.supremum_mass = 1+J0eff*Gtilde;
out.min_mesh_x_mass = min(1+J(:)*Gtilde);
out.min_mesh_medium_mass = min(1+(J(:)-K0)*Gstat);
end

function out = empty_record()
out = struct('B',NaN,'seed_B',NaN,'nH',NaN,'status',"", ...
    'converged',false,'hmf',NaN,'D_uni',NaN,'Fprime',NaN, ...
    'final_resid',NaN,'component_nodes',0,'component_edge',NaN, ...
    'root_bracket_count',0,'support_lo',NaN,'support_hi',NaN, ...
    'Gtilde',NaN,'supremum_mass',NaN,'min_mesh_x_mass',NaN, ...
    'min_mesh_medium_mass',NaN);
end

function out = empty_agreement()
out = struct('B',NaN,'nH',NaN,'both_converged',false,'agrees',false, ...
    'hmf_delta',NaN,'Sigma_delta',NaN,'D_uni_delta',NaN, ...
    'component_nodes_delta',NaN,'support_lo_delta',NaN, ...
    'support_hi_delta',NaN);
end
