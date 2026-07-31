function result = invzp_highfield_progressive_field_continuation()
%INVZP_HIGHFIELD_PROGRESSIVE_FIELD_CONTINUATION Trace the Bc-side HMF root.
% Diagnostic only.  Two independently accepted ordered states are advanced
% through the same fine field ladder.  An accepted target seeds only the next
% target on its own path.  This deliberately tests whether the direct-seed
% failure near 4.70 T is a continuation-distance problem; recovered states
% remain ineligible as production retry seeds.

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
source_B = [4.50 4.59];
target_B = [4.68 4.69 4.695 4.700 4.705 4.710 4.7125 ...
    4.715 4.717 4.718 4.7185 4.7188];
nH_values = [129 257];

[source_grid,resolution_grid] = ndgrid(source_B,nH_values);
source_points = cell(size(source_grid));
parfor k = 1:numel(source_points)
    opts = make_opts(J0eff,Jxx0,resolution_grid(k));
    source_points{k} = invz_solve_point_ordered( ...
        ion,T,source_grid(k),J,opts);
end
assert(all(cellfun(@accepted,source_points),'all'), ...
    'Every path must start from an independently accepted ordered state.');

% Each path is sequential in B, but the two sources and two resolutions are
% mutually independent and may run in parallel.
path_points = cell(size(source_points));
parfor ip = 1:numel(path_points)
    source = source_points{ip};
    points = cell(1,numel(target_B));
    carrier = source;
    for ib = 1:numel(target_B)
        opts = make_opts(J0eff,Jxx0,resolution_grid(ip));
        opts.hmf_profile_state_seed = state_seed(carrier);
        opts.hmf_profile_sweep_direction = 'descending';
        points{ib} = invz_solve_point_ordered( ...
            ion,T,target_B(ib),J,opts);
        if ~accepted(points{ib})
            break;
        end
        carrier = points{ib};
    end
    path_points{ip} = points;
end

records = repmat(empty_record(), ...
    numel(source_B)*numel(nH_values)*numel(target_B),1);
ir = 0;
for in = 1:numel(nH_values)
    for is = 1:numel(source_B)
        points = path_points{is,in};
        for ib = 1:numel(target_B)
            ir = ir+1;
            if isempty(points{ib})
                records(ir) = not_attempted_record( ...
                    target_B(ib),source_B(is),nH_values(in));
            else
                records(ir) = summarize(target_B(ib),source_B(is), ...
                    nH_values(in),points{ib},J,J0eff);
            end
        end
    end
end
summary = struct2table(records);

source_agreement = repmat(empty_agreement(), ...
    numel(target_B)*numel(nH_values),1);
ia = 0;
for in = 1:numel(nH_values)
    for ib = 1:numel(target_B)
        ia = ia+1;
        p1 = path_points{1,in}{ib};
        p2 = path_points{2,in}{ib};
        source_agreement(ia) = compare_points( ...
            target_B(ib),nH_values(in),"source",p1,p2);
    end
end
source_agreement = struct2table(source_agreement);

resolution_agreement = repmat(empty_agreement(), ...
    numel(target_B)*numel(source_B),1);
ia = 0;
for is = 1:numel(source_B)
    for ib = 1:numel(target_B)
        ia = ia+1;
        p1 = path_points{is,1}{ib};
        p2 = path_points{is,2}{ib};
        resolution_agreement(ia) = compare_points( ...
            target_B(ib),source_B(is),"resolution",p1,p2);
    end
end
resolution_agreement = struct2table(resolution_agreement);

accepted_rows = summary.converged;
assert(all(summary.D_uni(accepted_rows) > 0) && ...
    all(summary.supremum_mass(accepted_rows) > 0) && ...
    all(summary.min_mesh_x_mass(accepted_rows) > 0) && ...
    all(summary.min_mesh_medium_mass(accepted_rows) > 0), ...
    'Every accepted traced point must retain positive static masses.');
assert(all(summary.final_resid(accepted_rows) < 1e-8), ...
    'Every accepted traced point must pass the final residual gate.');

result = struct('T',T,'area_factor',1,'source_B',source_B, ...
    'target_B',target_B,'nH_values',nH_values,'summary',summary, ...
    'source_agreement',source_agreement, ...
    'resolution_agreement',resolution_agreement, ...
    'known_refined_pm_mass_zero_T',4.71897990927, ...
    'uses_recovered_targets_as_seeds',true, ...
    'production_changed',false, ...
    'interpretation',['Fine-step field tracing diagnoses numerical reach. ' ...
    'Agreement of these paths is evidence for a connected branch, but does ' ...
    'not authorize production reuse of recovered target states.']);
save(fullfile(here,'wp6_highfield_progressive_field_continuation.mat'), ...
    'result','-v7');

disp(summary(:,{'B','source_B','nH','status','refinement_status', ...
    'converged','hmf','D_uni','Fprime','final_resid','component_nodes'}));
disp(source_agreement);
disp(resolution_agreement);
end

function opts = make_opts(J0eff,Jxx0,nH)
opts = struct('J0eff',J0eff,'Jxx0',Jxx0,'mix_outer',0.35, ...
    'max_outer',2000,'nH',nH, ...
    'hmf_integral_mode','missing_area_approx', ...
    'hmf_missing_area_factor',1, ...
    'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
end

function tf = accepted(point)
tf = point.converged && point.is_ordered && isfinite(point.hmf);
end

function seed = state_seed(point)
seed = struct('Sigma',point.Sigma,'lambda',point.lambda,'K0',point.K(1));
end

function out = summarize(B,source_B,nH,point,J,J0eff)
out = empty_record();
out.B = B;
out.source_B = source_B;
out.nH = nH;
out.status = string(point.hmf_status);
out.converged = accepted(point);
out.hmf = point.hmf;
out.D_uni = point.D_uni;
out.final_resid = point.final_resid;
p = point.hmf_prof;
out.refinement_status = string(p.refinement_failure_status);
if isfield(p,'missing_area_integral') && ...
        isfield(p.missing_area_integral,'node_count')
    m = p.missing_area_integral;
    out.component_nodes = m.node_count;
    out.component_edge = m.component_edge;
    if m.node_count >= 2 && isfinite(m.component_edge) && ...
            isfinite(m.edge_r) && m.edge_r > 0
        ix = m.selected_indices(:).';
        upper_integral = cumtrapz(p.hgrid(ix),p.r(ix));
        required_factor = (J0eff*p.m(ix)-upper_integral)/ ...
            (m.component_edge*m.edge_r);
        out.support_lo = min(required_factor);
        out.support_hi = max(required_factor);
    end
end
if ~(out.converged && isfinite(p.r_star) && ...
        isfinite(p.Gstat_star) && isfinite(point.K(1)))
    return;
end
Gtilde = p.Gstat_star/(1-point.K(1)*p.Gstat_star);
out.Fprime = p.r_star*(1+J0eff*Gtilde);
out.supremum_mass = 1+J0eff*Gtilde;
out.min_mesh_x_mass = min(1+J(:)*Gtilde);
out.min_mesh_medium_mass = min(1+(J(:)-point.K(1))*p.Gstat_star);
end

function out = not_attempted_record(B,source_B,nH)
out = empty_record();
out.B = B;
out.source_B = source_B;
out.nH = nH;
out.status = "not_attempted_after_failure";
end

function out = compare_points(B,label,comparison,p1,p2)
out = empty_agreement();
out.B = B;
out.label = label;
out.comparison = comparison;
if isempty(p1) || isempty(p2)
    return;
end
out.both_converged = accepted(p1) && accepted(p2);
if ~out.both_converged
    return;
end
out.hmf_delta = abs(p1.hmf-p2.hmf);
out.Sigma_delta = max(abs(p1.Sigma-p2.Sigma));
out.D_uni_delta = abs(p1.D_uni-p2.D_uni);
out.component_edge_delta = abs( ...
    p1.hmf_prof.missing_area_integral.component_edge- ...
    p2.hmf_prof.missing_area_integral.component_edge);
out.agrees = out.hmf_delta <= 1e-6 && out.Sigma_delta <= 1e-7 && ...
    out.D_uni_delta <= 1e-6 && out.component_edge_delta <= 1e-9;
end

function out = empty_record()
out = struct('B',NaN,'source_B',NaN,'nH',NaN,'status',"", ...
    'refinement_status',"not_evaluated",'converged',false,'hmf',NaN, ...
    'D_uni',NaN,'Fprime',NaN,'final_resid',NaN,'component_nodes',0, ...
    'component_edge',NaN,'support_lo',NaN,'support_hi',NaN, ...
    'supremum_mass',NaN,'min_mesh_x_mass',NaN, ...
    'min_mesh_medium_mass',NaN);
end

function out = empty_agreement()
out = struct('B',NaN,'label',NaN,'comparison',"", ...
    'both_converged',false,'agrees',false,'hmf_delta',NaN, ...
    'Sigma_delta',NaN,'D_uni_delta',NaN,'component_edge_delta',NaN);
end
