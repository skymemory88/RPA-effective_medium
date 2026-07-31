function result = invzp_highfield_one_sided_seed_audit()
%INVZP_HIGHFIELD_ONE_SIDED_SEED_AUDIT Ordered-side continuation near Bc.
% Diagnostic only. Two independently accepted lower-field states seed the
% same target fields. This tests numerical reach and agreement, not whether
% a one-sided retry is an authorized branch/area selector.

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
target_B = [driver_fields(target_index) 4.68 4.70];
base = make_opts(J0eff,Jxx0);

sources = cell(size(source_B));
parfor k = 1:numel(source_B)
    sources{k} = invz_solve_point_ordered(ion,T,source_B(k),J,base);
end
assert(all(cellfun(@(p)p.converged && p.is_ordered,sources)), ...
    'Every same-side continuation source must be independently accepted.');

[target_grid,source_grid] = ndgrid(target_B,source_B);
points = cell(size(target_grid));
seed_state = cell(size(target_grid));
for k = 1:numel(target_grid)
    [~,source_index] = min(abs(source_B-source_grid(k)));
    source = sources{source_index};
    seed_state{k} = struct('Sigma',source.Sigma, ...
        'lambda',source.lambda,'K0',source.K(1));
end
parfor k = 1:numel(target_grid)
    opts = base;
    opts.hmf_profile_state_seed = seed_state{k};
    opts.hmf_profile_sweep_direction = 'descending';
    points{k} = invz_solve_point_ordered( ...
        ion,T,target_grid(k),J,opts);
end

records = repmat(empty_record(),numel(target_grid),1);
for k = 1:numel(target_grid)
    records(k) = summarize(target_grid(k),source_grid(k), ...
        points{k},J0eff);
end
summary = struct2table(records);

agreement = repmat(struct('B',NaN,'both_converged',false, ...
    'hmf_delta',NaN,'Sigma_delta',NaN,'D_uni_delta',NaN, ...
    'component_edge_delta',NaN,'support_lo_delta',NaN, ...
    'support_hi_delta',NaN),numel(target_B),1);
for k = 1:numel(target_B)
    a = points{k,1};
    b = points{k,2};
    ra = records(k);
    rb = records(k+numel(target_B));
    both_converged = a.converged && b.converged;
    hmf_delta = NaN;
    Sigma_delta = NaN;
    D_uni_delta = NaN;
    if both_converged && isfield(a,'Sigma') && isfield(b,'Sigma')
        hmf_delta = abs(a.hmf-b.hmf);
        Sigma_delta = max(abs(a.Sigma-b.Sigma));
        D_uni_delta = abs(a.D_uni-b.D_uni);
    end
    agreement(k) = struct('B',target_B(k), ...
        'both_converged',both_converged, ...
        'hmf_delta',hmf_delta,'Sigma_delta',Sigma_delta, ...
        'D_uni_delta',D_uni_delta, ...
        'component_edge_delta',abs(ra.component_edge-rb.component_edge), ...
        'support_lo_delta',abs(ra.support_lo-rb.support_lo), ...
        'support_hi_delta',abs(ra.support_hi-rb.support_hi));
end
agreement_table = struct2table(agreement);

pm_mass = nan(size(target_B));
pm_converged = false(size(target_B));
parfor k = 1:numel(target_B)
    popts = struct('J0eff',J0eff,'Jxx0',Jxx0, ...
        'mix_outer',0.25,'max_outer',1000);
    pm = invz_solve_point(ion,T,target_B(k),J,popts);
    pm_mass(k) = pm.crit;
    pm_converged(k) = pm.converged;
end
pm_table = table(target_B(:),pm_converged(:),pm_mass(:), ...
    'VariableNames',{'B_T','converged','pm_mass'});

result = struct('T',T,'nH',129,'area_factor',1, ...
    'source_B',source_B,'target_B',target_B,'summary',summary, ...
    'agreement',agreement_table,'pm_table',pm_table, ...
    'known_refined_pm_mass_zero_T',4.71897990927, ...
    'interpretation',['Agreement of two lower-field seeds establishes a ' ...
    'reproducible ordered-side numerical continuation. It does not by ' ...
    'itself select that branch or its missing area, particularly when the ' ...
    'root approaches a uniform-stability boundary.']);
save(fullfile(here,'wp6_highfield_one_sided_seed_audit.mat'), ...
    'result','-v7');
disp(summary);
disp(agreement_table);
disp(pm_table);
end

function opts = make_opts(J0eff,Jxx0)
opts = struct('J0eff',J0eff,'Jxx0',Jxx0,'mix_outer',0.25, ...
    'max_outer',1000,'nH',129, ...
    'hmf_integral_mode','missing_area_approx', ...
    'hmf_missing_area_factor',1, ...
    'hmf_approx_branch', ...
    'picard_attracting_contiguous_high_h_component');
end

function out = summarize(B,seed_B,point,J0eff)
out = empty_record();
out.B = B;
out.seed_B = seed_B;
out.status = string(point.hmf_status);
out.converged = point.converged;
out.hmf = point.hmf;
out.D_uni = point.D_uni;
out.final_resid = point.final_resid;
p = point.hmf_prof;
if isfinite(p.r_star) && isfinite(p.Gstat_star) && ...
        isfield(point,'K') && isfinite(point.K(1))
    Gtilde_star = p.Gstat_star/(1-point.K(1)*p.Gstat_star);
    out.Fprime_identity = p.r_star*(1+J0eff*Gtilde_star);
end
m = p.missing_area_integral;
out.component_nodes = m.node_count;
out.component_edge = m.component_edge;
out.root_bracket_count = p.root_bracket_count;
if all(isfinite(p.root_bracket_indices))
    ia = p.root_bracket_indices(1);
    ib = p.root_bracket_indices(2);
    out.bracket_F_lo = p.F(ia);
    out.bracket_F_hi = p.F(ib);
    out.bracket_secant = (p.F(ib)-p.F(ia))/ ...
        (p.hgrid(ib)-p.hgrid(ia));
end
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
out = struct('B',NaN,'seed_B',NaN,'status',"",'converged',false, ...
    'hmf',NaN,'D_uni',NaN,'final_resid',NaN, ...
    'Fprime_identity',NaN,'component_nodes',0,'component_edge',NaN, ...
    'root_bracket_count',0,'bracket_F_lo',NaN,'bracket_F_hi',NaN, ...
    'bracket_secant',NaN,'support_lo',NaN,'support_hi',NaN);
end
