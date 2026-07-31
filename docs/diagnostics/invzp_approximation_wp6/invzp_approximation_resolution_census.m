function result = invzp_approximation_resolution_census()
%INVZP_APPROXIMATION_RESOLUTION_CENSUS Spot-check the certified-edge grid error.
% Resolution is a separate numerical systematic from the missing-area shape
% ensemble.  Compare the central member at representative low/intermediate/
% high fields using the opt-in driver resolution and one doubled grid.
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
fields = [0.2 0.4 0.5 1.0 2.7];
expected_covered = [true true false true true];
node_counts = [129 257];
[B,nH] = ndgrid(fields,node_counts);
hmf = nan(size(B));
component_edge = nan(size(B));
missing_area = nan(size(B));
status = strings(size(B));
bracket_count = nan(size(B));
parfor k = 1:numel(B)
    opts = struct('J0eff',J0eff,'Jxx0',Jxx0, ...
        'mix_outer',0.25,'max_outer',1000,'nH',nH(k), ...
        'hmf_integral_mode','missing_area_approx', ...
        'hmf_missing_area_factor',1);
    p = invz_solve_point_ordered(ion,T,B(k),J,opts);
    hmf(k) = p.hmf;
    status(k) = string(p.hmf_status);
    component_edge(k) = ...
        p.hmf_prof.missing_area_integral.component_edge;
    missing_area(k) = p.hmf_prof.missing_area;
    bracket_count(k) = p.hmf_prof.root_bracket_count;
end

hmf_change = abs(hmf(:,2)-hmf(:,1));
edge_change = abs(component_edge(:,2)-component_edge(:,1));
assert(all(status(expected_covered,:) == "ok_missing_area_approx",'all'));
assert(all(status(~expected_covered,:) == ...
    "missing_area_no_certified_component",'all'));
assert(all(bracket_count(expected_covered,:) == 1,'all'));
assert(all(bracket_count(~expected_covered,:) == 0,'all'));
assert(all(hmf_change(expected_covered) < 1e-3));
summary = table(fields.',hmf(:,1),hmf(:,2),hmf_change, ...
    component_edge(:,1),component_edge(:,2),edge_change, ...
    'VariableNames',{'B','hmf_129','hmf_257','hmf_change', ...
    'edge_129','edge_257','edge_change'});
result = struct('T',T,'fields',fields,'node_counts',node_counts, ...
    'expected_covered',expected_covered, ...
    'hmf',hmf,'component_edge',component_edge, ...
    'missing_area',missing_area,'status',status, ...
    'bracket_count',bracket_count,'summary',summary, ...
    'interpretation',['Grid sensitivity is not part of the missing-area ' ...
    'shape band. A changed terminal accepted edge can dominate the residual ' ...
    'quadrature error and must remain separately reported.']);
save(fullfile(here,'wp6_approximation_resolution_census.mat'), ...
    'result','-v7');
disp(summary);
end
