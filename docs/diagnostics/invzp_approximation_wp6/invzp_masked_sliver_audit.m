function result = invzp_masked_sliver_audit()
%INVZP_MASKED_SLIVER_AUDIT Diagnose isolated masks in the production sweep.
% Separate component-coverage and missing-area-bracket failures from the
% ordered/paramagnetic critical boundary.  The three area factors are
% evaluated algebraically on the same certified profile; no failed node is
% inserted into a quadrature and no factor is interpreted statistically.

here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));
F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
J = F.J;
J0eff = F.info.Jcc0;
Jxx0 = F.info.Jaa0;

ion = invz_ion();
T = 0.1;
fields = [0 0.09 0.27 0.36 0.45 0.54 4.59 4.68 4.77];
node_counts = [129 257];
area_factors = [0.75 1 1.5];
B = [fields 0 0.36 0.45 4.68].';
nH = [129*ones(size(fields)) 257 257 257 257].';
records = repmat(empty_record(numel(area_factors)),size(B));

parfor k = 1:numel(B)
    opts = struct('J0eff',J0eff,'Jxx0',Jxx0, ...
        'mix_outer',0.25,'max_outer',1000,'nH',nH(k), ...
        'hmf_integral_mode','missing_area_approx', ...
        'hmf_missing_area_factor',1, ...
        'hmf_approx_branch', ...
        'picard_attracting_contiguous_high_h_component');
    pmopts = rmfield(opts,{'nH','hmf_integral_mode', ...
        'hmf_missing_area_factor','hmf_approx_branch'});
    ordered = [];
    pm = [];
    ordered_error = "";
    pm_error = "";
    try
        ordered = invz_solve_point_ordered(ion,T,B(k),J,opts);
    catch err
        if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
        ordered_error = string(err.identifier);
    end
    try
        pm = invz_solve_point(ion,T,B(k),J,pmopts);
    catch err
        if ~strncmp(err.identifier,'invz:',5), rethrow(err); end
        pm_error = string(err.identifier);
    end
    records(k) = make_record(B(k),nH(k),ordered,pm, ...
        ordered_error,pm_error,area_factors,J0eff);
end

ordered_status = strings(size(B));
ordered_hmf = nan(size(B));
ordered_D_uni = nan(size(B));
Fprime = nan(size(B));
component_nodes = zeros(size(B));
component_edge = nan(size(B));
pm_converged = false(size(B));
pm_crit = nan(size(B));
coarse_brackets = zeros(numel(B),numel(area_factors));
F_edge = nan(numel(B),numel(area_factors));
for k = 1:numel(B)
    r = records(k);
    ordered_status(k) = r.ordered_status;
    ordered_hmf(k) = r.ordered_hmf;
    ordered_D_uni(k) = r.ordered_D_uni;
    Fprime(k) = r.Fprime;
    component_nodes(k) = r.component_nodes;
    component_edge(k) = r.component_edge;
    pm_converged(k) = r.pm_converged;
    pm_crit(k) = r.pm_crit;
    coarse_brackets(k,:) = r.coarse_bracket_count;
    F_edge(k,:) = r.F_edge;
end

summary = table(B(:),nH(:),ordered_status(:),ordered_hmf(:), ...
    ordered_D_uni(:),Fprime(:),component_nodes(:),component_edge(:), ...
    pm_converged(:),pm_crit(:), ...
    'VariableNames',{'B','nH','ordered_status','hmf','D_uni','Fprime', ...
    'component_nodes','component_edge','pm_converged','pm_crit'});

result = struct('T',T,'fields',fields,'node_counts',node_counts, ...
    'area_factors',area_factors,'records',records,'summary',summary, ...
    'coarse_bracket_count',coarse_brackets,'F_edge',F_edge, ...
    'interpretation',['The profile-resolution ladder distinguishes a ' ...
    'missing certified component from an area-dependent bracket loss. The ' ...
    'PM critical mass and ordered F''/D_uni are independent diagnostics of ' ...
    'whether the high-field mask coincides with a continuous boundary.']);
save(fullfile(here,'wp6_masked_sliver_audit.mat'),'result','-v7');
disp(summary);
for k = 1:numel(B)
    fprintf(['B %.2f nH %d brackets [%d %d %d] Fedge ' ...
        '[%.6g %.6g %.6g]\n'],B(k),nH(k), ...
        coarse_brackets(k,1),coarse_brackets(k,2), ...
        coarse_brackets(k,3),F_edge(k,1),F_edge(k,2),F_edge(k,3));
end
end

function out = make_record(B,nH,ordered,pm,ordered_error,pm_error,factors,J0eff)
out = empty_record(numel(factors));
out.B = B;
out.nH = nH;
out.ordered_error = ordered_error;
out.pm_error = pm_error;
if isempty(ordered)
    out.ordered_status = "exception";
    if ~isempty(pm)
        out.pm_converged = pm.converged;
        out.pm_crit = pm.crit;
        out.pm_outer_iters = pm.outer_iters;
    end
    return;
end
p = ordered.hmf_prof;
component_nodes = 0;
component_edge = NaN;
edge_r = NaN;
area_scale = NaN;
factor_support = [NaN NaN];
bracket_count = zeros(1,numel(factors));
F_edge = nan(1,numel(factors));
F_min = nan(1,numel(factors));
F_max = nan(1,numel(factors));

if isfield(p,'missing_area_integral') && ...
        isfield(p.missing_area_integral,'node_count')
    component = p.missing_area_integral;
    component_nodes = component.node_count;
    component_edge = component.component_edge;
    edge_r = component.edge_r;
    if component_nodes >= 2 && isfinite(component_edge) && ...
            isfinite(edge_r) && edge_r > 0
        area_scale = component_edge*edge_r;
        selected = component.selected_mask;
        ix = component.selected_indices(:).';
        upper_integral = cumtrapz(p.hgrid(ix),p.r(ix));
        required_area = J0eff*p.m(ix)-upper_integral;
        factor_support = [min(required_area) max(required_area)]/area_scale;
        for j = 1:numel(factors)
            [h0,meta] = invz_missing_area_integral( ...
                p.hgrid,p.r,selected,factors(j)*area_scale);
            residual = h0-J0eff*p.m;
            valid = meta.selected_mask & isfinite(residual);
            pairs = valid(1:end-1) & valid(2:end) & ...
                residual(1:end-1) < 0 & residual(2:end) >= 0;
            bracket_count(j) = nnz(pairs);
            F_edge(j) = residual(ix(1));
            F_min(j) = min(residual(ix));
            F_max(j) = max(residual(ix));
        end
    end
end

Fprime = NaN;
if isfinite(p.r_star) && isfinite(p.Gstat_star) && ...
        isfield(ordered,'K') && isfinite(ordered.K(1))
    Gtilde_star = p.Gstat_star/(1-ordered.K(1)*p.Gstat_star);
    Fprime = p.r_star*(1+J0eff*Gtilde_star);
end
out = struct('B',B,'nH',nH, ...
    'ordered_status',string(ordered.hmf_status), ...
    'ordered_error',ordered_error,'pm_error',pm_error, ...
    'ordered_converged',ordered.converged, ...
    'ordered_hmf',ordered.hmf,'ordered_D_uni',ordered.D_uni, ...
    'Fprime',Fprime,'component_nodes',component_nodes, ...
    'component_edge',component_edge,'edge_r',edge_r, ...
    'area_scale',area_scale,'factor_support',factor_support, ...
    'coarse_bracket_count',bracket_count,'F_edge',F_edge, ...
    'F_min',F_min,'F_max',F_max, ...
    'predictor_converged',p.predictor_converged, ...
    'predictor_slope',p.slope0,'pm_converged',false, ...
    'pm_crit',NaN,'pm_outer_iters',NaN);
if ~isempty(pm)
    out.pm_converged = pm.converged;
    out.pm_crit = pm.crit;
    out.pm_outer_iters = pm.outer_iters;
end
end

function out = empty_record(nfactor)
out = struct('B',NaN,'nH',NaN,'ordered_status',"not_evaluated", ...
    'ordered_error',"",'pm_error',"", ...
    'ordered_converged',false,'ordered_hmf',NaN,'ordered_D_uni',NaN, ...
    'Fprime',NaN,'component_nodes',0,'component_edge',NaN, ...
    'edge_r',NaN,'area_scale',NaN,'factor_support',[NaN NaN], ...
    'coarse_bracket_count',zeros(1,nfactor), ...
    'F_edge',nan(1,nfactor),'F_min',nan(1,nfactor), ...
    'F_max',nan(1,nfactor),'predictor_converged',false, ...
    'predictor_slope',NaN,'pm_converged',false,'pm_crit',NaN, ...
    'pm_outer_iters',NaN);
end
