function result = invzp_hmf_node_resolution_census()
%INVZP_HMF_NODE_RESOLUTION_CENSUS Production ordered-profile nH ladder.
% Holds the coupling, outer solver, tolerances, and geometric h-range fixed
% while evaluating nH=33,65,129 at Bx=1,2,3 T. Reports both the final ordered
% mask and the underlying per-node coverage; the independent h=0 predictor is
% retained exactly as production uses it.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
F = load(fixture_path);
ion = invz_ion();
T = F.provenance.T;
Bx_values = [1;2;3];
nH_values = [33;65;129];
nrow = numel(Bx_values)*numel(nH_values);

Bx = nan(nrow,1);
nH = nan(nrow,1);
elapsed_seconds = nan(nrow,1);
point_converged = false(nrow,1);
point_is_ordered = false(nrow,1);
profile_status = strings(nrow,1);
predictor_status = strings(nrow,1);
predictor_converged = false(nrow,1);
n_profile_nodes = nan(nrow,1);
n_converged_nodes = nan(nrow,1);
converged_fraction = nan(nrow,1);
top_contiguous_converged = nan(nrow,1);
min_converged_h = nan(nrow,1);
max_failed_h = nan(nrow,1);
n_static_ok = nan(nrow,1);
n_static_no_root = nan(nrow,1);
n_static_multiple = nan(nrow,1);
n_static_unresolved = nan(nrow,1);
n_outer_failed = nan(nrow,1);
details = cell(nrow,1);

base_opts = struct('mix_outer',0.35,'max_outer',200,'tol_outer',1e-8, ...
    'Ecut',40,'J0eff',F.info.Jcc0,'Jxx0',F.info.Jaa0, ...
    'emt_static',struct('Jsup',F.info.Jcc0,'warn',false, ...
                        'scan_points',4097,'endpoint_margin',1e-10, ...
                        'resid_tol',1e-10));
row = 0;
for ib = 1:numel(Bx_values)
    for in = 1:numel(nH_values)
        row = row+1;
        Bx(row) = Bx_values(ib);
        nH(row) = nH_values(in);
        opts = base_opts;
        opts.nH = nH(row);
        started = tic;
        pt = invz_solve_point_ordered(ion,T,Bx(row),F.J,opts);
        elapsed_seconds(row) = toc(started);
        p = pt.hmf_prof;
        point_converged(row) = pt.converged;
        point_is_ordered(row) = pt.is_ordered;
        profile_status(row) = string(p.status);
        predictor_status(row) = string(p.predictor_static_status);
        predictor_converged(row) = p.predictor_converged;
        n_profile_nodes(row) = numel(p.hgrid);
        n_converged_nodes(row) = nnz(p.node_conv);
        converged_fraction(row) = mean(p.node_conv);
        top_contiguous_converged(row) = top_true_count(p.node_conv);
        if any(p.node_conv)
            min_converged_h(row) = min(p.hgrid(p.node_conv));
        end
        if any(~p.node_conv)
            max_failed_h(row) = max(p.hgrid(~p.node_conv));
        end
        st = string(p.static_status);
        n_static_ok(row) = nnz(st == "ok");
        n_static_no_root(row) = nnz(st == "no_admissible_static_root");
        n_static_multiple(row) = ...
            nnz(st == "multiple_admissible_static_roots");
        n_static_unresolved(row) = nnz(st == "static_search_unresolved");
        n_outer_failed(row) = nnz(st == "outer_iteration_failed");
        details{row} = struct('hgrid',p.hgrid,'node_conv',p.node_conv, ...
            'static_status',p.static_status,'D_uni',p.D_uni, ...
            'Sigma0',p.Sigma0,'K0',p.K0,'profile_status',p.status, ...
            'predictor_static_status',p.predictor_static_status);
        fprintf(['nH census Bx=%.1f nH=%d: point=%d predictor=%s, ' ...
            'nodes %d/%d, min converged h %.6g (%.2fs)\n'], ...
            Bx(row),nH(row),point_converged(row),predictor_status(row), ...
            n_converged_nodes(row),n_profile_nodes(row), ...
            min_converged_h(row),elapsed_seconds(row));
    end
end

tab = table(Bx,nH,elapsed_seconds,point_converged,point_is_ordered, ...
    profile_status,predictor_status,predictor_converged,n_profile_nodes, ...
    n_converged_nodes,converged_fraction,top_contiguous_converged, ...
    min_converged_h,max_failed_h,n_static_ok,n_static_no_root, ...
    n_static_multiple,n_static_unresolved,n_outer_failed);
result = struct('table',tab,'details',{details},'base_opts',base_opts, ...
    'provenance',struct('fixture',fixture_path, ...
                        'coupling_opts',F.provenance.coupling_opts), ...
    'note',['Production ordered-profile path with only nH changed. Final ' ...
            'mask and per-node movement are distinct observables because ' ...
            'the h=0 predictor is independent of nH.']);
save(fullfile(here,'wp2_hmf_node_resolution_census.mat'),'result','-v7');
disp(tab);
end

function n = top_true_count(v)
v = logical(v(:));
idx = find(~v,1,'last');
if isempty(idx)
    n = numel(v);
else
    n = numel(v)-idx;
end
end
