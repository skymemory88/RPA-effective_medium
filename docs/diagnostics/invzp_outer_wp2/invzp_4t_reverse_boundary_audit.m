function result = invzp_4t_reverse_boundary_audit()
%INVZP_4T_REVERSE_BOUNDARY_AUDIT Static scan/margin audit at the stop node.
% Re-evaluates the failed coarse continuation transfer without changing its
% Sigma seed or applying an outer iteration.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
reverse_path = fullfile(here,'wp2_4t_reverse_continuation.mat');
F = load(fixture_path);
R = load(reverse_path);
failed_node = R.result.summary.failed_node;
if ~(isscalar(failed_node) && isfinite(failed_node) && failed_node >= 1 && ...
        failed_node < height(R.result.table))
    error('invz:reverseBoundaryAudit', ...
        'reverse-continuation artifact has no auditable failed transfer node.');
end

ion = invz_ion();
T = R.result.config.T;
Bx = R.result.config.Bx;
h = R.result.table.h(failed_node);
[wn,wts,beta] = invz_matsubara(T,R.result.config.Ecut);
Sigma_seed = R.result.Sigma_solution(:,failed_node+1);
if any(~isfinite(Sigma_seed))
    error('invz:reverseBoundaryAudit', ...
        'the immediately higher node does not contain a certified Sigma seed.');
end
si = invz_single_ion(ion,T,[Bx 0 0], ...
    struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h));
tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
ctx = make_context(si,tl,T,wn,wts,beta,F.J,F.info);

scan_points = [2049;4097;8193];
endpoint_margin = [1e-8;1e-10;1e-12];
n = numel(scan_points);
status = strings(n,1);
defined = false(n,1);
n_roots = nan(n,1);
n_admissible_roots = nan(n,1);
unresolved_minima = nan(n,1);
discontinuity_count = nan(n,1);
max_grid_gap = nan(n,1);
min_sample_abs_residual = nan(n,1);
details = cell(n,1);
for k = 1:n
    opts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false, ...
        'scan_points',scan_points(k),'endpoint_margin',endpoint_margin(k), ...
        'resid_tol',1e-10));
    q = invz_ordered_outer_map(Sigma_seed,ctx,opts);
    details{k} = q;
    status(k) = string(q.status);
    defined(k) = q.defined;
    n_roots(k) = q.static.n_roots;
    n_admissible_roots(k) = q.static.n_admissible_roots;
    unresolved_minima(k) = q.static.search.unresolved_minima;
    discontinuity_count(k) = q.static.search.discontinuity_count;
    max_grid_gap(k) = q.static.search.max_grid_gap;
    finite_residual = abs(q.static.search.fgrid( ...
        isfinite(q.static.search.fgrid)));
    if ~isempty(finite_residual)
        min_sample_abs_residual(k) = min(finite_residual);
    end
end
tab = table(scan_points,endpoint_margin,status,defined,n_roots, ...
    n_admissible_roots,unresolved_minima,discontinuity_count,max_grid_gap, ...
    min_sample_abs_residual);
classification_stable = isscalar(unique(status)) && ...
    all(~defined) && all(n_roots == 0) && all(n_admissible_roots == 0) && ...
    all(unresolved_minima == 0) && all(discontinuity_count == 0);
result = struct('failed_node',failed_node,'h',h, ...
    'seed_node',failed_node+1,'table',tab,'classification_stable', ...
    classification_stable,'details',{details}, ...
    'provenance',struct('fixture',fixture_path, ...
                        'reverse_continuation',reverse_path), ...
    'note',['Frozen-seed finite-resolution audit only. Stable no-root ' ...
            'classification does not prove that no coupled branch exists ' ...
            'at this h or between the two coarse profile nodes.']);
save(fullfile(here,'wp2_4t_reverse_boundary_audit.mat'),'result','-v7');
disp(tab);
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
