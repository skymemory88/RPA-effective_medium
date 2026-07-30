function [K0, Gstat, out] = invz_emt_static_ordered( ...
    tl, lam, Sigma0, Jnu_flat, ~, beta, J0eff, G0inel0, G0el0, opts)
%INVZ_EMT_STATIC_ORDERED Bounded physical-x ordered static EMT solve.
% The solve coordinate is
%   x = Gstat/(1-K0*Gstat),  -1/Jsup < x < 0,
% where Jsup includes the lattice spectrum and the separately weighted
% uniform/directional mode. On this interval
%   Phi(x) = mean_q x/(1+J(q)x),
%   K0(x)  = 1/Phi(x)-1/x,
% and roots of
%   Rhat(x) = Phi(x)-invz_gstat_ordered(...,K0(x),...)
% are enumerated over the full configured scan interval.
%
% No sign-changing bracket is accepted unless refinement reaches resid_tol at
% a finite point. Pole/discontinuity brackets are reported separately. If more
% than one admissible root survives, no arbitrary branch is exported:
% out.status='multiple_admissible_static_roots' and out.converged=false.
% K0_seed is intentionally ignored, making the inner result seed-independent.
%
% opts:
%   Jsup             verified spectral/uniform supremum; default
%                    max(J0eff,max(Jnu_flat))
%   resid_tol        both Rhat and EMT closure tolerance (default 1e-10)
%   scan_points      endpoint-clustered scan size (default 4097)
%   endpoint_margin  normalized distance from x=0 and x=-1/Jsup (1e-10)
%   root_x_tol       normalized root-refinement tolerance (1e-13)
%   root_merge_tol   normalized candidate deduplication tolerance (1e-10)
%   root_maxit       refinement budget per candidate (120)
%   lattice_block    coupling-map work block (256)
%   mass_tol         strict positive-denominator floor (default 0)
%   lambda_affine    optional struct with vectors base/slope defining
%                    lambda_p(K0)=base_p+slope_p*K0 for a deterministic
%                    outer-map evaluation; absent preserves fixed input lam
%   warn             emit status warning (default true)
%
% A nonnegative finite xi with a positive xi denominator is the hard xi-domain
% gate. No universal xi upper bound is imposed: the production expression does
% not supply one, and inventing one would change the physical closure.
if nargin < 10, opts = struct(); end

if ~isvector(Jnu_flat) || isempty(Jnu_flat)
    error('invz:emtStaticDomain', 'Jnu_flat must be a nonempty static coupling vector.');
end
Jf = Jnu_flat(:);
if ~isreal(Jf) || any(~isfinite(Jf))
    error('invz:emtStaticDomain', 'Jnu_flat must contain finite real couplings.');
end
if ~(isnumeric(Sigma0) && isreal(Sigma0) && isscalar(Sigma0) && isfinite(Sigma0))
    error('invz:emtStaticDomain', 'Sigma0 must be a finite real scalar.');
end
if ~(isnumeric(beta) && isreal(beta) && isscalar(beta) && isfinite(beta) && beta > 0)
    error('invz:emtStaticDomain', 'beta must be a finite positive scalar.');
end
lambda_affine = getf(opts,'lambda_affine',[]);
affine_lambda = ~isempty(lambda_affine);
if affine_lambda
    if ~(isstruct(lambda_affine) && isscalar(lambda_affine) && ...
            all(isfield(lambda_affine,{'base','slope'})))
        error('invz:emtStaticDomain', ...
            'opts.lambda_affine must be a scalar struct with base and slope vectors.');
    end
    lambda_base = lambda_affine.base(:);
    lambda_slope = lambda_affine.slope(:);
    if numel(lambda_base) < 2 || numel(lambda_slope) ~= numel(lambda_base) || ...
            any(~isfinite(lambda_base)) || any(~isfinite(lambda_slope)) || ...
            ~isreal(lambda_base) || ~isreal(lambda_slope)
        error('invz:emtStaticDomain', ...
            'lambda_affine base/slope must be equal-length finite real vectors with >=2 entries.');
    end
else
    if numel(lam) < 2 || any(~isfinite(lam(1:2))) || ~isreal(lam(1:2))
        error('invz:emtStaticDomain', 'lam must provide two finite real static components.');
    end
    lambda_base = [];
    lambda_slope = [];
end
if ~(isfinite(G0inel0) && isreal(G0inel0) && isscalar(G0inel0) && ...
     isfinite(G0el0) && isreal(G0el0) && isscalar(G0el0))
    error('invz:emtStaticDomain', 'G0inel0 and G0el0 must be finite real scalars.');
end

Jsup = getf(opts, 'Jsup', max([J0eff; Jf]));
if ~(isnumeric(Jsup) && isreal(Jsup) && isscalar(Jsup) && isfinite(Jsup) && Jsup > 0)
    error('invz:emtStaticDomain', 'opts.Jsup must be a finite positive scalar.');
end
sup_tol = 64*eps(max(1, abs(Jsup)));
if Jsup + sup_tol < J0eff || Jsup + sup_tol < max(Jf)
    error('invz:emtStaticDomain', ...
        'opts.Jsup=%.17g is below J0eff or the sampled lattice maximum %.17g.', ...
        Jsup, max(Jf));
end

rtol = getf(opts, 'resid_tol', 1e-10);
nscan = getf(opts, 'scan_points', 4097);
emargin = getf(opts, 'endpoint_margin', 1e-10);
xtol = getf(opts, 'root_x_tol', 1e-13);
merge_tol = getf(opts, 'root_merge_tol', 1e-10);
root_maxit = getf(opts, 'root_maxit', 120);
lat_block = getf(opts, 'lattice_block', 256);
mass_tol = getf(opts, 'mass_tol', 0);
warn = getf(opts, 'warn', true);
if ~(isscalar(nscan) && nscan == round(nscan) && nscan >= 17)
    error('invz:emtStaticDomain', 'scan_points must be an integer >=17.');
end
if ~(isscalar(emargin) && emargin > 0 && emargin < 0.5)
    error('invz:emtStaticDomain', 'endpoint_margin must lie strictly between 0 and 0.5.');
end

map = local_lattice_map(Jf, Jsup, nscan, emargin, lat_block);
lam_scan = lambda_at_K(map.K0);
[Gs_scan, sd] = local_gstat_array(tl, lam_scan, map.K0, Sigma0, beta, G0inel0, G0el0);
Rscan = map.Phi - Gs_scan;
point_ok = isfinite(Rscan) & isfinite(sd.xi) & ...
           isfinite(sd.xi_denom) & sd.xi_denom ~= 0 & ...
           isfinite(sd.local_denom) & sd.local_denom ~= 0;
% A sign change of either rational denominator is a known discontinuity edge.
edge_ok = point_ok(1:end-1) & point_ok(2:end) & ...
          sign(sd.local_denom(1:end-1)) == sign(sd.local_denom(2:end)) & ...
          sign(sd.xi_denom(1:end-1)) == sign(sd.xi_denom(2:end));

bropts = struct('resid_tol', rtol, 'x_tol', xtol, 'maxit', root_maxit, ...
    'merge_tol', merge_tol, 'fgrid', Rscan, 'edge_ok', edge_ok);
[sroots, search] = invz_bounded_roots(@residual_at_s, map.s, bropts);

template = local_record_template();
records = repmat(template, numel(sroots), 1);
for k = 1:numel(sroots)
    records(k) = candidate_record(sroots(k));
end
if isempty(records)
    root_table = struct2table(repmat(template, 0, 1));
    accepted = false(0,1);
else
    root_table = struct2table(records);
    accepted = [records.accepted].';
end
naccepted = nnz(accepted);

K0 = NaN;
Gstat = NaN;
go = local_empty_go();
status = 'no_admissible_static_root';
selected_index = NaN;
if search.unresolved_minima > 0
    status = 'static_search_unresolved';
elseif naccepted == 1
    selected_index = find(accepted, 1);
    K0 = records(selected_index).K0;
    Gstat = records(selected_index).Gstat;
    selected_lam = lambda_at_K(K0);
    [~, go] = invz_gstat_ordered(tl, selected_lam, K0, Sigma0, beta, ...
        G0inel0, G0el0, struct('stable_form', true));
    status = 'ok';
elseif naccepted > 1
    status = 'multiple_admissible_static_roots';
end

out = go;
out.status = status;
out.converged = strcmp(status, 'ok');
out.resid = NaN;
out.root_resid = NaN;
out.D_uni = NaN;
if out.converged
    out.resid = records(selected_index).closure_resid;
    out.root_resid = records(selected_index).root_resid;
    out.D_uni = records(selected_index).D_uni;
end
out.iters = search.refine_evals;
out.Jsup = Jsup;
out.lambda_mode = 'fixed';
if affine_lambda, out.lambda_mode = 'affine_K0'; end
out.selected_lambda = [];
if out.converged, out.selected_lambda = lambda_at_K(K0); end
out.seed_ignored = true;
out.selected_index = selected_index;
out.n_roots = numel(records);
out.n_admissible_roots = naccepted;
out.root_table = root_table;
out.search = search;
out.search.sgrid = map.s;
out.search.xgrid = map.x;
out.search.point_ok = point_ok;
out.search.edge_ok = edge_ok;
out.search.discontinuity_count = size(search.discontinuity_brackets, 1);
out.lattice = struct('mesh_max', max(Jf), 'J0eff', J0eff, 'Jsup', Jsup, ...
    'scan_points', nscan, 'endpoint_margin', emargin);

if warn && ~out.converged
    switch status
        case 'no_admissible_static_root'
            warning('invz:emtStaticNoRoot', ...
                'no admissible ordered static root on -1/Jsup < x < 0.');
        case 'multiple_admissible_static_roots'
            warning('invz:emtStaticMultipleRoots', ...
                '%d admissible ordered static roots found; no branch exported.', naccepted);
        otherwise
            warning('invz:emtStaticSearch', ...
                'ordered static root search was numerically unresolved.');
    end
end

    function R = residual_at_s(s)
        [~, ~, R, ev] = eval_at_s(s);
        if ~ev.continuous, R = NaN; end
    end

    function [x, kval, R, ev] = eval_at_s(s)
        x = -s/Jsup;
        dx = 1 + Jf*x;
        if ~(isfinite(x) && x < 0 && 1 + Jsup*x > 0 && all(dx > 0))
            kval = NaN; R = NaN;
            ev = struct('Gs',NaN,'go',local_empty_go(),'Phi',NaN, ...
                'xi_physical',false,'continuous',false,'local_denom',NaN);
            return;
        end
        invdx = 1./dx;
        M = mean(invdx);
        Phi = x*M;
        kval = mean(Jf.*invdx)/M; % stable form of 1/Phi-1/x
        lam0 = lambda_at_K(kval);
        [Gs, go0] = invz_gstat_ordered(tl, lam0, kval, Sigma0, beta, ...
            G0inel0, G0el0, struct('stable_form', true));
        R = Phi - Gs;
        xi_physical = isfinite(go0.xi) && go0.xi >= 0 && ...
                      isfinite(go0.xi_denom) && go0.xi_denom > 0;
        continuous = isfinite(R) && isfinite(go0.xi) && ...
                     isfinite(go0.xi_denom) && go0.xi_denom ~= 0 && ...
                     isfinite(go0.gstat_local_denom) && go0.gstat_local_denom ~= 0;
        ev = struct('Gs',Gs,'go',go0,'Phi',Phi, ...
            'xi_physical',xi_physical,'continuous',continuous, ...
            'local_denom',go0.gstat_local_denom);
    end

    function rec = candidate_record(s)
        rec = template;
        [x, kval, R, ev] = eval_at_s(s);
        go0 = ev.go;
        Gs = ev.Gs;
        Dqx = 1 + Jf*x;
        Dq = 1 + (Jf-kval)*Gs;
        Duni = 1 + (J0eff-kval)*Gs;
        Gq = Gs./Dq;
        closure_resid = abs(mean(Gq)-Gs);
        xinel = G0inel0/(1+Sigma0);
        xprod = go0.Gtil0;

        reasons = strings(0,1);
        if ~(isfinite(x) && x < 0), reasons(end+1) = "x_not_negative"; end
        if ~(isfinite(1+Jsup*x) && 1+Jsup*x > mass_tol)
            reasons(end+1) = "nonpositive_supremum_mass";
        end
        if ~all(isfinite([kval,Gs,go0.U,go0.V,go0.xi]))
            reasons(end+1) = "nonfinite_static_value";
        end
        if ~ev.xi_physical, reasons(end+1) = "nonphysical_xi"; end
        if ~(isfinite(R) && abs(R) <= rtol), reasons(end+1) = "root_residual"; end
        if ~(isfinite(closure_resid) && closure_resid <= rtol)
            reasons(end+1) = "closure_residual";
        end
        if ~(isfinite(Gs) && Gs < 0), reasons(end+1) = "nonnegative_Gstat"; end
        if ~(isfinite(Duni) && Duni > mass_tol)
            reasons(end+1) = "nonpositive_uniform_mass";
        end
        if ~(all(isfinite(Dqx)) && min(Dqx) > mass_tol)
            reasons(end+1) = "nonpositive_mesh_x_mass";
        end
        if ~(all(isfinite(Dq)) && min(Dq) > mass_tol)
            reasons(end+1) = "nonpositive_mesh_medium_mass";
        end

        rec.G0_inel = G0inel0;
        rec.G0_el = G0el0;
        rec.xi = go0.xi;
        rec.xi_numer = go0.xi_numer;
        rec.xi_denom = go0.xi_denom;
        rec.U = go0.U;
        rec.V = go0.V;
        rec.Gstat = Gs;
        rec.K0 = kval;
        rec.Sigma0 = Sigma0;
        lam0 = lambda_at_K(kval);
        rec.lambda1 = lam0(1);
        rec.lambda2 = lam0(2);
        rec.x = x;
        rec.root_resid = abs(R);
        rec.closure_resid = closure_resid;
        rec.supremum_mass = 1 + Jsup*x;
        rec.D_uni = Duni;
        rec.min_mesh_x_abs = min(abs(Dqx));
        rec.min_mesh_x_signed = min(Dqx);
        rec.min_mesh_medium_signed = min(Dq);
        rec.x_inel = xinel;
        rec.x_inel_discrepancy = x-xinel;
        rec.x_production = xprod;
        rec.x_production_discrepancy = x-xprod;
        rec.J0eff_abs_G0inel = J0eff*abs(G0inel0);
        rec.Jsup_abs_x_inel = Jsup*abs(xinel);
        rec.gstat_local_denom = go0.gstat_local_denom;
        rec.accepted = isempty(reasons);
        if isempty(reasons), rec.reject_reason = "ok";
        else,                rec.reject_reason = strjoin(reasons, "|"); end
    end

    function lamk = lambda_at_K(kval)
        if ~affine_lambda
            lamk = lam(1:2);
        elseif isscalar(kval)
            lamk = lambda_base + lambda_slope*kval;
        else
            lamk = lambda_base(1:2) + lambda_slope(1:2).*kval(:).';
        end
    end
end

function map = local_lattice_map(Jf, Jsup, nscan, margin, block)
% Cache the outer-state-independent Phi/K0 scan for repeated Sigma iterations.
persistent oldJ oldJsup oldN oldMargin oldMap
if ~isempty(oldMap) && isequal(oldJ, Jf) && isequal(oldJsup, Jsup) && ...
        isequal(oldN, nscan) && isequal(oldMargin, margin)
    map = oldMap;
    return;
end
theta = linspace(0, pi, nscan);
s = margin + (1-2*margin)*(1-cos(theta))/2;
x = -s/Jsup;
Phi = nan(size(x));
K0 = nan(size(x));
for i0 = 1:block:nscan
    idx = i0:min(i0+block-1, nscan);
    den = 1 + Jf*x(idx);
    invden = 1./den;
    M = mean(invden, 1);
    Phi(idx) = x(idx).*M;
    K0(idx) = mean(Jf.*invden, 1)./M;
end
map = struct('s',s(:),'x',x(:),'Phi',Phi(:),'K0',K0(:));
oldJ = Jf; oldJsup = Jsup; oldN = nscan; oldMargin = margin; oldMap = map;
end

function [Gs, d] = local_gstat_array(tl, lam, K0, Sigma0, beta, G0i, G0e)
% The production closure itself is vector-capable; do not maintain a replica.
[Gs, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0i, G0e);
d = struct('xi',go.xi,'xi_denom',go.xi_denom, ...
    'local_denom',go.gstat_local_denom);
end

function rec = local_record_template()
rec = struct('G0_inel',NaN,'G0_el',NaN,'xi',NaN,'xi_numer',NaN, ...
    'xi_denom',NaN,'U',NaN,'V',NaN,'Gstat',NaN,'K0',NaN,'Sigma0',NaN, ...
    'lambda1',NaN,'lambda2',NaN, ...
    'x',NaN,'root_resid',NaN,'closure_resid',NaN,'supremum_mass',NaN, ...
    'D_uni',NaN,'min_mesh_x_abs',NaN,'min_mesh_x_signed',NaN, ...
    'min_mesh_medium_signed',NaN,'x_inel',NaN,'x_inel_discrepancy',NaN, ...
    'x_production',NaN,'x_production_discrepancy',NaN, ...
    'J0eff_abs_G0inel',NaN,'Jsup_abs_x_inel',NaN, ...
    'gstat_local_denom',NaN,'accepted',false,'reject_reason',"");
end

function go = local_empty_go()
go = struct('xi',NaN,'h0',NaN,'G0bare',NaN,'Gtil0',NaN,'r',NaN, ...
    'gstat_local_denom',NaN,'G0inel0',NaN,'G0el0',NaN, ...
    'xi_numer',NaN,'xi_denom',NaN,'U',NaN,'V',NaN);
end
