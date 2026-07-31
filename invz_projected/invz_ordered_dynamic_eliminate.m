function out = invz_ordered_dynamic_eliminate(lambda, K0, ctx, opts)
%INVZ_ORDERED_DYNAMIC_ELIMINATE Eliminate ordered n>0 K_n at fixed global data.
% Diagnostic-only exact reduction of the dynamic ordered equations. At fixed
% (lambda_1,lambda_2,lambda_3,K0), Jensen's ordered self-energy is affine:
%
%   Sigma_n = c_n - d_n K_n,
%   d_n = (M2/n01^2)(1-n01^2)g_n,                  n > 0.
%
% Combining this with equations (16)-(17) gives one bounded scalar equation
% per Matsubara frequency:
%
%   K_n = sum_q J_q/L_q / sum_q 1/L_q,
%   L_q = 1 + c_n - d_n K_n + J_q G0_n.
%
% A physical root has K_n in [min J_q,max J_q] and L_q>0 for every q. The
% routine proves uniqueness from a derivative bound when possible and otherwise
% performs finite-resolution root enumeration. It never selects one root from a
% multiple-root result.
%
% Required CTX fields: G0,g,tl,wts,beta,Jnu_flat.
% OPTS:
%   mass_tol          strict L_q/D_q floor (default 0)
%   resid_tol         scalar equation tolerance (default 1e-12)
%   root_x_tol        K refinement tolerance in meV (default 1e-14)
%   scan_points       fallback enumeration grid (default 129)
%   endpoint_margin   fraction of K span kept from a pole (default 1e-10)
%   uniqueness_slack  require derivative bound < 1-slack (default 1e-10)
if nargin < 4, opts = struct(); end
required = {'G0','g','tl','wts','beta','Jnu_flat'};
if ~(isstruct(ctx) && isscalar(ctx) && all(isfield(ctx,required)))
    error('invz:reducedDynamicContext', ...
        'ctx must be a scalar struct with fields: %s.',strjoin(required,', '));
end
lambda = lambda(:);
if numel(lambda) ~= 3 || any(~isfinite(lambda)) || ~isreal(lambda)
    error('invz:reducedDynamicLambda', 'lambda must contain three finite real values.');
end
if ~(isnumeric(K0) && isreal(K0) && isscalar(K0) && isfinite(K0))
    error('invz:reducedDynamicK0', 'K0 must be a finite real scalar.');
end

G0 = ctx.G0(:);
g = ctx.g(:);
wts = ctx.wts(:);
nw = numel(G0);
if numel(g) ~= nw || numel(wts) ~= nw
    error('invz:reducedDynamicShape', 'ctx.G0, ctx.g, and ctx.wts must have equal lengths.');
end
if any(~isfinite(G0)) || any(~isfinite(g)) || ~isreal(G0) || ~isreal(g)
    error('invz:reducedDynamicShape', 'ctx.G0 and ctx.g must be finite real vectors.');
end

retarded = ~isvector(ctx.Jnu_flat);
if retarded
    if size(ctx.Jnu_flat,2) ~= nw
        error('invz:reducedDynamicJ', ...
            'Matrix Jnu_flat must have one column per Matsubara frequency.');
    end
    if any(~isfinite(ctx.Jnu_flat),'all') || ~isreal(ctx.Jnu_flat)
        error('invz:reducedDynamicJ', 'Jnu_flat must contain finite real values.');
    end
else
    Jshared = ctx.Jnu_flat(:);
    if isempty(Jshared) || any(~isfinite(Jshared)) || ~isreal(Jshared)
        error('invz:reducedDynamicJ', 'Jnu_flat must contain finite real values.');
    end
end

mass_tol = getf(opts,'mass_tol',0);
rtol = getf(opts,'resid_tol',1e-12);
xtol = getf(opts,'root_x_tol',1e-14);
nscan = getf(opts,'scan_points',129);
emargin = getf(opts,'endpoint_margin',1e-10);
uslack = getf(opts,'uniqueness_slack',1e-10);
if ~(isscalar(nscan) && nscan == round(nscan) && nscan >= 17)
    error('invz:reducedDynamicOptions', 'scan_points must be an integer >=17.');
end

Kprobe = zeros(nw,1);
Kprobe(1) = K0;
sig_intercept = invz_sigma_ordered(ctx.tl,lambda,Kprobe,g,ctx.beta);
d = (ctx.tl.M2/ctx.tl.n01^2) * (1-ctx.tl.n01^2) * g;
d(1) = 0; % K0 is supplied globally and already included in sig_intercept(1).

K = nan(nw,1);
K(1) = K0;
root_count = zeros(nw,1);
root_count(1) = 1;
status = strings(nw,1);
status(1) = "static_supplied";
proof = strings(nw,1);
proof(1) = "static_supplied";
root_residual = nan(nw,1);
derivative = nan(nw,1);
derivative_bound = nan(nw,1);
min_lattice_mass = nan(nw,1);
min_medium_mass = nan(nw,1);
interval_lo = nan(nw,1);
interval_hi = nan(nw,1);
fallback_unresolved = false(nw,1);

for n = 2:nw
    if retarded, J = ctx.Jnu_flat(:,n); else, J = Jshared; end
    rec = solve_frequency(J,G0(n),1+sig_intercept.Sigma(n),d(n));
    K(n) = rec.K;
    root_count(n) = rec.root_count;
    status(n) = rec.status;
    proof(n) = rec.proof;
    root_residual(n) = rec.root_residual;
    derivative(n) = rec.derivative;
    derivative_bound(n) = rec.derivative_bound;
    min_lattice_mass(n) = rec.min_lattice_mass;
    min_medium_mass(n) = rec.min_medium_mass;
    interval_lo(n) = rec.interval(1);
    interval_hi(n) = rec.interval(2);
    fallback_unresolved(n) = rec.fallback_unresolved;
end

defined = all(root_count(2:end) == 1) && all(isfinite(K(2:end))) && ...
    ~any(fallback_unresolved(2:end));
Sigma = nan(nw,1);
G = nan(nw,1);
lambda_check = nan(3,1);
affine_error = NaN;
dynamic_min_lattice_mass = NaN;
dynamic_min_medium_mass = NaN;
if defined
    sig = invz_sigma_ordered(ctx.tl,lambda,K,g,ctx.beta);
    Sigma = sig.Sigma(:);
    affine_error = max(abs(Sigma-(sig_intercept.Sigma(:)-d.*K)));
    denom_local = 1+Sigma+K.*G0;
    G = G0./denom_local;
    lambda_check = invz_lambdas(K,g,wts,ctx.beta,[1 2 3]);
    dynamic_min_lattice_mass = min(min_lattice_mass(2:end));
    dynamic_min_medium_mass = min(min_medium_mass(2:end));
end

if defined
    overall_status = "ok";
elseif any(root_count(2:end) > 1)
    overall_status = "dynamic_multiple_roots";
elseif any(fallback_unresolved(2:end))
    overall_status = "dynamic_search_unresolved";
elseif any(status(2:end) == "no_physical_domain")
    overall_status = "dynamic_no_physical_domain";
else
    overall_status = "dynamic_no_root";
end

out = struct('status',char(overall_status),'defined',defined,'K',K, ...
    'Sigma',Sigma,'G',G,'lambda_check',lambda_check, ...
    'sigma_intercept',sig_intercept.Sigma(:),'K_coefficient',d, ...
    'affine_error',affine_error,'root_count',root_count,'frequency_status',status, ...
    'proof',proof,'root_residual',root_residual,'derivative',derivative, ...
    'derivative_bound',derivative_bound,'min_lattice_mass',min_lattice_mass, ...
    'min_medium_mass',min_medium_mass,'interval_lo',interval_lo, ...
    'interval_hi',interval_hi,'fallback_unresolved',fallback_unresolved, ...
    'dynamic_min_lattice_mass',dynamic_min_lattice_mass, ...
    'dynamic_min_medium_mass',dynamic_min_medium_mass);

    function rec = solve_frequency(J,G0n,an,dn)
        J = J(:);
        jlo = min(J);
        jhi = max(J);
        span = max(jhi-jlo,eps(max(1,max(abs([jlo,jhi])))));
        lo = jlo;
        hi = jhi;
        if dn == 0
            L = an+J*G0n;
            if any(~isfinite(L)) || min(L) <= mass_tol
                rec = empty_record("no_physical_domain",[lo hi]);
                return;
            end
            invL = 1./L;
            kr = mean(J.*invL)/mean(invL);
            fr = kr-mean(J.*invL)/mean(invL);
            Dlocal = an+G0n*kr;
            Dq = L/Dlocal;
            rec = struct('K',kr,'root_count',1,'status',"ok", ...
                'proof',"K_independent",'root_residual',abs(fr),'derivative',0, ...
                'derivative_bound',0,'min_lattice_mass',min(L), ...
                'min_medium_mass',min(Dq),'interval',[lo hi], ...
                'fallback_unresolved',false);
            return;
        end
        minJG = min(J*G0n);
        mass_at_lo = an-dn*lo+minJG;
        if ~(isfinite(mass_at_lo) && mass_at_lo > mass_tol)
            rec = empty_record("no_physical_domain",[lo hi]);
            return;
        end

        full_interval = true;
        mass_at_hi = an-dn*hi+minJG;
        if dn > 0 && mass_at_hi <= mass_tol
            pole = (an+minJG-mass_tol)/dn;
            hi = min(hi,pole-emargin*span);
            full_interval = false;
        end
        if ~(isfinite(hi) && hi > lo)
            rec = empty_record("no_physical_domain",[lo hi]);
            return;
        end

        [~,~,Lhi] = scalar_fun(hi);
        Lmin = min(Lhi);
        Lmax = max(Lhi);
        if abs(G0n) <= realmin
            dbound = 0;
        else
            dbound = (dn/abs(G0n))*0.25*(Lmax/Lmin-1)^2;
        end
        proved_unique = full_interval && isfinite(dbound) && dbound < 1-uslack;

        if proved_unique
            [flo,~,~] = scalar_fun(lo);
            [fhi,~,~] = scalar_fun(hi);
            ftol = 64*eps(max(1,max(abs([jlo,jhi]))));
            if flo > ftol || fhi < -ftol
                rec = empty_record("bracket_inconsistent",[lo hi]);
                rec.derivative_bound = dbound;
                return;
            end
            if dbound < 0.5
                [kr,fr] = fixed_point_root(min(max(mean(J),lo),hi));
                mode = "contraction_bound";
            else
                [kr,fr] = bisect_root(lo,hi,flo,fhi);
                mode = "derivative_bound";
            end
            roots = kr;
            residuals = abs(fr);
            search_unresolved = false;
        else
            grid = linspace(lo,hi,nscan).';
            fgrid = nan(size(grid));
            for ig = 1:numel(grid), fgrid(ig) = scalar_fun(grid(ig)); end
            bropts = struct('resid_tol',rtol,'x_tol',xtol,'maxit',120, ...
                'merge_tol',max(10*xtol,1e-12*span),'fgrid',fgrid);
            [roots,search] = invz_bounded_roots(@scalar_fun,grid,bropts);
            residuals = search.root_residual;
            search_unresolved = search.unresolved_minima > 0;
            mode = "finite_enumeration";
        end

        if isscalar(roots) && ...
                (~isfinite(residuals(1)) || residuals(1) > rtol)
            roots = zeros(0,1);
            residuals = zeros(0,1);
            search_unresolved = true;
        end

        if numel(roots) ~= 1
            rec = empty_record("no_root",[lo hi]);
            rec.root_count = numel(roots);
            if numel(roots) > 1, rec.status = "multiple_roots"; end
            if search_unresolved, rec.status = "search_unresolved"; end
            rec.proof = mode;
            rec.derivative_bound = dbound;
            rec.fallback_unresolved = search_unresolved;
            return;
        end

        kr = roots(1);
        [fr,psip,L] = scalar_fun(kr);
        Dlocal = an+(G0n-dn)*kr;
        Dq = L/Dlocal;
        rec = struct('K',kr,'root_count',1,'status',"ok",'proof',mode, ...
            'root_residual',max(abs(fr),residuals(1)), ...
            'derivative',psip,'derivative_bound',dbound, ...
            'min_lattice_mass',min(L),'min_medium_mass',min(Dq), ...
            'interval',[lo hi],'fallback_unresolved',search_unresolved);

        function [f,psip,L] = scalar_fun(k)
            L = an-dn*k+J*G0n;
            if any(~isfinite(L)) || min(L) <= mass_tol
                f = NaN; psip = NaN; return;
            end
            invL = 1./L;
            S0 = mean(invL);
            S1 = mean(J.*invL);
            psi = S1/S0;
            f = k-psi;
            if nargout > 1
                S0p = dn*mean(invL.^2);
                S1p = dn*mean(J.*invL.^2);
                psip = (S1p*S0-S1*S0p)/(S0^2);
            end
        end

        function [kr,fr] = bisect_root(a,b,fa,fb)
            if abs(fa) <= rtol, kr=a; fr=fa; return; end
            if abs(fb) <= rtol, kr=b; fr=fb; return; end
            kr = a; fr = fa;
            for ib = 1:160
                m = a+(b-a)/2;
                fm = scalar_fun(m);
                if ~isfinite(fm), break; end
                if abs(fm) < abs(fr), kr=m; fr=fm; end
                if sign(fa) ~= sign(fm)
                    b=m; fb=fm; %#ok<NASGU>
                else
                    a=m; fa=fm;
                end
                if abs(fr) <= rtol || b-a <= xtol, break; end
            end
        end

        function [kr,fr] = fixed_point_root(k)
            kr = k;
            fcur = scalar_fun(k);
            fr = fcur;
            for ip = 1:80
                knew = k-fcur; % f(k)=k-Psi(k)
                fnew = scalar_fun(knew);
                if ~isfinite(fnew), break; end
                if abs(fnew) < abs(fr), kr=knew; fr=fnew; end
                if abs(fr) <= rtol, break; end
                k = knew;
                fcur = fnew;
            end
        end
    end

    function rec = empty_record(st,interval)
        rec = struct('K',NaN,'root_count',0,'status',st,'proof',"none", ...
            'root_residual',NaN,'derivative',NaN,'derivative_bound',NaN, ...
            'min_lattice_mass',NaN,'min_medium_mass',NaN,'interval',interval, ...
            'fallback_unresolved',false);
    end
end
