function sol = invzf_stationary_scalar(model, lattice, H, wn, opts)
%INVZF_STATIONARY_SCALAR Enumerate and compare strict scalar stationary states.
%
%   SOL = INVZF_STATIONARY_SCALAR(MODEL,LATTICE,H,WN,OPTS) searches the
%   preregistered interval opts.h_bounds for every resolved root of
%
%     R(h) = h-H-J0*(m_loc(h)-d f_ring/dh).
%
%   Sign-changing roots and resolved tangent roots are retained.  Every root
%   is evaluated with INVZF_SCALAR_FUNCTIONAL.  Local stability is the Schur
%   curvature after eliminating h,
%
%     f_red,mm = f_mm - f_mh^2/f_hh.
%
%   selected contains every stable root tied for the lowest functional value
%   within opts.energy_tol; symmetry-related degeneracy is not broken.
%
%   Required option:
%     h_bounds  finite increasing [hmin hmax]
%
%   Optional:
%     n_scan       initial uniform nodes (default 801)
%     residual_tol accepted |R| (default 1e-10)
%     gradient_tol accepted functional gradient norm (default residual_tol)
%     root_tol     root clustering/solver tolerance (default 1e-10)
%     energy_tol   minimum tie tolerance (default 1e-10)
%     stability_tol reduced-curvature floor (default 1e-10)

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'h_bounds') || numel(opts.h_bounds) ~= 2
    error('invzf:hBounds', 'opts.h_bounds=[hmin hmax] is required.');
end
bounds = opts.h_bounds(:).';
if any(~isfinite(bounds)) || bounds(2) <= bounds(1)
    error('invzf:hBounds', 'opts.h_bounds must be finite and increasing.');
end
nscan = get_opt(opts,'n_scan',801);
rtol = get_opt(opts,'residual_tol',1e-10);
gtol = get_opt(opts,'gradient_tol',rtol);
xtol = get_opt(opts,'root_tol',1e-10);
etol = get_opt(opts,'energy_tol',1e-10);
stol = get_opt(opts,'stability_tol',1e-10);
validateattributes(nscan, {'numeric'}, {'real','scalar','finite','integer','>=',5});
validateattributes([rtol,gtol,xtol,etol,stol], {'numeric'}, ...
    {'real','vector','finite','nonnegative'});
validateattributes(H, {'numeric'}, {'real','scalar','finite'});

if lattice.J0 == 0
    if H < bounds(1) || H > bounds(2)
        roots = repmat(empty_record(),0,1);
        sol = package(bounds,H,NaN,false,roots,[],rtol,xtol);
        return
    end
    ev = evaluate_h(H,model,lattice,H,wn);
    roots = repmat(empty_record(),0,1);
    if ev.valid
        rec = make_record(ev,model,lattice,H,wn,stol,gtol);
        if rec.stationary, roots = rec; end
    end
    selected = select_roots(roots,etol);
    sol = package(bounds,H,ev.residual,ev.valid,roots,selected,rtol,xtol);
    return
end

hscan = linspace(bounds(1),bounds(2),nscan).';
rscan = nan(size(hscan));
valid = false(size(hscan));
for k = 1:nscan
    ev = evaluate_h(hscan(k),model,lattice,H,wn);
    rscan(k) = ev.residual;
    valid(k) = ev.valid;
end

candidates = hscan(valid & abs(rscan) <= rtol);
for k = 1:nscan-1
    if ~(valid(k) && valid(k+1)), continue; end
    if rscan(k)*rscan(k+1) < 0
        fun = @(x) residual_only(x,model,lattice,H,wn);
        try
            candidates(end+1,1) = fzero(fun,[hscan(k),hscan(k+1)], ...
                optimset('TolX',xtol,'Display','off')); %#ok<AGROW>
        catch err
            if ~strcmp(err.identifier,'invzf:rootDomain'), rethrow(err); end
        end
    end
end

% A double/tangent stationary point has no sign change.  Search every
% resolved local minimum of |R| and accept it only if the original residual
% satisfies the same hard tolerance.
for k = 2:nscan-1
    if ~(valid(k-1) && valid(k) && valid(k+1)), continue; end
    if rscan(k-1)*rscan(k+1) > 0 ...
            && abs(rscan(k)) <= abs(rscan(k-1)) && abs(rscan(k)) <= abs(rscan(k+1))
        obj = @(x) abs(residual_or_large(x,model,lattice,H,wn));
        [hr,ar] = fminbnd(obj,hscan(k-1),hscan(k+1), ...
            optimset('TolX',xtol,'Display','off'));
        if ar <= rtol
            candidates(end+1,1) = hr; %#ok<AGROW>
        end
    end
end

candidates = cluster_roots(candidates,xtol);
records = repmat(empty_record(),0,1);
for k = 1:numel(candidates)
    ev = evaluate_h(candidates(k),model,lattice,H,wn);
    if ev.valid && abs(ev.residual) <= rtol
        rec = make_record(ev,model,lattice,H,wn,stol,gtol);
        if rec.stationary
            records(end+1,1) = rec; %#ok<AGROW>
        end
    end
end
selected = select_roots(records,etol);
sol = package(bounds,hscan,rscan,valid,records,selected,rtol,xtol);
end

function ev = evaluate_h(h,model,lattice,H,wn)
loc = invzf_twolevel_local(model.Delta,model.M,h,model.beta,wn);
if isfield(lattice,'mode_weights'), qw = lattice.mode_weights; else, qw = []; end
ring = invzf_ring_scalar(loc,lattice.Jmodes,qw);
ev = struct('h',h,'loc',loc,'ring',ring,'valid',strcmp(ring.status,'ok'), ...
    'residual',NaN);
if ev.valid
    ev.residual = h-H-lattice.J0*(loc.m-ring.dfdh);
end
end

function r = residual_only(h,model,lattice,H,wn)
ev = evaluate_h(h,model,lattice,H,wn);
if ~ev.valid
    error('invzf:rootDomain', 'Root trial left the ring domain.');
end
r = ev.residual;
end

function r = residual_or_large(h,model,lattice,H,wn)
ev = evaluate_h(h,model,lattice,H,wn);
if ev.valid, r = ev.residual; else, r = realmax('double')/4; end
end

function rec = make_record(ev,model,lattice,H,wn,stol,gtol)
h = ev.h;
m = ev.loc.m-ev.ring.dfdh;
fun = invzf_scalar_functional(model,lattice,m,h,H,wn);
fhh = fun.hessian(2,2);
if isfinite(fhh) && abs(fhh) > eps(max(1,abs(fhh)))
    redcurv = fun.hessian(1,1)-fun.hessian(1,2)*fun.hessian(2,1)/fhh;
else
    redcurv = NaN;
end
rec = empty_record();
rec.h = h;
rec.m = m;
rec.f = fun.f;
rec.residual = ev.residual;
rec.grad_norm = norm(fun.grad,inf);
rec.stationary = isfinite(rec.grad_norm) && rec.grad_norm <= gtol;
rec.reduced_curvature = redcurv;
rec.stable = rec.stationary && isfinite(redcurv) && redcurv > stol;
if rec.stable, rec.chi_uniform = 1/redcurv; else, rec.chi_uniform = NaN; end
rec.min_den = ev.ring.min_den;
rec.functional = fun;
end

function rec = empty_record
rec = struct('h',NaN,'m',NaN,'f',NaN,'residual',NaN,'grad_norm',NaN, ...
    'stationary',false,'reduced_curvature',NaN,'stable',false, ...
    'chi_uniform',NaN,'min_den',NaN,'functional',[]);
end

function idx = select_roots(records,etol)
idx = [];
if isempty(records), return; end
stable = find([records.stable] & isfinite([records.f]));
if isempty(stable), return; end
fmin = min([records(stable).f]);
idx = stable(abs([records(stable).f]-fmin) <= etol).';
end

function roots = cluster_roots(roots,tol)
roots = sort(roots(isfinite(roots)));
if isempty(roots), return; end
keep = roots(1);
for k = 2:numel(roots)
    scale = max([1,abs(roots(k)),abs(keep(end))]);
    if abs(roots(k)-keep(end)) > tol*scale
        keep(end+1,1) = roots(k); %#ok<AGROW>
    else
        keep(end) = 0.5*(keep(end)+roots(k));
    end
end
roots = keep;
end

function sol = package(bounds,hscan,rscan,valid,roots,selected,rtol,xtol)
sol = struct('status','ok','h_bounds',bounds,'h_scan',hscan, ...
    'residual_scan',rscan,'valid_scan',valid,'roots',roots, ...
    'selected',selected,'residual_tol',rtol,'root_tol',xtol);
if isempty(roots)
    sol.status = 'no_stationary_root';
elseif isempty(selected)
    sol.status = 'no_stable_root';
elseif numel(selected) > 1
    sol.status = 'degenerate_minima';
end
end

function v = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end
