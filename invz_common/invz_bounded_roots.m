function [roots, info] = invz_bounded_roots(fun, grid, opts)
%INVZ_BOUNDED_ROOTS Enumerate isolated scalar roots on a bounded sampled interval.
%   ROOTS = INVZ_BOUNDED_ROOTS(FUN, GRID, OPTS) searches every adjacent GRID
%   interval rather than stopping after the first sign change. It also refines
%   local minima of abs(FUN) so an isolated even-multiplicity root need not
%   change sign. A sign-changing bracket is accepted only when refinement
%   reaches opts.resid_tol at a finite point; otherwise it is reported as a
%   discontinuity bracket, never as a root.
%
%   GRID must be finite and strictly increasing. FUN accepts a scalar and
%   returns a scalar; NaN/Inf marks a point outside the searchable continuous
%   domain. OPTS:
%     resid_tol       accepted abs(FUN(root)); default 1e-10
%     x_tol           refinement interval tolerance; default 1e-12
%     maxit           bisection iterations; default 100
%     merge_tol       root deduplication distance; default 10*x_tol
%     fgrid           optional precomputed FUN(GRID) values
%     edge_ok         optional logical vector, length numel(GRID)-1. False
%                     prevents refinement across a known domain boundary.
%
%   INFO records fgrid, residuals, discontinuity brackets, tangency attempts,
%   refinement evaluation count, and the explicit grid resolution. This is a
%   finite-resolution enumerator, not an interval-arithmetic proof that an
%   arbitrarily narrow unresolved root cannot exist.
if nargin < 3, opts = struct(); end
if ~(isa(fun, 'function_handle'))
    error('invz:boundedRoots', 'fun must be a function handle.');
end
grid = grid(:);
if numel(grid) < 3 || any(~isfinite(grid)) || any(diff(grid) <= 0)
    error('invz:boundedRoots', 'grid must contain at least three finite, strictly increasing points.');
end
rtol = getf(opts, 'resid_tol', 1e-10);
xtol = getf(opts, 'x_tol', 1e-12);
maxit = getf(opts, 'maxit', 100);
merge_tol = getf(opts, 'merge_tol', 10*xtol);
refine_evals = 0;
if isfield(opts, 'fgrid') && ~isempty(opts.fgrid)
    fgrid = opts.fgrid(:);
    if numel(fgrid) ~= numel(grid)
        error('invz:boundedRoots', 'opts.fgrid must have one value per grid point.');
    end
else
    fgrid = nan(size(grid));
    for k = 1:numel(grid), fgrid(k) = safe_eval(grid(k)); end
end
if isfield(opts, 'edge_ok') && ~isempty(opts.edge_ok)
    edge_ok = logical(opts.edge_ok(:));
    if numel(edge_ok) ~= numel(grid)-1
        error('invz:boundedRoots', 'opts.edge_ok must have numel(grid)-1 elements.');
    end
else
    edge_ok = true(numel(grid)-1, 1);
end

cand_x = zeros(0,1);
cand_f = zeros(0,1);
disc = zeros(0,2);
sign_attempts = 0;
tangent_attempts = 0;
unresolved_minima = 0;

% Exact/near-exact sampled roots.
idx = find(isfinite(fgrid) & abs(fgrid) <= rtol);
for k = 1:numel(idx)
    j = idx(k);
    xs = grid(j); fs = fgrid(j);
    % A small residual can locate an even root only to O(sqrt(resid_tol)).
    % Polish interior sampled hits by minimizing |f| on their local cell.
    if j > 1 && j < numel(grid) && edge_ok(j-1) && edge_ok(j)
        oo = optimset('Display', 'off', 'TolX', xtol, 'MaxIter', maxit);
        [xp, ~, exitflag] = fminbnd(@abs_eval, grid(j-1), grid(j+1), oo);
        fp = safe_eval(xp);
        if exitflag > 0 && isfinite(fp) && abs(fp) < abs(fs)
            xs = xp; fs = fp;
        end
    end
    add_candidate(xs, fs);
end

% Refine every finite sign-changing interval that is not a known domain edge.
for k = 1:numel(grid)-1
    if ~edge_ok(k) || ~isfinite(fgrid(k)) || ~isfinite(fgrid(k+1)), continue; end
    if fgrid(k) == 0 || fgrid(k+1) == 0, continue; end
    if sign(fgrid(k)) == sign(fgrid(k+1)), continue; end
    sign_attempts = sign_attempts + 1;
    [xr, fr, ok] = refine_sign(grid(k), grid(k+1), fgrid(k), fgrid(k+1));
    if ok
        add_candidate(xr, fr);
    else
        disc(end+1,:) = [grid(k), grid(k+1)]; %#ok<AGROW>
    end
end

% Search every sampled local minimum of |f| for a touching root.
af = abs(fgrid);
for k = 2:numel(grid)-1
    if ~edge_ok(k-1) || ~edge_ok(k), continue; end
    if ~all(isfinite(fgrid(k-1:k+1))), continue; end
    % A root already bracketed by either adjacent sign-changing edge is not
    % a tangency candidate. Re-running fminbnd on |f| can produce a second,
    % less accurately polished representation of the same simple root.
    if sign(fgrid(k-1)) ~= sign(fgrid(k)) || ...
            sign(fgrid(k)) ~= sign(fgrid(k+1))
        continue;
    end
    if ~(af(k) <= af(k-1) && af(k) <= af(k+1)), continue; end
    if af(k) <= rtol, continue; end
    tangent_attempts = tangent_attempts + 1;
    oo = optimset('Display', 'off', 'TolX', xtol, 'MaxIter', maxit);
    [xm, ~, exitflag] = fminbnd(@abs_eval, grid(k-1), grid(k+1), oo);
    fm = safe_eval(xm);
    if exitflag > 0 && isfinite(fm) && abs(fm) <= rtol
        add_candidate(xm, fm);
    elseif exitflag <= 0
        unresolved_minima = unresolved_minima + 1;
    end
end

% Sort and deduplicate while retaining the lower-residual representative.
if isempty(cand_x)
    roots = zeros(0,1);
    root_resid = zeros(0,1);
else
    [cand_x, order] = sort(cand_x);
    cand_f = cand_f(order);
    keep_x = cand_x(1);
    keep_f = cand_f(1);
    for k = 2:numel(cand_x)
        if cand_x(k) - keep_x(end) <= merge_tol
            if abs(cand_f(k)) < abs(keep_f(end))
                keep_x(end) = cand_x(k);
                keep_f(end) = cand_f(k);
            end
        else
            keep_x(end+1,1) = cand_x(k); %#ok<AGROW>
            keep_f(end+1,1) = cand_f(k); %#ok<AGROW>
        end
    end
    roots = keep_x(:);
    root_resid = abs(keep_f(:));
end

info = struct('fgrid', fgrid, 'root_residual', root_resid, ...
    'discontinuity_brackets', disc, 'sign_attempts', sign_attempts, ...
    'tangent_attempts', tangent_attempts, 'unresolved_minima', unresolved_minima, ...
    'refine_evals', refine_evals, 'grid_points', numel(grid), ...
    'max_grid_gap', max(diff(grid)), 'resid_tol', rtol, 'x_tol', xtol, ...
    'merge_tol', merge_tol, 'finite_resolution', true);

    function y = safe_eval(x)
        refine_evals = refine_evals + 1;
        y = fun(x);
        if ~(isnumeric(y) && isreal(y) && isscalar(y))
            y = NaN;
        end
    end

    function y = abs_eval(x)
        fx = safe_eval(x);
        if isfinite(fx), y = abs(fx);
        else,            y = realmax/4;
        end
    end

    function add_candidate(x, fx)
        cand_x(end+1,1) = x;
        cand_f(end+1,1) = fx;
    end

    function [xr, fr, ok] = refine_sign(a, b, fa, fb)
        xr = NaN; fr = NaN; ok = false;
        best_x = a; best_f = fa;
        if abs(fb) < abs(best_f), best_x = b; best_f = fb; end
        for it = 1:maxit
            m = a + (b-a)/2;
            fm = safe_eval(m);
            if ~isfinite(fm), break; end
            if abs(fm) < abs(best_f), best_x = m; best_f = fm; end
            if sign(fa) ~= sign(fm)
                b = m; fb = fm; %#ok<NASGU>
            else
                a = m; fa = fm;
            end
            if b-a <= xtol, break; end
        end
        if isfinite(best_f) && abs(best_f) <= rtol
            xr = best_x; fr = best_f; ok = true;
        end
    end
end
