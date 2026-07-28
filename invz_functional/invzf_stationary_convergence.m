function cv = invzf_stationary_convergence(model, lattice, H, cutoffs, solver_opts, opts)
%INVZF_STATIONARY_CONVERGENCE Repeat the complete root solve versus cutoff.
%
%   CV = INVZF_STATIONARY_CONVERGENCE(...,CUTOFFS,SOLVER_OPTS,OPTS)
%   reruns root enumeration on wn=0:N for every increasing cutoff.  It
%   compares all selected stable minima after sorting them by h; it never
%   assumes that a warm-started branch is the thermodynamic state.
%
%   This complements invzf_cutoff_convergence, which holds (m,h) fixed.
%   The default acceptance requires the last two cutoffs to have the same
%   selected-root and total-root counts, componentwise convergence of
%   [h,m,f,reduced_curvature,min_den], and an accepted analytic free-energy
%   tail bound at every selected last-cutoff state.
%
%   opts.abs_tol (1e-8), opts.rel_tol (1e-6), opts.require_root_count (true).

if nargin < 6 || isempty(opts), opts = struct(); end
validateattributes(cutoffs, {'numeric'}, ...
    {'real','vector','finite','integer','positive'});
cutoffs = cutoffs(:);
if numel(cutoffs) < 2 || any(diff(cutoffs) <= 0)
    error('invzf:stationaryCutoffs', ...
        'cutoffs must contain at least two increasing positive integers.');
end
atol = get_opt(opts,'abs_tol',1e-8);
rtol = get_opt(opts,'rel_tol',1e-6);
require_count = get_opt(opts,'require_root_count',true);
validateattributes([atol rtol], {'numeric'}, ...
    {'real','vector','finite','nonnegative'});
valid_require_count = isscalar(require_count) && (islogical(require_count) || ...
    (isnumeric(require_count) && isreal(require_count) ...
     && isfinite(require_count) && any(require_count == [0 1])));
if ~valid_require_count
    error('invzf:stationaryConvergenceOpts', ...
        'require_root_count must be scalar logical or numeric 0/1.');
end
require_count = logical(require_count);

ncut = numel(cutoffs);
solutions = cell(ncut,1);
selected = cell(ncut,1);
root_counts = zeros(ncut,1);
selected_counts = zeros(ncut,1);
for k = 1:ncut
    solutions{k} = invzf_stationary_scalar( ...
        model,lattice,H,(0:cutoffs(k)).',solver_opts);
    root_counts(k) = numel(solutions{k}.roots);
    selected{k} = selected_table(solutions{k});
    selected_counts(k) = size(selected{k},1);
end

names = {'h','m','f','reduced_curvature','min_den'};
last_delta = nan(0,numel(names));
last_threshold = nan(0,numel(names));
same_selected_count = selected_counts(end) == selected_counts(end-1) ...
    && selected_counts(end) > 0;
if same_selected_count
    a = selected{end-1};
    b = selected{end};
    last_delta = abs(b-a);
    last_threshold = atol+rtol*max(1,abs(b));
end
value_ok = same_selected_count && all(last_delta <= last_threshold,'all');
root_count_ok = ~require_count || root_counts(end) == root_counts(end-1);

tail_bounds = nan(selected_counts(end),1);
if selected_counts(end) > 0
    rec = solutions{end}.roots(solutions{end}.selected);
    for k = 1:numel(rec)
        tail_bounds(k) = rec(k).functional.ring.f_tail_bound;
    end
end
fscale = 1;
if ~isempty(selected{end}), fscale = max(1,max(abs(selected{end}(:,3)))); end
tail_ok = ~isempty(tail_bounds) && all(isfinite(tail_bounds)) ...
    && all(tail_bounds <= atol+rtol*fscale);

cv = struct('status','converged','cutoffs',cutoffs, ...
    'names',{names},'solutions',{solutions},'selected',{selected}, ...
    'root_counts',root_counts,'selected_counts',selected_counts, ...
    'last_delta',last_delta,'last_threshold',last_threshold, ...
    'value_ok',value_ok,'root_count_ok',root_count_ok, ...
    'tail_bounds',tail_bounds,'tail_ok',tail_ok);
if ~(value_ok && root_count_ok && tail_ok)
    cv.status = 'not_converged';
end
end

function a = selected_table(sol)
a = nan(0,5);
if isempty(sol.selected), return; end
r = sol.roots(sol.selected);
a = [[r.h].',[r.m].',[r.f].',[r.reduced_curvature].',[r.min_den].'];
[~,ix] = sort(a(:,1));
a = a(ix,:);
end

function v = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end
