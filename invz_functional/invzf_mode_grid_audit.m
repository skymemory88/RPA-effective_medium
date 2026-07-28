function audit = invzf_mode_grid_audit(ion, T, B, grids, cutoff, opts)
%INVZF_MODE_GRID_AUDIT Compare strict stationary states across production BZ grids.
%
%   AUDIT = INVZF_MODE_GRID_AUDIT(ION,T,B,GRIDS,CUTOFF,OPTS) recomputes the
%   scalar production coupling spectrum on each cubic N in GRIDS, maps it
%   through invzf_projected_inputs, and repeats complete root enumeration.
%   Production solver files and dispatch are not modified.
%
%   opts.coupling is forwarded to invz_bz_couplings except that .grid is
%   owned by this audit.  opts.adapter is forwarded to
%   invzf_projected_inputs.  opts.solver is forwarded to the root solver;
%   absent h_bounds defaults to H +/-1.5*abs(J0)*ion.J independently at
%   each grid.  opts.abs_tol (1e-6), opts.rel_tol (1e-4),
%   opts.min_den_floor (0), and opts.require_root_count (true) grade only
%   the final two grids.  The extremal min_den is reported but is not
%   treated as a BZ quadrature observable: a mesh approaching an excluded
%   Gamma point need not make that extremum converge like a normalized
%   q-average.  It is instead subject to the explicit domain floor.

if nargin < 6 || isempty(opts), opts = struct(); end
validateattributes(grids, {'numeric'}, ...
    {'real','vector','finite','integer','>=',2});
grids = grids(:);
if numel(grids) < 2 || any(diff(grids) <= 0)
    error('invzf:modeGrids', 'grids must contain at least two increasing sizes.');
end
validateattributes(cutoff, {'numeric'}, ...
    {'real','scalar','finite','integer','positive'});
coupling_opts = get_opt(opts,'coupling',struct());
adapter_opts = get_opt(opts,'adapter',struct());
solver_base = get_opt(opts,'solver',struct());
if isfield(coupling_opts,'grid')
    error('invzf:modeGridOwner', 'opts.coupling.grid is owned by the grid audit.');
end
atol = get_opt(opts,'abs_tol',1e-6);
rtol = get_opt(opts,'rel_tol',1e-4);
min_den_floor = get_opt(opts,'min_den_floor',0);
require_count = get_opt(opts,'require_root_count',true);
validateattributes([atol rtol min_den_floor], {'numeric'}, ...
    {'real','vector','finite','nonnegative'});
valid_require_count = isscalar(require_count) && (islogical(require_count) || ...
    (isnumeric(require_count) && isreal(require_count) ...
     && isfinite(require_count) && any(require_count == [0 1])));
if ~valid_require_count
    error('invzf:modeGridOpts', ...
        'require_root_count must be scalar logical or numeric 0/1.');
end
require_count = logical(require_count);

ng = numel(grids);
entries = repmat(empty_entry(),ng,1);
selected = cell(ng,1);
root_counts = zeros(ng,1);
selected_counts = zeros(ng,1);
for k = 1:ng
    co = coupling_opts;
    co.grid = repmat(grids(k),1,3);
    [Jnu,info,Jaa0] = invz_bz_couplings(ion,co);
    pilot = invzf_projected_inputs(ion,T,B,Jnu,info,Jaa0,adapter_opts);
    so = solver_base;
    if ~isfield(so,'h_bounds') || isempty(so.h_bounds)
        hspan = 1.5*abs(pilot.lattice.J0)*ion.J;
        if hspan == 0
            hspan = max(1e-8,abs(pilot.H)*1e-6);
        end
        so.h_bounds = pilot.H+[-hspan hspan];
    end
    sol = invzf_stationary_scalar( ...
        pilot.model,pilot.lattice,pilot.H,(0:cutoff).',so);
    tab = selected_table(sol);
    entries(k) = struct('grid',grids(k),'n_modes',numel(Jnu), ...
        'Jmin',min(Jnu),'Jmax',max(Jnu),'Jmean',mean(Jnu), ...
        'Jmu2',mean(Jnu.^2),'J0',info.Jcc0,'Jaa0',Jaa0, ...
        'pilot',pilot,'solution',sol);
    selected{k} = tab;
    root_counts(k) = numel(sol.roots);
    selected_counts(k) = size(tab,1);
end

names = {'h','m','f','reduced_curvature','min_den'};
graded_names = names(1:4);
same_selected_count = selected_counts(end) == selected_counts(end-1) ...
    && selected_counts(end) > 0;
last_delta = nan(0,numel(names));
last_threshold = nan(0,numel(names));
if same_selected_count
    last_delta = abs(selected{end}-selected{end-1});
    last_threshold = atol+rtol*max(1,abs(selected{end}));
end
value_ok = same_selected_count ...
    && all(last_delta(:,1:numel(graded_names)) ...
           <= last_threshold(:,1:numel(graded_names)),'all');
root_count_ok = ~require_count || root_counts(end) == root_counts(end-1);
domain_ok = ~isempty(selected{end}) ...
    && all(selected{end}(:,5) > min_den_floor);

audit = struct('status','converged','T',T,'B',invz_field_vec(B), ...
    'grids',grids,'cutoff',cutoff,'names',{names}, ...
    'graded_names',{graded_names},'entries',entries, ...
    'selected',{selected},'root_counts',root_counts, ...
    'selected_counts',selected_counts,'last_delta',last_delta, ...
    'last_threshold',last_threshold,'value_ok',value_ok, ...
    'root_count_ok',root_count_ok,'min_den_floor',min_den_floor, ...
    'domain_ok',domain_ok);
if ~(value_ok && root_count_ok && domain_ok)
    audit.status = 'not_converged';
end
end

function e = empty_entry
e = struct('grid',NaN,'n_modes',NaN,'Jmin',NaN,'Jmax',NaN, ...
    'Jmean',NaN,'Jmu2',NaN,'J0',NaN,'Jaa0',NaN, ...
    'pilot',[],'solution',[]);
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
