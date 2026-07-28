function cv = invzf_cutoff_convergence(model,lattice,m,h,H,cutoffs,opts)
%INVZF_CUTOFF_CONVERGENCE Geometric Matsubara-cutoff audit at fixed (m,h).
%
%   CV = INVZF_CUTOFF_CONVERGENCE(...,CUTOFFS,OPTS) evaluates the same
%   functional on WN=0:N for each N.  CUTOFFS must contain at least three
%   successive doublings.  The expected scalar ring tail is O(N^-3).
%
%   The audit reports Richardson tail estimates for f, u, both gradient
%   components, and the three independent Hessian components, plus the
%   rigorous analytic bound on the omitted ring free-energy tail.
%
%   opts.abs_tol (default 1e-10), opts.rel_tol (default 1e-8), and
%   opts.min_order (default 2.5) define status='converged'.  This helper
%   holds (m,h) fixed; stationary calculations must also repeat the root
%   solve at the accepted cutoff.

if nargin < 8 || isempty(opts), opts = struct(); end
cutoffs = cutoffs(:);
validateattributes(cutoffs, {'numeric'}, ...
    {'real','vector','finite','integer','positive'});
if numel(cutoffs) < 3 || any(diff(cutoffs) <= 0)
    error('invzf:cutoffGrid', 'cutoffs must contain at least three increasing values.');
end
if any(cutoffs(2:end) ~= 2*cutoffs(1:end-1))
    error('invzf:cutoffGrid', 'cutoffs must be successive doublings.');
end
atol = get_opt(opts,'abs_tol',1e-10);
rtol = get_opt(opts,'rel_tol',1e-8);
pmin = get_opt(opts,'min_order',2.5);
validateattributes([atol,rtol,pmin], {'numeric'}, ...
    {'real','vector','finite','nonnegative'});

names = {'f','u','grad_m','grad_h','hess_mm','hess_mh','hess_hh'};
values = nan(numel(cutoffs),numel(names));
tail_bound = nan(numel(cutoffs),1);
states = cell(numel(cutoffs),1);
for k = 1:numel(cutoffs)
    out = invzf_scalar_functional(model,lattice,m,h,H,(0:cutoffs(k)).');
    if ~strcmp(out.status,'ok')
        cv = struct('status',out.status,'cutoffs',cutoffs,'names',{names}, ...
            'values',values,'states',{states});
        return
    end
    values(k,:) = [out.f,out.u,out.grad(1),out.grad(2), ...
        out.hessian(1,1),out.hessian(1,2),out.hessian(2,2)];
    tail_bound(k) = out.ring.f_tail_bound;
    states{k} = out;
end

increments = diff(values,1,1);
tail_estimate = abs(increments(end,:))/(2^3-1);
prev = abs(increments(end-1,:));
last = abs(increments(end,:));
observed_order = inf(size(last));
nz = last > 0 & prev > 0;
observed_order(nz) = log2(prev(nz)./last(nz));
observed_order(last > 0 & prev == 0) = -Inf;
scale = max(1,abs(values(end,:)));
within = tail_estimate <= atol+rtol*scale;
order_ok = observed_order >= pmin | tail_estimate <= atol;
f_bound_ok = tail_bound(end) <= atol+rtol*max(1,abs(values(end,1)));

cv = struct('status','converged','cutoffs',cutoffs,'names',{names}, ...
    'values',values,'increments',increments,'tail_estimate',tail_estimate, ...
    'observed_order',observed_order,'analytic_f_tail_bound',tail_bound, ...
    'within_tolerance',within,'order_ok',order_ok,'states',{states});
if ~(all(within) && all(order_ok) && f_bound_ok)
    cv.status = 'not_converged';
end
end

function v = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end
