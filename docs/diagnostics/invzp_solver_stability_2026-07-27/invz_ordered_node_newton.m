function [state, info] = invz_ordered_node_newton(node, seed, opts)
%INVZ_ORDERED_NODE_NEWTON Diagnostic residual-Newton corrector for one ordered node.
%   This is a solver-only alternative for the existing resummed equations. It does not
%   alter the closure, regularize lattice poles, or relax invz_ordered_residual's A-D
%   acceptance audit. It is intentionally outside the production MATLAB path.
%
%   IMPORTANT: acceptance proves only a fixed-h algebraic root. This helper does not establish
%   branch identity, direction independence, or a valid Jensen integration path. Never insert
%   its state into production HMF quadrature without a separate branch-tracked continuation
%   driver and the reproducibility gates in invzp_convg_fix_Claude.md.
%
%   Unknowns are u = [Sigma(:); K0]. Dynamic K(2:end) and lambda are derived directly.
%   The static equation is the defactored closure K0 = sum(Jq*Gq)/sum(Gq), evaluated in
%   the reciprocal chart z = 1/Gstat. This removes the coordinate degeneracy that makes
%   the factored residual vanish identically at Gstat = 0, without changing the equation
%   away from that limit.
%
%   The dense analytic Jacobian is intended for the current diagnostic/experimental
%   temperature range. A structured solve can be added later if profiling at low T shows
%   it is necessary.
if nargin < 2, seed = []; end
if nargin < 3 || isempty(opts), opts = struct(); end

if ~isvector(node.Jnu_flat)
    error('invz:orderedNewtonCoupling', ...
        'Residual Newton currently requires a static coupling vector.');
end
if ~strcmp(getf(node.eso, 'static_medium', 'resummed'), 'resummed')
    error('invz:orderedNewtonScheme', ...
        'Residual Newton currently supports static_medium = ''resummed'' only.');
end

maxit = getf(opts, 'maxit', 12);
max_linesearch = getf(opts, 'max_linesearch', 12);
tol_outer = getf(opts, 'tol_outer', 1e-8);
residual_trigger = getf(opts, 'residual_trigger', tol_outer);
pole_margin_min = getf(opts, 'pole_margin_min', 1e-10);
mean_margin_min = getf(opts, 'mean_margin_min', 1e-10);
rcond_min = getf(opts, 'rcond_min', 1e-12);
if ~(isscalar(maxit) && isfinite(maxit) && maxit >= 1 && maxit == floor(maxit))
    error('invz:orderedNewtonOpts', 'opts.maxit must be a positive integer.');
end
if ~(isscalar(max_linesearch) && isfinite(max_linesearch) && ...
        max_linesearch >= 1 && max_linesearch == floor(max_linesearch))
    error('invz:orderedNewtonOpts', ...
        'opts.max_linesearch must be a positive integer.');
end
validate_positive_finite(residual_trigger, 'residual_trigger');
validate_positive_finite(pole_margin_min, 'pole_margin_min');
validate_positive_finite(mean_margin_min, 'mean_margin_min');
validate_positive_finite(rcond_min, 'rcond_min');

nw = numel(node.wn);
if isstruct(seed) && isfield(seed, 'Sigma') && isfield(seed, 'K0s') && ...
        ~isempty(seed.Sigma)
    u = [seed.Sigma(:); seed.K0s];
else
    u = [zeros(nw, 1); 0];
end
if numel(u) ~= nw+1
    error('invz:orderedNewtonSeed', ...
        'Newton seed has %d unknowns; expected %d.', numel(u), nw+1);
end

reason = 'max_newton';
audit = [];
history = repmat(blank_history(), 1, 0);
for iter = 1:maxit
    [R, diag, state] = invz_ordered_node_equations(node,u);
    resid_inf = norm(R, Inf);
    reason = event_reason(R, diag, pole_margin_min, mean_margin_min);
    if ~isempty(reason)
        hrec = history_record(iter, resid_inf, diag, NaN, NaN);
        hrec.reason = reason;
        history(end+1) = hrec; %#ok<AGROW>
        break;
    end

    if resid_inf < residual_trigger
        audit = invz_ordered_residual(node, state, struct('tol_outer', tol_outer));
        if audit.accepted
            [~,~,~,J] = invz_ordered_node_equations(node,u);
            jac_rcond = rcond(J);
            history(end+1) = history_record( ...
                iter, resid_inf, diag, jac_rcond, NaN); %#ok<AGROW>
            if ~(isfinite(jac_rcond) && jac_rcond > rcond_min)
                reason = 'rank_loss';
                history(end).reason = reason;
                break;
            end
            history(end).reason = 'accepted';
            info = accepted_info(node, state, audit, diag, iter, ...
                resid_inf, jac_rcond, history);
            return;
        end
        % The defactored residual and the independent factored A-D audit have different
        % conditioning near a cancellation. Continue correcting; never promote a vector
        % residual solely because it crossed residual_trigger.
    end

    [~,~,~,J] = invz_ordered_node_equations(node,u);
    jac_rcond = rcond(J);
    hrec = history_record(iter, resid_inf, diag, jac_rcond, NaN);
    if ~(isfinite(jac_rcond) && jac_rcond > rcond_min)
        hrec.reason = 'rank_loss';
        history(end+1) = hrec; %#ok<AGROW>
        reason = 'rank_loss';
        break;
    end
    step = -(J\R);
    if ~all(isfinite(step))
        hrec.reason = 'nonfinite_step';
        history(end+1) = hrec; %#ok<AGROW>
        reason = 'nonfinite_step';
        break;
    end
    base_norm = norm(R, 2);
    alpha = 1;
    accepted_step = false;
    for ls = 1:max_linesearch
        trial = u + alpha*step;
        [Rt,dt] = invz_ordered_node_equations(node,trial);
        if isempty(event_reason(Rt, dt, pole_margin_min, mean_margin_min)) && ...
                norm(Rt, 2) < base_norm
            u = trial;
            accepted_step = true;
            break;
        end
        alpha = alpha/2;
    end
    hrec.alpha = alpha;
    if accepted_step
        hrec.reason = 'step_accepted';
    else
        hrec.reason = 'line_search_failed';
    end
    history(end+1) = hrec; %#ok<AGROW>
    if ~accepted_step
        reason = 'line_search_failed';
        break;
    end
end

[R,diag,state] = invz_ordered_node_equations(node,u);
resid_inf = norm(R, Inf);
final_rcond = NaN;
if ~isempty(history), final_rcond = history(end).rcond; end
info = struct('accepted', false, 'reason', reason, 'iterations', iter, ...
    'residual_inf', resid_inf, 'pole_margin', diag.pole_margin, ...
    'mean_margin', diag.mean_margin, 'rcond', final_rcond, ...
    'q_cancel', diag.q_cancel, 'Q', diag.Q, 'z', diag.z, ...
    'r', diag.r, 'Gtil0', diag.Gtil0, 'audit', audit, ...
    'history', history);
end

function reason = event_reason(R, diag, pole_margin_min, mean_margin_min)
reason = '';
if ~all(isfinite(R))
    reason = 'nonfinite_residual';
elseif ~isfinite(diag.r)
    reason = 'unbounded_integrand';
elseif ~isfinite(diag.Gtil0)
    reason = 'nonfinite_Gtil0';
elseif ~isfinite(diag.pole_margin) || diag.pole_margin <= pole_margin_min
    reason = 'pole_event';
elseif ~isfinite(diag.mean_margin) || diag.mean_margin <= mean_margin_min
    reason = 'mean_event';
end
end

function info = accepted_info(node, state, audit, diag, iterations, resid_inf, ...
        jac_rcond, history)
med = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts);
Jf = node.Jnu_flat(:);
Dq = 1+(Jf-state.K0s)*diag.Gstat;
medium = struct('scheme', 'resummed', 'status', 'not_applicable', ...
    'ref', [], 'closure', []);
so = struct('xi', diag.xi, 'h0', diag.h0, 'G0bare', diag.G0bare, ...
    'Gtil0', diag.Gtil0, 'r', diag.r, ...
    'gstat_local_denom', diag.gstat_local_denom, ...
    'D_uni', audit.D_uni, 'Dq_min', min(Dq), 'Dq_max', max(Dq), ...
    'Dq_neg_count', nnz(Dq <= 0), 'iters', 0, 'medium', medium, ...
    'medium_status', 'not_applicable', 'omit_mu3', NaN, ...
    'omit_cubic', NaN, 'omit_max', NaN, ...
    'resid', audit.blockB.resid, 'converged', audit.blockB.pass, ...
    'Gstat', diag.Gstat);
info = struct('res', audit, 'loop_converged', false, 'so', so, ...
    'med', med, 'outer_iters', iterations, 'term_reason', 'accepted_newton', ...
    'iters', struct([]), 'medium', medium, ...
    'medium_status', 'not_applicable', 'accepted', true, ...
    'reason', 'accepted', 'iterations', iterations, ...
    'residual_inf', resid_inf, 'pole_margin', diag.pole_margin, ...
    'mean_margin', diag.mean_margin, 'rcond', jac_rcond, ...
    'q_cancel', diag.q_cancel, 'Q', diag.Q, 'z', diag.z, ...
    'r', diag.r, 'Gtil0', diag.Gtil0, 'audit', audit, ...
    'history', history);
end

function rec = history_record(iter, resid_inf, diag, jac_rcond, alpha)
rec = blank_history();
rec.iter = iter;
rec.residual_inf = resid_inf;
rec.pole_margin = diag.pole_margin;
rec.mean_margin = diag.mean_margin;
rec.q_cancel = diag.q_cancel;
rec.Q = diag.Q;
rec.z = diag.z;
rec.r = diag.r;
rec.Gtil0 = diag.Gtil0;
rec.gstat_local_denom = diag.gstat_local_denom;
rec.rcond = jac_rcond;
rec.alpha = alpha;
end

function rec = blank_history()
rec = struct('iter', 0, 'residual_inf', NaN, 'pole_margin', NaN, ...
    'mean_margin', NaN, 'q_cancel', NaN, 'Q', NaN, 'z', NaN, ...
    'r', NaN, 'Gtil0', NaN, 'gstat_local_denom', NaN, ...
    'rcond', NaN, 'alpha', NaN, 'reason', '');
end

function validate_positive_finite(value, name)
if ~(isscalar(value) && isfinite(value) && value > 0)
    error('invz:orderedNewtonOpts', ...
        'opts.%s must be a positive finite scalar.', name);
end
end

function value = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
