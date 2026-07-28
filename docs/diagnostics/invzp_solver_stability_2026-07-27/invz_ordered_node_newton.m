function [state, info] = invz_ordered_node_newton(node, seed, opts)
%INVZ_ORDERED_NODE_NEWTON Diagnostic residual-Newton corrector for one ordered node.
%   This is a solver-only alternative for the existing resummed equations. It does not
%   alter the closure, regularize lattice poles, or relax invz_ordered_residual's A-D
%   acceptance audit. It is intentionally outside the production MATLAB path.
%
%   IMPORTANT: acceptance proves only a fixed-h algebraic root. This helper does not establish
%   branch identity, direction independence, or a valid Jensen integration path. Never insert
%   its state into production HMF quadrature without a separate branch-tracked continuation
%   driver and the reproducibility gates in invzp_convg_fix.md.
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
%
%   opts.row_equilibrate (default false) scales each residual row by its
%   Jacobian infinity norm only for the Newton solve and rcond gate. Raw
%   residuals and the independent A--D audit remain unchanged.
%   opts.static_polish (default false) permits one physical-K0 Newton
%   proposal at each eligible Newton iterate after every Sigma residual has
%   passed. The rounded proposal is capped by
%   opts.static_polish_max_ulps (default 4096) and is retained only if the
%   complete raw residual passes without changing a tolerance. Acceptance
%   returns immediately, so at most one proposal can be retained.
%   opts.audit_formulation ('nested' default | 'coupled') selects the
%   independent audit's Block-A definition. The simultaneous solver may
%   use 'coupled' only with node.eso.audit_coordinate='defactored'.
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
row_equilibrate = logical_scalar( ...
    getf(opts, 'row_equilibrate', false), 'row_equilibrate');
static_polish = logical_scalar( ...
    getf(opts, 'static_polish', false), 'static_polish');
static_polish_max_ulps = positive_integer( ...
    getf(opts, 'static_polish_max_ulps', 4096), ...
    'static_polish_max_ulps');
audit_formulation = text_choice( ...
    getf(opts,'audit_formulation','nested'), ...
    {'nested','coupled'},'audit_formulation');
if strcmp(audit_formulation,'coupled') && ...
        ~strcmp(getf(node.eso,'audit_coordinate','raw_closure'),'defactored')
    error('invz:orderedNewtonOpts', ...
        ['opts.audit_formulation=''coupled'' requires ', ...
         'node.eso.audit_coordinate=''defactored''.']);
end
if ~(isscalar(maxit) && isfinite(maxit) && maxit >= 1 && maxit == floor(maxit))
    error('invz:orderedNewtonOpts', 'opts.maxit must be a positive integer.');
end
if ~(isscalar(max_linesearch) && isfinite(max_linesearch) && ...
        max_linesearch >= 1 && max_linesearch == floor(max_linesearch))
    error('invz:orderedNewtonOpts', ...
        'opts.max_linesearch must be a positive integer.');
end
validate_positive_finite(residual_trigger, 'residual_trigger');
validate_positive_finite(tol_outer, 'tol_outer');
validate_positive_finite(pole_margin_min, 'pole_margin_min');
validate_positive_finite(mean_margin_min, 'mean_margin_min');
validate_positive_finite(rcond_min, 'rcond_min');
config = struct( ...
    'schema','invz_ordered_node_newton_config/v1', ...
    'maxit',maxit,'max_linesearch',max_linesearch, ...
    'tol_outer',tol_outer,'residual_trigger',residual_trigger, ...
    'pole_margin_min',pole_margin_min, ...
    'mean_margin_min',mean_margin_min,'rcond_min',rcond_min, ...
    'row_equilibrate',row_equilibrate, ...
    'audit_formulation',audit_formulation, ...
    'static_polish',static_polish, ...
    'static_polish_max_ulps',static_polish_max_ulps);

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
polish_count = 0;
for iter = 1:maxit
    [R, diag, state] = invz_ordered_node_equations(node,u);
    resid_inf = norm(R, Inf);
    currentEvent = event_reason(R, diag, pole_margin_min, mean_margin_min);
    if ~isempty(currentEvent)
        reason = currentEvent;
        hrec = history_record(iter,resid_inf,diag,NaN,NaN,NaN,NaN);
        hrec.reason = reason;
        history(end+1) = hrec; %#ok<AGROW>
        break;
    end

    if resid_inf < residual_trigger
        audit = invz_ordered_residual(node, state, struct( ...
            'tol_outer',tol_outer,'formulation',audit_formulation));
        if audit.accepted
            [~,~,~,J] = invz_ordered_node_equations(node,u);
            [~,~,jac_rcond,raw_rcond,equilibrated_rcond,scale_valid] = ...
                linear_system(J,R,row_equilibrate);
            history(end+1) = history_record( ...
                iter,resid_inf,diag,jac_rcond,raw_rcond, ...
                equilibrated_rcond,NaN); %#ok<AGROW>
            if ~scale_valid
                reason = 'invalid_row_scale';
                history(end).reason = reason;
                break;
            end
            if ~(isfinite(jac_rcond) && jac_rcond > rcond_min)
                reason = 'rank_loss';
                history(end).reason = reason;
                break;
            end
            history(end).reason = 'accepted';
            info = accepted_info(node, state, audit, diag, iter, ...
                resid_inf,jac_rcond,raw_rcond,equilibrated_rcond, ...
                row_equilibrate,polish_count,history,config);
            return;
        end
        % The defactored residual and the independent factored A-D audit have different
        % conditioning near a cancellation. Continue correcting; never promote a vector
        % residual solely because it crossed residual_trigger.
    end

    [~,~,~,J] = invz_ordered_node_equations(node,u);
    [Jsolve,Rsolve,jac_rcond,raw_rcond,equilibrated_rcond,scale_valid] = ...
        linear_system(J,R,row_equilibrate);
    hrec = history_record(iter,resid_inf,diag,jac_rcond,raw_rcond, ...
        equilibrated_rcond,NaN);
    if ~scale_valid
        hrec.reason = 'invalid_row_scale';
        history(end+1) = hrec; %#ok<AGROW>
        reason = 'invalid_row_scale';
        break;
    end
    if static_polish && norm(R(1:end-1),Inf) <= residual_trigger && ...
            abs(R(end)) > residual_trigger
        K0 = u(end);
        K1 = K0-R(end)/J(end,end);
        physical_ulps = abs(K1-K0)/eps(abs(K0));
        if isfinite(K1) && isfinite(physical_ulps) && ...
                physical_ulps <= static_polish_max_ulps
            trial = u;
            trial(end) = K1;
            [Rt,dt,trial_state,Jt] = invz_ordered_node_equations(node,trial);
            if isempty(event_reason( ...
                    Rt,dt,pole_margin_min,mean_margin_min)) && ...
                    norm(Rt,Inf) < residual_trigger
                trial_audit = invz_ordered_residual( ...
                    node,trial_state,struct( ...
                    'tol_outer',tol_outer, ...
                    'formulation',audit_formulation));
                [~,~,trial_rcond,trial_raw_rcond,trial_eq_rcond, ...
                    trial_scale_valid] = ...
                    linear_system(Jt,Rt,row_equilibrate);
                if trial_audit.accepted && trial_scale_valid && ...
                        isfinite(trial_rcond) && trial_rcond > rcond_min
                    polish_count = polish_count+1;
                    history(end+1) = history_record( ...
                        iter,norm(Rt,Inf),dt,trial_rcond, ...
                        trial_raw_rcond,trial_eq_rcond,1); %#ok<AGROW>
                    history(end).reason = 'accepted_static_K0_polish';
                    state = trial_state;
                    info = accepted_info(node,state,trial_audit,dt,iter, ...
                        norm(Rt,Inf),trial_rcond,trial_raw_rcond, ...
                        trial_eq_rcond,row_equilibrate, ...
                        polish_count,history,config);
                    return
                end
            end
        end
    end
    if ~(isfinite(jac_rcond) && jac_rcond > rcond_min)
        hrec.reason = 'rank_loss';
        history(end+1) = hrec; %#ok<AGROW>
        reason = 'rank_loss';
        break;
    end
    step = -(Jsolve\Rsolve);
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
final_raw_rcond = NaN;
final_equilibrated_rcond = NaN;
if ~isempty(history)
    final_rcond = history(end).rcond;
    final_raw_rcond = history(end).raw_rcond;
    final_equilibrated_rcond = history(end).equilibrated_rcond;
end
info = struct('accepted', false, 'reason', reason, 'iterations', iter, ...
    'residual_inf', resid_inf, 'pole_margin', diag.pole_margin, ...
    'mean_margin', diag.mean_margin, 'rcond', final_rcond, ...
    'raw_rcond',final_raw_rcond, ...
    'equilibrated_rcond',final_equilibrated_rcond, ...
    'row_equilibrate',row_equilibrate, ...
    'static_polish_count',polish_count, ...
    'q_cancel', diag.q_cancel, 'Q', diag.Q, 'z', diag.z, ...
    'r', diag.r, 'Gtil0', diag.Gtil0, 'audit', audit, ...
    'history', history,'config',config);
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
        jac_rcond,raw_rcond,equilibrated_rcond,row_equilibrate, ...
        polish_count,history,config)
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
    'raw_rcond',raw_rcond, ...
    'equilibrated_rcond',equilibrated_rcond, ...
    'row_equilibrate',row_equilibrate, ...
    'static_polish_count',polish_count, ...
    'q_cancel', diag.q_cancel, 'Q', diag.Q, 'z', diag.z, ...
    'r', diag.r, 'Gtil0', diag.Gtil0, 'audit', audit, ...
    'history', history,'config',config);
end

function rec = history_record(iter,resid_inf,diag,jac_rcond,raw_rcond, ...
        equilibrated_rcond,alpha)
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
rec.raw_rcond = raw_rcond;
rec.equilibrated_rcond = equilibrated_rcond;
rec.alpha = alpha;
end

function rec = blank_history()
rec = struct('iter', 0, 'residual_inf', NaN, 'pole_margin', NaN, ...
    'mean_margin', NaN, 'q_cancel', NaN, 'Q', NaN, 'z', NaN, ...
    'r', NaN, 'Gtil0', NaN, 'gstat_local_denom', NaN, ...
    'rcond', NaN, 'raw_rcond',NaN,'equilibrated_rcond',NaN, ...
    'alpha', NaN, 'reason', '');
end

function [Jsolve,Rsolve,gate_rcond,raw_rcond,equilibrated_rcond,valid] = ...
        linear_system(J,R,row_equilibrate)
raw_rcond = rcond(J);
equilibrated_rcond = NaN;
valid = true;
if ~row_equilibrate
    Jsolve = J;
    Rsolve = R;
    gate_rcond = raw_rcond;
    return
end
row_norm = max(abs(J),[],2);
valid = all(isfinite(row_norm),'all') && all(row_norm > 0);
if ~valid
    Jsolve = nan(size(J));
    Rsolve = nan(size(R));
    gate_rcond = NaN;
    return
end
Jsolve = J./row_norm;
Rsolve = R./row_norm;
equilibrated_rcond = rcond(Jsolve);
gate_rcond = equilibrated_rcond;
end

function validate_positive_finite(value, name)
if ~(isscalar(value) && isfinite(value) && value > 0)
    error('invz:orderedNewtonOpts', ...
        'opts.%s must be a positive finite scalar.', name);
end
end

function value = positive_integer(value,name)
validate_positive_finite(value,name);
if value ~= floor(value)
    error('invz:orderedNewtonOpts', ...
        'opts.%s must be a positive integer.',name);
end
end

function value = logical_scalar(value,name)
if ~islogical(value) || ~isscalar(value)
    error('invz:orderedNewtonOpts', ...
        'opts.%s must be a scalar logical.',name);
end
end

function value = text_choice(value,choices,name)
if isstring(value)
    if ~isscalar(value) || ismissing(value)
        error('invz:orderedNewtonOpts', ...
            'opts.%s must be one of: %s.',name,strjoin(choices,', '));
    end
    value = char(value);
end
if ~(ischar(value) && isrow(value) && any(strcmp(value,choices)))
    error('invz:orderedNewtonOpts', ...
        'opts.%s must be one of: %s.',name,strjoin(choices,', '));
end
end

function value = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
