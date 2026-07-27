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
    [R, diag, state] = reciprocal_residual(u, node);
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
            J = analytic_jacobian(u, node);
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

    J = analytic_jacobian(u, node);
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
        [Rt, dt] = reciprocal_residual(trial, node);
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

[R, diag, state] = reciprocal_residual(u, node);
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

function [R, diag, state] = reciprocal_residual(u, node)
nw = numel(node.wn);
Sigma = u(1:nw);
K0 = u(end);
med = invz_emt_scalar(node.G0, Sigma, node.Jnu_flat, node.eopts);
K = [K0; med.K(2:end)];
lam = invz_lambdas(K, node.g, node.wts, node.beta, [1 2 3]);
sg = invz_sigma_ordered(node.tl, lam, K, node.g, node.beta);
[Gstat, go] = invz_gstat_ordered(node.tl, lam, K0, Sigma(1), node.beta, ...
    node.G0inel0, node.G0el0, struct('stable_form', true));

d0 = go.gstat_local_denom;
Hz = node.G0inel0 + go.xi*node.G0el0*d0;
z = d0/Hz;
Gtil0 = 1/(z-K0);
r = go.G0bare*(z-K0);

Jf = node.Jnu_flat(:);
Jscale = max(abs(Jf));
if ~(isfinite(Jscale) && Jscale > 0)
    error('invz:orderedNewtonCoupling', ...
        'The static coupling scale must be positive and finite.');
end
if isinf(z)
    Gbar = 0;
    Jloc = mean(Jf);
    pole_margin = 1;
    mean_margin = 0;
    q_cancel = 1;
    Q = 1;
else
    scale = max(abs(z), Jscale);
    E = z+Jf-K0;
    weights = scale./E;
    mean_weights = mean(weights);
    Gbar = mean_weights/scale;
    Jloc = mean(Jf.*weights)/mean_weights;
    pole_margin = min(abs(E))/scale;
    mean_margin = abs(Gbar)*Jscale;
    q_cancel = abs(mean_weights)/mean(abs(weights));
    Q = z*Gbar;
end

R = [sg.Sigma-Sigma; (K0-Jloc)/Jscale];
state = struct('Sigma', Sigma, 'K', K, 'lam', lam, 'K0s', K0);
diag = struct('z', z, 'Gstat', Gstat, 'Gtil0', Gtil0, 'r', r, ...
    'Gbar', Gbar, 'pole_margin', pole_margin, ...
    'mean_margin', mean_margin, 'gstat_local_denom', d0, ...
    'xi', go.xi, 'h0', go.h0, 'G0bare', go.G0bare, ...
    'q_cancel', q_cancel, 'Q', Q);
end

function J = analytic_jacobian(u, node)
nw = numel(node.wn);
Sigma = u(1:nw);
K0 = u(end);
G0 = node.G0(:);
Jf = node.Jnu_flat(:);
D = 1+Sigma;

den = D.'+Jf*G0.';
A = mean(Jf./den, 1).';
Ap = -mean(Jf./(den.^2), 1).';
Hk = 1-A.*G0;
Nk = A.*D;
kp = ((Ap.*D+A).*Hk+Nk.*Ap.*G0)./(Hk.^2);

nvar = nw+1;
dK = zeros(nw, nvar);
dK(2:nw, 2:nw) = diag(kp(2:nw));
dK(1,end) = 1;
dK0 = zeros(1, nvar);
dK0(end) = 1;

g = node.g(:);
wts = node.wts(:);
L = zeros(3, nvar);
for p = 1:3
    L(p,:) = ((wts.*g.^p).'/node.beta)*dK;
end

tl = node.tl;
pref = tl.M2/tl.n01^2;
c0 = 1-tl.n01^2;
c1 = 0.5*(tl.g0+node.beta*c0);
am_pref = tl.m^2/tl.n01^2;
am_K = c0*(1+0.5*node.beta*tl.Delta*tl.n01)*tl.g0;
mfac = 2*tl.m^2/tl.M2;

dalpha = pref*(L(2,:)-c1*L(1,:));
dalpha_m = am_pref*(L(2,:)-tl.g0*L(1,:)+(4/tl.g0)*L(3,:)-am_K*dK0);
dgamma = pref*(repmat(L(1,:), nw, 1)-c0*dK);
dgamma0 = pref*(L(1,:)-c0*dK0);
dSigma_map = repmat(dalpha-dalpha_m, nw, 1)+ ...
    g.*(dgamma-mfac*repmat(dgamma0, nw, 1));
J_sigma = dSigma_map-[eye(nw), zeros(nw,1)];

Kdyn = A.*D./Hk;
lam = invz_lambdas([K0; Kdyn(2:end)], g, wts, node.beta, [1 2 3]);
a = node.G0inel0;
b = node.G0el0;
d0 = 1+Sigma(1)+K0*a;
dd0 = zeros(1, nvar);
dd0(1) = 1;
dd0(end) = a;
t = tl.m^2*tl.n01^2*node.beta*K0-tl.M2*node.beta*lam(1);
num_xi = 1+tanh(t);
den_xi = 1+(4*tl.n01^2*K0*tl.g0+2*lam(2)+tl.g0*lam(1))*tl.M2/tl.n01^2;
dt = tl.m^2*tl.n01^2*node.beta*dK0-tl.M2*node.beta*L(1,:);
dnum_xi = (1-tanh(t)^2)*dt;
dden_xi = (tl.M2/tl.n01^2)*(4*tl.n01^2*tl.g0*dK0+2*L(2,:)+tl.g0*L(1,:));
xi = num_xi/den_xi;
dxi = (dnum_xi*den_xi-num_xi*dden_xi)/(den_xi^2);
Hz = a+xi*b*d0;
dHz = b*(d0*dxi+xi*dd0);
z = d0/Hz;
dz = (dd0*Hz-d0*dHz)/(Hz^2);

Jscale = max(abs(Jf));
if isinf(z) || abs(z) > Jscale/sqrt(eps)
    rho = Hz/d0;
    drho = (dHz*d0-Hz*dd0)/(d0^2);
    c = Jf-K0;
    F = 1./(1+rho*c);
    Qrho = mean(F);
    Nrho = mean(Jf.*F);
    dF = -(F.^2).*(c*drho-rho*dK0);
    dQ = mean(dF, 1);
    dN = mean(Jf.*dF, 1);
    dJloc = (Qrho*dN-Nrho*dQ)/(Qrho^2);
else
    E = z+Jf-K0;
    Gbar = mean(1./E);
    N = mean(Jf./E);
    dE = dz-dK0;
    dGbar = -mean(1./(E.^2))*dE;
    dN = -mean(Jf./(E.^2))*dE;
    dJloc = (Gbar*dN-N*dGbar)/(Gbar^2);
end
J = [J_sigma; (dK0-dJloc)/Jscale];
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
