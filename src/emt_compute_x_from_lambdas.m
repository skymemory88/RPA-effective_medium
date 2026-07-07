function [X, aux] = emt_compute_x_from_lambdas(g_scalar, lambda, track_params, k_scalar)
% EMT_COMPUTE_X_FROM_LAMBDAS Compute scalar X(iw) from lambda moments.
%
% Supports two models:
% 1) moment_ratio (default):
%      X = a + gamma * g
%      a = Re(lambda_1)
%      gamma = (Re(lambda_2) - Re(lambda_1)^2) / (g0 + eps)
% 2) linear_lambda:
%      a = c0 + sum_k c_k * Re(lambda_k)
%      gamma = d0 + sum_k d_k * Re(lambda_k)
%      X = a + gamma * g
% 3) jensen_216_219:
%      X(iw) = alpha + gamma(iw)*g(iw), with
%      alpha = (M2/n01^2) * [lambda2 - 0.5*(g0 + beta*(1-n01^2))*lambda1]
%      gamma(iw) = (M2/n01^2) * [lambda1 - (1-n01^2)K(iw)]
%      lambda_p = (1/beta) * sum_{n'} K(iwn') [g(iwn')]^p
% 4) custom_x_function:
%      user callback handles X construction exactly as desired.

if nargin < 3
    track_params = struct();
end
if nargin < 4
    k_scalar = [];
end

if ~isfield(track_params, 'model')
    track_params.model = 'moment_ratio';
end
if ~isfield(track_params, 'clamp_scale')
    track_params.clamp_scale = 10.0;
end
if ~isfield(track_params, 'linear_a')
    track_params.linear_a = [0 1 0 0 0];
end
if ~isfield(track_params, 'linear_gamma')
    track_params.linear_gamma = [0 0 1 0 0];
end
if ~isfield(track_params, 'custom_x_function')
    track_params.custom_x_function = [];
end
if ~isfield(track_params, 'M2')
    track_params.M2 = 1.0;
end
if ~isfield(track_params, 'beta')
    track_params.beta = 1.0;
end
if ~isfield(track_params, 'n01')
    track_params.n01 = 1.0;
end
if ~isfield(track_params, 'n0')
    track_params.n0 = [];
end
if ~isfield(track_params, 'n1')
    track_params.n1 = [];
end
g_scalar = g_scalar(:);
lambda_real = real(lambda(:)).';
g0 = g_scalar(1);
k_scalar = k_scalar(:);

n01 = track_params.n01;
if isempty(n01) || ~isfinite(n01)
    n01 = 1.0;
end
if (~isempty(track_params.n0) && ~isempty(track_params.n1))
    n01 = track_params.n0 - track_params.n1;
end

if ~isempty(track_params.custom_x_function)
    [X, custom_aux] = track_params.custom_x_function(g_scalar, lambda, track_params, k_scalar);
    X = X(:);
    if ~isfield(custom_aux, 'a')
        custom_aux.a = NaN;
    end
    if ~isfield(custom_aux, 'gamma')
        custom_aux.gamma = NaN;
    end
    model_name = 'custom';
else
    switch lower(track_params.model)
        case 'moment_ratio'
            l1 = get_lambda(lambda_real, 1);
            l2 = get_lambda(lambda_real, 2);
            a = l1;
            gamma = (l2 - l1^2) / (g0 + eps);
            X = a + gamma * g_scalar;
            custom_aux = struct('a', a, 'gamma', gamma);
            model_name = track_params.model;

        case 'linear_lambda'
            a_coeff = pad_coeff(track_params.linear_a, numel(lambda_real) + 1);
            g_coeff = pad_coeff(track_params.linear_gamma, numel(lambda_real) + 1);

            basis = [1, lambda_real];
            a = a_coeff * basis.';
            gamma = g_coeff * basis.';
            X = a + gamma * g_scalar;
            custom_aux = struct('a', a, 'gamma', gamma);
            model_name = track_params.model;

        case 'jensen_216_219'
            if isempty(k_scalar) || numel(k_scalar) ~= numel(g_scalar)
                error('emt_compute_x_from_lambdas:missingK', ...
                    'jensen_216_219 model requires k_scalar with same length as g_scalar.');
            end

            if abs(n01) < eps
                error('emt_compute_x_from_lambdas:badN01', ...
                    'n01 is too small for stable Jensen Track-A evaluation.');
            end

            lambda1 = get_lambda(lambda_real, 1);
            lambda2 = get_lambda(lambda_real, 2);

            alpha = (track_params.M2 / (n01^2)) * ...
                (lambda2 - 0.5 * (g0 + track_params.beta * (1 - n01^2)) * lambda1);

            gamma_vec = (track_params.M2 / (n01^2)) * ...
                (lambda1 - (1 - n01^2) * k_scalar);

            X = alpha + gamma_vec .* g_scalar;

            custom_aux = struct();
            custom_aux.a = alpha;
            custom_aux.gamma = gamma_vec;
            custom_aux.lambda1 = lambda1;
            custom_aux.lambda2 = lambda2;
            model_name = track_params.model;

        otherwise
            error('emt_compute_x_from_lambdas:badModel', ...
                'Unknown Track-A model: %s', track_params.model);
    end
end

% Clamp multiplicative scale to avoid denominator blow-up.
scale = 1 ./ (1 + X);
abs_scale = abs(scale);
mask = abs_scale > track_params.clamp_scale;
if any(mask)
    scale(mask) = track_params.clamp_scale * exp(1i * angle(scale(mask)));
    X(mask) = 1 ./ scale(mask) - 1;
end

aux = struct();
aux.lambda = lambda(:).';
aux.a = custom_aux.a;
aux.gamma = custom_aux.gamma;
if isfield(custom_aux, 'lambda1')
    aux.lambda1 = custom_aux.lambda1;
end
if isfield(custom_aux, 'lambda2')
    aux.lambda2 = custom_aux.lambda2;
end
aux.scale = scale;
aux.model = model_name;

end

function value = get_lambda(lambda_real, idx)
    if numel(lambda_real) >= idx
        value = lambda_real(idx);
    else
        value = 0;
    end
end

function coeff = pad_coeff(coeff_in, n)
    coeff = zeros(1, n);
    m = min(numel(coeff_in), n);
    coeff(1:m) = coeff_in(1:m);
end
