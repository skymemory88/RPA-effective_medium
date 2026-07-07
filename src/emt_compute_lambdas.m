function lambda = emt_compute_lambdas(k_scalar, g_scalar, beta, p_list, weights)
% EMT_COMPUTE_LAMBDAS Compute lambda_p moments for Track-A renormalization.
%
% lambda_p = (1 / beta) * sum_n w_n * K_n * (g_n ^ p)

if nargin < 5 || isempty(weights)
    weights = ones(size(g_scalar));
end

k_scalar = k_scalar(:);
g_scalar = g_scalar(:);
weights = weights(:);

if numel(weights) ~= numel(g_scalar)
    error('emt_compute_lambdas:badWeights', ...
        'weights and g_scalar must have same length.');
end

lambda = zeros(size(p_list));
for ip = 1:numel(p_list)
    p = p_list(ip);
    lambda(ip) = (1 / beta) * sum(weights .* k_scalar .* (g_scalar .^ p));
end

end
