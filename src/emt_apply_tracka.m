function [G_mod, track_state] = emt_apply_tracka(G_in, K_in, params)
% EMT_APPLY_TRACKA Apply optional Track-A scalar renormalization to local G.

G_mod = G_in;
track_state = struct('enabled', false);

if ~isfield(params, 'track_a') || ~params.track_a.enabled
    return;
end

comp = params.track_a.component;
if comp < 1 || comp > size(G_in, 1)
    error('emt_apply_tracka:badComponent', ...
        'track_a.component must be between 1 and %d', size(G_in, 1));
end

n_omega = size(G_in, 3);
g_scalar = squeeze(G_in(comp, comp, :));
k_scalar = squeeze(K_in(comp, comp, :));

weights = params.track_a.frequency_weights;
if isempty(weights)
    weights = ones(n_omega, 1);
end

lambda = emt_compute_lambdas( ...
    k_scalar, g_scalar, params.track_a.beta, params.track_a.p_list, weights);
[X, aux] = emt_compute_x_from_lambdas(g_scalar, lambda, params.track_a, k_scalar);

apply_mode = lower(params.track_a.apply_mode);
if strcmp(apply_mode, 'isotropic')
    for iw = 1:n_omega
        G_mod(:,:,iw) = aux.scale(iw) * G_mod(:,:,iw);
    end
elseif strcmp(apply_mode, 'channel')
    for iw = 1:n_omega
        G_mod(comp, comp, iw) = aux.scale(iw) * G_mod(comp, comp, iw);
    end
else
    error('emt_apply_tracka:badApplyMode', ...
        'Unknown track_a.apply_mode: %s', params.track_a.apply_mode);
end

track_state.enabled = true;
track_state.component = comp;
track_state.lambda = lambda;
track_state.X = X;
track_state.a = aux.a;
track_state.gamma = aux.gamma;
track_state.scale = aux.scale;
track_state.model = aux.model;

end
