function point = emt_solve_point(seed_local, J_q, params, state_label, seed_guess)
% EMT_SOLVE_POINT Solve EMT SCF for one field/temperature point.

if nargin < 5
    seed_guess = struct();
end

if isfield(seed_guess, 'G') && ~isempty(seed_guess.G)
    G_local = seed_guess.G;
else
    G_local = seed_local;
end

if isfield(seed_guess, 'K') && ~isempty(seed_guess.K)
    K = seed_guess.K;
else
    K = emt_initial_k_from_seed(G_local, J_q, params);
end

alpha = params.mix_alpha;
residual_history = inf(params.max_iter, 1);
track_history = cell(params.max_iter, 1);

converged = false;
iter_count = params.max_iter;
residual = inf;
G_q_last = [];

for iter = 1:params.max_iter
    K_old = K;
    G_old = G_local;

    [G_q, G_avg, track_state] = emt_compute_gq_and_average(G_old, K_old, J_q, params);
    G_local = (1 - alpha) * G_old + alpha * G_avg;

    K_new = emt_update_medium_from_gq(G_q, G_local, J_q, K_old, params);
    K = (1 - alpha) * K_old + alpha * K_new;

    residual = max(max(abs(K(:) - K_old(:))), max(abs(G_local(:) - G_old(:))));
    residual_history(iter) = residual;
    track_history{iter} = track_state;
    G_q_last = G_q;

    if params.verbose && (iter == 1 || mod(iter, 50) == 0)
        fprintf('    %s | iter=%d residual=%.2e alpha=%.3f\n', ...
            state_label, iter, residual, alpha);
    end

    if iter > 1 && residual > params.backoff_trigger * residual_history(iter - 1)
        alpha = max(alpha * params.backoff_factor, params.min_mix_alpha);
    end

    if residual < params.tol
        converged = true;
        iter_count = iter;
        break;
    end
end

closure = emt_closure_residual(K, G_local, J_q, params);

point = struct();
point.K = K;
point.G_local = G_local;
point.chi = -G_local;
point.converged = converged;
point.iter = iter_count;
point.residual = residual;
point.closure = closure;
if params.store_gq
    point.G_q = G_q_last;
else
    point.G_q = [];
end
point.track_history = track_history(1:max(iter_count, 1));

end
