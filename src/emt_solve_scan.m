function results = emt_solve_scan(inputs, params)
% EMT_SOLVE_SCAN Solve EMT over all scan points.

n_cvar = inputs.n_cVar;
n_omega = inputs.n_omega;
sample = inputs.chi_seed(1,1,1,1,1);

K_emt = zeros(3, 3, n_omega, n_cvar, 'like', sample);
G_local_emt = zeros(3, 3, n_omega, n_cvar, 'like', sample);
chi_emt = zeros(3, 3, n_omega, n_cvar, 'like', sample);
converged = false(n_cvar, 1);
iter = zeros(n_cvar, 1);
residual = inf(n_cvar, 1);
closure = inf(n_cvar, 1);

if params.store_gq
    G_q_emt = zeros(3, 3, inputs.n_q, n_omega, n_cvar, 'like', sample);
else
    G_q_emt = [];
end

fprintf('\n=== EMT scan solve (%d points) ===\n', n_cvar);

parfor ii = 1:n_cvar
    seed_slice = squeeze(inputs.chi_seed(:,:,:,ii,:));
    if ndims(seed_slice) == 4
        seed_local = mean(seed_slice, 4);
    else
        seed_local = seed_slice;
    end
    seed_local = -reshape(seed_local, 3, 3, []);

    J_q = emt_select_j_slice(inputs.Jq, ii, n_cvar);
    label = emt_state_label(inputs.scan_mode, inputs.cVar(ii), inputs.dscrt_value);

    point = emt_solve_point(seed_local, J_q, params, label);

    K_emt(:,:,:,ii) = point.K;
    G_local_emt(:,:,:,ii) = point.G_local;
    chi_emt(:,:,:,ii) = point.chi;
    converged(ii) = point.converged;
    iter(ii) = point.iter;
    residual(ii) = point.residual;
    closure(ii) = point.closure;

    if params.store_gq && ~isempty(point.G_q)
        G_q_emt(:,:,:,:,ii) = point.G_q;
    end
end

% Neighbor-seeded retry for non-converged points.
if params.use_neighbor_seed
    converged_idx = find(converged);
    if ~isempty(converged_idx) && numel(converged_idx) < n_cvar
        retry_list = find(~converged);
        fprintf('Retrying %d non-converged points with neighbor seeding...\n', ...
            numel(retry_list));

        for jj = 1:numel(retry_list)
            ii = retry_list(jj);

            seed_slice = squeeze(inputs.chi_seed(:,:,:,ii,:));
            if ndims(seed_slice) == 4
                seed_local = mean(seed_slice, 4);
            else
                seed_local = seed_slice;
            end
            seed_local = -reshape(seed_local, 3, 3, []);

            J_q = emt_select_j_slice(inputs.Jq, ii, n_cvar);
            label = emt_state_label(inputs.scan_mode, inputs.cVar(ii), inputs.dscrt_value);

            [~, order] = sort(abs(converged_idx - ii));
            idx1 = converged_idx(order(1));
            seed_guess = struct('K', K_emt(:,:,:,idx1), 'G', G_local_emt(:,:,:,idx1));
            point_try = emt_solve_point(seed_local, J_q, params, label, seed_guess);

            accept = point_try.converged && point_try.closure < params.closure_tol;

            if ~accept && numel(converged_idx) >= 2
                idx2 = converged_idx(order(2));
                d1 = abs(ii - idx1);
                d2 = abs(ii - idx2);
                w1 = d2 / max(d1 + d2, eps);
                w2 = d1 / max(d1 + d2, eps);

                seed_guess.K = w1 * K_emt(:,:,:,idx1) + w2 * K_emt(:,:,:,idx2);
                seed_guess.G = w1 * G_local_emt(:,:,:,idx1) + w2 * G_local_emt(:,:,:,idx2);
                point_try2 = emt_solve_point(seed_local, J_q, params, label, seed_guess);

                accept2 = point_try2.converged && point_try2.closure < params.closure_tol;
                if accept2 && (~accept || point_try2.closure < point_try.closure)
                    point_try = point_try2;
                    accept = true;
                end
            end

            if accept
                K_emt(:,:,:,ii) = point_try.K;
                G_local_emt(:,:,:,ii) = point_try.G_local;
                chi_emt(:,:,:,ii) = point_try.chi;
                converged(ii) = point_try.converged;
                iter(ii) = point_try.iter;
                residual(ii) = point_try.residual;
                closure(ii) = point_try.closure;
                if params.store_gq && ~isempty(point_try.G_q)
                    G_q_emt(:,:,:,:,ii) = point_try.G_q;
                end
                converged_idx = find(converged);
            end
        end
    end
end

results = struct();
results.K_emt = K_emt;
results.G_local_emt = G_local_emt;
results.chi_emt = chi_emt;
results.G_q_emt = G_q_emt;
results.converged = converged;

info = struct();
info.iter = iter;
info.residual = residual;
info.closure = closure;
info.scan_mode = inputs.scan_mode;
info.cVar = inputs.cVar(:);
info.freq = inputs.freq(:);
info.qvec = inputs.qvec;
info.params = params;
info.n_converged = sum(converged);
info.n_total = n_cvar;
results.info = info;

end
