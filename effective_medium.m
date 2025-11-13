% Effctive_medium.m
% Corrected implementation of effective medium theory for tensor susceptibility
% Now supports concurrent computation over cVar (temperature/field)

close all;

%% Step 1: Load existing RPA results from MF_RPA_Yikai.m
% Expected inputs (passed from main script):
% - cVar: continuous variable (temperature or field) [1 x n_cVar]
% - freq_total: frequency range [GHz] [1 x n_freq]
% - chi0: non-interacting susceptibility [3 x 3 x n_freq x n_cVar x n_q]
% - chi_ini: initial susceptibility [3 x 3 x n_freq x n_cVar x n_q]
% - qvec: q-points [n_q x 3]
% - Jq: RPA interaction [3 x 3 x n_cVar x n_q]
% - dscrt_var: discrete variable (single value)

% Extract dimensions
% chi_ini = chiq; % use RPA suscpetibility as the initial guess
chi_ini = chi0; % use single-ion susceptibility as the initial guess
n_omega = size(chi_ini, 3);  % Number of frequencies
n_cVar = size(chi_ini, 4);   % Number of continuous variable points
n_q = size(qvec,1);      % Number of q-points

fprintf('\n=== Effective Medium Theory with cVar dependence ===\n');
fprintf('Computing for %d %s points concurrently...\n', n_cVar, scanMode);
fprintf('Frequencies: %d points from %.2f to %.2f GHz\n', n_omega, min(freq_total), max(freq_total));
fprintf('Q-points: %d\n', n_q);
fprintf('----------------------------------------------------\n');

% Jq carries J(q) for each cVar: [3 x 3 x n_cVar x n_q]

%% Step 2: Initialize storage for effective medium results
% Pre-allocate arrays for all cVar points
K_emt = zeros(3, 3, n_omega, n_cVar);           % Effective medium K
G_local_emt = zeros(3, 3, n_omega, n_cVar);     % Local Green's function
chi_emt = zeros(3, 3, n_omega, n_cVar);         % Effective medium susceptibility
converged_flags = false(n_cVar, 1);             % Convergence status

%% Step 3: Parallel computation over cVar
% Setup convergence parameters (same for all cVar points)
scf_params_base = struct();
scf_params_base.max_iter = 1e3;
scf_params_base.tol = 1e-5;
scf_params_base.mixing_alpha = 0.2;
scf_params_base.G_damp = 0.2;

fprintf('\n=== Starting parallel self-consistent calculations ===\n');

% Loop over continuous variable (temperature/field)
parfor ii = 1:n_cVar
    cVar_val = cVar(ii);
    [beta_local, var_str] = describe_state(scanMode, cVar_val, dscrt_var);
    scf_params = prepare_scf_params(scf_params_base, beta_local, n_omega, n_q, chi0, ii, Jq);

    [K_local, G_local_local, converged] = compute_effective_medium(scf_params, var_str);

    % Store results
    K_emt(:,:,:,ii) = K_local;
    G_local_emt(:,:,:,ii) = G_local_local;
    converged_flags(ii) = converged;

    % Compute susceptibility: χ = -G/β
    for iw = 1:n_omega
        chi_emt(:,:,iw,ii) = -G_local_local(:,:,iw) / beta_local;
    end
end

% Check convergence status
n_converged = sum(converged_flags);
fprintf('\n=== Convergence Summary (Initial) ===\n');
fprintf('Converged: %d / %d points (%.1f%%)\n', n_converged, n_cVar, 100*n_converged/n_cVar);

% Neighbor-seeding retry for non-converged points
if n_converged < n_cVar && n_converged > 0
    fprintf('\n=== Retrying non-converged points with neighbor seeding ===\n');
    non_converged_indices = find(~converged_flags);
    retry_success = 0;

    for jj = 1:length(non_converged_indices)
        ii = non_converged_indices(jj);
        cVar_val = cVar(ii);

        [beta_local, var_str] = describe_state(scanMode, cVar_val, dscrt_var);

        % Find nearest converged neighbor(s)
        converged_indices = find(converged_flags);
        if isempty(converged_indices)
            continue;
        end

        scf_params = prepare_scf_params(scf_params_base, beta_local, n_omega, n_q, chi0, ii, Jq);

        % Try multiple neighbor strategies
        converged_this_point = false;

        % Strategy 1: Single nearest neighbor
        [distances, sort_idx] = sort(abs(converged_indices - ii));
        neighbor_idx = converged_indices(sort_idx(1));

        fprintf('  Retry %s using nearest neighbor (idx=%d)...\n', var_str, neighbor_idx);
        K_init = K_emt(:,:,:,neighbor_idx);
        G_init = G_local_emt(:,:,:,neighbor_idx);

        [K_local, G_local_local, converged] = compute_effective_medium_seeded(...
            scf_params, var_str, K_init, G_init);

        if converged
            [converged_this_point, closure_val] = seed_acceptable(K_local, G_local_local, scf_params.J_q, ...
                'Seeded solution');
            if converged_this_point
                fprintf('    Accepted (closure %.2e)\n', closure_val);
            end
        else
            % Strategy 2: Try second-nearest neighbor if available
            if length(converged_indices) >= 2
                neighbor_idx2 = converged_indices(sort_idx(2));
                fprintf('  Retry %s using 2nd-nearest neighbor (idx=%d)...\n', var_str, neighbor_idx2);
                K_init = K_emt(:,:,:,neighbor_idx2);
                G_init = G_local_emt(:,:,:,neighbor_idx2);

                [K_local, G_local_local, converged] = compute_effective_medium_seeded(...
                    scf_params, var_str, K_init, G_init);

                if converged
                    [converged_this_point, closure_val] = seed_acceptable(K_local, G_local_local, scf_params.J_q, ...
                        'Seeded solution (2nd)');
                    if converged_this_point
                        fprintf('    2nd neighbor accepted (closure %.2e)\n', closure_val);
                    end
                end
            end

            % Strategy 3: Interpolate between two nearest neighbors
            if ~converged_this_point && length(converged_indices) >= 2
                neighbor_idx1 = converged_indices(sort_idx(1));
                neighbor_idx2 = converged_indices(sort_idx(2));

                % Weighted interpolation based on distance
                d1 = abs(neighbor_idx1 - ii);
                d2 = abs(neighbor_idx2 - ii);
                w1 = d2 / (d1 + d2);  % Closer neighbor gets more weight
                w2 = d1 / (d1 + d2);

                fprintf('  Retry %s using interpolation (idx=%d,%.2f + idx=%d,%.2f)...\n', ...
                    var_str, neighbor_idx1, w1, neighbor_idx2, w2);

                K_init = w1 * K_emt(:,:,:,neighbor_idx1) + w2 * K_emt(:,:,:,neighbor_idx2);
                G_init = w1 * G_local_emt(:,:,:,neighbor_idx1) + w2 * G_local_emt(:,:,:,neighbor_idx2);

                [K_local, G_local_local, converged] = compute_effective_medium_seeded(...
                    scf_params, var_str, K_init, G_init);

                if converged
                    [converged_this_point, closure_val] = seed_acceptable(K_local, G_local_local, scf_params.J_q, ...
                        'Seeded solution (interp)');
                    if converged_this_point
                        fprintf('    Interpolation accepted (closure %.2e)\n', closure_val);
                    end
                end
            end
        end

        % Update results if any retry path succeeded
        if converged_this_point
            K_emt(:,:,:,ii) = K_local;
            G_local_emt(:,:,:,ii) = G_local_local;
            converged_flags(ii) = true;
            retry_success = retry_success + 1;

            % Compute susceptibility
            for iw = 1:n_omega
                chi_emt(:,:,iw,ii) = -G_local_local(:,:,iw) / beta_local;
            end
        end
    end

    n_converged = sum(converged_flags);
    fprintf('  Neighbor seeding: %d additional points converged\n', retry_success);
end

fprintf('\n=== Final Convergence Summary ===\n');
fprintf('Converged: %d / %d points (%.1f%%)\n', n_converged, n_cVar, 100*n_converged/n_cVar);
if n_converged < n_cVar
    warning('%d points did not converge!', n_cVar - n_converged);
end

plot_comparison(cVar, freq_total, chi_ini, chi_emt, scanMode, dscrt_var);

%% Supporting Functions

%% Utility: compute max closure residual across frequencies for a candidate solution
function max_res = closure_residual_max(K_in, G_local_in, J_q_in)
    I3 = eye(3);
    [~,~,n_omega] = size(G_local_in);
    n_q_loc = size(J_q_in, 3);
    max_res = 0;
    for iw = 1:n_omega
        Gl = G_local_in(:,:,iw);
        K  = K_in(:,:,iw);
        sumJG = zeros(3,3);
        for iq = 1:n_q_loc
            D = I3 + (J_q_in(:,:,iq) - K) * Gl;
            % Regularize if needed
            rc = rcond(D);
            if isnan(rc) || isinf(rc) || rc < 1e-10
                D = D + 1e-8 * I3;
            end
            Gq = D \ Gl;
            sumJG = sumJG + J_q_in(:,:,iq) * Gq;
        end
        sumJG = sumJG / n_q_loc;
        denom = norm(K*Gl, 'fro') + eps;
        res = norm(sumJG - K*Gl, 'fro') / denom;
        if res > max_res
            max_res = res;
        end
    end
end

%% Function 1: Self-consistent effective medium calculation (entry point)
function [K, G_local, converged] = compute_effective_medium(scf_params, var_str)
    % Entry point - calls internal function with no initial seed
    [K, G_local, converged] = compute_effective_medium_seeded(scf_params, var_str, [], []);
end

%% Function 2: Self-consistent effective medium calculation with optional seeding
function [K, G_local, converged] = compute_effective_medium_seeded(scf_params, var_str, K_seed, G_seed)
    % Extract parameters
    n_omega = scf_params.n_omega;
    n_q = scf_params.n_q;
    G0_RPA = scf_params.G0;
    J_q_RPA = scf_params.J_q;

    % Initialize with provided seed if available, otherwise use default
    if ~isempty(K_seed) && ~isempty(G_seed)
        K_init = K_seed;
        G_local_init = G_seed;
    else
        K_init = [];
        G_local_init = [];
    end

    opts = struct('mixing_alpha', scf_params.mixing_alpha, ...
                  'G_damp', scf_params.G_damp, ...
                  'max_iter', scf_params.max_iter);

    [K, G_local, converged, final_iter, residual] = ...
        attempt_convergence(scf_params, G0_RPA, J_q_RPA, n_omega, n_q, opts, K_init, G_local_init);

    if converged
        fprintf('  %s: Converged in %d iterations (residual = %.2e)\n', var_str, final_iter, residual);
    else
        fprintf('  %s: Not converged after %d iterations (residual = %.2e)\n', var_str, final_iter, residual);
    end
end

%% Attempt convergence with adaptive step
function [K, G_local, converged, final_iter, final_residual] = ...
    attempt_convergence(scf_params, G0_RPA, J_q_RPA, n_omega, n_q, opts, K_init, G_local_init)

    % Initialize K and G_local (use provided initial guess if available)
    if isempty(K_init) || isempty(G_local_init)
        K = zeros(3, 3, n_omega);
        G_local = G0_RPA;  % Default initial guess
    else
        K = K_init;
        G_local = G_local_init;
    end

    % Adaptive mixing parameters
    mixing_alpha = opts.mixing_alpha;
    G_damp = opts.G_damp;

    residual_history = zeros(opts.max_iter, 1);

    % Main SCF loop
    converged = false;
    final_iter = opts.max_iter;
    final_residual = inf;

    % Suppress warnings for nearly singular matrices (we handle them explicitly)
    warning_state = warning('query', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning_state2 = warning('query', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:singularMatrix');

    for iter = 1:opts.max_iter
        K_old = K;
        G_local_old = G_local;

        % Step A: Compute G(q,iω) for all q using CURRENT G_local
        G_q = zeros(3, 3, n_q, n_omega);

        for iq = 1:n_q
            for iw = 1:n_omega
                G_local_iw = G_local(:,:,iw);
                J_q_iq = J_q_RPA(:,:,iq);
                K_iw = K(:,:,iw);

                % Eq. 2.12: G(q,iω) = [1 + (J(q) - K(iω))G_local(iω)]^(-1) * G_local(iω)
                denom = eye(3) + (J_q_iq - K_iw) * G_local_iw;

                % Check for NaN/Inf in matrix before computing condition number
                has_nan_inf = any(isnan(denom(:))) || any(isinf(denom(:)));

                % Adaptive regularization based on condition number
                rc = rcond(denom);
                if has_nan_inf || isnan(rc) || isinf(rc) || rc < 1e-10
                    denom = denom + 1e-5 * eye(3);
                end

                % Solve denom * G_q = G_local using robust method
                G_q(:,:,iq,iw) = denom \ G_local_iw;

                % Verify result doesn't contain NaN/Inf
                G_temp = G_q(:,:,iq,iw);
                if any(isnan(G_temp(:))) || any(isinf(G_temp(:)))
                    denom = eye(3) + (J_q_iq - K_iw) * G_local_iw + 1e-4 * eye(3);
                    G_q(:,:,iq,iw) = denom \ G_local_iw;
                end
            end
        end

        % Step B: Update local Green's function (q-averaged)
        G_local_new = zeros(3, 3, n_omega);
        for iw = 1:n_omega
            for iq = 1:n_q
                G_local_new(:,:,iw) = G_local_new(:,:,iw) + G_q(:,:,iq,iw);
            end
            G_local_new(:,:,iw) = G_local_new(:,:,iw) / n_q;
        end

        % Apply adaptive damped mixing to G_local
        % Note: G_damp acts as (1-alpha) where larger G_damp means more old value
        G_local_mixed = G_damp * G_local + (1 - G_damp) * G_local_new;

        % Safety check: if mixed result contains NaN/Inf, reduce mixing and retry
        if any(isnan(G_local_mixed(:))) || any(isinf(G_local_mixed(:)))
            G_local = G_local_old;
        else
            G_local = G_local_mixed;
        end

        % Step C: New K from self-consistency equation (2.11)
        K_new = zeros(3, 3, n_omega);

        for iw = 1:n_omega
            % Sum: Σ_q J(q) * G(q,iω)
            sum_JG = zeros(3, 3);
            for iq = 1:n_q
                sum_JG = sum_JG + J_q_RPA(:,:,iq) * G_q(:,:,iq,iw);
            end
            sum_JG = sum_JG / n_q;

            % K = sum_JG * inv(G_local)
            G_local_iw = G_local(:,:,iw);
            if any(isnan(G_local_iw(:))) || any(isinf(G_local_iw(:))) || ...
                    any(isnan(sum_JG(:))) || any(isinf(sum_JG(:)))
                K_new(:,:,iw) = K(:,:,iw);
                continue;
            end

            rc_G = rcond(G_local_iw);
            if ~isnan(rc_G) && ~isinf(rc_G) && rc_G > 1e-12
                K_tmp = sum_JG / G_local_iw;
                K_tmp = (K_tmp + K_tmp') / 2;
                if any(isnan(K_tmp(:))) || any(isinf(K_tmp(:)))
                    K_new(:,:,iw) = K(:,:,iw);
                else
                    K_new(:,:,iw) = K_tmp;
                end
            else
                K_new(:,:,iw) = K(:,:,iw);
            end
        end

        % Step D: Apply simple mixing to K
        K_update = K_new - K;
        K = K + mixing_alpha * K_update;
        if any(isnan(K(:))) || any(isinf(K(:)))
            K = K_old;
        end

        % Step E: Compute residuals
        residual_K = 0;
        residual_G = 0;
        for iw = 1:n_omega
            residual_K = max(residual_K, max(max(abs(K(:,:,iw) - K_old(:,:,iw)))));
            residual_G = max(residual_G, max(max(abs(G_local(:,:,iw) - G_local_old(:,:,iw)))));
        end

        residual = max(residual_K, residual_G);
        residual_history(iter) = residual;
        final_residual = residual;

        % Option B: residual-based backoff to prevent overshoot
        if iter > 1 && residual > 1.05 * residual_history(iter-1)
            mixing_alpha = max(mixing_alpha * 0.5, 1e-3);
            G_damp = min(G_damp * 1.1, 0.98);
        end

        % Check convergence
        if residual < scf_params.tol
            converged = true;
            final_iter = iter;
            break;
        end
    end

    % Restore warning states
    warning(warning_state.state, 'MATLAB:nearlySingularMatrix');
    warning(warning_state2.state, 'MATLAB:singularMatrix');
end

%% Helper: build SCF parameter struct for a given cVar index
function scf_params = prepare_scf_params(base_params, beta_local, n_omega, n_q, chi0, idx, Jq)
    scf_params = base_params;
    scf_params.beta = beta_local;
    scf_params.n_omega = n_omega;
    scf_params.n_q = n_q;
    scf_params.omega_n = bosonic_frequencies(n_omega, beta_local);
    scf_params.G0 = extract_G0(chi0, idx);
    scf_params.J_q = Jq;
end

function omega = bosonic_frequencies(n_omega, beta_local)
    omega = 2*(0:n_omega-1) * pi / beta_local;
end

function G0 = extract_G0(chi0, idx)
    G0 = squeeze(-mean(chi0(:,:, :, idx, :), 5));
end

function [beta_local, descriptor] = describe_state(scanMode, cVar_val, dscrt_var)
    kB = 8.617e-5; % eV/K
    switch scanMode
        case 'field'
            T_local = dscrt_var;
            beta_local = 1 / (kB * T_local);
            descriptor = sprintf('B = %.3f T, T = %.3f K', cVar_val, T_local);
        case 'temp'
            T_local = cVar_val;
            beta_local = 1 / (kB * T_local);
            descriptor = sprintf('T = %.3f K, B = %.3f T', T_local, dscrt_var);
        otherwise
            error('Unknown scan mode: %s', scanMode);
    end
end

function [ok, max_closure] = seed_acceptable(K_candidate, G_candidate, J_q_slice, label)
    max_closure = closure_residual_max(K_candidate, G_candidate, J_q_slice);
    ok = max_closure < 1e-3;
    if ~ok
        fprintf('    %s rejected: closure residual %.2e > 1e-3\n', label, max_closure);
    end
end

function plot_comparison(cVar, freq_total, chi_ini, chi_emt, scanMode, dscrt_var)
    fprintf('\n=== Creating comparison plots ===\n');

    figure('Position', [100, 100, 1800, 1000]);
    components = [1, 2, 3];
    comp_labels = {'xx', 'yy', 'zz'};

    cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
        [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
        [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
    cmap = flip(cmap,1);

    switch scanMode
        case 'field'
            xlab = 'Magnetic Field (T)';
            title_prefix = sprintf('T = %.3f K', dscrt_var);
        case 'temp'
            xlab = 'Temperature (K)';
            title_prefix = sprintf('B = %.3f T', dscrt_var);
        otherwise
            xlab = 'cVar';
            title_prefix = '';
    end

    for jj = 1:3
        comp_idx = components(jj);

        subplot(2, 3, jj);
        chi_rpa_comp = squeeze(chi_ini(comp_idx, comp_idx, :, :, :));
        chi_rpa_avg = mean(chi_rpa_comp, 3);
        data_rpa = mag2db(abs(imag(chi_rpa_avg)));
        hp = pcolor(cVar, freq_total, data_rpa);
        set(hp, 'EdgeColor', 'none');
        colormap(gca, cmap);
        colorbar;
        ylabel('Frequency (GHz)');
        title(sprintf('RPA: Im[\\chi^{%s}] (%s)', comp_labels{jj}, title_prefix), 'Interpreter', 'tex');
        set(gca, 'FontSize', 10);

        subplot(2, 3, jj + 3);
        chi_emt_comp = squeeze(chi_emt(comp_idx, comp_idx, :, :));
        data_emt = mag2db(abs(imag(chi_emt_comp)));
        hp = pcolor(cVar, freq_total, data_emt);
        set(hp, 'EdgeColor', 'none');
        colormap(gca, cmap);
        colorbar;
        xlabel(xlab);
        ylabel('Frequency (GHz)');
        title(sprintf('Eff. Medium: Im[\\chi^{%s}]', comp_labels{jj}), 'Interpreter', 'tex');
        set(gca, 'FontSize', 10);
    end

    sgtitle('Comparison: RPA vs Effective Medium Theory', 'FontSize', 14, 'FontWeight', 'bold');
end
