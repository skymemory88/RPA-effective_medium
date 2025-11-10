% Effctive_medium.m
% Corrected implementation of effective medium theory for tensor susceptibility
% Now supports concurrent computation over cVar (temperature/field)

close all;

%% Step 1: Load existing RPA results from MF_RPA_Yikai.m
% Expected inputs (passed from main script):
% - cVar: continuous variable (temperature or field) [1 x n_cVar]
% - freq_total: frequency range [GHz] [1 x n_freq]
% - chi0: non-interacting susceptibility [3 x 3 x n_freq x n_cVar x n_q]
% - chiq: RPA susceptibility [3 x 3 x n_freq x n_cVar x n_q]
% - qvec: q-points [n_q x 3]
% - Jq_RPA: RPA interaction [3 x 3 x n_cVar x n_q]
% - dscrt_var: discrete variable (single value)

% Extract dimensions
n_omega = size(chiq, 3);  % Number of frequencies
n_cVar = size(chiq, 4);   % Number of continuous variable points
n_q = size(chiq, 5);      % Number of q-points

fprintf('\n=== Effective Medium Theory with cVar dependence ===\n');
fprintf('Computing for %d %s points concurrently...\n', n_cVar, scanMode);
fprintf('Frequencies: %d points from %.2f to %.2f GHz\n', n_omega, min(freq_total), max(freq_total));
fprintf('Q-points: %d\n', n_q);
fprintf('----------------------------------------------------\n');

% Generate LiHoF4 parameters
params_lihof4 = struct();
params_lihof4.Delta = 0.7e-3;  % 0.7 meV singlet splitting (in eV)
params_lihof4.M = 6.551;        % Dipole matrix element
params_lihof4.m = 0.1;          % Small compared to M
params_lihof4.g_factor = 1.25;  % Ho³⁺ g-factor
params_lihof4.a = 5.175;        % Lattice constant (Angstrom)
params_lihof4.c = 10.75;        % c-axis (tetragonal)

% Extract J(q) - frequency independent interaction
J_q_RPA = squeeze(Jq_RPA(:,:,1,:));  % Dimensions: [3, 3, n_q]

%% Step 2: Initialize storage for effective medium results
% Pre-allocate arrays for all cVar points
K_emt = zeros(3, 3, n_omega, n_cVar);           % Effective medium K
G_local_emt = zeros(3, 3, n_omega, n_cVar);     % Local Green's function
chi_emt = zeros(3, 3, n_omega, n_cVar);         % Effective medium susceptibility
converged_flags = false(n_cVar, 1);             % Convergence status

%% Step 3: Parallel computation over cVar
% Setup convergence parameters (same for all cVar points)
scf_params_base = struct();
scf_params_base.max_iter = 100;
scf_params_base.tol = 1e-5;
scf_params_base.mixing_alpha = 0.2;  % Start conservative
scf_params_base.use_anderson = true;
scf_params_base.mixing_depth = 7;
scf_params_base.anderson_beta = 0.7;

fprintf('\n=== Starting parallel self-consistent calculations ===\n');

% PARFOR loop over continuous variable (temperature/field)
parfor icVar = 1:n_cVar
    % Extract data for this specific cVar point
    cVar_val = cVar(icVar);

    % Determine temperature for this point
    switch scanMode
        case 'field'
            T_local = dscrt_var;  % Fixed temperature
            beta_local = 1/(8.617e-5 * T_local);
            var_str = sprintf('B = %.3f T, T = %.3f K', cVar_val, T_local);
        case 'temp'
            T_local = cVar_val;  % Variable temperature
            beta_local = 1/(8.617e-5 * T_local);
            var_str = sprintf('T = %.3f K, B = %.3f T', T_local, dscrt_var);
        otherwise
            error('Unknown scan mode: %s', scanMode);
    end

    % Generate Matsubara frequencies for this temperature
    omega_n_local = 2*(0:n_omega-1) * pi / beta_local;

    % Extract G0 from chi0 (non-interacting susceptibility)
    G0_local = zeros(3, 3, n_omega);
    for iw = 1:n_omega
        % Average over all q points for G0
        G0_local(:,:,iw) = -mean(chi0(:,:,iw,icVar,:), 5);
    end

    % Setup local SCF parameters
    scf_params = scf_params_base;
    scf_params.beta = beta_local;
    scf_params.n_omega = n_omega;
    scf_params.n_q = n_q;
    scf_params.omega_n = omega_n_local;
    scf_params.G0 = G0_local;
    scf_params.J_q = J_q_RPA;

    % Run self-consistent calculation for this cVar point
    [K_local, G_local_local, converged] = compute_effective_medium(scf_params, var_str);

    % Store results
    K_emt(:,:,:,icVar) = K_local;
    G_local_emt(:,:,:,icVar) = G_local_local;
    converged_flags(icVar) = converged;

    % Compute susceptibility: χ = -G/β
    for iw = 1:n_omega
        chi_emt(:,:,iw,icVar) = -G_local_local(:,:,iw) / beta_local;
    end
end

% Check convergence status
n_converged = sum(converged_flags);
fprintf('\n=== Convergence Summary ===\n');
fprintf('Converged: %d / %d points (%.1f%%)\n', n_converged, n_cVar, 100*n_converged/n_cVar);
if n_converged < n_cVar
    warning('%d points did not converge!', n_cVar - n_converged);
end

%% Step 4: Create 6-panel comparison figure (RPA vs Effective Medium)
fprintf('\n=== Creating comparison plots ===\n');

% Create figure with 2x3 panels
fig_comparison = figure('Position', [100, 100, 1800, 1000]);

% Diagonal components to plot
components = [1, 2, 3];  % xx, yy, zz
comp_labels = {'xx', 'yy', 'zz'};

% Custom colormap (from MF_RPA_Yikai.m)
cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
[linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
[ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
cmap = flip(cmap,1);

% Determine axis labels based on scan mode
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

% Plot each diagonal component
for ii = 1:3
    comp_idx = components(ii);

    % Top panels: RPA susceptibility (averaged over q)
    subplot(2, 3, ii);

    % Extract diagonal component and average over q-points
    chi_rpa_comp = squeeze(chiq(comp_idx, comp_idx, :, :, :));  % [n_freq x n_cVar x n_q]
    chi_rpa_avg = mean(chi_rpa_comp, 3);  % Average over q: [n_freq x n_cVar]

    % Plot as color map (log scale of absolute value)
    data_rpa = mag2db(abs(imag(chi_rpa_avg)));
    hp = pcolor(cVar, freq_total, data_rpa);
    set(hp, 'EdgeColor', 'none');
    colormap(gca, cmap);
    caxis('auto');
    colorbar;
    ylabel('Frequency (GHz)');
    title(sprintf('RPA: Im[\\chi^{%s}] (%s)', comp_labels{ii}, title_prefix), 'Interpreter', 'tex');
    set(gca, 'FontSize', 10);

    % Bottom panels: Effective Medium susceptibility
    subplot(2, 3, ii + 3);

    % Extract diagonal component from effective medium results
    chi_emt_comp = squeeze(chi_emt(comp_idx, comp_idx, :, :));  % [n_freq x n_cVar]

    % Plot as color map (log scale of absolute value)
    data_emt = mag2db(abs(imag(chi_emt_comp)));
    hp = pcolor(cVar, freq_total, data_emt);
    set(hp, 'EdgeColor', 'none');
    colormap(gca, cmap);
    caxis('auto');
    colorbar;
    xlabel(xlab);
    ylabel('Frequency (GHz)');
    title(sprintf('Eff. Medium: Im[\\chi^{%s}]', comp_labels{ii}), 'Interpreter', 'tex');
    set(gca, 'FontSize', 10);
end

% Add overall title
sgtitle('Comparison: RPA vs Effective Medium Theory', 'FontSize', 14, 'FontWeight', 'bold');

%% Supporting Functions

%% Function 1: Self-consistent effective medium calculation for single cVar point
function [K, G_local, converged] = compute_effective_medium(scf_params, var_str)
    % Extract parameters
    n_omega = scf_params.n_omega;
    n_q = scf_params.n_q;
    G0_RPA = scf_params.G0;
    J_q_RPA = scf_params.J_q;
    beta = scf_params.beta;

    % Initialize K and G_local
    K = zeros(3, 3, n_omega);
    G_local = G0_RPA;  % Initial guess

    % For Anderson mixing history
    if scf_params.use_anderson
        K_history = cell(scf_params.mixing_depth, 1);
        F_history = cell(scf_params.mixing_depth, 1);
    end

    % Main SCF loop
    converged = false;
    for iter = 1:scf_params.max_iter
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

                % Regularize if needed
                if rcond(denom) < 1e-12
                    denom = denom + 1e-12 * eye(3);
                end

                % Solve denom * G_q = G_local
                G_q(:,:,iq,iw) = denom \ G_local_iw;
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

        % Apply mixing to G_local
        G_local_mixing = 0.5;
        G_local = (1 - G_local_mixing) * G_local + G_local_mixing * G_local_new;

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
            if rcond(G_local(:,:,iw)) > 1e-12
                K_new(:,:,iw) = sum_JG / G_local(:,:,iw);
            else
                K_new(:,:,iw) = K(:,:,iw);
            end

            % Enforce Hermiticity
            K_new(:,:,iw) = (K_new(:,:,iw) + K_new(:,:,iw)') / 2;
        end

        % Step D: Apply mixing to K
        if scf_params.use_anderson && iter > 2
            idx = mod(iter-1, scf_params.mixing_depth) + 1;
            K_history{idx} = K;
            F_history{idx} = K_new;
            K = AndersonMix(K_history, F_history, iter, scf_params);
        else
            K = (1 - scf_params.mixing_alpha) * K + scf_params.mixing_alpha * K_new;
        end

        % Step E: Compute residuals
        residual_K = 0;
        residual_G = 0;
        for iw = 1:n_omega
            residual_K = max(residual_K, max(max(abs(K(:,:,iw) - K_old(:,:,iw)))));
            residual_G = max(residual_G, max(max(abs(G_local(:,:,iw) - G_local_old(:,:,iw)))));
        end

        residual = max(residual_K, residual_G);

        % Check convergence
        if residual < scf_params.tol
            converged = true;
            fprintf('  %s: Converged in %d iterations (residual = %.2e)\n', var_str, iter, residual);
            break;
        end

        % Adaptive mixing
        if iter == scf_params.max_iter
            fprintf('  %s: Not converged after %d iterations (residual = %.2e)\n', var_str, iter, residual);
        end
    end
end

%% Function 2: Improved mixing for tensor case
function K_new = AndersonMix(K_history, F_history, iter, params)
    % Anderson acceleration adapted for 3x3 tensor mixing
    % K_history: cell array of previous K tensors (circular buffer)
    % F_history: cell array of F(K) = K_new from self-consistency
    % iter: current iteration number

    m = min(params.mixing_depth, iter-1);  % Depth of history to use

    if m < 2
        % Not enough history: simple mixing
        % Get current index in circular buffer
        idx = mod(iter-1, params.mixing_depth) + 1;
        alpha = params.mixing_alpha;
        K_new = (1-alpha) * K_history{idx} + alpha * F_history{idx};
        return;
    end

    % Get indices of stored history (handling circular buffer correctly)
    % Most recent is at idx_current, go backwards m steps
    idx_current = mod(iter-1, params.mixing_depth) + 1;

    % Build indices going backwards from current iteration
    % Only access cells that have been filled
    indices = zeros(m, 1);
    for i = 1:m
        % Map to actual iteration numbers that were stored
        actual_iter = iter - (m - i + 1);  % Go back from current
        indices(i) = mod(actual_iter-1, params.mixing_depth) + 1;
    end

    % Flatten tensors for linear algebra
    % Check that the current cell is filled
    if isempty(K_history{idx_current})
        % Fallback to simple mixing
        alpha_simple = params.mixing_alpha;
        K_new = (1-alpha_simple) * K_history{idx_current} + alpha_simple * F_history{idx_current};
        return;
    end

    [n1, n2, n3] = size(K_history{idx_current});
    n_total = n1 * n2 * n3;

    % Build residual matrices using correct indices
    % Filter out any empty cells
    R = zeros(n_total, m);
    valid_count = 0;
    for i = 1:m
        idx = indices(i);
        % Check if this cell is filled
        if ~isempty(K_history{idx}) && ~isempty(F_history{idx})
            valid_count = valid_count + 1;
            K_flat = reshape(K_history{idx}, [], 1);
            F_flat = reshape(F_history{idx}, [], 1);
            R(:,valid_count) = F_flat - K_flat;
        end
    end

    % Trim R to only valid columns
    if valid_count < m
        R = R(:, 1:valid_count);
        m = valid_count;
    end

    % If not enough valid history, fall back to simple mixing
    if m < 2
        alpha_simple = params.mixing_alpha;
        K_new = (1-alpha_simple) * K_history{idx_current} + alpha_simple * F_history{idx_current};
        return;
    end

    % Compute optimal linear combination
    % Minimize ||R * alpha||² subject to sum(alpha) = 1

    % Gram matrix
    G = R' * R;

    % Add regularization for stability
    G = G + 1e-10 * eye(m);

    % Solve for coefficients
    ones_vec = ones(m, 1);
    try
        % Solve [G, 1; 1', 0] * [alpha; lambda] = [0; 1]
        A = [G, ones_vec; ones_vec', 0];
        b = [zeros(m,1); 1];
        x = A \ b;
        alpha = x(1:m);
    catch
        % Fallback to simple mixing if singular
        alpha_simple = params.mixing_alpha;
        K_new = (1-alpha_simple) * K_history{idx_current} + alpha_simple * F_history{idx_current};
        return;
    end

    % Construct new K from linear combination using only valid indices
    K_new = zeros(size(K_history{idx_current}));
    F_avg = zeros(size(F_history{idx_current}));

    valid_idx = 0;
    for i = 1:length(indices)
        idx = indices(i);
        % Only use filled cells (same check as before)
        if ~isempty(K_history{idx}) && ~isempty(F_history{idx})
            valid_idx = valid_idx + 1;
            K_new = K_new + alpha(valid_idx) * K_history{idx};
            F_avg = F_avg + alpha(valid_idx) * F_history{idx};
        end
    end

    % Add relaxation
    beta = params.anderson_beta;  % typically 0.5-1.0
    K_new = (1-beta) * K_new + beta * F_avg;
end
