% Effctive_medium.m
% Corrected implementation of effective medium theory for tensor susceptibility
% Key fix: Using self-consistent G_local in equation 2.12 instead of G0

clear; close all;

%% Step 1: Load your existing RPA results
% Assuming you have these from your RPA calculation:
% - G0_RPA: [3 x 3 x n_omega] non-interacting Green's function
% - J_q_RPA: [3 x 3 x n_q] dipole-dipole interaction in q-space  
% - qpoints: [3 x n_q] your q-point mesh
% - omega_n: Matsubara frequencies

% Example loading (replace with your actual data)
% load('my_RPA_results.mat', 'G0_RPA', 'J_q_RPA', 'qpoints', 'omega_n');

% For testing, let's create realistic dummy data
% REPLACE THIS SECTION WITH YOUR ACTUAL DATA
n_omega = 100;
n_q = 8^3;
T = 1.0;  % K
beta = 1/T;
omega_n = (2*(0:n_omega-1) + 1) * pi / beta;

% Generate realistic LiHoF4 parameters
params_lihof4 = struct();
params_lihof4.Delta = 0.7e-3;  % 0.7 meV singlet splitting (in eV)
params_lihof4.M = 6.551;        % Dipole matrix element
params_lihof4.m = 0.1;          % Small compared to M
params_lihof4.g_factor = 1.25;  % Ho³⁺ g-factor
params_lihof4.a = 5.175;        % Lattice constant (Angstrom)
params_lihof4.c = 10.75;        % c-axis (tetragonal)
params_lihof4.beta = beta;

%% Step 2: Generate G0 for LiHoF4 singlet-singlet system
fprintf('Generating G0 for LiHoF4 singlet-singlet system...\n');
G0_RPA = generate_G0_singlet(omega_n, params_lihof4);

%% Step 3: Generate dipole-dipole interaction J(q)
fprintf('Computing dipole-dipole interaction tensor...\n');

% Generate q-points (replace with your actual mesh)
% Simple cubic mesh for testing
qx = linspace(-pi, pi, 8);
qy = linspace(-pi, pi, 8);
qz = linspace(-pi, pi, 8);
[QX, QY, QZ] = meshgrid(qx, qy, qz);
qpoints = [QX(:)'; QY(:)'; QZ(:)'];

params_lihof4.use_ewald = false;  % Set true for accurate long-range
J_q_RPA = compute_dipolar_tensor(qpoints, params_lihof4);

%% Step 4: Setup self-consistent calculation parameters
scf_params = struct();
scf_params.beta = beta;
scf_params.n_omega = n_omega;
scf_params.n_q = n_q;
scf_params.omega_n = omega_n;
scf_params.qpoints = qpoints;
scf_params.G0 = G0_RPA;
scf_params.J_q = J_q_RPA;

% Convergence parameters
scf_params.max_iter = 100;
scf_params.tol = 1e-5;
scf_params.mixing_alpha = 0.2;  % Start conservative
scf_params.use_anderson = true;
scf_params.mixing_depth = 7;
scf_params.anderson_beta = 0.7;

%% Step 5: Run self-consistent calculation
fprintf('\n=== Starting self-consistent K(iω_n) calculation ===\n');
fprintf('System: LiHoF4 at T = %.2f K\n', T);
fprintf('Mesh: %d q-points, %d Matsubara frequencies\n', n_q, n_omega);
fprintf('------------------------------------------------\n');

% Initialize K and G_local
K = zeros(3, 3, n_omega);

% CRITICAL INITIALIZATION: Start with G_local = G0 for first iteration
G_local = G0_RPA;  % This is the key initialization!

% For Anderson mixing history
if scf_params.use_anderson
    K_history = cell(scf_params.mixing_depth, 1);
    F_history = cell(scf_params.mixing_depth, 1);
end

% Convergence history
history = struct();
history.residual = [];
history.energy = [];
history.G_local_change = [];  % Track G_local convergence

% Main SCF loop
for iter = 1:scf_params.max_iter
    K_old = K;
    G_local_old = G_local;  % Store previous G_local for convergence check
    
    % Step A: Compute G(q,iω) for all q using CURRENT G_local
    G_q = zeros(3, 3, n_q, n_omega);
    
    for iq = 1:n_q
        for iw = 1:n_omega
            % CORRECTED: Use G_local, not G0!
            G_local_iw = G_local(:,:,iw);  % Self-consistent local Green's function
            J_q_iq = J_q_RPA(:,:,iq);
            K_iw = K(:,:,iw);
            
            % Eq. 2.12 in matrix form with correct G_local
            denom = eye(3) + (J_q_iq - K_iw) * G_local_iw;
            
            % Regularize if needed
            if rcond(denom) < 1e-12
                denom = denom + 1e-12 * eye(3);
            end
            
            G_q(:,:,iq,iw) = G_local_iw / denom;
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
    
    % Apply mixing to G_local for stability
    G_local_mixing = 0.5;  % Can adjust this separately from K mixing
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
        
        % Divide by G_local (equation 2.11)
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
        % Store in history
        idx = mod(iter-1, scf_params.mixing_depth) + 1;
        K_history{idx} = K;
        F_history{idx} = K_new;
        
        % Anderson mixing
        K = tensor_anderson_mixing(K_history, F_history, iter, scf_params);
    else
        % Simple linear mixing
        K = (1 - scf_params.mixing_alpha) * K + ...
            scf_params.mixing_alpha * K_new;
    end
    
    % Step E: Compute residuals and energy
    residual_K = 0;
    residual_G = 0;
    for iw = 1:n_omega
        residual_K = max(residual_K, max(max(abs(K(:,:,iw) - K_old(:,:,iw)))));
        residual_G = max(residual_G, max(max(abs(G_local(:,:,iw) - G_local_old(:,:,iw)))));
    end
    
    % Use combined residual for convergence
    residual = max(residual_K, residual_G);
    
    % Energy (tr[K*G])
    energy = 0;
    for iw = 1:n_omega
        energy = energy + real(trace(K(:,:,iw) * G_local(:,:,iw)));
    end
    energy = -energy / beta;
    
    % Store history
    history.residual(iter) = residual;
    history.energy(iter) = energy;
    history.G_local_change(iter) = residual_G;
    
    % Print progress
    if mod(iter, 10) == 1 || residual < scf_params.tol
        fprintf('Iter %3d: residual_K = %.3e, residual_G = %.3e, E = %.6f\n', ...
                iter, residual_K, residual_G, energy);
        fprintf('         K₁₁(ω₁) = %.4f, K₃₃(ω₁) = %.4f, ', ...
                real(K(1,1,1)), real(K(3,3,1)));
        fprintf('G₁₁(ω₁) = %.4f\n', real(G_local(1,1,1)));
    end
    
    % Check convergence
    if residual < scf_params.tol
        fprintf('\n✓ Converged in %d iterations!\n', iter);
        break;
    end
    
    % Adaptive mixing: reduce alpha if oscillating
    if iter > 10 && history.residual(end) > history.residual(end-5)
        scf_params.mixing_alpha = scf_params.mixing_alpha * 0.8;
        fprintf('  Reducing K mixing to α = %.3f\n', scf_params.mixing_alpha);
    end
end

if iter == scf_params.max_iter
    warning('Not fully converged after %d iterations', scf_params.max_iter);
end

%% Step 6: Compare with RPA
fprintf('\n=== Comparing with RPA ===\n');

% Compute RPA Green's function for comparison
G_RPA_local = zeros(3, 3, n_omega);
for iw = 1:n_omega
    G0_iw = G0_RPA(:,:,iw);
    
    % RPA: sum over q
    for iq = 1:n_q
        J_q_iq = J_q_RPA(:,:,iq);
        G_RPA_q = G0_iw / (eye(3) + J_q_iq * G0_iw);
        G_RPA_local(:,:,iw) = G_RPA_local(:,:,iw) + G_RPA_q;
    end
    G_RPA_local(:,:,iw) = G_RPA_local(:,:,iw) / n_q;
end

%% Step 7: Calculate physical observables
fprintf('\n=== Physical Observables ===\n');

% Compute observables
params_lihof4.J_q0 = J_q_RPA(:,:,1);  % J(q=0)
obs_1z = compute_observables(K, G_local, omega_n, params_lihof4);
obs_RPA = compute_observables(zeros(size(K)), G_RPA_local, omega_n, params_lihof4);

fprintf('\nComparison (1/z vs RPA):\n');
fprintf('  χ_zz/χ_xx: %.2f vs %.2f\n', obs_1z.anisotropy_ratio, obs_RPA.anisotropy_ratio);
fprintf('  Energy shift: %.4f vs %.4f meV\n', obs_1z.energy_shift*1000, obs_RPA.energy_shift*1000);

%% Step 8: Visualize comparison
figure('Position', [100, 100, 1400, 600]);

% Plot diagonal Green's functions
subplot(1,3,1);
n_plot = min(50, n_omega);
plot(omega_n(1:n_plot), real(squeeze(G_local(1,1,1:n_plot))), 'b-', 'LineWidth', 2);
hold on;
plot(omega_n(1:n_plot), real(squeeze(G_local(3,3,1:n_plot))), 'r-', 'LineWidth', 2);
plot(omega_n(1:n_plot), real(squeeze(G_RPA_local(1,1,1:n_plot))), 'b--', 'LineWidth', 1);
plot(omega_n(1:n_plot), real(squeeze(G_RPA_local(3,3,1:n_plot))), 'r--', 'LineWidth', 1);
xlabel('iω_n');
ylabel('Re[G(iω_n)]');
title('Local Green''s Function');
legend('G_{xx} (1/z)', 'G_{zz} (1/z)', 'G_{xx} (RPA)', 'G_{zz} (RPA)');
grid on;

% Plot effective medium K
subplot(1,3,2);
plot(omega_n(1:n_plot), real(squeeze(K(1,1,1:n_plot))), 'b-', 'LineWidth', 2);
hold on;
plot(omega_n(1:n_plot), real(squeeze(K(3,3,1:n_plot))), 'r-', 'LineWidth', 2);
xlabel('iω_n');
ylabel('Re[K(iω_n)]');
title('Effective Medium K(iω_n)');
legend('K_{xx}', 'K_{zz}');
grid on;

% Plot convergence
subplot(1,3,3);
semilogy(history.residual, 'k-', 'LineWidth', 2);
hold on;
semilogy(history.G_local_change, 'b-', 'LineWidth', 1);
xlabel('Iteration');
ylabel('Residual');
title('Convergence History');
legend('Total residual', 'G_{local} change');
grid on;

%% Supporting Functions

%% Function 1: Generate G0 for singlet-singlet system
function G0 = generate_G0_singlet(omega_n, params)
    % Generate non-interacting Green's function for singlet-singlet system
    % Based on equation 2.7-2.8 in the paper
    
    n_omega = length(omega_n);
    G0 = zeros(3, 3, n_omega);
    
    % Parameters
    Delta = params.Delta;  % Singlet splitting
    M = params.M;          % Dipole matrix element
    m = params.m;          % Small magnetic moment within singlets
    beta = params.beta;
    
    % Population factors
    n0 = 1 / (1 + exp(-beta * Delta));
    n1 = 1 - n0;
    n01 = n0 - n1;
    
    % Response functions (from eq. 2.8)
    for iw = 1:n_omega
        g_iw = 2 * n01 * Delta / (omega_n(iw)^2 + Delta^2);
        h_iw = beta * (1 - n01) * (iw == 1);  % Delta function at ω=0
        
        % Build G0 tensor
        % For singlet-singlet: only x-component has off-diagonal matrix elements
        % z-component may have diagonal elements
        G0(1,1,iw) = M^2 * g_iw - m^2 * h_iw;  % xx component
        G0(2,2,iw) = M^2 * g_iw - m^2 * h_iw;  % yy component (same as xx for tetragonal)
        G0(3,3,iw) = m^2 * g_iw;                % zz component (smaller)
        
        % Off-diagonal elements (if any)
        % Usually zero for singlet-singlet in simple cases
    end
end

%% Function 2: Compute dipole-dipole interaction tensor
function J_q = compute_dipolar_tensor(qpoints, params)
    % Compute dipole-dipole interaction in momentum space
    % J_ij(q) = D * [δ_ij - 3*q_i*q_j/q²] * exp(-q²ξ²)
    
    n_q = size(qpoints, 2);
    J_q = zeros(3, 3, n_q);
    
    % Dipolar coupling strength
    % For magnetic dipoles: D = (μ₀/4π) * g²μ_B² / a³
    mu_B = 5.788e-5;  % Bohr magneton in eV/T
    mu_0_4pi = 1e-7;  % μ₀/4π in SI
    g = params.g_factor;
    a = params.a * 1e-10;  % Convert Angstrom to meters
    
    % Convert to energy units (meV)
    D0 = mu_0_4pi * (g * mu_B)^2 / a^3 * 1000;
    
    % For LiHoF4 specifically
    if isfield(params, 'D0')
        D0 = params.D0;  % Use specified value
    else
        D0 = 0.04;  % Typical value in meV for Ho³⁺ in LiHoF4
    end
    
    for iq = 1:n_q
        q = qpoints(:, iq);
        q_norm = norm(q);
        
        if q_norm < 1e-10
            % Special case for q=0: use trace-preserving form
            J_q(:,:,iq) = D0 * diag([1, 1, -2]) / 3;
        else
            q_hat = q / q_norm;
            
            % Dipolar tensor: δ_ij - 3*q_i*q_j/q²
            for i = 1:3
                for j = 1:3
                    if i == j
                        J_q(i,j,iq) = 1 - 3*q_hat(i)*q_hat(j);
                    else
                        J_q(i,j,iq) = -3*q_hat(i)*q_hat(j);
                    end
                end
            end
            
            % Multiply by coupling and form factor
            % Form factor accounts for finite size of dipoles
            form_factor = 1;  % Can add exp(-q²*r₀²) for finite size
            
            % Ewald summation contribution (for accurate long-range)
            if isfield(params, 'use_ewald') && params.use_ewald
                % Add Ewald correction for periodic boundaries
                alpha = params.ewald_alpha;
                form_factor = form_factor * exp(-(q_norm/(2*alpha))^2);
            end
            
            J_q(:,:,iq) = D0 * J_q(:,:,iq) * form_factor;
        end
        
        % Account for multiple Ho sites in unit cell
        if isfield(params, 'n_sites') && params.n_sites > 1
            % Structure factor for multiple sublattices
            structure_factor = compute_structure_factor(q, params);
            J_q(:,:,iq) = J_q(:,:,iq) * structure_factor;
        end
    end
end

%% Function 3: Improved mixing for tensor case
function K_new = tensor_anderson_mixing(K_history, F_history, iter, params)
    % Anderson acceleration adapted for 3x3 tensor mixing
    % K_history: cell array of previous K tensors
    % F_history: cell array of F(K) = K_new from self-consistency
    
    m = min(params.mixing_depth, iter-1);  % Depth of history to use
    
    if m < 2
        % Not enough history: simple mixing
        alpha = params.mixing_alpha;
        K_new = (1-alpha) * K_history{end} + alpha * F_history{end};
        return;
    end
    
    % Flatten tensors for linear algebra
    [n1, n2, n3] = size(K_history{end});
    n_total = n1 * n2 * n3;
    
    % Build residual matrices
    R = zeros(n_total, m);
    for i = 1:m
        K_flat = reshape(K_history{end-m+i}, [], 1);
        F_flat = reshape(F_history{end-m+i}, [], 1);
        R(:,i) = F_flat - K_flat;
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
        K_new = (1-alpha_simple) * K_history{end} + alpha_simple * F_history{end};
        return;
    end
    
    % Construct new K from linear combination
    K_new = zeros(size(K_history{end}));
    for i = 1:m
        K_new = K_new + alpha(i) * K_history{end-m+i};
    end
    
    % Add relaxation
    beta = params.anderson_beta;  % typically 0.5-1.0
    F_avg = zeros(size(F_history{end}));
    for i = 1:m
        F_avg = F_avg + alpha(i) * F_history{end-m+i};
    end
    
    K_new = (1-beta) * K_new + beta * F_avg;
end

%% Function 4: Check physical constraints for convergence
function [is_physical, diagnostics] = check_physical_constraints(K, G_local, params)
    % Check if solution satisfies physical constraints
    
    diagnostics = struct();
    is_physical = true;
    
    % 1. Hermiticity check
    max_herm = 0;
    for iw = 1:size(K,3)
        herm_error = norm(K(:,:,iw) - K(:,:,iw)');
        max_herm = max(max_herm, herm_error);
    end
    diagnostics.hermiticity_violation = max_herm;
    if max_herm > 1e-6
        warning('K not Hermitian: max violation = %.3e', max_herm);
        is_physical = false;
    end
    
    % 2. Causality check (spectral function should be positive)
    % A(ω) = -Im[G(ω+iδ)]/π should be positive
    for iw = 1:size(G_local,3)
        eigenvals = eig(G_local(:,:,iw));
        if any(real(eigenvals) > 0)  % For Matsubara, Re[G] should be negative
            diagnostics.causality_violation = true;
            is_physical = false;
            break;
        end
    end
    
    % 3. High-frequency behavior
    % K should decay as 1/ω for large ω
    n_check = min(10, size(K,3));
    omega_tail = (2*(size(K,3)-n_check:size(K,3)-1) + 1) * pi / params.beta;
    K_tail_11 = squeeze(K(1,1,end-n_check+1:end));
    
    % Fit power law in log-log
    p = polyfit(log(omega_tail), log(abs(K_tail_11)), 1);
    diagnostics.high_freq_exponent = p(1);
    
    if abs(p(1) + 1) > 0.5  % Should be close to -1
        warning('Incorrect high-frequency behavior: exponent = %.2f', p(1));
        is_physical = false;
    end
    
    % 4. Sum rule check
    % Tr[∫dω G(ω)] should equal something physical
    spectral_sum = zeros(3,3);
    for iw = 1:size(G_local,3)
        spectral_sum = spectral_sum + G_local(:,:,iw);
    end
    spectral_sum = spectral_sum / params.beta;
    
    diagnostics.spectral_sum_trace = trace(spectral_sum);
    
    % 5. Symmetry checks for tetragonal system
    if isfield(params, 'symmetry') && strcmp(params.symmetry, 'tetragonal')
        % x and y components should be equal
        xy_diff = 0;
        for iw = 1:size(K,3)
            xy_diff = max(xy_diff, abs(K(1,1,iw) - K(2,2,iw)));
        end
        diagnostics.xy_symmetry_breaking = xy_diff;
        
        if xy_diff > 1e-6
            warning('Tetragonal symmetry broken: |K_xx - K_yy| = %.3e', xy_diff);
        end
    end
    
    return;
end

%% Function 5: Compute physical observables from converged solution
function observables = compute_observables(K, G_local, omega_n, params)
    % Extract physical quantities from the converged Green's function
    
    observables = struct();
    
    % 1. Static susceptibility tensor χ(q=0,ω=0)
    % Need to analytically continue iω_n → 0
    % For lowest Matsubara frequency as approximation:
    chi_static = -params.beta * G_local(:,:,1);
    observables.chi_static = chi_static;
    observables.chi_xx = chi_static(1,1);
    observables.chi_zz = chi_static(3,3);
    observables.anisotropy_ratio = real(chi_static(3,3) / chi_static(1,1));
    
    % 2. Effective moment from Curie-Weiss fit
    % μ_eff² = 3k_B * T * χ
    k_B = 8.617e-5;  % eV/K
    T = 1/params.beta;
    mu_eff_z = sqrt(3 * k_B * T * real(chi_static(3,3)));
    observables.mu_eff = mu_eff_z;
    
    % 3. Energy shift from interactions
    % ΔE = -(1/β) * Σ_n Tr[K(iω_n) * G_local(iω_n)]
    energy_shift = 0;
    for iw = 1:length(omega_n)
        energy_shift = energy_shift - trace(K(:,:,iw) * G_local(:,:,iw));
    end
    energy_shift = real(energy_shift) / params.beta;
    observables.energy_shift = energy_shift;
    
    % 4. Quasiparticle weight (if applicable)
    % Z = [1 - ∂Σ/∂ω|_{ω=0}]^{-1}
    if length(omega_n) > 2
        % Finite difference for derivative
        dK_domega = (K(:,:,2) - K(:,:,1)) / (omega_n(2) - omega_n(1));
        Z_matrix = inv(eye(3) - 1i * dK_domega);
        observables.Z_matrix = Z_matrix;
        observables.Z_avg = real(trace(Z_matrix)) / 3;
    end
    
    % 5. Critical temperature estimate (mean-field)
    % T_c when largest eigenvalue of χ(0,0) * J(0) = 1
    if isfield(params, 'J_q0')
        chi_J = chi_static * params.J_q0;
        lambda_max = max(real(eig(chi_J)));
        observables.lambda_max = lambda_max;
        observables.T_c_estimate = T / lambda_max;  % Very rough estimate
    end
    
    % 6. Print summary
    fprintf('\n=== Physical Observables ===\n');
    fprintf('Static susceptibility:\n');
    fprintf('  χ_xx = %.4f, χ_zz = %.4f\n', real(chi_static(1,1)), real(chi_static(3,3)));
    fprintf('  Anisotropy ratio χ_zz/χ_xx = %.2f\n', observables.anisotropy_ratio);
    fprintf('Effective moment: μ_eff = %.2f μ_B\n', mu_eff_z);
    fprintf('Energy shift: ΔE = %.4f meV\n', energy_shift * 1000);
    if isfield(observables, 'T_c_estimate')
        fprintf('Estimated T_c (mean-field): %.3f K\n', observables.T_c_estimate);
    end
    
    return;
end

%% Function 6: Compute structure factor for multiple sublattices
function S_q = compute_structure_factor(q, params)
    % For systems with multiple sites per unit cell
    % Default to 1 if not needed
    S_q = 1;
    
    if isfield(params, 'site_positions')
        % Sum exp(i*q·r) over sites
        S_q = 0;
        for site = 1:size(params.site_positions, 2)
            r = params.site_positions(:, site);
            S_q = S_q + exp(1i * dot(q, r));
        end
        S_q = abs(S_q)^2 / params.n_sites;
    end
end