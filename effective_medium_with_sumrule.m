% Modification to effective_medium.m to include sum rule check (Eq. 2.23)
%
% ADD THIS SECTION to the attempt_convergence function after line 341
% (after "residual_history = zeros(opts.max_iter, 1);")
%
% This implements equation 2.23 from Jensen (1994):
% (1/β) Σ_n G(iω_n) = -M²

%% Additional SCF parameters for sum rule validation
% Extract M² (dipole matrix element squared) from system parameters
% For LiHoF4/LiErF4, this comes from the matrix elements
if isfield(scf_params, 'M_squared')
    M_squared = scf_params.M_squared;
else
    % Default: estimate from system if not provided
    % For singlet-singlet systems, M ≈ g_J * μ_B * J_matrix_element
    % This should be provided by the calling code based on eigenstate analysis
    M_squared = 1.0; % Placeholder - should be computed from eigenstates
    warning('M_squared not provided in scf_params. Using default value %.3e', M_squared);
end

% Get omega grid if available
if isfield(scf_params, 'omega_grid')
    omega_grid = scf_params.omega_grid; % [GHz]
else
    % Create a default grid based on number of frequency points
    omega_max = 50; % GHz, typical for rare-earth systems
    omega_grid = linspace(0, omega_max, n_omega);
end

% Sum rule check frequency (check every N iterations)
sum_rule_check_interval = 100; % Check every 100 iterations
sum_rule_history = zeros(opts.max_iter, 1); % Store sum rule error history

%% INSERT THIS SECTION inside the main SCF loop
% Place this AFTER "Step E: Compute residuals" (around line 460)
% and BEFORE "Check convergence"

% ===== SUM RULE CHECK (Eq. 2.23) =====
% Check sum rule periodically during iteration
if mod(iter, sum_rule_check_interval) == 0 || iter == 1
    [sum_ok, sum_val, exp_val, sum_err] = check_sum_rule_eq223(...
        G_local, scf_params.beta, M_squared, omega_grid, false);

    sum_rule_history(iter) = sum_err;

    if scf_params.verbose && mod(iter, sum_rule_check_interval*2) == 0
        fprintf('    [Iter %d] Sum rule error: %.3e ', iter, sum_err);
        if sum_ok
            fprintf('✓\n');
        else
            fprintf('(not satisfied)\n');
        end
    end
end
% ===== END SUM RULE CHECK =====

%% INSERT THIS SECTION after convergence is achieved
% Place this AFTER the "if residual < scf_params.tol" convergence check
% and BEFORE "break"

% ===== FINAL SUM RULE VALIDATION =====
% Perform detailed sum rule check on converged solution
if scf_params.verbose
    fprintf('\n--- Final Sum Rule Validation ---\n');
end

[sum_ok, sum_val, exp_val, sum_err] = check_sum_rule_eq223(...
    G_local, scf_params.beta, M_squared, omega_grid, scf_params.verbose);

if ~sum_ok && scf_params.verbose
    warning('Converged solution does not satisfy sum rule (Eq. 2.23) within tolerance!');
    fprintf('This may indicate:\n');
    fprintf('  1. Insufficient frequency resolution\n');
    fprintf('  2. Incorrect M² value\n');
    fprintf('  3. Numerical instability in Green''s function\n');
end
% ===== END FINAL VALIDATION =====

%% HELPER FUNCTION: Compute M² from eigenstates
% Add this function to the end of effective_medium.m file

function M_squared = compute_M_squared_from_eigenstates(eigenE, eigenW, ion, const)
% COMPUTE_M_SQUARED_FROM_EIGENSTATES Extract dipole matrix element
%
% For a singlet-singlet system, M is the matrix element of J_x between
% the two lowest singlet states.
%
% Inputs:
%   eigenE - Eigenvalues [meV]
%   eigenW - Eigenvectors (columns are eigenstates)
%   ion    - Ion structure with quantum numbers
%   const  - Physical constants structure
%
% Output:
%   M_squared - |<0|J_x|1>|² (scalar for isotropic, 3x3 for anisotropic)

    % Get J operators
    ionJ = ion.J(const.elem);
    Jz = diag(ionJ:-1:-ionJ);
    Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ)).*(ionJ+1+((ionJ-1):-1:-ionJ))),1);
    Jm = Jp';

    % Include hyperfine if present
    if ion.hyp(const.elem) > 0
        ionI = ion.I(const.elem);
        Iz = diag(ionI:-1:-ionI);
        IhT_z = kron(eye(2*ionJ+1), Iz);
        Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1);
        Im = Ip';
        Iph = kron(eye(2*ionJ+1), Ip);
        Imh = kron(eye(2*ionJ+1), Im);
        IhT_x = (Iph+Imh)/2;
        IhT_y = (Iph-Imh)/2i;

        Jph = kron(Jp, eye(2*ionI+1));
        Jmh = kron(Jm, eye(2*ionI+1));
        JhT_x = (Jph+Jmh)/2;
        JhT_y = (Jph-Jmh)/2i;
        JhT_z = kron(Jz, eye(2*ionI+1));
    else
        JhT_x = (Jp+Jm)/2;
        JhT_y = (Jp-Jm)/2i;
        JhT_z = Jz;
        IhT_x = 0;
        IhT_y = 0;
        IhT_z = 0;
    end

    % Hybridize electron-nuclear operators
    JhT_x = JhT_x + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT_x;
    JhT_y = JhT_y + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT_y;
    JhT_z = JhT_z + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT_z;

    % Matrix elements between ground state |0⟩ and first excited state |1⟩
    v0 = eigenW(:, 1); % Ground state
    v1 = eigenW(:, 2); % First excited state

    % Compute matrix elements for each component
    M_x = abs(v0' * JhT_x * v1)^2;
    M_y = abs(v0' * JhT_y * v1)^2;
    M_z = abs(v0' * JhT_z * v1)^2;

    % Return as diagonal tensor (anisotropic case)
    % or as maximum value (conservative estimate for isotropic case)
    M_squared = diag([M_x, M_y, M_z]);

    % Alternative: use trace or maximum for scalar
    % M_squared = max([M_x, M_y, M_z]); % Most conservative
    % M_squared = (M_x + M_y + M_z) / 3; % Average

    fprintf('Matrix elements: M_x² = %.4e, M_y² = %.4e, M_z² = %.4e\n', M_x, M_y, M_z);
end

%% USAGE EXAMPLE
% To use this enhancement in your main script:
%
% 1. Add M_squared calculation before calling effective_medium:
%    M_squared = compute_M_squared_from_eigenstates(eigenE(1,:)', eigenW(1,:,:), ion, const);
%    scf_params.M_squared = M_squared;
%    scf_params.omega_grid = freq_total; % Your frequency grid in GHz
%
% 2. The sum rule will be checked automatically during iteration
%
% 3. Results will be displayed at convergence if verbose=true
