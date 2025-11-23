%% Example: Run effective medium with Sum Rule Check (Eq. 2.23)
% This script demonstrates how to use the sum rule validation
% during effective medium iterations

clear; close all;

%% Setup parameters (modify for your system)
mion = 'Ho'; % or 'Er'
scanMode = 'temp';
temp = [0.15, 0.2, 0.3, 0.5, 1.0]; % K
fields = 0.0; % T
freq_total = linspace(0, 50, 201); % GHz
theta = 0; % degrees
phi = 0; % degrees
gama = 1e-3; % GHz
hyp = 1.0; % hyperfine proportion
RPA_mode = true;

%% Run MF/RPA to get initial susceptibilities
fprintf('Computing MF/RPA susceptibilities...\n');
[cVar, freq_total, chi0, chiq, qvec, Jq, dscrt_var] = MF_RPA_Yikai(...
    mion, scanMode, dscrt_var, freq_total, theta, phi, gama, hyp, RPA_mode);

%% Load eigenstate data for M² calculation
% This is needed to compute the dipole matrix element
if hyp > 0
    nZee_path = 'Hz_I=1';
else
    nZee_path = 'Hz_I=0';
end

if ispc
    location = ['C:\Users\skyme\OneDrive - Nexus365\Postdoc\Research projects\Li',mion,...
        'F4 project\Data\Simulations\mean field\eigen-states\', nZee_path, '\'];
else
    location = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Postdoc/Research projects/',...
        'Li', mion,'F4 project/Data/Simulations/mean field/eigen-states/', nZee_path, '/'];
end

switch scanMode
    case 'field'
        filename = strcat(['Hscan_Li',mion,'F4_'],...
            sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var, theta, phi, hyp),'.mat');
    case 'temp'
        filename = strcat(['Tscan_Li',mion,'F4_'],...
            sprintf('%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var, theta, phi, hyp),'.mat');
end

file = fullfile(location, filename);
load(file,'-mat','eee','vvv','ion');

eigenE = eee;
eigenW = vvv;

%% Setup constants
const.hbar = 1.05457E-34;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 1e9 * const.J2meV;
const.muB = 9.27401e-24;
const.muN = 5.05078e-27;
const.mu0 = 4*pi*1e-7;
const.dpRng = 100;
const.elem = find(ion.prop);
const.ELEf = ion.gLande(const.elem) * const.muB * const.J2meV;
const.NUCf = ion.nLande(const.elem) * const.muN * const.J2meV;
const.gfac = ion.gLande(const.elem)^2 * (const.muB)^2 * (const.mu0/4/pi) * const.J2meV * 10^30;

%% Compute M² for each temperature/field point
n_cVar = length(cVar);
M_squared_all = zeros(n_cVar, 1); % Store M² for each point

fprintf('\nComputing dipole matrix elements M²...\n');
for ii = 1:n_cVar
    eigenE_point = squeeze(eigenE(ii,:))';
    eigenW_point = squeeze(eigenW(ii,:,:));

    % Compute matrix element between ground and first excited state
    M_sq = compute_M_squared_from_eigenstates(eigenE_point, eigenW_point, ion, const);

    % For isotropic approximation, use trace/3
    if isequal(size(M_sq), [3,3])
        M_squared_all(ii) = trace(M_sq) / 3;
    else
        M_squared_all(ii) = M_sq;
    end

    fprintf('  %s = %.3f: M² = %.4e\n', scanMode, cVar(ii), M_squared_all(ii));
end

%% Setup for effective medium with sum rule validation
USE_SINGLE_ION_SEED = false; % Use RPA as seed

% Modified effective_medium.m parameters
kB = 8.61733e-2; % meV/K

% Get dimensions
n_omega = size(chi0, 3);
switch scanMode
    case 'field'
        n_cVar = length(fields);
        cVar_use = fields;
        dscrt_var_use = temp;
    case 'temperature'
        n_cVar = length(temp);
        cVar_use = temp;
        dscrt_var_use = fields;
end
n_q = size(qvec,1);

%% Run effective medium with sum rule checks
fprintf('\n=== Effective Medium with Sum Rule Validation ===\n');

% Storage
K_emt = zeros(3, 3, n_omega, n_cVar);
G_local_emt = zeros(3, 3, n_omega, n_cVar);
chi_emt = zeros(3, 3, n_omega, n_cVar);
sum_rule_errors = zeros(n_cVar, 1);
sum_rule_satisfied = false(n_cVar, 1);

% SCF parameters
scf_params_base = struct();
scf_params_base.max_iter = 5e3;
scf_params_base.tol = 1e-5;
scf_params_base.mixing_alpha = 0.05;
scf_params_base.G_damp = 0.1;
scf_params_base.verbose = true; % Enable verbose for sum rule output

% Loop over temperature/field
for ii = 1:n_cVar
    cVar_val = cVar_use(ii);

    switch scanMode
        case 'field'
            T_local = dscrt_var_use;
            beta_local = 1 / (kB * T_local);
            descriptor = sprintf('B = %.3f T, T = %.3f K', cVar_val, T_local);
        case 'temp'
            T_local = cVar_val;
            beta_local = 1 / (kB * T_local);
            descriptor = sprintf('T = %.3f K, B = %.3f T', T_local, dscrt_var_use);
    end

    fprintf('\n--- Processing %s ---\n', descriptor);

    % Prepare SCF parameters with sum rule info
    chi_seed_slice = chiq(:,:,:,ii,:);
    scf_params = scf_params_base;
    scf_params.beta = beta_local;
    scf_params.n_omega = n_omega;
    scf_params.n_q = n_q;
    scf_params.M_squared = M_squared_all(ii); % Add M² for sum rule check
    scf_params.omega_grid = freq_total; % Add frequency grid
    scf_params.G0 = -chiq(:,:,:,ii,:); % G = -χ for initial guess
    if ndims(scf_params.G0) == 5
        scf_params.G0 = mean(scf_params.G0, 5); % Average over q
    end
    scf_params.J_q = Jq(:,:,:,:);

    % Run effective medium (you'll need to modify this to include sum rule)
    % For now, this is a placeholder - the actual implementation requires
    % modifying the attempt_convergence function
    [K_sol, G_sol, converged] = compute_effective_medium_with_sumrule(scf_params, descriptor);

    % Store results
    K_emt(:,:,:,ii) = K_sol;
    G_local_emt(:,:,:,ii) = G_sol;

    % Compute susceptibility
    for iw = 1:n_omega
        chi_emt(:,:,iw,ii) = -G_sol(:,:,iw) / beta_local;
    end

    % Final sum rule check
    [sum_ok, ~, ~, sum_err] = check_sum_rule_eq223(...
        G_sol, beta_local, M_squared_all(ii), freq_total, true);

    sum_rule_errors(ii) = sum_err;
    sum_rule_satisfied(ii) = sum_ok;
end

%% Summary
fprintf('\n=== Sum Rule Validation Summary ===\n');
fprintf('%s points: %d\n', scanMode, n_cVar);
fprintf('Satisfied: %d / %d (%.1f%%)\n', ...
    sum(sum_rule_satisfied), n_cVar, 100*sum(sum_rule_satisfied)/n_cVar);
fprintf('Mean error: %.3e\n', mean(sum_rule_errors));
fprintf('Max error:  %.3e\n', max(sum_rule_errors));

% Plot sum rule error vs temperature/field
figure('Name','Sum Rule Error');
plot(cVar_use, sum_rule_errors, 'o-', 'LineWidth', 2);
xlabel(sprintf('%s', scanMode));
ylabel('Sum Rule Relative Error');
title('Equation 2.23 Sum Rule Validation');
grid on;
set(gca, 'YScale', 'log');

%% Helper function (simplified version - modify attempt_convergence for full integration)
function [K, G_local, converged] = compute_effective_medium_with_sumrule(scf_params, var_str)
    % This is a placeholder - in practice, you should modify
    % the attempt_convergence function in effective_medium.m
    % to include the sum rule checks as shown in effective_medium_with_sumrule.m

    fprintf('  Note: Full sum rule integration requires modifying attempt_convergence()\n');
    fprintf('  See effective_medium_with_sumrule.m for implementation details\n');

    % For now, return dummy results
    K = zeros(3, 3, scf_params.n_omega);
    G_local = scf_params.G0;
    converged = false;
end

function M_squared = compute_M_squared_from_eigenstates(eigenE, eigenW, ion, const)
    % Compute dipole matrix element M² between ground and first excited state

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

    % Hybridize
    JhT_x = JhT_x + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT_x;
    JhT_y = JhT_y + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT_y;
    JhT_z = JhT_z + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT_z;

    % Matrix elements
    v0 = eigenW(:, 1);
    v1 = eigenW(:, 2);

    M_x = abs(v0' * JhT_x * v1)^2;
    M_y = abs(v0' * JhT_y * v1)^2;
    M_z = abs(v0' * JhT_z * v1)^2;

    % Return diagonal tensor
    M_squared = diag([M_x, M_y, M_z]);
end
