function [sum_rule_satisfied, sum_value, expected_value, relative_error] = sum_rule_check(G_local, beta, M_squared, omega_grid, verbose)
% SUM_RULE_CHECK Validate the sum rule from equation 2.23
%
% Equation 2.23 from Jensen (1994): (1/β) Σ_n G(iω_n) = -M²
%
% Inputs:
%   G_local     - Local Green's function [3 x 3 x n_omega]
%   beta        - Inverse temperature 1/(k_B * T) [meV^-1]
%   M_squared   - Square of dipole matrix element M² [scalar or 3x3 for anisotropic]
%   omega_grid  - Frequency grid points [GHz] [1 x n_omega]
%   verbose     - (optional) Print diagnostic output (default: false)
%
% Outputs:
%   sum_rule_satisfied - Boolean: true if sum rule is satisfied within tolerance
%   sum_value         - Actual value of (1/β) * Σ_n G(iω_n) [3 x 3]
%   expected_value    - Expected value -M² [3 x 3 or scalar]
%   relative_error    - Relative error norm(sum_value + M²I) / norm(M²I)
%
% References:
%   Jensen, J. (1994). 1/z renormalization of the mean-field behavior of the
%   dipole-coupled singlet-singlet system HoF3. Physical Review B, 49(17), 11833.
%   Equation 2.23, page 11835.
%
% Example:
%   [ok, val, exp, err] = sum_rule_check(G_local, beta, M^2, freq, true);

if nargin < 5
    verbose = false;
end

% Get dimensions
[n_comp1, n_comp2, n_omega] = size(G_local);
assert(n_comp1 == 3 && n_comp2 == 3, 'G_local must be 3x3xN');

% CORRECTED: Use discrete Matsubara sum (Eq. 2.23)
% Equation 2.23: (1/β) Σ_n G(iω_n) = -M²
%
% This is a DISCRETE sum over Matsubara frequencies, NOT a continuous integral
% The Matsubara frequencies are: ω_n = 2πn/β where n is an integer
%
% In practice, we're evaluating G at a discrete set of frequencies,
% so we perform a simple discrete sum WITHOUT integration weights:

G_sum = zeros(3, 3);

for iw = 1:n_omega
    G_sum = G_sum + G_local(:,:,iw);  % Simple sum, NO domega weights
end

% Apply the 1/β factor to get the sum rule value
% This should equal -M² according to Eq. 2.23
sum_value = (1/beta) * G_sum;

% Expected value: -M²
if isscalar(M_squared)
    % Isotropic case: -M² * I
    expected_value = -M_squared * eye(3);
else
    % Anisotropic case: -M²_tensor
    expected_value = -M_squared;
end

% Compute relative error
residual = sum_value - expected_value;
relative_error = norm(residual, 'fro') / max(norm(expected_value, 'fro'), eps);

% Check if sum rule is satisfied (within 10% tolerance)
tolerance = 0.10; % 10% tolerance
sum_rule_satisfied = relative_error < tolerance;

% Verbose output
if verbose
    fprintf('\n=== Sum Rule Check (Eq. 2.23) ===\n');
    fprintf('Temperature: T = %.3f K (β = %.3e meV⁻¹)\n', 1/(8.61733e-2 * beta), beta);
    fprintf('Number of frequency points in discrete sum: %d\n', n_omega);
    fprintf('NOTE: Using discrete Matsubara sum (1/β)Σ_n G(iω_n), not continuous integral\n');

    if isscalar(M_squared)
        fprintf('M² = %.4e\n', M_squared);
        fprintf('\nExpected: -M² * I = \n');
        disp(expected_value);
    else
        fprintf('M² (tensor):\n');
        disp(M_squared);
        fprintf('\nExpected: -M² = \n');
        disp(expected_value);
    end

    fprintf('Computed: (1/β) Σ G(iω_n) = \n');
    disp(sum_value);

    fprintf('Residual: \n');
    disp(residual);

    fprintf('Relative error: %.4e (%.2f%%)\n', relative_error, relative_error*100);

    if sum_rule_satisfied
        fprintf('✓ Sum rule SATISFIED (error < %.0f%%)\n', tolerance*100);
    else
        fprintf('✗ Sum rule VIOLATED (error > %.0f%%)\n', tolerance*100);
        warning('Sum rule (Eq. 2.23) is not satisfied! Check convergence.');
    end
    fprintf('==================================\n\n');

    % Diagonal components comparison
    fprintf('Diagonal comparison:\n');
    fprintf('  Component    Expected       Computed       Error\n');
    fprintf('  ---------    --------       --------       -----\n');
    for ii = 1:3
        comp_name = char('x' + ii - 1);
        fprintf('  %s%s:          %+.4e     %+.4e     %.2e\n', ...
            comp_name, comp_name, ...
            expected_value(ii,ii), sum_value(ii,ii), ...
            abs(sum_value(ii,ii) - expected_value(ii,ii)));
    end
    fprintf('\n');
end

end
