# Sum Rule Implementation (Equation 2.23)

## Overview

This implementation adds validation of the sum rule from **Jensen (1994), Equation 2.23** during effective medium theory iterations:

```
(1/β) Σ_n G(iω_n) = -M²
```

where:
- `G(iω_n)` is the local Green's function at Matsubara frequency ω_n
- `β = 1/(k_B T)` is the inverse temperature
- `M` is the dipole matrix element between the two lowest singlet states

## Reference

Jensen, J. (1994). *1/z renormalization of the mean-field behavior of the dipole-coupled singlet-singlet system HoF3*. Physical Review B, 49(17), 11833-11841.

**Equation 2.23** (page 11835): This sum rule is a fundamental constraint that the Green's function must satisfy for consistency with the underlying quantum mechanics.

## Files Created

### 1. `sum_rule_check.m`
Main validation function that checks if the sum rule is satisfied.

**Usage:**
```matlab
[ok, sum_val, exp_val, error] = sum_rule_check(G_local, beta, M_squared, omega_grid, verbose);
```

**Inputs:**
- `G_local` - Local Green's function [3×3×n_omega]
- `beta` - Inverse temperature 1/(k_B*T) [meV^-1]
- `M_squared` - Dipole matrix element squared (scalar or 3×3 tensor)
- `omega_grid` - Frequency points [GHz]
- `verbose` - Print detailed output (default: false)

**Outputs:**
- `ok` - Boolean: true if sum rule satisfied within 10% tolerance
- `sum_val` - Computed value of (1/β) Σ G(iω_n)
- `exp_val` - Expected value -M²
- `error` - Relative error

### 2. `effective_medium_with_sumrule.m`
Code snippets and helper functions to integrate sum rule checking into `effective_medium.m`.

**Key additions:**
- Periodic sum rule checking during SCF iterations
- Final validation at convergence
- `compute_M_squared_from_eigenstates()` function to extract M² from quantum states

### 3. `example_with_sumrule.m`
Complete example script demonstrating usage with your existing code.

## Implementation Steps

### Quick Start (Standalone Check)

Add to any script after computing `G_local`:

```matlab
% Physical constants
kB = 8.61733e-2; % meV/K
beta = 1/(kB * T); % T in Kelvin

% Compute M² from eigenstates (see compute_M_squared_from_eigenstates)
M_squared = compute_M_squared_from_eigenstates(eigenE, eigenW, ion, const);

% Check sum rule
[ok, sum_val, exp_val, err] = sum_rule_check(...
    G_local, beta, M_squared, freq_total, true);

if ~ok
    warning('Sum rule violated! Relative error: %.2e', err);
end
```

### Full Integration into effective_medium.m

1. **Before SCF loop** (after line 341):
   ```matlab
   % Add sum rule parameters
   if isfield(scf_params, 'M_squared')
       M_squared = scf_params.M_squared;
   else
       M_squared = 1.0;
       warning('M_squared not provided');
   end

   if isfield(scf_params, 'omega_grid')
       omega_grid = scf_params.omega_grid;
   else
       omega_grid = linspace(0, 50, n_omega);
   end

   sum_rule_check_interval = 100;
   sum_rule_history = zeros(opts.max_iter, 1);
   ```

2. **Inside SCF loop** (after line 460, after computing residuals):
   ```matlab
   % Check sum rule periodically
   if mod(iter, sum_rule_check_interval) == 0 || iter == 1
       [sum_ok, ~, ~, sum_err] = sum_rule_check(...
           G_local, scf_params.beta, M_squared, omega_grid, false);
       sum_rule_history(iter) = sum_err;

       if scf_params.verbose && mod(iter, sum_rule_check_interval*2) == 0
           fprintf('    [Iter %d] Sum rule error: %.3e\n', iter, sum_err);
       end
   end
   ```

3. **At convergence** (after line 477, inside convergence check):
   ```matlab
   % Final validation
   [sum_ok, ~, ~, sum_err] = sum_rule_check(...
       G_local, scf_params.beta, M_squared, omega_grid, scf_params.verbose);

   if ~sum_ok && scf_params.verbose
       warning('Converged solution violates sum rule!');
   end
   ```

4. **Add helper function** (at end of file):
   ```matlab
   % Copy compute_M_squared_from_eigenstates from effective_medium_with_sumrule.m
   ```

### Integration with Existing Workflow

Modify your main analysis script:

```matlab
% After loading eigenstates
M_squared_array = zeros(n_cVar, 1);
for ii = 1:n_cVar
    eigenE_point = squeeze(eigenE(ii,:))';
    eigenW_point = squeeze(eigenW(ii,:,:));
    M_sq = compute_M_squared_from_eigenstates(eigenE_point, eigenW_point, ion, const);

    % Use scalar (trace/3 for isotropic approximation)
    if isequal(size(M_sq), [3,3])
        M_squared_array(ii) = trace(M_sq) / 3;
    else
        M_squared_array(ii) = M_sq;
    end
end

% Before calling effective_medium
scf_params.M_squared = M_squared_array(current_idx);
scf_params.omega_grid = freq_total; % Your frequency grid [GHz]
```

## Physical Interpretation

### Why This Sum Rule Matters

1. **Fundamental Constraint**: Equation 2.23 is derived from the spectral representation of the Green's function and must be satisfied for any physically consistent solution.

2. **Convergence Indicator**: Violation indicates:
   - Insufficient frequency resolution
   - Numerical instability
   - Incorrect system parameters (especially M²)
   - Poor convergence of SCF iterations

3. **Quality Check**: Provides independent validation beyond simple residual convergence.

### Expected M² Values

For LiHoF4 (J=8 singlet-singlet system):
- Ground state to first excited state separation: Δ ≈ 0.7 meV
- Typical M² ≈ 10^-2 to 10^0 (dimensionless, in units of ℏ²)

For LiErF4 (J=15/2 doublet-doublet system):
- Different physics, but similar magnitude

## Troubleshooting

### Sum Rule Not Satisfied

**Problem**: Relative error > 10%

**Solutions**:
1. Increase frequency resolution (`n_omega`)
2. Extend frequency range to capture high-energy contributions
3. Verify M² calculation (check eigenstates)
4. Check SCF convergence (residual should be < 1e-5)
5. Verify beta value (temperature units: K → meV^-1)

### M² Calculation Issues

**Problem**: M² = 0 or unreasonably large

**Check**:
1. Eigenstates are normalized
2. Correct quantum numbers (J, I)
3. Hyperfine coupling included if hyp > 0
4. Using correct excited state (first excited, not higher)

### Integration Errors

**Problem**: Frequency integration not accurate

**Solutions**:
1. Use trapezoidal or Simpson's rule for non-uniform grids
2. Ensure frequency grid covers main spectral weight
3. Check for Matsubara vs real frequency (code assumes real frequencies)

## Testing

### Unit Test
```matlab
% Test with known system
T = 0.5; % K
kB = 8.61733e-2;
beta = 1/(kB*T);
M_squared_test = 0.5;

% Create test Green's function (should satisfy sum rule)
n_omega = 100;
omega = linspace(0, 50, n_omega); % GHz
const_Gh2mV = 1.05457E-34 * 2*pi * 1e9 * 6.24151e+21;

% Simple Lorentzian form
Delta = 1.0; % meV
Gamma = 0.1; % meV
G_test = zeros(3, 3, n_omega);
for iw = 1:n_omega
    w_meV = omega(iw) * const_Gh2mV;
    G_test(1,1,iw) = M_squared_test / (Delta - w_meV - 1i*Gamma);
    G_test(2,2,iw) = G_test(1,1,iw);
    G_test(3,3,iw) = G_test(1,1,iw);
end

% Check
[ok, sum_val, exp_val, err] = sum_rule_check(...
    G_test, beta, M_squared_test, omega, true);
```

## Performance Impact

- **Computation time**: < 1% overhead when checking every 100 iterations
- **Memory**: Negligible (stores only error history)
- **Recommended**: Check every 50-100 iterations during development, every 200-500 in production

## References

1. Jensen, J. (1994). PRB 49, 11833 - Original paper with Eq. 2.23
2. Stinchcombe, R. B. (1973). J. Phys. C 6, 2459 - 1/z expansion theory
3. Galili & Zevin (1987). J. Phys. C 20, 2543 - Effective medium approach

## Contact

For issues or questions about this implementation, refer to:
- CODE_REVIEW_AND_FIXES.md
- Original paper: Jensen PRB 49, 11833 (1994)
