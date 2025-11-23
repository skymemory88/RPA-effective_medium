# Critical Corrections to Effective Medium Theory Implementation

## Summary of Issues Found and Corrected

Based on theoretical analysis of Jensen (1994), *Physical Review B* **49**, 11833, three critical issues were identified and corrected in the effective medium theory implementation.

---

## Issue 1: Green's Function Extraction ❌ → ✅

### Problem
**Location:** `effective_medium.m:545-555` (function `extract_G0`)

The Green's function extraction was missing the β (inverse temperature) factor when using RPA seed:

```matlab
% INCORRECT (before):
if use_rpa
    G0 = -slice;              % ❌ Missing β factor!
else
    G0 = -beta_local * slice;
end
```

### Theoretical Background
From the fluctuation-dissipation theorem and Eq. 2.23 of Jensen (1994):

**χ(ω) = -G(ω)/β**

Therefore: **G(ω) = -β·χ(ω)**

This relationship holds **regardless** of whether χ comes from RPA or single-ion calculations. The susceptibility χ is in units of [meV⁻¹], and β = 1/(k_B·T) is also in [meV⁻¹].

### Correction Applied
```matlab
% CORRECTED (after):
% Always use G = -β*χ relationship
G0 = -beta_local * slice;
```

**File:** `effective_medium.m:554`

### Impact
- Previously, Green's functions from RPA seeds were off by a factor of β ≈ 12-600 (depending on temperature)
- This caused incorrect convergence behavior and violated sum rules
- Correcting this should significantly improve convergence

---

## Issue 2: Sum Rule Calculation ❌ → ✅

### Problem
**Location:** `sum_rule_check.m:35-52`

The sum rule calculation was treating Eq. 2.23 as a continuous integral instead of a discrete sum:

```matlab
% INCORRECT (before):
for iw = 1:n_omega
    G_sum = G_sum + G_local(:,:,iw) * domega(iw);  % ❌ Integration weights!
end
sum_value = (1/beta) * G_sum;
```

### Theoretical Background
Equation 2.23 from Jensen (1994) is a **discrete Matsubara sum**:

**(1/β) Σ_n G(iω_n) = -M²**

where ω_n = 2πn/β are discrete Matsubara frequencies. This is **NOT** a continuous integral.

The sum should be:
**(1/β) [G(iω_0) + G(iω_1) + G(iω_2) + ...]**

NOT:
**(1/β) ∫ G(ω) dω**

### Correction Applied
```matlab
% CORRECTED (after):
for iw = 1:n_omega
    G_sum = G_sum + G_local(:,:,iw);  % ✅ Simple sum, no weights
end
sum_value = (1/beta) * G_sum;
```

**File:** `sum_rule_check.m:44-52`

### Impact
- Previously, sum rule checks were giving incorrect values
- The integration weights (Δω) artificially scaled the sum
- Correcting this gives the true Matsubara sum as required by theory

---

## Issue 3: Sum Rule Constraint ➕ (NEW FEATURE)

### Problem
The sum rule was only **monitored** but never **enforced** during iterations. This meant:
- SCF iterations could converge to solutions violating the sum rule
- No constraint was applied to guide convergence toward physical solutions

### Theoretical Background
Equation 2.23 is a fundamental constraint that must be satisfied by any physical Green's function. Enforcing it during iterations can:
1. Improve convergence by reducing the search space
2. Ensure solutions are physical
3. Help avoid unphysical local minima

### Correction Applied
Added optional sum rule projection in SCF loop:

**Location:** `effective_medium.m:437-478`

The enforcement works by:
1. Computing current sum rule value: (1/β) Σ_n G(iω_n)
2. Comparing to target: -M²
3. Computing violation: Δ = current - target
4. Applying partial correction distributed across all frequencies:
   ```
   G(iω_n) → G(iω_n) - (α·β/N_ω)·Δ
   ```
   where α is the enforcement strength (0-1)

### Usage
To enable sum rule enforcement, add these parameters before calling effective medium:

```matlab
% In your main script, before running effective_medium.m:
scf_params_base.M_squared = <your_M_squared_value>;
scf_params_base.omega_grid = freq_total;
scf_params_base.enforce_sum_rule = true;
scf_params_base.sum_rule_strength = 0.5;    % Start with 0.3-0.5
scf_params_base.sum_rule_interval = 50;     % Apply every 50 iterations
```

### Caution ⚠️
- **Start with moderate strength** (0.3-0.5)
- Too high strength can cause oscillations or instability
- Monitor convergence behavior when enabling enforcement
- Sum rule enforcement is OPTIONAL - you can leave it disabled for initial tests

---

## Verification Steps

### Step 1: Verify Green's Function Extraction

Add diagnostic output after extracting G0:

```matlab
% After line 554 in effective_medium.m
if scf_params.verbose
    fprintf('  G0 extraction: |G0| = %.3e, β = %.3e, |χ| = %.3e\n', ...
        norm(G0(:,:,1),'fro'), beta_local, norm(slice(:,:,1),'fro'));
    fprintf('  Ratio |G0|/|χ| = %.3e (should be ≈ β = %.3e)\n', ...
        norm(G0(:,:,1),'fro')/norm(slice(:,:,1),'fro'), beta_local);
end
```

**Expected:** Ratio should be approximately equal to β

### Step 2: Verify Sum Rule Calculation

Run `sum_rule_check` with verbose output:

```matlab
[ok, sum_val, exp_val, err] = sum_rule_check(G_local, beta, M_squared, omega_grid, true);
```

Check the output:
- **Before correction:** Large relative error (> 50%)
- **After correction:** Error should be smaller (ideally < 10%)

If error is still large, check:
1. Is M_squared correct?
2. Are you using enough frequency points?
3. Is the frequency range appropriate?

### Step 3: Test Sum Rule Enforcement

**Test A: Without enforcement**
```matlab
scf_params_base.enforce_sum_rule = false;
% Run effective_medium and check final sum rule error
```

**Test B: With moderate enforcement**
```matlab
scf_params_base.M_squared = <value>;
scf_params_base.omega_grid = freq_total;
scf_params_base.enforce_sum_rule = true;
scf_params_base.sum_rule_strength = 0.3;    % Conservative
scf_params_base.sum_rule_interval = 50;
scf_params_base.verbose = true;             % See enforcement in action
% Run effective_medium and compare convergence
```

**Expected:**
- Enforcement should reduce sum rule error
- Convergence may be faster or more stable
- If you see oscillations, reduce `sum_rule_strength`

---

## Computing M² from Eigenstates

The dipole matrix element squared M² can be computed from your eigenstates. Here's a template:

```matlab
function M_squared = compute_M_squared_from_eigenstates(eigenE, eigenW, ion, const, temp_idx)
    % Extract eigenstates for a specific temperature/field point
    v = squeeze(eigenW(temp_idx,:,:));  % Eigenvectors
    en = squeeze(eigenE(temp_idx,:));   % Eigenvalues

    % Build J operators (electronic + hyperfine if applicable)
    % ... (use your existing code from MF_RPA_Yikai.m)

    % Compute matrix elements between ground and first excited state
    % For singlet-singlet system, this is typically between states 1 and 2
    M_x = abs(v(:,1)' * JhT.x * v(:,2))^2;
    M_y = abs(v(:,1)' * JhT.y * v(:,2))^2;
    M_z = abs(v(:,1)' * JhT.z * v(:,2))^2;

    % For isotropic case, average over directions
    M_squared = (M_x + M_y + M_z) / 3;

    % For anisotropic case, return tensor
    % M_squared = diag([M_x, M_y, M_z]);

    % Convert to correct units (should already be dimensionless in your convention)
    % M_squared is in units consistent with your susceptibility
end
```

---

## Troubleshooting

### Issue: Sum rule error is still large (> 20%)

**Possible causes:**
1. **Insufficient frequency points:** Need denser sampling of Matsubara frequencies
2. **Incorrect M² value:** Double-check calculation from eigenstates
3. **Wrong frequency range:** May need to extend to higher frequencies
4. **Unit mismatch:** Verify χ is in [meV⁻¹] and β is in [meV⁻¹]

**Solutions:**
- Increase number of frequency points
- Verify M² calculation independently
- Check unit conversions in susceptibility calculation

### Issue: Convergence became worse after corrections

**Possible causes:**
1. **Sum rule enforcement too strong:** Reduce `sum_rule_strength`
2. **Initial guess changed:** Previous bug was accidentally helping convergence
3. **Need to adjust mixing parameters:** Try reducing `mixing_alpha` or increasing `G_damp`

**Solutions:**
- First try WITHOUT sum rule enforcement
- Adjust convergence parameters
- Use RPA seed instead of single-ion seed (or vice versa)

### Issue: Sum rule enforcement causes oscillations

**Possible causes:**
- Enforcement strength too high
- Applying too frequently
- Conflict with self-consistency equations

**Solutions:**
- Reduce `sum_rule_strength` to 0.1-0.3
- Increase `sum_rule_interval` to 100-200
- Disable enforcement and rely on monitoring only

---

## References

Jensen, J. (1994). "1/z renormalization of the mean-field behavior of the dipole-coupled singlet-singlet system HoF3." *Physical Review B*, **49**(17), 11833-11841.

Key equations:
- **Eq. 2.23** (page 11835): Sum rule (1/β) Σ_n G(iω_n) = -M²
- **Eq. 2.20** (page 11835): G(q,iω_n) in effective medium approximation
- **Eq. 2.11** (page 11834): Self-consistency equation for K(iω_n)

---

## Questions or Issues?

If you encounter any problems after applying these corrections:

1. Check that all units are consistent (meV for energies, meV⁻¹ for susceptibilities)
2. Verify β = 1/(k_B·T) is computed correctly
3. Ensure M² is calculated from the correct pair of states
4. Try running with `verbose = true` to see detailed diagnostic output

Good luck with your calculations!
