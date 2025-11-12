# Integration Guide: Using compute_qvectors with MF_RPA_Yikai.m

This guide shows you how to integrate the new `compute_qvectors()` function into your existing RPA code.

---

## Quick Start

### 1. Test the Function
```matlab
cd /path/to/RPA-effective_medium
run test_qvectors.m
```

This will verify that the function works correctly.

### 2. Run the Example
```matlab
run example_LiHoF4_qvectors.m
```

This generates and visualizes all q-vectors for LiHoF4 and saves them to `LiHoF4_qvectors.mat`.

---

## Modifying MF_RPA_Yikai.m

### Option 1: Minimal Changes (Quick Integration)

In `MF_RPA_Yikai.m`, around **line 24**, replace:

```matlab
% OLD CODE (lines 24-38):
if Options.Kplot == true
    qz = linspace(0,-1, 101)';
    qx = linspace(0, 1, 101)';
    qy = zeros(length(qx),1);
    qvec = [qx qy qz];
    cVar0 = [0.44]; % selected points of the continuous variable for k-plot
else
    qx = 0.0;
    qy = zeros(size(qx,1),1);
    qz = zeros(size(qx,1),1);
    qvec = [qx qy qz];
end
```

With:

```matlab
% NEW CODE - Using compute_qvectors for LiHoF4:

% Define LiHoF4 lattice parameters
a_LiHoF4 = 5.175;  % Angstrom
c_LiHoF4 = 10.75;  % Angstrom
lattice = [a_LiHoF4, 0, 0; 0, a_LiHoF4, 0; 0, 0, c_LiHoF4];

if Options.Kplot == true
    % Generate path for dispersion plot
    % Common path: Gamma -> X -> M -> Gamma -> Z -> A
    path_points = {'Gamma', 'X', 'M', 'Gamma', 'Z', 'A'};
    [qvec, ~, ~] = compute_qvectors(lattice, 'mode', 'path', ...
                                     'path', path_points, ...
                                     'npoints', 101);
    cVar0 = [0.44]; % selected points of the continuous variable for k-plot
else
    % Use high-symmetry points or grid for integration
    % Option A: Just Gamma point (original behavior)
    qvec = [0.0, 0.0, 0.0];

    % Option B: Use uniform grid for q-space integration
    % [qvec, ~, ~] = compute_qvectors(lattice, 'mode', 'grid', ...
    %                                  'grid', [11, 11, 11], ...
    %                                  'range', [-0.5, 0.5]);
end
```

### Option 2: Comprehensive Integration

Create a new section at the beginning of `MF_RPA_Yikai.m` (after line 100):

```matlab
%% Q-vector setup for LiHoF4
function [qvec, lattice, recip_lattice] = setup_qvectors_LiHoF4(Options)
    % Setup q-vectors based on LiHoF4 crystal structure

    % LiHoF4 lattice parameters (tetragonal scheelite structure)
    a = 5.175;  % Angstrom
    c = 10.75;  % Angstrom
    lattice = [a, 0, 0; 0, a, 0; 0, 0, c];

    if Options.Kplot
        % Path for plotting dispersion relations
        path_labels = {'Gamma', 'X', 'M', 'Gamma', 'Z', 'A', 'R', 'Gamma'};
        [qvec, ~, recip_lattice] = compute_qvectors(lattice, ...
                                                      'mode', 'path', ...
                                                      'path', path_labels, ...
                                                      'npoints', 101);
    elseif isfield(Options, 'q_integration') && Options.q_integration
        % Dense grid for q-space integration (EMT, RPA corrections)
        grid_density = Options.q_grid_size;  % e.g., 21
        [qvec, ~, recip_lattice] = compute_qvectors(lattice, ...
                                                      'mode', 'grid', ...
                                                      'grid', [grid_density, grid_density, grid_density], ...
                                                      'range', [-0.5, 0.5]);
    elseif isfield(Options, 'q_custom')
        % Custom q-vectors provided by user
        qvec = Options.q_custom;
        [~, ~, recip_lattice] = compute_qvectors(lattice);
    else
        % Default: just Gamma point
        qvec = [0.0, 0.0, 0.0];
        [~, ~, recip_lattice] = compute_qvectors(lattice);
    end
end
```

Then in the main code (around line 24), replace with:

```matlab
% Setup q-vectors
[qvec, lattice_LiHoF4, recip_lattice] = setup_qvectors_LiHoF4(Options);

% Store in ion structure if needed
if ~isfield(ion, 'abc')
    ion.abc = cell(1, length(ion.prop));
end
ion.abc{const.elem} = lattice_LiHoF4;
```

---

## Using with Effective Medium Theory

In `effective_medium.m`, when you need q-space integration:

```matlab
%% Setup q-grid for effective medium calculation
a = 5.175; c = 10.75;
lattice = [a, 0, 0; 0, a, 0; 0, 0, c];

% Choose grid density based on desired accuracy
% 11^3 = 1,331 points   (quick test)
% 21^3 = 9,261 points   (standard)
% 31^3 = 29,791 points  (high accuracy)
% 41^3 = 68,921 points  (publication quality)

grid_size = 21;  % Standard choice
[qvec, qvec_cart, recip_lattice] = compute_qvectors(lattice, ...
                                                     'mode', 'grid', ...
                                                     'grid', [grid_size, grid_size, grid_size], ...
                                                     'range', [-1, 1]);
Nq = size(qvec, 1);

fprintf('Q-space integration using %d points\n', Nq);

%% Compute q-dependent quantities
% ... your existing code ...

% When integrating over q-space:
% G_local(iω) = (1/Nq) Σ_q G(q, iω)
G_local = zeros(3, 3, Nfreq);
for iq = 1:Nq
    for iw = 1:Nfreq
        G_local(:,:,iw) = G_local(:,:,iw) + G_q(:,:,iw,iq);
    end
end
G_local = G_local / Nq;  % Average over q-points
```

---

## Common Q-Vector Configurations

### Configuration 1: Dispersion Along [001]
```matlab
% For plotting χ(q,ω) along c-axis
qz = linspace(0, 1, 201)';  % Reduced coordinates
qvec = [zeros(201,1), zeros(201,1), qz];

% Or using compute_qvectors:
[qvec, ~] = compute_qvectors(lattice, 'mode', 'path', ...
                              'path', [0,0,0; 0,0,1], ...
                              'npoints', 201);
```

### Configuration 2: Dispersion in ab-Plane
```matlab
% Path: Gamma -> X -> M -> Gamma
path = {'Gamma', 'X', 'M', 'Gamma'};
[qvec, ~] = compute_qvectors(lattice, 'mode', 'path', ...
                              'path', path, 'npoints', 100);
```

### Configuration 3: Full 3D Grid
```matlab
% For computing structure factor or full susceptibility
[qvec, qvec_cart] = compute_qvectors(lattice, 'mode', 'grid', ...
                                      'grid', [31, 31, 31], ...
                                      'range', [-1, 1]);

% The range [-1, 1] covers the first Brillouin zone plus some extension
```

### Configuration 4: High-Symmetry Points Only
```matlab
% Just compute at special q-points
[qvec, qvec_cart] = compute_qvectors(lattice, 'mode', 'highsym');

% Or specific points:
[qvec, qvec_cart] = compute_qvectors(lattice, 'mode', 'highsym', ...
                                      'sympoints', {'Gamma', 'X', 'Z', 'M'});
```

---

## Verifying Results

### Check 1: Gamma Point
The Gamma point should always be at (0,0,0):
```matlab
[qvec, ~] = compute_qvectors(lattice);
assert(all(qvec(1,:) == [0,0,0]), 'Gamma point should be first');
```

### Check 2: Reciprocal Lattice
For LiHoF4 with a=5.175 Å, c=10.75 Å:
```matlab
[~, ~, recip] = compute_qvectors(lattice);
expected_a_star = 2*pi/5.175;  % = 1.215 Å⁻¹
expected_c_star = 2*pi/10.75;  % = 0.584 Å⁻¹

assert(abs(norm(recip(1,:)) - expected_a_star) < 1e-3, 'a* incorrect');
assert(abs(norm(recip(3,:)) - expected_c_star) < 1e-3, 'c* incorrect');
```

### Check 3: Grid Coverage
Verify that grid covers the desired range:
```matlab
[qvec, ~] = compute_qvectors(lattice, 'mode', 'grid', ...
                              'grid', [11,11,11], 'range', [-0.5, 0.5]);

assert(min(qvec(:,1)) >= -0.5 && max(qvec(:,1)) <= 0.5, 'qx range');
assert(min(qvec(:,2)) >= -0.5 && max(qvec(:,2)) <= 0.5, 'qy range');
assert(min(qvec(:,3)) >= -0.5 && max(qvec(:,3)) <= 0.5, 'qz range');
```

---

## Performance Considerations

### Grid Size vs Computation Time

| Grid Size | # Points | RPA Time* | EMT Time* | Accuracy |
|-----------|----------|-----------|-----------|----------|
| 7³        | 343      | ~1 min    | ~5 min    | Low      |
| 11³       | 1,331    | ~5 min    | ~20 min   | Medium   |
| 21³       | 9,261    | ~35 min   | ~2 hrs    | Good     |
| 31³       | 29,791   | ~2 hrs    | ~6 hrs    | High     |
| 41³       | 68,921   | ~5 hrs    | ~15 hrs   | Very High|

*Approximate times on typical workstation, depends on frequency grid size

### Recommendations:
- **Development/Testing**: Use 7³ or 11³ grid
- **Standard Calculations**: Use 21³ grid
- **Publication Quality**: Use 31³ or 41³ grid
- **Convergence Test**: Run with increasing grid sizes until results change by < 1%

---

## Troubleshooting

### Issue 1: "Undefined function 'compute_qvectors'"
**Solution**: Make sure `compute_qvectors.m` is in your MATLAB path:
```matlab
addpath('/path/to/RPA-effective_medium');
```

### Issue 2: "Unknown symmetry point: XYZ"
**Solution**: Check the lattice type. Available points depend on detected lattice:
- Tetragonal (LiHoF4): Gamma, X, Y, M, Z, A, R
- Cubic: Gamma, X, M, R
- See `QVECTORS_LiHoF4.md` for full list

### Issue 3: Memory error with large grid
**Solution**: Reduce grid size or process q-points in batches:
```matlab
% Process in batches
batch_size = 1000;
for ibatch = 1:ceil(Nq/batch_size)
    q_start = (ibatch-1)*batch_size + 1;
    q_end = min(ibatch*batch_size, Nq);
    q_batch = qvec(q_start:q_end, :);

    % Process q_batch
    % ...
end
```

### Issue 4: RPA not converging at certain q-points
**Solution**: Check if q is at zone boundary where interactions diverge:
```matlab
% Add small damping for problematic q-points
q_mag = sqrt(sum(qvec_cart.^2, 2));
damping = 1e-6 * (q_mag > 0.9*max(q_mag));  % Damp near zone boundary
```

---

## Advanced Usage

### Custom Path Definition
```matlab
% Define custom path in reduced coordinates
custom_path = [
    0.0, 0.0, 0.0;   % Gamma
    0.3, 0.0, 0.0;   % Intermediate point
    0.5, 0.0, 0.0;   % X
    0.5, 0.5, 0.0;   % M
    0.0, 0.0, 0.5    % Z
];

[qvec, ~] = compute_qvectors(lattice, 'mode', 'path', ...
                              'path', custom_path, ...
                              'npoints', 50);
```

### Weighted Integration
```matlab
% For non-uniform grids, use proper weights
[qvec, qvec_cart] = compute_qvectors(lattice, 'mode', 'grid', ...
                                      'grid', [21,21,21]);

% Simple average (uniform grid)
G_avg = mean(G_q, 4);  % Average over 4th dimension (q-points)

% Or explicit sum with weights
Nq = size(qvec, 1);
weights = ones(Nq, 1) / Nq;  % Uniform weights
G_avg = sum(G_q .* reshape(weights, [1,1,1,Nq]), 4);
```

### Symmetry Reduction
```matlab
% For tetragonal lattice, use symmetry to reduce computation
% qx and qy are equivalent, so only need to compute half the grid

[qvec_full, ~] = compute_qvectors(lattice, 'mode', 'grid', ...
                                   'grid', [21,21,21]);

% Select irreducible q-points (qx >= qy >= 0, qz >= 0)
irreducible = (qvec_full(:,1) >= qvec_full(:,2)) & ...
              (qvec_full(:,2) >= 0) & ...
              (qvec_full(:,3) >= 0);
qvec_irred = qvec_full(irreducible, :);

fprintf('Reduced from %d to %d q-points using symmetry\n', ...
        size(qvec_full,1), size(qvec_irred,1));
% Typically reduces by factor of ~16 for tetragonal
```

---

## Summary Checklist

Before running calculations:

- [ ] Test `compute_qvectors` with `test_qvectors.m`
- [ ] Run `example_LiHoF4_qvectors.m` to understand output
- [ ] Choose appropriate q-grid size for your accuracy needs
- [ ] Modify `MF_RPA_Yikai.m` to use new q-vector generation
- [ ] Verify lattice parameters (a=5.175 Å, c=10.75 Å for LiHoF4)
- [ ] Check that high-symmetry points match expectations
- [ ] Run convergence test with increasing grid sizes
- [ ] Save q-vectors and results for documentation

---

## Files Created

1. **`compute_qvectors.m`** - Main function for q-vector generation
2. **`example_LiHoF4_qvectors.m`** - Comprehensive example and visualization
3. **`test_qvectors.m`** - Unit tests
4. **`QVECTORS_LiHoF4.md`** - Detailed documentation
5. **`INTEGRATION_GUIDE.md`** - This file

---

For questions or issues, refer to:
- `QVECTORS_LiHoF4.md` for theory and background
- Function help: `help compute_qvectors`
- Example script: `example_LiHoF4_qvectors.m`

*Last updated: 2025-11-12*
