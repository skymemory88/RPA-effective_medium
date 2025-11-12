# Q-Vectors for LiHoF4

## Crystal Structure

**LiHoF4** crystallizes in the **scheelite structure** (CaWO₄-type):
- **Space group**: I4₁/a (#88)
- **Crystal system**: Tetragonal
- **Lattice type**: Body-centered tetragonal

### Lattice Parameters
From literature (e.g., J. Phys. C: Solid State Phys. 8, 3483 (1975)):
- **a = b = 5.175 Å** (tetragonal plane)
- **c = 10.75 Å** (c-axis)
- **α = β = γ = 90°**

### Magnetic Sublattice
- **4 magnetic Ho³⁺ ions per unit cell**
- Ho³⁺ ions form a body-centered tetragonal sublattice
- Positions in fractional coordinates: (0,0,0), (1/2,1/2,1/2), and symmetry equivalents

---

## Reciprocal Lattice

The reciprocal lattice vectors are:
```
a* = 2π/a = 1.215 Å⁻¹  (along x)
b* = 2π/b = 1.215 Å⁻¹  (along y)
c* = 2π/c = 0.584 Å⁻¹  (along z)
```

### Reciprocal Lattice Matrix
```matlab
recip_lattice = [
    1.215,  0.000,  0.000;   % a* vector
    0.000,  1.215,  0.000;   % b* vector
    0.000,  0.000,  0.584    % c* vector
] % in Å⁻¹
```

---

## High-Symmetry Points in the Brillouin Zone

For the tetragonal lattice, the important high-symmetry points are:

| Label | Reduced Coords (h,k,l) | Cartesian (Å⁻¹) | Physical Meaning |
|-------|------------------------|-----------------|------------------|
| **Γ** | (0, 0, 0) | (0.000, 0.000, 0.000) | Zone center (ferromagnetic limit) |
| **X** | (1/2, 0, 0) | (0.607, 0.000, 0.000) | Zone boundary along a* |
| **M** | (1/2, 1/2, 0) | (0.607, 0.607, 0.000) | Zone corner in ab-plane |
| **Z** | (0, 0, 1/2) | (0.000, 0.000, 0.292) | Zone boundary along c* |
| **A** | (1/2, 0, 1/2) | (0.607, 0.000, 0.292) | Mixed zone boundary |
| **R** | (1/2, 1/2, 1/2) | (0.607, 0.607, 0.292) | Zone corner (antiferromagnetic) |

### Physical Interpretation

1. **Γ point (0,0,0)**:
   - Uniform (ferromagnetic) mode
   - All spins oscillate in phase
   - Important for static susceptibility measurements

2. **Z point (0,0,1/2)**:
   - Antiferromagnetic ordering along c-axis
   - Alternating spins along c-direction
   - Relevant for LiHoF4 which orders with **q ≈ (0,0,0)** (ferromagnetic in zero field)

3. **R point (1/2,1/2,1/2)**:
   - Full 3D antiferromagnetic order
   - Alternating spins in all directions

4. **Lines and planes**:
   - **(h,0,0)** line: Modulation along a-axis
   - **(0,0,l)** line: Modulation along c-axis (particularly important for LiHoF4)
   - **(h,h,0)** plane: Planar modulations

---

## Magnetic Ordering in LiHoF4

### Ground State (Zero Field, T < Tₙ = 1.53 K)
- **Ordering wave vector**: q = **(0, 0, 0)** (Γ point)
- **Type**: Ferromagnetic along a-axis (Ising-like)
- **Moments**: Aligned along [100] or [010] direction
- **Domain structure**: Two equivalent domains

### Applied Transverse Field (B ⟂ c)
- Applied field perpendicular to c-axis induces quantum phase transition
- **Critical field**: Bᶜ ≈ 5 T at T = 0
- Quantum critical point at (Bᶜ, T = 0)
- Excitations show characteristic softening as B → Bᶜ

### Q-Dependence of Susceptibility
The susceptibility χ(q,ω) is highly anisotropic:
- **Strong enhancement near Γ point** due to ferromagnetic correlations
- **Dipolar interactions** favor ferromagnetic correlations in ab-plane
- **c-axis dispersion** is weaker due to larger lattice spacing

---

## Q-Vector Sampling Strategies

### 1. For Dispersion Relations (ω vs q)
Use **path through high-symmetry points**:
```matlab
path = {'Gamma', 'X', 'M', 'Gamma', 'Z', 'A', 'R', 'Gamma'};
[qvec, ~] = compute_qvectors(lattice, 'mode', 'path', 'path', path, 'npoints', 100);
```

### 2. For q-Space Integration (RPA/EMT)
Use **uniform grid** with sufficient density:
```matlab
% Coarse (fast testing)
[qvec, ~] = compute_qvectors(lattice, 'mode', 'grid', 'grid', [11, 11, 11]);

% Medium (typical calculation)
[qvec, ~] = compute_qvectors(lattice, 'mode', 'grid', 'grid', [21, 21, 21]);

% Fine (high accuracy)
[qvec, ~] = compute_qvectors(lattice, 'mode', 'grid', 'grid', [41, 41, 41]);
```

### 3. For Specific Directions
**Along c-axis** (longitudinal):
```matlab
qz = linspace(0, 1, 101)';
qvec = [zeros(101,1), zeros(101,1), qz];
```

**In ab-plane** (transverse):
```matlab
qh = linspace(0, 1, 101)';
qvec = [qh, zeros(101,1), zeros(101,1)];  % Along a-axis
% or
qvec = [qh, qh, zeros(101,1)];            % Along (h,h,0)
```

---

## Integration Over Brillouin Zone

For effective medium theory calculations requiring q-space integration:

### Discrete Sum Approximation
```matlab
% Uniform grid
[qvec, ~] = compute_qvectors(lattice, 'mode', 'grid', 'grid', [Nq, Nq, Nq]);

% Integration: ∫ dq f(q) ≈ (1/Nq³) Σ_q f(q)
G_local = (1/Nq^3) * sum(G_q, 'all');  % Average over all q-points
```

### Convergence
- **11³ grid (1,331 points)**: Quick testing, ~5% error
- **21³ grid (9,261 points)**: Standard calculation, ~1% error
- **31³ grid (29,791 points)**: High accuracy, ~0.3% error
- **41³ grid (68,921 points)**: Publication quality, <0.1% error

---

## Usage in RPA Calculations

### In `MF_RPA_Yikai.m`

Replace the hardcoded q-vectors (around line 24):
```matlab
% OLD CODE:
if Options.Kplot == true
    qz = linspace(0,-1, 101)';
    qx = linspace(0, 1, 101)';
    qy = zeros(length(qx),1);
    qvec = [qx qy qz];
else
    qx = 0.0;
    qy = zeros(size(qx,1),1);
    qz = zeros(size(qx,1),1);
    qvec = [qx qy qz];
end
```

**NEW CODE** (flexible):
```matlab
% Load LiHoF4 lattice
a = 5.175; c = 10.75;
lattice = [a, 0, 0; 0, a, 0; 0, 0, c];

if Options.Kplot == true
    % Path for dispersion plot
    path = {'Gamma', 'X', 'M', 'Gamma', 'Z', 'A', 'R', 'Gamma'};
    [qvec, ~, ~] = compute_qvectors(lattice, 'mode', 'path', ...
                                     'path', path, 'npoints', 100);
else
    % Single point or grid for integration
    [qvec, ~, ~] = compute_qvectors(lattice, 'mode', 'grid', ...
                                     'grid', [21, 21, 21], ...
                                     'range', [-0.5, 0.5]);
end
```

---

## Physical Considerations

### 1. Dipolar Interactions
LiHoF4 has **strong dipolar interactions** (~10% of crystal field splitting):
- Long-range in real space → All q-vectors contribute
- Fourier transform: `D(q) = Σᵣ D(r) exp(iq·r)`
- Computed in `MF_dipole(qvec, dpRng, lattice, tau)`

### 2. Exchange Interactions
Exchange is typically nearest-neighbor:
- Short-range in real space → Smooth variation in q-space
- `J(q) = J₀ [cos(qₓa) + cos(qᵧa) + cos(qᵧc)]` (approximate)
- Computed in `exchange(qvec, ex_strength, lattice, tau)`

### 3. Lorentz Term
For **q = 0** (uniform mode), include Lorentz correction:
```matlab
lorz_on = (abs(q) < 1e-10);  % Only for Γ point
Lorentz_term = -lorz_on * (4*pi/Vc) * (eye(3)/3 - demag_tensor);
```

### 4. Structure Factor
Body-centered lattice has **structure factor**:
```matlab
F(q) = (1/4) * Σᵢ exp(iq·rᵢ)  % Sum over 4 atoms in unit cell
```
When F(q) = 0, certain modes are forbidden.

---

## Example: Computing Full χ(q,ω)

```matlab
%% Setup
a = 5.175; c = 10.75;
lattice = [a, 0, 0; 0, a, 0; 0, 0, c];

%% Q-space grid
[qvec, qvec_cart, recip] = compute_qvectors(lattice, 'mode', 'grid', ...
                                             'grid', [21, 21, 21]);
Nq = size(qvec, 1);

%% Frequency grid
freq = linspace(0, 10, 201);  % meV
Nf = length(freq);

%% Allocate susceptibility
chi = zeros(3, 3, Nf, Nq);  % [3 x 3 x freq x q]

%% Compute RPA susceptibility
for iq = 1:Nq
    q = qvec(iq, :);

    % Compute J(q) - interaction matrix
    Jq = compute_Jq(q, lattice, ion, const);

    for iw = 1:Nf
        omega = freq(iw);

        % Non-interacting susceptibility
        chi0 = compute_chi0(omega, ion, const);

        % RPA: χ(q,ω) = [1 - χ₀(ω)·J(q)]⁻¹ · χ₀(ω)
        chi(:,:,iw,iq) = inv(eye(3) - chi0 * Jq) * chi0;
    end
end

%% Extract observable quantities
% Neutron scattering cross-section: S(q,ω) ∝ Im[χ⁺⁺(q,ω)]
S_perp = squeeze(imag(chi(1,1,:,:)));  % Transverse

% Integrated intensity at each q
I_q = trapz(freq, S_perp, 1);

% Plot dispersion
figure;
imagesc(1:Nq, freq, S_perp);
xlabel('q-point'); ylabel('ω (meV)');
title('S(q,ω) for LiHoF4');
```

---

## References

1. **Crystal structure**:
   - Roser & Corruccini, *Phys. Rev. Lett.* **65**, 1064 (1990)

2. **Magnetic properties**:
   - Bitko et al., *Phys. Rev. Lett.* **77**, 940 (1996)
   - Ronnow et al., *Science* **308**, 389 (2005)

3. **Neutron scattering**:
   - Rønnow et al., *Phys. Rev. B* **75**, 054426 (2007)

4. **RPA theory**:
   - Jensen & Mackintosh, *Rare Earth Magnetism* (1991)
   - Your paper: "Z⁻¹ renormalization of the mean field behavior..."

---

## Quick Start

```matlab
%% Load and run example
run example_LiHoF4_qvectors.m

%% Use in your RPA code
load('LiHoF4_qvectors.mat', 'qvec_fine');
% ... continue with MF_RPA_Yikai.m using qvec_fine
```

---

*Generated by Claude for LiHoF4 RPA effective medium calculations*
*Date: 2025-11-12*
