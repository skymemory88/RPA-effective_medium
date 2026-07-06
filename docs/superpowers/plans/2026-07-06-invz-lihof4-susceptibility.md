# 1/z-Renormalized Susceptibility for LiHoF4 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute the dynamic susceptibility χ(q,ω) of LiHoF4 in the paramagnetic phase (thermal and quantum paramagnet, transverse field along a) one order beyond MF-RPA, using Jensen's high-density 1/z effective-medium expansion, reproducing the published benchmarks of Rønnow et al., PRB 75, 054426 (2007).

**Architecture:** A new self-contained `invz/` MATLAB module, Matsubara-native. Four layers: (1) single-ion — crystal field + hyperfine + Zeeman + transverse mean field, exact diagonalization (17 electronic / 136 electronuclear states), χ₀ tensor at arbitrary complex frequency from the spectral representation; (2) lattice — 4×4 sublattice cc-coupling matrix (recycled `MF_dipole`/`exchange`), reduced to 4 eigenvalue branches 𝒥ν(q); (3) scalar 1/z core — self-consistent loop G ↔ K ↔ (λ₁,λ₂) ↔ Σ = α + γg on the Matsubara axis; (4) observables — cc component of χ₀ divided by 1+Σ(ω), scalar mode-RPA on the real axis, phase boundary. Existing `src/` EMT solver is left untouched.

**Tech Stack:** MATLAB (function-based unit tests via `functiontests`/`runtests`, matching `tests/` conventions). No toolboxes required. Recycled from repo root: `MF_dipole.m`, `exchange.m`, `qVec_generator.m`.

## Global Constraints

- Units: energies in **meV**, temperatures in **K**, magnetic fields in **T**. `kB = 0.08617333 meV/K`, `muB = 0.05788382 meV/T`. β = 1/(kB·T) in meV⁻¹.
- Susceptibility convention: χ is the "spin" susceptibility built from bare J matrix elements (units meV⁻¹·J²); couplings 𝒥(q) in meV·J⁻². Green function **G = −χ**. RPA: χ_q = χ₀/(1−𝒥χ₀). Criticality: 𝒥_max(0)·χ₀cc(0) = 1+Σ(0).
- Matsubara frequencies bosonic: ωₙ = 2πn/β; all Green functions real and even in ωₙ; sums over n∈ℤ implemented as n≥0 with weights w₀=1, wₙ≥₁=2.
- Ion parameters (Rønnow 2007): J=8, gL=5/4, I=7/2, A=3.36e-3 meV; CF (meV): B20=−0.06, B40=0.35e-3, B44=3.6e-3, B60=4.0e-7, B64c=7.0e-5, B64s=0.98e-5; lattice a=b=5.175 Å, c=10.75 Å, 4 Ho/cell; 𝒥₁₂=−0.1e-3 meV (nn exchange); dipole prefactor gfac = μ0/4π·(gL·μB)² = 0.08388 meV·Å³.
- **Paramagnetic phase only** (⟨Jz⟩=0). Functions must error out, not silently continue, if called at an ordered-phase point. Ordered phase (α_m, ξ, H_MF machinery) is an explicit non-goal.
- Published benchmark targets: 𝒥_D·D_aa(0)=3.912 μeV, 𝒥_D·D_cc(0)=6.821 μeV, 𝒥(0)_eff=6.421 μeV; Σc(0) closed form = 0.3447 (fcc nn), 0.3004 (LiHoF4); Tc(B=0)=1.74 K; Hc(0.31 K)≈42.4–43 kOe with Σc(0)=0.0932; soft-mode minimum ≈0.19 meV; χ_J(0)/χ(0)≈0.77 at 0.31 K.
- Theory references by equation: `jensen_1z_framework.html` (HTML eq numbers) and Rønnow 2007 (R eq numbers). Core: λ_p (HTML 21), γ (HTML 23), α (HTML 27), Dyson G=G₀/(1+Σ+KG₀) (HTML 28 / R 9), effective medium (HTML 14–15 / R 7–8), closed-form critical Σ (R 10).
- Tests live in `invz/tests/`, run with `matlab -batch "results = runtests('invz/tests'); assertSuccess(results)"` from the repo root. Slow integration tests are guarded by `assumeTrue(strcmp(getenv('INVZ_SLOW'),'1'), ...)` so the default suite stays fast; the slow suite runs with `INVZ_SLOW=1`.
- Commit after every task (working tree already contains unrelated untracked files — `git add` only the files named in the task).

## Module Interface Contract (used across tasks)

- `C = invz_const()` → struct: `kB, muB, gfac, Gh2mV`.
- `ops = stevens_ops(J)` → struct of (2J+1)² matrices: `Jx,Jy,Jz,Jplus,Jminus,O20,O40,O44,O60,O64c,O64s`.
- `ion = invz_ion()` → parameter struct (fields listed in Task 2).
- `si = invz_single_ion(ion, T, B, opts)` → struct: `E` [n×1] meV (ascending, min subtracted), `V` [n×n], `P` [n×1] Boltzmann populations, `Mx,My,Mz` [n×n] matrix elements of electronic J in the eigenbasis, `Jexp` [3×1], `hx` (converged transverse mean field, meV), `JzJz_fluct` scalar ⟨J̃z²⟩. `opts.hyp` logical; n = 17 (false) or 136 (true).
- `chi = invz_chi0z(si, T, z, opts)` → [3,3,numel(z)] complex; z any complex row/col vector (Matsubara: `1i*wn`; real axis: `w + 1i*eta`); `opts.elastic` (default true) adds the degenerate/elastic term to entries with |z|<`opts.ztol`.
- `[wn, wts, beta] = invz_matsubara(T, Ecut)` → ωₙ column [nw×1], weights column, β.
- `tl = invz_twolevel(ion, T, Bx)` → struct: `Delta, M2, m, n01, g0` (electronic two-level parameters).
- `g = invz_g(tl, z)` → two-level response 2·n01·Δ/(Δ²−z²), same shape as z.
- `[Jnu, info] = invz_jq_modes(ion, qvec, opts)` → `Jnu` [n_q×4] real (meV), eigenvalues of the 4×4 cc sublattice coupling matrix per q; `info.Jcc0` = uniform-mode coupling at Γ.
- `Sc = invz_sigma_crit(J0, Jnu_flat)` → scalar: mean(Jnu/(J0−Jnu)) excluding |J0−Jnu|<1e-12 (R eq 10).
- `med = invz_emt_scalar(G0, Sigma, Jnu_flat, opts)` → struct: `G, K` [nw×1] real, `closure` scalar, `converged`, `iters`. Solves the K fixed point at fixed Σ.
- `lam = invz_lambdas(K, g, wts, beta, plist)` → column, λ_p = (1/β)Σₙ wₙ·Kₙ·gₙ^p.
- `sig = invz_sigma(tl, lam, K, g, beta)` → struct: `alpha` scalar, `gamma` [nw×1], `Sigma` [nw×1].
- `pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts)` → struct: `Sigma0, alpha, lambda, K, G, Sigma, tl, chi0cc0, crit, sumrule_rel, converged`.
- `bx = invz_critical(ion, T, Jnu_flat, opts)` → critical transverse field at temperature T (bisection on `pt.crit`); `invz_critical_T0field(ion, Sc, J0eff)` → zero-field Tc from the closed-form Σc.
- `out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)` → struct: `Sigma_w` [nw_real×1] complex, `chi_cc_q` [n_q_sel×4×nw_real], `chi0cc_w`.

---

### Task 1: Module scaffold, constants, Stevens operators

**Files:**
- Create: `invz/invz_const.m`
- Create: `invz/stevens_ops.m`
- Test: `invz/tests/test_invz_stevens.m`

**Interfaces:**
- Consumes: nothing.
- Produces: `invz_const()` and `stevens_ops(J)` exactly as in the contract above. Every later task calls both.

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_stevens
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_constants(testCase)
C = invz_const();
verifyEqual(testCase, C.kB,  0.08617333, 'RelTol', 1e-6);
verifyEqual(testCase, C.muB, 0.05788382, 'RelTol', 1e-6);
verifyEqual(testCase, C.gfac, 0.08388, 'RelTol', 2e-3);   % mu0/4pi*(gL*muB)^2 in meV*Ang^3
end

function test_angular_momentum_algebra(testCase)
ops = stevens_ops(8);
verifyEqual(testCase, size(ops.Jz), [17 17]);
comm = ops.Jx*ops.Jy - ops.Jy*ops.Jx;
verifyLessThan(testCase, max(abs(comm - 1i*ops.Jz), [], 'all'), 1e-12);
Jsq = ops.Jx^2 + ops.Jy^2 + ops.Jz^2;
verifyLessThan(testCase, max(abs(Jsq - 8*9*eye(17)), [], 'all'), 1e-10);
end

function test_stevens_J2_explicit(testCase)
% For J=2 (X=6): O20 diagonal = 3m^2-6 for m=2..-2; O44(1,5)=12 (=(1/2)*24)
ops = stevens_ops(2);
verifyEqual(testCase, diag(ops.O20), [6;-3;-6;-3;6], 'AbsTol', 1e-12);
verifyEqual(testCase, real(ops.O44(1,5)), 12, 'AbsTol', 1e-10);
for f = {'O20','O40','O44','O60','O64c','O64s'}
    M = ops.(f{1});
    verifyLessThan(testCase, max(abs(M - M'), [], 'all'), 1e-10);  % Hermitian
end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_stevens.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_const'".

- [ ] **Step 3: Implement `invz/invz_const.m`**

```matlab
function C = invz_const()
%INVZ_CONST Physical constants for the invz module (meV / K / T units).
C.kB     = 0.08617333;    % meV/K
C.muB    = 0.05788382;    % meV/T
C.Gh2mV  = 4.135667e-3;   % GHz -> meV
gL = 5/4;
C.gfac   = 0.08388;       % mu0/4pi*(gL*muB)^2 [meV*Ang^3]; check: gfac*4/Vc = 1.1654e-3 meV = J_D
end
```

- [ ] **Step 4: Implement `invz/stevens_ops.m`**

```matlab
function ops = stevens_ops(J)
%STEVENS_OPS Angular momentum and Stevens operator matrices for given J.
% Basis ordered mJ = J..-J. Conventions: Hutchings; O64s carries the 1/(4i) factor.
dim = round(2*J+1);
mJ  = (J:-1:-J).';
X   = J*(J+1);
Jz  = diag(mJ);
jp  = sqrt(X - mJ(2:end).*(mJ(2:end)+1));  % <m+1|J+|m>, m = mJ(2:end)
Jplus  = diag(jp, 1);
Jminus = Jplus';
E = eye(dim);
ops.Jz = Jz; ops.Jplus = Jplus; ops.Jminus = Jminus;
ops.Jx = (Jplus + Jminus)/2;
ops.Jy = (Jplus - Jminus)/(2i);
ops.O20 = 3*Jz^2 - X*E;
ops.O40 = 35*Jz^4 - (30*X - 25)*Jz^2 + (3*X^2 - 6*X)*E;
ops.O44 = 0.5*(Jplus^4 + Jminus^4);
ops.O60 = 231*Jz^6 - (315*X - 735)*Jz^4 + (105*X^2 - 525*X + 294)*Jz^2 ...
          + (-5*X^3 + 40*X^2 - 60*X)*E;
Cm = 11*Jz^2 - (X + 38)*E;
P4 = Jplus^4 + Jminus^4;  M4 = Jplus^4 - Jminus^4;
ops.O64c = 0.25 *(P4*Cm + Cm*P4);
ops.O64s = -0.25i*(M4*Cm + Cm*M4);
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_stevens.m'); assertSuccess(results)"`
Expected: PASS (3 passed).

- [ ] **Step 6: Commit**

```bash
git add invz/invz_const.m invz/stevens_ops.m invz/tests/test_invz_stevens.m
git commit -m "feat(invz): module scaffold, constants, Stevens operators"
```

---

### Task 2: Single-ion Hamiltonian (electronic + hyperfine + transverse mean field)

**Files:**
- Create: `invz/invz_ion.m`
- Create: `invz/invz_single_ion.m`
- Test: `invz/tests/test_invz_single_ion.m`

**Interfaces:**
- Consumes: `invz_const()`, `stevens_ops(J)`.
- Produces: `ion = invz_ion()` with fields `J, gL, I, A, B20, B40, B44, B60, B64c, B64s, a (3×3 Å), tau (4×3 fractional), J12, Jxx0, J0eff, Vc`; `si = invz_single_ion(ion, T, B, opts)` returning the struct defined in the contract (fields `E,V,P,Mx,My,Mz,Jexp,hx,JzJz_fluct`).

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_single_ion
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_zero_field_electronic_levels(testCase)
% Rønnow 2007 Table II: ground doublet, first excited singlet at 11 K.
ion = invz_ion(); C = invz_const();
si = invz_single_ion(ion, 1.0, [0 0 0], struct('hyp', false));
verifyEqual(testCase, numel(si.E), 17);
verifyLessThan(testCase, si.E(2) - si.E(1), 1e-9);                  % degenerate doublet
gap_K = (si.E(3) - si.E(1)) / C.kB;
verifyGreaterThan(testCase, gap_K, 9);  verifyLessThan(testCase, gap_K, 13);
% Ising moment of the doublet: eigvals of the 2x2 Jz block = ±<Jz>, g|| = 2*gL*<Jz> = 13.78
blk = si.Mz(1:2,1:2);
mu  = sort(abs(eig(blk)));
verifyEqual(testCase, mu(2), 5.51, 'RelTol', 0.10);
end

function test_hyperfine_manifold(testCase)
% 136 states; lowest 16 (doublet x nuclear octet) spread ~ 2*A*I*<Jz> ≈ 0.13 meV (1.5 K)
ion = invz_ion();
si = invz_single_ion(ion, 0.5, [0 0 0], struct('hyp', true));
verifyEqual(testCase, numel(si.E), 136);
spread = si.E(16) - si.E(1);
verifyGreaterThan(testCase, spread, 0.10);  verifyLessThan(testCase, spread, 0.16);
verifyGreaterThan(testCase, si.E(17) - si.E(16), 0.3);   % gap to next manifold (11 K = 0.95 meV scale)
end

function test_transverse_field_splitting_and_para(testCase)
ion = invz_ion();
d = zeros(1,3); bx = [1 3 6];
for k = 1:3
    si = invz_single_ion(ion, 0.31, [bx(k) 0 0], struct('hyp', false));
    d(k) = si.E(2) - si.E(1);
    verifyLessThan(testCase, abs(si.Jexp(3)), 1e-8);      % paramagnetic: <Jz>=0
    verifyLessThan(testCase, abs(si.Mz(1,1)), 1e-6);      % m ≈ 0
end
verifyGreaterThan(testCase, d(2), d(1));  verifyGreaterThan(testCase, d(3), d(2));
verifyGreaterThan(testCase, d(3), 0.1);   verifyLessThan(testCase, d(3), 1.0);
% mean field converged: hx nonzero under transverse field
si = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false));
verifyGreaterThan(testCase, abs(si.hx), 1e-6);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_single_ion.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_ion'".

- [ ] **Step 3: Implement `invz/invz_ion.m`**

```matlab
function ion = invz_ion()
%INVZ_ION LiHoF4 parameters (Rønnow et al., PRB 75, 054426 (2007), Table I last row).
ion.J   = 8;    ion.gL = 5/4;   ion.I = 3.5;
ion.A   = 3.36e-3;                                 % meV
ion.B20 = -0.06;      ion.B40 = 0.35e-3;  ion.B44 = 3.6e-3;    % meV
ion.B60 = 4.0e-7;     ion.B64c = 7.0e-5;  ion.B64s = 0.98e-5;  % meV
ion.a   = [5.175 0 0; 0 5.175 0; 0 0 10.75];       % Å (rows = lattice vectors)
ion.tau = [0 0 0; 0 0.5 0.25; 0.5 0.5 0.5; 0.5 0 0.75];  % Ho sites, fractional
ion.Vc  = abs(det(ion.a));                          % 287.9 Å^3
ion.J12 = -0.1e-3;                                  % meV, nn exchange (4 neighbours)
% Uniform-mode couplings (validated in Task 5): J_D*D_cc(0)=6.821e-3? NO: 6.821 μeV = 6.821e-3 meV
ion.J0eff = 6.421e-3;   % meV: J_D*D_cc(0) + 4*J12  (R 2007, after eq 11)
ion.Jxx0  = 3.512e-3;   % meV: J_D*D_aa(0) + 4*J12, transverse MF channel
end
```

- [ ] **Step 4: Implement `invz/invz_single_ion.m`**

```matlab
function si = invz_single_ion(ion, T, B, opts)
%INVZ_SINGLE_ION Exact single-ion diagonalization with transverse mean field.
% B = [Bx By Bz] in T. opts.hyp: include nuclear I=7/2 (dim 136). Paramagnetic phase:
% no <Jz> mean field is applied; caller must check si.Jexp(3) ≈ 0.
if nargin < 4, opts = struct(); end
hyp = isfield(opts,'hyp') && opts.hyp;
Jxx0 = ion.Jxx0; if isfield(opts,'Jxx0'), Jxx0 = opts.Jxx0; end
C   = invz_const();
oJ  = stevens_ops(ion.J);
Hcf = ion.B20*oJ.O20 + ion.B40*oJ.O40 + ion.B44*oJ.O44 ...
    + ion.B60*oJ.O60 + ion.B64c*oJ.O64c + ion.B64s*oJ.O64s;
if hyp
    oI = stevens_ops(ion.I);
    nI = size(oI.Jz,1);  nJ = size(oJ.Jz,1);
    kJ = @(M) kron(M, eye(nI));
    Hhf = ion.A*(kron(oJ.Jx,oI.Jx) + kron(oJ.Jy,oI.Jy) + kron(oJ.Jz,oI.Jz));
else
    kJ = @(M) M;  Hhf = 0;
end
Jx = kJ(oJ.Jx);  Jy = kJ(oJ.Jy);  Jz = kJ(oJ.Jz);
H0 = kJ(Hcf) + Hhf - ion.gL*C.muB*(B(1)*Jx + B(2)*Jy + B(3)*Jz);
beta = 1/(C.kB*T);
hx = 0;
for it = 1:200                                   % transverse mean-field fixed point
    H = H0 - hx*Jx;
    H = (H + H')/2;
    [V, D] = eig(H, 'vector');
    [E, ix] = sort(real(D));  V = V(:, ix);
    p = exp(-beta*(E - E(1)));  p = p/sum(p);
    jx = real(diag(V'*Jx*V)).'*p;
    hx_new = Jxx0*jx;
    if abs(hx_new - hx) < 1e-12, hx = hx_new; break; end
    hx = hx_new;
end
si.E  = E - E(1);
si.V  = V;
si.P  = p;
si.Mx = V'*Jx*V;  si.My = V'*Jy*V;  si.Mz = V'*Jz*V;
si.Jexp = [real(diag(si.Mx)).'*p; real(diag(si.My)).'*p; real(diag(si.Mz)).'*p];
si.hx = hx;
jz2 = real(diag(V'*(Jz*Jz)*V)).'*p;
si.JzJz_fluct = jz2 - si.Jexp(3)^2;
end
```

- [ ] **Step 5: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_single_ion.m'); assertSuccess(results)"`
Expected: PASS (3 passed). If `test_zero_field_electronic_levels` fails on the 11 K gap, first suspect the B64s sign/magnitude (both signs are physical per the paper — try `-0.98e-5`) before touching anything else.

- [ ] **Step 6: Commit**

```bash
git add invz/invz_ion.m invz/invz_single_ion.m invz/tests/test_invz_single_ion.m
git commit -m "feat(invz): LiHoF4 single-ion Hamiltonian with hyperfine and transverse mean field"
```

---

### Task 3: χ₀ at arbitrary complex frequency (spectral representation)

**Files:**
- Create: `invz/invz_chi0z.m`
- Test: `invz/tests/test_invz_chi0z.m`

**Interfaces:**
- Consumes: `invz_single_ion` output struct `si` (fields `E,P,Mx,My,Mz,Jexp`), `invz_const`.
- Produces: `chi = invz_chi0z(si, T, z, opts)` → [3,3,nz] complex. Kubo form χμν(z) = Σ_{a,b: |ΔE|>degtol} (p_a−p_b)·Mμ(a,b)·Mν(b,a)/(E_b−E_a−z), plus for grid entries with |z|<ztol the elastic term β·[Σ_{|ΔE|≤degtol} p_a·Mμ(a,b)·Mν(b,a) − ⟨Jμ⟩⟨Jν⟩]. Options: `degtol` (default 1e-8 meV), `ztol` (default 1e-12), `elastic` (default true).

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_chi0z
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_matsubara_equal_time_sum_rule(testCase)
% (1/beta)*sum_n chi_cc(iwn) = <Jz~^2>  (HTML eq 31 with chi = -G), electronuclear ion.
ion = invz_ion(); C = invz_const();
T = 4.0; beta = 1/(C.kB*T);
si = invz_single_ion(ion, T, [2 0 0], struct('hyp', true));
Ecut = 150;                                    % meV; CF spectrum tops out near 40 meV
nmax = ceil(Ecut*beta/(2*pi));
wn  = 2*pi*(0:nmax).'/beta;  wts = [1; 2*ones(nmax,1)];
chi = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
s   = real(squeeze(chi(3,3,:)));
lhs = sum(wts.*s)/beta;
verifyEqual(testCase, lhs, si.JzJz_fluct, 'RelTol', 0.03);
end

function test_even_real_on_matsubara_axis(testCase)
ion = invz_ion(); 
si = invz_single_ion(ion, 1.0, [3 0 0], struct('hyp', false));
wn = [0.3 0.9];
cp = invz_chi0z(si, 1.0,  1i*wn, struct());
cm = invz_chi0z(si, 1.0, -1i*wn, struct());
verifyLessThan(testCase, max(abs(imag(cp(3,3,:)))), 1e-12);
verifyEqual(testCase, cp(3,3,:), cm(3,3,:), 'AbsTol', 1e-12);
end

function test_static_positive_and_hyp_enhancement(testCase)
% chi_J(0)/chi(0) ≈ 0.77 at 0.31 K near the critical field (R 2007, Sec IIC).
ion = invz_ion();
T = 0.31; Bx = 4.24;
sJ = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false));
sF = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));
cJ = real(invz_chi0z(sJ, T, 0, struct()));  cJ = cJ(3,3);
cF = real(invz_chi0z(sF, T, 0, struct()));  cF = cF(3,3);
verifyGreaterThan(testCase, cJ, 0);
r = cJ/cF;
verifyGreaterThan(testCase, r, 0.68);  verifyLessThan(testCase, r, 0.86);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_chi0z.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_chi0z'".

- [ ] **Step 3: Implement `invz/invz_chi0z.m`**

```matlab
function chi = invz_chi0z(si, T, z, opts)
%INVZ_CHI0Z Single-ion susceptibility tensor at arbitrary complex frequencies.
% chi(mu,nu,iz) = sum_{a,b inelastic} (p_a-p_b) M_mu(a,b) M_nu(b,a) / (E_b-E_a - z(iz))
%              [+ elastic beta-term on entries with |z|<ztol].
if nargin < 4, opts = struct(); end
degtol = 1e-8;  if isfield(opts,'degtol'), degtol = opts.degtol; end
ztol   = 1e-12; if isfield(opts,'ztol'),   ztol   = opts.ztol;   end
elast  = true;  if isfield(opts,'elastic'), elast = opts.elastic; end
C = invz_const();  beta = 1/(C.kB*T);
E = si.E;  p = si.P;  n = numel(E);
dE = E.' - E;                    % dE(a,b) = E(b)-E(a)
dp = p - p.';                    % dp(a,b) = p(a)-p(b)
inel = abs(dE) > degtol;
M = {si.Mx, si.My, si.Mz};
z = z(:);  nz = numel(z);
chi = zeros(3,3,nz);
for mu = 1:3
    for nu = 1:3
        Nmn = M{mu} .* (M{nu}.');            % M_mu(a,b)*M_nu(b,a)
        w   = Nmn .* dp;  w(~inel) = 0;
        dEi = dE(inel);  wi = w(inel);
        for iz = 1:nz
            chi(mu,nu,iz) = sum(wi ./ (dEi - z(iz)));
        end
        if elast
            P2 = repmat(p, 1, n);                % P2(a,b) = p(a)
            el = beta*(sum(Nmn(~inel).*P2(~inel)) - si.Jexp(mu)*si.Jexp(nu));
            idx0 = abs(z) < ztol;
            chi(mu,nu,idx0) = chi(mu,nu,idx0) + el;
        end
    end
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_chi0z.m'); assertSuccess(results)"`
Expected: PASS (3 passed). The equal-time test carries a 3% tolerance for Matsubara truncation; if it misses, raise `Ecut`, and if it misses badly (>10%) the elastic term or the dp sign is wrong.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_chi0z.m invz/tests/test_invz_chi0z.m
git commit -m "feat(invz): spectral chi0 at complex frequency with elastic term"
```

---

### Task 4: Matsubara grid, two-level extraction, g identities

**Files:**
- Create: `invz/invz_matsubara.m`
- Create: `invz/invz_twolevel.m`
- Create: `invz/invz_g.m`
- Test: `invz/tests/test_invz_matsubara_twolevel.m`

**Interfaces:**
- Consumes: `invz_single_ion`, `invz_const`.
- Produces: `[wn,wts,beta] = invz_matsubara(T,Ecut)`; `tl = invz_twolevel(ion,T,Bx)` with `Delta,M2,m,n01,g0`; `g = invz_g(tl,z)`.

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_matsubara_twolevel
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_grid(testCase)
[wn, wts, beta] = invz_matsubara(2.0, 30);
verifyEqual(testCase, wn(1), 0);
verifyEqual(testCase, wn(2), 2*pi/beta, 'RelTol', 1e-12);
verifyEqual(testCase, wts(1), 1);  verifyEqual(testCase, wts(2), 2);
verifyGreaterThan(testCase, wn(end), 30);
end

function test_g_sum_identities(testCase)
% (1/beta)*sum_n g = 1  and  (1/beta)*sum_n g^2 = (1/2)[g(0)+beta(1-n01^2)]  (HTML eqs 26, Sec 7)
tl.Delta = 0.4;  T = 1.3;  C = invz_const();  beta = 1/(C.kB*T);
tl.n01 = tanh(beta*tl.Delta/2);  tl.M2 = 1;  tl.m = 0;  tl.g0 = 2*tl.n01/tl.Delta;
[wn, wts, ~] = invz_matsubara(T, 400*tl.Delta);
g = real(invz_g(tl, 1i*wn));
s1 = sum(wts.*g)/beta;
s2 = sum(wts.*g.^2)/beta;
verifyEqual(testCase, s1, 1, 'RelTol', 1e-3);
verifyEqual(testCase, s2, 0.5*(tl.g0 + beta*(1-tl.n01^2)), 'RelTol', 1e-6);
end

function test_twolevel_from_ion(testCase)
ion = invz_ion();  C = invz_const();
T = 0.31;  Bx = 4.0;
tl = invz_twolevel(ion, T, Bx);
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false));
verifyEqual(testCase, tl.Delta, si.E(2)-si.E(1), 'AbsTol', 1e-12);
verifyEqual(testCase, tl.M2, abs(si.Mz(1,2))^2, 'RelTol', 1e-12);
verifyLessThan(testCase, abs(tl.m), 1e-6);
verifyEqual(testCase, tl.n01, tanh(tl.Delta/(2*C.kB*T)), 'RelTol', 1e-12);
verifyGreaterThan(testCase, tl.M2, 1);       % sizeable Ising matrix element
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_matsubara_twolevel.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_matsubara'".

- [ ] **Step 3: Implement the three functions**

`invz/invz_matsubara.m`:
```matlab
function [wn, wts, beta] = invz_matsubara(T, Ecut)
%INVZ_MATSUBARA Bosonic Matsubara grid up to Ecut (meV), n>=0 with doubling weights.
C = invz_const();  beta = 1/(C.kB*T);
nmax = max(8, ceil(Ecut*beta/(2*pi)));
wn  = 2*pi*(0:nmax).'/beta;
wts = [1; 2*ones(nmax,1)];
end
```

`invz/invz_twolevel.m`:
```matlab
function tl = invz_twolevel(ion, T, Bx)
%INVZ_TWOLEVEL Electronic two-level (split doublet) parameters for the Jensen self-energy.
C  = invz_const();
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false));
tl.Delta = si.E(2) - si.E(1);
if tl.Delta < 1e-4
    error('invz:degenerateDoublet', ...
        'Doublet splitting %.2e meV too small: Bx=0 limit needs the closed-form Sigma_c (invz_sigma_crit).', tl.Delta);
end
tl.M2  = abs(si.Mz(1,2))^2;
tl.m   = real(si.Mz(1,1));
if abs(tl.m) > 1e-3
    error('invz:orderedPhase', 'Nonzero diagonal moment m=%.3g: outside paramagnetic-phase scope.', tl.m);
end
tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.g0  = 2*tl.n01/tl.Delta;
end
```

`invz/invz_g.m`:
```matlab
function g = invz_g(tl, z)
%INVZ_G Two-level inelastic response g(z) = 2*n01*Delta/(Delta^2 - z^2)  (HTML eq 11).
g = 2*tl.n01*tl.Delta ./ (tl.Delta^2 - z.^2);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_matsubara_twolevel.m'); assertSuccess(results)"`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add invz/invz_matsubara.m invz/invz_twolevel.m invz/invz_g.m invz/tests/test_invz_matsubara_twolevel.m
git commit -m "feat(invz): Matsubara grid, two-level parameters, g(z) with identity tests"
```

---

### Task 5: Coupling branches 𝒥ν(q) (recycled dipole/exchange sums)

**Files:**
- Create: `invz/invz_jq_modes.m`
- Create: `invz/cache/.gitignore` (content: `*` `!.gitignore`)
- Test: `invz/tests/test_invz_jq_modes.m`

**Interfaces:**
- Consumes: repo-root `MF_dipole(q, N, a, tau)` → [3,3,4,4] (Å⁻³) and `exchange(q, Jex, a, tau)` → [3,3,4,4]; `invz_const`, `invz_ion`.
- Produces: `[Jnu, info] = invz_jq_modes(ion, qvec, opts)` — `qvec` [n_q×3] reduced (h,k,l); `Jnu` [n_q×4] real meV (ascending per q); `info.Jcc0` scalar (uniform cc coupling at Γ incl. Lorentz term, no demag — sample-shape term cancels in the critical condition per R 2007); `opts.dpRng` (default 30), `opts.cache` (default true, key = sprintf of grid hash into `invz/cache/`).
- Before coding, copy the exact working call pattern of `MF_dipole`/`exchange` from `tests/test_dipole_batch.m` / `tests/test_exchange_batch.m` (they contain executable examples of the argument conventions).

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_jq_modes
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));   % repo root: MF_dipole, exchange
end

function test_gamma_point_constants(testCase)
% R 2007 eq (4): J_D*D_aa(0) = 3.912 ueV, J_D*D_cc(0) = 6.821 ueV;
% uniform cc mode with exchange: J(0) = 6.821 + 4*(-0.1) = 6.421 ueV.
ion = invz_ion();
[~, info] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 30, 'cache', false));
verifyEqual(testCase, info.Jcc0_dipole*1e3, 6.821e-3*1e3, 'RelTol', 0.03);  % meV -> ueV displayed
verifyEqual(testCase, info.Jaa0_dipole*1e3, 3.912e-3*1e3, 'RelTol', 0.03);
verifyEqual(testCase, info.Jcc0*1e3, 6.421e-3*1e3, 'RelTol', 0.03);
end

function test_modes_real_and_bounded(testCase)
ion = invz_ion();
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09];
[Jnu, info] = invz_jq_modes(ion, q, struct('dpRng', 20, 'cache', false));
verifyEqual(testCase, size(Jnu), [3 4]);
verifyLessThan(testCase, max(abs(imag(Jnu(:)))), 1e-12);
verifyLessThan(testCase, max(Jnu(:)), info.Jcc0 + 1e-4);   % Γ uniform mode is the max coupling
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_jq_modes.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_jq_modes'".

- [ ] **Step 3: Implement `invz/invz_jq_modes.m`**

```matlab
function [Jnu, info] = invz_jq_modes(ion, qvec, opts)
%INVZ_JQ_MODES Eigenvalue branches of the 4x4 cc sublattice coupling matrix (meV).
% J_cc(q)_{rs} = gfac*dip_cc_{rs}(q) [+ Lorentz at q≡0] + J12*ex_cc_{rs}(q)
% Convention: ferromagnetic-positive; criticality when J(0)*chi0 = 1+Sigma(0).
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng'), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;
C = invz_const();
cacheDir = fullfile(fileparts(mfilename('fullpath')), 'cache');
key = sprintf('jq_%d_%s.mat', dpRng, hash_qvec(qvec));
cacheFile = fullfile(cacheDir, key);
if useCache && exist(cacheFile, 'file')
    S = load(cacheFile);  Jnu = S.Jnu;  info = S.info;  return;
end
nq = size(qvec,1);
Jnu = zeros(nq, 4);
lorz = 4*pi/(3*ion.Vc)*C.gfac;   % Lorentz term per site pair block, meV (sphere, demag cancelled)
for iq = 1:nq
    q = qvec(iq,:);
    dip = MF_dipole(q, dpRng, ion.a, ion.tau);       % [3,3,4,4], Å^-3
    ex  = exchange(q, abs(ion.J12), ion.a, ion.tau); % [3,3,4,4], carries |J12|
    Jcc = squeeze(C.gfac*dip(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:));
    if is_gamma_equiv(q, ion.tau)
        Jcc = Jcc + lorz;                            % add Lorentz cavity term at Γ
    end
    Jcc = (Jcc + Jcc')/2;
    Jnu(iq,:) = sort(real(eig(Jcc))).';
end
% Γ-point info block (always computed at q=[0 0 0]):
dip0 = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
Jcc0d = squeeze(C.gfac*dip0(3,3,:,:)) + lorz;
Jaa0d = squeeze(C.gfac*dip0(1,1,:,:)) + lorz;
info.Jcc0_dipole = max(real(eig((Jcc0d+Jcc0d')/2)));
info.Jaa0_dipole = max(real(eig((Jaa0d+Jaa0d')/2)));
info.Jcc0 = info.Jcc0_dipole + 4*ion.J12;
info.dpRng = dpRng;
if useCache
    if ~exist(cacheDir,'dir'), mkdir(cacheDir); end
    save(cacheFile, 'Jnu', 'info');
end
end

function tf = is_gamma_equiv(q, tau)
tf = abs(real(sum(exp(2i*pi*(tau*q.'))))/size(tau,1) - 1) < 1e-9;
end

function h = hash_qvec(qvec)
h = sprintf('%dq_%08x', size(qvec,1), ...
    typecast(single(sum(qvec(:).*(1:numel(qvec))')), 'uint32'));
end
```

- [ ] **Step 4: Run, then resolve the Lorentz-term convention checkpoint**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_jq_modes.m'); assertSuccess(results)"`

This step is a **known physics checkpoint**, not a formality. The target is `info.Jcc0_dipole = 6.821e-3 meV`. Diagnosis table if it fails:
- Result ≈ `6.821e-3 − 4.88e-3 = 1.94e-3` → the Lorentz term is missing entirely (`lorz` not added, or `is_gamma_equiv` returning false at Γ): 𝒥_D·4π/3 = 4.881 μeV.
- Result ≈ `6.821e-3 + 4.88e-3` → Lorentz double-counted (brute-force sphere sum already contains it); drop `lorz`.
- Result off by an overall factor ≈ 0.084/0.0839·4 or /4 → per-cell vs per-site normalization of the 4×4 block sum; compare with `MF_RPA_Yikai.m` lines 823–827 (`J_avg = sum(...)/unitN`).
- `Jaa0` and `Jcc0` swapped → `tau` row/column or a-vs-c axis convention; c must be the third axis.
Record the resolution as a comment in `invz_jq_modes.m`. Expected after resolution: PASS (2 passed).

- [ ] **Step 5: Grid-convergence spot check (manual, no assert)**

Run in MATLAB:
```matlab
addpath invz; addpath .
ion = invz_ion();
for R = [10 20 40]
    [~, info] = invz_jq_modes(ion, [0 0 0], struct('dpRng', R, 'cache', false));
    fprintf('dpRng=%2d : Jcc0_dipole = %.4f ueV\n', R, info.Jcc0_dipole*1e3);
end
```
Expected: values converging to 6.82 μeV within ~1% by R=30. Note the chosen production `dpRng` in the file header comment.

- [ ] **Step 6: Commit**

```bash
git add invz/invz_jq_modes.m invz/cache/.gitignore invz/tests/test_invz_jq_modes.m
git commit -m "feat(invz): 4-branch cc coupling J_nu(q) with Gamma-point benchmarks"
```

---

### Task 6: BZ averages and the closed-form critical self-energy (R eq 10)

**Files:**
- Create: `invz/invz_sigma_crit.m`
- Test: `invz/tests/test_invz_sigma_crit.m`

**Interfaces:**
- Consumes: `invz_jq_modes`, `qVec_generator` (repo root).
- Produces: `Sc = invz_sigma_crit(J0, Jnu_flat)` — Σc(0) = mean over all (q,ν) of 𝒥/(J0−𝒥), entries with `J0−Jnu < 1e-12` excluded (uniform mode at Γ).

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_sigma_crit
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function test_fcc_watson(testCase)
% R 2007 below eq (10): fcc nearest-neighbour value 0.3447 (= Watson integral 1.34466 - 1).
n = 40;  f = (0:n-1)/n;
[F1,F2,F3] = ndgrid(f,f,f);
b = 2*pi*[-1 1 1; 1 -1 1; 1 1 -1];          % fcc primitive reciprocal vectors (a=1)
qx = F1*b(1,1)+F2*b(2,1)+F3*b(3,1);
qy = F1*b(1,2)+F2*b(2,2)+F3*b(3,2);
qz = F1*b(1,3)+F2*b(2,3)+F3*b(3,3);
Jq = 4*(cos(qx/2).*cos(qy/2) + cos(qy/2).*cos(qz/2) + cos(qz/2).*cos(qx/2));
Sc = invz_sigma_crit(12, Jq(:));
verifyEqual(testCase, Sc, 0.3447, 'AbsTol', 0.004);
end

function test_lihof4_sigma_crit(testCase)
% R 2007: Sigma_c(0; H=0) = 0.3004 with J12 = -0.1 ueV.  SLOW (dipole sums on a 3D grid).
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'size', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);          % drop Γ
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Sc = invz_sigma_crit(info.Jcc0, Jnu(:));
verifyEqual(testCase, Sc, 0.3004, 'AbsTol', 0.015);
end
```

- [ ] **Step 2: Run fast test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_sigma_crit.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_sigma_crit'" (the LiHoF4 test is filtered by the INVZ_SLOW assumption).

- [ ] **Step 3: Implement `invz/invz_sigma_crit.m`**

```matlab
function Sc = invz_sigma_crit(J0, Jnu_flat)
%INVZ_SIGMA_CRIT Closed-form critical self-energy, Rønnow 2007 eq (10):
%   Sigma_c(0) = (1/N) sum_{q,nu} J_nu(q) / (J(0) - J_nu(q)),
% valid at the zero-frequency critical point of the degenerate (elastic) doublet.
x = Jnu_flat(:);
keep = (J0 - x) > 1e-12;
if ~all(keep)
    warning('invz:sigmaCritExcluded', 'Excluded %d uniform-mode entries.', sum(~keep));
end
Sc = mean(x(keep) ./ (J0 - x(keep)));
end
```

- [ ] **Step 4: Run fast suite to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_sigma_crit.m'); assertSuccess(results)"`
Expected: PASS (1 passed, 1 filtered by assumption).

- [ ] **Step 5: Run the slow LiHoF4 test once**

Run: `INVZ_SLOW=1 matlab -batch "results = runtests('invz/tests/test_invz_sigma_crit.m'); assertSuccess(results)"`
Expected: PASS (2 passed) — this validates the entire 𝒥ν(q) machinery against 0.3004. Grid caching makes reruns fast. If it lands at ~0.30±0.03 but outside tolerance, increase the grid to 20³ before suspecting physics.

- [ ] **Step 6: Commit**

```bash
git add invz/invz_sigma_crit.m invz/tests/test_invz_sigma_crit.m
git commit -m "feat(invz): closed-form critical self-energy with fcc Watson and LiHoF4 benchmarks"
```

---

### Task 7: Scalar effective-medium K-loop (fixed Σ)

**Files:**
- Create: `invz/invz_emt_scalar.m`
- Test: `invz/tests/test_invz_emt_scalar.m`

**Interfaces:**
- Consumes: nothing invz-specific (pure numerics; test uses `invz_g`, `invz_matsubara`).
- Produces: `med = invz_emt_scalar(G0, Sigma, Jnu_flat, opts)`; `G0,Sigma` [nw×1] real (Matsubara) or complex (real axis); `Jnu_flat` [nJ×1] meV; returns `med.G, med.K` [nw×1], `med.closure`, `med.converged`, `med.iters`. Equations: G = G₀/(1+Σ+K·G₀) (HTML 28/R 9); G_q = G/(1+(𝒥−K)G) (HTML 15/R 7); K = (1/N)Σ 𝒥·G_q / G (HTML 14/R 8). `opts`: `tol` (1e-10), `max_iter` (500), `mix` (0.5), `K0` (initial K, default zeros).

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_emt_scalar
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function [G0, Jf, wn, wts, beta] = synthetic_case()
% Two-level G0 on a simple-cubic nn lattice, subcritical.
C = invz_const();  T = 1.0;
tl.Delta = 0.3;  tl.n01 = tanh(tl.Delta/(2*C.kB*T));  tl.M2 = 30;
[wn, wts, beta] = invz_matsubara(T, 60*tl.Delta);
G0 = -tl.M2*real(invz_g(tl, 1i*wn));
n = 8;  f = 2*pi*(0:n-1)/n;
[qx,qy,qz] = ndgrid(f,f,f);
J0 = 0.8/(tl.M2*2*tl.n01/tl.Delta);                 % J(0)*chi0(0) = 0.8
Jf = (J0/6)*2*(cos(qx(:))+cos(qy(:))+cos(qz(:)));
end

function test_rpa_recovery_at_sigma_zero(testCase)
% HTML Sec 5: with Sigma=0 the effective-medium G(q) equals the RPA G0/(1+J*G0) exactly.
[G0, Jf, wn, ~, ~] = synthetic_case();
med = invz_emt_scalar(G0, zeros(size(G0)), Jf, struct());
verifyTrue(testCase, med.converged);
Grpa_avg = mean(G0.'./(1 + Jf.*G0.'), 1).';        % (1/N)sum_q RPA
verifyEqual(testCase, med.G, Grpa_avg, 'RelTol', 1e-8);
verifyLessThan(testCase, med.closure, 1e-8);
% K decays like wn^-2: check last point much smaller than K(0)
verifyLessThan(testCase, abs(med.K(end)), 0.05*max(abs(med.K(1)), eps));
end

function test_dyson_shift_with_constant_sigma(testCase)
% A constant Sigma acts exactly like G0 -> G0/(1+Sigma) (HTML eq 30).
[G0, Jf, ~, ~, ~] = synthetic_case();
s0 = 0.2;
medS = invz_emt_scalar(G0, s0*ones(size(G0)), Jf, struct());
medR = invz_emt_scalar(G0/(1+s0), zeros(size(G0)), Jf, struct());
verifyEqual(testCase, medS.G, medR.G, 'RelTol', 1e-8);
verifyEqual(testCase, medS.K, medR.K, 'RelTol', 1e-6);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_emt_scalar.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_emt_scalar'".

- [ ] **Step 3: Implement `invz/invz_emt_scalar.m`**

```matlab
function med = invz_emt_scalar(G0, Sigma, Jnu_flat, opts)
%INVZ_EMT_SCALAR Scalar effective-medium fixed point at fixed self-energy.
% G  = G0./(1 + Sigma + K.*G0)            (Dyson, R eq 9)
% Gq = G ./ (1 + (J_nu - K).*G)           (R eq 7)
% K  = (1/N) sum_{q,nu} J_nu.*Gq ./ G     (R eq 8)
if nargin < 4, opts = struct(); end
tol  = 1e-10; if isfield(opts,'tol'),      tol  = opts.tol;      end
mit  = 500;   if isfield(opts,'max_iter'), mit  = opts.max_iter; end
mix  = 0.5;   if isfield(opts,'mix'),      mix  = opts.mix;      end
G0 = G0(:);  Sigma = Sigma(:);  Jf = Jnu_flat(:);
K = zeros(size(G0)); if isfield(opts,'K0') && ~isempty(opts.K0), K = opts.K0(:); end
nJ = numel(Jf);
converged = false;
for it = 1:mit
    G  = G0 ./ (1 + Sigma + K.*G0);
    Gq = G.' ./ (1 + (Jf - K.').*G.');           % [nJ, nw]
    Kn = (Jf.'*Gq).'/nJ ./ G;
    res = max(abs(Kn - K) ./ max(abs(Kn), 1e-14));
    K = K + mix*(Kn - K);
    if res < tol, converged = true; break; end
end
G  = G0 ./ (1 + Sigma + K.*G0);
Gq = G.' ./ (1 + (Jf - K.').*G.');
med.G = G;  med.K = K;
med.closure = max(abs(mean(Gq,1).' - G) ./ max(abs(G), 1e-14));
med.converged = converged;  med.iters = it;
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_emt_scalar.m'); assertSuccess(results)"`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
git add invz/invz_emt_scalar.m invz/tests/test_invz_emt_scalar.m
git commit -m "feat(invz): scalar Matsubara effective-medium K loop with RPA-recovery test"
```

---

### Task 8: λ_p and Σ = α + γ·g, verified by the sum-rule cancellation

**Files:**
- Create: `invz/invz_lambdas.m`
- Create: `invz/invz_sigma.m`
- Test: `invz/tests/test_invz_sigma.m`

**Interfaces:**
- Consumes: `invz_matsubara`, `invz_g`, two-level struct `tl`.
- Produces: `lam = invz_lambdas(K, g, wts, beta, plist)` (column over plist); `sig = invz_sigma(tl, lam, K, g, beta)` with `alpha` (HTML eq 27), `gamma` [nw×1] (HTML eq 23), `Sigma = alpha + gamma.*g`. `lam(1)`,`lam(2)` must correspond to p=1,2.

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_sigma
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function [tl, wn, wts, beta, g, K] = fixture(Ecut_factor)
C = invz_const();  T = 1.0;
tl.Delta = 0.3;  tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.M2 = 30;  tl.m = 0;  tl.g0 = 2*tl.n01/tl.Delta;
[wn, wts, beta] = invz_matsubara(T, Ecut_factor*tl.Delta);
g = real(invz_g(tl, 1i*wn));
K = 2e-3 ./ (1 + (wn/0.5).^2);          % synthetic even K decaying like wn^-2
end

function test_lambda_analytic(testCase)
% With K = const = c: lambda_1 = c*(1/beta)*sum g = c;  lambda_2 = c*(1/2)[g0+beta(1-n01^2)].
[tl, ~, wts, beta, g, ~] = fixture(400);
c = 3.7e-3;
lam = invz_lambdas(c*ones(size(g)), g, wts, beta, [1 2]);
verifyEqual(testCase, lam(1), c, 'RelTol', 1e-3);
verifyEqual(testCase, lam(2), c*0.5*(tl.g0 + beta*(1-tl.n01^2)), 'RelTol', 1e-6);
end

function test_sumrule_cancellation(testCase)
% THE identity that fixes alpha (HTML eq 24): (1/beta) sum_n G0.*(K.*G0 + Sigma) = 0.
[tl, ~, wts, beta, g, K] = fixture(400);
G0 = -tl.M2*g;
lam = invz_lambdas(K, g, wts, beta, [1 2]);
sig = invz_sigma(tl, lam, K, g, beta);
lhs   = sum(wts .* G0 .* (K.*G0 + sig.Sigma)) / beta;
scale = sum(wts .* abs(G0 .* K .* G0)) / beta;
verifyLessThan(testCase, abs(lhs)/scale, 1e-3);
% convergence order: residual falls when the grid doubles
[tl2, ~, wts2, beta2, g2, K2] = fixture(800);
G02 = -tl2.M2*g2;
lam2 = invz_lambdas(K2, g2, wts2, beta2, [1 2]);
sig2 = invz_sigma(tl2, lam2, K2, g2, beta2);
lhs2 = sum(wts2 .* G02 .* (K2.*G02 + sig2.Sigma)) / beta2;
verifyLessThan(testCase, abs(lhs2), abs(lhs));
end

function test_sigma_zero_when_K_zero(testCase)
[tl, ~, wts, beta, g, ~] = fixture(400);
lam = invz_lambdas(zeros(size(g)), g, wts, beta, [1 2]);
sig = invz_sigma(tl, lam, zeros(size(g)), g, beta);
verifyEqual(testCase, sig.alpha, 0, 'AbsTol', 1e-15);
verifyLessThan(testCase, max(abs(sig.Sigma)), 1e-15);
end

function test_matches_existing_src_formula(testCase)
% Cross-check against the already unit-tested src/emt_compute_x_from_lambdas 'jensen_216_219'.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','..','src'));
[tl, ~, wts, beta, g, K] = fixture(400);
lam = invz_lambdas(K, g, wts, beta, [1 2]);
sig = invz_sigma(tl, lam, K, g, beta);
tp = struct('model','jensen_216_219','beta',beta,'M2',tl.M2, ...
            'n0',(1+tl.n01)/2,'n1',(1-tl.n01)/2,'clamp_scale',1e9);
[X, ~] = emt_compute_x_from_lambdas(g, lam, tp, K);
verifyEqual(testCase, sig.Sigma, X(:), 'RelTol', 1e-9);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_sigma.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_lambdas'".

- [ ] **Step 3: Implement the two functions**

`invz/invz_lambdas.m`:
```matlab
function lam = invz_lambdas(K, g, wts, beta, plist)
%INVZ_LAMBDAS lambda_p = (1/beta) sum_n K(iwn) g(iwn)^p   (HTML eq 21 / J 2.19).
K = K(:); g = g(:); wts = wts(:);
lam = zeros(numel(plist),1);
for ip = 1:numel(plist)
    lam(ip) = real(sum(wts .* K .* g.^plist(ip))) / beta;
end
end
```

`invz/invz_sigma.m`:
```matlab
function sig = invz_sigma(tl, lam, K, g, beta)
%INVZ_SIGMA Jensen self-energy Sigma = alpha + gamma.*g  (HTML eqs 22, 23, 27).
pref = tl.M2 / tl.n01^2;
sig.alpha = pref * (lam(2) - 0.5*(tl.g0 + beta*(1 - tl.n01^2))*lam(1));
sig.gamma = pref * (lam(1) - (1 - tl.n01^2)*K(:));
sig.Sigma = sig.alpha + sig.gamma .* g(:);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_sigma.m'); assertSuccess(results)"`
Expected: PASS (4 passed). If `test_matches_existing_src_formula` fails on interface details of `emt_compute_x_from_lambdas` (e.g. it expects `n01` directly, or `g` as row), adapt the *test call*, not `invz_sigma` — read `src/emt_compute_x_from_lambdas.m` lines 60–130 and `tests/test_tracka_moments.m` for the exact contract.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_lambdas.m invz/invz_sigma.m invz/tests/test_invz_sigma.m
git commit -m "feat(invz): lambda_p and Jensen self-energy with sum-rule cancellation test"
```

---

### Task 9: Self-consistent solve at one (T, Bx) point

**Files:**
- Create: `invz/invz_solve_point.m`
- Test: `invz/tests/test_invz_solve_point.m`

**Interfaces:**
- Consumes: everything from Tasks 2–8.
- Produces: `pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts)` with fields `Sigma0, alpha, lambda [2×1], K [nw×1], G [nw×1], Sigma [nw×1], tl, chi0cc0, crit, sumrule_rel, converged, outer_iters`. `opts`: `Ecut` (default 40 meV), `hyp` (default true), `J0eff` (default `ion.J0eff`), `mix_outer` (0.7), `tol_outer` (1e-8), `max_outer` (60), `emt` (opts struct passed to `invz_emt_scalar`). Definitions: `crit = 1 + Sigma0 - J0eff*chi0cc0` (paramagnet stable while > 0); `sumrule_rel = |(1/β)Σ wts·G − (−⟨J̃z²⟩)| / ⟨J̃z²⟩`.

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_solve_point
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function Jf = toy_couplings(J0)
% Cheap stand-in lattice (sc-like) with uniform-mode coupling J0, so the test
% does not depend on slow dipole sums.
n = 10;  f = 2*pi*(0:n-1)/n;
[qx,qy,qz] = ndgrid(f,f,f);
Jf = (J0/6)*2*(cos(qx(:))+cos(qy(:))+cos(qz(:)));
Jf = Jf(2:end);                                   % drop uniform mode (measure zero)
end

function test_converges_and_is_physical(testCase)
ion = invz_ion();
Jf = toy_couplings(ion.J0eff);
pt = invz_solve_point(ion, 1.0, 4.0, Jf, struct('hyp', true));
verifyTrue(testCase, pt.converged);
verifyGreaterThan(testCase, pt.Sigma0, 0);        % fluctuations suppress the response
verifyLessThan(testCase, pt.Sigma0, 1.0);
verifyGreaterThan(testCase, pt.crit, 0);          % paramagnetic at 1 K, 4 T
verifyLessThan(testCase, pt.sumrule_rel, 0.10);   % Jensen: resummed G obeys sum rule approximately
end

function test_sigma_grows_toward_criticality(testCase)
% Cooling at fixed field: K(0) grows, so Sigma(0) must grow monotonically.
ion = invz_ion();
Jf = toy_couplings(ion.J0eff);
s = zeros(1,3);  Ts = [2.0 1.0 0.5];
for k = 1:3
    pt = invz_solve_point(ion, Ts(k), 4.0, Jf, struct('hyp', true));
    s(k) = pt.Sigma0;
end
verifyGreaterThan(testCase, s(2), s(1));
verifyGreaterThan(testCase, s(3), s(2));
end

function test_interaction_off_gives_zero_sigma(testCase)
ion = invz_ion();
pt = invz_solve_point(ion, 1.0, 4.0, 1e-12*ones(100,1), struct('hyp', false));
verifyLessThan(testCase, abs(pt.Sigma0), 1e-8);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_solve_point.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_solve_point'".

- [ ] **Step 3: Implement `invz/invz_solve_point.m`**

```matlab
function pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_POINT Self-consistent 1/z solution at one paramagnetic (T, Bx) point.
% Outer loop: EMT K (Task 7) <-> lambda_p (Task 8) <-> Sigma, at fixed single-ion input.
if nargin < 5, opts = struct(); end
Ecut  = 40;   if isfield(opts,'Ecut'),      Ecut  = opts.Ecut;      end
hyp   = true; if isfield(opts,'hyp'),       hyp   = opts.hyp;       end
J0eff = ion.J0eff; if isfield(opts,'J0eff'), J0eff = opts.J0eff;    end
mixo  = 0.7;  if isfield(opts,'mix_outer'), mixo  = opts.mix_outer; end
tolo  = 1e-8; if isfield(opts,'tol_outer'), tolo  = opts.tol_outer; end
maxo  = 60;   if isfield(opts,'max_outer'), maxo  = opts.max_outer; end
eopts = struct(); if isfield(opts,'emt'), eopts = opts.emt; end

[wn, wts, beta] = invz_matsubara(T, Ecut);
si  = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', hyp));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));                 % full (electro)nuclear cc Green function
tl  = invz_twolevel(ion, T, Bx);                 % electronic two-level params for Sigma
g   = real(invz_g(tl, 1i*wn));

Sigma = zeros(size(wn));  K = zeros(size(wn));
converged = false;
for outer = 1:maxo
    eopts.K0 = K;
    med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
    K   = med.K;
    lam = invz_lambdas(K, g, wts, beta, [1 2]);
    sg  = invz_sigma(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + mixo*(sg.Sigma - Sigma);
    if dS < tolo, converged = true; break; end
end
pt.Sigma0 = Sigma(1);  pt.alpha = sg.alpha;  pt.lambda = lam;
pt.K = K;  pt.G = med.G;  pt.Sigma = Sigma;  pt.tl = tl;
pt.chi0cc0 = -G0(1);
pt.crit = 1 + pt.Sigma0 - J0eff*pt.chi0cc0;
pt.sumrule_rel = abs(sum(wts.*med.G)/beta + si.JzJz_fluct) / si.JzJz_fluct;
pt.converged = converged && med.converged;
pt.outer_iters = outer;
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_solve_point.m'); assertSuccess(results)"`
Expected: PASS (3 passed). Typical runtime: seconds per point (nw ≈ 100–500, nJ = 999).

- [ ] **Step 5: Commit**

```bash
git add invz/invz_solve_point.m invz/tests/test_invz_solve_point.m
git commit -m "feat(invz): self-consistent 1/z solve at a (T,Bx) point with physical sanity tests"
```

---

### Task 10: Phase boundary against published numbers

**Files:**
- Create: `invz/invz_critical.m`
- Create: `invz/invz_critical_T0field.m`
- Test: `invz/tests/test_invz_critical.m`

**Interfaces:**
- Consumes: `invz_solve_point`, `invz_sigma_crit`, `invz_chi0z`, `invz_single_ion`, `invz_jq_modes`, `qVec_generator`.
- Produces: `bx = invz_critical(ion, T, Jnu_flat, opts)` — bisection of `pt.crit` over Bx in `opts.window` (default [2 7] T, tol 0.02 T); `Tc = invz_critical_T0field(ion, Sc, J0eff)` — solves 1+Sc = J0eff·χ₀cc(0;T,Bx=0) by bisection over T∈[0.8, 3] K (full electronuclear χ₀ with elastic term; this is the Bx=0 route where the closed-form Σc replaces the α+γg machinery, HTML Sec 9 caveat).

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_critical
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function [Jf, J0] = lihof4_couplings()
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'size', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Jf = Jnu(:);  J0 = info.Jcc0;
end

function test_zero_field_tc(testCase)
% R 2007: Tc(B=0) = 1.74 K in the 1/z theory (experiment 1.53 K; MF higher). SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
Sc = invz_sigma_crit(J0, Jf);
Tc  = invz_critical_T0field(ion, Sc, J0);
TcMF = invz_critical_T0field(ion, 0,  J0);
verifyEqual(testCase, Tc, 1.74, 'AbsTol', 0.08);
verifyGreaterThan(testCase, TcMF, Tc);            % fluctuations suppress Tc
end

function test_critical_field_at_310mK(testCase)
% R 2007: Hc(0.31 K) = 42.4-43.0 kOe with Sigma_c(0) = 0.0932. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[Jf, J0] = lihof4_couplings();
bx = invz_critical(ion, 0.31, Jf, struct('J0eff', J0));
verifyGreaterThan(testCase, bx, 4.0);  verifyLessThan(testCase, bx, 4.6);   % Tesla
pt = invz_solve_point(ion, 0.31, bx, Jf, struct('hyp', true, 'J0eff', J0));
verifyEqual(testCase, pt.Sigma0, 0.0932, 'AbsTol', 0.02);
end
```

- [ ] **Step 2: Run fast suite to verify it fails cleanly**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_critical.m'); assertSuccess(results)"`
Expected: both tests filtered (assumption), suite passes trivially — then with `INVZ_SLOW=1` it must FAIL with "Unrecognized function or variable 'invz_critical_T0field'".

- [ ] **Step 3: Implement the two functions**

`invz/invz_critical_T0field.m`:
```matlab
function Tc = invz_critical_T0field(ion, Sc, J0eff)
%INVZ_CRITICAL_T0FIELD Zero-field Tc from 1 + Sigma_c = J(0)*chi0_cc(0;T).
% Sc is the closed-form critical self-energy (invz_sigma_crit); Sc=0 gives the MF-RPA Tc.
f = @(T) J0eff*static_chi_cc(ion, T) - (1 + Sc);
Tlo = 0.8;  Thi = 3.0;
assert(f(Tlo) > 0 && f(Thi) < 0, 'invz:bracket', 'Tc not bracketed in [0.8, 3.0] K');
for it = 1:60
    Tm = 0.5*(Tlo + Thi);
    if f(Tm) > 0, Tlo = Tm; else, Thi = Tm; end
end
Tc = 0.5*(Tlo + Thi);
end

function c = static_chi_cc(ion, T)
si = invz_single_ion(ion, T, [0 0 0], struct('hyp', true));
cc = real(invz_chi0z(si, T, 0, struct('elastic', true)));
c  = cc(3,3);
end
```

`invz/invz_critical.m`:
```matlab
function bx = invz_critical(ion, T, Jnu_flat, opts)
%INVZ_CRITICAL Critical transverse field at temperature T (paramagnetic side).
% Bisection on pt.crit = 1 + Sigma(0) - J(0)*chi0_cc(0), which is positive in the
% paramagnet and crosses zero at the boundary.
if nargin < 4, opts = struct(); end
win = [2 7];  if isfield(opts,'window'), win = opts.window; end
tol = 0.02;   if isfield(opts,'tol'),    tol = opts.tol;    end
f = @(B) crit_of(ion, T, B, Jnu_flat, opts);
flo = f(win(1));  fhi = f(win(2));
assert(flo < 0 && fhi > 0, 'invz:bracket', ...
    'Boundary not bracketed: crit(%.2fT)=%.3g, crit(%.2fT)=%.3g', win(1), flo, win(2), fhi);
lo = win(1);  hi = win(2);
while hi - lo > tol
    mid = 0.5*(lo + hi);
    if f(mid) < 0, lo = mid; else, hi = mid; end
end
bx = 0.5*(lo + hi);
end

function c = crit_of(ion, T, B, Jf, opts)
pt = invz_solve_point(ion, T, B, Jf, opts);
c = pt.crit;
end
```

- [ ] **Step 4: Run the slow suite**

Run: `INVZ_SLOW=1 matlab -batch "results = runtests('invz/tests/test_invz_critical.m'); assertSuccess(results)"`
Expected: PASS (2 passed). Runtime: minutes (bisection ≈ 8 solve-point calls; the Jq grid comes from cache). Two known failure modes to check in order: (i) Tc lands near the MF value → Σ not being applied (check `pt.Sigma0` > 0 at the boundary); (ii) Hc off by >0.5 T → `J0eff` inconsistent between `invz_jq_modes` info and `invz_solve_point` opts — always pass `info.Jcc0` explicitly as in the test.

- [ ] **Step 5: Commit**

```bash
git add invz/invz_critical.m invz/invz_critical_T0field.m invz/tests/test_invz_critical.m
git commit -m "feat(invz): phase boundary solvers benchmarked against Ronnow 2007 (Tc=1.74K, Hc(0.31K))"
```

---

### Task 11: Real-axis χ(q,ω), soft mode, and driver scripts

**Files:**
- Create: `invz/invz_chi_realaxis.m`
- Create: `invz/invz_run_phase_diagram.m` (script)
- Create: `invz/invz_run_spectra.m` (script)
- Create: `invz/README.md`
- Test: `invz/tests/test_invz_chi_observable.m`

**Interfaces:**
- Consumes: converged `pt` from `invz_solve_point` (uses `pt.alpha`, `pt.lambda(1)`, `pt.tl`), `invz_chi0z`, `invz_emt_scalar`.
- Produces: `out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)` where `w` is a real frequency grid (meV): `out.Sigma_w` (complex, Σ(ω) = α + (M²/n₀₁²)[λ₁ − (1−n₀₁²)K(ω)]·g(ω), with K(ω) from `opts.npass` (default 3) alternating real-axis EMT passes at `z = w + 1i*eta`, `eta` default 5e-3 meV); `out.chi_cc_q` [nJ_sel × nw] mode-resolved χcc(q,ν,ω) = χ̃₀/(1−𝒥ν·χ̃₀) with χ̃₀(ω) = −G₀(ω)/(1+Σ(ω)); `out.chi0cc_w`.

- [ ] **Step 1: Write the failing test**

```matlab
function tests = test_invz_chi_observable
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function test_rpa_limit_and_positivity(testCase)
% With Sigma forced to zero the output must equal scalar RPA built directly from chi0.
ion = invz_ion();
T = 1.0;  Bx = 4.0;  w = (0.02:0.02:0.8).';  eta = 5e-3;
Jsel = [5e-3; 6.4e-3];                       % two representative couplings (meV)
pt = struct('alpha', 0, 'lambda', [0;0], 'tl', invz_twolevel(ion, T, Bx));
out = invz_chi_realaxis(ion, T, Bx, pt, w, struct('Jsel', Jsel, 'eta', eta, 'npass', 1));
si  = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));
c0  = squeeze(invz_chi0z(si, T, w + 1i*eta, struct('elastic', false)));
c0cc = squeeze(c0(3,3,:));
rpa = c0cc ./ (1 - Jsel(2)*c0cc);
verifyEqual(testCase, out.chi_cc_q(2,:).', rpa, 'RelTol', 1e-6);
verifyGreaterThan(testCase, min(imag(out.chi_cc_q(2,:))), -1e-10);  % chi'' >= 0 for w > 0
end

function test_soft_mode_near_criticality(testCase)
% R 2007 Fig 2: at T=0.31 K, H=Hc the lowest mode bottoms at ~0.19 meV (calc), never zero. SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'size', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
T = 0.31;
bx = invz_critical(ion, T, Jnu(:), struct('J0eff', info.Jcc0));
pt = invz_solve_point(ion, T, bx, Jnu(:), struct('hyp', true, 'J0eff', info.Jcc0));
w  = (0.01:0.005:0.6).';
out = invz_chi_realaxis(ion, T, bx, pt, w, struct('Jsel', info.Jcc0, 'eta', 5e-3));
[~, ipk] = max(imag(out.chi_cc_q(1,:)));
Epk = w(ipk);
verifyGreaterThan(testCase, Epk, 0.10);  verifyLessThan(testCase, Epk, 0.28);
end
```

- [ ] **Step 2: Run fast test to verify it fails**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_chi_observable.m'); assertSuccess(results)"`
Expected: FAIL — "Unrecognized function or variable 'invz_chi_realaxis'".

- [ ] **Step 3: Implement `invz/invz_chi_realaxis.m`**

```matlab
function out = invz_chi_realaxis(ion, T, Bx, pt, w, opts)
%INVZ_CHI_REALAXIS 1/z-renormalized cc susceptibility on the real axis.
% Sigma(w) = alpha + (M2/n01^2)*[lambda1 - (1-n01^2)*K(w)] * g(w)   (HTML eqs 22-23),
% alpha and lambda1 fixed by the converged Matsubara solve (pt).
% chi~0(w) = chi0_cc(w)/(1+Sigma(w));  chi(q,nu,w) = chi~0/(1 - J_nu*chi~0)   (HTML eq 29-30).
if nargin < 6, opts = struct(); end
eta   = 5e-3; if isfield(opts,'eta'),   eta   = opts.eta;   end
npass = 3;    if isfield(opts,'npass'), npass = opts.npass; end
Jsel  = ion.J0eff; if isfield(opts,'Jsel'), Jsel = opts.Jsel; end
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));
z  = w(:) + 1i*eta;
c0 = invz_chi0z(si, T, z, struct('elastic', false));
G0 = -squeeze(c0(3,3,:));
tl = pt.tl;
g  = invz_g(tl, z);
pref = tl.M2/tl.n01^2;
if isfield(pt,'K') && ~isempty(pt.K)
    Kw = pt.K(1)*ones(size(z));                    % seed with static Matsubara K
else
    Kw = zeros(size(z));                           % RPA-limit callers may omit K
end
for pass = 1:npass
    Sw  = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw).*g;
    if isfield(opts,'Jfull') && ~isempty(opts.Jfull)
        med = invz_emt_scalar(G0, Sw, opts.Jfull, struct('max_iter', 100, 'tol', 1e-8));
        Kw = med.K;
    else
        break;                                     % no lattice pass requested: single shot
    end
end
Sw = pt.alpha + pref*(pt.lambda(1) - (1 - tl.n01^2)*Kw).*g;
chit = (-G0) ./ (1 + Sw);
out.Sigma_w  = Sw;
out.chi0cc_w = -G0;
out.chi_cc_q = zeros(numel(Jsel), numel(z));
for k = 1:numel(Jsel)
    out.chi_cc_q(k,:) = (chit ./ (1 - Jsel(k)*chit)).';
end
end
```

- [ ] **Step 4: Run fast test, then slow test**

Run: `matlab -batch "results = runtests('invz/tests/test_invz_chi_observable.m'); assertSuccess(results)"`
Expected: PASS (1 passed, 1 filtered). Then `INVZ_SLOW=1 matlab -batch "results = runtests('invz/tests/test_invz_chi_observable.m'); assertSuccess(results)"` — expected PASS (2 passed).

- [ ] **Step 5: Write the driver scripts and README**

`invz/invz_run_phase_diagram.m`:
```matlab
%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'size', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Ts = [0.15 0.31 0.6 0.9 1.2 1.4 1.6];
Bc = nan(size(Ts));
for k = 1:numel(Ts)
    try
        Bc(k) = invz_critical(ion, Ts(k), Jnu(:), struct('J0eff', info.Jcc0));
    catch err
        warning('T=%.2f K: %s', Ts(k), err.message);
    end
    fprintf('T = %.2f K  ->  Bc = %.2f T\n', Ts(k), Bc(k));
end
Tc0 = invz_critical_T0field(ion, invz_sigma_crit(info.Jcc0, Jnu(:)), info.Jcc0);
fprintf('Zero-field Tc (1/z) = %.3f K  [target 1.74 K]\n', Tc0);
figure; plot(Ts, Bc*10, 'o-'); hold on; plot(Tc0, 0, 'ks');
xlabel('T (K)'); ylabel('B_c (kOe)'); title('LiHoF_4 1/z phase boundary (paramagnetic side)');
```

`invz/invz_run_spectra.m`:
```matlab
%INVZ_RUN_SPECTRA chi''_cc(omega) at the uniform mode vs field at T = 0.31 K (cf. R 2007 Fig 2,
% Kovacevic 2016 Fig 3d). RPA vs 1/z overlay.
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'size', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
T = 0.31;  w = (0.01:0.005:0.7).';
fields = [3.6 4.2 4.8 5.4 6.0];
figure; hold on;
for B = fields
    pt = invz_solve_point(ion, T, B, Jnu(:), struct('hyp', true, 'J0eff', info.Jcc0));
    out  = invz_chi_realaxis(ion, T, B, pt, w, struct('Jsel', info.Jcc0));
    pt0  = struct('alpha', 0, 'lambda', [0;0], 'tl', pt.tl, 'K', 0*pt.K);
    out0 = invz_chi_realaxis(ion, T, B, pt0, w, struct('Jsel', info.Jcc0, 'npass', 1));
    plot(w, imag(out.chi_cc_q(1,:)), '-', 'DisplayName', sprintf('1/z, %.1f T', B));
    plot(w, imag(out0.chi_cc_q(1,:)), '--', 'DisplayName', sprintf('RPA, %.1f T', B));
    fprintf('B = %.1f T : Sigma(0) = %.4f, 1+Sigma(0) = %.4f\n', B, pt.Sigma0, 1+pt.Sigma0);
end
xlabel('\omega (meV)'); ylabel('\chi''''_{cc}'); legend show;
```

`invz/README.md`: one page — module purpose, layer diagram, the interface contract table from this plan, how to run fast/slow tests, the published benchmark numbers with their source equations, and the scope statement (paramagnetic phase only; ordered phase / α_m / ξ / H_MF and the 4×4 tensor RPA observable layer are future work; small-Bx region where Δ ≲ hyperfine width is expected to overestimate the Tc suppression, per R 2007's own comparison with QMC).

- [ ] **Step 6: Run both drivers manually and sanity-check output**

Run: `matlab -batch "invz_run_phase_diagram"` then `matlab -batch "invz_run_spectra"`.
Expected: boundary curve decreasing from ≈5.3 T at 0.15 K to 0 near 1.74 K; spectra show the 1/z curves shifted/damped relative to RPA, with the soft-mode minimum near the critical field at ≈0.19 meV. Paste the two printed summaries into the commit message body.

- [ ] **Step 7: Commit**

```bash
git add invz/invz_chi_realaxis.m invz/invz_run_phase_diagram.m invz/invz_run_spectra.m invz/README.md invz/tests/test_invz_chi_observable.m
git commit -m "feat(invz): real-axis 1/z susceptibility, soft-mode benchmark, drivers and README"
```

---

## Non-goals (explicitly out of scope for this plan)

- Ordered phase and zero-field low-T region: the elastic-sector self-energy (HTML eqs 37–40: α_m, ξ, λ₃) and the modified mean field H_MF (HTML eqs 41–43). The Bx=0 boundary point is handled only via the closed-form Σc route.
- Full 3×3 tensor / 4×4 sublattice RPA observable layer (neutron cross-sections, transverse components) — the existing `MF_RPA_Yikai.m` remains the tool for that; hook point documented in README.
- Free energy, δU, heat capacity (HTML eqs 34–35, 44).
- Ewald acceleration of the dipole sums (brute-force + cache is sufficient at 16³).
- Any modification of `src/` (the existing EMT solver stays as is; `invz` only *reads* `emt_compute_x_from_lambdas` in one cross-check test).

## Self-Review Notes

- Spec coverage: single-ion (T2), χ₀ complex-z (T3), two-level+grid (T4), 𝒥ν(q) (T5), closed-form Σc + both published lattice benchmarks (T6), EMT core (T7), λ/Σ + sum-rule cancellation (T8), self-consistent point (T9), phase boundary vs 1.74 K / 42.4–43 kOe / Σc=0.0932 (T10), real-axis spectra vs 0.19 meV soft mode + hyperfine ratio 0.77 (T3/T11). All five gaps identified in the code-map synthesis are closed: Matsubara path (T3/T4), Δ/M²/n₀₁ extraction (T4), sublattice branches (T5), wiring (T9), convention fix (Global Constraints).
- Type consistency: `tl` fields (`Delta,M2,m,n01,g0`) used identically in T4/T8/T9/T11; `pt` fields defined in T9 and consumed in T10/T11; `info.Jcc0` produced in T5, consumed in T6/T10/T11; `wts/beta` conventions identical throughout.
- Known-risk checkpoints are explicit steps, not hidden assumptions: Lorentz-term convention (T5 step 4), B64s sign (T2 step 5), `emt_compute_x_from_lambdas` call contract (T8 step 4).
