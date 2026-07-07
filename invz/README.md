# `invz` — 1/z effective-medium expansion for LiHoF4

## Purpose

This module implements the Jensen 1/z effective-medium (EMT) expansion for the
dipolar-coupled, hyperfine-split Ising ferromagnet LiHoF4, following Rønnow
*et al.*, PRB **75**, 054426 (2007) ("R 2007"). It computes the paramagnetic-phase
transverse-field phase boundary and the real-frequency electronic + nuclear
susceptibility χ_cc(q,ν,ω), going beyond mean-field/RPA by resumming the leading
1/z fluctuation correction to the single-ion self-energy self-consistently with
a scalar effective-medium lattice sum. The full derivation (all numbered
equations referenced in the code as "HTML eq N") lives in
`jensen_1z_framework.html` at the repository root.

Everything here is **paramagnetic-phase, scalar (cc-channel) only**; see
Scope below.

## Layer diagram

Each layer consumes only the outputs of the layers above it (T# = task number
in the implementation plan):

```
T1  invz_ion, invz_const, stevens_ops              ion parameters, CF/Stevens operators
T2  invz_single_ion                                 exact diagonalization + transverse MF (si)
T3  invz_chi0z                                      single-ion chi0(mu,nu; z) at any complex z
T4  invz_matsubara, invz_twolevel, invz_g            Matsubara grid; two-level (Delta,M2,n01,g0); g(z)
T5  invz_jq_modes                                   dipole+exchange coupling branches J_nu(q) (4x4 cc sublattice eigs)
T6  invz_sigma_crit, invz_critical_T0field           closed-form critical self-energy; zero-field Tc
T7  invz_emt_scalar                                 scalar EMT K-loop (Dyson G, cavity Gq, self-consistent K)
T8  invz_lambdas, invz_sigma                        lambda_p moments; Sigma = alpha + gamma*g
T9  invz_solve_point                                self-consistent (T,Bx) solve: outer loop T7<->T8
T10 invz_critical                                   bisection for Bc(T) using pt.crit sign
T11 invz_chi_realaxis                               real-axis analytic continuation -> chi(q,nu,w)
    invz_run_phase_diagram.m, invz_run_spectra.m     interactive driver scripts (this task)
```

Dependency direction is strictly top-to-bottom; T11 is the only consumer of
`pt` (the struct produced by T9/`invz_solve_point`).

## Interface contract table

| Function | Signature | Produces |
|---|---|---|
| `invz_ion` | `ion = invz_ion()` | LiHoF4 CF/exchange/lattice parameters (R 2007 Table I) |
| `invz_single_ion` | `si = invz_single_ion(ion,T,B,opts)` | `si.E,V,P,Mx,My,Mz,Jexp,hx,JzJz_fluct` |
| `invz_chi0z` | `chi = invz_chi0z(si,T,z,opts)` | `[3,3,nz]` single-ion tensor susceptibility at complex `z` |
| `invz_twolevel` | `tl = invz_twolevel(ion,T,Bx)` | `tl.Delta,M2,m,n01,g0` (doublet two-level params) |
| `invz_g` | `g = invz_g(tl,z)` | two-level response `g(z)` |
| `invz_jq_modes` | `[Jnu,info] = invz_jq_modes(ion,qvec,opts)` | `Jnu` `[nq,4]` cc branch eigenvalues; `info.Jcc0` uniform-mode coupling |
| `invz_sigma_crit` | `Sc = invz_sigma_crit(J0,Jnu_flat)` | closed-form critical self-energy Σc(0) |
| `invz_critical_T0field` | `Tc = invz_critical_T0field(ion,Sc,J0eff)` | zero-field Tc from `1+Sc = J(0)*chi0_cc(0;T)` |
| `invz_emt_scalar` | `med = invz_emt_scalar(G0,Sigma,Jnu_flat,opts)` | `med.G,K,closure,converged` |
| `invz_lambdas` | `lam = invz_lambdas(K,g,wts,beta,plist)` | moments λ_p |
| `invz_sigma` | `sig = invz_sigma(tl,lam,K,g,beta)` | `sig.alpha,gamma,Sigma` |
| `invz_solve_point` | `pt = invz_solve_point(ion,T,Bx,Jnu_flat,opts)` | `pt.alpha,lambda[2x1],K[nw],Sigma,tl,Sigma0,crit,converged` |
| `invz_critical` | `bx = invz_critical(ion,T,Jnu_flat,opts)` | critical transverse field Bc(T), bisection on `pt.crit` |
| `invz_chi_realaxis` | `out = invz_chi_realaxis(ion,T,Bx,pt,w,opts)` | `out.Sigma_w,chi0cc_w,chi_cc_q[nJsel,nw]` real-axis spectra |

## Prerequisites

Repo-root `MF_dipole.m`, `exchange.m`, `qVec_generator.m` (tracked). `src/` is
optional — one cross-check test (`test_matches_existing_src_formula` in
`invz/tests/test_invz_sigma.m`) auto-skips (via `assumeTrue`) if it is absent.

## Running the tests

Fast suite (seconds; the default/CI gate):

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

Full suite including the slow, physics-benchmark tests (Jq caches for 12/16/24³
grids at `dpRng=30` are pre-warmed in `invz/cache/` in this working tree, so
even the "slow" tests typically run in well under a minute once warm;
`invz/cache/` is gitignored, so a fresh clone starts cold — cold-cache runs of
the dipole-sum-heavy tests can then take ~10-15 min on first run):

```bash
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz/tests'); assertSuccess(results)"
```

`qVec_generator` prints lattice diagnostics to stdout on every call; inside the
test suite these calls are wrapped in `evalc(...)` to keep test output clean.
The two driver scripts below call it directly (unwrapped) since their
printout is expected/useful for an interactive tool.

## Driver scripts

- `invz_run_phase_diagram.m` — reproduces the R 2007 Fig. 1 paramagnetic-side
  boundary Bc(T) plus the zero-field Tc. **Full run: ~1-2 hours** (bisection
  over full lattice EMT solves at each of 7 temperatures). Prints one progress
  line per temperature as it runs.
- `invz_run_spectra.m` — χ''_cc(ω) at the uniform mode vs field at T=0.31 K,
  1/z vs bare-RPA overlay (cf. R 2007 Fig. 2, Kovačević 2016 Fig. 3d).

## Published benchmarks reproduced by this module

| Quantity | Published | This module | Source |
|---|---|---|---|
| Σc, fcc nearest-neighbour (Watson integral) | 0.3447 | 0.3447 (Richardson 40³/80³) | R 2007, text below eq. (10) |
| Σc(0), LiHoF4, H=0 | 0.3004 | 0.2980 (Richardson 12³/24³, `dpRng=30`) | R 2007 eq. (10) |
| 𝒥_D·D_cc(0) | 6.821 μeV | 6.821-6.824 μeV (`dpRng`-converged) | R 2007 eq. (4) |
| 𝒥_D·D_aa(0) | 3.912 μeV | 3.910-3.912 μeV (`dpRng`-converged) | R 2007 eq. (4) |
| Tc(H=0) | 1.74 K | ≈1.74 K (`invz_critical_T0field`, AbsTol 0.08) | R 2007 |
| Hc(0.31 K) | 42.4-43 kOe | within [4.0,4.6] T bisection window | R 2007 Fig. 1/2 |
| Σ(0) at Hc(0.31 K) | 0.0932 | AbsTol 0.02 | R 2007 |
| Soft-mode energy at H≈Hc, T=0.31 K | ≈0.19 meV (calc.) | computed 0.22 meV (within the [0.10, 0.28] acceptance band; published calculation ≈0.19 meV, R 2007 Fig. 2) | R 2007 Fig. 2 |

## Scope

**In scope:** the paramagnetic phase only — single-ion response, the scalar
(cc-channel) 1/z self-energy and effective-medium lattice sum, the
paramagnetic-side phase boundary, and the real-axis cc-susceptibility
χ(q,ν,ω) (this task).

**Explicitly out of scope / future work:**
- The ordered phase and the associated elastic-sector self-energy machinery
  (α_m, ξ, λ₃ — HTML eqs 37-40) and the modified mean field H_MF (HTML eqs
  41-43). The Bx=0 boundary point is handled only via the closed-form Σc
  route (`invz_sigma_crit` / `invz_critical_T0field`), not by solving inside
  the ordered phase.
- The full 3×3 tensor / 4×4 sublattice RPA observable layer (neutron
  cross-sections, transverse components): `MF_RPA_Yikai.m` (repository root)
  remains the tool for that; `invz_jq_modes`'s per-sublattice `Jcc(q)` matrix
  is the natural hook point for extending this module to it.
- Free energy, δU, and heat capacity (HTML eqs 34-35, 44).
- Ewald acceleration of the dipole lattice sums — brute-force summation with
  `dpRng` real-space cutoff plus on-disk caching (`invz/cache/`) is sufficient
  at the 12-24³ grids used here.
- **Small-Bx accuracy caveat:** near zero field the electronic doublet
  splitting Δ can become comparable to or smaller than the ⁶⁷Ho hyperfine
  coupling scale (A ≈ 3.36 μeV). In that region the two-level (electronic
  doublet only) treatment underlying `invz_twolevel`/`invz_sigma` is expected
  to overestimate the fluctuation suppression of Tc, mirroring R 2007's own
  observation that the 1/z theory over-suppresses Tc relative to QMC in this
  regime (R 2007, comparison with quantum Monte Carlo). Results near the
  degenerate-doublet limit should be treated as directional only; this is why
  `invz_twolevel` raises `invz:degenerateDoublet` rather than silently
  returning a value when Δ < 1e-4 meV.
