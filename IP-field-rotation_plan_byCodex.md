# In-plane field-rotation investigation and implementation plan

Date: 2026-07-16  
Scope: rotation of an applied magnetic field within the crystallographic `ab` plane in the current scalar-cc `invz` codebase.

## Executive conclusion

The leading angular dependence of the LiHoF4 longitudinal response under

```matlab
Bvec = B*[cosd(phi_ab), sind(phi_ab), 0];
```

is generated locally by the non-axial crystal-field terms together with the Zeeman interaction. Bare crystal field + Zeeman reproduces nearly all of the angular shape seen after a physically consistent two-component transverse mean-field solve.

It is not, however, the complete quantitative theory:

- an isotropic, self-consistent `(hx,hy)` transverse mean field renormalizes the angular modulation, especially at lower fields;
- the current x-only feedback breaks C4 equivalence and substantially distorts an in-plane scan;
- finite-frequency `xz` and `yz` response channels are nonzero even when the field is exactly in-plane, so tensor RPA gives a small additional correction to the cc peak;
- ordered and near-critical spectra cannot be inferred from the bare single-ion calculation alone because the longitudinal molecular field and critical amplification remain essential.

The recommended implementation is therefore staged:

1. Preserve a crystal-field + Zeeman-only diagnostic mode.
2. Add an opt-in, always-active vector transverse mean field `(hx,hy)` for quantitative in-plane scalar-1/z spectra.
3. Add `phi_ab` to the spectra driver and require the vector transverse-MF mode for nonzero azimuth.
4. Validate the scalar cc result against a full 3x3 `Sigma=0` RPA reference over field, angle, and frequency.
5. Defer tensorized 1/z unless that reference shows material cc-mode errors.

## 1. Existing implementations and the relevant distinction

The external mean-field stack under

```text
/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/Mean Field/LiReF4
```

already handles a general field vector:

```matlab
Hx = Bfield*cos(phi)*cos(theta);
Hy = Bfield*sin(phi)*cos(theta);
Hz = Bfield*sin(theta);
```

`LiIonsF4.m`/`remf.m` iterate all moment components through the full dipolar/exchange mean field, `MF_Chi.m` builds the full 3x3 local susceptibility, and `RPA.m` solves

```matlab
chiq = (eye(3) - chi0*Jq) \ chi0;
```

with a 3x3 `chi0`. Its interaction is unit-cell averaged and diagonal-projected, but off-diagonal local susceptibility still participates in the matrix inversion.

The current `invz` route is different:

- `invz_single_ion` accepts `[Bx By Bz]` in the Zeeman Hamiltonian;
- its transverse molecular-field feedback remains only `hx = Jxx0*<Jx>`;
- `invz_chi0z` already returns the full 3x3 local susceptibility;
- the effective medium, self-energy, and final lattice response retain only the cc channel.

Consequently, simply passing `[B*cos(phi), B*sin(phi), 0]` through the present `invz` chain includes the local Zeeman rotation but does not give a C4-consistent transverse mean field.

## 2. Diagnostic question

The first question was whether the leading in-plane rotation effect could be accounted for entirely by the crystal-field and Zeeman terms.

Three single-ion models were compared at `T = 0.31 K`:

1. **Bare:** crystal field + Zeeman, with transverse MF disabled (`Jxx0 = 0`).
2. **Legacy:** current `invz_single_ion`, which iterates `hx` only.
3. **Vector MF:** diagnostic self-consistency in both transverse channels,

   ```matlab
   hx_new = Jaa0*<Jx>;
   hy_new = Jaa0*<Jy>;
   ```

   with both channels active for every angle, including the crystallographic axes.

The electronic doublet parameters and full electronuclear static `chi0_cc` were evaluated for `phi_ab = 0:5:90 deg` (10-degree spacing in the electronuclear vector-MF diagnostic).

## 3. Numerical findings

### 3.1 Bare CF + Zeeman already produces a sizeable fourfold effect

Relative angular span `(max-min)/mean` over one C4 sector:

| `B` | Bare electronic `Delta` | Bare electronic `M2` | Bare electronuclear `chi0_cc(0)` |
|---:|---:|---:|---:|
| 2 T | 5.27% | 0.089% | 4.60% |
| 4 T | 6.97% | 0.407% | 7.45% |
| 6 T | 8.24% | 0.770% | 8.48% |

Thus the leading angular sensitivity is mainly a modulation of the doublet splitting and the associated cc susceptibility; the longitudinal doublet matrix element `M2` is much less angle-sensitive.

At 4 T, the bare electronic splitting was

```text
Delta_min = 0.345693898281 meV at phi_ab about 15 deg
Delta_max = 0.370662538778 meV at phi_ab about 60 deg
```

The extrema are not tied to `0/45 deg`. Their shift is a physical consequence of the complete non-axial crystal field, including the sine Stevens coefficient `B64s`.

### 3.2 The non-axial crystal field is the origin of the angular pattern

At 4 T, bare CF + Zeeman diagnostics gave:

| Crystal-field model | Result |
|---|---|
| Full `B44`, `B64c`, `B64s` | 6.97% splitting modulation; extrema near 15/60 deg |
| `B64s = 0` | 1.52% modulation; extrema return to symmetry directions |
| `B44` only among the non-axial terms | Very large modulation, showing strong cancellation/interference with sixth-order terms in the physical parameter set |
| All non-axial `m=4` terms removed | In-plane angular dependence vanishes to numerical accuracy |

The in-plane effect is therefore not a new many-body anisotropy generated by the 1/z machinery. Its leading harmonic and orientation are already encoded in the local crystal-field Hamiltonian.

### 3.3 A proper vector MF preserves the bare angular shape but renormalizes it

Comparison of bare and vector-MF angular spans:

| `B` | Bare `Delta` | Vector-MF `Delta` | Bare `chi0_cc(0)` | Vector-MF `chi0_cc(0)` |
|---:|---:|---:|---:|---:|
| 2 T | 5.27% | 5.60% | 4.60% | 5.18% |
| 4 T | 6.97% | 7.28% | 7.45% | 7.81% |
| 6 T | 8.24% | 8.45% | 8.48% | 8.72% |

After fitting an offset and overall scale, the vector-MF angular curve differs from the bare curve by only approximately 0.5-1.3% of the total angular span. The local CF + Zeeman terms therefore account for essentially all of the *shape*.

The scale is not identical. For the doublet splitting, the fitted vector-MF amplification of the bare angular variation was approximately:

```text
2 T: 1.18
4 T: 1.10
6 T: 1.05
```

This is small-to-moderate at high field but material at lower field. Near a critical point, even a modest change in the local splitting can be amplified strongly by the lattice denominator.

### 3.4 The current x-only MF is not acceptable for an in-plane scan

The legacy x-only feedback does not merely rescale the bare result. Its residual angular-shape error after an affine fit is about 9-28% of the angular span across the tested fields and observables.

At 4 T and `T = 0.31 K`, for example:

```text
legacy [4 0 0]: E2-E1 = 0.369235620278 meV
legacy [0 4 0]: E2-E1 = 0.352382559983 meV
```

These two fields are related by a 90-degree C4 rotation and should be equivalent. With both `hx` and `hy` iterated, their splitting is identical to numerical precision:

```text
vector MF x/y axes: E2-E1 = 0.369276647675 meV
```

The reason is visible even on the x axis: because `B64s` is nonzero, `[4 0 0]` induces

```text
<Jy> = -0.0689991117463
```

in the current single-ion state. Suppressing `hy` discards feedback from a real transverse moment. Activating `hy` only when `By` is nonzero would remain discontinuous and C4-inconsistent; both transverse channels must be active in the vector-MF mode for every field direction.

### 3.5 Tensor RPA is a secondary but nonzero correction

For an exactly in-plane field, the diagnostics found:

- `<Jz>` is zero to numerical precision;
- static `chi_xz(0)` and `chi_yz(0)` are zero to numerical precision;
- finite-frequency `chi_xz(omega)` and `chi_yz(omega)` are not zero.

Therefore a scalar static argument does not prove that the full dynamical cc response is tensor-decoupled.

A representative high-field comparison was performed at

```text
T = 0.31 K
B = 6 T
eta = 1e-3 meV
Sigma = 0
J = diag(Jaa0, Jaa0, Jcc0)
```

using the vector-MF single-ion state:

| Result | Scalar cc RPA | Full 3x3 RPA |
|---|---:|---:|
| Angular span of cc peak energy | 13.43% | 13.66% |
| Maximum tensor-induced peak shift | — | about 0.00245 meV |

The full tensor calculation preserved the locations of the angular extrema. This supports treating tensor mixing as a secondary correction for the high-field cc peak position, but it is only one paramagnetic, `Sigma=0`, uniform-mode diagnostic. It is not sufficient evidence for the critical band, ordered phase, transverse observables, or full spectral lineshapes.

Pointwise or unnormalized spectral differences are not a useful metric for narrow peaks: a small energy shift can produce a large L2 difference. Validation should separately track peak energy, integrated spectral weight, and a broadened/peak-normalized shape metric.

## 4. Answer to the diagnostic question

**Leading effect:** yes. Crystal field + Zeeman accounts for nearly all of the in-plane angular pattern in the tested local quantities and high-field cc peak.

**Entire quantitative effect:** no. A vector transverse MF changes the modulation scale, finite-frequency tensor mixing shifts the RPA pole, and ordered/near-critical feedback can amplify small local differences. Bare CF + Zeeman should be retained as a diagnostic/reference model, not used as the final production theory across the phase diagram.

## 5. Recommended code design

### 5.1 Add an explicit transverse-MF mode

Preserve current benchmarks by making the new behavior opt-in:

```matlab
opts.transverse_mf = 'legacy_x';  % default: current behavior
% 'none'       : CF + Zeeman only
% 'vector_ab'  : self-consistent hx and hy
```

`legacy_x` must remain the default until the project explicitly elects to rebaseline all historical x-field results.

For `vector_ab`, define the current tetragonal transverse coupling as

```matlab
Jperp0 = Jaa0*eye(2);
hperp_new = Jperp0*si.Jexp(1:2);
```

or equivalently

```matlab
hx_new = Jaa0*jx;
hy_new = Jaa0*jy;
```

The matrix form is preferable because it makes the physical approximation explicit and leaves a clean extension point for a future full uniform transverse coupling matrix.

In `invz_single_ion`:

```matlab
H = H0 - hx*Jx - hy*Jy - hz*Jz;
dmf = max([abs(hx_new-hx), abs(hy_new-hy), abs(hz_new-hz)]);
```

Requirements:

- seed `hx = hy = 0`;
- mix `hy` identically to `hx`;
- iterate both components in paramagnetic, ordered, and `hz_fixed` modes;
- never guard `hy` on `By ~= 0`;
- return `si.hy` and include it in MF diagnostics/free-energy bookkeeping;
- retain the existing sign-aware longitudinal branch handling unchanged.

### 5.2 Thread the mode through the entire local-state chain

The electronuclear state used in `chi0_cc` and the electronic state used for the two-level self-energy must use the same transverse-MF model.

Propagate `transverse_mf` and, if exposed, `Jperp0` through:

- `invz_single_ion.m`
- `invz_twolevel.m`
- `invz_twolevel_ordered.m`
- `invz_solve_point.m`
- `invz_solve_point_ordered.m`
- `invz_solve_auto.m`
- the fallback construction in `invz_chi_realaxis.m`
- `invz_spectra_map.m` and `invz_spectra_qpath.m`

The existing `solve_opts` route in the spectra functions is the natural carrier for `transverse_mf`. Live couplings remain driver-owned: `Jaa0` should be derived from `info.Jaa0` and used for both `a` and `b` in the current spheroidal/tetragonal model.

### 5.3 Add the azimuthal driver knob

Extend the spectra driver direction to

```matlab
theta_c = 0;   % out-of-plane tilt toward c
phi_ab  = 0;   % in-plane rotation a -> b

dhat = [cosd(theta_c)*cosd(phi_ab), ...
        cosd(theta_c)*sind(phi_ab), ...
        sind(theta_c)];
```

For `phi_ab ~= 0`, the driver should either:

- automatically select `transverse_mf = 'vector_ab'` and print/store that model choice; or
- reject `legacy_x` with a clear error requiring the user to opt in.

Silently running a rotated field with `legacy_x` should not be allowed because the result is not C4-consistent.

For `theta_c = 0`, the existing transverse ordered-first/paramagnetic-fallback routing remains valid. For simultaneous `theta_c ~= 0` and `phi_ab ~= 0`, combine the already implemented forced-moment longitudinal route with `vector_ab` transverse MF.

### 5.4 Store enough metadata to make results reproducible

Add to spectra results:

```matlab
S.theta_c
S.phi_ab
S.field_dir
S.Bvec
S.transverse_mf
S.Jperp0
```

Titles and logs should state both angles and the MF model. This avoids comparing a legacy-x result with a vector-MF result as though only the angle differed.

## 6. Tensor-reference strategy

The first tensor step need not be tensorized 1/z. Build a local/reference 3x3 `Sigma=0` RPA from the already available `invz_chi0z`:

```matlab
chi_full = (eye(3) - chi0*diag([Jaa0 Jaa0 Jsel])) \ chi0;
chi_scalar_cc = chi0(3,3)/(1 - Jsel*chi0(3,3));
```

Use the exact same converged vector-MF single-ion state on both sides. Compare over:

- `phi_ab = 0:5:90 deg`;
- representative fields below, near, and above `Bc`;
- positive and negative frequency grids or the existing positive real-axis convention;
- both electronic-only and electronuclear states where useful;
- `q=0` and selected q-path couplings if the reference is extended beyond the uniform mode.

Report separately:

1. peak-energy difference;
2. integrated cc spectral-weight difference;
3. peak-normalized spectral-shape difference;
4. angular harmonic coefficients;
5. angle of the extrema.

This reference answers whether tensor mixing changes the *in-plane anisotropy*, not merely whether scalar and tensor spectra differ at one angle.

If the tensor correction remains small for cc peak positions across the target regime, the scalar 1/z route with vector transverse MF is a justified staged approximation. If it becomes material near the crossover or ordered phase, the next step is the previously identified A0/A1-style tensor propagation rather than an ad hoc correction to scalar `chi0_cc`.

## 7. Acceptance tests

### 7.1 Backward compatibility

1. Default `transverse_mf = 'legacy_x'` reproduces all existing theta/phi-zero tests and benchmarks bit-for-bit.
2. Explicit `'legacy_x'` equals the omitted/default option.
3. `'none'` with `Jaa0` present is identical to the current `Jxx0 = 0` bare calculation.

### 7.2 Vector-MF correctness

1. C4 axes:

   ```text
   spectrum([B 0 0]) == spectrum([0 B 0])
   <Jx>([B 0 0]) == <Jy>([0 B 0])
   <Jx>([0 B 0]) == -<Jy>([B 0 0])
   ```

2. Periodicity:

   ```text
   observable(phi) == observable(phi + 90 deg)
   ```

3. `Jz = 0` for a purely in-plane field on the paramagnetic single-ion path.
4. `si.hy = Jaa0*si.Jexp(2)` at convergence.
5. The returned MF residual includes all active `hx`, `hy`, and `hz` channels.
6. The variational MF free energy includes `0.5*hy*<Jy>` in vector mode.

### 7.3 End-to-end consistency

1. `si` and `tl` record/use the same transverse-MF mode.
2. Scalar and `[B 0 0]` inputs remain identical in legacy mode.
3. A nonzero `phi_ab` cannot silently enter legacy-x mode.
4. Ordered and `hz_fixed` electronic rebuilds reproduce the same `hx/hy` convention as the electronuclear state.
5. `phi_ab -> phi_ab + 90 deg` leaves the scalar 1/z cc spectrum invariant within numerical tolerance.

### 7.4 Diagnostic physics tests

1. With all non-axial `m=4` Stevens terms zeroed, the bare CF+Zeeman spectrum is invariant under continuous in-plane rotation within numerical tolerance.
2. With the physical crystal field, fit

   ```text
   O(phi) = A0 + A4c*cos(4*phi) + A4s*sin(4*phi) + A8c*cos(8*phi) + ...
   ```

   and report the dominant amplitude and principal-axis phase.
3. Reproduce the recorded bare/vector-MF spans above within stated numerical tolerances.
4. Run the full-3x3 `Sigma=0` RPA comparison and log the chosen peak/weight/shape metrics rather than relying on pointwise relative differences at spectral zeros.

## 8. Recommended delivery sequence

### Phase IP0 — diagnostic utility

- Add no production behavior change.
- Produce bare CF+Zeeman angular scans for `Delta`, `M2`, `chi0_cc`, and peak energy.
- Fit C4 harmonics and record extrema.
- Add the 3x3 `Sigma=0` RPA reference at high field, then near the target crossover.

### Phase IP1 — opt-in vector transverse MF

- Add `transverse_mf = 'legacy_x'|'none'|'vector_ab'`.
- Implement `hy` with no angle/component guard.
- Thread the option through all local-state builders.
- Add C4, convergence, free-energy, and compatibility tests.

### Phase IP2 — driver and spectra integration

- Add `phi_ab` and the full direction formula.
- Require/activate vector MF for nonzero azimuth.
- Store model/angle metadata and update titles/logs.
- Test field-sweep and q-path views in paramagnetic and ordered regimes.

### Phase IP3 — quantitative scalar-versus-tensor validation

- Sweep field, angle, and frequency with the common vector-MF state.
- State a supported regime for scalar cc peak predictions using explicit error thresholds.
- Decide whether A0/A1 tensor propagation is warranted.

## 9. Final recommendation

Use crystal field + Zeeman as the first diagnostic and as the explanation of the leading fourfold anisotropy. Do not treat it as the complete production calculation, and do not rotate the field through the current x-only mean-field path.

The smallest physically consistent production extension to the current scalar 1/z code is an opt-in, always-active `(hx,hy)` transverse mean field with the same `Jaa0` in both tetragonal directions. It preserves the established scalar self-energy machinery, restores exact C4 equivalence, and can be validated directly against the existing full-tensor RPA capability before committing to a research-scale tensor 1/z implementation.
