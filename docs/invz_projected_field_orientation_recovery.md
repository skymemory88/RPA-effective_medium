# Projected-spin field-orientation recovery note

This note preserves the validated field-orientation work removed from the
production projected-spin stream during the 2026-07-29 reduction. The
production calculation is deliberately restricted to a purely transverse
field along the crystallographic `a` axis. The material below is retained so
that in-plane rotation can be restored after the spectra stream is stable.

## Physical model

For an in-plane field

```matlab
Bvec = B * [cosd(phi_ab), sind(phi_ab), 0];
```

the correct transverse mean-field model is `transverse_mf = 'vector_ab'`.
It solves

```text
hx = Jaa(0) <Jx>
hy = Jaa(0) <Jy>
```

self-consistently. The historical `legacy_x` model solved only `hx` and
forced `hy=0`. That restriction is not C4-consistent: the non-axial
crystal-field coefficient `B64s` permits a nonzero `<Jy>` even when the
applied field points along `a`. The `none` model (`hx=hy=0`) remains a useful
bare crystal-field-plus-Zeeman control, but it is not the interacting
transverse mean field.

The implementation formerly lived primarily in:

- `invz_common/invz_single_ion.m` (`transverse_mf`)
- `invz_projected/invz_run_spectra.m` (`phi_ab` and field direction)
- `invz_projected/invz_run_ip_scan.m` and `invz_c4fit.m` (validation)
- `invz_projected/invz_chi_tensor_ref.m` (scalar/tensor reference)

The last committed pre-removal versions and full measurement notes remain
recoverable from Git:

```text
git show 2a8d5b7^:docs/SESSION-2026-07-16-inplane-rotation.md
git show 2a8d5b7^:docs/SESSION-2026-07-16-field-angle.md
```

## Validation anchors

The single-ion angular scan used `T=0.31 K`, hyperfine off,
`B={2,4,6} T`, and angles `0:5:90` plus `11` and `79` degrees.
`vector_ab` produced a stable C4 principal axis near -29 to -30 degrees and
a monotonic angular splitting span:

| B (T) | splitting span | C4 amplitude | principal axis |
|---:|---:|---:|---:|
| 2 | 5.608% | 2.810% | -28.94° |
| 4 | 7.289% | 3.643% | -30.04° |
| 6 | 8.464% | 4.239% | -29.63° |

By contrast, `legacy_x` moved from +29.37° at 2 T to -36.24° at 4 T,
including a sign flip, and was rejected for rotated fields.

The Sigma=0 scalar-versus-tensor comparison used `T=0.1 K`,
`B={2,4.95,6} T`, angles `{0,5,11,15,30,45,60,75,79,90}` degrees,
`w=0:0.005:0.6 meV`, and `eta=0.02 meV`. Every tested point passed

```text
relative peak-amplitude error <= 0.10
peak-energy error <= max(0.02 * tensor peak energy, eta)
```

The worst amplitude error was 2.66%; the worst peak-energy error was
7.87 micro-eV. These bounds apply to peak observables, not full spectral
lineshape fidelity.

The external crystal-field rotation convention was pinned numerically:

```text
ion.cfRot = -11 degrees  <=>  phi_ab = -11 degrees
```

The same-sign mapping agreed to `1.8e-14 meV`; the opposite sign differed by
`0.435 meV`.

## Restoration constraints

1. Restore `vector_ab` end to end; do not restore an angle knob that feeds a
   single-axis (`legacy_x`) mean field.
2. Re-run the angular C4 scan and scalar/tensor peak-observable gate.
3. Do not compare `legacy_x` and `vector_ab` runs as though only the angle
   changed; the two are different mean-field models even at `phi_ab=0`.
4. Out-of-plane tilt and in-plane rotation were validated only on separate
   one-dimensional slices. Their simultaneous use was never validated.
5. The projected scalar self-energy still dresses only the cc channel.
   Arbitrary-orientation fidelity ultimately requires the full-tensor route.
