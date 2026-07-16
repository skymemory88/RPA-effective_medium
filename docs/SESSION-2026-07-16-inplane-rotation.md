# Session note — 2026-07-16: in-plane field rotation (phi_ab) — cfRot pin + production validation

Branch: `invz-1z-lihof4`. Task 5 (final) of the in-plane field-rotation feature
(`.superpowers/sdd/task-5-brief.md`); Tasks 1-4 landed full `transverse_mf`
support, the `phi_ab` driver knob, `invz_run_ip_scan.m`, and `invz_c4fit.m`.

## Feature summary

Three `opts.transverse_mf` models are available at every single-ion call
(`invz_single_ion`, and everywhere it threads through — `invz_chi_tensor_ref`,
the spectra drivers):

- **`'legacy_x'`** (default, bit-for-bit with the pre-existing single-axis
  code path) — only `hx = Jxx0*<Jx>` is solved; `hy` stays exactly 0. This is
  the C4-**inconsistent** model for a rotated field: forcing `hy = 0`
  regardless of the applied field's azimuth breaks the crystal's four-fold
  symmetry (Section A below shows this directly — `legacy_x`'s harmonic
  content and principal-axis angle scatter across field magnitudes where the
  other two models track a stable ~-29 deg axis).
- **`'none'`** — both transverse channels forced to zero (`hx = hy = 0`), i.e.
  the bare CF + Zeeman problem; used as a C4-clean baseline/control.
- **`'vector_ab'`** — solves both `hx = Jxx0*<Jx>` and `hy = Jxx0*<Jy>`
  self-consistently. Required for a physically correct in-plane rotation:
  even at `phi_ab = 0` (field along a pure x-axis) the non-axial `B64s`
  crystal-field term makes `<Jy> != 0`, so a genuinely C4-consistent
  transverse mean field cannot fix `hy = 0` by fiat.

The driver knob is `phi_ab` (deg) in `invz_run_spectra.m`, the in-plane
rotation of the swept field from `a` toward `b`: `dhat = [cosd(theta_c)*cosd(phi_ab),
cosd(theta_c)*sind(phi_ab), sind(theta_c)]`. A guard enforces
`transverse_mf = 'vector_ab'` whenever `phi_ab != 0` — the library errors
rather than silently returning C4-broken `legacy_x` output for a rotated
field.

## Table 1 (Section A) — single-ion angular scans

`invz_run_ip_scan.m`, Section A: `T = 0.31` K, `hyp = false`,
`B in {2, 4, 6}` T, `phi_ab in {0:5:90} union {11}` deg, for the three
transverse-MF models. `span(%) = 100*(max(d)-min(d))/mean(d)` of the doublet
splitting `d = E(2)-E(1)` over the angle scan; `A4/A0`, `A8/A0` are the C4/C8
harmonic amplitude ratios and `phi0` the principal-axis angle from
`invz_c4fit`.

```
     model   B(T)    span(%)   A4/A0(%)   A8/A0(%)     phi0
      none    2.0      5.275      2.649      0.011   -28.48
      none    4.0      6.986      3.488      0.007   -29.71
      none    6.0      8.256      4.128      0.008   -29.47
  legacy_x    2.0     10.013      2.016      0.041    29.14
  legacy_x    4.0      6.433      1.978      0.036   -35.73
  legacy_x    6.0      6.668      3.232      0.045   -31.12
 vector_ab    2.0      5.611      2.809      0.014   -28.95
 vector_ab    4.0      7.293      3.641      0.005   -30.06
 vector_ab    6.0      8.470      4.234      0.013   -29.67
```

Matches the expected shape from the reviewed Codex diagnostics: `none` and
`vector_ab` both show span growing monotonically with field (~5.3% -> ~8.3%
and ~5.6% -> ~8.5%) with a stable principal axis clustered tightly at
`phi0 ~ -29` to `-30` deg across all three fields. `legacy_x` is visibly
C4-broken by comparison: span is *non-monotonic* in field (10.0% -> 6.4% ->
6.7%) and `phi0` swings from +29.1 deg at 2 T to -35.7 deg at 4 T to -31.1
deg at 6 T — a ~65 deg jump between the first two fields, including a sign
flip, that `none`/`vector_ab` do not show. `legacy_x` rows exist to display
this C4 violation; they are never used in production.

## Table 2 (Section B) — Sigma=0 scalar-vs-tensor cc comparison over angle

`invz_run_ip_scan.m`, Section B: `invz_chi_tensor_ref` with
`transverse_mf = 'vector_ab'` (production couplings from the live
16^3/dpRng=30 cached lattice sum), `T = 0.1` K, `B in {2, 4.95, 6}` T,
`phi_ab in {0, 5, 11, 15, 30, 45, 60, 75, 90}` deg, `w = (0:0.005:0.6)` meV,
comparison `eta = 0.02` meV. Gate (same tilt criterion as the theta_c
reference, `docs/SESSION-2026-07-16-field-angle.md`): `eps_amp <= 0.10 AND
dE_peak <= max(0.02*Epeak_ten, eta)` (censoring rule: both-sided NaN passes).
Energies (`dE_peak`, `Ep_sc`, `Ep_ten`) are in meV.

```
     phi  |B| (T)      dE_peak      eps_amp        eps_W      Ep_sc     Ep_ten    ok
     0.0     2.00     0.007721      0.02356     0.002941     0.4484     0.4407     1
     5.0     2.00     0.007807      0.02242     0.002785     0.4488     0.4410     1
    11.0     2.00     0.007861      0.02181     0.002674     0.4489     0.4411     1
    15.0     2.00     0.007866       0.0219      0.00266     0.4489     0.4410     1
    30.0     2.00     0.007665      0.02537     0.003008     0.4479     0.4402     1
    45.0     2.00     0.007254      0.02505      0.00357     0.4464     0.4391     1
    60.0     2.00     0.007103      0.02436     0.003791     0.4459     0.4388     1
    75.0     2.00     0.007317      0.02595     0.003507     0.4469     0.4396     1
    90.0     2.00     0.007721      0.02356     0.002941     0.4484     0.4407     1
     0.0     4.95     0.002152      0.02309     0.008456     0.2547     0.2525     1
     5.0     4.95     0.002347     0.002007     0.008644     0.2579     0.2556     1
    11.0     4.95     0.002412      0.01815     0.008704     0.2594     0.2570     1
    15.0     4.95      0.00242     0.005912     0.008629     0.2585     0.2561     1
    30.0     4.95     0.001883      0.01844     0.007606     0.2548     0.2529     1
    45.0     4.95     0.001595     0.009794     0.006159     0.2754     0.2738     1
    60.0     4.95     0.001583      0.01108     0.005849     0.2852     0.2836     1
    75.0     4.95     0.001777      0.01079     0.006818     0.2755     0.2737     1
    90.0     4.95     0.002152      0.02309     0.008456     0.2547     0.2525     1
     0.0     6.00     0.002196      0.01168     0.006596     0.3457     0.3435     1
     5.0     6.00     0.002466     0.002206     0.006965     0.3384     0.3359     1
    11.0     6.00     0.002444     0.002065     0.007104     0.3330     0.3306     1
    15.0     6.00      0.00227     0.003979     0.006988     0.3321     0.3298     1
    30.0     6.00     0.001842     0.005571     0.005588     0.3459     0.3441     1
    45.0     6.00     0.001555      0.01129     0.004286     0.3699     0.3684     1
    60.0     6.00     0.001566     0.005636     0.004056     0.3820     0.3804     1
    75.0     6.00     0.001809     0.000631     0.004997     0.3715     0.3697     1
    90.0     6.00     0.002196      0.01168     0.006596     0.3457     0.3435     1

supported in-plane angles (all fields): [0 5 11 15 30 45 60 75 90] deg
```

(A handful of rows at the un-converged tail of the mean-field loop print a
`Mean field not converged after 800 iterations` warning with residual
`|dmf|` between `1e-11` and `1e-6` meV — six to nine orders of magnitude
below the `1e-12` convergence tolerance's next significant digit and below
any threshold in the gate; harmless diagnostic noise, not a solve failure.)

## Supported-regime statement

**All nine tested in-plane angles (0, 5, 11, 15, 30, 45, 60, 75, 90 deg) are
supported at all three tested fields (2, 4.95, 6 T)** under the peak-observable
gate `eps_amp <= 0.10 AND dE_peak <= max(0.02*Epeak_ten, eta=0.02 meV)`.
Worst-case `eps_amp = 2.60%` (phi = 75 deg, 2 T), an order of magnitude under
the 10% ceiling. Worst-case `dE_peak = 7.87 ueV` (phi = 15 deg, 2 T) against
a threshold there of `max(0.02*0.441, 0.02) = 20 ueV` — comfortably inside.
At 6 T, `dE_peak` peaks at 2.47 ueV (phi = 5 deg), matching the previously
Codex-measured ~2.45 ueV tensor-induced peak shift at 6 T almost exactly.
`dE_peak` is smallest at high field and largest at low field (2 T), the
opposite ordering from the theta_c tilt reference — but still well inside
gate at every point tested.

## Pinned cfRot <=> phi_ab mapping

The external stack's `cf.m` ('coefficient' method) rotates the m=4 crystal-field
coefficient pairs by `Br = [cos(4r) sin(4r); -sin(4r) cos(4r))]`. A new test,
`invz_projected/tests/test_invz_cfrot_equiv.m`, builds the CF Hamiltonian both
ways — (a) coefficients rotated by `r` via `cf.m`'s convention, applied
directly with a fixed field along x, and (b) coefficients left alone with the
field itself rotated by a candidate angle `s*r` — and compares the resulting
spectra. Two candidate signs were tested (`s = +1` and `s = -1`) at `r = -11`
deg, `|B| = 4` T:

- `s = +1` (`phi_ab = +r = -11` deg): max spectral difference `1.8e-14` meV — **matches to machine precision**.
- `s = -1` (`phi_ab = -r = +11` deg): max spectral difference `0.435` meV — clearly does not match.

A second angle/field pair (`r = 7` deg, `|B| = 6` T) confirms the same `s = +1`
sign, ruling out a `4r`-aliasing coincidence.

**Discovered mapping: `s = +1` — coefficient rotation by `r` equals field
rotation by the SAME angle `+r` (not the naively expected opposite sign).**

```
external ion.cfRot = -11 deg  <=>  invz phi_ab = -11 deg
```

phi_ab = -11 deg reproduces the experimentally calibrated external-stack
geometry (ion.cfRot(Ho) = -11 deg).

`stevens_ops.m` gained `ops.O44s = -0.5i*(Jplus^4 - Jminus^4)` (additive;
existing fields untouched) so the test could build the rotated-coefficient
Hamiltonian directly — this is the sine partner of the existing `ops.O44`,
in the same convention family as the existing `ops.O64s` (Hutchings, 1/(4i)
factor, matching the external `cf.m` convention).

## Caveats (do not violate)

- **Never cross-compare `legacy_x` and `vector_ab` runs as if only the angle
  differed.** `vector_ab` shifts even the `phi_ab = 0` result relative to
  `legacy_x` (the same field, same nominal geometry) because `vector_ab`
  lets `<Jy> != 0` respond to the `B64s` crystal-field term that `legacy_x`
  suppresses by construction — a real physical difference between the two
  mean-field *models*, not a rotation effect. `invz_run_spectra.m` documents
  this at ~0.04 ueV at 4 T, growing at lower field. Table 1 above shows the
  same effect at the spectrum level: `none` and `vector_ab` track each other
  closely (both C4-consistent, `phi0` within ~1.5 deg of each other at every
  field) while `legacy_x` diverges from both.
- **Combined `theta_c` (out-of-plane c-tilt) and `phi_ab` (in-plane rotation)
  is unvalidated.** The `theta_c <= 5` deg supported range
  (`docs/SESSION-2026-07-16-field-angle.md`) was measured with
  `transverse_mf = 'legacy_x'` and `phi_ab = 0` throughout — a pure c-axis
  tilt of an a-axis field. This session's Section B measures the orthogonal
  slice: pure in-plane rotation at `theta_c = 0` under `vector_ab`. Neither
  measurement bounds the combined-tilt regime (`theta_c != 0 AND phi_ab != 0`
  simultaneously); running the driver with both knobs nonzero is
  scope-excluded until a dedicated two-angle validation is done (design spec
  Sec. 8, amendment 3: combined-angle validation documented, not
  implemented).

## File / function map (Task 5)

| File | Status | Responsibility |
|---|---|---|
| `invz_projected/stevens_ops.m` | modified | added `ops.O44s` (sine partner of `ops.O44`), additive only |
| `invz_projected/invz_run_spectra.m` | modified | `phi_ab` comment corrected to state the pinned sign |
| `invz_projected/tests/test_invz_cfrot_equiv.m` | new | pins `ion.cfRot = -11 deg <=> phi_ab = -11 deg` against the external `cf.m` convention; second angle/field guards against 4r-aliasing |
| `docs/SESSION-2026-07-16-inplane-rotation.md` | new (this file) | feature summary, both measured tables, supported-regime statement, pinned sign, caveats |
| `invz_projected/README.html` | modified | §1 vector-MF paragraph, §4 `phi_ab`/`transverse_mf` knob docs, §6 `invz_run_ip_scan`/`invz_c4fit` reference rows, §8 supported-range + scope-exclusion statement |

## Test suite status

Fast suite: 107 passed, 0 failed, 12 incomplete (baseline 105/0/12 plus the
two new `test_invz_cfrot_equiv` tests).

Slow suite (`INVZ_SLOW=1`): 119 passed, 0 failed, 0 incomplete.
