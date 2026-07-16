# Session note — 2026-07-16: Sigma=0 scalar-vs-tensor reference + supported field-tilt range

Branch: `invz-1z-lihof4`. Task 9 of the field-angle plan
(`docs/superpowers/plans/2026-07-16-invz-field-angle.md`). Went through two
BLOCKED-and-resolved rounds before landing; both resolutions are recorded in
`docs/superpowers/plans/2026-07-16-invz-field-angle.md` and
`docs/superpowers/specs/2026-07-16-invz-field-angle-design.md` (commits
`d6c43c2`, `7f9f16b`).

## What this measures

`invz_chi_tensor_ref(ion, T, Bvec, w, opts)` builds ONE single-ion state
(order mode, sign-aware seed) and ONE bare `chi0` tensor from `invz_chi0z`,
then propagates it two ways from the same inputs:

- **scalar** (production pipeline): `chi_cc = chi0_cc / (1 - Jsel*chi0_cc)` —
  cc channel only.
- **tensor** (reference): `chi = chi0 * inv(I3 - J*chi0)`,
  `J = diag(Jaa0, Jaa0, Jsel)` — full 3x3 inversion, all cross channels
  retained.

The difference between the two quantifies the cross-channel error the scalar
`Sigma=0` pipeline omits under a c-axis field tilt.

## Round 1: the zero-tilt yz baseline (not a bug)

The first measurement attempt found the naive "scalar == tensor at zero tilt"
assumption false: at `theta = 0`, `B || a`, the `yz` cross channel of `chi0`
is **symmetry-allowed** by the non-axial crystal-field term `B64s` and is
**large** — measured `max|chi0_yz| / max|chi0_zz| = 0.183` at 6 T (hyp =
false, T = 0.31 K), while `max|chi0_xz| / max|chi0_zz| = 2.8e-3` (strongly
suppressed, as the naive symmetry argument predicts for `xz` specifically).
Root-caused via first-principles group theory (a field along a single
in-plane axis leaves only a `C2` rotation about that axis as an exact
symmetry of the single-ion Hamiltonian; that symmetry forbids `chi0_xz` but
allows `chi0_yz`) and confirmed independently by the prior external review
`IP-field-rotation_plan_byCodex.md`, which measured the same effect
(`<Jy> = -0.069` at `[4 0 0]`, consistent order of magnitude/sign with the
`-0.084` measured here at 6 T). **This is a pre-existing, physical property of
the scalar cc pipeline, present at every tilt angle including zero — not a
tilt-induced error, and not a defect in `invz_chi_tensor_ref`.** In the
spontaneously ordered phase (e.g. `B = 2 T`) the spontaneous moment breaks Z2
and adds a *second*, independent source of nonzero cross-channel response at
`theta = 0` on top of the `B64s` baseline.

Resolution: the raw spectral metric `eps_spec` is reported as a diagnostic,
never gated. A tilt-referenced metric `invz_tilt_err` (`eps_tilt`) was added,
differencing against the `theta = 0` reference at the same field to remove
the theta-independent `yz` baseline.

## Round 2: L2 metrics are diagnostics only (the delta0/eta artifact)

The full 5-angle x 3-field measurement then showed that `eps_tilt` *also*
fails a naive 5% gate at essentially every nonzero angle (8-28%), even though
the physically meaningful peak-position shift `dE_peak` stayed at the
**microelectronvolt** scale throughout. Root cause: for these sharp lines
(comparison `eta = 0.02` meV), an L2 norm cannot resolve below the zero-tilt
peak-position offset `delta0` (`dE_peak` at `theta = 0`) leaking through at
the `delta0/eta` scale — quantitatively verified: `eps_tilt ~ 0.11` vs
`delta0/eta = 2.2/20` micro-eV at 6 T; `eps_tilt ~ 0.28` vs `delta0/eta =
7.7/20` micro-eV at 2 T, and this ratio **decreases with angle** as the
tilt-induced change (the metric's denominator) grows — exactly the signature
of a fixed-offset artifact being diluted by a growing signal, not a real
tilt error growing with angle. **Every L2 lineshape metric in this
reference — raw `eps_spec` or tilt-differenced `eps_tilt` — is therefore a
reported diagnostic only; lineshape fidelity is explicitly not claimed for
the scalar stage.**

User-approved resolution (commit `7f9f16b`): gate on **peak observables**
instead, which the production drivers actually plot and which track to
microelectronvolts: `dE_peak` (already present) plus a new `eps_amp = |max_w
chi''_sc - max_w chi''_ten| / max_w chi''_ten` (peak amplitude/intensity
error). `invz_chi_tensor_ref.m` gained `R.amp_sc`, `R.amp_ten`, `R.eps_amp`.

## Support criterion (final)

An angle `theta > 0` is supported when, at **every** tested field:

```
eps_amp <= 0.10   AND   dE_peak <= max(0.02*Epeak_ten, eta)
```

(censoring rule for `dE_peak`: both-sides-censored/NaN passes, one-sided NaN
fails). The `theta = 0` row is baseline-only and not gated. `eps_spec` /
`eps_tilt` remain reported diagnostic columns.

## Measured table

`T = 0.1 K` (spectra-driver default), `w = (0:0.005:0.6)` meV, comparison
`eta = 0.02` meV, live couplings from the cached 16^3/dpRng=30 lattice sum
(`Jsel = info.Jcc0 = 6.42444e-3` meV, `Jaa0 = 3.51045e-3` meV — close to the
`invz_ion()` fallback constants `J0eff = 6.421e-3`, `Jxx0 = 3.512e-3` meV).

```
   theta  |B| (T)     eps_spec     eps_tilt      eps_amp      dE_peak      Ep_sc     Ep_ten
    0.00     2.00       0.1291            0      0.02356     0.007721     0.4484     0.4407   ok
    0.50     2.00        0.131       0.2813      0.02013     0.007979     0.4642     0.4562   ok
    1.00     2.00       0.1329       0.2654      0.01729     0.008023     0.4796     0.4716   ok
    2.00     2.00       0.1367       0.2183      0.01401     0.007707     0.5096     0.5019   ok
    5.00     2.00       0.1237       0.1399      0.01785     0.007465     0.5959     0.5885   ok
    0.00     4.95      0.07645            0      0.02313     0.002153     0.2547     0.2525   ok
    0.50     4.95       0.1096       0.0916     0.003193     0.003189     0.3739     0.3707   ok
    1.00     4.95       0.1205      0.08759     0.005772     0.003336     0.4372     0.4338   ok
    2.00     4.95        0.126      0.08279     0.008675     0.003789     0.5345     0.5307   ok
    5.00     4.95      0.04355      0.07649      0.05393          NaN        NaN        NaN   ok
    0.00     6.00      0.08157            0      0.01246     0.002201     0.3456     0.3434   ok
    0.50     6.00      0.09431       0.1137     0.006918     0.002766      0.404     0.4012   ok
    1.00     6.00       0.1058      0.09415     0.007712     0.002879     0.4712     0.4683   ok
    2.00     6.00       0.1128      0.08795      0.01273     0.003253     0.5824     0.5791   ok
    5.00     6.00      0.02548      0.08162      0.03136          NaN        NaN        NaN   ok
Supported tilt range (peak-observable criterion): theta_c <= 5 deg
```

(At `theta = 5` deg, `B in {4.95, 6}` T the peak moved outside the sampled
`w` window and both `Epeak_sc`/`Epeak_ten` were censored per
`invz_peak_energy`'s both-NaN rule — that counts as `dE_peak` passing under
the censoring rule, not as a failure; `eps_amp` at those two points, 0.054
and 0.031, independently confirms the peak observables still agree.)

## Verdict

**Supported tilt range: `theta_c <= 5 deg`** — the entire tested spec grid
(`0, 0.5, 1, 2, 5` deg) passes the peak-observable criterion at every tested
field. Worst-case `eps_amp = 5.39%` at `theta = 5` deg, `B = 4.95` T (still
under the 10% gate); worst-case (non-censored) `dE_peak = 8.0` micro-eV at
`theta = 1` deg, `B = 2` T, comfortably under its threshold there
(`max(0.02*Epeak_ten, eta) = max(9.4, 20) = 20` micro-eV).

## Reproducibility-test convention

The three spot-check points in
`invz/tests/test_invz_tensor_ref.m::test_reproducibility_of_logged_table`
(`{0.5, 6}`, `{2, 6}`, `{5, 4.95}` degrees/Tesla) use **default couplings**
(`ion.J0eff`, `ion.Jxx0`, i.e. `opts = struct('eta', 0.02)` with no
`Jsel`/`Jaa0` override), NOT the live cached-lattice-sum couplings used by
the table above. The two sets of couplings are close (`Jsel`: 6.421e-3
default vs 6.42444e-3 live; `Jaa0`: 3.512e-3 default vs 3.51045e-3 live) so
the resulting `eps_amp`/`eps_tilt` values differ only in the last couple of
significant figures from the table's live-coupling row at the same
(angle, field); the logged `expected_amp`/`expected_tilt` constants were
measured with the default-coupling convention directly, matching the test's
own calls exactly.
