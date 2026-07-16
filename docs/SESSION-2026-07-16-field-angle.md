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

> Provenance: this table was captured at commit 4938c8f, before (a) the wmin mask on the amp fields (ea7ca89) and (b) the "base" label for the ungated theta = 0 rows. The mask is verified immaterial here -- the slow reproducibility test re-measured the extreme cells against the masked code to 1%, and masked == unmasked max was confirmed directly at several points (the electronic peak dominates everywhere tested); theta = 0 rows in this capture read "ok" but are report-only baselines.

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
the table above. The two sets of couplings are close in absolute terms
(`Jsel`: 6.421e-3 default vs 6.42444e-3 live; `Jaa0`: 3.512e-3 default vs
3.51045e-3 live), but `eps_amp` and `eps_tilt` do **not** inherit that
closeness equally: `eps_tilt` (an L2 norm dominated by the tilt-induced
signal) shifts by only ~0.1% between the two conventions, while `eps_amp` is
a *small-difference* metric (`|amp_sc - amp_ten|` is only ~1% of `amp_ten`),
so the same tiny per-spectrum coupling shift is amplified into a much larger
RELATIVE move -- at `theta = 2` deg, `B = 6` T the live-coupling row's
`eps_amp = 0.01273` vs. the default-coupling `expected_amp(2) = 0.01135`, a
~12% relative difference, even though the underlying spectra are nearly
identical. The logged `expected_amp`/`expected_tilt` constants were measured
with the default-coupling convention directly, matching the test's own calls
exactly -- so the test's 1% `RelTol` reproducibility check is sound (it never
crosses conventions). Cross-convention agreement between this paragraph's
default-coupling spot checks and the live-coupling table above should be
read loosely, and only for `eps_tilt`; `eps_amp` is not expected to agree
closely across the two conventions.

## File / function map (Task 10)

Full plan: `docs/superpowers/plans/2026-07-16-invz-field-angle.md` ("File
Structure"). Design spec: `docs/superpowers/specs/2026-07-16-invz-field-angle-design.md`.

| File | Status | Responsibility |
|---|---|---|
| `invz/invz_field_vec.m` | new | scalar/3-vector field normalization (single validation point, `invz:fieldVec`); every solver boundary calls this once |
| `invz/invz_single_ion.m` | modified | sign-aware `mz_seed`, `mf_converged`/`mf_iters`/`mf_residual`, `E0`, `F_mf` |
| `invz/invz_twolevel.m`, `invz/invz_twolevel_ordered.m` | modified | accept scalar-or-vector `B` |
| `invz/invz_solve_point.m` | modified | normalizes `B` via `invz_field_vec` once, passes the vector down |
| `invz/invz_solve_point_ordered.m` | modified | `opts.forced_moment`, sign-aware seed + mirrored retry, complete early returns, `pt.moment_branch` |
| `invz/invz_solve_auto.m` | modified | `opts.bz_tol` dead band + longitudinal routing (exclusively through `invz_solve_point_ordered` with `forced_moment = true` once `|Bz| > bz_tol`) |
| `invz/invz_chi_realaxis.m` | modified | vector-capable paramagnet fallback |
| `invz/invz_spectra_map.m` | modified | `opts.field_dir`/`bz_tol`/`solve_opts` API, dead-band resolved once before the `parfor`, longitudinal failure contract, `S.field_dir`/`S.Bvec` metadata |
| `invz/invz_spectra_qpath.m` | modified | vector `B`, `S.Bvec`/`S.Bmag`, `mat2str`-based vector error message, same `opts` contract as `invz_spectra_map` |
| `invz/invz_run_spectra.m` | modified | `theta_c` knob, `dhat` field direction, tilt-aware plot labels |
| `invz/invz_plot_spectra_map.m` | modified | axis label `\|B\| (T)` (magnitude, not signed `Bx`) |
| `invz/invz_chi_tensor_ref.m` | new | Sigma=0 scalar-vs-tensor cross-channel reference at one `(T, Bvec)` |
| `invz/invz_tilt_err.m` | new | tilt-only metric `eps_tilt`: differences `invz_chi_tensor_ref` against the `theta=0` reference at the same field, to remove the theta-independent `yz` baseline (diagnostic only, see "Round 2" above) |
| `invz/invz_run_tensor_ref.m` | new | measurement driver: prints the theta/\|B\| table of `eps_spec`/`eps_tilt`/`eps_amp`/`dE_peak` that sets the supported tilt range |
| `invz/invz_peak_energy.m` | new (factored out) | shared censored, parabolic-refined peak-energy helper used by `invz_spectra_map`, `invz_spectra_qpath`, and `invz_chi_tensor_ref` |
| `invz/tests/test_invz_field_vec.m` | new | Task 1 tests |
| `invz/tests/test_invz_field_angle.m` | new | Tasks 2-4 solver-level tests |
| `invz/tests/test_invz_field_angle_spectra.m` | new | Tasks 5-7 spectra/label tests + fast mirror test |
| `invz/tests/test_invz_field_angle_slow.m` | new | Task 8 slow tests (`INVZ_SLOW`-gated) |
| `invz/tests/test_invz_tensor_ref.m` | new | Task 9 reference tests (3 fast structural + 1 slow reproducibility) |
| `docs/SESSION-2026-07-16-field-angle.md` | new (this file) | measured tensor-reference table + supported angle range |
| `invz/README.html` | modified | function reference (`invz_field_vec`, `invz_chi_tensor_ref`, `invz_run_tensor_ref`) + `theta_c` knob documentation (`§4`) + scope entries (`§8`) |

## New option / metadata names

| Name | Where | Meaning |
|---|---|---|
| `opts.field_dir` | `invz_spectra_map` | nonzero finite real 3-vector, normalized internally; sets the sweep direction of `fields` (`\|B\|`); default `[1 0 0]` (pure transverse, unchanged behaviour); invalid input -> `invz:fieldDir` |
| `opts.bz_tol` | `invz_spectra_map`, `invz_spectra_qpath`, `invz_solve_auto` | Tesla, default `1e-9`; dead band on the longitudinal component `Bz` -- resolved once per driver call (before the `parfor`/solve), zeroing `Bz` below tolerance so `theta_c = 0` reproduces the pure-transverse benchmark exactly |
| `opts.solve_opts` | `invz_spectra_map`, `invz_spectra_qpath` | struct merged into the per-field `invz_solve_auto` opts; fields `J0eff`/`Jxx0`/`hyp` are reserved (driver-owned) -> `invz:solveOpts` if set |
| `opts.forced_moment` | `invz_solve_point_ordered` (set internally by `invz_solve_auto`'s longitudinal route) | when true, bypasses the spontaneous `\|m0\| > m_tol` gate and treats the moment as field-induced; enforces sign-aware seed alignment with the applied `Bz` (one mirrored retry); a non-converged mean-field loop is itself an early-return condition |
| `pt.moment_branch` | `invz_solve_point_ordered` | `'spontaneous'` \| `'field_induced'` \| `'none'` (`'none'` only on the spontaneous-mode paramagnetic early return); present on every return path, `tl = []` on early returns |
| `S.field_dir` | `invz_spectra_map` | the normalized field direction actually used (echoes `opts.field_dir`) |
| `S.Bvec` | `invz_spectra_map` (`[nB x 3]`), `invz_spectra_qpath` (`[1x3]`) | the field vector(s) actually used, dead-band applied |
| `S.Bmag` | `invz_spectra_qpath` | `norm(S.Bvec)` |
| `S.phase = 1` | `invz_spectra_map`, `invz_spectra_qpath` | now means "moment-form self-energy" generically -- spontaneous FM below \(B_c\) at `theta_c = 0`, OR field-induced (forced-moment) under a nonzero tilt; `2` = strict paramagnet, `0` = masked/no solution |

## Deferred follow-ups (design spec §8; not implemented here)

1. **Azimuthal field + two-channel transverse MF.** Requires iterating both
   `hx` and `hy` for every field direction (including the legacy `[Bx 0 0]`,
   where `B64s` induces `<Jy> != 0`), which changes the default x-field model
   and demands benchmark revalidation. Do not resurrect the rejected `By != 0`
   guard (angle-discontinuous, C4-inconsistent) -- see the plan's "Review
   resolutions" table, finding 1.
2. **Full-tensor A0(+A1)** (`odd_implementation_plan.html` Appendix A). **A0**:
   rigorous tensor-RPA layer, `[12,12,nq]` Cartesian x sublattice `J(q)`
   against the full 3x3 `chi0` from `invz_chi0z`; captures the `chi0_xz`/`yz`
   cross-channel exactly, superseding this stage's `theta_c <= 5 deg` scalar
   limit for arbitrary tilt. **A1**: projected 1/z bridge (`chi_tilde0 =
   chi_dom/(1+Sigma_c) + chi_rest`), still a dominant-transition approximation,
   not a rigorous tensor 1/z.

Both are out of scope for this feature (scalar stage, c-tilt only); see
`docs/superpowers/specs/2026-07-16-invz-field-angle-design.md` §8 for the full
discussion and the design decisions that led to descoping them.
