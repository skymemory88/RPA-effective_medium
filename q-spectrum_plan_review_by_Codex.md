# Review of the Demag Knob + q-Path Spectra Plan

**Reviewed plan:** `docs/superpowers/plans/2026-07-14-demag-knob-qpath-spectra.md`  
**Review date:** 2026-07-14  
**Reviewer:** Codex

## Executive assessment

The proposed demagnetization architecture is coherent: the intrinsic ordering coupling (`info.Jcc0` and `Jnu`) is separated from the strict-uniform observable correction (`info.Jshape_cc`), while the transverse applied/internal-field relation is carried by the demag-aware `info.Jaa0`.

The plan is not ready to execute as written, however. The q-path implementation has two numerical blockers that would produce misleading Rønnow et al. Fig. 3 curves:

1. the finite real-space dipole sum is not converged near the Γ-equivalent endpoint `(2,0,0)`; and
2. the proposed frequency window truncates a substantial part of the published dispersion while the peak picker silently reports the grid boundary as a peak.

The plan also overstates the quantitative status of the resulting spectra while known Brillouin-zone, real-axis, and phase-selection problems remain unresolved.

## Findings

### 1. Blocker: the q-path dipole sum is not converged near `(2,0,0)`

Task 5 evaluates every requested path point directly with:

```matlab
Jpath = invz_jq_modes(ion, qpath, struct('dpRng', dpRng, 'cache', true));
Jq = Jpath(:, branch).';
```

This uses the sharply truncated real-space sum in `MF_dipole`. Near a reciprocal-lattice/Γ-equivalent point, the required real-space range grows as the distance to Γ decreases. A fixed `dpRng = 30` therefore does not provide a converged approach to `(2,0,0)`.

A direct check of the largest coupling branch, in meV, gave:

| `h` in `(h,0,0)` | `dpRng=10` | `dpRng=20` | `dpRng=30` | `dpRng=40` |
|---:|---:|---:|---:|---:|
| 1.90  | 0.006676 | 0.006387 | 0.006335 | 0.006317 |
| 1.96  | 0.003995 | 0.006693 | 0.006448 | 0.006293 |
| 1.98  | 0.002279 | 0.004000 | 0.005745 | 0.006705 |
| 1.99  | 0.001746 | 0.002269 | 0.003069 | 0.004003 |
| 1.999 | 0.001560 | 0.001549 | 0.001560 | 0.001571 |
| 2.0   | 0.006439 | 0.006423 | 0.006424 | 0.006422 |

Thus the computed branch can collapse immediately before Γ and then jump at the exact endpoint, where `is_gamma_equiv` adds the Lorentz term explicitly. The proposed three-point unit test (`h = 1, 1.5, 2`) and nine-point smoke test are too sparse to detect this artifact.

#### Required revision

Use a convergent long-range dipolar treatment for path calculations, preferably an Ewald-based Fourier sum with the appropriate direction-dependent nonanalytic Γ limit. If that is deferred, the plan must at least:

- stop claiming general or quantitative q-path support;
- define an explicit direction-aware Γ-limit convention;
- establish `dpRng` convergence at dense points close to Γ;
- reject or mask path points that have not converged; and
- add a regression test over points such as `h = 1.90, 1.96, 1.98, 1.99, 1.999, 2.0`, comparing multiple cutoffs and checking for an artificial endpoint jump.

This issue must be resolved before calling the result a reproduction of Fig. 3.

### 2. Blocker: the driver frequency window clips Fig. 3

The rewritten driver retains:

```matlab
w = (0:0.002:0.45).';
```

Rønnow et al. Fig. 3 extends to about `0.75 meV` near `h = 1` for the 60 kOe curve after the published factor-1.15 scaling. Consequently, the suggested workflow with:

```matlab
Bq = [3.6 4.24 6.0];
```

cannot capture the full dispersion with a `0.45 meV` upper bound.

The proposed `peak_energy` helper compounds the problem:

```matlab
[~, i] = max(c);
E(k) = wm(i);
```

If a resonance lies above the sampled range, this can return the final frequency bin as though it were a valid peak, producing a clipped or artificially flat dispersion without warning.

#### Required revision

- Give q-path mode its own frequency grid, extending to at least `0.8 meV` for the Fig. 3 workflow.
- Treat a maximum in the first or last usable frequency bin as censored, returning `NaN` or issuing a diagnostic rather than reporting it as a peak.
- Require a positive finite spectral maximum and preferably a genuine local maximum.
- Add tests for upper-bound clipping and for `peak_wmin >= max(w)`.
- Consider local interpolation around the peak so displayed dispersions are not quantized at the frequency-grid spacing.

### 3. High: “reproducing Fig. 3” overstates the numerical scope

The plan explicitly leaves the following existing high-impact issues out of scope:

- biased Brillouin-zone quadrature from the closed `[-0.5,0.5]` grid;
- non-causal/negative-weight real-axis spectra; and
- the FM/PM handoff at the bare mean-field boundary rather than the implemented 1/z boundary.

The new q-path engine inherits all three. In particular, `invz_solve_auto` explicitly retains the bare-MF phase selector, and `invz_chi_realaxis` retains the real-axis continuation already known to produce negative spectral weight.

In addition, selecting one sorted coupling eigenvalue and applying the scalar pole formula produces a branch-resolved susceptibility. It does not reproduce neutron scattering intensity, which also depends on eigenvectors, sublattice interference/structure factors, polarization, and the magnetic form factor. Rønnow Fig. 3 itself is an energy-dispersion comparison, so an energy-only workflow is reasonable, but the output should not be presented as the measured neutron intensity.

#### Required revision

Until the inherited blockers are resolved, rename the goal to something like:

> Add an exploratory branch-resolved q-path energy/susceptibility workflow for comparison with the trends in Rønnow et al. Fig. 3.

Documentation and plot titles should distinguish branch susceptibility from neutron scattering intensity and label ordered/near-critical results exploratory.

### 4. High: Task 0 would commit unrelated user work

Task 0 instructs:

```bash
git add -A
git reset -- invz/invz_run_spectra.asv
git commit ...
```

The current working tree contains 32 changed or untracked entries, including unrelated documentation, deleted files, PDFs, MATLAB figures/data, ordered-phase work, and verification scripts. A repository-wide `git add -A` would combine all of that user work into an implementation-plan checkpoint. The plan also commits before running its baseline test, so discovering a red baseline would occur only after mutating history.

#### Required revision

- Run the baseline test suite before any commit.
- Do not automatically commit pre-existing user changes.
- If a checkpoint is genuinely wanted, require explicit user approval of an exact path list and review the staged diff.
- Prefer implementing on top of the dirty tree with explicit-path commits limited to files changed by each task.
- Remove automated `Co-Authored-By` attribution unless it accurately reflects the desired repository attribution policy.

### 5. Medium: the proposed demag documentation contradicts itself

The context correctly states that demag-on `Bc(T)` may shift through the transverse `Jaa0` channel while `Tc(B=0)` remains invariant. The replacement `invz_run_spectra` comment instead says that `info.Jcc0/Jnu`—“and with it Bc and Tc”—is demag-invariant, immediately before noting that q-path spectra still see demag through `Jaa0`.

Those statements cannot both describe the full applied-field boundary.

#### Required wording

- `info.Jcc0`, `Jnu`, and the ordering-channel contribution to criticality are demag-invariant.
- `Tc(B=0)` is exactly demag-invariant because the transverse moment vanishes there.
- `Bc(T)` expressed versus applied transverse field can shift through demag-aware `info.Jaa0`.
- q-path calculations omit the longitudinal strict-uniform `Jshape_cc` transform, but can still change through the transverse applied/internal-field relation.

The `invz_ion` and README wording should use the same distinction. Avoid saying simply that q-path spectra have “no demag” when `Jaa0` remains demag-aware.

### 6. Medium: the Task 2 test snippet is missing `end`

The proposed function:

```matlab
function test_jxx0_override(testCase)
...
verifyGreaterThan(testCase, abs(t1.Delta - t0.Delta), 1e-9);
```

does not include a closing `end`. Because the existing test file consistently terminates its local functions with `end`, the first test run would fail during parsing rather than with the expected “Too many input arguments” error.

Add the missing `end` before executing Task 2.

The test also covers only the paramagnetic two-level helper even though Task 2 changes both `invz_twolevel` and `invz_twolevel_ordered`. Add an ordered-path assertion proving that `Jxx0` reaches `invz_twolevel_ordered` and `invz_solve_point_ordered`.

### 7. Medium: `.hyp=false` remains internally inconsistent

`invz_spectra_qpath` advertises `.hyp` as a supported option and passes it to `invz_solve_auto`. On a paramagnetic real-axis evaluation, however, the proposed `copts` contains `Jsel`, `eta`, and `Jxx0` but not `hyp`; `invz_chi_realaxis` therefore reconstructs a single-ion state with `hyp=true` unconditionally.

The result for `.hyp=false` would use an electronic-only Matsubara medium and an electronuclear real-axis response.

Although the plan lists this older mismatch as out of scope, a new API should not advertise a knowingly inconsistent option.

#### Required revision

Either:

- thread `hyp` through `invz_chi_realaxis` and use it when constructing the paramagnetic single-ion state; or
- remove `.hyp` from the new q-path API and document that only the default electronuclear calculation is supported.

### 8. Medium: broad exception swallowing will conceal implementation defects

`invz_solve_auto` catches every exception from both solvers and discards it. The replacement `one_field` function likewise retains broad, empty catches. A typo, missing field, dimension error, or invalid option would therefore be misreported as “no solution” or a column of `NaN` values.

#### Required revision

- Catch only explicitly expected numerical exceptions, such as the known degenerate-doublet condition.
- Preserve the exception identifier/message in a returned diagnostic.
- Rethrow unexpected programming and API errors.
- Add a test proving that an intentionally invalid option or malformed input is not silently converted into `phase = 0`.

### 9. Low: the “any r.l.u. path” contract is too broad

The proposed path coordinate is:

```matlab
S.s = [0 cumsum(vecnorm(diff(qpath, 1, 1), 2, 2)).'];
```

This is Euclidean distance in raw Miller-index coordinates, not reciprocal-space distance. For tetragonal LiHoF4, equal changes in `h`, `k`, and `l` do not represent equal Cartesian reciprocal distances because `a*` and `c*` differ. The formula is harmless for the specific `(h,0,0)` Fig. 3 path but misleading for the advertised general API.

Either label it explicitly as cumulative distance in index coordinates or compute physical distance from `diff(qpath) * recip_lattice`, optionally returning both coordinates.

Sorted branch number is also not a persistent mode identity through branch crossings. If the API is meant to support arbitrary paths, branch tracking should use eigenvector overlap between adjacent points rather than only the sorted eigenvalue index.

## What is sound in the plan

The following parts are well motivated and can be retained after the blockers above are addressed:

- removing `dm_cc` from `Jnu` and `info.Jcc0`;
- exporting `info.Jshape_cc = 4*dm_cc` for strict-uniform observable conversion;
- exporting demag-aware `info.Jaa0` and passing it as `opts.Jxx0`;
- preserving `ion.Jxx0` as a backward-compatible fallback;
- bumping the q-coupling cache schema;
- testing the algebraic identity
  `chi_int/(1+Jshape*chi_int) = chit/(1-(Jcc0-Jshape)*chit)`;
- testing exact ordering-channel demag invariance with transverse coupling pinned;
- solving the 1/z medium once per `(T,B)` and vectorizing the pole evaluation over path couplings; and
- excluding the low-frequency hyperfine pole when identifying the crystal-field-doublet excitation, provided peak-boundary validation is added.

## Recommended revision order

1. Remove or redesign Task 0 so it does not commit unrelated work.
2. Keep Tasks 1–4 for the demag separation and option threading, fixing the missing test `end`, ordered-path coverage, broad catches, and `hyp` consistency.
3. Before Task 5, implement or select a convergent long-range q-dependent dipolar evaluation and define the directional Γ limit.
4. Add dense q-convergence tests near `(2,0,0)`.
5. Add a q-path-specific frequency grid and peak-censoring tests.
6. Limit the feature claim to branch-resolved exploratory dispersion unless the known BZ, real-axis, and phase-selection blockers are also resolved.
7. Reconcile all demag wording before updating the README and drivers.
8. Run the fast and slow suites, plus explicit cutoff, endpoint, peak-boundary, and Fig. 3 range checks.

## Verification performed during this review

- Current fast MATLAB suite: **32 passed, 0 failed, 8 filtered as slow**.
- Direct `invz_jq_modes` cutoff study at `h = 1.90` through `2.0`, summarized above.
- Comparison with the local primary reference:
  `References/Magnetic excitations near the quantum phase transition in the Ising ferromagnet LiHoF4.pdf`, especially Fig. 3 and its energy range.
- Inspection of the current solver, real-axis, q-coupling, driver, and test call chains.

No implementation files were modified as part of the review.
