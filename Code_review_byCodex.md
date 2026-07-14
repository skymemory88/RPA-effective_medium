# Code Review by Codex: Jensen 1/z Framework and LiHoF4 Implementation

Date: 2026-07-10

## Scope

This review covers:

- `jensen_1z_framework.html`
- `invz/README.html`
- the current LiHoF4 implementation in `invz/`
- the MATLAB tests in `invz/tests/`
- comparison with the local primary references:
  - Jensen, Phys. Rev. B 49, 11833 (1994)
  - Ronnow et al., Phys. Rev. B 75, 054426 (2007)

The review was performed against the dirty working tree as found. No existing
changes were reverted or modified.

## Executive Assessment

The paramagnetic Matsubara core, single-ion diagonalization, and algebraic
self-energy formulas form a credible research prototype. Most equations in
`jensen_1z_framework.html` are transcribed correctly, including the correction
of the printed `m^2 n01` error in Jensen Eq. 2.24.

The quantitative phase-diagram and ordered real-frequency outputs are not yet
reliable production results. The main blockers are:

1. biased Brillouin-zone quadrature;
2. incorrect use of demagnetization in the critical condition;
3. non-causal real-axis spectra whose negative weight is hidden by plotting;
4. a phase selector that switches at the bare mean-field boundary rather than
   the implemented 1/z boundary;
5. production defaults that do not reproduce the headline LiHoF4 benchmarks.

## Findings

### 1. High: Demagnetization is applied incorrectly to criticality

Relevant code:

- `invz/invz_jq_modes.m:69-80`
- `invz/invz_jq_modes.m:103-119`
- `invz/invz_ion.m:14-19`
- `invz/README.html:130-136`

The implementation folds the sample-shape demagnetization field into
`info.Jcc0`, which is then used by the ordering mean field, RPA denominator,
and 1/z critical condition. Ronnow et al. explicitly state that the
demagnetization field cancels from the critical condition because ordering
occurs at wave vectors infinitesimally different from zero, rather than in the
strict uniform `q = 0` mode.

Consequences:

- setting `ion.demag ~= 0` changes `Bc` and `Tc`, contrary to the cited theory;
- the feature test in `test_invz_jq_modes.m:30-53` currently enforces the wrong
  criticality behavior;
- sample shape should be applied to the measured uniform response or external
  field relation, not to the intrinsic phase-instability coupling.

The displayed equation in `invz/README.html:132-134` also omits the factor of
four arising from projection of a scalar-broadcast term onto the four-site
uniform sublattice mode. The code's own comments and test use that factor of
four.

### 2. High: Production Brillouin-zone grids are biased

Relevant code:

- `qVec_generator.m:181-198`
- `invz/invz_run_phase_diagram.m:44-46`
- `invz/invz_spectra_map.m:47-50`
- `invz/tests/test_invz_sigma_crit.m:32-47`

The grid generator uses:

```matlab
qx = linspace(-0.5, 0.5, n);
```

This includes both reciprocal-equivalent zone faces. Equal-weight averaging
therefore double-counts the boundary. For even `n`, it also omits Gamma
entirely. The usual line that removes Gamma has no effect on the default 16^3
grid.

Measured for the current 16^3 production grid:

```text
nq       = 4096
hasGamma = 0
Sigma_c  = 0.272077
Tc(0)    = 1.77945 K
```

The documented extrapolated values are approximately `Sigma_c = 0.2980` and
`Tc(0) = 1.74 K`. The observed `O(1/n)` bias is therefore not solely the
integrable Gamma singularity: it also contains endpoint-weighting error. The
current `2*S(2n)-S(n)` extrapolation is not the exact Richardson combination
for a grid whose interval count is `n-1`.

Use a half-open periodic grid or a midpoint grid, and define explicitly whether
Gamma is sampled or excluded. Recompute all J(q) caches and convergence studies
after this change.

### 3. High: Real-axis susceptibility has negative spectral weight

Relevant code:

- `invz/invz_chi_realaxis.m:32-48`
- `invz/invz_spectra_map.m:85-119`
- `invz/invz_plot_spectra_map.m:15-35`

Production calls to `invz_chi_realaxis` do not pass `Jfull`. Consequently,
`K(w)` is initialized to the static Matsubara value `K(0)` and remains static
at every real frequency. The documented dynamic `K(w)` continuation is not the
default behavior.

More importantly, the resulting diagonal retarded susceptibility violates
spectral positivity. Direct probes gave:

```text
LiHoF4, 16^3 grid, T=0.31 K, Bx=5.5 T:
  9 of 400 positive-frequency samples had chi'' < 0
  minimum chi'' = -27.4796

Near-critical paramagnetic calculation, T=0.31 K, Bx=4.3 T:
  4 negative samples
  minimum chi'' = -40.5829
```

Passing `Jfull` and running the currently untested dynamic-K branch changed the
spectrum but did not eliminate the negative weight in the sampled case.

`invz_plot_spectra_map` converts negative values to `realmin` before taking the
logarithm, visually hiding the causality violation. A plausible colormap is
therefore not evidence that the calculated retarded response is physical.

Required validation should include:

- `chi''(w) >= 0` for every positive frequency in a diagonal channel;
- convergence with respect to `eta`, frequency spacing, and dynamic-K passes;
- a spectral sum rule or Kramers-Kronig consistency check;
- failure rather than plotting when the response violates the chosen tolerance.

### 4. High: Spectra switch phase at the wrong boundary

Relevant code:

- `invz/invz_spectra_map.m:79-100`
- `invz/invz_solve_point_ordered.m:13-17`
- `invz/README.html:121-128`
- `invz/README.html:297-304`

At each field, `invz_spectra_map` accepts any converged bare mean-field ordered
solution before attempting the paramagnetic 1/z solution. This switches the
spectrum at the bare mean-field boundary, not at the boundary defined by

```text
1 + Sigma(0) - J(0)*chi0(0) = 0.
```

Measured on the current default lattice at `T = 0.31 K`:

```text
1/z critical field from invz_critical: Bc = 4.01872 T

At Bx = 4.3 T:
  paramagnetic solve converged
  paramagnetic crit = +0.0686408
  invz_spectra_map selected phase = 1 (ferromagnet)
```

The bare ordered moment remains above the `m_tol` threshold through roughly
4.3 T. The actual branch mismatch is therefore about 0.3-0.4 T in the current
code, not merely a negligible handoff. Statements that the map is continuous
and is FM below the 1/z `Bc` are incorrect.

The correct resolution is Jensen's modified-field/order-parameter
self-consistency (HTML Eqs. 41-43), or an explicitly documented exploratory
branch-selection rule that does not claim to represent the 1/z phase boundary.

### 5. High: Headline LiHoF4 benchmarks are accepted too loosely

Relevant code:

- `invz/tests/test_invz_critical.m:47-55`
- `invz/tests/test_invz_chi_observable.m:26-49`
- `invz/invz_run_phase_diagram.m:49-52`
- `invz/README.html:348-366`

Current measured results include:

```text
Published calculated Hc(0.31 K): 4.24-4.30 T
Current invz_critical result:     4.01872 T

Published Sigma(0) near Hc:      0.0932
Current Sigma(0) at returned Hc: 0.100413

Documented zero-field Tc:        approximately 1.74 K
Production raw-16^3 endpoint:    1.77945 K
```

The critical-field test accepts the range 4.0-4.6 T, so the current result
passes only 0.019 T above the lower limit despite lying outside the published
range. The soft-mode test fixes the field at 4.3 T and accepts 0.10-0.28 meV;
it does not verify that 4.3 T is the implementation's critical field or that
the spectrum is causal.

The 0.22 meV peak is reproducible, but at that point the paramagnetic solver's
critical denominator is positive and the production map chooses the ordered
branch. Peak-position agreement alone does not validate the phase logic.

### 6. Medium: The phase-diagram driver does not run its documented two-regime method

Relevant code:

- `invz/invz_run_phase_diagram.m:1-40`
- `invz/invz_run_phase_diagram.m:54-68`
- `invz/README.html:306-319`

The driver explains that:

- fixed-T field cuts become ill-conditioned near `Tc(0)`;
- the high-temperature boundary should use fixed-B temperature cuts;
- `Ts` should preferably stay below about 1.6 K.

Its current defaults instead are:

```matlab
Ts = linspace(0.05, 1.85, 26);
Bs = [];
```

Thus it runs only the fixed-T method, including the region its own header says
is ill-conditioned. The README still describes a two-regime production run and
also contains stale text about a fixed `[1.0, 2.0] K` temperature window, while
`invz_critical_T` now uses an adaptive window.

The exported `phase_boundary` also omits the separately plotted `[Tc0, 0]`
endpoint.

### 7. Medium: Broad exception handling can convert bugs into physics output

Relevant code:

- `invz/invz_critical.m:49-69`
- `invz/invz_critical.m:78-86`
- `invz/invz_critical.m:110-120`
- `invz/invz_critical_T.m:100-108`
- `invz/invz_spectra_map.m:85-120`

The critical and spectra paths use unqualified `catch` blocks. Dimension
errors, missing functions, invalid parameters, and programmer defects are
therefore treated the same as known numerical conditions such as a degenerate
doublet.

`invz_critical` may additionally call `para_edge`, which returns the lower edge
of numerical paramagnetic convergence as a critical field when no converged
ordered-side value exists. That is a solver diagnostic, not a root of the
physical critical equation.

Catch only expected error identifiers. Unexpected errors should be rethrown.
Fallback estimates should carry an explicit status/result structure and should
not be returned through the same scalar interface as a bracketed physical root.

### 8. Medium: Ordered-phase implementation is incomplete and poorly conditioned in limits

Relevant code:

- `invz/invz_sigma_ordered.m`
- `invz/invz_solve_point_ordered.m:13-17`
- `invz/invz_twolevel_ordered.m`
- `jensen_1z_framework.html:397-440`
- `invz/README.html:369-378`

The README correctly discloses that the implementation omits:

- the full static elastic `xi*h` expression of Jensen Eqs. 39-40;
- the `H_MF <-> H` self-consistency of Eqs. 41-43;
- ordered-phase thermodynamic consistency.

This means the ordered solver is an Option-A approximation, not Jensen's full
ordered-phase theory. In addition, the claim that the doublet is necessarily
degenerate at `Bx -> 0` in both phases is false. The ordered longitudinal field
lifts it, and the current full-lattice ordered solve converged at `Bx = 0`:

```text
is_ordered = true
converged  = true
m0         = 5.51986
M2         = approximately 1.45e-28
Sigma0     = -0.0570171
```

The real concern is the `M2 -> 0` conditioning of expressions containing
`m^2/M2`, not an unavoidable `degenerateDoublet` exception.

`invz_single_ion` also computes a local `converged` flag but does not return it.
Consequently, `invz_solve_point_ordered` can report `pt.converged = true` even
after the mean-field subsolve emitted `invz:mfNotConverged`.

### 9. Medium: Jensen resummation and validity are overstated in the derivation document

Relevant text:

- `jensen_1z_framework.html:329-340`
- `jensen_1z_framework.html:453-470`

The document labels the Dyson resummation as "Derived in full" and states that
higher single-site terms repeat the same pattern. Jensen's paper is more
cautious: it says the unconditional expansion does not estimate sixth and
higher orders well, and that replacing the leading `G0` by `G` is a
Dyson-like prescription found to account for sixth-order terms more
effectively. It is a motivated resummation, not a proof that all higher terms
repeat.

The final claim that the method is quantitatively accurate "essentially
everywhere" is also too strong. Ronnow states that the approach misses the
logarithmic critical corrections and, compared with QMC, overestimates the
small-field suppression of the transition temperature. The ordered formulas
also contain explicit finite-frequency and elastic-sector approximations.

The equation transcription is otherwise strong. In particular, the following
were checked against the local Jensen PDF and agree with the source:

- paramagnetic self-energy Eqs. 2.17-2.19;
- ordered corrections Eqs. 2.26-2.29;
- modified-field relations Eqs. 2.30-2.33;
- free-energy consistency Eq. 2.34;
- the identification of printed Eq. 2.24 as dimensionally and algebraically
  inconsistent with `<Jx> = m*n01`.

### 10. Medium: Documented options and current defaults have drifted

Examples:

- `invz_spectra_map` documents `.hyp=false`, but the paramagnetic real-axis
  path reconstructs `invz_single_ion(..., 'hyp', true)` unconditionally in
  `invz_chi_realaxis.m:24`. The Matsubara and real-axis halves then use
  different Hilbert spaces.
- `invz/README.html` describes a 0.31 K quick-start spectrum, while
  `invz_run_spectra.m` currently defaults to 0.1 K and 201 fields from 0 to 9 T.
- `invz_run_spectra.m:28` sets `eta = 5e-5 meV`, while its comment discusses
  `5e-3 meV`; the selected `eta` is also smaller than the 0.0001 meV frequency
  spacing despite the instruction to keep it at least comparable to the grid
  spacing.
- the README function table mentions an `ordered` spectra-map output, while
  the actual field is `S.phase`.
- the README refers to `info.Jaa0`, which is not produced; only
  `info.Jaa0_dipole` exists.

## Verification Performed

### MATLAB tests

Fast suite:

```text
32 passed, 0 failed, 8 filtered
```

Full suite with `INVZ_SLOW=1`:

```text
40 passed, 0 failed, 0 incomplete
361.6661 seconds
```

The green suite does not cover:

- half-open periodic BZ sampling;
- phase-boundary independence from demagnetization;
- positivity or Kramers-Kronig consistency of the 1/z real-axis response;
- selection of FM/PM branches at the 1/z boundary;
- propagation of mean-field convergence status;
- `hyp=false` consistency on the real axis;
- strict reproduction of the published critical field.

### MATLAB Code Analyzer

Four non-blocking warnings were reported:

- unused `gL` in `invz_const.m`;
- one unused test output in `test_invz_emt_scalar.m`;
- one unused test output in `test_invz_jq_modes.m`;
- one unused variable in `test_invz_sigma_crit.m`.

### Primary-reference comparison

The local Jensen and Ronnow PDFs were converted to text and the relevant
equations were also visually inspected in the rendered PDF pages. This
confirmed both the strong equation transcription and the discrepancies noted
above regarding resummation strength, demagnetization, and quantitative
validity.

## Recommended Remediation Order

1. Replace the BZ generator used by `invz` with a half-open or midpoint periodic
   quadrature and invalidate all J(q) caches.
2. Separate intrinsic `q -> 0` critical coupling from sample-shape correction
   to the strict uniform observable.
3. Establish a causal real-axis continuation with positivity and sum-rule
   regression tests; do not hide negative spectral weight.
4. Make spectra branch selection consistent with the 1/z critical boundary,
   or complete the modified-field self-consistency.
5. Recompute `Tc(0)`, `Hc(0.31 K)`, `Sigma(0)`, and the soft-mode spectra with
   tightened published-value tests.
6. Repair driver defaults and synchronize `invz/README.html` with the actual
   APIs and run configurations.
7. Replace broad catches with identifier-specific numerical handling and
   propagate convergence status through every solver layer.

Until items 1-5 are resolved, phase diagrams and ordered spectral maps should
be labeled exploratory rather than quantitative LiHoF4 predictions.
