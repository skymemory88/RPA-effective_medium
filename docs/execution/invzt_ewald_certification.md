# Full-tensor Ewald default certification

Status: `certified_default`

Branch: `invzt-multilevel-1z`

Base: `29d4117e55ba0a2e25543f5a2a19dea195169741`

Primary temperature for solver anchors: `0.1 K`

## Scope and decision

This work package decides only whether the full-tensor lattice builder may use
the conducting, exact-`k=0`-omitted Ewald sum by default. It does not claim that
Ewald repairs the scalar-A1 fixed point, certify a multilevel closure, or change
the projected-spin implementation.

The candidate production controls are inherited from the frozen primitive
calibration:

```text
alpha0 = sqrt(pi) / Vc^(1/3)
r_cut = 5.5 / alpha0
g_cut = 11 * alpha0
boundary = conducting_k0_omitted
```

For `invz_ion()` these evaluate to `alpha0=0.268430511437721 1/Angstrom`,
`r_cut=20.4894740562161 Angstrom`, and `g_cut=2.95273562581493 1/Angstrom`.
The exploratory `(0.3,16,3)` controls are not a production-default candidate.

Brute force must remain explicitly selectable as a diagnostic reference.

## Evidence inherited from the certified primitive

The primitive was implemented and reviewed in commits `4c86c83` through
`adf2d63`. Its frozen Gate-A calibration used:

- alpha multipliers `{0.6,0.8,1.0,1.2,1.5}` with matched generous cutoffs;
- real-cutoff ladder `Cr={4.5,5.0,5.5}` at `Cg=13`;
- reciprocal-cutoff ladder `Cg={9,10,11}` at `Cr=6.5`;
- candidate default `(Cr,Cg)=(5.5,11)` versus joint refinement `(6,12)`;
- complete-tensor metric
  `|A-B| <= 1e-8*Tscale + 1e-8*max(|A|,|B|)`;
- coupling metric
  `|JA-JB| <= 1e-8*Jref + 1e-8*max(|JA|,|JB|)`, with
  `Jref=0.006424435656 meV`.

The current primitive may inherit those numerical results only if its
executable statements are identical to the reviewed version. Comment-only
differences are allowed and must be demonstrated mechanically.

## Pre-registered gates

### G1 — primitive identity and safety

1. Executable-statement comparison with `adf2d63` is identical after removing
   comment-only and blank lines.
2. The primitive's boundary guard, alpha/cutoff guards, geometry fingerprint,
   candidate-count caps, and 4-GiB pre-allocation memory cap remain present.
3. Any executable difference stops inherited certification and requires the
   original primitive suite to be rerun.

### G2 — tensor assembly

The focused tensor contracts must pass:

- all nine Cartesian blocks match direct primitive-plus-exchange assembly;
- every full `12x12` page is Hermitian within Frobenius residual `2e-12`;
- Ewald Gamma receives no second wrapper Lorentz term;
- `full-dipole` is exchange-only;
- reciprocal-shifted pages have sorted-eigenvalue difference at most `2e-11`;
- absent and explicit brute-force requests remain bit-identical before the
  default flip.

### G3 — tensor parameter convergence

On the frozen 30-point primitive sample and the `16^3` half-open production
grid, the candidate default and joint refinement `(Cr,Cg)=(6,12)` must pass
both frozen metrics above for:

- every component of the full `12x12` coupling pages;
- sorted page eigenvalues;
- uniform `aa` and `cc` Gamma projections.

The full-grid comparison is cache-disabled and records maximum absolute errors,
metric ratios, wall time, candidate counts, and estimated peak bytes.

### G4 — backend identity and cache safety

1. Cache-disabled and cache-miss/hit outputs are bit-identical.
2. Cache identity includes backend, exact Ewald controls, q-grid and weights,
   lattice/basis/couplings, convention, and output shape.
3. Explicit-q real-axis reconstruction uses the converged point's exact backend
   and controls. Backend mixing must fail closed.

### G5 — converged solver anchor

At `T=0.1 K`, `B=[6,0,0] T`, and the `16^3` half-open grid, both the candidate
and joint-refined lattices must produce converged PM states. The comparison
must satisfy:

- `abs(Sigma0_default-Sigma0_refined) <= 1e-7`;
- `abs(crit_default-crit_refined) <= 1e-7`;
- maximum absolute difference in the returned `Sigma`, `K`, and local-medium
  `G` arrays `<= 1e-7`.

This is a parameter-stability gate, not a brute/Ewald-equivalence gate and not
a low-field convergence claim.

### G6 — default flip and regression

Only after G1--G5 pass:

1. absent backend and explicit `dipole='ewald'` both resolve to the calibrated
   Ewald controls;
2. the production spectra driver requests the same calibrated default;
3. explicit `dipole='bruteforce'` retains the historical `dpRng` behavior;
4. provenance and cache namespace make the two backends impossible to mix;
5. focused tests, MATLAB Code Analyzer, and `git diff --check` pass.

## Diagnostic, non-gating comparisons

A brute-force `dpRng={10,20,30}` ladder will be recorded for Gamma projections
and selected finite-q eigenvalues. Because dipolar sums and their Gamma limits
are boundary-convention sensitive, proximity to a finite brute cutoff is not a
scientific acceptance test for Ewald. The comparison is retained to expose the
size and direction of the production-default change.

Grid convergence of a future multilevel physical observable is a separate
closure-level obligation. It is not inferred from convergence of the lattice
sum at each grid point.

## Stop conditions

Any failed G1--G5 gate leaves brute force as the full-tensor default and records
the first failure. Tolerances and samples are not changed after results are
observed. Passing all gates authorizes the Ewald default only for the isolated
full-tensor branch.

## Executed result — 2026-08-02

Decision: **all gates passed; calibrated Ewald is the full-tensor default.**

### G1

After removal of comment-only and blank lines, both the reviewed `adf2d63`
primitive and the current `invz_dipole_ewald.m` have executable SHA-256
`62051b3ec54eb63986cec8517e872b551a8277bced2d14c7039e061df1f8021d`.
The safety constants and failure paths are unchanged.

### G2 and G4

Nine focused MATLAB contracts pass. They cover direct all-block assembly,
Hermiticity, Gamma/no-double-Lorentz behavior, exchange-only subtraction,
reciprocal-shift spectra, option validation, implicit-default identity,
explicit brute reachability, cache miss/hit identity, corrupted-cache
recomputation, exact real-axis backend forwarding, and the legacy-point
brute-force fallback.

### G3

Candidate versus joint refinement results:

| Sample | max page error (meV) | tensor tolerance ratio | max eigenvalue error (meV) | coupling tolerance ratio |
|---|---:|---:|---:|---:|
| frozen 30 q | `2.283e-15` | `6.46e-5` | `6.224e-15` | `8.50e-5` |
| half-open 16^3 | `1.860e-15` | `3.58e-5` | `3.875e-15` | `5.38e-5` |

The two cache-disabled production-grid builds took `1.524 s` and `1.525 s`.
Candidate/refined preflight estimates were `5,912,832` and `7,999,960` bytes,
well below the primitive's 4-GiB cap.

### G5

Both 6 T PM anchors converged in five iterations. Differences were:

```text
Sigma0       5.14e-14
crit         2.77e-13
max Sigma    5.14e-14
max K        5.33e-15
max G        1.01e-11
```

The implicit production default independently reproduced Ewald provenance,
`Jcc0=0.00642166181275534 meV`, `Sigma0=0.053348176938097`, and
`crit=0.23581987140837`.

### G6

`invzt_jq_tensor` now resolves an absent backend to calibrated Ewald;
`invzt_run_spectra`, `invzt_run_phase_diagram`, and `invzt_run_ladder` request
and record it explicitly. An Ewald request carrying `dpRng` fails rather than
silently ignoring the brute-force control. Explicit `dipole='bruteforce'`
retains the historical path. The real-axis continuation reconstructs from the
point's provenance; a pre-provenance legacy point explicitly falls back to
brute force.

MATLAB Code Analyzer reports zero findings on all nine modified or added MATLAB
files. `git diff --check` passes.

### Non-gating brute-force diagnostic

| dpRng | Jcc0 (meV) | Jaa0 (meV) | max selected-q eigenvalue difference from Ewald (meV) |
|---:|---:|---:|---:|
| 10 | `0.006439382606` | `0.003502972730` | `1.114e-4` |
| 20 | `0.006423040501` | `0.003511143782` | `2.294e-5` |
| 30 | `0.006424435656` | `0.003510446205` | `8.009e-6` |

These finite-cutoff differences are retained as diagnostics and were not used
to weaken or replace the Ewald self-convergence gates.

### Artifact

`docs/diagnostics/invzt_ewald/certification.mat` is 1,672 bytes with SHA-256
`0ff0f30340639690f9345570e4cb55f8cd90245b1cf65fed8b455e7d6d278e5d`.
