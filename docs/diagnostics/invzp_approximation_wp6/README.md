# Controlled missing-area approximation evidence

This directory contains reproducible diagnostics and retained MATLAB outputs
for the explicitly opt-in missing-area approximation. The strict
`full_profile` default and the rigorous-route conclusions are unchanged.

## Current production evidence

- `invzp_approximation_census.m` inventories certified terminal components and
  algebraic missing-area support across 0.2--1.2 T and 2.5--3.0 T.
- `invzp_branch_area_ensemble_validation.m` checks the declared numerical
  branch against independently continued exact sheets.
- `invzp_approximation_production_validation.m` propagates the factor
  0.75/1/1.5 ensemble through the real-axis susceptibility and checks the
  strict-default oracle.
- `invzp_approximation_resolution_census.m` keeps profile-node uncertainty
  separate from missing-area sensitivity.
- `invzp_adjacent_retry_map_smoke.m` verifies the two-sided 129-node recovery
  of the 0.36 and 0.45 T cold masks for every ensemble member.
- `invzp_adjacent_retry_highfield_smoke.m` verifies that the two-sided retry
  cannot bridge the ordered/paramagnetic boundary.
- `invzp_ordered_boundary_retry_smoke.m` verifies the separately gated
  two-lower-source recovery of the central 4.68 T member and preserves the
  incomplete noncentral members and the 4.77 T PM labels.

## Active high-field sliver evidence

- `invzp_highfield_sliver_convergence_audit.m` reproduces the moving
  4.663636 T mask on the 111-point driver grid. It separates cold component
  support from node resolution, damping, missing-area factor, and seeded
  continuation effects.
- `invzp_highfield_one_sided_seed_audit.m` compares two independent lower-field
  seeds at 4.663636, 4.68, and 4.70 T and records root, component, PM-mass,
  stability, derivative, and seed-agreement diagnostics.
- `invzp_highfield_derivative_section_audit.m` reconstructs narrow fixed-field
  sections around the 4.663636 and 4.68 T continued roots and checks direct
  moment derivatives, exact analytic identities, and local/coarse secants.
- `invzp_highfield_boundary_continuation_ladder.m` repeats direct two-seed
  continuation at 129 and 257 nodes over 4.64--4.71 T without reusing a
  recovered target as a seed.
- `invzp_highfield_progressive_field_continuation.m` advances two independent
  source histories on a fine ladder through 4.7188 T at both resolutions.
- `invzp_highfield_boundary_retry_contract_audit.m` checks the current 101-grid
  4.68 T target against direct two-source, independent-PM, stability,
  resolution, and real-axis gates.

The earlier 4.70/4.71 loss was a continuation-distance failure, not an observed
branch endpoint. With accepted intermediate states used only inside a
diagnostic trace, both independent histories and both resolutions reach
4.7188 T with positive margins and residuals below (6.3\times10^{-9}).
At the production 4.68 T target, direct untouched 4.50 and 4.59 T sources agree,
the independently damped PM solve converges with negative mass, and the
129/257 peak shift is 0.00110 GHz with a 0.21% integrated-weight change. A
distinct opt-in ordered-boundary retry is therefore enabled in the spectra
driver. It requires two frozen lower ordered sources, an accepted PM field
above, a separately converged negative PM mass, single-component agreement,
and (D_{\rm uni},F'>10^{-3}). It never retries `phase==2` or uses a recovered
state as a seed. Only factor one passes at 4.68 T; the 0.75/1.5 sensitivity
members remain honestly incomplete.

## Exact-zero endpoint evidence

- `invzp_zero_field_limit_ladder.m` compares the raw electronic doublet,
  basis-invariant full-ion observables, coupled states, physical margins, and
  real-axis spectra on a 19-point logarithmic/linear ladder to exact zero.
- `invzp_zero_field_ensemble_validation.m` evaluates all three missing-area
  members end-to-end at exact zero and four positive control fields, including
  the production 0--6 GHz real-axis grid and the strict-mode comparator.
- `invzp_zero_field_resolution_audit.m` retains the 129/257-node profile
  systematic separately from the missing-area member spread.

The exact-zero approximation is now complete without copying a positive-field
spectrum or inserting a transverse epsilon. The singular predictor is rejected
with `twolevel_domain_invalid`, while the selected positive-\(h\) component and
final ordered state pass every ordinary gate. All ensemble spectra are finite
and continuous to \(10^{-6}\) T; no member has an interior peak in the retained
0--6 GHz window, so peak energy is honestly censored. Strict `full_profile`
remains masked. The former high-field central sliver is handled only by the
separately labeled ordered-boundary retry; its area-sensitivity interval is
still incomplete at that field.
