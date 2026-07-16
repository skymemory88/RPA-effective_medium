# invz_tensor — full-tensor 1/z (ODD inclusion)

Home of the full-tensor implementation of Jensen's 1/z effective-medium theory
for LiHoF4, carrying the calculation with the complete 3x3 Cartesian (x 4
sublattice) coupling and response tensors — in contrast to `../invz_projected/`,
which projects onto the scalar cc (Ising) channel with a two-level self-energy.

Implementation plan: `../odd_implementation_plan.html` — the staged route is
Appendix A (A0 tensor-RPA parity layer -> A1 projected-1/z bridge -> A2 matrix
effective medium -> A3 tensor 1/z self-energy), with the off-diagonal dipolar
(ODD) blocks of Tiers 1-2 as the physics target. The scalar-stage field-tilt
validation that motivates the tensor route is documented in
`../docs/SESSION-2026-07-16-field-angle.md` (cross-channel yz baseline; scalar
stage supported to theta_c <= 5 deg for peak observables only).
