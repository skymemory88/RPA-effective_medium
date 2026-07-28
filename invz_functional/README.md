# Strict scalar functional prototype

This directory is an isolated, non-production implementation of the strict `O(1/z)` scalar
ring functional derived in `../invzp_functional_wp0_ring_derivation.md`.

It deliberately contains no call from `invz_projected/`, `invz_common/`, or `invz_tensor/`.
The retained objects are:

- an exact source-biased scalar two-level local reference;
- its positive connected Matsubara `C2`, source derivatives, and beta derivative;
- the finite-mode tadpole-subtracted ring functional and its analytic source derivatives;
- a rigorous scalar free-energy tail bound and successive-cutoff convergence audit;
- the independently varied `(m,h)` common functional;
- an interval-bounded stationary-root enumerator that compares stable roots by that functional;
- a stationary cutoff audit that reruns complete root enumeration rather than only re-evaluating a
  warm-started branch;
- dense exact finite-cluster and two-site oracles (currently capped at eight scalar sites).

`invzf_centered_pair_exact` and `invzf_centered_cluster_exact` hold the uncoupled local moment fixed
while scanning centred bonds.  The pair cubic coefficient isolates the leading `C3-C3` skeleton
topology; the three-site mixed quartic coefficient isolates ring plus one connected `C4` and grades
the cancellation of reducible `C3*C2*C3` content after 1PI conversion.  These are exact-cluster
diagnostics, not physical source prescriptions.

`invzf_projected_inputs` and `invzf_mode_grid_audit` form a read-only diagnostic bridge to the
production transverse doublet and BZ coupling spectrum.  They do not install a production dispatch:
the bridge uses a fixed electronic doublet, treats `Bz` only as its projected longitudinal source,
and explicitly excludes hyperfine, ODD/retarded modes, and the ordered `tanh/xi` replacement.
For the future nonlocal skeleton, the optional detailed output of `invz_bz_couplings` supplies
q-resolved Hermitian `Jcc` pages and normalized row weights; ordinary production callers retain their
existing scalar-output path.

`invzf_electronuclear_local` is the next isolated local oracle.  It uses the full source-biased
electronuclear Hamiltonian and a stable Lehmann correlator, with nested source/beta stencils for the
ring derivatives.  Its first version deliberately requires `transverse_mf='none'`: admitting a
self-consistent transverse molecular field requires adding its moment and conjugate field to the
common functional, not patching only the static susceptibility.
`invzf_electronuclear_inputs` safely separates an applied `Bz` into the common-functional external
source while passing only the transverse field to that local oracle.

`invzf_local_1pi_static` performs the first WP4 amputation gate for the static scalar
`gamma2/gamma3/gamma4`.  `invzf_twolevel_cumulant` and `invzf_twolevel_1pi_vertex` extend that gate
to exact frequency-labelled two-level `C2/C3/C4` and `gamma2/gamma3/gamma4`, using ordered-simplex
matrix exponentials for analytic KMS/Hermite limits.  `invzf_twolevel_vertex_table` freezes those
objects on a signed convolution grid.

The first varied-covariance ansatz in `../invzp_functional_wp4_skeleton_spec.md` was implemented
temporarily and rejected by its immutable mixed-chain coefficient before commit.  Its Gaussian local
trace plus 1PI `gamma4` double bubble uses the wrong local bilocal kernel and changes the retained
`a^2*b^2` coefficient under stationarity.  `invzf_local_bilocal_hessian_static` is the next, smaller
oracle: it derives the exact `Q=0` local Hessian from the connected `C2/C3/C4` response block and
recovers the immutable zero-source static coefficient.
`invzf_local_bilocal_hessian_modes` extends that response block to nonnegative Matsubara
representatives, with the distinct zero- and nonzero-frequency pair-disconnected factors and the
nonzero-source `C3-C3` Schur term.  It remains a finite-mode local oracle: the Matsubara-tail
completion, lattice trace, and double counting are not derived.  No skeleton solver remains in this
directory; a corrected route must pass those gates rather than patching the failed ansatz.

It does not contain EMT, a dressed self-energy, the ordered `tanh/xi` replacement, spectral
continuation, or any production branch solver.
