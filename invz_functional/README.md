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
- dense exact finite-cluster and two-site oracles (currently capped at eight scalar sites).

It does not contain EMT, a dressed self-energy, the ordered `tanh/xi` replacement, spectral
continuation, or any production branch solver.
