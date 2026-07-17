# invz_common/

Branch-shared single-ion engine, factored out of `invz_projected/` (pure
move, zero logic changes) so `invz_projected/` and `invz_tensor/` depend on
one copy instead of two diverging ones.

## Scope

16 functions: single-ion CF+Zeeman diagonalization and susceptibility
(`invz_ion`, `invz_single_ion`, `invz_chi0z`, `invz_chiperp`,
`invz_check_transverse_mf`, `invz_cfrot`, `invz_field_vec`, `stevens_ops`),
Matsubara/two-level scalar-Sigma machinery (`invz_matsubara`, `invz_twolevel`,
`invz_g`, `invz_lambdas`, `invz_sigma`), and utilities (`getf`, `invz_const`,
`invz_is_gamma_equiv`) -- signatures unchanged from `invz_projected`.

`invz_cache_key.m` is separate: a weak-hash cache-filename helper for
`invz_tensor` ONLY. Never call it from `invz_projected` code, which keeps its
own local `hash_vec`.

## Tests-of-record

No tests live here. `invz_projected/tests` remains the acceptance suite for
these 16 functions, unchanged function-for-function before/after the move.

## Path requirement

Any caller (test, driver, exploratory script) must `addpath` this folder
explicitly -- it is not implied by having `invz_projected` or `invz_tensor`
on the path.
