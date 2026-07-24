# Ewald reformulation — convergence / acceptance PRE-REGISTRATION (FROZEN)

**Status: FROZEN 2026-07-24.** This is the pre-registration the reviewed design
(`docs/invzp_ewald_design.md`, §5/§6/§8/§10) requires as Step 3. Every value, metric, and gate below
is now fixed **before** any Ewald-based physical result is inspected. One controller correction was
applied at freeze — the §3 Gate-A test-5 cell-phase covariance conjugation direction
(`U'*Dcell*U`, numerically verified; see that line). Any later change requires a dated prospective
re-registration, not an in-place edit. Nothing here authorizes production code, default changes,
ODD/tensor adoption, or Phase-3/3A/3B/4 work. It governs Step-4/5/6/7 implementation and acceptance
only.

**Provenance of the frozen numbers:** the primitive convention is `docs/invzp_ewald_design.md` §3
(re-derived in `docs/invzp_ewald_derivation.md` §A/§B, finite-difference-confirmed). The parameter
values are the scratch calibration in that note's Verification/Calibration A–E. **Per design §4.3/§5.2
those scratch results are calibration only; the retained repository tests (Gates A–C) rerun every
check after this freeze.** J_ref, gfac, Vc, lorz below are the live values used throughout.

Frozen reference constants (live `invz_ion()`/`invz_const()`): `gfac = 0.08388`,
`Vc = 287.8917187 Å³` (`Vc^{1/3} = 6.6030`), `4π/Vc = 0.0436496425424 Å⁻³`,
`lorz = 4π/(3Vc)·gfac = 1.22044400549e-3 meV`, `J12 = −1e-4`, `J_ref = |Jcc0| = 0.006424435656 meV`,
physical anchors `Jcc0 = 0.00642444 meV`, `Jaa0 = 0.00351045 meV`.

---

## 1. Frozen primitive convention (fixed; not re-opened here)

The primitive `[dip, counts, geom] = invz_dipole_ewald(q, a, tau, eopts, geom)` returns, in
`MF_dipole`'s quantity/sign/units (`dip = −Σ'_R T(R+d) e^{−iq·(R+d)}`, Å⁻³, `[3,3,4,4]`,
`dip(:,:,m,n)=conj(dip(:,:,n,m))`), the sum of the three convergent parts with the **exact `k=0`
reciprocal term omitted** (`boundary='conducting_k0_omitted'`, no primitive surface term):

```
T_ab(x) = (3 x_a x_b − |x|²δ_ab)/|x|⁵ = +∂_a∂_b(1/|x|)          x = R + d_nm, d_nm = τ_m−τ_n
dip^(r)  = −Σ'_{|x|≤r_cut} g_ab(x) e^{−iq·x},
           g_ab = P(r)(3x_ax_b−r²δ)/r⁵ + Q(r)x_ax_b,
           P(r)=erfc(αr)+(2αr/√π)e^{−α²r²},  Q(r)=(4α³/√π)e^{−α²r²}/r²
dip^(G)  = +(4π/Vc) Σ_{k≠0,|k|≤g_cut} (k_a k_b/|k|²) e^{−|k|²/4α²} e^{+iG·d_nm},   k = q+G
dip^(self) = −δ_nm δ_ab · 4α³/(3√π)
```

Q-domain contract (design §2.1): reduce every input row to the half-open reciprocal cube
`qbar = q − K`, `K = floor(q+1/2)`, `qbar_i ∈ [−1/2,1/2)`; evaluate at `qbar`; restore the gauge
`dip_nm(q) = e^{−i 2π K·(τ_m−τ_n)} dip_nm(qbar)`. The cached reciprocal candidate set is the union
of all `G` with `|qbar+G| ≤ g_cut` over the **entire** half-open cube. The `geom` fingerprint
includes `a`, `tau`, `alpha`, `r_cut`, `g_cut`, the boundary label, the q-reduction convention, and a
schema version.

`counts` is a diagnostic struct with this frozen schema:

```text
counts.real_pair                         [ntau x ntau] actual retained real vectors
counts.recip_candidates                  scalar cached reciprocal-candidate count
counts.recip_used                        [nq x 1] q-specific retained reciprocal counts
counts.preflight.real_cube_bound         [ntau x ntau] conservative integer-mesh bounds
counts.preflight.recip_cube_bound        scalar conservative integer-mesh bound
counts.preflight.estimated_peak_bytes    scalar, including the 25% margin
counts.preflight.array_manifest          struct array
```

Each manifest row contains `name`, `class`, `is_complex`, `size`,
`numel`, and `bytes`; it covers retained and temporary arrays included in the
preflight estimate. `counts` is never collapsed to a legacy scalar.

## 2. Frozen parameters (α, cutoffs, brackets, guards, resource caps)

- **Default splitting parameter:** `α₀ = √π / Vc^{1/3} = 0.268431 Å⁻¹` (deterministic rule
  `α₀ = √π / Vc^{1/3}`).
- **Real / reciprocal cutoffs:** `r_cut = C_r/α₀`, `g_cut = C_g·α₀`. **Production defaults `C_r = 5.5`,
  `C_g = 11`** ⇒ `r_cut = 20.49 Å`, `g_cut = 2.953 Å⁻¹` (default self-convergence residual: raw-tensor
  3.6e-13, coupling 2.1e-13 — §3 metric).
- **Cutoff ladders (≥3 rungs, tested on SEPARATE axes per design §5.1):**
  - real axis: `C_r ∈ {4.5,5.0,5.5}` at fixed generous `C_g=13`;
  - reciprocal axis: `C_g ∈ {9,10,11}` at fixed generous `C_r=6.5`; and
  - joint refinement: production default `(C_r,C_g)=(5.5,11)` versus `(6.0,12)`.
  Both separate ladders must bracket and clear the §3 convergence tolerance; the joint comparison
  must also pass it.
- **α bracket:** multipliers `{0.6,0.8,1.0,1.2,1.5}·α₀`
  (`[0.161,0.403] Å⁻¹`). For α-independence, use α-matched generous cutoffs
  `r_cut=6.5/α`, `g_cut=13α`, so the dimensionless truncation bounds stay fixed.
- **Accuracy guards (hard):** reject/raise unless `α·r_cut ≥ 4.5` **and** `g_cut/(2α) ≥ 4.5` (the
  coarsest calibrated ladder rung; below these the truncation tails exceed 1e-8).
- **Deterministic resource caps (hard, frozen; evaluated BEFORE any large allocation):** max
  real-space pair vectors `≤ 3.0e6` per sublattice pair, max reciprocal candidates `≤ 3.0e6`, and
  estimated peak working bytes `≤ 4 GiB`. Exceeding one is a hard failure. The estimate must include
  retained real/reciprocal geometry, q-specific work arrays, output tensors, and a 25% allocator
  margin. These conservative caps are intended to catch a runaway enumeration before allocation.
  Before constructing an integer mesh, use the following conservative cube bounds. With direct
  lattice matrix `a`, pair displacement `d`, and
  `nmax_r=ceil((r_cut+norm(d))/svd_min(a))`, the preflight real candidate bound is
  `(2*nmax_r+1)^3` for that pair. With Cartesian reciprocal row matrix
  `B=2π*inv(a)'`, let `qmax` be the largest `norm(qbar*B)` over the eight canonical-domain
  corners, set `nmax_G=ceil((g_cut+qmax)/svd_min(B))`, and bound the reciprocal integer mesh by
  `(2*nmax_G+1)^3`. These conservative bounds, not a partially allocated candidate array, are
  checked against the vector-count caps. The byte estimate is the sum, over every planned retained
  and temporary numeric/logical array, of `number_of_elements × MATLAB-class-byte-width`, multiplied
  by `1.25`; before enumeration, any unknown retained count is replaced by its conservative cube
  bound. This initial conservative estimate is the value retained in `counts` and used for the
  hard gate; it is not replaced by a post-allocation measurement.
- **Operational timing target (not a scientific gate):** report elapsed wall time for every
  configuration. `120 s` per single-q primitive configuration on the development machine is a
  scheduling target only. Because wall time cannot be preflighted exactly and depends on hardware,
  load, MATLAB version, and thread count, exceeding it pauses/reports the run but does not invalidate
  a numerically converged result. Any enforced timeout must record that environment and is
  classified separately from Gate-A convergence.

## 3. Frozen primitive-validation metrics + tolerances (Gate A) — explicit definitions

Per design §5.2, every metric is defined here (norm, absolute floor, relative denominator, direction,
near-zero rule). All comparisons are **symmetric** in the two operands.

- **Tensor scale (per q and operand pair):**
  `T_scale(q;A,B) = max over all 3·3·ntau·ntau components of max(|A(q)|,|B(q)|)`. For a comparison
  to a converged reference, `B=dip_ref`. This symmetric scale is used only for complete tensors.
- **Tensor agreement `M_T`:** two tensors `A, B` agree iff, componentwise,
  `|A−B| ≤ AbsTol_T + RelTol_T·max(|A|,|B|)` with
  **`AbsTol_T = 1e-8·T_scale(q;A,B)`**,
  **`RelTol_T = 1e-8`**. Report the single worst component's `|A−B| − (AbsTol_T+RelTol_T·max)`; ≤0
  everywhere = PASS. The `AbsTol_T` floor is the frozen near-zero-component rule.
- **Screened-Hessian finite-difference metric `M_HFD`:** for one real-space sample `x`, let
  `H_scale(x;A,B)=max(|x|⁻³, max_ab|A_ab|, max_ab|B_ab|)`. The Richardson-extrapolated numerical
  Hessian `B` agrees with the closed tensor `A` iff
  `max_ab |A_ab−B_ab|/H_scale ≤ 1e-7`. A component is sign-gated only when
  `max(|A_ab|,|B_ab|) ≥ 1e-8·H_scale`; every such component must have the same sign.
- **Scalar-oracle finite-difference metric `M_FD`:** a complete oracle tensor `B_FD` agrees with the
  primitive tensor `A` iff, componentwise,
  `|A−B_FD| ≤ AbsTol_FD + RelTol_FD·max(|A|,|B_FD|)`, with
  **`AbsTol_FD=2e-8 Å⁻³`** and **`RelTol_FD=1e-7`**. This deliberately looser, dimensioned floor
  covers the calibrated finite-difference error; it is not `M_T`.
- **Coupling agreement `M_J`** (for `Jnu[nq×4]`, `Jcc0`, `Jaa0`, `Juni`): agree iff
  `|J_A−J_B| ≤ AbsTol_J + RelTol_J·max(|J_A|,|J_B|)`, **`AbsTol_J = 1e-8·J_ref`**, **`RelTol_J=1e-8`**;
  worst over all branches/q reported.
- **Exact-identity tolerance `M_id`** (machine-level structural identities):
  `max_components|A−B| ≤ 1e-12·T_scale(q;A,B)`.

### Frozen deterministic samples

Unless a test explicitly names the hand-built fixture below, use the production
`invz_ion()` lattice and basis.

- Generic reduced-coordinate wavevectors:
  \[
  q_{\rm int}\in\{
  (0.137,0.291,0.453),\
  (-0.311,0.173,-0.227),\
  (0.25,0,0.1)\}.
  \]
- Candidate-boundary probe: set \(\delta_q=2^{-40}\), and test all eight
  Cartesian-product corners whose components are independently chosen from
  \(\{-0.5,\ 0.5-\delta_q\}\).
- Face probes: the six points obtained by fixing one component in turn to
  either \(-0.5\) or \(0.5-\delta_q\), with the other two components equal to
  \(0.137\) and \(0.291\), respectively.
- Edge probes: the twelve points obtained by choosing one component in turn
  as the free component \(0.173\), with the other two independently chosen
  from \(\{-0.5,\ 0.5-\delta_q\}\).
- Define \(Q_A\) as exact Γ together with all generic, face, edge, and corner
  probes above.
- Reciprocal translations:
  \[
  K\in\{(\pm1,0,0),(0,\pm1,0),(0,0,\pm1),
  (1,-2,3),(-2,1,-1)\}.
  \]
- Common fractional-origin shift: \(s=(0.137,-0.211,0.089)\).
- Representative per-basis Bravais shifts: for every basis index in turn,
  repeat the calculation after adding each of \((1,0,0)\) and \((0,-1,1)\)
  to that representative while leaving all other representatives unchanged.
- Every identity test covers all basis pairs and all nine tensor components.

Generate the screened-Hessian sample once in MATLAB by
`rng(20260724,'twister')`, followed by `U=randn(200,3)`, row normalization of
`U`, `p=randperm(200)'`, and
`r=exp(linspace(log(0.4),log(6),200))'` in Angstrom. The sample vectors are
`x=U.*r(p)`. Use these same 200 vectors at every frozen positive-\(\alpha\)
value. The \(\alpha=0\) closed-form test is separate.

The explicit small-shell sign/phase fixture is:

- \(a=\operatorname{diag}(3.1,4.2,5.3)\) Angstrom;
- \(\tau=\{(0,0,0),(0.23,0.31,0.17)\}\) in fractional coordinates;
- \(\alpha=0.27\ {\rm Angstrom}^{-1}\);
- \(r_{\rm cut}=18.0\) Angstrom; and
- \(q=(0.17,-0.23,0.29)\) in reduced reciprocal coordinates.

For this fixture, directly enumerate every integer-cell vector whose physical
displacement satisfies the stated real-space cutoff. This direct finite shell,
not production candidate generation, is the reference for the Gate-A
real-space sign/phase comparison.

Gate-A tests (all retained, rerun post-freeze):
1. **Screened Hessian:** `g_ab` vs central-difference `∂_a∂_b[erfc(αr)/r]` at Cartesian steps
   `h ∈ {4e-3,2e-3,1e-3} Å`, with O(h²) Richardson removal from both adjacent step pairs, over
   the 200 frozen `x` samples per α across the α bracket. The two Richardson estimates must differ by
   `≤1e-7·H_scale`; the finer estimate must pass `M_HFD`, with 0 sign mismatches on sign-gated
   components. A pre-freeze double-precision arithmetic check on this exact sample gave worst
   normalized analytic and adjacent-Richardson errors `6.08e-9` and `5.99e-9`, respectively;
   these are calibration only and the retained test reruns them after freeze. Evaluate the closed
   formula separately at `α=0` and require the bare tensor exactly.
2. **Signs vs explicit small screened real-space sum:** use the frozen hand-enumerated fixture above
   and require `M_T`.
3. **Reciprocal phase + gauge** (`M_id`):
   `dip_nm(q+K) = exp[−i2πK·(tau_m−tau_n)]dip_nm(q)` for every frozen
   `q in Q_A` and `K`,
   and sorted-eigenvalue periodicity at the **Phase-1 item-2** tolerance (`AbsTol_J=1e-10`,
   `RelTol_J=1e-8`) — NOT raw matrix equality.
4. **Canonical half-open reduction + extended-zone gauge restoration + reciprocal-candidate
   completeness** at all eight frozen BZ corners and every frozen integer shift (`M_id` on the
   restored dip; the candidate set must contain every needed `G`).
5. **Structural identities and basis representation:** at every `q in Q_A`, require sublattice-pair
   conjugation/Hermiticity at `M_id`; at exact Γ also require realness at `M_id`. Apply the frozen
   common origin shift and require raw invariance at `M_id`. For each frozen individual
   representative shift, explicitly re-index the real-space cells and require raw invariance in
   this primitive's total-displacement gauge. As a separate convention check, convert to the
   cell-phase gauge
   `Dcell_nm=exp[+i2πq·(tau_m−tau_n)]dip_nm` and require
   `Dcell'=U'*Dcell*U` at `M_id`, where
   `U_nn=exp(+i2πq·L_n)` and `L_n` is that representative's Bravais shift.
   (Conjugation direction verified numerically against `MF_dipole`: a Bravais
   shift of representative n multiplies cell-phase row n by `exp(−i2πq·L_n)`, so
   the reverse `U*Dcell*U'` fails by O(1) — 8.2e-2 vs 1.0e-6 at the finite-N
   floor.) This is the precise interpretation of design §5.2's basis-gauge
   requirement: covariance is nontrivial in the cell-phase gauge, while the
   re-indexed total-displacement tensor itself is invariant.
6. **Self term** `−δ_nm δ_ab·4α³/(3√π)` on same-sublattice blocks; **absent** off-sublattice (`M_T`).
7. **Boundary-label guard:** every label except `conducting_k0_omitted` errors; assert no
   surface/demag term is present in the primitive.
8. **`counts` validation** against explicit enumeration; **resource caps enforced before allocation**.
9. **α- and cutoff-independence:** at every `q in Q_A`, the §2 α bracket,
   separate-axis cutoff ladders, and joint refinement must pass `M_T` (raw)
   and `M_J` (couplings). The coarsest rung must bracket convergence, and the
   default must retain the calibrated ≥3-orders margin relative to the
   controlling tolerance.

## 4. Frozen Gate-B independent numerical reference (retained deliverable)

The oracle is a **separately-implemented scalar-Coulomb Ewald potential** (phase on `R` only, a
genuinely distinct code path from the tensor primitive), from which the dipolar tensor is obtained by
`dip_ab,nm(q) = −e^{−iq·d_nm}·∂_a∂_b φ_lat|_{x=d_nm}`. Use the frozen Cartesian step ladder
`h ∈ {4e-3,2e-3,1e-3} Å`, form O(h²)-removed Richardson estimates from both adjacent pairs, and
retain every raw and extrapolated tensor as machine-readable test output. For
`n=m`, define the differentiated scalar potential as the regularized non-self
potential: omit the `R=0` Coulomb image at every displaced finite-difference
evaluation before taking derivatives. For `n≠m`, omit no image.

Frozen acceptance:

- the two Richardson estimates differ componentwise by no more than `M_FD`;
- off-diagonal (`n≠m`) blocks from the finer Richardson estimate agree with the primitive at
  `M_FD` (calibration worst `9.2e-9 Å⁻³`);
- on `n=m`, the residual between the primitive and the scalar-oracle non-self contribution equals
  the analytic self term at `M_FD` (calibration worst `4.1e-9 Å⁻³`); and
- the primitive's α/cutoff self-convergence still independently passes the tighter `M_T`. Passing
  `M_FD` is never reported as passing `M_T`.

A brute-force `MF_dipole` spherical extrapolation at generic q bounded away from Γ is a
**secondary** cross-check only (state the summation shape and extrapolation model); it cannot
validate the Γ boundary term.

Run the retained scalar oracle at every `q_int`, for every basis pair and all
nine Cartesian tensor components, at every frozen positive-\(\alpha\) value
using that α value's matched generous cutoffs from §2. Store the primitive
tensor, each raw finite-difference tensor, both Richardson tensors, the
analytic self term, the componentwise residuals, and the pass/fail masks.
A pre-freeze double-precision arithmetic check over this complete frozen
sample gave worst primitive-reference and adjacent-Richardson absolute errors
`7.51e-10 Å⁻³` and `6.31e-10 Å⁻³` (`0.034` and `0.029` of `M_FD`);
these are calibration only and must be reproduced by the retained oracle.

## 5. Frozen Gate-C — Γ / Lorentz / demag reconciliation (formal rerun)

The §4.2 algebraic decision (**Ewald adds `0` at Γ**; the isotropic Lorentz is already inside
`dip_reg`) is fixed. Gate C reruns the reconciliation as retained tests with frozen samples/tolerances
and tests `Jpath_base_cc` and `Jgamma_cc` **separately**:

```
Jgamma_cc  (exact-Γ production, backend-agnostic physical Γ matrix):
  bruteforce: −gfac·dip_sphere_cc(0) + lorz·ones(4) + Jex_cc(0)
  ewald:      −gfac·dip_reg_cc(0)              + Jex_cc(0)          (adds 0)
Jpath_base_cc (q-path reconstruction base, normalized so ONE formula serves both):
  bruteforce: −gfac·dip_sphere_cc(0) + Jex_cc(0)
  ewald:      −gfac·dip_reg_cc(0)    + Jex_cc(0) − lorz·ones(4)
  reconstruct: J(q→0,q̂) = Jpath_base_cc + gfac(4π/Vc)(1/3 − q̂_c²)·ones(4)   [both backends]
```

Frozen checks (retained; rerun post-freeze):
1. `k=0` boundary term derived in the exact `dip=−ΣT` total-displacement convention, confirmed against
   the Gate-B oracle at `M_FD`; primitive self-convergence remains separately gated at `M_T`.
2. `dip_reg(0) = dip^(r)(0)+Σ_{G≠0}dip^(G)(0)+dip^(self)` (direction-independent), distinct from the
   directional `q→0` limit.
3. **Isolated projector:** for Γ approaches along `a*`, `b*`, `c*`, and the two r.l.u. rays
   `[1 1 1]`, `[2 1 −1]`, at `s ∈ {±1e-3,±1e-4}`, convert each ray to Cartesian `q` and verify the
   isolated reciprocal `G=0` contribution equals
   `P0_ab(q)=(4π/Vc)q̂_aq̂_b exp[−|q|²/(4α²)]` on every sublattice pair at `M_id`.
4. **Full-tensor analytic remainder:** define `R(q)=dip(q)−P0(q)·ones(ntau)` for nonzero q and
   `R(0)=dip_reg(0)`. For each positive ray/magnitude form
   `R_even(q)=[R(q)+R(−q)]/2` and `R_odd(q)=[R(q)−R(−q)]/2`.
   Raw non-Bravais off-diagonal blocks are allowed an O(q), odd-in-q imaginary term.
   Define `E_even(s)=max_rays,components|R_even(s)−R(0)|`,
   `E_odd(s)=max_rays,components|R_odd(s)|`, and
   `A_T(s)=1e-8·max_rays T_scale(0;R_even(s),R(0))`.
   At `s=1e-4`, require
   `E_even(s)≤1e-6·max_rays T_scale(0;R_even(s),R(0))`, and require contraction
   `E_even(1e-4)≤0.02·E_even(1e-3)+A_T(1e-4)`.
   Report `R_odd`; require only the O(q)-consistent bound
   `E_odd(1e-4)≤0.2·E_odd(1e-3)+A_T(1e-4)`, not that it vanish.
   At integration level verify the complete Cartesian
   `δ_ab/3−q̂_aq̂_b` q-path reconstruction for all components, not only cc.
5. Let `v=ones(4,1)/2` and let `Vperp` be any machine-orthonormal basis for its orthogonal
   complement. On every frozen Γ ray, the macroscopic term couples **only** to the uniform mode:
   `v'Δv = 4·(4π/Vc)q̂_c²` and the non-uniform block `V⊥'ΔV⊥` carries no
   `4π/Vc`-scale term (maximum absolute component `≤1e-4·(4π/Vc)`).
6. Both `Jcc0` and `Jaa0` are gated in exactly three caller-level demagnetization cases:
   **off** (`demag=0`), **sphere** (`demag=1`, top-level `opts.alpha=1`), and a
   **c-axis needle** (`demag=1`, top-level `opts.alpha=0`).
   This top-level aspect-ratio control is distinct from `opts.ewald.alpha`.
   `Jcc0`/zero-field criticality must be demag-invariant; `Jshape_cc` and demag-aware `Jaa0`
   must retain the existing observable semantics; the primitive contains no surface/demag term.
7. **Integration, cache, and provenance regression:**
   - before adoption, an absent `opts.dipole` and an explicit
     `opts.dipole='bruteforce'` must each reproduce all three legacy
     `invz_jq_modes` numerical outputs and every pre-existing `info` field
     bit-for-bit; additive provenance fields are compared separately;
   - brute-force and Ewald cache keys must differ. Cold and warm calls for
     each backend must be `isequaln` in all outputs and metadata. A cached
     payload is accepted only after exact validation of q array, lattice,
     basis, exchange, demag/aspect ratio, backend, Ewald controls, grid/offset,
     Γ policy, and schema version;
   - `opts.ewald` without `opts.dipole='ewald'`, Ewald without all four
     controls `{alpha,r_cut,g_cut,boundary}`, brute force with an Ewald
     control, an unknown backend/control/boundary, and an active ODD request
     (`opts.odd` neither absent nor false) combined with Ewald must each raise a stable
     namespaced error;
   - `invz_bz_couplings`, `invz_jq_path`, `invz_spectra_map`,
     `invz_spectra_qpath`, and `invz_run_spectra` must forward one exact
     backend/grid provenance set. A mixed BZ/path backend or precomputed
     `Jnu`/`info` that conflicts with an explicit request must error. Retained
     cold/warm field-map and q-path tests exercise both Γ policies.

Gate-C PASS requires all seven checks plus removal of the brute-force `dpRng`
trust-radius from the Ewald q-path branch: every nonzero q is evaluated
directly, while exact Γ uses the local path direction and analytic limit.

## 6. Frozen Gate-D — BZ numerical gates (half-open, both Γ policies)

Rerun the **Phase-1** structural/multiset gates on the Ewald couplings, half-open grid, **production
offset `[0 0 0]`** (the other seven `{0,½}³` offsets are translation-sensitivity diagnostics), both
**P-complete** and **P-drop** reported symmetrically. The Ewald α/cutoff ladder (§2) replaces the
`dpRng` sub-ladder; it does **not** replace the N or offset ladders.

For an `N^3` grid, an offset component `o_i in {0,1/2}` means
`q_i=-1/2+(j_i+o_i)/N`, `j_i=0,...,N-1`, followed by the canonical half-open
wrap. Thus `{0,1/2}^3` denotes half-grid-step phase offsets, not a translation
by one half of the Brillouin zone. P-complete assigns uniform weight over all
`N^3` rows. P-drop removes every exact Γ-equivalent row and assigns uniform
weight over the remaining rows.

The exact required configuration set is:

1. grid ladder `N={12,16,20}`, production offset `[0 0 0]`, frozen default
   Ewald parameters;
2. all eight offsets at every `N={12,16,20}`, frozen default Ewald parameters;
3. at `N=16`, offset `[0 0 0]`, the full real-cutoff ladder at fixed
   `C_g=13`, the full reciprocal-cutoff ladder at fixed `C_r=6.5`, the full
   α bracket with its α-matched generous cutoffs, and the frozen default versus
   joint-refinement pair from §2.

Run this set for both Γ policies. No unlisted cross-product combination is
required, and no further economization of this declared set is allowed after
results are inspected.

- **Item 1/3/4** (cheap regression assertions): points `q_i,q_j` are
  duplicates iff
  `max(abs(mod(q_i-q_j+0.5,1)-0.5))<tol_uniq`, with
  `tol_uniq=1e-12`. Every offset has `N^3` distinct points. P-complete has
  `N^3` rows and exactly one Γ row only at offset `[0 0 0]`; P-drop has
  `N^3-1` rows there and `N^3` elsewhere, with zero Γ rows after filtering.
  Require `|sum(w)-1|≤1e-12`. The P-complete `invz_sigma_crit`
  singular-pole guard on the unshifted grid must remove **exactly one**
  uniform Γ eigenvalue (zero under P-drop); any other count is a hard
  failure.
- **Item 2** (reciprocal periodicity): for every frozen `q_int` and `K`,
  require the sorted-branch spectrum at `q+K` to equal that at `q`, with
  `AbsTol_J=1e-10 meV`, `RelTol_J=1e-8`.
- **Item 5a** (multiset grid refinement): N ladder `N ∈ {12,16,20}`;
  normalized mean/var/min/max + quantiles {0.05,0.25,0.5,0.75,0.95} ÷`J_ref` (var ÷`J_ref²`), energy
  scalars `J0eff`/`Jcc0`/`max(Jnu)`; **shape gate** `|Δs| ≤ 1e-6 + 1e-3·max(|s₁|,|s₂|)`, **energy
  gate** `|ΔJ| ≤ 1e-6·J_ref + 1e-4·max(|J₁|,|J₂|)`; gate on the finest step
  `N=16→20`, with its absolute spread no larger than the `N=12→16` spread.
- **Item 5b** (full-BZ Ewald-parameter invariance): compute the same summaries and
  scalars for every configuration in declared set 3. On each three-rung
  one-axis cutoff ladder, apply the item-5 shape/energy gates to the finest
  step and require its spread not to exceed the preceding step. Every
  α-bracket member must agree pairwise with every other member at those same
  gates. The default and joint-refinement configurations must agree at those
  gates. Gate-A's tighter `M_T`/`M_J` primitive checks remain independently
  required; these multiset gates do not replace them.
- **Item 6** (offset sensitivity): for each item-5 statistic and scalar,
  define its offset spread as `max_over_offsets(value)-min_over_offsets(value)`
  in the same normalized-shape or energy units used by its gate. Report every
  spread at every N. At **`N=20`** with frozen default cutoffs, the two
  extreme-offset values must satisfy the corresponding item-5 gate, and the
  spread must not exceed its `N=16` value. **This is the gate Phase 1 failed
  at 143×; Gate-D PASS requires it to now hold.**

Gate-D produces a separate PASS/FAIL verdict for each Γ policy; no policy has
yet been selected. Eligibility requires items 2, 5a, 5b, and 6 to pass and
items 1/3/4 to be asserted for that policy. Failure means the blocker is
unresolved for that policy—do **not** proceed to physical adoption merely
because Gate-A α-convergence passed.

## 7. Frozen Gate-E — physics anchors (legacy + Jensen/HMF target)

Parameterized wrappers injecting the Ewald candidate (`opts.dipole='ewald'`, half-open, per Γ policy);
legacy tests unchanged; shared helper `invz_anchor_couplings.m` legacy-identical when its optional
fields are absent. Each anchor at **its own documented tolerance** (do not substitute a blanket bar):

Every candidate wrapper starts from `ion=invz_ion()` (`demag=0`), uses
offset `[0 0 0]`, uniform weights, the frozen §2 Ewald defaults, and
`cache=true`; the only two variants are P-complete and P-drop. Grid sizes and
physics inputs are exactly those listed below. The legacy-only `dpRng=30`
argument may remain in mirrored call signatures for parity but is ignored by
the Ewald backend and is not part of Ewald convergence.

- LiHoF₄ `Σ_c` / `Jcc0` — mirror only `test_lihof4_sigma_crit` from
  `invz_projected/tests/test_invz_sigma_crit.m`: use `N={12,24}`,
  `S_c=2*S_24-S_12`; require `|S_c-0.3004|≤0.006` and
  `|Jcc0-6.421e-3|≤0.03*|6.421e-3|`.
- Zero-field `Tc` — mirror only `test_zero_field_tc` from
  `invz_projected/tests/test_invz_critical.m`: use the same `N={12,24}`
  construction and `S_c=2*S_24-S_12`; require `|Tc-1.74 K|≤0.08 K` and
  `TcMF>Tc`.
- 310 mK critical field — mirror only `test_critical_field_at_310mK` from
  `invz_projected/tests/test_invz_critical.m`: use `N=16`, require
  `4.0 T<Bc<4.6 T`, then solve at `Bc` with `hyp=true` and require
  `|Sigma0-0.0932|≤0.02`.
- Legacy qualitative **bare-MF** ordered state + soft mode — mirror
  `test_ordered_solve_and_soft_mode` from
  `invz_projected/tests/test_invz_ordered_phase.m`: use `N=16`,
  `T=0.31 K`, `hyp=true`, `B={2,4,6} T`, and
  `w=(0.005:0.005:0.7)' meV`. Require ordered and converged states at
  `2,4 T`; positive `crit` at both with `crit(2)>crit(4)`; `m0(2)>m0(4)`;
  `E2>E4>0.05 meV`; `E2<0.6 meV`; and no ordered state at `6 T`.
  (The bare-MF anchor requires a finite ordered-side gap and CANNOT stand in for the user-facing
  symptom.)
- **NEW — target-path Jensen/HMF ordered acceptance** (`ordered_1z='jensen'`, same route as
  `invz_spectra_map`). The complete frozen input is:

  ```text
  ion              = invz_ion(), demag=0
  T                = 0.1 K
  fields           = [1.0, 2.5, 3.5, 6.0] T
  field_dir        = [1 0 0]
  w_GHz            = (0:0.01:5.5)' GHz
  C                = invz_const()
  w_internal       = w_GHz * C.Gh2mV
  eta              = 5e-5 meV
  hyp              = true
  peak_wmin        = 0
  bz_tol           = 1e-9 T
  ordered_1z       = 'jensen'
  solve_opts       = struct('transverse_mf','legacy_x')
  parallel/verbose = false/false
  BZ grid          = [16 16 16], half-open, offset [0 0 0]
  dipole backend   = Ewald with the frozen §2 defaults
  Gamma policy     = the wrapper's P-complete or P-drop variant
  ```

  **PASS criteria, applied independently to both Γ-policy variants:**

  1. At `B={1.0,2.5,3.5} T`, `phase_1z=1`; `Sigma0`, `m_1z`, and `D_ord` are finite;
     `m_1z>0`, `D_ord>0`; and each `S.chiz` column contains finite data with a strictly positive
     spectral maximum.
  2. The ordered static quantities soften toward the boundary:
     `m_1z(1.0)>m_1z(2.5)>m_1z(3.5)>0` and
     `D_ord(1.0)>D_ord(2.5)>D_ord(3.5)>0`. Report `Epeak` at these points when the censored peak
     picker returns a finite value, but do not use it in place of the static-mass gate.
  3. At `B=6.0 T`, `phase_1z=2`, `Sigma0` and `crit_pm` are finite, `crit_pm>0`, and the `S.chiz`
     column contains finite data with a strictly positive spectral maximum. By branch contract,
     `m_1z` and `D_ord` are NaN here.
  4. No declared point has `phase_1z=0`, an all-NaN `S.chiz` column, or a non-finite state variable
     required by its branch.

  No monotonic relation is imposed between the ordered `3.5 T` spectrum and the far-PM `6.0 T`
  spectrum: the mode may harden again beyond the transition and `D_ord` is undefined in the PM
  branch. These four sparse fields test branch coverage and the ordered static trend; they do not
  claim a solver-grade continuity measurement at the boundary. Failure leaves the user-facing
  spectra defect unresolved even if the legacy bare anchor passes.

## 8. Frozen Gate-F — Γ-policy decision tree (no scoring-by-closeness)

1. A policy (P-complete / P-drop) is eligible only if it independently passes Gates D **and** E.
2. If neither is eligible → stop, no production selection.
3. If exactly one is eligible → select it.
4. If both are eligible:
   - **Common-support implementation check:** for every declared Gate-D
     configuration, remove the exact-Γ row, when present, from P-complete and
     apply the P-drop normalization. Corresponding sorted branches must agree
     at `M_J`, as must corresponding `Juni`, `Jcc0`, and `Jaa0`; their item-5
     summaries and `J0eff`/`Jcc0` must agree at the item-5 shape/energy gates.
     This checks implementation identity on the same mathematical support.
   - **Direct integrated-policy check:** at the production `N=16`, offset
     `[0 0 0]`, frozen-default configuration, directly compare P-complete and
     P-drop normalized mean, variance, and the five frozen quantiles at the
     item-5 shape gates. Do **not** compare raw min, max, `max(Jnu)`, or another
     extremum whose support deliberately differs at Γ.
   - **Physical-output check:** require
     `|Delta Sc|≤0.006`,
     `|Delta Jcc0|≤0.03*6.421e-3 meV`,
     `|Delta Tc|≤0.08 K`,
     `|Delta Bc|≤0.10 T`, and
     `|Delta Sigma0(Bc)|≤0.02`.
     For the legacy bare-MF and Jensen/HMF anchors, whose frozen gates are
     categorical/ordinal rather than proximity tolerances, both policies must
     independently satisfy every predicate in §7; report their numeric
     differences, but do not invent an after-the-fact distance tolerance.
   If any of these three checks fails, stop for a Γ-theory/finite-size audit.
   If all pass, select **P-drop** (predeclared measure-zero critical-pole
   regularization and production-compatible tie-breaker).
5. An unexpected singular-entry count or a Gate-E failure makes that policy
   ineligible. The surviving policy may be selected only through steps 2–4;
   there is no majority vote and no post-hoc tolerance change.

## 9. Frozen adoption + `Bc_PM` re-freeze (Step-7 procedure)

Only after Gates A–C pass globally and Gate F selects a policy that passed
Gates D and E: freeze the production `[16 16 16]` half-open
grid, uniform weights, offset `[0 0 0]`, selected Γ policy,
`conducting_k0_omitted` boundary, `α₀`, and cutoffs; flip the supported non-ODD
projected-path default to `ewald`; repopulate caches under the new schema; then
resume Phase 3.

Before any resumed ordered audit, recompute `Bc_PM` by the same independent
PM-mass route used by Task 2a:

```text
ion = invz_ion()
[Jf,info,Jaa0] = invz_bz_couplings(ion, selected frozen production options)
Bc_PM = invz_critical(ion, 0.1, Jf, struct('J0eff',info.Jcc0,'Jxx0',Jaa0))
F_new = [0.25,0.55,0.80] * Bc_PM
```

Re-freeze all and only the enumerated derived artifacts:

- the `Bc_PM` provenance comments and three derived field literals in
  `invz_projected/invz_task2_matrix_enumerate.m` (the independent `2.85 T`
  defect anchor remains unchanged);
- the `Bc_PM` assertion, the three derived-field assertions, and any matching
  enumerator field table in
  `invz_projected/tests/test_invz_task2_matrix.m`;
- a superseding Task-2 freeze record and `docs/invzp_task2_report.md`; and
- every coupling/Task-2 cache whose provenance identifies the old backend,
  grid convention, Γ policy, α/cutoffs, or old derived field table.

Record the new root and full-precision derived values before rounding the
three six-decimal field literals, and retain the existing `1e-3 T` regression
tolerance unless a separate prospective re-registration changes it. The
brute-force default and `MF_dipole` remain reachable indefinitely.

## 10. Selection / escalation / stop (frozen)

- **Select** only per the Gate-F tree, after Gates A–E.
- **Escalation:** if the Ewald primitive itself cannot meet the §2 accuracy guards within the frozen
  resource caps, that is a hard failure of this construction (not a silent cap relaxation) — report and
  stop.
- **Stop / hard-fail:** a failure of global Gate A, B, or C, or failure to
  select any eligible policy through Gate F, keeps the brute-force default,
  reports the first failed gate and diagnostic, and remains short of any
  3A/3B/4 path. A Gate-D/E failure disqualifies only that Γ policy as stated
  in Gate F. If the Jensen target fails for every otherwise viable policy,
  the user-facing symptom is not resolved by this route.

## 11. Execution order (frozen) + what this does NOT authorize

1. Freeze this pre-registration (this file). 2. Implement the opt-in primitive + retained Gate-A/B
tests + primitive-level Gate C (no integration, no default change). 3. Integrate the opt-in
`invz_jq_modes`/BZ/q-path/spectra backend + cache schema + provenance; run remaining Gate C +
legacy-regression (bruteforce byte-identical). 4. Gate D (both Γ). 5. Gate E (legacy + Jensen). 6.
Gate F. 7. On success only: flip default, re-freeze `Bc_PM` (§9), resume Phase 3.

This document does **not** authorize any default change, ODD/tensor-path adoption, or Phase-3/3A/3B/4
work until the gates pass in that order. The frozen Task-2 §D observable tolerance is unchanged.

## 12. Errata (dated prospective re-registration; per the freeze rule, not in-place edits)

- **E1 (2026-07-24) — §5 Gate-C check 5 target.** Check 5's written target `v'·Δ_cc·v =
  4·(4π/Vc)q̂_c²` is the **q→0 limit**. The exact target, consistent with check 3's isolated projector
  `P0_ab(q)=(4π/Vc)q̂_aq̂_b·exp(−|q|²/(4α²))`, retains the Gaussian factor:
  `v'·Δ_cc·v = 4·(4π/Vc)q̂_c²·exp(−|q|²/(4α²))`. This is because `Δ_cc` (the isolated `G=0` summand) is
  uniform over sublattice pairs, so `v'·Δ_cc·v = 4·P0_cc(q)` by linear algebra. The Step-4 Gate-C C5
  test (`invz_projected/tests/test_invz_dipole_ewald_gammaC.m`) implements the exact form and passes at
  `M_id` margin 0.0; the literal no-`exp` reading is unpassable at `M_id` (deviates ~2.7e-7 at the frozen
  ray magnitudes). Read check 5's target with the `exp(−|q|²/(4α²))` factor; this is a prose
  clarification — no tolerance, sample, or scope change. Confirmed by the Step-4 whole-branch review.
  **Relevant to Step 5:** when the caller-level Gate-C q-path reconstruction is wired, do not
  re-introduce the small-q limit as an exact target.
