# Ewald Dipolar Primitive (Step 4) — Corrected Implementation Plan

**Status:** review-corrected and ready to execute under the frozen
`docs/invzp_ewald_prereg.md` dated 2026-07-24.

**Goal:** add the opt-in Ewald dipolar-tensor primitive
`invz_dipole_ewald.m` and retained tests for Gate A, Gate B, and the
primitive-computable parts of Gate C. Production callers, defaults, caches,
ODD/tensor paths, and physical anchors remain untouched.

**Authority order:** if this plan conflicts with another document, follow:

1. `docs/invzp_ewald_prereg.md` (frozen acceptance contract);
2. `docs/invzp_ewald_design.md` and `docs/invzp_ewald_derivation.md`
   (derivation and architecture);
3. this plan.

Do not loosen or reinterpret a frozen tolerance after seeing a result. A
failure is evidence to diagnose, not permission to alter the preregistration.

## Scope and deliverables

Add exactly these seven files:

- `invz_dipole_ewald.m`
- `invz_projected/tests/invz_ewald_fixtures.m`
- `invz_projected/tests/invz_ewald_metrics.m`
- `invz_projected/tests/invz_scalar_ewald_ref.m`
- `invz_projected/tests/test_invz_dipole_ewald.m`
- `invz_projected/tests/test_invz_dipole_ewald_ref.m`
- `invz_projected/tests/test_invz_dipole_ewald_gammaC.m`

Do not modify:

- `MF_dipole.m`;
- any existing `invz_projected/*.m` production file;
- `invz_common/*`;
- the frozen design, derivation, or preregistration;
- existing tests or physical-anchor expectations.

This step does not authorize production integration, a default change,
Gate D/E/F, Bc re-freezing, ODD/tensor adoption, or Phase 3/3A/3B/4.

## Fixed primitive contract

The public signature is:

```matlab
[dip, counts, geom] = invz_dipole_ewald(q, a, tau, eopts, geom)
```

The fifth argument is optional and is a reusable, q-independent geometry
payload. For one q row, `dip` is `[3,3,ntau,ntau]`; for `nq>1`, it is
`[3,3,ntau,ntau,nq]`. `counts.recip_used` is always `[nq,1]`.

Use the following indexing consistently:

```text
output block       dip(:,:,n,m)
pair displacement  d_nm = tau(m,:)*a - tau(n,:)*a
pair identity      dip(:,:,m,n) = conj(dip(:,:,n,m))
```

The frozen expression is:

```text
dip = dip_real + dip_recip + dip_self

dip_real_nm =
  -sum_{R+d_nm, nonself, |R+d_nm|<=r_cut}
     g_ab(R+d_nm) exp[-i q_cart.(R+d_nm)]

dip_recip_nm =
  +(4*pi/Vc) sum_{k=q_cart+G, k~=0, |k|<=g_cut}
     (k_a*k_b/|k|^2) exp[-|k|^2/(4*alpha^2)] exp[+i G.d_nm]

dip_self_nm =
  -delta_nm delta_ab 4*alpha^3/(3*sqrt(pi))
```

The primitive contains no Lorentz, surface, shape, demagnetization, exchange,
or P-complete/P-drop term. Its only boundary label is
`conducting_k0_omitted`.

### Important implementation invariants

- Compute `B = 2*pi*inv(a)'`; do not use a rounded reciprocal basis.
- The deterministic default, where a test needs it, is
  `alpha0 = sqrt(pi)/abs(det(a))^(1/3)`. Do not substitute the displayed
  rounded value `0.268431`.
- Reduce each q row by `K=floor(q+0.5)` into `[-0.5,0.5)`, evaluate there,
  and restore
  `exp[-i*2*pi*K.(tau(m,:)-tau(n,:))]`.
- The cached reciprocal candidate set must cover the exact union over the
  whole canonical cube. It must contain `G=[0,0,0]`: only exact
  `k=qbar+G=0` is omitted per q. Never discard `G=0` globally.
- Build every `(n,m)` pair independently in the implementation. Do not fill
  half the tensor by conjugation; otherwise the retained Hermiticity test
  becomes tautological and cannot catch a displacement/phase indexing error.
- Geometry reuse must be exact-fingerprint checked. The fingerprint contains
  `a`, `tau`, `alpha`, `r_cut`, `g_cut`, boundary, q reduction convention,
  and schema.
- Input validation is strict: `q` is finite real `[nq,3]`, `a` is finite real
  nonsingular `[3,3]`, `tau` is finite real `[ntau,3]`, and `eopts` has
  exactly the four finite controls `alpha`, `r_cut`, `g_cut`, and `boundary`.
  Unknown or missing fields raise stable `invz:*` errors.

## Preflight is a pre-allocation gate

This ordering is mandatory:

1. Validate small inputs and controls.
2. From `ntau` and `nq`, first check a dimension-only lower bound for the
   unavoidable pair-indexed metadata and output. This must happen before
   allocating an `ntau^2` displacement/count array.
3. Compute only `Vc`, `B`, singular values, one pair displacement at a time,
   and scalar conservative bounds.
4. Build the conservative manifest and apply all three hard caps.
5. Only after all caps pass may `ndgrid`, candidate vectors, screened tensors,
   reciprocal phases, q work arrays, or output arrays be allocated.

Use the frozen bounds:

```matlab
nmax_r = ceil((r_cut + norm(d_nm))/min(svd(a)));
real_cube_bound(n,m) = (2*nmax_r + 1)^3;

qcorner = all eight rows of {-0.5,+0.5}^3;
qmax = max(vecnorm(qcorner*B,2,2));
nmax_G = ceil((g_cut + qmax)/min(svd(B)));
recip_cube_bound = (2*nmax_G + 1)^3;
```

Before enumeration require:

```text
every real_cube_bound(n,m) <= 3.0e6
recip_cube_bound           <= 3.0e6
estimated_peak_bytes       <= 4*2^30
```

The byte manifest uses conservative cube bounds wherever an actual retained
count is not yet known. It includes all retained and temporary numeric and
logical arrays: integer meshes, Cartesian meshes, radii/masks, retained
`x/g_ab`, reciprocal `Ghkl/Gcart`, pair phases, q-specific work arrays,
and complex output. Apply the frozen `1.25` allocator margin. Retain this
initial conservative estimate in `counts`; do not replace it with a smaller
post-allocation estimate.

When a geometry is reused with a different `nq`, rerun the q-work/output part
of the preflight before allocating the new output. A manifest captured from a
previous one-q call is not sufficient for a later multi-q call.

Construct the reciprocal candidate union only after preflight. Start from the
preflight integer cube and retain every integer `G_hkl` for which

```text
min_{qbar in [-0.5,0.5)^3} norm((qbar + G_hkl)*B) <= g_cut.
```

Implement the three-dimensional box-constrained quadratic minimization
deterministically (for example, enumerate the 27 free/lower/upper active
sets); it must not require Optimization Toolbox. Treat the upper face with
the frozen half-open convention and include numerically ambiguous boundary
candidates conservatively. Per-q filtering still applies the requested
`|qbar*B+Gcart|<=g_cut` and omits only the exact zero k vector.

## Common test helpers

`invz_ewald_fixtures()` returns exactly:

- the three frozen `q_int`;
- `delta_q=2^-40`;
- `Q_A = Gamma + q_int + six faces + twelve edges + eight corners`;
- all eight frozen reciprocal translations `K`;
- origin shift `[0.137,-0.211,0.089]`;
- Bravais shifts `[1,0,0]` and `[0,-1,1]`;
- the 200 deterministic Hessian vectors generated in the exact frozen RNG
  order;
- the explicit two-site small-shell fixture;
- `alpha0` computed from the live `invz_ion()` cell.

The fixture constructor asserts its own cardinalities and uniqueness so a
malformed sample cannot silently reduce test coverage.

`invz_ewald_metrics()` returns symmetric helpers for `M_T`, `M_J`, `M_id`,
`M_HFD`, and `M_FD`. Each result contains at least:

```text
pass
worst_margin
worst_ratio
max_abs_error
```

The helpers reject NaN/Inf. `M_T` and `M_id` derive one scale from the complete
tensor being compared, not a different scale per sublattice block.
`M_FD` remains dimensioned in inverse cubic Angstrom and is never labelled
`M_T`.

For primitive-level coupling convergence, add a local test-only mapper in
`test_invz_dipole_ewald.m` that exactly mirrors the non-ODD Ewald algebra
without calling or modifying `invz_jq_modes`:

```matlab
ex = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX);
Jcc = -C.gfac*squeeze(dip(3,3,:,:)) ...
       + sign(ion.J12)*squeeze(ex(3,3,:,:));
Jaa = -C.gfac*squeeze(dip(1,1,:,:)) ...
       + sign(ion.J12)*squeeze(ex(1,1,:,:));
```

Hermitize these matrices exactly as the future caller will. Return sorted cc
eigenvalues (`Jnu`), the uniform projection (`Juni`), and at exact Gamma the
uniform `Jcc0` and `Jaa0`. Add no Lorentz term to the Ewald branch.

## Task 0 — authorization, baseline, and workspace fence

- [ ] Confirm the preregistration header still says
      `Status: FROZEN 2026-07-24`.
- [ ] Confirm none of the seven deliverables already contains user work. If
      one exists, stop and inspect rather than overwrite it.
- [ ] Record `git status --short`; preserve all unrelated changes and
      untracked files.
- [ ] Inventory pre-existing cache/MAT artifacts before the baseline because
      some legacy tests may legitimately exercise their own caches.
- [ ] Run the existing suite before adding files:

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch \
  "addpath('.'); addpath('invz_common'); r=runtests('invz_projected/tests'); \
   fprintf('passed=%d failed=%d incomplete=%d\n',sum([r.Passed]),sum([r.Failed]),sum([r.Incomplete]))"
```

- [ ] The previously observed baseline was 225 passed / 0 failed /
      23 incomplete. Treat the live run as authoritative. Any new failure is
      investigated before implementation; do not merely rewrite this number.

No commit is made in Task 0.

## Task 1 — strict input contract, geometry, exact candidate union, preflight

**Files:** create `invz_dipole_ewald.m` and
`invz_projected/tests/test_invz_dipole_ewald.m`.

### Tests first

Add tests for:

- [ ] exact output/geometry/counts schema;
- [ ] exact `alpha0` rule rather than the rounded display value;
- [ ] real vectors within cutoff and exact self-cell omission;
- [ ] deterministic ordering of retained real and reciprocal vectors;
- [ ] `G_hkl=[0,0,0]` is present in the cached candidate union;
- [ ] all actual retained counts are no greater than conservative bounds;
- [ ] every manifest row has
      `name,class,is_complex,size,numel,bytes`;
- [ ] the retained preflight estimate equals `1.25*sum(manifest.bytes)` and
      is based on conservative bounds;
- [ ] one-q and multi-q calls update the q-work/output estimate correctly;
- [ ] identical geometry reuse is identical, while changing each fingerprint
      field separately raises `invz:ewaldGeomReuse`;
- [ ] missing/unknown controls, invalid shapes, nonfinite values, nonpositive
      controls, singular `a`, and an unsupported boundary raise stable errors;
- [ ] both hard accuracy guards reject values just below 4.5 and accept the
      exact boundary.

The first run must fail because `invz_dipole_ewald` is absent.

### Implement

- [ ] Implement the public validator and q-independent geometry builder.
- [ ] Compute and enforce the conservative preflight before any large
      allocation.
- [ ] Implement the exact candidate-union filter; retain `G=0`.
- [ ] Enumerate full real geometry for every ordered basis pair.
- [ ] Precompute `g_ab(x)` and `exp(+i*G.d_nm)` for every ordered pair.
- [ ] Store the frozen fingerprint and counts schema.
- [ ] It is acceptable for `dip` to be an internal `[]` placeholder while
      running this task's geometry-only red test, but do not commit or hand
      off that incomplete public state. Proceed directly through Task 2
      before the first implementation checkpoint.

Run the new geometry tests, then proceed directly to Task 2. Do not commit an
incomplete primitive.

## Task 2 — real, reciprocal, self assembly and gauge restoration

**Files:** modify `invz_dipole_ewald.m` and
`test_invz_dipole_ewald.m`.

### Tests first

- [ ] One q returns `[3,3,ntau,ntau]`; multiple rows return the documented
      five-dimensional result and `[nq,1]` reciprocal counts.
- [ ] A test-only reconstruction independently forms real, full reciprocal,
      isolated `G=0`, and self contributions from `geom`. Their sum agrees
      with the returned complete tensor at `M_id`.
- [ ] Every ordered basis pair and all nine components are reconstructed;
      no test helper fills conjugate blocks by assignment.
- [ ] Returned pair conjugation/Hermiticity passes `M_id`.
- [ ] `counts.recip_used(i)` equals explicit q-specific enumeration.
- [ ] At nonzero reduced q, the `G=0` term is used. At exact Gamma it is the
      one exact k vector omitted.
- [ ] A benign generic-q `MF_dipole` spherical sequence is recorded as the
      frozen secondary cross-check only. Do not promote an unregistered
      brute-force truncation threshold to Gate-A acceptance.

### Implement

For each q row:

1. reduce to `qbar` and retain `K`;
2. use `qcart=qbar*B`;
3. filter candidates by `|qcart+Gcart|<=g_cut`;
4. omit only rows with all Cartesian k components exactly zero;
5. assemble real, reciprocal, and self terms for all ordered `(n,m)`;
6. apply the extended-zone gauge phase;
7. preserve the documented output rank.

Do not allocate both an unused four-dimensional `dip` and a second full
five-dimensional output. The preflight manifest must match the arrays that
are actually allocated.

Run the new test file and the full suite, then checkpoint the two named files.

## Task 3 — frozen fixtures, metrics, and Gate-A test 1

**Files:** create `invz_ewald_fixtures.m` and
`invz_ewald_metrics.m`; modify `test_invz_dipole_ewald.m`.

### Finite-difference implementation

For each Hessian step `h`, use:

```text
diagonal:
  [f(x+h e_a) - 2f(x) + f(x-h e_a)] / h^2

mixed a~=b:
  [f(x+h e_a+h e_b) - f(x+h e_a-h e_b)
   -f(x-h e_a+h e_b) + f(x-h e_a-h e_b)] / (4h^2)
```

Do not use the mixed formula with `a=b`; that silently turns the labelled
step h into a `2h` diagonal stencil.

From `h1=4e-3`, `h2=2e-3`, `h3=1e-3` form both:

```matlab
R12 = (4*H_h2 - H_h1)/3;
R23 = (4*H_h3 - H_h2)/3;
```

### Retained tests

- [ ] Generate the exact frozen samples and assert their counts.
- [ ] At each of the five positive alpha values, test all 200 vectors.
- [ ] `R12` and `R23` differ by at most `1e-7*H_scale`.
- [ ] `R23` passes `M_HFD` against the independent closed formula with zero
      sign mismatches on gated components.
- [ ] Evaluate the closed formula at `alpha=0` exactly and require the bare
      tensor, rather than approximating zero by `1e-8`.
- [ ] Bridge the independent formula to production: for retained real
      vectors from every ordered pair, compare every
      `geom.real{n,m}.gab` with the independently evaluated closed expression
      at `M_id`. This
      prevents a duplicated correct formula in the test from passing while
      the private production implementation is wrong.
- [ ] Log the retained worst analytic and adjacent-Richardson ratios.

Fix formula/stencil defects, never the frozen tolerances. Run the new file and
the full suite, then checkpoint the three named files.

## Task 4 — Gate-A tests 2, 6, 7, and 8

**Files:** modify `test_invz_dipole_ewald.m`.

### A2: explicit small-shell real-space signs and phase

- [ ] Independently enumerate integer cells for the exact frozen two-site
      fixture.
- [ ] Build the complete `[3,3,2,2]` screened real-space tensor directly
      from the defining formula, including the leading minus sign and
      total-displacement phase.
- [ ] Independently reconstruct the primitive's real contribution from its
      retained geometry.
- [ ] Compare the two complete tensors at `M_T`.

Do not compare the full Ewald tensor to a “reference” that quietly adds a
second copy of the reciprocal implementation. Gate A2 is specifically the
finite screened real shell.

### A6: self term

- [ ] Reconstruct real plus reciprocal nonself contributions from `geom`.
- [ ] Subtract them from the returned primitive.
- [ ] Build one complete expected self tensor with
      `-4*alpha^3/(3*sqrt(pi))*eye(3)` on every same-site block and zero on
      every off-site block; compare the complete actual and expected tensors
      once at `M_T`.
- [ ] Report same-site and off-site maxima separately without replacing the
      complete-tensor frozen scale.
- [ ] Repeat at every frozen positive alpha.

An imaginary-trace or Gamma-realness placeholder is not a self-term test.

### A7: boundary and absence of a surface term

- [ ] Every boundary string except `conducting_k0_omitted` errors.
- [ ] Unknown surface/demag controls error as unknown primitive controls.
- [ ] Reconstruction proves the output contains exactly real + reciprocal +
      self and no fourth macroscopic term.

### A8: counts and caps

- [ ] Compare every `counts.real_pair` entry with independent integer-cell
      enumeration.
- [ ] Compare `recip_candidates` with an independent candidate-union
      enumeration and each `recip_used` with independent per-q filtering.
- [ ] Exercise all three hard caps separately with synthetic small-input
      configurations and stable errors:
      `invz:ewaldRealCap`, `invz:ewaldRecipCap`, and
      `invz:ewaldMemoryCap`.
- [ ] These cases must return the cap error without attempting the prohibited
      large allocation. Keep timing diagnostic only; do not turn elapsed time
      into a scientific assertion.

Run the file and full suite, then checkpoint the test file.

## Task 5 — Gate-A tests 3 and 4

**Files:** modify `test_invz_dipole_ewald.m`.

### A3: reciprocal gauge and periodic spectrum

For every `q in Q_A` and every frozen `K`:

- [ ] construct the full gauge-predicted tensor and compare it to
      `dip(q+K)` once at complete-tensor `M_id`;
- [ ] use the test-only coupling mapper, including exchange but no Ewald
      Lorentz term, and compare all four sorted cc branches with
      `AbsTol=1e-10 meV`, `RelTol=1e-8`.

Do not apply a meV tolerance directly to an inverse-cubic-Angstrom dipolar
eigenvalue.

### A4: half-open reduction and candidate completeness

- [ ] Test exact `+0.5` reduction to `-0.5` as well as all eight frozen
      near-upper-face corners.
- [ ] For every frozen corner and every frozen K, compare the restored full
      tensor at `M_id`.
- [ ] Independently enumerate every q-specific G satisfying the cutoff and
      assert it is present in the cached candidate set.
- [ ] Explicitly cover `G=0`, face, edge, and corner cases.

The earlier scalar loops over only three x-axis shifts are insufficient.
Run the file and full suite, then checkpoint the test file.

## Task 6 — Gate-A test 5 (A5)

**Files:** modify `test_invz_dipole_ewald.m`.

For every `q in Q_A`:

- [ ] require full-tensor pair conjugation/Hermiticity at `M_id`;
- [ ] at exact Gamma require the complete tensor, including off-site blocks,
      to be real at `M_id`;
- [ ] apply the frozen common fractional-origin shift and require raw
      complete-tensor invariance at `M_id`;
- [ ] for every basis representative and both frozen Bravais shifts, rebuild
      geometry (thereby explicitly reindexing real cells) and require raw
      complete-tensor invariance at `M_id`;
- [ ] convert both results to cell-phase gauge and require the complete
      transformed tensor to obey
      `Dcell' = U' * Dcell * U` at `M_id`, with
      `U_nn=exp(+i*2*pi*q.L_n)`.

Build a full predicted four-dimensional tensor before calling `M_id`; do not
replace the frozen metric with an ad hoc block-relative `1e-10` tolerance.
Vectorize the q rows and reuse each compatible geometry so the exhaustive
sample remains practical.

Run the file and full suite, then checkpoint the test file.

## Task 7 — Gate-A test 9

**Files:** modify `test_invz_dipole_ewald.m`.

Evaluate every configuration below at every `q in Q_A`, reusing one geometry
per parameter set:

- [ ] alpha multipliers `{0.6,0.8,1.0,1.2,1.5}` with
      `r_cut=6.5/alpha`, `g_cut=13*alpha`;
- [ ] real ladder `C_r={4.5,5.0,5.5}` at fixed `C_g=13`;
- [ ] reciprocal ladder `C_g={9,10,11}` at fixed `C_r=6.5`;
- [ ] production default `(5.5,11)` versus joint refinement `(6.0,12)`.

Acceptance:

- [ ] all required raw complete-tensor comparisons pass `M_T`;
- [ ] alpha-bracket members agree independently of alpha; retain pairwise
      diagnostics so cancellation is visible;
- [ ] both one-axis ladders demonstrate the frozen coarsest-to-finest
      convergence bracket without combining the two cutoff axes;
- [ ] the default/joint comparison passes;
- [ ] the test-only mapper's `Jnu`, `Juni`, `Jcc0`, and `Jaa0` pass `M_J`;
- [ ] the default-versus-joint raw and coupling utilization of the
      controlling tolerances is at most `1e-3`, preserving the frozen
      three-orders margin;
- [ ] elapsed time is logged per configuration but never used as a
      numerical PASS/FAIL condition.

Testing only `q_int(1,:)`, omitting `M_J`, or omitting joint refinement does
not satisfy Gate A9.

Run the file and full suite, then checkpoint the test file.

## Task 8 — Gate-B independent scalar-Coulomb oracle

**Files:** create `invz_scalar_ewald_ref.m` and
`test_invz_dipole_ewald_ref.m`.

The oracle must not call the tensor primitive, its local functions, or its
geometry builder. It returns a machine-readable struct containing:

```text
h
raw                    tensor at each of the three h values
richardson_coarse
richardson_fine
self_analytic
adjacent_residual
```

### Correct scalar potential

For canonical q, `qcart=q*B`, and displacement x:

```text
phi_lat(x) =
  sum_{R, |R+x|<=r_cut}
      exp(-i*qcart.R) erfc(alpha*|R+x|)/|R+x|
  + (4*pi/Vc) sum_{k=qcart+G, k~=0, |k|<=g_cut}
      exp(-|k|^2/(4*alpha^2))/|k|^2 * exp(+i*k.x)
```

Then:

```text
dip_nm = -exp(-i*qcart.d_nm) * Hessian(phi_lat) at x=d_nm.
```

The reciprocal denominator, Gaussian, cutoff, and phase use `k=q+G`, not
bare `G`. For `n=m`, exclude the integer cell `R_hkl=[0,0,0]` from every
displaced finite-difference evaluation by its cell index, not by testing
whether `|R+x|` is near zero. For `n~=m`, exclude no cell.

Use the same correct diagonal/mixed finite-difference stencils as Task 3.
Return all three raw tensors and both adjacent-pair Richardson tensors.

### Retained Gate-B test

For all three `q_int`, all five positive alpha values with matched generous
cutoffs, every ordered basis pair, and all nine components:

- [ ] the two Richardson tensors agree at `M_FD`;
- [ ] every off-site finer-Richardson block agrees with the primitive at
      `M_FD`;
- [ ] on each same-site block, primitive minus regularized oracle equals the
      analytic self tensor at `M_FD`;
- [ ] all raw, extrapolated, residual, and pass-mask arrays are finite and
      exposed in the returned machine-readable result;
- [ ] the tighter primitive self-convergence remains independently covered
      by Gate A9 and is reported as `M_T`, never inferred from `M_FD`.

Compute all ordered pairs independently; do not force oracle Hermiticity by
copying conjugate blocks. Batch/reuse the oracle's own scalar candidate
meshes where possible, but keep its code path independent of the tensor
primitive.

Run the Gate-B file, then all three new Ewald test files, then the full suite.
Checkpoint only the oracle and its new test.

## Task 9 — primitive-level Gate C, checks 1–5

**Files:** create `test_invz_dipole_ewald_gammaC.m`.

Gate-C checks 6 and 7 are caller integration/cache/provenance tests and remain
Step 5. The complete Cartesian q-path reconstruction in check 4 also remains
an integration test. This task covers the boundary/oracle evidence from
Task 8 and the primitive-computable checks below.

Add a local independent parts reconstruction using public `geom` data. It
returns real, reciprocal-total, isolated reciprocal `G=0`, and self tensors
and first proves their sum equals the primitive at `M_id`.

Use the five frozen reciprocal-space rays and
`s={+1e-3,-1e-3,+1e-4,-1e-4}`, with `alpha=alpha0`,
`r_cut=6.5/alpha0`, and `g_cut=13*alpha0`.

### C1/C2: boundary convention and regular Gamma tensor

- [ ] Cite/consume the retained Gate-B result as the independent
      scalar-oracle confirmation of the sign and self convention.
- [ ] At exact Gamma prove the per-q filter omits only `G=0`.
- [ ] Reconstruct
      `dip_reg(0)=dip_real(0)+sum_{G~=0}dip_recip(0)+dip_self`
      and compare the complete tensor to the returned Gamma primitive at
      `M_id`.
- [ ] Assert there is no directional or surface term at exact Gamma.

### C3: isolated projector

For every signed ray/magnitude, every ordered sublattice pair, and every
Cartesian component:

- [ ] extract the implementation's isolated `G=0` summand;
- [ ] compare it at `M_id` to
      `(4*pi/Vc)*qhat*qhat.'*exp(-|q|^2/(4*alpha^2))`;
- [ ] prove the `G=0` candidate participates for nonzero q.

Subtracting the analytic projector and checking only `isfinite` is not this
gate.

### C4: even/odd analytic remainder

Define complete tensors:

```text
R(q)      = dip(q) - P0(q) replicated over every (n,m) pair
R(0)      = dip_reg(0)
R_even(q) = (R(q)+R(-q))/2
R_odd(q)  = (R(q)-R(-q))/2
```

Aggregate over all rays, pairs, and components exactly as preregistered.
Retain and assert:

```text
E_even(1e-4) <= 1e-6*Tscale_even(1e-4)
E_even(1e-4) <= 0.02*E_even(1e-3) + A_T(1e-4)
E_odd (1e-4) <= 0.20*E_odd (1e-3) + A_T(1e-4)
A_T(s)        = 1e-8*Tscale_even(s)
```

Report the odd part; do not require it to vanish.

### C5: uniform-mode support

For every frozen ray, use the isolated implementation `G=0` contribution,
`v=ones(4,1)/2`, and a deterministic machine-orthonormal `Vperp`:

- [ ] `v'*Delta_cc*v` agrees with
      `4*(4*pi/Vc)*qhat_c^2` at the exact-identity scale;
- [ ] every component of `Vperp'*Delta_cc*Vperp` is at most
      `1e-4*(4*pi/Vc)`.

Do not use `dip(q)-dip(0)` as though it were purely macroscopic; that
difference also contains the finite-q analytic remainder.

Run Gate C, all Ewald tests, and the full suite. Checkpoint the new test file.

## Task 10 — closure and handoff

- [ ] Run each new test file independently so failures map to one gate.
- [ ] Run the full `invz_projected/tests` suite and report live
      passed/failed/incomplete counts.
- [ ] Run `git diff --check`; before each scoped commit also run
      `git diff --cached --check` after staging only that task's files.
- [ ] Verify the implementation delta is limited to the seven scoped new
      files. Leave this reviewed plan and all other pre-existing user changes
      untouched.
- [ ] Compare the cache/MAT inventory with Task 0. No new artifact
      attributable to the Ewald primitive tests may remain; preserve
      pre-existing user/legacy caches. Verify no physical report, production
      file, or frozen document was modified.
- [ ] Record Gate A1–A9, Gate B, and primitive Gate-C C1–C5 as separate
      verdicts with their worst frozen metrics. Do not collapse them into one
      “tests pass” statement.
- [ ] If any global Gate-A/B/primitive-C item fails, keep the brute-force
      default and stop before Step 5.
- [ ] If all pass, hand off Step 5 as the next separately reviewed plan:
      opt-in caller integration, cache/provenance schema, and remaining
      Gate-C checks. The default still remains brute force.

When committing during execution, stage only the files named by the current
task. Never stage the user's unrelated uncommitted or untracked files, and do
not add a `Co-Authored-By` trailer.
