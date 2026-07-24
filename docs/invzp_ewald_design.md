# Ewald dipolar-summation reformulation — DESIGN (PREREG-REVIEWED)

> **STATUS ADDENDUM (2026-07-24 — supersedes the pre-freeze status language below).** The five §10
> amendments were incorporated and the preregistration `docs/invzp_ewald_prereg.md` was **FROZEN
> 2026-07-24** (with §12 Errata E1). The Step-4 primitive `invz_dipole_ewald.m` was subsequently
> **built, reviewed, and committed** (`086d102..fcb3031`, additive; production default still
> `bruteforce`; suite 269/0/23). Consequently the "not freeze-ready", "No Ewald implementation …
> should begin", "Step 4 implementation remains unauthorized", and "does not authorize
> implementation" statements in the header below and in §10 are **historical snapshots of the
> pre-freeze review and are no longer operative**. Current operative status and roadmap:
> `docs/invzp_ordered_1z_state.md`. Per the freeze rule, the original review text is retained
> unedited below this note rather than rewritten.

**Status: PREREG-REVIEWED DESIGN — design closure is complete, but this is not an implementation
authorization.** This revision incorporates the independent review of the completed Phase-1 report,
the derivation calibration, every live user-facing projected-spectra path, and review of the draft
`docs/invzp_ewald_prereg.md`. That preregistration exists but is not freeze-ready until the
corrections recorded in §10 are incorporated. No Ewald implementation or physical-anchor run should
begin until the corrected preregistration is approved and frozen.

This design executes the frozen escalation clause in
`docs/invzp_phase1_quadrature_prereg.md`: failure of the corrected half-open construction under
real-space-cutoff refinement escalates the coupling calculation to an Ewald or equivalent
convergence-accelerated dipolar sum before Phase 3.

---

## 0. Determination and immediate recommendation

Phase 1 reached a valid hard stop:

- The half-open grid removes the duplicate-face defect and passes the structural grid checks.
- Both Γ policies still fail coupling-multiset grid refinement and offset agreement.
- The brute-force spherical real-space sum also fails the frozen `dpRng=40→50` convergence gate.
- `Jcc0` being comparatively stable does not rescue the construction: Phase 1 shows that a Γ-only
  scalar can look stable while the full coupling distribution remains materially unstable.

The immediate next step is therefore **amend and freeze the Ewald preregistration**, followed by an
additive, opt-in primitive. Do not implement against the current draft preregistration, continue
increasing `N` or `dpRng` in the same brute-force construction, choose a Γ policy yet, or begin the
ordered exact-`h` Phase-3 audit.

Phase 1 supports the Ewald escalation; it does not by itself prove that every observed offset effect
is caused by one particular surface term. Treat conditional convergence/truncation of `MF_dipole`
as the leading, testable diagnosis. The Ewald candidate must demonstrate that it removes the failed
numerical behavior.

## 1. Scope and non-negotiable constraints

### 1.1 Initial scope

- Add a new Ewald primitive for the **non-ODD projected path**.
- Leave `MF_dipole` untouched and permanently reachable for legacy parity and diagnostics.
- Leave the tensor path (`invzt_jq_tensor`), `MF_RPA_Yikai`, and `RPA_line` untouched.
- Keep the current production default on `bruteforce` until every gate in §6 passes.
- Do not modify or reinterpret the frozen Task-2 or Phase-1 pre-registrations.

The current ODD diversion is not automatically covered: `invz_jq_modes` sends `opts.odd` to
`invz_odd_blocks`, which calls `MF_dipole` directly for all Cartesian blocks. During this project:

- `opts.odd` continues to use the labeled brute-force backend;
- requesting `opts.dipole='ewald'` together with an active `opts.odd` must raise a clear
  unsupported-combination error rather than silently mixing backends; and
- adopting Ewald in `invz_odd_blocks` is a separate, later scope requiring the full `ca`, `cb`, and
  `cc` block validation.

### 1.2 Acceptance cannot be “physics anchors only”

Primitive α/cutoff convergence and BZ quadrature convergence test different things. An Ewald sum can
be converged for each sampled `q` while a finite BZ quadrature remains grid- or offset-sensitive.
The Phase-1 result itself demonstrates why Γ scalars and physical anchors alone are insufficient.

Production adoption therefore requires all of:

1. primitive formula and boundary-condition validation;
2. Ewald α/real-cutoff/reciprocal-cutoff convergence;
3. reciprocal-gauge periodicity;
4. the applicable frozen Phase-1 **grid-refinement and eight-offset gates** rerun with Ewald;
5. the frozen Phase-2 physical anchors; and
6. a target-path Jensen/HMF ordered-state acceptance check, because the legacy ordered anchor tests
   only the bare-MF branch; and
7. a predeclared Γ-policy decision tree consistent with the implemented boundary convention.

The existing Phase-1 tolerances remain the acceptance tolerances for the BZ multiset and offset
checks unless a new tolerance is justified and frozen **before** candidate results are inspected.

## 2. Architecture and interfaces

### 2.1 Primitive

Add a pure, deterministic primitive:

```matlab
[dip, counts, geom] = invz_dipole_ewald(q, a, tau, eopts, geom)
```

Contract:

- `q` is in reciprocal-lattice units, exactly as in `MF_dipole`.
- `dip` uses `MF_dipole`'s convention and units:
  `dip = -Σ'_R T(R+d) exp[-i q·(R+d)]`, in Å⁻³.
- The output has the same singleton/squeeze behavior as `MF_dipole`; for the production four-site
  basis it is `[3,3,4,4]`.
- Sublattice Hermiticity is
  `dip(:,:,m,n) = conj(dip(:,:,n,m))`.
- `eopts` contains `alpha`, `r_cut`, `g_cut`, and the only initially supported boundary label,
  `boundary='conducting_k0_omitted'`. This means: omit the exact reciprocal `k=0` term and add no
  Ewald surface term inside the primitive (conducting/tinfoil exterior). Lorentz reconciliation and
  the caller-owned ellipsoidal demagnetization term are separate operations described in §4.
- `geom` caches only q-independent data: validated lattice metadata, pair-specific real-space
  vectors/factors, a reciprocal candidate set covering the canonical q domain below, and index
  shells. Quantities involving `qbar+G`, including the Gaussian and the `|qbar+G|<g_cut` mask, are
  computed per `q`.
- Reusing an incompatible `geom` must error. Its fingerprint includes `a`, `tau`, all Ewald controls,
  the boundary label, the q-reduction convention, and a schema version.
- `counts` is diagnostic, not a legacy scalar: it is a struct with
  `real_pair` `[ntau x ntau]` retained-vector counts, `recip_candidates` (cached candidate count),
  and `recip_used` (q-specific count after the `|qbar+G|<=g_cut` mask).

The q-domain contract is fixed rather than implicit:

1. Reduce every input row `q` deterministically to the half-open reciprocal parallelepiped
   `qbar = q-K`, `qbar_i in [-1/2,1/2)`, with integer `K=floor(q+1/2)`.
2. Evaluate the Ewald sums at `qbar`. The cached reciprocal candidate set must contain every `G`
   that can satisfy `|qbar+G|<=g_cut` for any `qbar` in that full half-open cube.
3. Restore the exact total-displacement gauge before returning:
   `dip_nm(q)=exp[-i 2*pi*K·(tau_m-tau_n)] dip_nm(qbar)`.

This makes extended-zone paths and reciprocal-gauge tests complete without letting the first q
passed to a reusable geometry object determine its reciprocal shell.

`r_cut` and `g_cut` are physical Cartesian cutoffs (Å and Å⁻¹), not integer box extents. For each
sublattice pair, the real-space generator must include every `R` for which `|R+d|≤r_cut`; the
reciprocal generator must satisfy the canonical-domain completeness rule above.

### 2.2 Projected-path options

Use a namespaced integration option. The unchanged legacy request is:

```matlab
opts.dipole = 'bruteforce';       % current default until adoption
% opts.ewald is absent
```

An explicit Ewald candidate request is:

```matlab
opts.dipole = 'ewald';
opts.ewald  = struct( ...
    'alpha', ..., ...
    'r_cut', ..., ...
    'g_cut', ..., ...
    'boundary', 'conducting_k0_omitted');
```

Do **not** reuse top-level `opts.alpha`: it already denotes the ellipsoid aspect ratio used by
`ellipsoid_demagn`. Ewald α lives only at `opts.ewald.alpha`.

The Ewald backend and every control must be represented losslessly in the coupling-cache key.
Do not build a cache key from a weak weighted-sum hash alone for the new schema. The cached payload
must also store and validate the exact q array, lattice parameters, backend, Ewald options,
demagnetization options, exchange, and schema version before reuse.

Option combinations are strict rather than permissive:

- any `opts.ewald` field without `opts.dipole='ewald'` must error;
- `opts.dipole='bruteforce'` with an Ewald field must error rather than silently ignore
  them;
- `opts.dipole='ewald'` requires the complete frozen Ewald-control schema; and
- unknown backend or boundary labels must error.

An absent `opts.dipole` field remains the only implicit legacy-default request before adoption.

### 2.3 Integration surface

The change is larger than a one-line swap because the live code has backend-specific callers:

- `invz_projected/invz_jq_modes.m`: backend dispatch, Γ reconciliation, cache schema, metadata.
- `invz_projected/invz_bz_couplings.m`: forward `dipole`, `ewald`, grid convention/offset, and Γ
  policy so the parameterized acceptance wrappers can actually inject the candidate.
- `invz_projected/invz_jq_path.m`: forward the backend and replace the brute-force
  `dpRng`-dependent Γ snap logic. With converged Ewald, nonzero q points are evaluated directly;
  only an exactly Γ-equivalent endpoint needs an explicit path-direction prescription.
- `invz_projected/invz_spectra_map.m`: expose and forward the BZ backend, Ewald controls,
  `gridConvention`, `gridOffset`, and `gammaPolicy`; validate matching provenance when precomputed
  `Jnu`/`info` are supplied.
- `invz_projected/invz_spectra_qpath.m`: use one consistent backend for both its BZ medium and its
  path dispersion. It must forward the BZ grid policy to `invz_bz_couplings`, forward
  `dipole`/`ewald` to `invz_jq_path`, and reject precomputed `info` whose backend/controls conflict
  with an explicit request.
- `invz_projected/invz_run_spectra.m`: add user-facing `dipoleBackend`, Ewald-control,
  grid-convention, production-offset, and Γ-policy knobs and pass the same resolved configuration
  into both field-map and q-path routes. After any future default flip, an explicit
  `dipoleBackend='bruteforce'` remains supported.
- Phase-1 Ewald wrapper/driver: reuse the committed grid and summary machinery without altering the
  frozen brute-force report.

`info` must identify the backend and controls used. Backend-specific geometry must not be exposed
under an ambiguous contract that causes `invz_jq_path` to call `MF_dipole` on Ewald metadata.
The backend-agnostic q-path export is named `info.Jpath_base_cc`, not `Jreg_cc`:

```text
bruteforce:
  Jpath_base_cc = -gfac*dip_sphere_cc(0) + Jex_cc(0)

Ewald:
  Jpath_base_cc = -gfac*dip_reg_cc(0) + Jex_cc(0) - lorz*ones(4)

both:
  J(q->0,qhat) = Jpath_base_cc
                 + gfac*(4*pi/Vc)*(1/3-qhat_c^2)*ones(4).
```

The Ewald subtraction in the metadata export is a normalization for reusing the existing q-path
formula; it is **not** an extra term in the primitive or in the exact-Γ production matrix. Export
the latter separately as `info.Jgamma_cc` so neither object is called “regular Γ” ambiguously.

## 3. Mathematical convention — pin this before coding

Let the rows of `a` be real-space primitive vectors, `b` the reciprocal rows used by `MF_dipole`,
`Vc=|det(a)|`, and `τ_n` the Cartesian basis positions (`tau*a` in the live code). In the equations
below, `q` and `G` are Cartesian reciprocal vectors; the API input is converted as
`q_cart=q_rlu*b`. For the existing `(n,m)` matrix convention,

```text
d_nm = τ_m − τ_n
x    = R + d_nm
T_ab(x) = (3 x_a x_b − |x|² δ_ab)/|x|⁵
        = +∂_a∂_b (1/|x|)
dip_ab,nm(q) = −Σ'_R T_ab(R+d_nm) exp[−i q·(R+d_nm)] .
```

The sign `T=+∂∂(1/r)` is essential. Writing `T=-∂∂(1/r)` reverses the reciprocal and self
bookkeeping.

With Fourier transform
`f_hat(k)=∫f(x)exp(-i k·x)dx`, `k=q+G`, and the **total-displacement phase used by
`MF_dipole`**, a direct formula for the primitive's `dip` convention is:

```text
dip^(r)_ab,nm(q)
  = −Σ'_R ∂_a∂_b[erfc(alpha |x|)/|x|] exp(−i q·x),
    x=R+d_nm, |x|≤r_cut

dip^(G)_ab,nm(q)
  = +(4π/Vc) Σ_{G:k≠0}
      (k_a k_b/|k|²) exp[−|k|²/(4 alpha²)] exp(+i G·d_nm),
    |k|≤g_cut

dip^(self)_ab,nm
  = −δ_nm δ_ab 4 alpha³/(3 sqrt(π)).
```

The returned intrinsic candidate is
`dip = dip^(r) + dip^(G) + dip^(self)`, with the `k=0` reciprocal term omitted.

Important consequences:

- The reciprocal phase is `exp(+i G·d_nm)`, not
  `exp[-i(q+G)·d_nm]`, for `MF_dipole`'s total-displacement Fourier convention.
- Under a reciprocal vector `K`,
  `dip_nm(q+K)=exp(-i K·d_nm) dip_nm(q)`. Raw sublattice matrices are gauge-covariant; their
  eigenvalue multisets are periodic.
- The positive reciprocal projector and negative self term above apply to `dip=-ΣT`. If an
  implementation instead accumulates `D=+ΣT`, both signs reverse.
- The expanded screened real-space coefficients must be derived from the displayed Hessian and
  checked symbolically/numerically. A named paper alone is not enough because papers use differing
  Hamiltonian, Fourier, and tensor signs.

Primary references for the boundary and cutoff analysis are de Leeuw, Perram & Smith,
Proc. R. Soc. A **373**, 27–56 (1980), DOI
[`10.1098/rspa.1980.0135`](https://doi.org/10.1098/rspa.1980.0135), and Wang & Holm,
J. Chem. Phys. **115** (2001), DOI
[`10.1063/1.1398588`](https://doi.org/10.1063/1.1398588). The implementation plan must map their
symbols and signs explicitly to the equations above.

## 4. Γ, Lorentz, and demagnetization reconciliation

### 4.1 Keep two distinct “Γ removals” separate

- **Ewald `k=0` omission plus the explicit absence/presence of a surface term** defines the
  macroscopic boundary convention of the dipolar lattice sum. The initial
  `conducting_k0_omitted` branch has no primitive surface term.
- **BZ P-drop** removes an exact Γ row from a numerical Brillouin-zone average.

They are different operations and must never share a flag or be used to justify one another.
P-complete and P-drop remain symmetric candidate quadratures until Phase 2.

### 4.2 RESOLVED ON PAPER — the Ewald `conducting_k0_omitted` branch adds `0` at Γ

The anticipated directional structure is confirmed: at nonzero `q→0` the Ewald reciprocal G=0 term is
`+(4π/Vc)qhat_a qhat_b`, giving `J_cc(q→0,qhat) = J_reg,cc − gfac(4π/Vc)qhat_c²` — identical in its
directional part to the existing convention's `+gfac(4π/Vc)(1/3−qhat_c²)`. The two differ only by the
isotropic `lorz = gfac·4π/(3Vc)` in the regular part, and that one constant is now decided
in scratch calibration (`docs/invzp_ewald_derivation.md`, Calibration D):

- `v'·[−gfac·dip_reg,cc]·v + 4·J12 = 0.00642166` reproduces the physical
  `Jcc0 = 0.00642444` to **2.77e-6**
  — which is the brute-force dpRng=30 target's *own* conditional-convergence error (a converged sphere
  matches Ewald to **3.2e-7**) — while the "adds-`lorz`" alternative is off by **4.88e-3**, three
  orders larger. Unambiguous.
- Per element, `dip_reg,ab(Ewald) = dip_sphere,ab(brute) − (4π/Vc)(1/3)δ_ab` (verified over all 16
  sublattice pairs to 6.2e-6 cc / 7.8e-6 aa). **The isotropic Lorentz term is already inside the Ewald
  regular part.**

**DECISION: the Ewald branch adds `0` at Γ** (retaining the caller's `+lorz` would double-count), for
**both the cc and aa channels** (`Jaa0` verified identically, Δ=1.4e-6). Integration consequence: the
Ewald branch simply **omits** the `if Γ: Jcc += lorz` line; demag (`dm_cc`/`dm_aa` via
`ellipsoid_demagn`), `Jshape_cc`, and the demag-aware `Jaa0` ride on top exactly as today, so the
four-bullet semantics (`Jcc0`/`Tc(B=0)` demag-invariant, demag exported) are preserved. The
brute-force branch keeps its `+lorz` and remains bitwise identical in every pre-existing numerical
output and `info` field.

### 4.3 Γ reconciliation calibration — passed provisionally; formal rerun required

The six checks below were run in scratch calibration
(`docs/invzp_ewald_derivation.md`, Calibration D) and support the design decision. The scratch code
was not retained, so these results are **not** prospective Gate-C acceptance evidence. Step 3 must
freeze exact samples, tolerances, retained result artifacts, and the independent oracle; the
repository implementation must rerun every check after the freeze.

1. `k=0` boundary term derived in the exact `dip=−ΣT`, total-displacement convention (reciprocal phase
   `e^{+iG·d}`, self `−δ_nm δ_ab·4α³/(3√π)`), independently confirmed against a scalar-Coulomb Ewald
   reference to 9.2e-9.
2. `dip_reg(0)` = `dip^(r)(0)+Σ_{G≠0}dip^(G)(0)+dip^(self)` (direction-independent), distinct from the
   directional `q→0` limit.
3. The scratch run approached Γ along `a*`, `b*`, `c*`, and two generic rays and reported the same
   regular matrix after subtracting the analytic projector. The formal test must distinguish the
   isolated reciprocal `G=0` term from the full analytic remainder: raw non-Bravais off-diagonal
   blocks may contain an allowed O(q), odd-in-q imaginary part.
4. Formally test the isolated `G=0` contribution against
   `(4π/Vc)qhat_a qhat_b exp[−q²/(4alpha²)]` for every Cartesian component and sublattice pair. For
   the full tensor, subtract that exact finite-q term and either fit the allowed O(q) analytic term
   or use the paired even part `[R(q)+R(−q)]/2`, which must approach `R(0)` as O(q²). Report the odd
   part separately. At integration level, test the complete
   `δ_ab/3−qhat_a qhat_b` reconstruction, not only cc.
5. Macroscopic term couples only to the uniform mode (`v'Δv=1.81e-2` = predicted `4·(4π/Vc)qhat_c²`;
   the non-uniform block carries no `4π/Vc`-scale term, 1.2e-6).
6. Both `Jcc0` and `Jaa0` gated; demag-invariance and the `Jshape_cc`/`Jaa0` observable semantics
   retained.

The §4.2 algebraic convention is therefore fixed for preregistration: the Ewald branch adds `0` at
Γ; the brute-force branch stays numerically identical in all pre-existing outputs and fields.
Formal acceptance remains conditional on the post-freeze Gate-C rerun.

## 5. Pre-registration contract and current draft

The draft `docs/invzp_ewald_prereg.md` now exists. Amend it to satisfy this section and §10, then
freeze it before implementation or inspection of any Ewald-based physical anchor. Numerical
calibration used to choose safe cutoffs may inspect only primitive convergence residuals, not
`Σ_c`, `Tc`, `Bc`, ordered moments, or soft modes.

The pre-registration must fix:

### 5.1 Parameter ladder

- a deterministic rule for the default `alpha`, expressed in Å⁻¹;
- a predeclared α bracket around that default;
- at least three real cutoffs and three reciprocal cutoffs;
- separate one-axis cutoff tests: grow `r_cut` with a generous fixed `g_cut`, and grow `g_cut` with
  a generous fixed `r_cut`, so cancellation of truncation errors cannot fake convergence;
- a final joint-refinement comparison; and
- minimum accuracy guards no weaker than `alpha*r_cut>=4.5` and
  `g_cut/(2*alpha)>=4.5` (the coarsest calibrated ladder values);
- separate hard, deterministic resource bounds on real-pair vectors, reciprocal candidates,
  estimated working bytes, and other quantities that can actually be evaluated before allocation;
- any elapsed-time limit must name the reference environment and distinguish an operational timeout
  from a scientific convergence failure. Wall time cannot be preflighted exactly and must not make
  a mathematically valid construction backend- or machine-dependent; and
- a clear failure if either the requested tolerance or the deterministic resource bounds cannot be
  met.

### 5.2 Primitive metrics and q sample

Freeze raw-tensor and derived-coupling tolerances for:

- α independence;
- real- and reciprocal-cutoff independence;
- Hermiticity and realness at Γ;
- reciprocal-gauge covariance of raw matrices and periodicity of sorted eigenvalues;
- sublattice-pair conjugation;
- raw invariance under a common basis-origin shift, and gauge covariance (not raw equality) when an
  individual basis representative is shifted by a Bravais vector;
- Γ-ray reconciliation from §4; and
- agreement with an independent reference.

The q sample must enumerate exact Γ regular data, small-q magnitudes and rays, generic interior
points, BZ-face/edge/corner points, reciprocal-equivalent integer shifts, basis translations, and
all tensor components and sublattice pairs. “Raw-tensor REL” is not a sufficient metric name:
Step 3 must define every norm, absolute floor, relative denominator, comparison direction, and
zero/near-zero-component rule explicitly.

The screened-Hessian derivative check and the scalar-Coulomb finite-difference oracle also require
their own local scales and discretization-error contracts. A q-dependent whole-tensor scale is not
defined for a single screened real-space Hessian sample. A finite-difference comparison may use a
frozen FD floor only after demonstrating convergence across at least two step sizes (preferably
Richardson-extrapolated); it must not be described as passing the tighter analytic tensor metric
when its calibrated FD error is larger than that metric.

Use the frozen `J_ref` to express coupling-level tolerances. Do not validate only normalized
eigenvalues: a sign or sublattice phase error can be hidden by selected spectra.
The retained independent scalar-Coulomb oracle, its finite-difference steps, and exact inputs and
outputs are formal Gate-B deliverables; the unretained scratch calibration cannot substitute for
them.

### 5.3 BZ numerical gates

For the half-open candidate, rerun:

- reciprocal periodicity (Phase-1 item 2);
- multiset grid refinement at the pre-registered N ladder (item 5);
- all eight offsets at every N and the frozen finest rung (item 6);
- both P-complete and P-drop, reported symmetrically; and
- the original normalized summary and energy-scalar tolerances.

Items 1, 3, and 4 are structural grid/weight checks already established for the half-open builder,
but their cardinality, Γ count, and weight sum should still be emitted as cheap regression
assertions. The Ewald cutoff ladder replaces the `dpRng` ladder; it does **not** replace the N or
offset ladders.

Because Ewald is expected to be much cheaper per converged q sweep than the old large sphere, prefer
the complete declared cross-product. Any economization must be written into the pre-registration
before results are inspected.

The production offset is fixed now as `[0 0 0]`; the other seven offsets are translation-sensitivity
diagnostics. Production adoption is allowed only if the frozen offset-agreement gate passes, so
this canonical choice cannot select a materially different finite-grid answer.

For any P-complete `invz_sigma_crit` anchor at the unshifted grid, Step 3 must assert the number of
entries removed by its singular-pole guard. The expected contract is one excluded uniform Γ
eigenvalue under P-complete and zero under P-drop; any other count is a hard diagnostic failure,
not a warning to ignore.

## 6. Validation and adoption gates

Validation proceeds in this order:

### Gate A — algebra and primitive

- Test the Hessian coefficients against numerical derivatives of `erfc(alpha*r)/r`.
- Test signs against a small explicitly enumerated screened real-space sum.
- Test the reciprocal phase and gauge transformation, not raw periodic equality.
- Test canonical half-open q reduction, extended-zone gauge restoration, and reciprocal-candidate
  completeness at the worst BZ corners and declared integer shifts.
- Test the self term on same-sublattice blocks and its absence off-sublattice.
- Reject every boundary label except `conducting_k0_omitted`; verify that no surface or demag term
  is present in the primitive.
- Validate every `counts` field against explicit enumeration and enforce the frozen resource caps
  before allocation.
- Test α and both cutoffs independently, as frozen in §5.

### Gate B — independent numerical reference

The brute-force conditional sum is not the sole oracle. Require a retained independent reference
implementation and machine-readable results, such as:

- differentiating a separately implemented scalar Coulomb Ewald potential with respect to basis
  displacements; or
- a second trusted dipolar-Ewald implementation with an explicit convention map.

For a finite-difference oracle, freeze at least two derivative step sizes, the extrapolation rule,
and a separate FD acceptance floor supported by the observed step convergence. Report both the
oracle error and the tighter primitive self-convergence metric; do not conflate them.

A spherical brute-force extrapolation at generic q bounded away from Γ is useful as a secondary
cross-check only. It must state the summation shape and extrapolation model. It cannot validate the
Γ boundary term by itself.

### Gate C — Γ/Lorentz/demag reconciliation

Pass every proof-driven test in §4.3 and remove the current brute-force `dpRng` trust-radius behavior
from the Ewald q-path branch. Nonzero q must be evaluated directly; exactly Γ uses the path's local
direction and the analytic limit. Test `Jpath_base_cc` and `Jgamma_cc` separately, and verify that
the primitive contains no surface/demag term while caller-level demagnetization preserves the
documented `Jcc0`/`Jaa0`/`Jshape_cc` semantics.

### Gate D — corrected Phase-1 numerical acceptance

The Ewald half-open candidate must pass the frozen item-2, item-5, and item-6 gates and the cheap
item-1/3/4 assertions in §5.3. Failure here means the original blocker is not resolved; do not move
to physical adoption merely because α convergence passes.

### Gate E — frozen Phase-2 physics anchors

Run parameterized wrappers, leaving legacy tests unchanged:

- LiHoF4 `Σ_c` and `Jcc0` (`test_invz_sigma_crit.m`);
- zero-field `Tc` (`test_invz_critical.m`);
- 310 mK critical field `Bc` (`test_invz_critical.m`);
- legacy qualitative **bare-MF** ordered state and soft mode (`test_invz_ordered_phase.m`); and
- a new target-path **Jensen/HMF** ordered acceptance test using the same `ordered_1z='jensen'`
  route as `invz_spectra_map`.

Use each anchor's documented tolerance; do not replace them with a blanket “<0.3%” claim. Exact
grid-dependent values in `invz_odd_anchors.m` remain legacy/ODD diagnostics, not selection gates.
The bare-MF test explicitly requires a finite ordered-side gap and therefore cannot stand in for the
current symptom. Step 3 must freeze field points and criteria for the Jensen target test before
candidate results are seen. At minimum it must require accepted FM states (`phase_1z=1`, finite
`Sigma0`, `m_1z`, and `D_ord`) at predeclared ordered fields, stable PM states on the other side,
and no unexplained masked columns at every declared point. It must also freeze `T`, field direction,
field grid, frequency grid and units, `eta`, `hyp`, `peak_wmin`, `transverse_mf`, solver options,
grid/backend provenance, and the exact finite/positivity/ordering inequalities.

Do not require one monotonic “softening” sequence through a far-PM field: a mode normally softens
toward the transition and may harden again on the PM side, while `D_ord` is undefined there. Gate
ordered-side `m_1z`, `D_ord`, and `Epeak` only on ordered points; gate PM stability with
`phase_1z=2`, finite positive `crit_pm`, and its own spectral criterion. Any continuity claim must
use explicitly bracketed fields and static masses, not a qualitative peak trend across phases.
Failure leaves the user-facing spectra defect unresolved even if the legacy bare anchor passes.

### Gate F — Γ policy and production selection

After theory, numerical gates, and physical anchors are available, apply this decision tree without
scoring candidates by closeness to an anchor:

1. A policy is eligible only if it independently passes Gates D and E.
2. If neither is eligible, stop with no production selection.
3. If exactly one is eligible, select it.
4. If both are eligible, first compare numerical multisets on **common support** (remove the full
   exact-Γ row from P-complete and apply the same normalization used by P-drop) as an implementation
   consistency check. Do not require raw P-complete/P-drop extrema—especially `max(Jnu)`—to agree:
   one policy contains the exact Γ maximum and the other deliberately does not. Compare physical
   anchors and genuinely integrated summaries using separately frozen cross-policy tolerances. If
   those materially disagree, stop for a Γ-theory/finite-size audit. Otherwise select **P-drop** as
   the predeclared measure-zero critical-pole regularization and current production-compatible
   tie-breaker.
5. Any conflicting physical anchors, unexpected singular-entry count, or Jensen-target failure is
   a stop; there is no majority vote or post-hoc tolerance change.

Then:

- freeze the half-open grid, uniform weights, production offset `[0 0 0]`, selected Γ policy,
  `conducting_k0_omitted` boundary convention, α, and cutoffs;
- only then flip the supported non-ODD projected-path default to Ewald;
- recompute and re-freeze `Bc_PM` and the exact derived audit artifacts enumerated in Step 3;
- repopulate caches under the new schema; and
- resume Phase 3.

On any gate failure, keep the brute-force default, report the first failed gate and diagnostic, and
remain short of the 3A/3B/4 paths.

## 7. Files and regression coverage

Expected implementation surface:

- `invz_dipole_ewald.m` — primitive and private geometry builder.
- `invz_projected/invz_jq_modes.m` — opt-in backend, Γ mapping, cache schema, backend metadata.
- `invz_projected/invz_bz_couplings.m` — candidate-option and Γ-policy threading.
- `invz_projected/invz_jq_path.m` — backend threading and Ewald exact-Γ directional handling.
- `invz_projected/invz_spectra_map.m` and `invz_projected/invz_spectra_qpath.m` — consistent
  backend/grid provenance forwarding.
- `invz_projected/invz_run_spectra.m` — user-facing backend/grid controls.
- additive Ewald preregistration, driver/report, and parameterized anchor wrappers.

Required regression coverage:

- Existing calls with no `opts.dipole` remain bitwise identical in their numerical outputs and
  every pre-existing `info` field before the production flip. The complete `info` struct is not
  `isequaln` because additive backend/provenance fields are required.
- Explicit `opts.dipole='bruteforce'` remains numerically bitwise identical in those same
  pre-existing outputs and fields indefinitely.
- Brute-force and Ewald caches never collide.
- Cache cold/warm results are identical for each backend.
- `opts.odd` remains labeled brute force and rejects an Ewald request until separately supported.
- All three `invz_jq_modes` outputs and the relevant `info` semantics are covered.
- Field-map and q-path tests prove that BZ and path couplings cannot silently use different
  backends, including when precomputed `Jnu`/`info` are supplied.
- Existing legacy regression tests remain unchanged; candidate tests are additive.

## 8. Execution order (corrected — blockers close BEFORE the prereg is frozen)

This supersedes the earlier ordering, which froze the prereg before resolving §9 — the contradiction
the user flagged (§9 required blocker closure before approval, yet the old order froze the prereg
first). The Γ convention (old Gate C) is resolved on paper here, not after integration.

1. **Close the design blockers on paper** — screened Hessian, k=0/Lorentz mapping, independent
   reference, parameter ladders/tolerances, wrapper API, q-path contract. **[DONE — §9,
   `docs/invzp_ewald_derivation.md`, `docs/invzp_ewald_integration_map.md`.]**
2. **Amend and approve this design.** **[DONE — prereg-reviewed 2026-07-24; the design is approved,
   but the Step-3 draft still requires the §10 amendments before it can be frozen.]**
3. Amend the reviewed draft and freeze `docs/invzp_ewald_prereg.md` only after the §10 corrections
   are incorporated.
4. Implement the opt-in primitive plus retained Gate-A/B tests and the primitive-level parts of
   Gate C; no production integration and no default change.
5. Integrate the opt-in `invz_jq_modes`/BZ/q-path/spectra backend, cache schema, and provenance;
   run the remaining Gate-C and legacy-regression tests. Production default stays `bruteforce`.
6. Rerun the frozen Phase-1 grid/offset gates with both Γ policies on the Ewald couplings (Gate D).
7. Run the frozen legacy and Jensen-target physical anchors (Gate E), apply the deterministic Γ
   decision tree (Gate F), flip the default only on success, and recompute/re-freeze `Bc_PM` and
   the enumerated derived audit artifacts; then resume Phase 3.

## 9. Design blockers — RESOLVED ON PAPER (formal gates remain prospective)

The review-corrected contracts below are closed well enough to freeze a prospective preregistration.
The appended scratch numbers are calibration evidence only; “resolved on paper” does not mean that
Gates A–F have passed.

1. **Screened Hessian — DONE.** `g_ab = P(r)(3x_ax_b−r²δ)/r⁵ + Q(r)x_ax_b`,
   `P=erfc(αr)+(2αr/√π)e^{−α²r²}`, `Q=(4α³/√π)e^{−α²r²}/r²`; finite-difference-confirmed to 4.9e-9,
   0 sign mismatches, α→0 → bare tensor.
2. **Boundary/k=0/Γ mapping — DONE ON PAPER.** The only initial boundary is
   `conducting_k0_omitted`: no primitive surface term, Ewald adds **0** at Γ (§4.2), and
   caller-owned demagnetization remains separate. The isotropic Lorentz is already in `dip_reg`,
   for both cc and aa; brute force keeps its existing `+lorz`.
3. **q-domain/cache completeness — DONE ON PAPER.** Reduce to the canonical half-open reciprocal
   cube, cache a reciprocal candidate union complete over that cube, and restore the exact
   sublattice gauge for extended-zone q. The reduction convention is fingerprinted.
4. **α/cutoff ladders + tolerances — CALIBRATED.** Default
   `α₀=√π/Vc^{1/3}=0.2684 Å⁻¹`, `r_cut=C_r/α₀`
   (`C_r=5.5→20.49 Å`), `g_cut=C_g·α₀` (`C_g=11→2.95 Å⁻¹`); default residual raw 3.6e-13 / coupling
   2.1e-13; ladders `C_r∈{4.5,5.0,5.5}`, `C_g∈{9,10,11}` bracket ≤1e-8 on separate axes; α bracket
   `[0.161,0.403]`. Accuracy guards are corrected to 4.5 on both dimensionless axes; Step 3 must
   freeze deterministic enumeration/memory caps and treat any wall-clock limit as an
   environment-qualified operational bound, not a preflightable scientific gate.
5. **Independent reference design (Gate B) — DONE ON PAPER.** Scalar-Coulomb Ewald potential
   differentiated w.r.t. basis displacement (phase on R only — a genuinely separate code path)
   agreed in scratch calibration to 9.2e-9 (off-diagonal) and isolated the self term to 4.1e-9.
   The formal oracle and its outputs must be retained and rerun after preregistration.
6. **Physics-anchor wrapper API — DONE** (`invzp_ewald_integration_map.md`, Blocker 5). The anchors
   inline their own `qVec_generator`+`invz_jq_modes` (they do NOT call `invz_bz_couplings`) and today
   run the legacy endpoint-inclusive grid — so Ewald validation flips grid convention (→half-open) and
   backend (→ewald) together. Plan: a shared `invz_anchor_couplings.m` helper + sibling
   `test_invz_*_ewald.m` files with `_pcomplete`/`_pdrop` variants, each reusing its anchor's own
   documented tolerance. The helper omits optional grid/backend fields when absent, so its default
   path is genuinely legacy-identical. A separate Jensen/HMF target gate is mandatory.
7. **q-path Γ/metadata contract — DONE** (map, Blocker 6). Hoist the explicitly normalized
   `info.Jpath_base_cc` out of `invz_jq_modes` so `invz_jq_path` reads it instead of calling
   `MF_dipole(...,info.geomD)`; export `info.Jgamma_cc` separately. The `dpRng` trust radius becomes
   backend-conditional (exactly-Γ-only under Ewald).
8. **User-facing integration contract — DONE ON PAPER.** `invz_spectra_map`,
   `invz_spectra_qpath`, and `invz_run_spectra` forward one consistent provenance set; a mixed
   Ewald-BZ/brute-force-path calculation must error.
9. **`Bc_PM` re-freeze — adoption-time procedure** (§8 step 7), not a design blocker. After the Γ
   policy and parameters are frozen and the default flips, recompute `Bc_PM`; update
   `invz_task2_matrix_enumerate.m`'s `Bc_PM` comments and three derived field constants,
   `tests/test_invz_task2_matrix.m`'s frozen root/field assertions, a superseding Task-2 freeze
   record/report, and all coupling/Task-2 caches whose provenance contains the old construction.

**Two clarifications folded in:** (a) **exchange is out of scope** — `invz_jq_modes` keeps its
short-range `sign(J12)·|J12|·ex(3,3)` term via `exchange.m` (absolutely convergent, untouched by
Ewald); only the dipolar term changes backend. (b) the Ewald primitive returns the diagnostic
`counts` struct defined in §2.1, not an ambiguous legacy-style scalar `nN`.

## 10. Review disposition of `docs/invzp_ewald_prereg.md` (2026-07-24)

The draft preregistration is comprehensive but **must not be frozen as currently written**. Amend
these prospective rules without inspecting Ewald physical-anchor results:

1. **Separate analytic and finite-difference metrics.** Its `M_T` whole-tensor scale is undefined
   for a single screened-Hessian sample, and the calibrated scalar-oracle error (~`9.2e-9` Å⁻³) is
   larger than `M_T` at the observed tensor scale. Define a local Hessian scale and a
   step-converged FD floor; do not claim the FD oracle passes `M_T`.
2. **Split deterministic resource gates from timing.** Vector counts and estimated bytes may be
   hard pre-allocation gates. A `120 s` wall-clock limit cannot be evaluated before allocation and
   has no meaning without a reference machine/process configuration; make it an operational
   timeout/reporting rule or freeze the benchmark environment separately.
3. **Correct the Γ-ray test.** Test the isolated finite-q `G=0` projector directly. For the full
   non-Bravais tensor, allow the analytic odd-in-q imaginary term and apply the q² requirement only
   to the paired even remainder or to a justified extrapolation.
4. **Correct the Γ-policy comparison.** P-complete and P-drop cannot be required to have the same
   raw `max(Jnu)` at finite N because one includes the exact Γ row. Use common-support numerical
   comparisons and separately frozen integrated/physical cross-policy criteria.
5. **Fully specify the Jensen/HMF gate.** Freeze all calculation knobs and exact inequalities.
   Ordered-side softening must not be required to remain monotone through a far-PM point; gate
   ordered and PM observables on their respective branches and use static masses for any
   continuity statement.

After those five amendments and a consistency check, the preregistration can be considered for
freeze. Until then, Step 4 implementation remains unauthorized.

This document does not authorize implementation, default changes, ODD/tensor-path adoption, or any
Phase-3/3A/3B/4 work. This review authorizes only amendment and re-review of
`docs/invzp_ewald_prereg.md`.

> **[Superseded 2026-07-24 — see the STATUS ADDENDUM at the top of this file.]** The prereg was
> amended per the five points above and **FROZEN** (2026-07-24, §12 Errata E1); the Step-4 primitive
> was then built, reviewed, and committed (`086d102..fcb3031`, additive, production default still
> `bruteforce`). The "remains unauthorized" / "does not authorize implementation" language in this
> section is a historical pre-freeze snapshot and is no longer operative. Step 5 (opt-in integration,
> default unchanged) is the current authorized phase.
