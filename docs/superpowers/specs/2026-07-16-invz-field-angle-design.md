# 3D field-direction (misalignment) knob for `invz_run_spectra`

Date: 2026-07-16. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-16). Scope: **scalar stage** (leading-order,
small-misalignment), with the full-tensor rigorous route deferred and documented
(see §8).

## Problem

`invz_run_spectra.m` sweeps a scalar field magnitude and treats it as a purely
**transverse** field along the crystal **a**-axis: every solver in the chain
turns the scalar `Bx` into the literal vector `[Bx 0 0]` (ordering axis = **c** =
z). Real LiHoF4 experiments have a small **misalignment** of the nominally
transverse field toward **c**; the resulting longitudinal component `Bz` rounds
the sharp quantum phase transition into a crossover and is a known source of
experiment/theory discrepancy. The user wants to port the field-angle knob from
`LiReF4_MF_Yikai.m` (which parametrizes the direction via `theta`/`phi`) into the
1/z spectra driver so tilted-field spectra can be computed.

## Physics background — why longitudinal field is *not* trivial in 1/z

At the **MF** and **RPA** levels a longitudinal field is trivial (as
`LiReF4_MF_Yikai.m` + `MF_chi.m` + `RPA.m` already do): both are **full-tensor**
methods. MF just diagonalizes `H(B)` for any `B = [Bx By Bz]` and
self-consistifies `⟨Jx⟩,⟨Jy⟩,⟨Jz⟩`; RPA builds the full 3×3 `χ0(ω)` from all
transitions and inverts `χ = [χ0⁻¹ − 𝒥(q)]⁻¹`. A tilt just changes the numbers
in the same matrices.

The **1/z self-energy is a different kind of object.** Jensen's `Σ(ω)` is a
*scalar* built from one transition — the field-split Ising doublet — and it
renormalizes **only the cc (c-axis) component**: the code literally computes
`χ̃0_cc = χ0_cc/(1+Σ)` and never forms `χ0_xx`, `χ0_xz`, or any off-diagonal.
Every ingredient (`α`, `γ`, the `λp` Matsubara sums, the ordered `α_m`, the
elastic sector) is a projection onto the doublet parameters
`(Δ, M²=|⟨0|Jz|1⟩|², n01, m=⟨0|Jz|0⟩)`.

Consequences for a tilted field:

- **The moment-carrying machinery is already the general case.** A longitudinal
  `Bz` adds `−gL·µB·Bz·Jz` to `H0`, inducing `⟨Jz⟩ ≠ 0` at every field. This is
  structurally identical to how the *ordered* branch already carries a
  longitudinal molecular field `hz`. And `invz_sigma_ordered` (the m≠0
  self-energy) **reduces exactly to** `invz_sigma` as `m→0` (`alpha_m → 0`);
  `invz_chi0z` already folds the induced moment into the elastic term via
  `Jexp`; the sum rule already uses the *variance* `JzJz_fluct = ⟨Jz²⟩−⟨Jz⟩²`.
  So the induced-moment case flows through the existing ordered path with no new
  self-energy — the transverse-only case is its `m=0` special case.

- **What the scalar pipeline structurally cannot represent** is the
  **cross-channel (xz) dressing**. A tilt mixes the eigenstates, so the *true*
  `χ0` acquires `χ0_xz` cross terms; a full-tensor RPA inversion mixes them
  (even with diagonal couplings, since `χ0⁻¹` of a full matrix is full). The
  scalar cc pipeline discards `χ0_xz` entirely. Therefore the scalar treatment
  captures the **leading** longitudinal effect (transition rounding via the
  induced moment `m` and the re-split doublet), correct to **O(tilt²)** for a
  small misalignment, and **omits** the cross-channel correction, which grows
  with tilt. Genuine rigor for arbitrary tilt requires the full-tensor route
  (§8).

## Design decisions (resolved with the user)

1. **Scope of the field direction:** support arbitrary direction including the
   c-tilt (the transition-rounding misalignment), *staged* — see decision 4.
2. **Knob shape:** full 3D via **two angles** `theta_c` (tilt out of the ab-plane
   toward c) and `phi_ab` (in-plane azimuth a→b), **plus** adding the transverse
   `hy` mean-field self-consistency to `invz_single_ion` so the in-plane azimuth
   is treated rigorously (not just the x-component). Convention matches
   `LiReF4_MF_Yikai.m`.
3. **The ODD extension does not help here.** `odd_implementation_plan.html`
   Tiers 1–2 are *deliberately scalar-cc* (they fold off-diagonal *dipole*
   couplings back into the cc channel; orthogonal to an external longitudinal
   *field*). Only its deferred Appendix A (A0+A1, full-tensor) would make 1/z
   full-tensor and hence make arbitrary field direction natural.
4. **Staged delivery:** implement the **scalar port now** (this spec); document
   the full-tensor **A0+A1** bridge as the rigorous follow-up that would
   supersede it (§8).

## Solution overview

Thread a full field **vector** through the solve chain in place of the hardcoded
`[Bx 0 0]`, in a **backward-compatible** way: a tiny helper maps a scalar to
`[B 0 0]` (today's behavior) and passes a 3-vector through untouched. The driver
builds the vector from a fixed misalignment direction × the swept magnitude. When
the field has a longitudinal component, the solve is routed exclusively through
the moment-carrying (ordered) path with the "is it spontaneous?" gate disabled,
since the moment is now field-induced. The single-ion mean field gains a `hy`
channel for the in-plane azimuth. All changes are inert when the field is along
x, so every existing benchmark and test is preserved.

## Components

### 1. New helper `invz/invz_field_vec.m`

```matlab
function B = invz_field_vec(B)
%INVZ_FIELD_VEC Normalize a field argument to a 1x3 row [Bx By Bz] in Tesla.
%   Scalar b -> [b 0 0] (transverse along a; the historical convention).
%   3-element vector -> row passthrough. Anything else errors.
```

Used at every site that currently writes `[Bx 0 0]`. A scalar still maps to
`[Bx 0 0]`, so all scalar callers are unchanged.

### 2. `invz/invz_single_ion.m` — add the `hy` transverse mean field

Currently the MF loop self-consistifies only `hx = Jxx0·⟨Jx⟩`. Add a `hy =
Jyy0·⟨Jy⟩` channel (`Jyy0 = Jxx0` by tetragonal a≡b symmetry; optional
`opts.Jyy0` override), **guarded on `B(2) ≠ 0`** so a field with no y-component
leaves the loop bit-for-bit as today.

- `H = H0 − hx·Jx − hy·Jy − hz·Jz` inside the fixed-point loop and in the final
  recompute.
- `jy = ⟨Jy⟩`, `hy_new = Jyy0·jy`; extend the convergence test to
  `dmf = max(|hx_new−hx|, |hy_new−hy|, |hz_new−hz|)`.
- Seed `hy = 0`; mix identically to `hx`.
- Header note: `hy` active only for an in-plane-azimuthal field; `Jyy0 = Jxx0`
  (tetragonal).

### 3. Field-vector plumbing (the 5 `[Bx 0 0]` leaf sites)

Each becomes `invz_field_vec(B)` and accepts a scalar-or-vector `B`:

- `invz_twolevel.m:7` — keep the `m=0` gate (it protects the strict-paramagnet
  path; only reached when `Bz=0`).
- `invz_twolevel_ordered.m:13` — rebuild the doublet with the **full vector** +
  the fixed MF `hz`. No double-counting: external `Bz` lives in `H0`, `hz` is
  the MF piece only, exactly as in the `order` solve.
- `invz_solve_point.m:16,19` — pass the vector to `invz_single_ion` /
  `invz_twolevel`.
- `invz_solve_point_ordered.m:34,45` — pass the vector to `invz_single_ion`
  (`order`) / `invz_twolevel_ordered`; add `forced_moment` (see §4).
- `invz_chi_realaxis.m:37` — the paramagnet fallback vector (the ordered path
  already reuses `pt.si` built with the full vector).

`invz_critical.m`, `invz_critical_T.m`, `invz_critical_T0field.m`,
`invz_run_phase_diagram.m` keep passing scalars → `[Bx 0 0]` → unchanged. The
angle knob is **not** added to the phase-diagram driver (a longitudinal field has
no sharp `Bc`/`Tc`).

### 4. `invz/invz_solve_auto.m` + `invz_solve_point_ordered.m` — longitudinal routing

`invz_solve_auto`:

- Compute `Bvec = invz_field_vec(B)`.
- **`Bvec(3) == 0`** (transverse, incl. pure in-plane rotation): keep today's
  logic verbatim — ordered-first (spontaneous moment), paramagnet fallback.
  In-plane fields keep `⟨Jz⟩ = 0` (a transverse field splits the Ising doublet
  into `(|+⟩±|−⟩)/√2` combinations, each with zero diagonal `Jz`; only a
  longitudinal field biases it), so the paramagnet path and its `m=0` gate stay
  valid.
- **`Bvec(3) ≠ 0`**: route *exclusively* to `invz_solve_point_ordered` with
  `opts.forced_moment = true`; `phase = 1` when converged, else `phase = 0`. Do
  not attempt the strict-paramagnet solver (its `m=0` gate would reject the
  induced moment).

`invz_solve_point_ordered` gains `opts.forced_moment` (default false): when set,
`pt.is_ordered = pt.converged` (skip the `|m0| > mtol` spontaneous-moment test),
because with an explicit `Bz` any moment is physical and the ordered self-energy
reduces smoothly to the paramagnet form as `m→0`. Routing threshold: treat
`|Bvec(3)|` above a small absolute tolerance (proposal `1e-9` T) as longitudinal.

### 5. Drivers `invz_spectra_map.m` / `invz_spectra_qpath.m`

- Accept a field **vector** (or a direction + magnitude) and pass it down. In
  `invz_spectra_map`, `one_field` builds `Bvec = mag·dhat`; the `parfor` slices
  over magnitudes as today (x-axis stays the swept `|B|`).
- Phase labels: under `Bz≠0` the FM/PM "V" is a **rounded crossover** and
  `phase` is always the moment-carrying branch — update the display strings /
  header comments accordingly (e.g. "field-induced moment (crossover)" rather
  than "ferromagnet"). The RPA overlay (`Sigma=0`) already uses the ordered-style
  `pt0` for `phase==1`, so no overlay-path change is needed.
- `invz_spectra_qpath` solves the medium once at `Bvec` (magnitude `Bq` ×
  direction) and is otherwise unchanged.

### 6. Driver knob `invz/invz_run_spectra.m`

Two scalars near the existing knobs:

```matlab
theta_c = 0;   % deg: tilt of the field OUT of the transverse ab-plane toward c (Ising axis)
phi_ab  = 0;   % deg: in-plane azimuth from a (x) toward b (y)
```

Direction (matches `LiReF4_MF_Yikai.m`):

```matlab
tc = theta_c*pi/180;  pa = phi_ab*pi/180;
dhat = [cos(tc)*cos(pa), cos(tc)*sin(pa), sin(tc)];   % unit; fixed across the sweep
```

Fed to both the field-sweep view (each `fields(k)·dhat`) and the q-path view
(`Bq·dhat`). `theta_c = phi_ab = 0` reproduces today's `[|B| 0 0]` exactly.
Header documents: the convention; that the sweep x-axis is the total magnitude
`|B|` (c-component `|B|·sinθc`); and the O(tilt²) / cc-only validity note with a
pointer to the deferred full-tensor follow-up (§8).

## Backward compatibility

- A scalar field anywhere → `[B 0 0]` → identical Hamiltonian, identical path.
- `hy` guarded on `B(2)≠0` → field-along-x callers unchanged.
- `Bz=0` → `invz_solve_auto` keeps today's exact branch logic.
- Phase-diagram / critical drivers and all their tests untouched.

**Non-negotiable:** with `theta_c = phi_ab = 0` the full existing test suite and
every published benchmark (Σc, Tc(0), Hc, soft-mode minimum) must reproduce
within their current tolerances.

## Testing

New tests under `invz/tests/` (fast unless noted); mostly exact-symmetry:

1. **Regression:** full suite green with the knob at defaults; a direct check
   that `invz_single_ion(ion,T,[Bx 0 0])` is unchanged by the `hy` addition.
2. **hy mean field (C4):** `[0 B 0]` spectrum ≡ `[B 0 0]` spectrum;
   `⟨Jy⟩|[0 B 0]` ≡ `⟨Jx⟩|[B 0 0]` (a↔b equivalence).
3. **In-plane azimuth:** `⟨Jz⟩ = 0` for any `Bz=0` field; single-ion spectrum
   invariant under `phi_ab → phi_ab + 90°`.
4. **Longitudinal:** ±`Bz` mirror symmetry (Z2) of `χ''_cc`; no NaN masking
   across the former `Bc`; peak energy stays finite through the crossover;
   `pt.sumrule_rel` stays small; `theta_c → 0` continuously recovers the
   transverse result.
5. **forced_moment:** an induced-moment point at small `Bz` converges with
   `pt.is_ordered = true` and matches the transverse result as `Bz→0`.

## 8. Deferred — full-tensor rigorous follow-up (A0+A1)

The scalar stage omits the cross-channel (xz) dressing (see Physics background).
The rigorous route is the **A0+A1** bridge of
`odd_implementation_plan.html` Appendix A, with the ODD blocks set to zero
(ODD is orthogonal to the field-angle need):

- **A0** — full-tensor RPA parity layer: build `𝒥(q)` as `[12,12,nq]`
  (Cartesian⊗sublattice) and evaluate `χ = [χ0⁻¹ − 𝒥(q)]⁻¹` using the full 3×3
  `χ0` that `invz_chi0z` already returns. Under a tilt this captures the
  `χ0_xz` cross-channel response the scalar pipeline drops.
- **A1** — projected-1/z bridge: keep the full-tensor, all-136-state
  propagation but apply the scalar `Σ_c` only to the dominant longitudinal
  sector via a transition mask: `χ̃0 = χ_dom/(1+Σ_c) + χ_rest`, then the full
  inversion. This is the multilevel "dominant transition renormalized, weak
  transitions kept at RPA" approximation, and it **supersedes** the scalar port
  for arbitrary tilt (also subsuming retardation).

This is a research-scale build (new tensor-RPA path, Schur-complement RPA-parity
tests, convention/double-counting validation) and is intentionally out of scope
for this spec. When implemented it would make the scalar routing of §4 a special
case.

## Notes / open items

- The two-level identification uses the lowest two states `E(1),E(2)`. Under a
  strong tilt a third level can approach the doublet, degrading the two-level
  approximation — a validity limit of the scalar stage, documented, not
  code-enforced.
- `forced_moment` longitudinal threshold (`1e-9` T) is a proposal; final value
  set during implementation against the routing tests.
