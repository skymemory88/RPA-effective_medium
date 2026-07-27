# invZ design rationale (consolidated)

This file consolidates the durable design rationale from `docs/superpowers/specs/` (six
design documents, 2026-07-08 through 2026-07-25) into one standalone record, ahead of
deprecating the `docs/superpowers/` skill directory. It captures **decisions, reasoning,
rejected alternatives, binding cautions, and frozen numeric/formula conventions** — not
the step-by-step implementation or test-authoring process that got them built. Compiled
2026-07-27. The original spec files remain recoverable from git history after this
consolidation; nothing here supersedes them as the historical record, but this file is
the intended durable reference going forward.

**Repo-wide conventions** (stated once here; several specs restate them): `G = -chi`
(units meV⁻¹); ferromagnetic couplings are positive `J`.

---

## 1. (2026-07-08) Two-regime phase-boundary search for `invz_run_phase_diagram`

**Decision.** Split the paramagnetic-side boundary trace at a regime boundary into:
- **Low-T** (`T` at/below the split): unchanged fixed-`T` bisection of `Bc(T)` via the
  existing `invz_critical`.
- **High-T** (boundary above the split, approaching `Tc0`): invert the roles — bisect
  the critical *temperature* `Tc(B)` at fixed field with a new mirror function
  `invz_critical_T`.

**Why.** Near the classical critical point (`T → Tc0 ≈ 1.74 K, B → 0`) the boundary runs
nearly parallel to the field axis. A fixed-`T` (vertical) cut crosses it at a glancing
angle: tiny temperature errors become huge `Bc` errors, brackets become hard to
establish, and the last few high-`T` points are ill-conditioned or fail outright. A
fixed-`B` (horizontal) cut crosses the same near-vertical boundary transversally, so it
is well-conditioned exactly where the vertical cut is not. The ordered/paramagnet
classifier is identical in both directions: `pt.crit` non-finite or ≤ 0 ⇒ ordered side.

`invz_critical_T(ion, Bx, Jnu_flat, opts)` bisects `pt.crit` over temperature at fixed
transverse field; `opts.window = [Tlo Thi]` (default `[1.0 2.0]` K), `opts.tol` (default
`0.01` K). Bracket assert requires ordered at `Tlo`, paramagnetic at `Thi`; the error
names the likely cause (`Bx` exceeds `Bc(Tlo)`). **`invz_critical.m` itself was not
modified** (it is validated by the slow test suite).

### Amendments (same day, during implementation — these govern over the original design text)

1. **Measured defect, small fields:** the paramagnetic solve develops non-convergence
   patches near the boundary; the non-finite-⇒-ordered classifier then biases `Tc(B)`
   upward by **+0.04–0.05 K at 0.2–0.3 T** (gentler `max_outer`/`mix_outer` does **not**
   cure it). At `B ≥ 0.5 T` the method is clean: `Tc(0.5 T) = 1.777 K`, just below the
   closed-form `Tc0 = 1.7795 K` on the **same 16³ q-grid**. Binding clarification: the
   spec's original `1.74 K` baseline is the Richardson-*extrapolated* value — the wrong
   comparison for a 16³ computation; the undershoot direction was right once compared
   grid-consistently.
2. **`Bs` floored at 0.5 T**, default `[0.5 0.75 1.0 1.25 1.5]`. The `0 < B < 0.5 T`
   segment spans only ~4 mK and is represented on the plot by the closed-form `Tc0`
   endpoint alone.
3. **The `Tsplit` knob was REMOVED** (user decision: trimming `Ts` already achieves the
   same effect, so a separate knob is redundant). The `Tc(B)` bisection window is fixed
   inside the driver as constants `Tlo = 1.0` K / `Tmax = 2.0` K, with documented
   constraint `max(Bs) < Bc(Tlo) ≈ 2.8 T`. Default `Ts` now ends at 1.6 K.
4. Near-zero-field test renamed `test_tc_small_field`: `B = 0.5` T, window `[1.0 2.0]` K,
   bounds `(1.70, 1.79)` K.

**Frozen defaults** (as amended): `Ts = [0.05 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6]`
(trimmed at the split); `Bs = [0.5 0.75 1.0 1.25 1.5]`; `Tlo = 1.0` K, `Tmax = 2.0` K
constants (not knobs). One flat `parfor` over all `numel(Ts)+numel(Bs)` jobs; each keeps
the existing try/catch-to-NaN pattern.

**Why `Bs` never reaches near zero:** (1) `invz_twolevel` raises `invz:degenerateDoublet`
when the field-induced doublet splitting `Δ < 1e-4` meV (`Bx ≈ 0`); (2) R2007's
small-`Bx` caveat — the two-level `Σ` overestimates the `Tc` suppression when the
doublet splitting is comparable to or smaller than the hyperfine width, so `Tc(B → 0)`
slightly undershoots the closed-form `Tc0`. This is expected physics, not a bug.

**Rejected alternatives:** slanted/normal cuts or predictor–corrector boundary
continuation (serial, more code, overkill for this problem); generalising
`invz_critical` into a line-segment bisection (would refactor a validated function for
no current need). Neither was pursued; the solver layer (`invz_solve_point` and below)
was explicitly out of scope for this change.

**Scope limit / later status.** The projected `invz_critical_T` referenced here was
itself rewritten on **2026-07-09**, one day after this spec, after a "rugged-boundary"
failure (near the boundary, EMT self-consistency suffers critical slowing down and
naive bisection latches onto spurious sign flips). That rewritten algorithm — valid-
samples-only T-grid voting, highest-T ordered→para crossing, regula-falsi refinement,
adaptive `Tc0`-anchored window — is the version later transplanted to the tensor branch
in §5 below, which also documents and closes two classifier gaps this projected
rewrite still carries (exact-zero votes, re-entrant lower leg) as a flagged, un-actioned
follow-up on the projected side.

---

## 2. (2026-07-08) Display-only energy-unit knob (`eUnit`) for `invz_run_spectra`

**Decision.** Add a knob `eUnit = 'GHz' | 'meV'` near the top of `invz_run_spectra.m`.
Conversion is **display-only** and lives entirely in the driver, reusing the existing
(previously unused) `invz_const().Gh2mV` constant (GHz → meV). The frequency grid `w` fed
into `invz_spectra_map` for the actual solve is **never** touched — only plotted values
are rescaled.

**Why this shape.** The driver has two mutually exclusive plot branches selected by field
count (`numel(fields) <= sliceMax`): a line-slice overlay and a 2D colormap via
`invz_plot_spectra_map.m`. A unit knob covering only one branch would be a silent no-op
under the script's own default settings (`fields = linspace(0,9,201)`, `sliceMax = 6`
lands in the colormap branch) — so both branches needed the conversion.

**Frozen formula:**
```matlab
C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end
```
Line-slice branch plots `w*eScale` with `xlabel(eLabel)`. Colormap branch builds a
display-only shallow copy `Splot = S; Splot.w = S.w*eScale;` and calls
`invz_plot_spectra_map(ax, Splot, Splot.chiz, ttl, eUnit)` — a **new 5th argument**.
`invz_plot_spectra_map.m` uses it only to build `ylabel(ax, sprintf('\omega (%s)',
eUnit))`; it adds **no conversion arithmetic**, so callers passing numeric data in
whatever unit they like keep working unchanged.

**Binding caution:** invalid `eUnit` must error *immediately after the knob*, before the
potentially 15–25 minute sweep runs — not discovered only at plot time.

**Scope limit.** Explicitly excludes converting any non-plotted quantity (`eta`'s doc
comment, printed `Sigma0` diagnostics stay in meV), and does not touch
`invz_spectra_map.m` or any solver-layer file. No automated test was added — this is a
display-only labeling knob with no physics/numerics to unit-test, and
`invz_plot_spectra_map.m` had no existing test coverage at the time. Verification was
manual: axis numbers must scale by `1/C.Gh2mV` (~241.8×) with the label updating.

---

## 3. (2026-07-16) c-axis field-misalignment (tilt) knob for `invz_run_spectra`

Scope: **scalar stage, c-axis tilt only** (`phi_ab = 0`); azimuthal support and full
tensor propagation are deferred follow-ups (§8 of the spec, summarized below). Revised
twice the same day after external review; findings verified against code before
acceptance (table at the end of this section).

**Problem.** The driver treats a scalar swept field as purely transverse along **a**:
every solver turns scalar `Bx` into `[Bx 0 0]` (ordering axis = **c** = z). Real LiHoF4
experiments have a small misalignment of the nominally transverse field toward **c**;
the resulting longitudinal `Bz` rounds the sharp QPT into a crossover. Convention
follows `LiReF4_MF_Yikai.m` (external file, `.../Mean Field/LiReF4/LiReF4_MF_Yikai.m`,
lines 60–66): the **formula** `Hx = B·cosθ·cosφ, Hz = B·sinθ` is authoritative; that
file's line-60 comment ("angle from c-axis") contradicts its own formula and is to be
ignored.

**Why longitudinal field is not trivial in 1/z (physics background).** At MF and RPA,
a longitudinal field is trivial — both are full-tensor methods (MF diagonalizes `H(B)`
for any `B`; RPA inverts the full 3×3 `χ0(ω)`). The **1/z self-energy is a different
kind of object**: Jensen's `Σ(ω)` is a *scalar* built from one transition (the
field-split Ising doublet) and renormalizes **only** the cc component
(`χ̃0_cc = χ0_cc/(1+Σ)`); it never forms `χ0_xx`, `χ0_xz`, or any off-diagonal. Every
ingredient (`α`, `γ`, the `λp` Matsubara sums, the ordered `α_m`, the elastic sector) is
a projection onto doublet parameters `(Δ, M²=|⟨0|Jz|1⟩|², n01, m=⟨0|Jz|0⟩)`.

Consequences: the moment-carrying machinery is *already the general case* — a
longitudinal `Bz` induces `⟨Jz⟩ ≠ 0` structurally identically to the ordered branch's
molecular field `hz` (`invz_sigma_ordered` reduces exactly to `invz_sigma` as `m → 0`),
so the induced-moment case flows through the existing ordered path with **no new
self-energy**.

**Accuracy statement (corrected per review — binding).** `χ_cc` is even in the tilt
angle, so the entire tilt effect starts at `O(θ²)` — and the scalar port keeps only
*part* of the `O(θ²)` coefficient (nonperturbative doublet re-splitting + induced-moment
effects) while omitting another same-order part (the `χ_zx·χ_xz` cross-channel
contribution a full tensor inversion would mix in). The defensible claim is: **exact at
zero tilt; uncontrolled relative error in the tilt-induced change, beginning at
`O(θ²)`.** Near the former critical point the response to the conjugate field is
non-analytic, so Taylor counting does not apply there at all — the supported angle
range must be established **numerically** (Σ=0 tensor-reference test below), never by
power counting.

### Design decisions

1. **Scope to c-axis tilt (`theta_c`) only; `phi_ab` descoped.** With `B64s ≠ 0` an
   x-field already induces a small perpendicular `⟨Jy⟩ ≈ −0.069` at 4 T (verified) with
   no feedback channel — a "rigorous azimuth" would need a two-channel `(hx, hy)`
   transverse mean field that is exactly C4-consistent, which is **incompatible with a
   bit-for-bit x-field baseline**. The transverse MF stays the legacy x-only
   approximation, now documented as such. **Rejected:** resurrecting a `By ≠ 0` guard
   (angle-discontinuous, C4-inconsistent) — explicitly declined even as a deferred
   option.
2. **Longitudinal routing through the moment-form self-energy**, sign-aware branch
   selection (default `mz_seed = +1` converges to the *metastable* anti-aligned branch
   for `Bz < 0` — verified failure mode).
3. **The ODD extension does not help here.** `odd_implementation_plan.html` Tiers 1–2
   are deliberately scalar-cc (orthogonal to an external longitudinal field). Only its
   Appendix A route is relevant: A0 is the rigorous tensor-RPA layer; A1 remains a
   projected, dominant-transition 1/z approximation (not itself a rigorous tensor 1/z).
4. **Staged delivery:** scalar c-tilt now; azimuth/two-channel MF and A0(+A1) tensor
   propagation later (deferred, §8).

### Frozen conventions and components

`invz_field_vec(B)` normalizes to a `1×3` row `[Bx By Bz]` (Tesla): scalar → `[B 0 0]`
(historical convention) or exactly 3 elements; anything else (NaN/Inf, complex, empty,
wrong length) errors `invz:fieldVec`. **Sign-aware seed** in `invz_single_ion.m`:
`mz_seed = sign(B(3))` when `B(3) ≠ 0` (else `+1`; `opts.mz_seed` still overrides) —
verified digit-for-digit at `B = [2 0 −0.01]` T: seed `+1` gives `⟨Jz⟩ = +4.815`
(metastable), seed `−1` gives `−4.86686` (exact Z2 mirror). **Variational free energy**
`si.F_mf = −kT·ln Tr e^{−βH} + hx²/(2·Jxx0) + hz²/(2·J0z)` (double-counting correction
restoring `−½J⟨J⟩²` per channel, on the **unshifted** spectrum `si.E0`). **Binding
caution:** the naive shifted-spectrum comparison is *wrong* — at `T = 0.31 K`,
`B = [2 0 −0.01]` it mis-ranks the branches (`−5.15e-8` vs `−2.74e-8` meV), while
`F_mf` correctly ranks the aligned branch lower (`−21.47664574` vs `−21.46963936` meV).

**Routing rule**, `opts.bz_tol` (default `1e-9` T), single-sourced everywhere:
`|Bvec(3)| <= bz_tol` → zero the component, run today's transverse logic verbatim;
`|Bvec(3)| > bz_tol` → route **exclusively** to `invz_solve_point_ordered` with
`opts.forced_moment = true` (never `invz_solve_point`, whose two-level gate rejects
`m ≠ 0`). **`forced_moment` acceptance is explicitly non-circular**: bypass the early
paramagnetic return; require `si.mf_converged`; assert
`sign(si.Jexp(3)) == sign(Bvec(3))` (retry once with a mirrored seed on mismatch, else
return non-converged with a warning); run the outer EMT⇆Σ loop as today; only then set
`pt.is_ordered = true` plus machine-readable
`pt.moment_branch = 'spontaneous' | 'field_induced'`. **Every** return path of
`invz_solve_point_ordered` (including early/failed returns) populates the same field
set (`is_ordered, converged, Sigma0, crit, si, tl, m0, moment_branch`) — callers never
probe a missing struct member.

**Map API:** `opts.field_dir = [1 0 0]` (unit-normalized internally); `fields` stays a
vector of nonnegative magnitudes; `S.Bvec` carries the dead-band-normalized vectors
actually used. `opts.solve_opts` passes extra solver knobs through, with
`J0eff`/`Jxx0`/`hyp` **reserved** (error `invz:solveOpts` if present — the driver owns
them). **Driver knob:** `theta_c = 0` (deg, tilt out of the ab-plane toward c),
`dhat = [cosd(theta_c) 0 sind(theta_c)]`; axis label `'B_x (T)'` → `'|B| (T)'`.

**Σ=0 tensor-reference validation (the numerical support-range test, measurement-driven,
not assumed).** Compares the scalar-chain `Σ=0` result against a full 3×3 Cartesian RPA
built from `invz_chi0z` with `diag(Jaa0, Jaa0, Jsel)`, at `theta_c ∈ {0, 0.5, 1, 2, 5}`
deg, `ion.demag = 0`. **Measured baseline fact:** at `θ = 0` with `B ∥ a`, the `yz` cross
channel of `χ0` is symmetry-allowed and LARGE (`max|χ0_yz|/max|χ0_zz| = 0.183` at 6 T,
`hyp = false`, vs `xz/zz = 2.8e-3`) — this raw discrepancy is present at **all** angles
(a pre-existing property of the scalar cc pipeline) and must be reported as baseline,
**never** gated as tilt error. An L2 spectral metric at the production linewidth
(`eta = 5e-3`) is dominated by sub-linewidth peak misalignment (measured 44%
discrepancy at `θ = 0` coexisting with a physically negligible 3.5 µeV peak shift) — a
metric instability, not physics — and even a baseline-differenced L2 tilt metric stays
dominated by the zero-tilt `yz`-induced peak-position offset `δ₀` leaking through at the
`δ₀/η` scale (`≈ 0.11` vs `δ₀/η = 2.2/20 µeV` at 6 T; `≈ 0.28` vs `7.7/20 µeV` at 2 T).
**Final support criterion, on peak observables instead:** `θ > 0` is supported when, at
every tested field, `dE_peak <= max(0.02·Epeak_ten, eta)` AND
`eps_amp = |max_w χ''_sc − max_w χ''_ten| / max_w χ''_ten <= 10%` (both `max_w` over
`w >= peak_wmin`, excluding the hyperfine line). The L2 metrics remain reported
diagnostics with the artifact explanation; **lineshape fidelity is explicitly not
claimed** for the scalar stage.

**Backward compatibility (non-negotiable):** scalar field anywhere → `[B 0 0]` →
identical Hamiltonian. At `theta_c = 0` the routing dead band zeroes nothing (already
zero) and the full existing suite plus every published benchmark (Σc, Tc(0), Hc,
soft-mode minimum) must reproduce **bit-for-bit**.

**Deferred follow-ups (§8):** azimuthal field + two-channel transverse MF (requires
iterating `hx` AND `hy` for every direction including the default `[Bx 0 0]` case,
changing the default x-field model and demanding benchmark revalidation, in exchange
for exact C4 symmetry); full-tensor A0(+A1) — A0 the rigorous `[12,12,nq]`
Cartesian⊗sublattice tensor-RPA layer capturing `χ0_xz` exactly, A1 a projected 1/z
bridge (`χ̃0 = χ_dom/(1+Σ_c) + χ_rest`), still a dominant-transition approximation. Either
would supersede this stage's scalar routing for arbitrary tilt.

**Review-verified cautions:** the `one_field` overlay's longitudinal branch previously
called `invz_twolevel` *outside* its try block, so a failed longitudinal point
(`phase = 0`) would raise uncaught and abort a `parfor` — now an explicit contract
(masked/RPA-only column, never a crash). A vector `B` formatted with a bare `%.3f` is
silently mis-formatted by MATLAB's format recycling — fixed with `mat2str`-style
formatting everywhere a vector might reach a `sprintf`.

---

## 4. (2026-07-18) Tensor-branch drivers: `invzt_run_phase_diagram` / `invzt_run_spectra`

**Decisions.**
1. **Scope: tensor-only drivers, PM-side, as-is.** Ship `invzt_run_phase_diagram.m` and
   `invzt_run_spectra.m` in `invz_tensor/`, honestly scoped to what
   `invzt_solve_point`/`invzt_critical`/`invzt_chi_realaxis` supported *at the time*
   (no ordered-phase branch, no temperature-cut finder — both real physics/numerics
   gaps, explicitly deferred, and in fact closed one day later by §5 below). One
   narrow **non-physics** module fix is in scope: a one-line projection bug in
   `invzt_chi_realaxis`'s explicit-q branch, below.
2. **`invz_peak_energy.m` moves to `invz_common/`** (`git mv`, zero logic change) so
   both branches call one shared copy — matching the precedent already set moving 16
   single-ion functions there for the same reason.

**Problem.** `invz_projected` ships runnable drivers; `invz_tensor` had only point-solver
primitives, forcing users to hand-roll the sweep loop. Porting the projected drivers
verbatim was not possible: `invzt_solve_point` only ever returns a paramagnetic fixed
point (**LOCKED convention 7**), and there was no temperature-cut root finder yet.

### Verified facts (binding — checked against source across two verification passes + external review)

- **`imag(out.chi_uniform(3,3,:))` is the correct, already-positive `χ''`** — no sign
  flip anywhere in the tensor chain, matching the projected chain's convention.
  **`out.chi_cc_q` was real by construction — a real bug (blocker F1), confirmed:** the
  explicit-q accumulation did `acc = acc + real(Xq(...))` — a transplant of a
  Matsubara-axis noise-cleaning idiom onto the real axis, where it deletes the entire
  dissipative part. No production code consumed tensor `chi_cc_q` at the time; the only
  existing test asserted size/finiteness (zeros silently passed).
- **`invzt_tc_pm_extrap` does not guard its `critfun`** (a bare call whose exceptions
  propagate and abort the whole call); the driver's local `invzt_crit_at` supplies that
  missing selectivity. **`invzt_solve_point` throws on a realistic sweep**, no internal
  try/catch: sweep-plausible identifiers are `invzt:a1ZeroField` (transverse field
  below `1e-6` T), `invz:degenerateDoublet` (splitting `< 1e-4` meV), `invz:orderedPhase`
  (`|m| > 1e-3`) — a *failed* PM fixed point deep in the ordered phase does **not**
  throw, it returns `pt.converged = false` cleanly.
- **`invzt:bracket` doubles as `invzt_critical`'s argument-validation identifier**: a
  malformed/reversed `Brange` raises the *same* id as a genuine no-crossing window, so a
  driver absorbing `invzt:bracket` per point must preflight `Brange` itself, or a typo
  becomes a silent all-NaN sweep. **`invzt_critical`'s nested sampler absorbs every
  non-`invzt:*` exception as an invalid/ordered vote** (its committed "non-convergence =
  phase signal" policy); the drivers' fail-loud guarantee is scoped so *the driver*
  rethrows every non-bracket error that *escapes* `invzt_critical`, while the finder's
  own broader classification stands as designed (narrowing it is a module-behavior
  change with its own test surface, declined as out of scope for a driver task).
- **`invz_odd_zero_field` requires `ion.demag == 0`** for every mode, solves **seven**
  variants per grid at `mode='full'`, and always builds its own legacy
  endpoint-inclusive `[-0.5, 0.5]` mesh — a **different quadrature convention** from the
  tensor driver's `'halfopen'` grid even at equal nominal N.
- **Strict Γ in an explicit q-list is accepted, not rejected**, but `invzt_jq_tensor`
  assembles Γ-equivalent points with the Lorentz cavity — the strict-`q=0` *observable*,
  not the `q→0⁺` intrinsic limit a dispersion plot wants. **Convention decision:**
  dispersion q-paths in drivers/recipes/smoke tests are **Γ-excluded**; every example
  starts at finite q.

### Component 0 — prerequisite fix, `invzt_chi_realaxis` explicit-q complex contract

(a) **Keep the explicit-q response complex** (fixes the blocker): preallocate
`out.chi_cc_q = complex(zeros(nq, nw))` and drop the `real()` projection in the
accumulation, `acc = acc + Xq(3*(s-1)+3, 3*(s-1)+3, iq)`. Header updated: `chi_cc_q` is
COMPLEX; `imag()` is the positive `χ''` — contract parity with the projected
`invz_chi_realaxis`.

(b) **Enforce the documented A1-only scope** at the callee, protecting every caller:
```matlab
ptmode = getf(pt, 'mode', 'a1');
if ~strcmp(ptmode, 'a1')
    error('invzt:realaxisMode', ...);   % A2/A3 continuation is a separate open item
end
```

(c) **Honor `pt.odd` at explicit q** (added post-execution after a second review found
the explicit-q continuation ignoring `pt.odd = false` gave a **17.2% response error**):
the odd-off Cartesian-off-diagonal zeroing rule was extracted into a shared helper
`invzt_odd_mask.m` (byte-identical semantics) and applied to the explicit-q `latq.Jt`
(the load-bearing site) and the `'gamma'` `Jpage` (a C2 no-op at Γ, ~1e-19); the
Cartesian-diagonal `'gamma_uniform'` `Jd` needs no mask by construction.

**Binding caution preserved from the review record:** the physical (full-Σ) call has a
**pre-existing near-resonance negative-`χ''` artifact** (`−312.9` beside a `+652.5` peak
at the test point), present identically in the untouched `chi_uniform` path *and* in the
projected reference — inherited, not introduced by this fix. The exact
`χ''(ω>0)`-non-negativity gate therefore runs only on a separate `force_sigma0 = true`
(bare-RPA) call, where passivity is empirically established (`min = 0.319 ≥ 0`).

### Driver decisions and frozen conventions

- **`invzt_run_spectra.m` scope:** PM side only (cannot sweep across `Bc` — no
  ordered-phase branch). `solve_opts.mode` must be `'a1'`, enforced at both the driver
  (`invzt_run_spectra:mode`) and the callee (`invzt:realaxisMode`); no `dress` knob (A3
  is dead under forced a1).
- **Error policy, deliberately asymmetric by view:** the field sweep does per-field
  *selective masking* (physics signals + invalid samples → NaN column with a console
  note, never an error); the q-path view **fails loud by design**
  (`invzt:qpathNotPM`, with `converged`/`crit`/`Sigma0` in the message) because the
  whole product hinges on the single solved point — an empty map would be worse than a
  loud failure.
- **`opts.sigma_floor`** (default `−0.5`) is a single-sourced existing `invzt_critical`
  option, read identically by both drivers and the proxy helper — never re-defaulted
  independently.
- **`phi_ab ~= 0` under `transverse_mf = 'legacy_x'`** is rejected
  (`invzt_run_spectra:transverseMF`) — legacy_x is x-only and C4-inconsistent for a
  rotated field; needs `'vector_ab'`.
- **q-path preflight** (`invzt_run_spectra:qpath` / `invzt_run_spectra:qpathGamma`) uses
  the module's own `invz_is_gamma_equiv` gate (stricter than a naive `abs(q) <= 1e-12`,
  since it also catches Γ-equivalent points) and runs before any lattice/solve work.

**Rejected alternatives, explicitly recorded:** running the cross-model anchor
comparison on the tensor driver's own `'halfopen'` grid convention, to reduce the
mismatch with the projected comparator — declined, since the anchor is already a
cross-model comparator (different model, different ODD treatment) and quadrature is
the smallest of its three differences; and unifying the two branches' drivers behind a
shared `invz_common` entry point / model toggle — declined "for now," since the adapter
work is only worth it once the tensor branch has FM + T-cut support (§5 supplies the
T-cut half one day later; FM/ordered support remains open).

**Scope limits (explicitly out of scope):** `invzt_critical_T` and any ordered-phase
tensor solver (the actual blockers to full parity — see §5, which closes the first of
these the next day); any change to `invz_projected`'s own drivers beyond the one
`git mv`; A2/A3 real-axis continuation; warm-starting `Sigma_seed` *across* sweep points
(the projected drivers do not do this either; `invzt_critical` warm-starts only within
one bisection).

**Superseded note:** this spec's original `invzt_run_phase_diagram.m` (Component 1,
field-cuts only, `Ts = linspace(1.0, 2.2, 25)`) is **superseded** by §5's two-regime
rewiring, written and landed the same day/next day — see §5's Component 5 for the
governing version.

---

## 5. (2026-07-18/19) Tensor-branch temperature-cut finder: `invzt_critical_T` + two-regime driver

This spec extends §4 (chronologically after it: it explicitly supersedes §4's
phase-diagram driver Component 1, and its own revisions run through 2026-07-19). User
request: fixed-field `Tc(B)` search "just as in the projected branch," conditional on
the tensor/projected cores being fundamentally different enough to justify duplicating
the approach — condition affirmatively established (full 12×12 tensor RPA + min-eig
criticality vs. the projected scalar cc channel; measured **+0.016 K** boundary
difference at the identical grid). **Approach A** (full mirror of the projected
algorithm via shared helpers) was explicitly selected over alternatives.

### Decisions

1. **Mirror the projected algorithm's core idea (the 2026-07-09 rugged-boundary fix,
   §1 above), with its inherited classifier gaps closed here — not there.** Two gaps,
   confirmed by hand-simulation of the projected classifier: (i) an exactly-sampled
   root (`crit == 0`) is mis-read as two sign changes with no returnable crossing,
   producing a spurious `multipleCrossings` warning plus a bracket error; (ii) the
   lower leg of a re-entrant region (`[+,+,-,-]`) makes the classifier `break` instead
   of sliding toward the physical high-T paramagnet. Both are fixed in a new pure
   decision core, `invzt_tc_pick`, and are explicitly flagged as a possible **projected**
   -side follow-up — not touched there in this work (validated code, own test surface).
2. **Pure decision core:** `invzt_tc_pick(cv)` takes only the ascending-T vector of
   valid crit votes and returns `act ∈ {'zero','bracket','up','down'}` plus a crossing
   count `ncross` (strict sign flips + exactly-critical runs, a zero-run counted once).
   `'zero'`: the sampled root itself is the boundary. `'bracket'`: refine between the
   highest-T ordered→para pair. `'up'`: the top voter is ordered — the boundary (or the
   re-entrant lower leg) lies above the window. `'down'`: every voter is paramagnetic.
   No solves, no state — millisecond-testable on synthetic votes.
3. **An explicit `opts.window` is a HARD bound**: one grid pass, no sliding,
   `invzt:bracket` if it contains no returnable root. Only the adaptive
   (`Tc0`-anchored) mode slides, up to 9 window attempts total. (This is a deliberate
   improvement over the projected finder, whose explicit-window override still slides.)
4. **Shared helpers (Approach A):** `invz_refine_crossing.m` moved (`git mv`, zero
   content change) to `invz_common/`; the driver's local `invzt_crit_at` was promoted
   to a shared module file used by the finder, the driver, and the proxy.
5. **Namespaced tolerances, never one shared name in two units:** driver knobs `Btol`
   (tesla → `invzt_critical`'s `tol`) and `Ttol`/`Twidth`/`Tgridstep` (kelvin →
   `invzt_critical_T`), merged into per-finder opts only at the call boundary.
6. **Strict input validation via a hardened safe formatter, never raw `mat2str`.**
   MATLAB behavior confirmed: `mat2str(struct())` throws
   `MATLAB:mat2str:NumericInput` (masking the intended error identifier if used to
   *build* an error message); a vector on the right of `&&` throws
   `MATLAB:nonLogicalConditional`; `isfinite(1+1i)` is `true`. `invzt_str.m` (the
   module's existing safe formatter) itself fell through to raw `mat2str` and so
   inherited the same struct-throwing bug — hardened here with a class/size
   placeholder branch (e.g. `'<1x1 struct>'`) for anything `mat2str` cannot format.
7. **A unified cross-model driver stays declined** — the tensor side still has no
   ordered/FM solve; this work removes only one of the two blockers noted in §4
   (temperature-cut search). Revisit unification when FM lands.

### Frozen algorithm and conventions (`invzt_critical_T`)

```
crit(h)-style sample validity (single-sourced with the field-cut finder):
  ok = converged && isfinite(crit) && Sigma0 >= getf(opts,'sigma_floor', -0.5)
```
validity-only, deliberately with **no** `crit > 0` term — each consumer applies its own
phase logic. Selective catch on the sampler: `{invz:degenerateDoublet,
invz:orderedPhase, invzt:a1ZeroField}` → absorbed as an invalid/ordered vote; every
other identifier (misconfiguration) rethrows. Options (all finite positive real
scalars, else `invzt:tcOpts`): `opts.window` — hard `[Tlo Thi]` K bound, degenerate
spans (`Thi-Tlo <= 1e-9`) rejected, no sliding; `opts.Tc0` — **required** when `window`
is absent (the tensor branch has **no zero-field closed form** to fall back on, unlike
the projected branch's `invz:oddTc0` case) → else `invzt:tcAnchor`; `opts.width`
default `0.5` K, `opts.gridstep` default `1/30` K, `opts.tol` default `0.005` K.
Re-entrance indicators (`ncross`) **accumulate across adaptive window attempts**; more
than one across the whole scan warns `invzt:multipleCrossings` (never masked —
candidate re-entrance is physics to report) and returns the highest-T root. Reported
window (`out.window`) is always the **last actually-sampled** grid, captured only
*after* the floor-collapse guard — never a slide-mutated pair the user didn't actually
get sampled. Errors: `invzt:tcWindow` (malformed/below-floor `opts.window`),
`invzt:tcAnchor` (adaptive mode without a usable `Tc0`), `invzt:bracket` (no returnable
root — hard window without a crossing, adaptive attempts exhausted, or the window
collapsed at the `0.02` K single-ion solve floor).

**Execution finding E1 (a real algorithm gap none of the three code reviews caught,
caught instead by integration testing):** with every synthetic test green and every
code transcription byte-verified, the sampled `crit(T)` at fixed field flipped sign
**five times** in one round-trip test (returning 1.596 K vs. the field-cut's 1.4 K) —
cold-started tensor solves hop between A1 fixed-point branches near the boundary, and
the field-cut finder only avoids this because it warm-starts every sample. **Fix:** the
T-grid is sampled **descending** in T (paramagnet-first) with a rolling `Sigma` seed
re-interpolated onto each new temperature's Matsubara grid; state advances only on
valid samples, so an absorbed solve never poisons the chain. A caller-supplied
`opts.Sigma_seed` is still stripped (T-dependent Matsubara length fits at most one
sampled temperature) — the finder's own internal rolling seed is a separate mechanism.

**Execution finding E2 (a genuine, reported-not-fixed scope limit):** a controller
probe found **cross-finder path hysteresis** in a shallow-crit corner: at 8³/`dpRng`-15,
`odd` on, the T-cut's warm-tracked profile at `B = 1.5` T is smooth with a single clean
crossing near 1.35–1.36 K, yet the field-cut computes `Bc(1.4 K) = 1.916` T — the two
scan directions disagree by **~50 mK / ~0.4 T** in a near-degenerate corner at this
coarse grid (`|crit| ≲ 0.01` throughout — the A1 fixed point near the boundary depends
on the warm-start path taken to reach it). Per repo policy this was **reported, never
tuned**: the odd-off round-trip stayed exact, gates were left unchanged, and
re-examination at production 16³/`dpRng`-30 (where the shallow region should shrink)
was left as a follow-up.

**Operational gotcha:** the small-`Bx` proxy `Tc0` samples are **cold** (not
branch-tracked) — the proxy grid (`Ts_proxy`) must reach *well above* `Tc0`, or a
too-low anchor silently pushes the T-cuts onto the adaptive slide-up path. The odd-on
slow-gate floor was de-flaked `1.35` → `1.30` K (the original margin was flake-prone;
this is E2's physics, not a bug).

### `invzt_run_phase_diagram.m` — governing two-regime version (supersedes §4 Component 1)

Frozen default knobs: `Ts = linspace(0.4, 1.4, 11)` K (low-T field-cut branch);
`Bs = linspace(0.25, 1.5, 6)` T (T-cut fields); `Twin = []` (adaptive, anchored to the
small-`Bx` proxy `Tc0`) or an explicit hard `[Tlo Thi]`; `Brange = [0.05 6.0]` T;
`Btol = 0.02` T; `Ttol = 0.005` K; `Twidth = 0.5` K; `Tgridstep = 1/30` K; `gridN = 16`,
`gridConv = 'halfopen'`; `dpRng = 30`; `Bproxy = 0.05` T;
`Ts_proxy = 1.40:1/30:2.00` K. One flat `parfor` over both cut kinds; only
`invzt:bracket` is absorbed per point, everything else rethrows. The `show_proxy_anchor`
marker and the opt-in `show_projected_anchor` cross-model comparator carry the same
four caveats as in §4 (puts `invz_projected` on the path, requires `ion.demag == 0`,
solves 7 variants per grid, is a *different* ODD treatment — comparator, not the same
quantity).

**Scope limit, still standing:** PM side only — there is still no ordered-phase tensor
solve, so nothing below the boundary is computed. Unlike the projected T-cut finder
(which cannot bracket at all with ODD on, for lack of a metastable window), the tensor
A1 solver *does* converge metastable PM fixed points inside the ordered phase, so this
T-cut brackets with `odd` on or off alike.

**Out of scope (unchanged from §4):** ordered-phase tensor solve (still the FM
blocker to full driver unification); fixing the projected `invz_critical_T`'s own
inherited exact-zero/re-entrant-lower-leg gaps (flagged to the user as a candidate
follow-up, not touched); warm-starting *across* different T's Matsubara grids in
general (documented as impossible without an interpolating seed scheme — this spec's
own rolling seed is the interpolating scheme, scoped to one finder's internal use).

---

## 6. (2026-07-25) invzp ordered-solver reformulation: strict-order static medium

> **STATUS: implemented as an opt-in scheme; candidate FALSIFIED at Gate 0 on
> 2026-07-27.** The preregistered negative-outcome predicate's clauses **(a), (c), and
> (e)** fired on the real production coupling multiset (clauses (b) and (d) held); per
> the predicate's own binding rule, the run **stopped at diagnosis** — no scheme
> switch, tolerance widening, or regularization was applied, and the `'resummed'`
> production default was never flipped. See `invzp_convg_fix_attmpt_claude.md`
> for the failure narrative and `invzp_strict_medium_gate0_report.md` for
> the measured verdict; neither is restated here. **The theory basis (§A/§B below —
> the `crit(h) = F'(h)` Landau formulation and the moment expansion showing both legs
> truncate identically) is NOT refuted by this failure.** What failed is the specific
> truncated `K0` medium (one-shot, `mu2`-only, `Gref` = Dyson reference) on the real
> 16384-point BZ multiset — not the Landau stationarity argument, and not the identity
> that the PM and ordered legs' exact closures are the same function of the local
> dressed propagator.

**Task/scope of the design.** Un-mask the 1/z ordered spectra panel below `Bc_1z`, keep
the window `Bc_1z < B < Bc_bare` reported as PM, with a defensible 1/z stability
criterion rather than node-convergence luck. Production phase dispatch is already
**PM-mass-first** (`invz_spectra_map.m:288-333`): a successful PM probe with
`crit_pm > 0` owns the window; a *failed* PM probe is **unknown**, not evidence for
order. The defect: an unexplained `phase_1z = 0` region below `Bc_1z`, without letting
solver availability become a phase criterion.

### Decisions of record

1. **Keep Jensen's path-integral self-consistency; reformulate only the static
   medium.** An earlier proposal to replace the full Jensen condition with a new local
   order-parameter fixed point is **superseded**: this removes a beyond-order static
   resummation, it does not invent a new ordered-state stationarity condition.
2. **The pole-free static closure applies to BOTH legs (PM and ordered), not the
   ordered leg alone. Rejected alternative: ordered-leg-only truncation** — a formally
   `O(1/z²)` PM/FM mismatch is dangerous exactly at a continuous boundary (target mass
   exactly zero) and can become the leading residual gap; the shared approximation
   preserves the exact `m → 0` PM/ordered identity, which an ordered-leg-only variant
   would abandon; and that variant does not even have a one-function blast radius —
   `invz_ordered_residual.m:109` independently re-demands the full resummed closure,
   and Gate C1 requires exact equality with the PM medium, so both contracts need
   revision either way. The resummed PM calculation is retained as a **clearly labelled
   legacy comparator**, not the production path — joining two different resummation
   schemes across one transition was judged not defensible merely to preserve old PM
   numbers.
3. **One-shot, no `K0`-denominator feedback:** `K0 = Jbar − mu2·Gref`, evaluated once
   on a `K0`/`λ`-independent `O(1)` reference propagator, in both legs, at `ω_n = 0`
   only. Frozen convention: the **Dyson reference**, `Gref = G0bare/(1+Sigma0)`.
   Documented artifact (not assumed negligible, must be measured — Gate G7): `K(1)`
   (strict) and `K(2)` (exact resummed) differ while `ω_1 → 0` as `T` falls, and both
   feed the same `invz_lambdas` sums — a scheme discontinuity across adjacent grid
   points, not a physical one.
4. **Domain-switched prescriptions are rejected** — the result would depend on where
   the switch occurs. A separately derived local `Phi_1z(m)` is unnecessary unless the
   strict-order construction fails; §A already supplies the Landau-like functional.

**The `D_q` correction is absorbed, not discarded.** `D_q⁻¹` is **not** an `O(1/z²)`
artifact — it carries the leading collective/RPA stability physics. What begins beyond
retained order is the *feedback* of the fully resummed denominator into `K0`. `D_q` is
kept and reported as the susceptibility/stability observable; only its use *as the
definition of the medium* is dropped.

### §A Theory basis (preserved verbatim — not refuted by the Gate-0 failure)

Jensen's spontaneous condition (framework §9.3, J 2.31–2.33) is
`h0(hmf) = ∫₀^hmf r dh'` with `r = G0bare/Gtil0`, and `F = h0 − J0eff·m = 0`. Using
`F(0) = 0` and J 2.31 (`dm/dh = −G0bare`):

```
F'(h) = r(h) + J0eff*G0bare(h) =: crit(h)          [ dimensionless mass ]
F(h)  = integral_0^h crit(h') dh'
```

`crit(0)` is exactly `invz_hmf_ordered.m:138`'s `slope_pred = r0n + J0eff*Gb0` (the
existing `prof.slope0`), generalised to finite `h`.

**Consequence 1 — Jensen's root supplies a one-dimensional Landau-like stationarity
condition.** Setting `Phi_path(m) = ∫₀^m F dm'` gives `dPhi_path/dm = F` and
`dF/dm = crit/chi_path`. Since `chi_path = −G0bare > 0`, the ordered root is a local
minimum **iff `crit(h*) > 0`**. `crit` is dimensionless; the dimensional
inverse-susceptibility/curvature is `crit/chi_path`, so `chi_1z,uni = chi_path/crit`.
The existing bisection (`invz_hmf_ordered.m:216`) is *intended* to select such a
locally stable crossing, but the endpoint derivative must be evaluated, not inferred
from a finite-grid sign change; if multiple stable roots exist, compare `Phi_path`.

Before calling `Phi_path` a thermodynamic free energy, its normalization must be
cross-checked against Jensen's independent relation (J 2.34); **that cross-check is not
available from existing code.** `invz_deltaF_ordered.m` is, by its own docstring, a
"PARTIAL HYBRID DIAGNOSTIC" — explicitly not `dF(m=0)`, meaningful only at a common
cutoff, carrying a documented ~13.7% same-retained-order static-elastic residual —
validation-only, and **cannot be promoted into a "two routes agree" gate**; a genuine
J 2.34 normalization check is separate derivation work, out of scope.

**Consequence 2 — the obstruction localises to the static elastic correction, not to
`Sigma0` alone:** `crit = [1 - J0eff*chi_path(h')] + [r(h') - 1]`. The first bracket is
a direct diagonalization, well defined in the unstable window. At `m = 0` the exact
hybrid identity gives `r = 1 + Sigma0`; at finite ordered moment it does **not** — `r`
also contains the elastic factor `ξ` and depends on `K0` and `λ(1:2)`. The ill-defined
content of the current construction is accordingly `∫[r(h)-1] dh` over the unstable
window; `∫Sigma0 dh` is a useful component diagnostic but not the whole correction
except at onset.

**Binding caution.** The ~0.3% PM boundary shift (`Bc_1z 3.025` vs. bare `3.033`) does
**not** bound either `∫Sigma0(h) dh` or the full `∫[r(h)-1] dh` deep in the ordered
phase — both must be measured (Gate G5), never inferred from the boundary shift.

**Consequence 3 — re-anchoring the integral above the unstable window is impossible.**
`h0(h*)` genuinely requires `r` on all of `[0, h*]`: `H0 = H + J0eff·m` is a state
function, so the integration constant is fixed at `Hmf = 0`. What needs a pole-free
definition is the finite-order static correction `r − 1` throughout that path — not a
replacement for Jensen's condition.

### §B The moment expansion, and why both legs truncate identically (preserved — not refuted)

The ordered static closure (`invz_emt_static_ordered.m`) is
`Gq = Gs/(1 + (J-K0)*Gs)`, `mean_q Gq = Gs`, `K0 = mean_q(J*Gq)/mean_q(Gq)`. Expanding
both closure conditions independently in `Gs` (`J = Jbar + d`, `mu_n = E[d^n]`,
`mu_1 = 0`) — they agree identically:
```
K0 = Jbar - mu2*Gs + mu3*Gs^2 + (2*mu2^2 - mu4)*Gs^3 + (mu5 - 5*mu2*mu3)*Gs^4 + O(Gs^5)
```
Every coefficient is polynomial in the moments: **no q-denominator remains after formal
expansion.** Under high-density power counting `mu2 ~ 1/z`, so the strict `O(1/z)`
medium is `K0 = Jbar - mu2*Gref`, with `Gref` evaluated without feeding the resulting
`K0` back into itself. This is a **formal order statement**, not a proof that the
numerical moment series converges near the critical path.

**The PM leg has the same pole, and the same expansion.** `invz_emt_scalar.m`
eliminates `K` in closed form, but its `A = mean_q J/(D + J*G0)` carries the identical
q-denominator (already flagged in that file as a "vanishing RPA denominator").
Re-expressed in terms of its own local dressed `G`, verified term-by-term to `O(G^4)`
(difference exactly 0):
```
K_PM(G_local) = Jbar - mu2*G + mu3*G^2 + (2*mu2^2 - mu4)*G^3 + (mu5 - 5*mu2*mu3)*G^4
              = K_ordered(G_local)                                   <-- SAME FUNCTION
```
So the two legs' exact closures are the same function of `(G_local; mu_n)`. Truncating
both at the same order **and using the same explicitly defined `Gref` at `m → 0`**
makes the cross-phase static-medium identity exact within the numerical scheme — **but
only if the shared primitive receives an already-constructed `Gref`** from each caller;
sharing the polynomial alone is not enough if the two callers feed it differently
(hence Gate C1 tests the complete caller wiring, not just the shared polynomial).

### Measured multiset (provenance-tagged — **not universal**, valid only for this exact tuple)

Provenance: production default grid `[16 16 16]`, `dpRng 30`, `bruteforce` dipole
backend, no grid-policy fields; `n = 16384` entries (4096 q × 4 branches, no exact
Γ point in this symmetric even-count grid).

| quantity | value |
|---|---|
| `J0eff` (`info.Jcc0`) | `6.42444e-3` meV |
| `Jbar = mean_q J` | `1.20766e-4` |
| `sqrt(mu2)` (rms spread) | `2.3415e-3` |
| `mu2` | `5.48264e-6` |
| `mu3` | `-3.42228e-11` (skewness `-0.0027`) |
| `mu4` | `2.3894 * mu2^2` |
| `J_max`, `J_min` | `5.98514e-3`, `-6.7631e-3` |
| `(J_max - Jbar)/sqrt(mu2)` | `2.5045` sigma |

On this multiset `Jmax < J0eff`; with the physical static sign `Gstat < 0` this gives
`D_q − D_uni = (J(q) − J0eff)·Gstat > 0`, hence `D_q ≥ D_uni` (from those two facts, not
from omitting Γ alone).

**Binding caution — these moments support truncation diagnostics, not a convergence
claim.** The `mu3` contribution is negligible *on this multiset only* — generalising
that is precisely the inference error that produced the original synthetic-`Jnu`
defect this design was meant to repair. Both omitted ratios (`omit_mu3`, `omit_cubic`)
are therefore *always* reported, never assumed small: the first omitted cubic
coefficient `(2*mu2^2 - mu4) = -0.3894*mu2^2` is, relative to the retained `mu2*Gref`
term, approximately 5.2%, 6.0%, 9.3%, and 37.1% of it at `|G| = 155.7, 167.1, 208.2,
416.2`, reaching parity only near `|G| ~ 684`.

**Breakdown thresholds** (`D = 1+Sigma ~ 1`, `|G0|` = bare local static weight) —
monotone in how much beyond-`O(1/z)` denominator feedback is retained, and the direct
justification for decision 3: the uniform instability (real physics, kept as an
observable) precedes the first resummed q-average pole by ~7%.

| `\|G0\|` | event | kind |
|---|---|---|
| **155.7** | `D_uni = 1 + (J0eff-K0)*G <= 0` — physical uniform instability (Γ) | real physics |
| **167.1** | first q-average pole, `D + J_max*G0 = 0` — resummed closure dies | arithmetic |
| **208.2** | `mu2` closure, self-consistent PM/Dyson branch point `\|D+Jbar*G0\| = 2√mu2·\|G0\|` | arithmetic |
| **416.2** | `mu2` closure, one-shot PM/Dyson pole `1+Jbar*G0-mu2*G0^2=0` at `D=1` | arithmetic |

**Binding caution:** the last two numbers are PM ordinary-Dyson estimates only — **not**
universal ordered-elastic thresholds, since `invz_gstat_ordered` has a different local
function; they motivate removing feedback but do not replace an actual ordered-path
domain map (that mapping is Gate G0).

**Rejected alternative, for the record — the self-consistent variant.** Merging
`K0 = Jbar - mu2*G` with `G = G0/(D + K0*G0)` gives a quadratic with a physical branch,
real only while `|D+Jbar*G0| > 2√mu2·|G0|`. It extends the domain over the resummed
closure (measured gain `1.2523`, vs. `2.49` for the one-shot form) but needs different
inputs (`G0`, `D`, a branch choice) and is **not** a normal `static_medium` scheme — it
was placed in a separately named diagnostic comparator, unselectable from the
production resolver, rather than adopted.

### §1 The frozen construction

Jensen's condition is kept **verbatim**; only the `omega_n = 0` medium changes.
```
chi_path = -G0bare = dm/dh > 0
crit(h)  = r(h) + J0eff*G0bare(h)                                [ = F'(h), dimensionless ]
chi_1z,uni(h) = chi_path(h)/crit(h)                              [ where stable ]
F(h)     = integral_0^h crit dh'   = h0(h) - J0eff*m(h)          [ Jensen, unchanged ]
h*       : F(h*) = 0,  crit(h*) > crit_tol      [ Landau minimum ]

Gref     = G0bare/(1 + Sigma0)                                   [ selected Dyson reference ]
K0       = Jbar - mu2 * Gref                                     [ ONE-SHOT, closed form ]
r        = G0bare * (1/Gstat - K0)                               [ arrangement fix, below ]
```
`Gref` is deliberately **not** the physical ordered `Gstat` — the latter still uses
Jensen's hybrid elastic expression with the solved `Sigma`, `lambda`, and the one-shot
`K0`. `D_q = 1 + (J(q)-K0)*Gstat` and `D_uni = 1 + (J0eff-K0)*Gstat` remain built from
`Gstat` and are reported in full as the collective/RPA observables.

**Two-tier verdict — do not collapse these:** **path-node consistency** (the complete
outer `Sigma`/`λ` residual contract passes, every exported quantity `r` needs is
finite; intermediate nodes may have `crit`/`D_uni`/`Dq_min` non-positive — they are
explicitly the analytic continuation through the unstable Landau interval) vs.
**physical endpoint stability** (a consistent root additionally needs `crit > crit_tol`,
`D_uni > D_tol`, `Dq_min > Dq_tol`, tolerances **frozen before the first strict-mode
run** — choosing tolerances after seeing output is tolerance-tuning through the back
door, forbidden).

**Invariant that survives by construction, not tuning:** `r = 1 + Sigma0` at `m = 0`,
for **any** `K0` — `1/Gtil0 = 1/Gstat - K0 = (1+Sigma0)/G0inel0 + K0 - K0` cancels the
`K0` dependence identically, so the pinned `invz_gstat_ordered` identity cannot be
broken by the truncation. The static `K0` also **stops being an inner iteration** —
no fixed point, no branch margin, in the one-shot default.

**Arrangement note — the `Gstat` singularity is removable in the integrand, not merely
in principle.** Because `1/Gtil0 = 1/Gstat - K0`, as `Gstat → ±∞`:
`r → -G0bare*K0` (same limit both sides, continuous through the pole),
`Gtil0 → -1/K0`, `crit → G0bare*(J0eff - K0)` — the divergence **cancels** in the
quantity Jensen's condition integrates, surviving only in `D_uni`/`D_q` (diagnostics).
But the as-written code (`Gtil0 = Gstat/(1-K0*Gstat)`, `r = G0bare/Gtil0`) evaluates
`Inf/(-Inf) = NaN` at `Gstat = Inf` and loses precision approaching it; the
algebraically identical `r = G0bare*(1/Gstat - K0)` is exact to `|Gstat| ~ 1e18` where
the as-written form already goes indeterminate. **This is a reassociation, not a
regulariser** — identical value, no broadening, no tolerance — required for the
removability to actually hold in floating point.

### Components, scheme ids, and error taxonomy (condensed)

New pure leaves in `invz_common/`: `invz_coupling_moments` (moments from `Jnu_flat`);
`invz_static_medium_reference` (constructs `Gref` **and owns its denominator
metadata** — a closure leaf handed only the quotient cannot reconstruct `1+Sigma0` or
its margin); `invz_medium_moment_closure` (`[K0,info]`, both legs call the same map
with an explicit caller-supplied `Gref`, carrying both omitted-order ratios);
`invz_check_static_medium` (sole authority for the scheme id, stamped into both legs so
they cannot disagree); `invz_is_recoverable_solver_error` (below). Existing files
change *in place*, not duplicated, each scheme-gated at `resummed`\|`strict` with the
legacy path bitwise untouched: `invz_emt_scalar.m` (PM medium), `invz_emt_static_ordered.m`
(ordered medium, strict mode removes the damped `K0` iteration entirely),
`invz_gstat_ordered.m` (the `r`/`Gtil0` arrangement fix), `invz_ordered_residual.m`
(block B revised in place, load-bearing residual becomes `\|K0-Kstrict(Gref)\|`),
`invz_hmf_ordered.m`/`invz_ordered_node_solve.m`/`invz_solve_point[_ordered].m` (thread
scheme + provenance), `invz_spectra_map.m` (three-way dispatch, next).

**Scheme ids (frozen):** `'resummed'` (DEFAULT, bit-identical); `'strict_1z_dyson_ref'`
(the candidate that failed Gate 0); `'strict_1z_bare_ref'` (comparator: `Gref` without
`(1+Sigma0)`); unknown id → `invz:staticMedium`. `opts.static_medium` is the sole
public authority (forbidden inside `opts.solve_opts` at driver level); an unstamped
legacy call stays bit-identical.

**Three-way dispatcher** (replaces a fallthrough that let a PM *failure* manufacture an
ordered label): `crit_pm > crit_tol` → **PM**; `crit_pm < -crit_tol` → ordered solve
**phase-eligible**; PM non-convergence/`|crit_pm| <= crit_tol` →
**unknown/boundary-indeterminate** (may run the ordered solver diagnostically but
**cannot** emit a production phase label without a separately validated free-energy
rule). `S.phase_1z_reason` categorical vocabulary: `'ordered'|'pm'|'unstable_endpoint'|
'medium_out_of_domain'|'degenerate_doublet'|'solver_failed'|'pm_probe_unknown'|
'boundary_indeterminate'|'not_attempted_longitudinal'|'bare_not_ordered'|
'bare_escape_hatch'`. `phase_1z` keeps its existing `{0,1,2}` enum (widening would
silently change every saved figure); a stability failure stays `phase_1z = 0` with the
new reason field attached — masked, no longer unexplained. Both the dispatcher and any
un-nesting of the 1/z leg from the auto/overlay branch are **gated on a strict scheme
being active** — under `'resummed'` the historical dispatch is preserved exactly
(applying the rule there would classify the same-pole PM failure as unknown, i.e.
**more** masking).

**Error taxonomy (frozen), one shared classifier for every catch:** wiring/programming
(`invz:staticMedium`, missing `node.Jmom`, `Gref` not supplied, `invz:emtJnu`,
`invz:hzFixed`) → throw loudly, never absorb; physics non-convergence (Picard
`max_iter`, block residual fail) → non-accepted, never throw; out of domain, new class
(reference denominator invalid, `Delta < 1e-4` meV) → distinct status.
`invz_is_recoverable_solver_error(id)` uses an explicit whitelist (initially only
`invz:orderedPhase`, `invz:degenerateDoublet`); unclassified ids rethrow. **Binding
rationale:** at least four broad catches sit on the strict path, and a wiring error
silently downgraded to "node not accepted" at any one is *exactly* the failure mode
that let the original defect hide for a whole stage — adding to the whitelist is a
reviewed contract change, never a convenience fix. Domain contract:
`Delta(Bx, h) >= 1e-4` meV (the constructor's existing threshold, not a new number); a
`domain_policy = 'return'` mode on `invz_twolevel_ordered` records `Delta`/`valid=false`
and returns before any division assuming a resolved doublet; at exactly `B = 0` the
column is labelled `'degenerate_doublet'` and the driver does **not** substitute a
small-field proxy.

### Gate 0 — the frozen negative-outcome predicate (quoted; this is what failed)

> Promotion FAILS if any of: **(a)** any required solved-path node has an invalid
> **reference** denominator status; **(b)** any skipped or invalid node is unaccounted
> for in coverage counters; **(c)** `max(omit_max)` over the solved path exceeds the
> frozen `omit_promote` control bound; **(d)** any **local** `Gstat` denominator
> crossing at which `r` or `crit` is non-finite or discontinuous (a crossing that stays
> finite/continuous — the removable case of §1 — does *not* fail promotion); **(e)**
> any required ordered field fails to return `status='ok'`, a finite nonzero root, and
> a stable endpoint under the frozen `crit`/`D_uni`/`Dq` margins, **or** either
> designated PM control field fails to return a converged, finite, **positive-mass** PM
> state. **On failure the run stops at diagnosis.** Carrying another moment, changing
> `Gref`, or truncating other Matsubara sectors is a **new theory candidate requiring a
> new spec and fresh preregistration** — never an in-run fallback. Regularisation,
> broadening, and tolerance widening remain forbidden.

Clause (e) exists because, without it, (a)–(d) could all pass while every field stayed
masked — a Gate 0 that PASSes without the ordered leg producing even one stable root
cannot serve its purpose. **Preregistration was binding at spec-approval time:** the
spec states it is "design-approved... **not executable**" until `crit_tol`, `D_tol`,
`Dq_tol`, `K_atol`, `K_rtol`, the reference-denominator margin, the `omit_mu3`/
`omit_cubic` bounds, and the full convergence/grid set are frozen with actual numbers
(subsequently frozen; see `docs/invzp_strict_medium_gate0_report.md` for the values
used in the run that failed clauses a/c/e). The spec is explicit Gate 0 "can falsify
the candidate but cannot demonstrate the masking fix" — exactly what happened;
end-to-end validation and any default flip were deferred past it.

### Scope limits the spec itself flagged as unsettled

`invz_static_domain_scan` must accept an explicit `hgrid`, not recreate
`invz_hmf_ordered`'s redensification logic independently; G2 floating-point equality is
exact at `m=0` but the two callers need not reach it bitwise (use `K_atol`/`K_rtol`
unless a fixture proves intentional bitwise identity); whether the `one_field`
un-nesting (lifting the 1/z leg out of the auto/overlay branch) is in scope was
recommended **deferred** — a genuine robustness improvement, but not on the critical
path and unaffected by the Gate-0 outcome; the pre-existing `[nJ,nw]` static-flattening
defect in `invz_emt_static_ordered.m:43` (`Jf = Jnu_flat(:)` flattens *all* frequency
columns into the static q-average) is **not introduced by this design** — recommended
the resolver reject the ordered-retarded combination outright rather than silently
average.

**What does not change, regardless of scheme outcome:** the v5 coupling cache/
`cacheMeta` (moments derived at call time, no schema change); the Ewald backend
(`invz_jq_modes`/`invz_bz_couplings`/`invz_jq_path`, production dipole stays
`bruteforce`); `K(2:end)` in both legs; `pt.crit`'s historical ordinary-Dyson
definition (a legacy diagnostic, not the ordered pole mass).
