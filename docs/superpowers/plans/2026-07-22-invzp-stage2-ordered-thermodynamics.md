# Projected Stage-2: Full Ordered 1/z Thermodynamics — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the bare-MF diagnostic ordered leg of the projected 1/z spectra with Jensen's thermodynamically consistent ordered state (elastic static sector J 2.28–2.29 + nonlinear applied-field/H_MF relation J 2.31–2.33, validated by the free-energy test J 2.34), so the 1/z FM and PM soft-mode branches close continuously at the same critical field `Bc_1z`.

**Vehicle decision (user, 2026-07-22):** the projected path owns the ordered completion. The observable being fixed is the projected spectra map; the tensor `a3d` route cannot supply it (the two routes are never stacked, per `invz_projected/README.html`). The boundary-linearized cosmetic handoff is explicitly rejected — the outcome must be physical. This resolves the open decision recorded in `invzp_QCP_diagnosis.md` Stage 2.

**Revision (2026-07-22, round 1):** amended per `invzp_QCP_stage2_review_codex.md` round 1 (all findings independently re-verified against the source): fixed-field node solves no longer pass `order=true` (P0-1: `invz_single_ion` overrides `hz_fixed` when ordering — verified at `invz_common/invz_single_ion.m:92-97`); a static-sector EMT closure is now its own task before the H_MF integral (P0-2: `invz_emt_scalar` is a direct solve of the ordinary Dyson form, `opts.K0` unused); the free-energy temperature route's sign is corrected to `δF = +T∫_T^∞ δU/T′² dT′` and its tests repaired (P0-3); the H_MF grid is geometric with independent root refinement and a grid-convergence test, and jensen-mode acceptance is gated on root existence rather than `m_tol` (P1-4); the closure regression asserts *pole* closure via the static inverse response, not argmax peaks (P1-5); option threading, the `B(3)` guard, and `ordered_mode` reservation are made explicit (P1-6).

**Revision (2026-07-22, round 2):** amended per the second review (numeric diagnostics re-verified): the static closure now operates on a **boundary-preserving hybrid static propagator** — the J 2.28 *structure* applied to the full-electronuclear static weights (inelastic/elastic split from `invz_chi0z`), with only the ξ resummation factor two-level — because closing the pure two-level propagator approaches a *different* m→0 fixed point than the full-electronuclear PM solver (measured: K0 23% off, crit −0.011 vs −0.062 at 3.0 T), invalidating the onset-coincidence claim (P0-A; §3/§5/§6 rewritten, load-bearing PM-limit identity test added); the route-B Matsubara grid is reconstructed executably (P0-B); the near-zero H_MF bracket extends adaptively below `hmin_frac` when the slope predicts ordering, tested just below a refined `Bc_1z` (P1-C); the final jensen state writes the closed elastic `Gstat` into `pt.G(1)` before diagnostics/sum-rule (P1-D); the regression's remaining argmax subtest is replaced by a minimum-inverse-response tracker (P1-E); static-closure and bare-bracket options are fully threaded (P1-F); the `pt.crit` vs `pt.D_uni` and truncation-estimate contracts are explicit (P2-G).

**Revision (2026-07-22, round 3):** amended per the third review (repo claims re-verified: both solvers now default `max_outer = 200`; `invz_matsubara(T,Ecut)` exists; `invz_chi0z` returns the full 3×3 tensor): the J 2.31 transcription gate uses a CLOSED analytic 2×2 fixture — the repository's `invz_twolevel_ordered` re-embeds the doublet in the full CF manifold, so the two-level identity fails against it by up to 31% (measured; P0-1); the H_MF integrand's `G0bare` is now the derivative of the ACTUAL node path — transverse-MF feedback included via the static-tensor chain rule, with a production finite-field FD gate (P0-2, measured 0.1–0.7% path mismatch with `Jxx0 ≠ 0`); the adaptive lower extension's predicate was logically unreachable and is now driven by an independent h=0 PM predictor, with an explicit `unresolved` status, exposed extension count, and a non-circular near-boundary test refined from the PM `crit = 0` side (P0-3); the final-state `Gstat` is exposed and tested discriminatively — the old ordinary-Dyson value must FAIL the test (P1-4); the H_MF `max_outer` default is aligned to the solvers' 200 (P1-5); the refinement subtest tracks the lowest LOCAL minimum of the inverse response in a bounded basin, with the finite-m non-identity to `D_uni` documented (P1-6); route B's tail is `tail_est_B` with a higher-`Tmax` convergence check (P2-7).

**Revision (2026-07-22, round 4 — review verdict: "nearly implementation-ready", architecture confirmed):** localized fixes per the fourth review: the chain rule now branches explicitly on all three `transverse_mf` modes — `'none'` has NO feedback (`invz_single_ion` forces `hx = hy = 0`), and the FD gate is parameterized over all three (P1-A); the final static refresh in both the H_MF nodes and the jensen final loop KEEPS its newly closed `K0` and exports the self-consistent set (P1-B); route A integrates only `status == 'ok'` profiles with all nodes converged, and a `'node_failed'` status replaces the untruthful `'ok'`-then-reject path (P1-C); the adaptive slow test constructs its coarse lower limit provably above a fine-grid reference root so `n_extend ≥ 1` is guaranteed by correctness, not luck (P2-D); the two stale round-2 instructions are gone (P2-E); and the transverse-feedback bucketing is recorded as a second finite-m modeling choice alongside the two-level ξ (P2-F).

**Revision (2026-07-22, round 5 — FINAL; review verdict: "implementation-ready" after these hardenings, no further physics/architecture revision indicated):** a failed bisection or final-evaluation node now TERMINATES the root solve with `'node_failed'` (never a root from a partial bracket), an unrefined bracket returns `'unresolved'`, and the failure contract is exercised by a dedicated test hook (P1-A); the post-loop static refresh is gated on its own `so.converged` and `so.resid < resid_tol` (a documented `emt_static` contract field) in both the H_MF nodes and the jensen final loop, folded into `pt.converged`, with a residual-only λ/Σ revalidation (`pt.final_resid < tolo`) so the exported tuple is closure-consistent to the stated outer tolerance (P1-B/P2); the Task-4 test regains the exported-tuple lattice-closure assertion (formula consistency alone does not prove the exported `K` closes); the status enumeration documents `'node_failed'` and requires the solver to mask both failure states (P2).

**Architecture:** Four new/modified layers, each gated by numeric identity tests before the next builds on it: (1) a pure two-level function for the ordered elastic static single-site response (J 2.28–2.29); (2) the static-sector EMT closure that makes `K(0)` self-consistent with that elastic propagator; (3) the H_MF integral relation and spontaneous-root solver (J 2.31–2.33), on fixed-field single-ion states; (4) a `'jensen'` mode in `invz_solve_point_ordered` wired into `invz_spectra_map`'s 1/z ordered leg, with a pole-closure regression. The auto/overlay leg, the PM leg, and the longitudinal route are untouched. The free-energy consistency test (J 2.34) is a standalone validation, not a production dependency.

**Tech Stack:** MATLAB R2025a, matlab.unittest function-based tests, existing `invz_projected`/`invz_common` solver stack.

---

## Physics and algebra background (for the independent reviewer)

This section is the review-ready specification. Every formula is transcribed from the repository's framework document `jensen_1z_framework.html` (Jensen's 1/z expansion; J x.y = equation numbers in Jensen's published paper, which are revision-stable). The HTML's own equation numbers are volatile — the working-tree revision renumbers Section 9 by +2 relative to `HEAD` — so citations below give the J-number first and the working-tree HTML number in parentheses.

### 0. The problem being fixed

Stage 1 (`docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md`, executed and merged) gave the 1/z spectrum its own stability-gated phase: above `Bc_1z` (where the PM 1/z mass `crit = 1 + Σ(0) − J(0)χ₀^cc(0)` crosses zero) it uses the strict-PM state; below, only the **bare-MF** ordered state existed, whose moment persists up to the *bare* boundary `Bc_auto > Bc_1z`. Result: at `Bc_1z` the plotted 1/z map jumps from a still-gapped ordered branch to a PM branch rising from zero — the ordered side "does not know" the renormalized critical point. Diagnosis regression 2 ("1/z FM and PM poles meet at `Bc_1z`") is waived until the consistent ordered state exists. This plan builds that state.

### 1. Conventions and code dictionary

| Symbol (framework) | Meaning | Code |
|---|---|---|
| `G(ωn)` | site Green function, **G = −χ** (units meV⁻¹) | `chi = -G` throughout |
| `β` | 1/k_B T in meV⁻¹ | solvers' `beta` (mirrors `invz_single_ion`'s Boltzmann factors) |
| `h ≡ gμ_B H` | fields expressed as energies (meV) | `si.hz`, `hz_fixed`, `J0eff*m` |
| `m, M², n01, Δ, g(0)` | ordered two-level params: diagonal moment, transverse matrix element², population diff `tanh(βΔ/2)`, splitting, static inelastic response | `tl.m, tl.M2, tl.n01, tl.Delta, tl.g0` from `invz_twolevel_ordered(ion,T,Bx,hz,opts)` |
| `g(z) = 2 n01 Δ/(Δ² − z²)` | two-level inelastic response (J/HTML eq 11) | `invz_g(tl, z)` |
| `K(ωn)` | effective-medium (cavity) coupling, J 2.11 (HTML 16) | `med.K` from `invz_emt_scalar`; bosonic grid, `K(1)` = ω=0 slot |
| `λ_p = (1/β)Σ_n K(ωn) g(ωn)^p` | J 2.19 (HTML 21) | `invz_lambdas(K, g, wts, beta, plist)` |
| `Σ(ωn)` ordered | J 2.26–2.27 (HTML 39–40) | `invz_sigma_ordered(tl, lam, K, g, beta)` — **already implemented** |
| `G₀(ωn)` | bare single-site propagator | full electronuclear: `G0 = -real(squeeze(c0(3,3,:)))`, `c0 = invz_chi0z(si, T, 1i*wn, struct('elastic',true))`; static `pt.chi0cc0 = -G0(1)` |
| `⟨Jz⟩` | ordered moment (full manifold) | `si.Jexp(3)`; solver `pt.m0` |
| `J(0)` | uniform ordering coupling (demag/ODD-shifted) | `J0eff` |

**`hz_fixed` semantics (P0-1, binding):** in `invz_single_ion`, `opts.hz_fixed` holds the longitudinal molecular field ONLY when `opts.order` is absent/false — the loop's else-branch comment says "para (0) or held-fixed hz_fixed" (`invz_single_ion.m:97`). With `order=true`, `hz_new = J0z*jz` overrides the imposed field every iteration and the solve converges to the *bare* fixed point regardless of `hz_fixed` (verified numerically: `order+hz_fixed` at `hz_fixed=0.001` returns `hz=0.0237`, identical to `order` alone; `hz_fixed` without `order` holds `hz=0.001`). Therefore: **every H_MF node solve and the final jensen state use `hz_fixed` WITHOUT `order`**; only the separate bare bracketing solve uses `order=true`. Each fixed-field call is followed by a hard assert `si.hz == hz_fixed` (error id `invz:hzFixed`).

### 2. The resummation backbone (already in the code)

Jensen's single-site resummation (J 2.30, worktree HTML eq 30):

```
G(ωn) = G₀(ωn) / [1 + Σ(ωn) + K(ωn) G₀(ωn)]
```

Inserting into the effective-medium form (J 2.12, HTML 17)

```
G(q, ωn) = G(ωn) / [1 + (J(q) − K(ωn)) G(ωn)]
```

makes `K` cancel identically (J 2.20, HTML 31): `G(q,ωn) = G₀/(1 + Σ + J(q)G₀)`. Equivalently the renormalized single-ion propagator (HTML 32) is

```
G̃₀(ωn) ≡ G(ωn)/[1 − K(ωn)G(ωn)] = G₀(ωn)/(1 + Σ(ωn)),
```

which in χ-convention is exactly the code's `chi~0 = chi0/(1+Sigma)`, `chi(q) = chi~0/(1 − J_ν chi~0)` (`invz_chi_realaxis.m:8`). **Nothing in this plan changes the ωn ≠ 0 pipeline.** The EMT solver `invz_emt_scalar` solves this ordinary form in closed form (it is a direct solve — its header states `tol/max_iter/mix/K0` are "accepted for backward compatibility but unused"), enforcing the closure that the q-average of `G(q,ωn)` equals the site `G(ωn)` ("R eq 8" diagnostic in its source). That closed form is exact for the ordinary Dyson structure — and therefore does NOT apply to the elastic static sector below, which breaks that structure (P0-2).

### 3. What is missing, piece 1: the ordered elastic static sector (J 2.28–2.29; worktree HTML 41–42)

When `m ≠ 0`, `G₀` acquires the elastic term and the ω=0 single-site function must be written (J 2.28):

```
G(0) = −M² g(0) / [1 + Σ(0) − M² K(0) g(0)]  −  m² ξ h(0),        h(0) = β(1 − n01²)
```

with the resummed elastic weight (J 2.29):

```
ξ = [1 + tanh( m² n01² β K(0) − M² β λ₁ )]
    / [1 + (4 n01² K(0) g(0) + 2 λ₂ + g(0) λ₁) M² / n01²]
```

Structure to note (reviewer checkpoints):
- The first term's denominator is precisely the J 2.30 resummation at ω=0 for the **inelastic** part, since `G₀,inel(0) = −M²g(0)`: `1 + Σ(0) + K(0)·G₀,inel(0) = 1 + Σ(0) − M²K(0)g(0)`.
- The `tanh` is Jensen's bounded resummation of the elastic weight's leading 1/z term (the linear term × `h(0)` diverges at low T; the resummed form is bounded, and the exponential vanishing of `h(0)` keeps elastic corrections finite at all T). Framework §9.2, transcription verified there against the published paper term-by-term.
- At ω=0 this `G(0)` is inserted **directly into the effective-medium form (HTML 17)** — the K-cancellation does *not* go through for the elastic sector: `G(q,0) = G(0)/[1 + (J(q) − K(0))G(0)]`.
- **Consequence (P0-2): `K(0)` needs its own closure.** The medium coupling is *defined* by the EMT consistency (J 2.10–2.11, HTML 16): the q-average of `G(q,0)` must equal the site `G(0)`, and `K(0) = ⟨J(q)G(q,0)⟩_q / ⟨G(q,0)⟩_q`. For the ordinary form this has the closed-form solution `invz_emt_scalar` implements; for the elastic form it does not, because `G(0)` itself depends on `K(0)` (through both the inelastic denominator and ξ). The static sector is therefore a scalar fixed-point problem solved per outer pass (Task 2), and the resulting `K(0)` feeds back into `λ_p` and `Σ` on the next pass. The ordinary direct solve remains exact — and is the m→0 identity gate — for every `ωn ≠ 0` and for the m=0 limit of `ωn = 0`.

**The boundary-preserving hybrid static propagator (P0-A — the round-2 blocker and its resolution).** Jensen's J 2.28 is written in the two-level model, but the repository's EMT runs on the FULL electronuclear propagator: at m=0 the PM solver's `K(0)` is the direct solve on `G₀,full(0)` from `invz_chi0z`, not on `−M²g(0)`. Closing the pure two-level `G(0)` therefore approaches a *different* m→0 fixed point (measured on the synthetic model at 3.0 T: `G₀(0)` −251.4 full vs −202.6 two-level, `K0` 0.00412 vs 0.00317, and after the Σ loop `crit` −0.0109 vs −0.0620) — the onset would NOT coincide with the PM boundary and the discontinuity would reappear. The plan therefore applies the **J 2.28 structure to the full-electronuclear static weights**, keeping only the ξ resummation factor two-level:

```
Gstat = G₀ᶦⁿᵉˡ_full(0) / [1 + Σ(0) + K(0)·G₀ᶦⁿᵉˡ_full(0)]  +  ξ · G₀ᵉˡ_full(0)
```

where `G₀ᶦⁿᵉˡ_full(0)` is the static full-electronuclear propagator WITHOUT its elastic term (`invz_chi0z(..., 'elastic', false)` at ωn=0), `G₀ᵉˡ_full(0)` is the elastic weight (`elastic:true` minus `elastic:false`, static slot), and ξ is the two-level J 2.29 formula unchanged. Fidelity checks: (a) substituting the two-level weights `G₀ᶦⁿᵉˡ = −M²g(0)`, `G₀ᵉˡ = −m²h(0)` recovers J 2.28 verbatim — the transcription gates in Task 1 run in exactly this parametrization; (b) at m→0 the elastic weight vanishes by the Jz → −Jz symmetry of the hz=0 state (in BOTH sectors), so `Gstat → G₀,full(0)/(1+Σ+K·G₀,full(0))` — the ordinary form on the FULL propagator, whose closure fixed point IS the PM solver's `K(0)` (load-bearing identity test, Task 3 Gate 6b); (c) the bare limit gives `Gstat = G₀,full(0)` and `r = 1`, recovering bare MF with full-manifold weights consistently on both sides of the root equation.
- The framework's erratum note (§9.1) fixes the bare static identity (J 2.31 requirement): `−G₀(0) = d⟨Jz⟩/d h_MF`. **Two distinct realizations, gated separately (round-3 P0-1/P0-2):**
  - *Transcription gate (closed two-level model only):* `−G₀(0) = M²g(0) + m²β(1−n01²) = d(m·n01)/dh_MF` holds EXACTLY for a genuinely closed 2×2 Hamiltonian. It does NOT hold against the repository's `invz_twolevel_ordered`, which re-diagonalizes the full CF Hamiltonian at each `hz` — the doublet subspace itself rotates into higher levels, adding subspace-drift terms (measured mismatch: 31% at `hz = 0.15` meV, 0.6–2% at small `hz`). Gate 1 therefore uses an ANALYTIC 2×2 fixture built inside the test, never the embedded doublet.
  - *Production path gate (full sector):* the H_MF derivation requires the numerator of `r` to be the derivative of the SAME moment curve `m(h)` that enters `F(h)`. The node path re-solves the transverse mean field `hx = Jxx0⟨Jx⟩` at every `hz`, while `invz_chi0z` is the fixed-Hamiltonian response — measured path mismatch 0.1–0.7% with the repo's `Jxx0 ≠ 0` (1e-8 with `Jxx0 = 0`). The production `G0bare` is therefore the PATH derivative, computed analytically from the full static tensor via the transverse-MF chain rule (§4a below) and gated against a centered FD of the node moments at finite `hz`.
- The uniform static inverse response of the ordered lattice is `D_uni = 1 + (J(0) − K(0))·Gstat`: the q=0, ω=0 pole condition. At m→0 the K-cancellation gives `D_uni → [1+Σ(0)+J(0)G₀,full(0)]/[1+Σ(0)+K(0)G₀,full(0)]`, whose zero is the PM `crit = 0` on the SAME full-sector quantities — so `D_uni` (from below) and `crit` (from above) vanish at the same field. `D_uni` is the *pole-based* closure observable the regression in Task 5 uses (P1-5).

### 4. What is missing, piece 2: the applied-field/H_MF relation (J 2.31–2.33; worktree HTML 43–45)

The MF moment is generated by the single-ion problem in the molecular field `h_MF`; its differential response is `−G₀(0; h_MF)` (J 2.31). The exact response to the *applied* field is the bulk q=0 susceptibility. Equating the two expressions for `d⟨Jz⟩` (J 2.32) and using the effective-medium structure, the equation separates into the boxed applied-field integral (J 2.33):

```
h₀ = ∫₀^{h_MF} [ G₀(0; h′) / G̃₀(0; h′) ] dh′ ,       h₀ ≡ h_applied + J(0)·⟨Jz⟩
```

where at ω=0 the ordered `G̃₀(0) = G(0)/[1 − K(0)G(0)]` uses the elastic `G(0)` of piece 1 with the closed `K(0)` of Task 2, and `G₀(0; h′)` is the bare static function at the state imposed by molecular field `h′` (fixed-field solve, §1). For the **spontaneous transverse-field problem** (`h_applied = 0` along the ordering axis c/z), the self-consistency becomes a scalar root problem in `h_MF`:

```
F(h_MF) ≡ h₀(h_MF) − J(0)·⟨Jz⟩(h_MF) = 0,   with the nonzero root taken when it exists.
```

### 4a. The path-consistent static derivative (round-3 P0-2)

The bare numerator `G₀(0;h′)` of `r` must be `−dm/dh′` along the path the nodes actually trace — with the transverse mean field re-solved at every `h′`. From the full static susceptibility tensor `X` (χ-convention, `X(a,b) = −G₀(a,b;0)` from `invz_chi0z(..., 'elastic', true)` at ωn=0), the chain rule through the transverse channel(s) `t` (x under `legacy_x`; {x,y} under `vector_ab`) gives

```
dm/dh′ = X_zz + X_zt · Jxx0 (I − Jxx0·X_tt)⁻¹ · X_tz        ⇒   G0bare = −dm/dh′
```

(1×1 or 2×2 transverse block). The sector split then is: `G0inel` stays the FIXED-Hamiltonian inelastic static (`invz_chi0z` `elastic:false`, ωn=0 slot) — consistent with the EMT's dynamic propagator, which also carries no transverse feedback (the feedback is an adiabatic, strictly-static effect); `G0el = G0bare − G0inel` is the elastic weight PLUS the small transverse-feedback correction, both of which vanish at the boundary (the cross response `X_zt → 0` by the Jz → −Jz symmetry at hz = 0, and the elastic weight vanishes quadratically — both measured). Consequences: the bare limit `r ≡ 1` now holds exactly ON THE PATH (so the bare root reproduces the bare solver *including* its transverse feedback), and the m→0 limit is unchanged. A production gate compares the chain-rule `G0bare` against a centered FD of the node moments at a finite-moment `hz` (two extra fixed-field single-ion solves, test-only).

Define the integrand ratio `r(h′) ≡ G₀(0;h′)/G̃₀(0;h′)` with this path-consistent numerator. Its two limits are the plan's load-bearing checks:

- **Bare limit** (`Σ = K = λ = 0`): `ξ = 1`, `G(0) = G₀(0)`, `G̃₀(0) = G₀(0)` ⇒ `r ≡ 1` ⇒ `h₀ = h_MF` ⇒ `h_MF = J(0)⟨Jz⟩` — the current bare-MF code exactly.
- **m → 0 limit**: the elastic term vanishes as m², the K-cancellation goes through, `G̃₀(0) → G₀,inel(0)/(1+Σ(0))` ⇒ `r → 1 + Σ(0)`.

### 5. The onset-coincidence result (why the discontinuity closes *by construction*)

Linearize the root condition at m → 0. The single-ion moment responds as `⟨Jz⟩ = χ₀^cc(0)·h_MF` (full electronuclear static, exactly what `invz_single_ion` + `invz_chi0z` deliver), and `h₀ = (1+Σ(0))·h_MF` by the m→0 limit of `r`. A nonzero root of `F` exists iff

```
J(0)·χ₀^cc(0) > 1 + Σ(0)    ⇔    crit = 1 + Σ(0) − J(0)χ₀^cc(0) < 0.
```

**Round-2 correction (P0-A):** this conclusion identifies the onset with the code's PM boundary ONLY if the `Σ(0)` reached by the ordered loop's m→0 limit is the PM solver's `Σ(0)` — and `Σ(0)` depends on `K(0)` indirectly through every `λ_p`. The K-cancellation inside `r` removes the *direct* static-propagator factor but NOT this indirect dependence. With the pure two-level static closure the m→0 `K(0)` differs from the PM solver's by ~23% (measured, §3) and the masses disagree. The **hybrid static propagator of §3 restores the identity**: at m→0 its closure degenerates to the ordinary direct solve on the full-electronuclear `G₀,full(0)` — the PM solver's own fixed point — so `K(0)`, `λ_p`, `Σ(0)`, and hence the onset condition all coincide with the PM leg's `crit = 0`. **The ordered moment then onsets exactly at `Bc_1z`**, the ordered state approaches the PM state continuously (`m → 0`, `h_MF → 0`), and the 1/z spectra join continuously. This also fixes the *analytic small-field slope*: `F(h)/h → crit` as `h → 0⁺` — the root-existence detector Task 3 uses (P1-4). This is a theorem about the hybrid construction, verified numerically by Gate 6 (onset bracketing) and Gate 6b (the h→0 loop reproduces the ACTUAL `invz_solve_point` fixed point: `K(1)`, `Sigma0`, `crit` — the load-bearing PM-limit identity test the round-2 review demanded); it is not assumed.

### 6. The hybrid sector mapping (a deliberate, documented approximation — revised in round 2)

The projected code has an established hybrid: `Σ, α_m, λ_p, g` live in the **two-level** (electronic doublet) sector, while the propagator `G₀` fed to the EMT is the **full electronuclear** `invz_chi0z` (hyperfine manifold included, renormalized "in one piece"). Round 1 of this plan put the static closure entirely in the two-level sector; the round-2 review showed numerically that this breaks the m→0 identity with the PM solver (§3, §5). The corrected mapping:

- the **static propagator weights** (`G₀ᶦⁿᵉˡ_full(0)`, `G₀ᵉˡ_full(0)`) and therefore the static closure of `K(0)`, the integrand `r(h′)`, and the pole observable `D_uni` live in the **full electronuclear** sector — matching what the EMT actually propagates;
- only the **ξ resummation factor** (J 2.29) and the dynamic machinery (`Σ(ωn)`, `α_m`, `λ_p`, `g`) remain two-level — where Jensen's ordered algebra is defined;
- the moment `⟨Jz⟩(h_MF)` and the criticality mass use the full electronuclear quantities, as the code already does.

This hybrid is *exact at the boundary by construction* (§3 fidelity check b, §5, Gate 6b) and *exact in the bare limit with full-manifold weights* (§3 check c). At finite m there are TWO residual modeling choices (round-4 P2-F), both exact in the two load-bearing limits, both ACCEPTED as documented scope limitations — the closed-model J 2.34 test of §7b does NOT validate them (third §7 review); they are probed only indirectly by the boundary gates and the closure regression:

1. **ξ stays two-level**: the resummation factor's fourth-order semi-invariant content (J 2.29's λ-groupings) exists only in the two-level algebra; it is applied multiplicatively to the full elastic weight.
2. **Transverse-feedback bucketing**: the chain rule (§4a) rigorously establishes the path derivative `G0bare`, but placing the feedback correction inside `G0el` — where it is multiplied by ξ and enters the local static EMT closure — is an additional finite-m choice, not Jensen's original elastic spectral weight. It is exact at m→0 (the cross response `X_zt` vanishes) and in the bare limit (ξ = 1 and the closure is trivial); between them it is a modeling decision.

An independent reviewer should confirm they accept both scopings; the alternatives (an electronuclear ξ, or a first-principles ordered treatment of the transverse channel) are beyond the framework document.

### 7. Free-energy consistency test (J 2.34; worktree HTML 46), the global validation

Because `⟨Jz⟩ = ⟨Jz⟩₀`, comparing the 1/z and pure-MF theories at the same single-ion state gives `δh = h₀ − h_MF`, and integrating the Legendre-transformed relation from saturation down to zero moment (framework §9.4):

```
δF(m=0) = − ∫₀^{M₀} δh  d⟨Jz⟩          (per site, meV; route A — probes the FIELD dependence of Σ)
```

The same quantity follows independently by temperature-integrating the internal-energy correction δU (J 2.22; worktree HTML 37):

```
δU = ½ { α n01 Δ/(1+α)  −  M² λ₁  +  (1/β) Σ_n K(ωn)[G(ωn) − G₀(ωn)] }     (per site)
```

**Sign derivation (P0-3, corrected):** from `d(δF/T)/dT = −δU/T²` and `δF/T → 0` as `T → ∞`, integrating from `T` to `∞`:

```
δF(T) = + T ∫_T^∞ (δU/T′²) dT′          (route B — probes the FREQUENCY structure)
```

(The framework's compressed phrase "δF = −T∫(δU/T²)dT from T=∞" integrates *downward* from ∞, which is the same statement; an earlier draft of this plan transcribed it with the wrong sign.) Agreement of the two routes is the framework's own "stringent global check on a numerical implementation." **Execution outcome:** the check as originally scoped (production hybrid, full-moment saturation) FAILED for domain reasons diagnosed in §7a/§7b — the gate now runs in the closed 2×2 reference model (§7b, Task 6b), which exercises the same production functions with an unbroken one-moment conjugacy chain; the production hybrid's absolute δF is a documented domain limitation, not a gate.

### 7a. Route A resolution — SUPERSEDED DRAFT (kept for the record; see §7b)

The first route-A amendment below (two-level measure `m·n01` + quench-detected endpoint on the PRODUCTION profile) was **rejected by independent review** (`invzp_7a_plan_review_Codex.md`, 2026-07-22) on two load-bearing grounds, both empirically verified: (1) substituting `dm2` only in the final Legendre integral mixes thermodynamic coordinates — the production `h0` is built from the FULL-path derivative and `F` keeps the full moment, and the measured `d(m2)/dM` ratio swings 0.84–1.24; (2) the production `tl` comes from `invz_twolevel_ordered`, which re-selects two states of the full CF solve at every `hz`, so `m2` keeps acquiring CF admixture and NO quench is observed on any tested grid (`|r−1|` = 0.29/0.18/0.16 at `hmax_fac` 4/8/16 against a 0.0067 threshold). The superseding design is §7b. Original draft follows unchanged for provenance:

#### (superseded) Route A corrected: the two-level measure and the quench-detected endpoint (execution amendment 4)

**Measured failure of the original route A (Task-6 execution, 2026-07-22).** The as-planned route A — `δF = −∫ δh dm` over the FULL-manifold moment `m = si.Jexp(3)` toward saturation `M0 = ion.J = 8` — does not converge: `dFA` = −0.0094 at `hmax_fac = 4`, −0.0032 at 8 (197% apart), and +0.49 at 300 (sign flip); `tail_est` exceeded its 5% bound by 172×. Route B is internally clean at the same point (`dFB` = −0.0376, grid-stable to 0.69%, tail 1.7e-4).

**Why the construction, not the algebra, is at fault.** Jensen's §9.4 anchor — "infinite field, fluctuations quenched, δΦ = 0" — is a statement inside the TWO-LEVEL model: there, doublet saturation and the quenching of every fluctuation weight (`h(0) = β(1−n01²) → 0`, `g(0) = 2n01/Δ → 0`, hence λ, Σ, α_m → 0 and `r → 1`) happen together at `h_MF` of a few splittings. The hybrid's FULL moment instead approaches `J = 8` only at crystal-field-mixing field scales, where (a) the doublet parametrization underlying ξ, Σ, and the closure loses meaning, and (b) the static-closure fixed point is far outside the regime every Task 1–5 gate validated (the fac-300 sign flip is consistent with a spurious closure branch there — a domain boundary, not a defect in the boundary-anchored algebra). Integrating the full measure to `J = 8` therefore drags the validation through territory the model does not cover.

**Corrected construction.** Evaluate route A in the model's own terms:

```
δF_A = − ∫₀^{m₂(h_q)} δh d m₂ ,     m₂(h′) ≡ tl.m · tl.n01   (Jensen's ⟨Jz⟩₀, per node)
```

with `δh(h) = h₀(h) − h` exactly as before, and the endpoint `h_q` QUENCH-DETECTED rather than assumed: the first grid point beyond which the integrand has died, `|r(h) − 1| < q_frac · |r(0⁺) − 1|` (default `q_frac = 0.01`, i.e. the fluctuation corrections have decayed to 1% of their h→0 value `Σ(0)`), sustained to the end of the profile. Contracts: (i) if no quench is observed before `hmax`, the routine returns status `'no_quench'` and the test FAILS with the measured `|r−1|` profile — never a silent tail guess; (ii) `tail_est_A = |δh(h_q)| · (m₂_plateau − m₂(h_q))` with `m₂_plateau = max(m₂)` over the profile — genuinely small at quench because `m₂` saturates WITHIN the doublet regime; (iii) the `hmax_fac` 4-vs-8 convergence assert (5%) is retained and is now expected to pass because beyond `h_q` the two-level measure has plateaued, so extending the grid adds ~nothing. The full-moment `m` stays in `F(h)` and everywhere else — ONLY the δF measure changes, because the §9.4 Legendre bookkeeping is the two-level model's identity.

**Honest scope note (reviewer must weigh):** route B's δU (J 2.22) is likewise the model's internal energy evaluated in the hybrid (two-level α, n01, Δ, M², λ₁; the `K(G−G₀)` term with the hybrid's full-manifold propagator). Both routes thus target the MODEL's free-energy correction as realized by the hybrid; residual sector-mixing differences between the two routes are part of the §6-documented finite-m approximations, and the 10% two-route gate is precisely their global probe. The corrected route A requires one additive change to `invz_hmf_ordered` (expose `prof.m2` per node) and a rewrite of `invz_deltaF_ordered` + its tests.

**For the independent reviewer of this amendment:** (1) is the two-level measure `m·n01` the correct work-conjugate for §9.4's `d⟨Jx⟩`, given the framework's own `⟨Jx⟩ = ⟨Jx⟩₀ = m·n01` (§9.1)? (2) is the quench criterion (`|r−1|` decayed to 1% of its boundary value, sustained) the right operationalization of "fluctuations quenched, δΦ = 0"? (3) does capping the domain at `h_q` legitimately excise the CF-mixing regime, or does the excised tail carry model-physics the check should see? (4) are the `'no_quench'`/tail/convergence contracts sufficient to keep a false pass impossible?

### 7b. Route A resolution — ADOPTED: closed-model J 2.34 validation + production scope limitation (per the independent review's recommended option)

**Scope decision.** The J 2.34 two-route identity is validated in a **genuinely closed 2×2 model** in which every object — the moment, its field derivative, the field relation, both free-energy routes, and the saturation anchor — belongs to ONE theory, exercising the SHARED LOW-LEVEL KERNELS (`invz_g`, `invz_matsubara`, `invz_lambdas`, `invz_sigma`, `invz_sigma_ordered`, `invz_gstat_ordered`, `invz_emt_static_ordered`, `invz_emt_scalar`) through an independent local harness — it does NOT execute `invz_hmf_ordered` or the production hybrid end to end (second §7 review, P2 wording correction). The production hybrid's saturation-normalized absolute `δF(m=0)` is declared **outside its validated domain**: its acceptance rests on the boundary-anchored gates and the pole-closure regression (which is where its physics claims live), the original failed route-A numbers are PRESERVED in the execution record as evidence delimiting the hybrid's domain, and `invz_deltaF_ordered` is retained only as an explicitly labeled **partial, cutoff-dependent hybrid diagnostic** (`dF_partial`) that is never called `δF(m=0)` and never gates anything.

**The closed validation model** (frozen 2×2 Hilbert space — never re-selected from a CF solve):

```
H(h) = diag(+Δ0/2, −Δ0/2) − h·Jz2,      Jz2 = [0 M0op; M0op 0]   (FIXED operator)
```

**Lattice fixture identity (second §7 review, P0-1 — load-bearing):** the coupling spectrum MUST satisfy Jensen's no-self-site identity `⟨J(q)⟩_q = J(ii) = 0`. The earlier synthetic `linspace(−2e-3, 6e-3, 24)` has mean 2e-3 and is NOT a valid closed Jensen fixture (measured consequence: the two routes come out with opposite signs, +9.6e-4 vs −9.9e-3). The fixture is `Jnu = linspace(−4e-3, 4e-3, 24).'` with a load-bearing assertion `abs(mean(Jnu)) < 1e-15`, commented as enforcing `J(ii) = 0` — not as numerical centering. (The stage-1/2 *dispatch* tests may keep the old spectrum: they test labels and consistency identities, not the thermodynamic identity.)

**OPEN closed-model discrepancy (second §7 review, P0-2 — must be explained, never tolerated away):** with the corrected fixture the reviewer measures route A = −5.23e-4 vs route B = −6.26e-4 (16.5%; 12.9% after 4× denser T-spacing) — stable under `hmax`/`nH` doubling (0.06/0.08%), coupling-scale halving (~12%), and `Ecut` 40→80. Route A is internally converged; the mismatch is therefore either a residual loop/transcription inconsistency or a structural consequence of comparing the self-consistently Dyson-resummed implementation with the order retained in J 2.22's δU. **Task 6b therefore mandates a localization investigation (δU via J 2.21 vs J 2.22; same-final-tuple audit; T-density convergence; coupling-scale exponent scan) and ends in a BLOCKED escalation with the evidence — the 10% gate is finalized only after the discrepancy is explained (fixed if a transcription error; derived and re-gated if a controlled resummation-order limitation). The rough scale-independence under coupling halving currently points toward a transcription/order inconsistency rather than a higher-order-in-J effect.**

Choosing the PURELY OFF-DIAGONAL `Jz2` makes the h=0 state a genuine paramagnet (`m(0) = 0` — the eigenstates carry no diagonal moment), so the `δF(m=0)` anchor exists in-model, and the saturation moment is ANALYTIC: `M_sat = M0op` exactly (the eigenvalue of `Jz2`), with `n01 → 1` and all fluctuation weights (`g0 = 2n01/Δ ~ 1/h`, `h(0) = β(1−n01²) → 0` exponentially) genuinely dying at large `h` — the §9.4 quenching premise HOLDS here by construction, unlike in the production hybrid. Per node `h`: diagonalize the 2×2, read `m, M², n01, Δ` (Task 1's `twolevel_fix` pattern with `a = 0`), and run the ordered Σ↔EMT loop with TWO-LEVEL static weights (`G0inel = −M²g0`, `G0el = −m²h0`) — the exact parametrization in which J 2.28 is verbatim (§3 check a). One moment `M(h) = m·n01` serves J 2.31 (`−G0 = dM/dh`, exact — Gate 1), J 2.33 (`h0 = ∫ r dh`), the work term, and J 2.34's measure — the conjugacy chain is unbroken.

- **Route A (closed):** `δF_A = −∫₀^{h_end} δh · (dM/dh) dh` on a geometric grid to `h_end` where saturation is VERIFIED, not assumed: require `|r−1| < max(q_abs, q_rel·|r(0⁺)−1|)` (defaults `q_abs = 1e-4`, `q_rel = 0.01`) sustained over the last ≥ 5 nodes AND the primitive corrections `|Σ(0)|` AND `|λ₁·g0|` small over the SAME tail window (both exported per node — second-review P1) AND `M_sat − M(h_end) < 0.01·M_sat`. `tail_est_A = |δh(h_end)|·(M_sat − M(h_end))` uses the ANALYTIC `M_sat` and is an ESTIMATE, not a proved bound (second-review P2 — the `h_end`-doubling result carries the convergence claim). Convergence: `h_end`-doubling and grid-doubling each move `δF_A` < 2%.
- **Route B (closed):** `δU(T′)` per J 2.22 at `h = 0` (PM, m = 0, elastic weight zero): `α` and `λ₁` from the PM machinery (`invz_sigma`), `K` from `invz_emt_scalar` on `G0(iωn) = −M²g(iωn)`, `G = G0/(1+Σ+K·G0)`; `δF_B = +T∫_T^∞ δU/T′² dT′` on a geometric T-grid with BOTH convergence checks (second-review P1): `Tmax` extension (5%) AND temperature-DENSITY (ratio 1.35 vs √1.35 over the same endpoints, 5%; the DENSER value enters the gate).
- **Gate (PENDING the P0-2 investigation):** `|δF_A − δF_B|` at a tolerance finalized only after the measured 13–16% closed-model discrepancy is localized and explained (transcription fix → 10%; derived resummation-order limitation → a gate derived from that analysis). Until then Task 6b ends in the mandated BLOCKED escalation with the localization evidence. Cost: trivial (2×2 diagonalizations) — FAST, not slow-gated.

**What this validates and what it does not (corrected per the Task-6b blocker review — no overclaim):** the closed-model harness validates EXACT LIMITS (zero-coupling, m→0 boundary, bare limit), CONVERGENCE (field grid/ceiling, temperature density/ceiling, final tuples), and LEADING-ORDER SCALING (`|δF_A|, |δF_B|, |δF_A−δF_B| ~ s²`). It does **NOT** validate cross-route thermodynamic closure: the measured ~13.65% route residual is a **same-retained-order static-elastic approximation residual** (its absolute size is O(s²) — the SAME leading coupling order as both corrections; it is NOT a "higher-order" effect and must never be labeled one). Nor does it validate the hybrid's sector-mixing choices at finite m.

**The published-paper evidence (Task-6b blocker review, checked against J. Jensen, PRB 49, 11833 (1994)):** Jensen's own Eq. 2.34 check used the FULL ordered machinery — Eq. 2.33 solved fully self-consistently, no paramagnetic shortcut — and was satisfied to **2–3% at low temperature** on the physical HoF3 lattice; Jensen himself flags that at higher temperature the elastic contribution "requires more careful treatment." That is direct published evidence for a static-elastic same-order consistency limitation, AND it means: (a) a nonzero route residual exists even in Jensen's own implementation — Eq. 2.34 is a stringent numerical check, not a machine-precision identity; (b) our 13.65% cannot be excused as testing something Jensen didn't; (c) the 2–3% is regime/lattice-specific, not a universal bound, and 13.65% justifies no universal 15% allowance. The reviewer's low-temperature cross-check on OUR synthetic fixture (T = 0.10 K → 14.48%) shows the synthetic residual does not approach Jensen's 2–3% by cooling — it is fixture-dependent approximation behavior or an unexposed implementation issue, presently undecided.

**Three-outcome resolution structure (supersedes the earlier binary):**
1. a DEMONSTRATED implementation defect is found and fixed → rerun the equality gate;
2. a faithful reproduction of Jensen's published low-temperature HoF3 Eq.-2.34 check (full ordered machinery, the paper's physical lattice/parameters, quadrature error budgeted separately) reproduces the 2–3% → the synthetic 13.65% is recorded as fixture-dependent static-approximation behavior and the implementation gains its scientific clearance;
3. absent either → J 2.34 remains a DIAGNOSTIC and the thermodynamic closure remains a documented, scientifically unvalidated approximation limitation — even while engineering work proceeds.
The HoF3 discriminator (outcome 2's experiment) — or, alternatively, restoring the omitted frequency-dependent elastic terms / deriving both routes from one common approximate functional — is a substantial follow-up physics task OUTSIDE this plan's scope.

**The engineering gate that IS derivable now (replaces the invalid percentage-equality gate):** keep every exact-limit and convergence gate, add the ORDER-CONSISTENCY gate — under `Jnu → s·Jnu` require the fitted exponents `p_A, p_B, p_Δ ∈ [1.9, 2.1]` and `|p_ε| ≤ 0.15` (engineering TEST TOLERANCES encoding the observed leading-order structure, explicitly NOT analytic error bounds) — and pin the route values themselves as a NON-GATING **approximation regression** (`dFA`, `dFB` at T = 0.31 K and the T = 0.10 K cross-check point pinned to their recorded values within 2%, named an approximation fingerprint against code drift, never a two-route validation). Neither 10% nor 15% equality may be committed as a physics gate. The 2.68% J 2.21/J 2.22 residual is a DIFFERENT-order effect (vanishes under the coupling scan) and is no part of any 13.65% budget. The probe ablations are counterfactual trajectory changes — sensitivities, NOT additive shares or error bounds.

**Reviewer answers incorporated:** the quench criterion is absolute+relative with primitive-correction checks and a minimum tail span, singular-free in the bare limit (where `r ≡ 1` makes route A trivially zero — asserted, not `no_quench`); the tail uses the analytic asymptote; the production `dF_partial` diagnostic reports its cutoff dependence and is non-gating; the original failure numbers stay in the record.

### 8. What deliberately does NOT change

- The **auto/overlay leg** (`S.phase`, `S.chirpa`, `S.Bc_auto`): stays bare-MF — it is the Σ=0/RPA proxy and the input to the future independent-RPA-evaluator follow-up (diagnosis regression 4, still open, unaffected by this plan). `ordered_mode` becomes a RESERVED `solve_opts` field (P1-6) so a caller cannot flip the auto leg into jensen mode.
- The **PM leg** and its `crit` definition.
- The **longitudinal (field-induced) route**: rounded crossover, no sharp `Bc`; keeps the bare mode. `'jensen'` mode errors on `forced_moment` AND on `|B(3)| > bz_tol` directly (P1-6).
- The ωn ≠ 0 spectral pipeline (`invz_chi_realaxis`), demag handling, ODD machinery, and everything in `invz_tensor/`.

### 9. Escalation policy (binding on every implementer)

The identity gates in Tasks 1–3 pin all sign conventions numerically. If any gate fails: STOP, report BLOCKED with the numeric mismatch. **Never** flip a sign, insert an `abs()`, or loosen a tolerance to make a gate pass — a failing gate means the transcription or the mapping is wrong, and that is a controller/human decision.

---

## Global Constraints

- Run every MATLAB command from the repository root; the path contains spaces — keep the binary path quoted: `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "..."`.
- Fast-suite gate after each task: `runtests('invz_projected/tests')` reports **0 failed**. Baseline entering this plan: **145 passed / 0 failed / 19 incomplete**. Expected counts per task are stated in each task.
- TDD throughout: failing test first, then implementation; identity-gate failures are BLOCKED escalations (§9 above), not iteration fodder.
- Error policy: absorb only `invz:*` identifiers; anything else rethrows.
- Fixed-field single-ion calls NEVER combine `order=true` with `hz_fixed` (§1); every fixed-field call asserts `si.hz == hz_fixed` (`invz:hzFixed`).
- `'jensen'` mode is TRANSVERSE-ONLY: it raises `invz:orderedMode` on `opts.forced_moment` AND on `|Bx(3)| > bz_tol` (default 1e-9 T).
- Back-compat is binding: `opts.ordered_mode` defaults to `'bare'` in `invz_solve_point_ordered` and the bare path must remain byte-identical (the auto leg, `Bc_auto`, and tensor-parity work depend on it). The map's 1/z leg opts in to `'jensen'` explicitly; `invz_check_solve_opts` reserves `ordered_mode` as driver-owned (P1-6).
- Option threading (P1-6): the H_MF machinery must run under the SAME numerical settings as the final solve — `Ecut`, `mf_maxit`, `mf_mix`, `mix_outer`, `tol_outer`, `max_outer`, and any EMT options are forwarded from the caller's opts into every node solve; the jensen branch passes the caller's full opts through (with `J0eff` replaced by the ODD-shifted value and mode fields stripped), never a hand-picked subset.
- Stage-1 knobs still apply near the boundary: `solve_opts.max_outer` / `mf_maxit` (critical slowing).
- Stage only each task's files — the worktree may carry the user's unrelated uncommitted edits; never `git add -A`. (`invz_projected/README.html` may be dirty: use the established stash choreography — `git stash push -m "<task>-user-readme" -- <file>` before editing, pop after committing, never touching pre-existing user stashes.)
- New physics code lives beside its sector: two-level algebra in `invz_common/`, projected-path solvers/drivers in `invz_projected/`.

---

### Task 0: Stage-1 closeout — stale-file cleanup and execution records

**Files:**
- Delete: `.superpowers/sdd/task-*-brief.md`, `.superpowers/sdd/task-*-report.md`, `.superpowers/sdd/review-*.diff`, `.superpowers/sdd/final-fix-report.md` (git-ignored scratch; KEEP `progress.md`, the ledger)
- Modify: `invzp_QCP_diagnosis.md` (record the vehicle decision)
- Commit (new, untracked): `docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md`, `invzp_QCP_diagnosis.md`

**Interfaces:** none — housekeeping. Repo convention: executed plans are committed as execution records (precedent: commit `1a27d54`).

- [ ] **Step 1: Remove the Stage-1 SDD scratch files**

```bash
cd "<repo root>" && rm .superpowers/sdd/task-*-brief.md .superpowers/sdd/task-*-report.md .superpowers/sdd/review-*.diff .superpowers/sdd/final-fix-report.md && ls .superpowers/sdd/
```
Expected: only `progress.md` remains.

- [ ] **Step 2: Record the vehicle decision in the diagnosis**

In `invzp_QCP_diagnosis.md`, Stage 2 section, replace the sentence beginning `Decision to make before starting: which route owns the ordered completion.` and its paragraph's final sentence `The projected README forbids stacking the two routes, so choose one vehicle explicitly rather than duplicating the ordered-side work.` — keep the middle of the paragraph intact — so the paragraph now opens:

```
Decision (2026-07-22, user): the PROJECTED path owns the ordered completion -- the observable being fixed is the projected spectra map, which the tensor route cannot supply (never stack the two routes). Implementation plan: docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md. The boundary-linearized cosmetic handoff is rejected. For context, the tensor a3d machinery remains NOT a complete vehicle by itself:
```
and closes:
```
The projected README forbids stacking the two routes; the tensor branch keeps its own scope.
```

- [ ] **Step 3: Commit the execution records**

```bash
git add docs/superpowers/plans/2026-07-21-invzp-qcp-stage1-split-overlays.md invzp_QCP_diagnosis.md
git commit -m "docs(plans): commit the stage-1 split plan + QCP diagnosis as execution records; record the stage-2 vehicle decision (projected path) -- repo convention"
```

---

### Task 1: `invz_gstat_ordered` — the elastic static sector (J 2.28–2.29)

**Files:**
- Create: `invz_common/invz_gstat_ordered.m`
- Test: `invz_projected/tests/test_invz_gstat_ordered.m` (new file)

**Interfaces:**
- Consumes: `tl` from `invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0',..,'transverse_mf',..))` (fields `m, M2, n01, Delta, g0` — used ONLY for ξ); `lam = invz_lambdas(K, g, wts, beta, [1 2])` (λ₁, λ₂ suffice — ξ does not use λ₃); scalars `K0 = K(1)`, `Sigma0 = Sigma(1)`, `beta`; static propagator weights `G0inel0`, `G0el0` (G-convention scalars, meV⁻¹) supplied by the CALLER — the hybrid design of §3: the production loops pass full-electronuclear weights, the transcription gates pass the two-level weights `−M2*g0` / `−m²h0` to pin J 2.28 verbatim.
- Produces (Tasks 2–3 rely on): `[Gstat, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0)` with `out.xi`, `out.h0`, `out.G0bare` (= `G0inel0 + G0el0`), `out.Gtil0`, `out.r` (all scalars).

- [ ] **Step 1: Write the failing identity tests**

Create `invz_projected/tests/test_invz_gstat_ordered.m` (note the helper builds `tl` directly from `invz_twolevel_ordered` — no `invz_single_ion` call, and in particular no `order`+`hz_fixed` combination, per §1):

```matlab
function tests = test_invz_gstat_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function tl = ordered_tl(hz)
% Ordered two-level params at an imposed molecular field.
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
end

function [gi, ge] = twolevel_weights(tl, beta)
% Two-level static weights: the parametrization in which J 2.28 is EXACT (SS3 check a).
gi = -tl.M2*tl.g0;
ge = -tl.m^2 * beta*(1 - tl.n01^2);
end

function tl = twolevel_fix(h, D0, a, M0, beta)
% CLOSED analytic 2x2 fixture (round-3 P0-1): H = diag([D0/2, -D0/2]) - h*Jz2 with the
% FIXED traceless operator Jz2 = [a M0; M0 -a]. Unlike invz_twolevel_ordered -- whose
% doublet re-embeds in the full CF manifold, adding subspace-drift terms that break the
% J 2.31 identity by up to 31% (review-measured) -- this model is closed, so the
% identity is EXACT. Diagonals of V'*Jz2*V stay +/-m (traceless), so <Jz> = m*n01.
Jz2 = [a, M0; M0, -a];
H = [D0/2, 0; 0, -D0/2] - h*Jz2;
[V, E] = eig((H+H')/2, 'vector');  [E, ix] = sort(real(E));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
Mz = V'*Jz2*V;
tl = struct('m', real(Mz(1,1)), 'M2', abs(Mz(1,2))^2, 'n01', p(1) - p(2), ...
            'Delta', E(2) - E(1), 'g0', 0);
tl.g0 = 2*tl.n01/tl.Delta;                              % g(0) = 2*n01*Delta/Delta^2
end

function test_bare_limit_and_fd_sign_anchor(testCase)
% GATE 1 (sign anchor, J 2.31; round-3 P0-1): the identity
%   -G0bare = M2*g0 + m^2*beta*(1-n01^2) = d(m*n01)/d(hmf)
% is exact ONLY for a closed two-level model -- built analytically here, NEVER from
% invz_twolevel_ordered (see twolevel_fix header).
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
D0 = 0.2;  a = 2.0;  M0 = 3.0;  hz = 0.15;  d = 1e-6;
tlf = @(h) twolevel_fix(h, D0, a, M0, beta);
tl  = tlf(hz);  tlp = tlf(hz + d);  tlm = tlf(hz - d);
[gi, ge] = twolevel_weights(tl, beta);
[Gs, out] = invz_gstat_ordered(tl, [0; 0], 0, 0, beta, gi, ge);
verifyEqual(testCase, out.xi, 1, 'AbsTol', 1e-14);                 % tanh(0)=0, denom 1
verifyEqual(testCase, Gs, out.G0bare, 'RelTol', 1e-14);            % elastic + inelastic bare
fd = (tlp.m*tlp.n01 - tlm.m*tlm.n01) / (2*d);                      % d<Jz>_closed/d hmf
verifyEqual(testCase, -out.G0bare, fd, 'RelTol', 1e-7);            % EXACT for the closed model
verifyEqual(testCase, out.r, 1, 'RelTol', 1e-14);                  % bare integrand ratio
end

function test_m_zero_recovers_resummed_pm_static(testCase)
% GATE 2 (J 2.30 at w=0): at m = 0 the elastic weight vanishes and
% Gstat = G0inel/(1+Sigma0+K0*G0inel);  GATE 3 (K-cancellation, HTML 32):
% Gtil0 = Gstat/(1-K0*Gstat) = G0inel/(1+Sigma0);  hence r = 1+Sigma0.
% Run in BOTH parametrizations: two-level weights AND a generic full-sector weight,
% since the m->0 algebra must hold for any G0inel0 (the hybrid's whole point, SS3).
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(1e-9);  tl.m = 0;                       % PM limit of the two-level params
Sigma0 = 0.31;  K0 = 0.018;  lam = [0.012; 0.004];      % generic nonzero test values
for G0inel = [-tl.M2*tl.g0, -251.4]                     % two-level and full-scale weights
    [Gs, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel, 0);
    verifyEqual(testCase, Gs, G0inel/(1 + Sigma0 + K0*G0inel), 'RelTol', 1e-13);
    verifyEqual(testCase, out.Gtil0, G0inel/(1 + Sigma0), 'RelTol', 1e-13);
    verifyEqual(testCase, out.r, 1 + Sigma0, 'RelTol', 1e-13);
end
end

function test_xi_formula_direct(testCase)
% GATE 4: xi transcription check against the closed formula (J 2.29) with
% hand-computed values -- independent arithmetic, not a copy of the source.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(0.15);
Sigma0 = 0.2;  K0 = 0.02;  lam = [0.01; 0.003];
[gi, ge] = twolevel_weights(tl, beta);
[~, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, gi, ge);
num = 1 + tanh(tl.m^2*tl.n01^2*beta*K0 - tl.M2*beta*lam(1));
den = 1 + (4*tl.n01^2*K0*tl.g0 + 2*lam(2) + tl.g0*lam(1))*tl.M2/tl.n01^2;
verifyEqual(testCase, out.xi, num/den, 'RelTol', 1e-14);
verifyGreaterThan(testCase, out.xi, 0);                 % bounded resummation stays positive
end
```

- [ ] **Step 2: Run to verify failure**

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests/test_invz_gstat_ordered.m'); disp(results); assertSuccess(results)"
```
Expected: all tests error — `invz_gstat_ordered` undefined.

- [ ] **Step 3: Implement**

Create `invz_common/invz_gstat_ordered.m`:

```matlab
function [Gstat, out] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0)
%INVZ_GSTAT_ORDERED Jensen ordered elastic static single-site function (framework SS9.2, J 2.28-2.29),
% in the BOUNDARY-PRESERVING HYBRID parametrization (stage-2 plan SS3, round-2 P0-A):
%   Gstat = G0inel0/(1 + Sigma0 + K0*G0inel0)  +  xi*G0el0
%   xi    = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1)))
%           / (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2)         (J 2.29)
% G-convention (G = -chi, meV^-1). The CALLER supplies the static weights: the production
% loops pass the FULL-electronuclear split (invz_chi0z elastic:false static, and
% elastic:true minus elastic:false), so the m -> 0 closure fixed point is the PM solver's
% own (round-2 P0-A); the transcription tests pass the two-level weights
% G0inel0 = -M2*g0, G0el0 = -m^2*h0, under which the first term's denominator becomes
% 1 + Sigma0 - M2*K0*g0 and Gstat is J 2.28 VERBATIM. Only xi is two-level (SS6).
% At w = 0 the K-cancellation does NOT go through -- Gstat inserts directly into the
% effective-medium form G(q,0) = Gstat/[1 + (J(q)-K0)*Gstat], and K0 itself requires the
% static closure of INVZ_EMT_STATIC_ORDERED (invz_emt_scalar's direct solve is
% ordinary-Dyson only).
% Identities pinned by test_invz_gstat_ordered (NEVER adjust signs to pass them, SS9):
%   bare (Sigma0=K0=lam=0):  xi = 1,  Gstat = G0bare = G0inel0 + G0el0,
%                            -G0bare = d<Jz>/d(hmf)                          (J 2.31)
%   m = 0 (G0el0 = 0):  Gstat = G0inel0/(1+Sigma0+K0*G0inel0);
%                       Gtil0 = G0inel0/(1+Sigma0);  r = 1+Sigma0   (any G0inel0)
% out: xi, h0 = beta*(1-n01^2), G0bare, Gtil0 = Gstat/(1-K0*Gstat), r = G0bare/Gtil0
% (the H_MF integrand of J 2.33).
m = tl.m;  M2 = tl.M2;  n01 = tl.n01;  g0 = tl.g0;
h0 = beta*(1 - n01^2);
xi = (1 + tanh(m^2*n01^2*beta*K0 - M2*beta*lam(1))) / ...
     (1 + (4*n01^2*K0*g0 + 2*lam(2) + g0*lam(1))*M2/n01^2);
Gstat  = G0inel0/(1 + Sigma0 + K0*G0inel0) + xi*G0el0;
G0bare = G0inel0 + G0el0;
Gtil0  = Gstat/(1 - K0*Gstat);
out = struct('xi', xi, 'h0', h0, 'G0bare', G0bare, 'Gtil0', Gtil0, 'r', G0bare/Gtil0);
end
```

Note for the implementer: check the repo's Boltzmann-constant usage before trusting the test's `kB` literal — grep `8.617` in `invz_common/`; if the codebase defines its own constant, use that in the tests instead (the FD gate's `beta` must match the one `invz_single_ion` uses).

- [ ] **Step 4: Run tests to verify pass; then the full fast suite**

Expected: 3/3 in the file (four gates across three test functions), then **148 passed / 0 failed / 19 incomplete** (145 + 3).

- [ ] **Step 5: Commit**

```bash
git add invz_common/invz_gstat_ordered.m invz_projected/tests/test_invz_gstat_ordered.m
git commit -m "feat(invz): ordered elastic static single-site function (J 2.28-2.29) with sign-anchor identity gates (stage-2 task 1)"
```

---

### Task 2: `invz_emt_static_ordered` — the static-sector EMT closure (P0-2)

**Files:**
- Create: `invz_projected/invz_emt_static_ordered.m`
- Test: `invz_projected/tests/test_invz_emt_static_ordered.m` (new file)

**Interfaces:**
- Consumes: Task 1's `invz_gstat_ordered` (hybrid signature); the branch-eigenvalue list `Jnu_flat` (same flattened `[nq*4]` vector `invz_emt_scalar` averages over); `invz_emt_scalar` (m→0 reference, in tests).
- Produces (Tasks 3–4 rely on): `[K0, Gstat, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu_flat, K0_seed, beta, J0eff, G0inel0, G0el0, opts)`:
  - `K0`: the closed static medium coupling (scalar, meV);
  - `Gstat`: the elastic static single-site function at the closed `K0` (meV⁻¹), built on the caller's static weights (production: full-electronuclear split — §3);
  - `out`: `xi`, `h0`, `G0bare`, `Gtil0`, `r` (forwarded from `invz_gstat_ordered` at the closed `K0`), plus `D_uni = 1 + (J0eff − K0)*Gstat` (the uniform static inverse response, the pole observable of §3), `resid` (the closure residual `|⟨G(q,0)⟩_q − Gstat|` at exit), `iters`, `converged` (logical).
  - opts (threaded from the caller as ONE NAMED STRUCT `emt_static` — P1-F, never a bare `struct()` at call sites): `resid_tol` (default 1e-10, meV⁻¹ — the PRIMARY convergence criterion AND the acceptance threshold callers gate on: the iteration stops, and `converged` is true, only when the closure residual `|⟨G(q,0)⟩_q − Gstat| < resid_tol`), `tol` (default 0 — an OPTIONAL absolute `|ΔK0|` stall floor; the built-in stall test is machine resolution `4·eps(|K0|)`, because `resid ≈ 1.8e5·|ΔK0|` at measured conditioning means any coarser floor preempts the residual criterion), `maxit` (default 200), `mix` (default 0.5). `converged` is measured on the EXPORTED tuple (post-loop recompute), so the flag and `out.resid` can never disagree. *(Execution amendments, 2026-07-22: (1) the original ΔK0-based stopping under-delivered closure by ~5 orders; (2) the first residual-based fix left a 1e-12 dK floor that structurally preempted the 1e-10 residual target and evaluated `converged` one half-step behind the export — both found by the Task-2 implementer's escalations, both verified numerically before this final form.)*

- [ ] **Step 1: Write the failing tests**

Create `invz_projected/tests/test_invz_emt_static_ordered.m`:

```matlab
function tests = test_invz_emt_static_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function tl = ordered_tl(hz)
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
end

function test_m_zero_reproduces_direct_solve(testCase)
% GATE C1 (exactness anchor, round-2 P0-A form): at m = 0 (elastic weight 0) the hybrid
% degenerates to the ordinary Dyson structure ON WHATEVER G0inel0 the caller supplies --
% invz_emt_scalar's direct solve is exact there. Check at BOTH a two-level-scale and a
% full-electronuclear-scale weight (the measured PM value is ~ -251 meV^-1): the closure
% must reproduce the direct solve's K(0) and G(0) for the SAME (G0, Sigma) inputs.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(1e-9);  tl.m = 0;
Sigma0 = 0.27;  lam = [0.012; 0.004];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
for G0inel = [-tl.M2*tl.g0, -251.4]
    med = invz_emt_scalar(G0inel, Sigma0, Jnu, struct());     % one-frequency direct solve
    [K0, Gs, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, 0, beta, 6.4e-3, ...
                                            G0inel, 0, struct());
    verifyTrue(testCase, out.converged);
    verifyEqual(testCase, K0, med.K(1), 'RelTol', 1e-9);
    verifyEqual(testCase, Gs, med.G(1), 'RelTol', 1e-9);
end
end

function test_closure_residual_is_enforced(testCase)
% GATE C2: at the returned K0 the EMT consistency holds -- the q-average of the
% elastic-inserted G(q,0) equals the site Gstat (the "R eq 8" convention) -- with
% genuinely ordered (m > 0) hybrid weights of full-electronuclear scale.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
tl = ordered_tl(0.15);                                    % genuinely ordered (m > 0)
Sigma0 = 0.2;  lam = [0.01; 0.003];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
G0inel = -230.0;  G0el = -12.0;                           % full-scale test weights
[K0, Gs, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu, 0, beta, 6.4e-3, ...
                                        G0inel, G0el, struct());
verifyTrue(testCase, out.converged);
Gq = Gs ./ (1 + (Jnu - K0).*Gs);
verifyEqual(testCase, mean(Gq), Gs, 'RelTol', 1e-9);      % independent recomputation
verifyLessThan(testCase, out.resid, 1e-10);
verifyEqual(testCase, out.D_uni, 1 + (6.4e-3 - K0)*Gs, 'RelTol', 1e-14);
end

function test_bare_limit_closure(testCase)
% GATE C3: with Sigma0 = 0, lam = 0 the closure agrees with the direct solve at the
% m = 0 degenerate limit and stays finite/converged at finite m; opts threading works.
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
Jnu = linspace(-2e-3, 6.0e-3, 24).';
tl = ordered_tl(0.15);
[gi, ge] = deal(-tl.M2*tl.g0, -tl.m^2*beta*(1 - tl.n01^2));
[K0m, Gsm, outm] = invz_emt_static_ordered(tl, [0;0], 0, Jnu, 0, beta, 6.4e-3, gi, ge, ...
                                           struct('resid_tol', 1e-12, 'maxit', 400));
verifyTrue(testCase, outm.converged && isfinite(K0m) && isfinite(Gsm));
verifyLessThan(testCase, outm.resid, 1e-11);              % tighter resid_tol was honored
tl0 = ordered_tl(1e-9);  tl0.m = 0;
med = invz_emt_scalar(-tl0.M2*tl0.g0, 0, Jnu, struct());
[K00, ~, ~] = invz_emt_static_ordered(tl0, [0;0], 0, Jnu, 0, beta, 6.4e-3, ...
                                      -tl0.M2*tl0.g0, 0, struct());
verifyEqual(testCase, K00, med.K(1), 'RelTol', 1e-9);
end
```

- [ ] **Step 2: Run to verify failure** (function undefined).

- [ ] **Step 3: Implement**

Create `invz_projected/invz_emt_static_ordered.m`:

```matlab
function [K0, Gstat, out] = invz_emt_static_ordered(tl, lam, Sigma0, Jnu_flat, K0_seed, beta, J0eff, G0inel0, G0el0, opts)
%INVZ_EMT_STATIC_ORDERED Static-sector EMT closure for the ordered elastic propagator (SS3, P0-2/P0-A).
% The elastic G(0) of J 2.28-2.29 breaks the ordinary Dyson structure, so the closed-form
% direct solve of INVZ_EMT_SCALAR does not apply at w = 0 for m ~= 0. This function solves
% the scalar fixed point demanded by the EMT definitions (J 2.10-2.11 / HTML 16):
%   Gstat(K0)  = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0)
%   G(q,0)     = Gstat ./ (1 + (J(q) - K0).*Gstat)                       (HTML 17 insertion)
%   closure:   mean_q G(q,0) = Gstat   and   K0 = mean_q(J.*Gq)/mean_q(Gq)
% by damped iteration on K0. The static weights are the caller's: production passes the
% FULL-electronuclear split (round-2 P0-A), so at m = 0 (G0el0 -> 0) the fixed point
% coincides with invz_emt_scalar's direct solve ON THE FULL PROPAGATOR -- the PM solver's
% own K(0) (gate C1; the load-bearing loop-level identity is Task 3's Gate 6b). lam and
% Sigma0 are FIXED inputs here; the caller's outer loop refreshes them with the closed K0
% written into the K vector's w = 0 slot (see invz_hmf_ordered / the jensen solver mode).
% out: xi/h0/G0bare/Gtil0/r at the closed K0 (from invz_gstat_ordered), D_uni = the uniform
% static inverse response 1 + (J0eff - K0)*Gstat (pole observable), resid, iters, converged.
if nargin < 10, opts = struct(); end
rtol  = getf(opts, 'resid_tol', 1e-10);  % PRIMARY: closure-residual convergence (meV^-1)
tol   = getf(opts, 'tol', 0);            % optional ABSOLUTE |dK0| stall floor (0 = machine-only)
maxit = getf(opts, 'maxit', 200);
mix   = getf(opts, 'mix', 0.5);
Jf = Jnu_flat(:);
K0 = K0_seed;
for it = 1:maxit
    Gs = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
    Gq = Gs ./ (1 + (Jf - K0).*Gs);
    Gbar = mean(Gq);
    if abs(Gbar - Gs) < rtol, break; end % closed at the CURRENT K0 -- exported as-is
    K0_new = mean(Jf .* Gq) / Gbar;
    dK = abs(K0_new - K0);
    if dK < max(tol, 4*eps(abs(K0)))     % TRUE stall: no representable progress possible
        break;                           % (resid ~ 1.8e5*dK at measured conditioning, so an
    end                                  % arbitrary dK floor must NOT preempt the residual)
    K0 = K0 + mix*(K0_new - K0);
end
[Gstat, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
Gq = Gstat ./ (1 + (Jf - K0).*Gstat);
out = go;
out.D_uni = 1 + (J0eff - K0)*Gstat;
out.resid = abs(mean(Gq) - Gstat);
out.iters = it;
out.converged = out.resid < rtol;        % measured on the EXPORTED tuple -- the flag and the
if ~out.converged                        % residual can never disagree (execution amendment 2)
    warning('invz:emtStatic', 'static closure not converged after %d iterations: resid = %.3g', it, out.resid);
end
end
```

Reviewer note on the update rule: at the fixed point, `mean_q Gq = Gstat` (self-average closure) together with `K0 = mean(J.*Gq)/mean(Gq)` reproduces the definition J 2.11 with the site function as denominator. At m=0 substituting the ordinary `Gs(K0)` shows the direct solve's `(K, G)` pair satisfies both — gate C1 checks this to 1e-9.

- [ ] **Step 4: Run the tests (3/3), then the fast suite: expected 151 passed / 0 failed / 19 incomplete** (148 + 3).

- [ ] **Step 5: Commit**

```bash
git add invz_projected/invz_emt_static_ordered.m invz_projected/tests/test_invz_emt_static_ordered.m
git commit -m "feat(invz): static-sector EMT closure for the ordered elastic propagator -- direct-solve degeneracy and closure-residual gates (stage-2 task 2)"
```

---

### Task 3: `invz_hmf_ordered` — the applied-field/H_MF relation and spontaneous root (J 2.31–2.33)

**Files:**
- Create: `invz_projected/invz_hmf_ordered.m`
- Test: `invz_projected/tests/test_invz_hmf_ordered.m` (new file)

**Interfaces:**
- Consumes: Tasks 1–2 (`invz_gstat_ordered` via `invz_emt_static_ordered`); `invz_single_ion(..., 'hz_fixed', h')` — **WITHOUT `order`** (§1); `invz_twolevel_ordered`; `invz_chi0z`; `invz_g`; `invz_emt_scalar`; `invz_lambdas`; `invz_sigma_ordered`. The Matsubara grid/weights/beta construction must MIRROR `invz_projected/invz_solve_point_ordered.m` — copy its setup block (the lines building `wn`, `wts`, `beta`, `eopts` from opts; currently the block preceding line 126) verbatim, honoring the caller's `Ecut` and EMT options (P1-6).
- Per-node solve loop (the ordered Σ↔EMT loop with the static-sector correction): each outer pass runs, in order: (1) `med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts); K = med.K;` (dynamic sector, all ωn); (2) `[K0s, Gs, sout] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, beta, J0eff, G0inel0, G0el0, eso);` then `K(1) = K0s;` (static sector on the FULL-electronuclear weights — P0-A; `eso` is the threaded `emt_static` opts struct — P1-F; on the first pass use `lam = [0;0;0]`); (3) `lam = invz_lambdas(K, g, wts, beta, [1 2 3]);` (4) `sg = invz_sigma_ordered(tl, lam, K, g, beta);` (5) damped mix of `Sigma` (`mixo`), converge on `max(abs(sg.Sigma − Sigma)) < tolo` AND the static closure converged. The static weights per node (round-4 P2-E — this supersedes the round-2 `G0(1)` difference): `G0inel0 = -real(c0i(3,3,1))` with `c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false))` (one-frequency call), and `G0el0 = G0bare0 − G0inel0` where `G0bare0` is the PATH-CONSISTENT chain-rule derivative of §4a (mode-switched: `'none'` has no feedback, `'legacy_x'` the x channel, `'vector_ab'` the 2×2 block — round-4 P1-A). The final post-loop static refresh KEEPS its newly closed `K0` (first output) and writes it to `K(1)` — the reported `K0k` must be the value the returned `Gstat/r/D_uni` were computed with, never the pre-refresh seed (round-4 P1-B). Warm-start `Sigma` (and `K0s`) from the previous node.
- Produces (Tasks 4/6 rely on): `[hmf_star, prof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, opts)`:
  - `hmf_star`: scalar meV, the nonzero spontaneous root of `F(h) = h0(h) − J0eff·m(h)`, refined by direct evaluation (below); `NaN` when no nonzero root exists;
  - `prof` struct: `hgrid` (meV, GEOMETRIC, clustered at 0 — P1-4, possibly extended below `hmin_initial` by the adaptive path), `r`, `h0`, `m` (full moments `si.Jexp(3)`), `Sigma0`, `K0` (closed static values), `D_uni`, `G0bare` (path-consistent, §4a), `Gstat` (closed elastic static per node — P1-4 discrimination), `node_conv`, `F`; the h=0 predictor node's `slope0` (= the §5 mass, ≈ PM `crit` — Gate 6b comparator), `Sigma0_pm0`, `K0_pm0`; the adaptivity record `n_extend`, `hmin_initial`, `status` (`'ok' | 'unresolved' | 'node_failed' | 'no_bare_order'` — round-5 P2: **`'unresolved'` means ordering was predicted but no bracket was found above `hmin_abs` OR the bracket did not refine to `tol_root`; `'node_failed'` means a profile, bisection, or final-evaluation node failed to converge/close. The jensen solver MUST map BOTH to `converged = false`, never to a PM label** — round-3 P0-3 / round-5 P1-A); and at the refined root `m_star`, `D_uni_star`, `r_star`, `Gstat_star`;
  - opts: `nH` (default 33), `hmax_fac` (default 1.25 — grid ceiling = `hmax_fac ×` the BARE solver's converged `si.hz`, obtained by one SEPARATE `order=true` solve BUILT FROM THE SAME MF OPTION BASE plus `order=true, J0z=J0eff` — P1-F; if bare does not order, return `NaN` immediately), `hmin_frac` (default 1e-3 — `hgrid = hmax·ratio.^((nH-1):-1:0)` with `ratio = hmin_frac^(1/(nH-1))`), `hmin_abs` (default `1e-10*hmax` — the ADAPTIVE lower-extension floor, P1-C: whenever `slope0 < 0` predicts ordering but no negative `F` sample exists, the grid extends geometrically downward until a negative sample is found or the floor is hit), `tol_root` (default 1e-3 relative, bisection by direct node evaluation, ≤ 12 evaluations), `J0eff` (required), `emt_static` (opts struct threaded verbatim to every `invz_emt_static_ordered` call — P1-F), plus the full numerical context forwarded to every node: `Jxx0`, `hyp`, `transverse_mf`, `Ecut`, `mf_maxit`, `mf_mix`, `mix_outer`, `tol_outer`, `max_outer`, EMT options (P1-6), and `force_bare` (test hook: `Σ=0, K=0, λ=0, r≡1`).

- [ ] **Step 1: Write the failing tests**

Create `invz_projected/tests/test_invz_hmf_ordered.m`:

```matlab
function tests = test_invz_hmf_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function o = syn_opts(extra)
o = struct('J0eff', 6.4e-3, 'Jxx0', [], 'hyp', true, 'transverse_mf', 'legacy_x');
ion = invz_ion();  o.Jxx0 = ion.Jxx0;
if nargin, f = fieldnames(extra); for k = 1:numel(f), o.(f{k}) = extra.(f{k}); end, end
end

function test_bare_hook_reproduces_bare_mf_root(testCase)
% GATE 5 (wiring): with force_bare the relation collapses to h0 = hmf, so the root
% must satisfy hmf = J0eff*m(hmf) -- the bare MF fixed point invz_single_ion finds.
% Also pins the fixed-field contract: every node held its imposed hz (P0-1).
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[hstar, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, syn_opts(struct('force_bare', true)));
verifyTrue(testCase, isfinite(hstar) && hstar > 0);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'order', true, ...
        'J0z', 6.4e-3, 'Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyEqual(testCase, hstar, si.hz, 'RelTol', 1e-3);        % same bare fixed point
verifyEqual(testCase, prof.r, ones(size(prof.r)), 'AbsTol', 1e-14);
end

function test_onset_coincides_with_pm_crit_zero(testCase)
% GATE 6 (THE closure theorem, SS5): the nonzero root exists at 2.85 T (PM crit
% -0.058) and NOT at 3.30 T (PM crit +0.087) -- the same Bc_1z the PM mass defines.
% Anchors from the stage-1 window probe (re-verify if the synthetic model changes).
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[h_lo, p_lo] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts());
[h_hi, ~]    = invz_hmf_ordered(ion, T, [3.30 0 0], Jnu, syn_opts());
verifyTrue(testCase, isfinite(h_lo) && h_lo > 0, 'ordered root missing below Bc_1z');
verifyTrue(testCase, isnan(h_hi), 'spurious ordered root above Bc_1z');
verifyTrue(testCase, all(p_lo.node_conv), 'non-converged H_MF nodes at 2.85 T');
% small-field slope diagnostic ~ crit: negative below Bc_1z (SS5)
verifyLessThan(testCase, p_lo.slope0, 0);
% fluctuation suppression: the 1/z root sits BELOW the bare root, with a smaller moment
si = invz_single_ion(ion, T, [2.85 0 0], struct('hyp', true, 'order', true, ...
        'J0z', 6.4e-3, 'Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyLessThan(testCase, h_lo, si.hz);
verifyLessThan(testCase, p_lo.m_star, si.Jexp(3));
% pole observable at the root: ordered static inverse response is small and positive
% approaching closure from the ordered side (vanishes AT Bc_1z)
verifyTrue(testCase, isfinite(p_lo.D_uni_star));
end

function test_moment_vanishes_continuously_toward_onset(testCase)
% GATE 7: m(B) decreases toward the onset -- the continuity stage 1 lacked.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
Bs = [2.85 2.95 3.05];  ms = nan(size(Bs));
for k = 1:numel(Bs)
    [h, p] = invz_hmf_ordered(ion, T, [Bs(k) 0 0], Jnu, syn_opts());
    if isfinite(h), ms(k) = p.m_star; end
end
fin = isfinite(ms);
verifyGreaterThanOrEqual(testCase, nnz(fin), 2);
verifyTrue(testCase, all(diff(ms(fin)) < 0));               % monotone decrease toward Bc_1z
end

function test_grid_convergence(testCase)
% GATE 8 (P1-4): the refined root is grid-converged -- doubling nH moves hmf_star by
% less than 1e-2 relative. This is the promised nH-refinement probe, as an actual test.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[h1, ~] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('nH', 21)));
[h2, ~] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, syn_opts(struct('nH', 42)));
verifyEqual(testCase, h1, h2, 'RelTol', 1e-2);
end

function test_pm_limit_reproduces_solve_point(testCase)
% GATE 6b (round-2 P0-A, LOAD-BEARING): with hyp = true, the ordered machinery's
% DEDICATED h = 0 predictor node must land on the ACTUAL PM solver's fixed point --
% K(0), Sigma0, and the inverse mass -- not an isolated two-level reference. This is
% what makes the ordered onset coincide with the code's Bc_1z (SS5). PM-side field.
ion = invz_ion();  T = 0.31;  Bx = [3.30 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
ptpm = invz_solve_point(ion, T, Bx, Jnu, o);
verifyTrue(testCase, ptpm.converged);
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, ...
    struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'nH', 5));
verifyEqual(testCase, prof.K0_pm0,     ptpm.K(1),   'RelTol', 1e-4);   % same medium
verifyEqual(testCase, prof.Sigma0_pm0, ptpm.Sigma0, 'RelTol', 1e-4);   % same self-energy
% same inverse mass: the h = 0 predictor slope equals the PM crit (SS5). NOTE: the
% cross response X_zt vanishes at hz = 0, so the path correction is inert here --
% this gate is sector-identity, not path-convention (that is the FD gate below).
verifyEqual(testCase, prof.slope0, ptpm.crit, 'RelTol', 1e-3);
end

function test_path_derivative_gate(testCase)
% Production FD gate (round-3 P0-2, SS4a; round-4 P1-A: ALL THREE transverse-MF modes):
% at a FINITE-moment node the chain-rule G0bare must equal the centered FD of the
% ACTUAL node-path moments under the SAME mode -- 'none' has NO feedback (fb = 0),
% 'legacy_x' the x channel, 'vector_ab' the {x,y} block. Small nH keeps this fast.
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
for tmfc = {'none', 'legacy_x', 'vector_ab'}
    tmf = tmfc{1};
    [~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, ...
        syn_opts(struct('transverse_mf', tmf, 'nH', 7)));
    k = find(prof.hgrid > 0.5*max(prof.hgrid), 1);    % a well-ordered node
    hp = prof.hgrid(k);  d = 1e-6;
    sib = struct('hyp', true, 'Jxx0', ion.Jxx0, 'transverse_mf', tmf);
    sp = sib;  sp.hz_fixed = hp + d;  sm = sib;  sm.hz_fixed = hp - d;
    mp = invz_single_ion(ion, T, Bx, sp);  mm = invz_single_ion(ion, T, Bx, sm);
    fd = (mp.Jexp(3) - mm.Jexp(3)) / (2*d);           % path derivative by FD
    verifyEqual(testCase, -prof.G0bare(k), fd, 'RelTol', 1e-3, ...
        sprintf('path-derivative mismatch under transverse_mf = %s', tmf));
end
end

function test_failure_contract(testCase)
% Round-5 P1-A hook: the failure paths must be EXERCISED, not only documented. With an
% absurdly small outer budget the Sigma loop cannot converge, so every node fails and
% the profile must report 'node_failed' with hmf_star = NaN -- never an 'ok' status or
% a root built from unconverged nodes.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[h, p] = invz_hmf_ordered(ion, T, [2.85 0 0], Jnu, ...
    syn_opts(struct('max_outer', 1, 'nH', 5)));
verifyTrue(testCase, isnan(h));
verifyEqual(testCase, p.status, 'node_failed');
end

function test_near_boundary_root_not_missed(testCase)
% P1-C / round-3 P0-3 (blind interval, NON-CIRCULAR): refine Bc_1z from the PM
% crit = 0 side (INDEPENDENT of the root detector under test), pick a field just
% below it, and prove the ADAPTIVE EXTENSION actually ran: a deliberately coarse
% hmin_frac must miss the tiny root in its initial grid (n_extend >= 1) and still
% find it after extension. Several profiles -> INVZ_SLOW-gated.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (multi-profile).');
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
% independent Bc_1z: PM crit interpolated to zero on a small field grid (crit is the
% PM solver's own mass -- NOT invz_hmf_ordered)
Bg = 2.90:0.05:3.20;  cg = nan(size(Bg));
for k = 1:numel(Bg)
    pt = invz_solve_point(ion, T, [Bg(k) 0 0], Jnu, o);
    if pt.converged, cg(k) = pt.crit; end
end
ok = isfinite(cg);
Bc_pm = interp1(cg(ok), Bg(ok), 0, 'linear', 'extrap');
Btest = Bc_pm - 0.01;                                   % just below the PM-refined boundary
% round-4 P2-D (DETERMINISTIC forcing): first obtain a fine-grid REFERENCE root, then
% CONSTRUCT the coarse lower limit to sit provably above it, so n_extend >= 1 is a
% guaranteed consequence of correct adaptivity, not luck.
[h_ref, p_ref] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, syn_opts());
verifyTrue(testCase, isfinite(h_ref) && h_ref > 0, 'no reference root just below Bc_pm');
hmax_ref = max(p_ref.hgrid);
frac = min(0.5, 4*h_ref/hmax_ref);                      % initial lower node > 4x the root
[h, p] = invz_hmf_ordered(ion, T, [Btest 0 0], Jnu, ...
    syn_opts(struct('hmin_frac', frac, 'nH', 17)));
verifyLessThan(testCase, p.slope0, 0);                  % ordering predicted (independent)
verifyGreaterThan(testCase, p.hmin_initial, h_ref);     % coarse limit provably above the root
verifyGreaterThanOrEqual(testCase, p.n_extend, 1);      % the adaptive path actually ran
verifyTrue(testCase, isfinite(h) && h > 0 && h < p.hmin_initial, ...
    sprintf('root missed at B=%.3f (Bc_pm=%.3f): h=%.3g, hmin_initial=%.3g, n_extend=%d', ...
            Btest, Bc_pm, h, p.hmin_initial, p.n_extend));
verifyEqual(testCase, h, h_ref, 'RelTol', 5e-2);        % same root as the fine reference
end
```

- [ ] **Step 2: Run to verify failure** (function undefined).

- [ ] **Step 3: Implement**

Create `invz_projected/invz_hmf_ordered.m` (structure below; the two `MIRROR` markers are instructions to copy the named source blocks byte-for-byte — the repo's literal-copy convention — with the Task-3 Interfaces bullet spelling out the loop's exact statement order including the static-closure insertion):

```matlab
function [hmf_star, prof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_HMF_ORDERED Jensen applied-field/H_MF self-consistency, spontaneous root (SS9.3, J 2.31-2.33).
%   h0(hmf) = int_0^hmf r(h') dh',   r = G0(0;h')/Gtil0(0;h')
% with Gtil0 built on the STATIC-CLOSURE K0 (invz_emt_static_ordered, P0-2), evaluated on
% fixed-field single-ion states (hz_fixed WITHOUT order -- P0-1, invz:hzFixed asserted).
% Spontaneous condition (zero applied longitudinal field): h0(hmf) = J0eff*<Jz>(hmf); the
% nonzero root is bracketed on a GEOMETRIC profile clustered at 0 (P1-4) and refined by
% bisection with DIRECT node evaluations to opts.tol_root. F(h)/h -> crit as h -> 0+
% (SS5), returned as prof.slope0. Returns NaN when no nonzero root exists, or when the
% separate bare (order=true) bracketing solve does not order.
if nargin < 5, opts = struct(); end
J0eff = opts.J0eff;                                  % required, no default (caller-owned)
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
hyp   = getf(opts, 'hyp', true);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
nH    = getf(opts, 'nH', 33);
hfac  = getf(opts, 'hmax_fac', 1.25);
hfrac = getf(opts, 'hmin_frac', 1e-3);
trt   = getf(opts, 'tol_root', 1e-3);
fbare = getf(opts, 'force_bare', false);
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 200);                % ALIGNED with both solvers' default
                                                     % (round-3 P1-5; verify at implementation
                                                     % time and keep identical)
eso   = getf(opts, 'emt_static', struct());          % static-closure opts, threaded (P1-F)

% single-ion opts for FIXED-FIELD nodes: NO 'order' (P0-1); forward mf knobs (P1-6)
sibase = struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), sibase.(f{1}) = opts.(f{1}); end
end
hmin_abs = getf(opts, 'hmin_abs', NaN);              % resolved after hmax below (P1-C)

prof = struct('hgrid', [], 'r', [], 'h0', [], 'm', [], 'Sigma0', [], 'K0', [], ...
              'D_uni', [], 'G0bare', [], 'Gstat', [], 'node_conv', [], 'F', [], ...
              'slope0', NaN, 'Sigma0_pm0', NaN, 'K0_pm0', NaN, 'J0eff', J0eff, ...
              'n_extend', 0, 'hmin_initial', NaN, 'status', 'no_bare_order', ...
              'redensified', false, ...
              'm_star', NaN, 'D_uni_star', NaN, 'r_star', NaN, 'Gstat_star', NaN);
hmf_star = NaN;
% Implementation note (amendment 3): factor the per-node nH sweep into a nested helper
%   [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s)
% used by BOTH the initial sweep and the re-densify pass (behavior-neutral refactor of
% the inline loop below; the extension's 3-node prepends may stay inline).

% Bracket ceiling from the BARE ordered fixed point: SAME MF option base plus
% order=true and J0z (P1-F -- the bracket runs under the caller's MF knobs too).
sibo = sibase;  sibo.order = true;  sibo.J0z = J0eff;
sib = invz_single_ion(ion, T, Bx, sibo);
if ~(sib.mf_converged && abs(sib.Jexp(3)) > 1e-6), return; end     % bare does not order
hmax = hfac * abs(sib.hz);
if isnan(hmin_abs), hmin_abs = 1e-10*hmax; end

% --- Matsubara grid, weights, beta: MIRROR invz_solve_point_ordered's setup block
% verbatim (wn, wts, beta, eopts from opts -- honors Ecut and EMT options, P1-6).
% [transcribe here]

% Independent h = 0 PM predictor node (round-3 P0-3; doubles as Gate 6b's comparator):
% ONE node solve at hz_fixed = 0 gives THIS machinery's PM fixed point. Its mass
%   slope_pred = r(0) + J0eff*G0bare(0) = 1 + Sigma0(0) - J0eff*chi_path(0)   (= crit, SS5)
% predicts root existence INDEPENDENTLY of any sampled profile value.
Sigma = [];  K0s = 0;                                % warm-start carriers across nodes
[r0n, ~, S0pm, K0pm, ~, Gb0, ~, ok0, Sigma, K0s] = eval_node(0, Sigma, K0s);
if ~ok0, prof.status = 'node_failed'; return; end     % a failed PREDICTOR node is a node
                                                      % failure, not an unrefined bracket
                                                      % (as-built; test_failure_contract)
slope_pred = r0n + J0eff*Gb0;
prof.Sigma0_pm0 = S0pm;  prof.K0_pm0 = K0pm;  prof.slope0 = slope_pred;

ratio = hfrac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);                 % geometric, clustered at 0 (P1-4)
prof.hmin_initial = hgrid(1);

[rv, mv, S0v, K0v, Dv, Gbv, Gsv] = deal(nan(1, nH));  cnv = false(1, nH);
for k = 1:nH
    [rv(k), mv(k), S0v(k), K0v(k), Dv(k), Gbv(k), Gsv(k), cnv(k), Sigma, K0s] = ...
        eval_node(hgrid(k), Sigma, K0s);
end

h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);  % first panel seeded with the EXACT r(0)
F  = h0 - J0eff*mv;

% ADAPTIVE lower extension (round-3 P0-3): predictor-driven, NOT self-referential.
% slope_pred < 0 predicts an ordered root; extend geometrically downward until a
% negative F sample appears or the absolute floor is reached.
n_extend = 0;
while slope_pred < 0 && all(F >= 0) && hgrid(1) > hmin_abs
    n_extend = n_extend + 1;
    hext = hgrid(1) * ratio.^(3:-1:1);                % three more decades-fraction nodes
    [re, me, S0e, K0e, De, Gbe, Gse] = deal(nan(1, 3));  ce = false(1, 3);
    for k = 1:3
        [re(k), me(k), S0e(k), K0e(k), De(k), Gbe(k), Gse(k), ce(k), Sigma, K0s] = ...
            eval_node(hext(k), Sigma, K0s);
    end
    hgrid = [hext hgrid];  rv = [re rv];  mv = [me mv];  cnv = [ce cnv];
    S0v = [S0e S0v];  K0v = [K0e K0v];  Dv = [De Dv];  Gbv = [Gbe Gbv];  Gsv = [Gse Gsv];
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
end
if n_extend > 0 && any(F < 0)
    % RE-DENSIFY (execution amendment 3, 2026-07-22): the extension's sparse geometric
    % panels feed O(coarse-grid) quadrature error into h0 exactly where F is a small
    % difference of large terms (measured: 11% root error at Bc_1z - 0.01 on a
    % deliberately coarse grid vs the fine default). Rebuild the profile at FULL nH
    % resolution anchored to the discovered bracket scale, so adaptive-path roots match
    % default-path quality. Cost: one extra nH-sweep, only when extension fired.
    idx0 = find(F < 0, 1, 'first');
    hfrac_eff = max(hmin_abs/hmax, 0.25*hgrid(idx0)/hmax);
    ratio2 = hfrac_eff^(1/(nH-1));
    hgrid = hmax * ratio2.^((nH-1):-1:0);
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s);
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
    prof.redensified = true;
end
prof.n_extend = n_extend;

prof.hgrid = hgrid;  prof.r = rv;  prof.h0 = h0;  prof.m = mv;
prof.Sigma0 = S0v;   prof.K0 = K0v;  prof.D_uni = Dv;  prof.node_conv = cnv;  prof.F = F;
prof.G0bare = Gbv;   prof.Gstat = Gsv;

if slope_pred < 0 && all(F >= 0)                      % floor hit without a bracket:
    prof.status = 'unresolved';                       % NEVER silently PM (round-3 P0-3)
    warning('invz:hmfUnresolved', ...
        'ordering predicted (slope_pred = %.3g) but no negative F above hmin_abs = %.3g', ...
        slope_pred, hmin_abs);
    return;                                           % hmf_star stays NaN; the jensen solver
end                                                   % must return converged = false here
if any(~cnv)                                          % round-4 P1-C: status must be truthful
    prof.status = 'node_failed';                      % on node failure -- never 'ok'
    return;
end
prof.status = 'ok';
s = sign(F);  idx = find(s(1:end-1) < 0 & s(2:end) >= 0, 1, 'last');
if isempty(idx), return; end                          % no nonzero root: PM side

% --- Root refinement by DIRECT evaluation (P1-4): bisection on F between the
% bracketing nodes, fresh node solve per iterate, cumulative h0 via local trapezoid
% panel from the bracket's left node.
a = hgrid(idx);  b = hgrid(idx+1);  Fa = F(idx);  h0a = h0(idx);  ra = rv(idx);
for it = 1:12
    c = 0.5*(a + b);
    [rc, mc, ~, ~, Dc, ~, ~, okc, Sigma, K0s] = eval_node(c, Sigma, K0s);
    if ~okc                                           % round-5 P1-A: a failed bisection node
        prof.status = 'node_failed';  hmf_star = NaN; % TERMINATES the solve -- never a root
        return;                                       % from a partial bracket
    end
    h0c = h0a + 0.5*(ra + rc)*(c - a);
    Fc  = h0c - J0eff*mc;
    if sign(Fc) == sign(Fa), a = c; Fa = Fc; h0a = h0c; ra = rc; else, b = c; end
    if (b - a) < trt*b, break; end
end
if (b - a) >= trt*b                                   % round-5 P1-A: tol_root not reached --
    prof.status = 'unresolved';  hmf_star = NaN;      % a distinct refinement failure
    warning('invz:hmfUnresolved', 'root bracket not refined to tol_root: (b-a)/b = %.3g', (b-a)/b);
    return;
end
hmf_star = 0.5*(a + b);
[r_s, m_s, ~, ~, D_s, ~, Gs_s, ok_s] = eval_node(hmf_star, Sigma, K0s);
if ~ok_s
    prof.status = 'node_failed';  hmf_star = NaN;  return;
end
prof.m_star = m_s;  prof.D_uni_star = D_s;  prof.r_star = r_s;  prof.Gstat_star = Gs_s;

    function [rk, mk, S0k, K0k, Dk, Gbk, Gsk, ok, Sigma, K0s] = eval_node(hp, Sigma, K0s)
    % One fixed-field node: si (hz_fixed, NO order), tl, c0/G0, then the ordered
    % Sigma<->EMT loop WITH the static-sector closure each pass (Interfaces bullet).
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bx, sio);
    if abs(si.hz - hp) > 1e-12
        error('invz:hzFixed', 'hz_fixed not held: si.hz = %.6g vs %.6g', si.hz, hp);
    end
    tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
    mk = si.Jexp(3);
    if fbare
        rk = 1;  S0k = 0;  K0k = 0;  Dk = NaN;  Gbk = NaN;  Gsk = NaN;  ok = true;  return;
    end
    c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
    G0 = -real(squeeze(c0(3,3,:)));
    c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));   % static inelastic only
    G0inel0 = -real(c0i(3,3,1));                                   % fixed-Hamiltonian slot
    % PATH-CONSISTENT bare static (SS4a, round-3 P0-2; round-4 P1-A): -dm/dh ALONG the
    % node path. Transverse-MF feedback by mode: 'none' forces hx = hy = 0 in
    % invz_single_ion, so the path derivative is X_zz alone; 'legacy_x' feeds back the
    % x channel; 'vector_ab' the 2x2 {x,y} block.
    X = real(c0(:, :, 1));                                         % static chi tensor (chi = -G)
    switch tmf
        case 'none'
            fb = 0;
        case 'legacy_x'
            fb = X(3, 1) * (Jxx0 / (1 - Jxx0*X(1, 1))) * X(1, 3);
        case 'vector_ab'
            t = [1 2];
            fb = X(3, t) * (Jxx0 * ((eye(2) - Jxx0*X(t, t)) \ X(t, 3)));
        otherwise
            error('invz:transverseMF', 'unknown transverse_mf ''%s''', tmf);
    end
    G0bare0 = -(X(3, 3) + fb);
    G0el0   = G0bare0 - G0inel0;                                   % elastic + feedback (SS4a)
    g  = real(invz_g(tl, 1i*wn));
    if isempty(Sigma), Sigma = zeros(size(wn)); end
    K = zeros(size(wn));  lam = [0; 0; 0];  ok = false;
    for outer = 1:maxo
        % (1) dynamic sector -- MIRROR invz_solve_point_ordered's emt call verbatim
        % [transcribe here: eopts.K0 = K; med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts); K = med.K;]
        % (2) static sector (P0-2/P0-A), threaded opts (P1-F):
        [K0s, ~, sout] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
                                                 beta, J0eff, G0inel0, G0el0, eso);
        K(1) = K0s;
        % (3)-(5) lambdas, ordered Sigma, damped mix -- MIRROR the solver's statements
        % [transcribe here: lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
        %  sg = invz_sigma_ordered(tl, lam, K, g, beta);
        %  dS = max(abs(sg.Sigma - Sigma)); Sigma = Sigma + mixo*(sg.Sigma - Sigma);
        %  if dS < tolo && sout.converged, ok = true; break; end]
    end
    [K0s, Gsk, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
                                             beta, J0eff, G0inel0, G0el0, eso);
    K(1) = K0s;                                       % round-4 P1-B: the final refresh runs
    % on the newly mixed Sigma(1), so its closed K0 differs from the seed -- KEEP it, and
    % report the SAME value the returned Gstat/r/D_uni were computed with.
    ctol = getf(eso, 'resid_tol', 1e-10);             % documented closure tolerance (meV^-1)
    ok = ok && so.converged && isfinite(so.resid) && so.resid < ctol;
    % round-5 P1-B: the final refresh must ITSELF converge and close -- an unconverged
    % refresh makes this node not-ok (callers then mark node_failed), never silent export.
    rk = so.r;  S0k = Sigma(1);  K0k = K0s;  Dk = so.D_uni;  Gbk = G0bare0;
    end
end
```

Note the ξ formula takes `lam(1:2)` — `invz_emt_static_ordered`/`invz_gstat_ordered` require only λ₁, λ₂ (the loop's `lam` is the 3-vector for `invz_sigma_ordered`).

- [ ] **Step 4: Run the new tests — 7 fast must pass; run the near-boundary test once with `INVZ_SLOW=1` (must pass); then the fast suite: 158 passed / 0 failed / 20 incomplete** (151 + 7 fast: bare hook, onset, continuity, grid convergence, Gate 6b, path-derivative gate, failure contract; the near-boundary test registers Incomplete without the env var).

- [ ] **Step 5: Commit**

```bash
git add invz_projected/invz_hmf_ordered.m invz_projected/tests/test_invz_hmf_ordered.m
git commit -m "feat(invz): applied-field/H_MF self-consistency + spontaneous root (J 2.31-2.33) -- fixed-field nodes, static-closure integrand, geometric grid with refined root and grid-convergence gate (stage-2 task 3)"
```

---

### Task 4: `'jensen'` mode in `invz_solve_point_ordered`

**Files:**
- Modify: `invz_projected/invz_solve_point_ordered.m` (mode gate; bare path byte-identical)
- Test: `invz_projected/tests/test_invz_ordered_jensen.m` (new file; or appended to `test_invz_ordered_phase.m` if its setup fits — implementer picks and states which)

**Interfaces:**
- Consumes: Task 3's `invz_hmf_ordered` (full-opts pass-through).
- Produces (Task 5 relies on): `opts.ordered_mode = 'bare'` (default) `| 'jensen'`. In jensen mode `pt` additionally carries `pt.ordered_mode = 'jensen'`, `pt.hmf` (the refined root, meV), `pt.D_uni` (the ordered static inverse response at the FINAL state — the pole observable), `pt.hmf_prof` (the Task-3 `prof`), and `pt.m0 = si.Jexp(3)` at `hz_fixed = pt.hmf`; `pt.is_ordered` is gated on ROOT EXISTENCE (`isfinite(hmf_star)`), NOT on `m_tol` (P1-4) — a paramagnetic early return (existing struct shape) when the root is absent. Guards: `invz:orderedMode` on `forced_moment` OR `|Bx(3)| > bz_tol` (P1-6). The final state's Σ↔EMT loop includes the SAME static-sector closure insertion as Task 3's nodes (statement order in Task 3's Interfaces bullet, full-sector weights, threaded `emt_static` opts), and **before packing diagnostics the ω=0 slot of the medium propagator is replaced by the closed elastic value: `med.G(1) = Gstat;` (P1-D)** — so `pt.G`, the sum rule, and every downstream equal-time/energy quantity see the elastic static function, and `pt.D_uni = 1 + (J0eff − K(1))·pt.G(1)` is evaluated at the final state (not copied from the profile). **Contract (P2-G):** `pt.crit` KEEPS its historical ordinary-Dyson definition (`1 + Sigma0 − J0eff·chi0cc0`) as a legacy diagnostic and is NOT the ordered pole mass below the boundary — `pt.D_uni` is; the solver docstring must state this explicitly so callers cannot conflate them.

- [ ] **Step 1: Write the failing tests**

```matlab
function test_jensen_mode_root_and_suppression(testCase)
% jensen mode at 2.85 T: ordered, converged, moment suppressed below bare, hmf below
% the bare fixed point, si.hz == hmf (fixed-field contract); at 3.30 T: PM early return.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
pt = invz_solve_point_ordered(ion, T, [2.85 0 0], Jnu, o);
verifyTrue(testCase, pt.is_ordered && pt.converged);
verifyEqual(testCase, pt.ordered_mode, 'jensen');
verifyEqual(testCase, pt.si.hz, pt.hmf, 'AbsTol', 1e-12);      % P0-1 contract on the final state
ptb = invz_solve_point_ordered(ion, T, [2.85 0 0], Jnu, rmfield(o, 'ordered_mode'));
verifyLessThan(testCase, pt.m0, ptb.m0);                       % fluctuation-suppressed moment
verifyLessThan(testCase, pt.hmf, ptb.si.hz);                   % root below the bare fixed point
verifyTrue(testCase, isfinite(pt.D_uni));                      % pole observable exposed
% P1-D (round-3 DISCRIMINATIVE form): the closure and D_uni identities alone cannot
% distinguish the elastic Gstat from the ordinary-Dyson value (both satisfy them), so:
% (a) pt.G(1) must equal the profile's closed Gstat_star;
% (b) pt.G(1) must equal an INDEPENDENT recomputation of Gstat from the final state;
% (c) pt.G(1) must DIFFER from the ordinary-Dyson static value at this finite moment.
verifyEqual(testCase, pt.G(1), pt.hmf_prof.Gstat_star, 'RelTol', 1e-8);
C = invz_const();  beta = 1/(C.kB*T);
c0e = invz_chi0z(pt.si, T, 0, struct('elastic', true));
c0i = invz_chi0z(pt.si, T, 0, struct('elastic', false));
X = real(c0e(:, :, 1));
fb = X(3,1) * (ion.Jxx0 / (1 - ion.Jxx0*X(1,1))) * X(1,3);     % legacy_x: 1x1 block (SS4a)
G0bare0 = -(X(3,3) + fb);  G0inel0 = -real(c0i(3,3,1));
Gs_ind = invz_gstat_ordered(pt.tl, pt.lambda(1:2), pt.K(1), pt.Sigma(1), beta, ...
                            G0inel0, G0bare0 - G0inel0);
verifyEqual(testCase, pt.G(1), Gs_ind, 'RelTol', 1e-8);        % independent recomputation
G_dyson = (-X(3,3)) / (1 + pt.Sigma(1) + pt.K(1)*(-X(3,3)));   % the OLD ordinary value
verifyGreaterThan(testCase, abs(pt.G(1) - G_dyson)/abs(G_dyson), 1e-6);  % must differ
verifyEqual(testCase, pt.D_uni, 1 + (6.4e-3 - pt.K(1))*pt.G(1), 'RelTol', 1e-10);
% round-5 P1-B: the EXPORTED tuple must satisfy the EMT lattice closure -- the formula
% consistency above proves what Gstat(K) is, this proves the exported K actually closes
Gq = pt.G(1) ./ (1 + (Jnu - pt.K(1)).*pt.G(1));
verifyEqual(testCase, mean(Gq), pt.G(1), 'RelTol', 1e-8);
verifyLessThan(testCase, pt.final_resid, 1e-8);                % outer-tol revalidation held
pth = invz_solve_point_ordered(ion, T, [3.30 0 0], Jnu, o);
verifyFalse(testCase, pth.is_ordered);                         % no ordered state above Bc_1z
end

function test_jensen_mode_guards(testCase)
ion = invz_ion();  Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
o1 = base;  o1.forced_moment = true;
verifyError(testCase, @() invz_solve_point_ordered(ion, 0.31, [2 0 0.5], Jnu, o1), ...
            'invz:orderedMode');                               % forced_moment forbidden
verifyError(testCase, @() invz_solve_point_ordered(ion, 0.31, [2 0 0.5], Jnu, base), ...
            'invz:orderedMode');                               % longitudinal field forbidden (P1-6)
end
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement the mode gate**

In `invz_solve_point_ordered.m`, immediately BEFORE the `si = invz_single_ion(ion, T, Bx, siopts);` call (so `siopts` is already built, opts parsed, and the ODD `-d` shift of `J0eff` applied — the jensen branch must see the shifted coupling), insert:

```matlab
omode = getf(opts, 'ordered_mode', 'bare');
if ~any(strcmp(omode, {'bare', 'jensen'}))
    error('invz:orderedMode', 'ordered_mode must be ''bare'' or ''jensen''.');
end
if strcmp(omode, 'jensen')
    if fmom || abs(Bx(3)) > getf(opts, 'bz_tol', 1e-9)
        error('invz:orderedMode', 'ordered_mode ''jensen'' is transverse/spontaneous only.');
    end
    hopts = opts;                                    % FULL numerical context (P1-6) ...
    hopts.J0eff = J0eff;                             % ... with the ODD-shifted coupling
    for f = {'ordered_mode', 'forced_moment'}        % ... and mode fields stripped
        if isfield(hopts, f{1}), hopts = rmfield(hopts, f{1}); end
    end
    [hstar, hprof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, hopts);
    if ~isfinite(hstar)
        si = invz_single_ion(ion, T, Bx, struct('hyp', hyp, 'hz_fixed', 0, ...
                                                'Jxx0', Jxx0, 'transverse_mf', tmf));
        pt = early_return(0, si, 'none');            % paramagnetic: PM leg owns this field
        pt.ordered_mode = omode;  pt.hmf_status = hprof.status;
        if strcmp(hprof.status, 'unresolved')
            pt.converged = false;                    % round-3 P0-3: ordering was PREDICTED
        end                                          % but unbracketed -- NOT a PM verdict;
        return;                                      % the map masks this column
    end
    siopts.hz_fixed = hstar;                         % impose the jensen molecular field ...
    siopts.order = false;                            % ... WITHOUT the ordering update (P0-1)
end
```
then, in the jensen case only: gate `pt.is_ordered` on the root (bypass the `m_tol` test — P1-4); make the final Σ↔EMT loop include the static-sector closure insertion (same statement order as Task 3's `eval_node`: full-sector weights from the final `si` with the mode-switched chain rule — round-4 P1-A; threaded `eso = getf(opts,'emt_static',struct())` — P1-F); the final post-loop static refresh KEEPS its newly closed `K0`, writes it to `K(1)` before packing, and is GATED on its own `so.converged` and `so.resid < resid_tol` — folded into `pt.converged` (round-4 P1-B / round-5 P1-B). After the refresh, perform a residual-only revalidation: recompute `lam_check = invz_lambdas(K, g, wts, beta, [1 2 3])` and `Sigma_check = invz_sigma_ordered(tl, lam_check, K, g, beta)`, store `pt.final_resid = max(abs(Sigma_check.Sigma − Sigma))`, and require `pt.final_resid < tolo` for `pt.converged` (round-5 P2: the exported tuple is thereby closure-consistent TO THE STATED OUTER TOLERANCE — the docstring must say exactly that, not "exactly self-consistent"); the discriminative test below verifies the exported set because its independent recomputation uses the exported `pt.K(1)`; and BEFORE the existing diagnostics packing (`pt.G = med.G` and the sum rule) write the closed elastic value into the static slot (P1-D):
```matlab
med.G(1) = Gstat;                                    % elastic static function, not ordinary Dyson
```
then after `pt` assembly add:
```matlab
if strcmp(omode, 'jensen')
    pt.ordered_mode = omode;  pt.hmf = hstar;  pt.hmf_prof = hprof;
    pt.D_uni = 1 + (J0eff - K(1))*med.G(1);          % pole observable AT THE FINAL STATE
    % CONTRACT (P2-G): pt.crit keeps its historical ordinary-Dyson definition and is NOT
    % the ordered pole mass below the boundary -- pt.D_uni is (docstring must say so).
    if abs(pt.si.hz - hstar) > 1e-12
        error('invz:hzFixed', 'jensen final state did not hold hmf: %.6g vs %.6g', pt.si.hz, hstar);
    end
end
```
The `'bare'` path must not gain, lose, or reorder a single line.

- [ ] **Step 4: Run the new tests (2/2) and the fast suite: expected 160 passed / 0 failed / 20 incomplete** (158 + 2). Any drift in pre-existing ordered-solver tests means the bare path was NOT left byte-identical — fix before proceeding.

- [ ] **Step 5: Commit**

```bash
git add invz_projected/invz_solve_point_ordered.m invz_projected/tests/<test file used>
git commit -m "feat(invz): 'jensen' ordered mode -- H_MF root imposed via order-free hz_fixed, static-closure-consistent final state, root-existence acceptance, transverse guards (stage-2 task 4)"
```

---

### Task 5: Map integration, pole-closure regression (diagnosis regression 2, un-waived), docs

**Files:**
- Modify: `invz_projected/invz_spectra_map.m` (`one_field` else-branch + docstring + pack: new `S.m_1z`, `S.D_ord`)
- Modify: `invz_projected/invz_check_solve_opts.m` (reserve `ordered_mode` — P1-6)
- Modify: `invz_projected/invz_run_spectra.m` (panel label)
- Modify: `invz_projected/README.html` (callout sentence; stash choreography if dirty)
- Modify: `invzp_QCP_diagnosis.md` (Stage-2 status; regression 2 un-waived)
- Test: `invz_projected/tests/test_invz_spectra_map.m` (add the reserved-field assert) + new `invz_projected/tests/test_invz_qcp_closure.m` (INVZ_SLOW)

**Interfaces:**
- Consumes: Task 4's jensen mode.
- Produces: `one_field` signature grows to `[chiz, chirpa, Sigma0, phase, phase_1z, crit_pm, m_1z, D_ord]` (new parfor-sliced outputs, NaN except on jensen-ordered columns); `S.m_1z`, `S.D_ord [1 x nB]`; map option `opts.ordered_1z = 'jensen'` (default) `| 'bare'` (Stage-1 diagnostic escape hatch), threaded to `one_field` as an argument; `invz_check_solve_opts` errors (`invz:solveOpts`) on `solve_opts.ordered_mode` — the mode is driver-owned, so the auto/overlay leg can never be flipped (P1-6). `phase_1z = 1` under `'jensen'` means the CONSISTENT ordered 1/z state; jensen failure with an invalid PM probe leaves `phase_1z = 0` (masked), never a silent bare fallback.

- [ ] **Step 1: Reserve `ordered_mode` and assert it**

In `invz_check_solve_opts.m`, add `'ordered_mode'` to the reserved-fields list (same error id `invz:solveOpts`, same message pattern as `J0eff/Jxx0/hyp`), and update its header comment: "driver-owned: the map's 1/z leg sets it internally; the auto/overlay leg must never see it." Add to `test_invz_spectra_map.m` (in `test_split_phase_diagnostics`, after the existing asserts):

```matlab
% ordered_mode is driver-owned (stage-2 P1-6): a caller-supplied value must error
verifyError(testCase, @() invz_spectra_map(ion, T, fields, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, ...
           'solve_opts', struct('ordered_mode', 'jensen'))), 'invz:solveOpts');
```

- [ ] **Step 2: Replace the `one_field` moment-form else-branch**

Current branch (the Stage-1 diagnostic; retained verbatim under `'bare'`):

```matlab
    else                                     % --- below Bc_1z: ordered 1/z, DIAGNOSTIC ---
        phase_1z = 1;  Sigma0 = pt.Sigma0;
        o  = invz_chi_realaxis(ion, T, B, pt, w, copts);
        chiz = imag(o.chi_cc_q(1, :)).';
        pt0 = invz_zero_sigma_overlay(pt);
        c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;
        o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
        chirpa = imag(o0.chi_cc_q(1, :)).';
    end
```

New (with `o1z` the threaded map option):

```matlab
    else                                     % --- below Bc_1z: ordered 1/z leg ---
        if strcmp(o1z, 'jensen') && B(3) == 0   % Stage 2: consistent ordered state (H_MF).
            % TRANSVERSE ONLY (as-built amendment, task 5): jensen mode hard-errors on a
            % longitudinal tilt (plan SS8 -- the tilted route keeps the bare mode, a
            % rounded crossover with no sharp Bc), so tilted fields fall through to the
            % byte-identical bare diagnostic branch below: out of scope, not "failed".
            so2 = sopts;  so2.ordered_mode = 'jensen';
            ptj = [];
            try
                ptj = invz_solve_point_ordered(ion, T, B, Jnu, so2);
            catch err
                if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
            end
            if ~isempty(ptj) && ptj.is_ordered && ptj.converged && isfinite(ptj.Sigma0)
                phase_1z = 1;  Sigma0 = ptj.Sigma0;  m_1z = ptj.m0;  D_ord = ptj.D_uni;
                oj = invz_chi_realaxis(ion, T, B, ptj, w, copts);   % jensen si differs from
                chiz = imag(oj.chi_cc_q(1, :)).';                   % the auto pt's -- no sharing
            end                              % else phase_1z stays 0 -> chiz column masked
            pt0 = invz_zero_sigma_overlay(pt);                      % overlay: UNCHANGED auto state
            c0opts = copts;  c0opts.npass = 1;
            o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        else                                 % 'bare': Stage-1 diagnostic escape hatch
            phase_1z = 1;  Sigma0 = pt.Sigma0;
            o  = invz_chi_realaxis(ion, T, B, pt, w, copts);
            chiz = imag(o.chi_cc_q(1, :)).';
            pt0 = invz_zero_sigma_overlay(pt);
            c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;
            o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
            chirpa = imag(o0.chi_cc_q(1, :)).';
        end
    end
```

Plus: `m_1z = NaN; D_ord = NaN;` in the initializer; `o1z` appended to `one_field`'s arguments, read in the map body as `o1z = getf(opts, 'ordered_1z', 'jensen');` (validate against `{'jensen','bare'}`, error `invz:ordered1z`); sliced `m1zM`/`DordM` init + parfor call site + pack (`S.m_1z`, `S.D_ord`, `S.ordered_1z = o1z;`). Docstring: under `'jensen'`, `phase_1z = 1` is the consistent ordered 1/z state (J 2.28–2.33; diagnosis Stage 2 delivered), the FM/PM branches meet at `Bc_1z`, and `S.D_ord` is the ordered static inverse response (→ 0 at `Bc_1z` from below, the pole observable); `'bare'` reproduces the Stage-1 diagnostic.

- [ ] **Step 3: Driver + README + diagnosis text updates**

- `invz_run_spectra.m` panel title: `'1/z (own phase; FM side diagnostic below B_c^{1/z}), T = %.2f K%s'` → `'1/z (own phase; Stage-2 ordered below B_c^{1/z}), T = %.2f K%s'`.
- `README.html` callout: replace the fragment `the bare-MF ordered curve below \(B_c^{1/z}\) as a flagged <em>diagnostic</em> — its ordered-side pole closure at \(B_c^{1/z}\) is a Stage-2 requirement, not a Stage-1 claim` with `the Jensen-consistent ordered 1/z state below \(B_c^{1/z}\) (elastic static sector J 2.28–2.29 + static-closure medium + H<sub>MF</sub> relation J 2.31–2.33), whose moment and static inverse response <code>S.D_ord</code> vanish at \(B_c^{1/z}\) so the FM and PM 1/z branches close at the same field; <code>ordered_1z = 'bare'</code> reproduces the retired Stage-1 diagnostic` (stash choreography if the file carries user edits).
- `invzp_QCP_diagnosis.md`: in "Required regression tests", replace regression 2's parenthetical `(Blocked on Stage 2: whichever route is chosen must include J 2.28--2.29 and J 2.31--2.34; until then the ordered-side 1/z branch is labelled diagnostic and this test is waived, not weakened.)` with `(Delivered by the Stage-2 plan on the projected path: enforced by test_invz_qcp_closure, INVZ_SLOW-gated, via the static inverse response -- pole-based, per Stage 3's argmax warning.)`; in Stage 2, prepend `STATUS: implemented on the projected path (docs/superpowers/plans/2026-07-22-invzp-stage2-ordered-thermodynamics.md).` to the section body.

- [ ] **Step 4: Write the pole-closure regression (INVZ_SLOW; P1-5)**

Create `invz_projected/tests/test_invz_qcp_closure.m`:

```matlab
function tests = test_invz_qcp_closure
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_fm_pm_poles_close_at_bc1z(testCase)
% Diagnosis regression 2, un-waived -- POLE-BASED (P1-5): the ordered static inverse
% response D_ord -> 0 from below and the PM mass crit -> 0 from above must vanish at
% the SAME field, the moment vanishing continuously; plus broadening/mesh refinement
% stability of the near-boundary spectrum. Minutes-scale -> INVZ_SLOW-gated.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (jensen sweep).');
ion = invz_ion();  T = 0.31;
fields = linspace(2.70, 3.50, 17);                   % brackets Bc_1z ~ 3.0 (synthetic)
w = (0.005:0.005:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
S = invz_spectra_map(ion, T, fields, w, ...
        struct('Jnu', Jnu, 'info', info, 'verbose', false));

fm = S.phase_1z == 1;  pm = S.phase_1z == 2;
verifyGreaterThanOrEqual(testCase, nnz(fm), 3);
verifyGreaterThanOrEqual(testCase, nnz(pm), 3);
step = fields(2) - fields(1);
% (i) the label flip agrees with the PM mass zero-crossing within one grid step
cp = S.crit_pm;  okp = isfinite(cp) & pm;
Bc_crit = interp1(cp(okp), fields(okp), 0, 'linear', 'extrap');
verifyLessThan(testCase, abs(S.Bc_1z - Bc_crit), 2*step);
% (ii) POLE closure from BOTH sides: D_ord (ordered static inverse response) monotone
% decreasing toward the boundary and extrapolating to zero at the same field as crit
okd = isfinite(S.D_ord) & fm;
verifyGreaterThanOrEqual(testCase, nnz(okd), 3);
Dd = S.D_ord(okd);  Bd = fields(okd);
verifyTrue(testCase, all(Dd > 0) && all(diff(Dd) < 0));    % positive, closing from below
B_dn = interp1(Dd, Bd, 0, 'linear', 'extrap');
verifyLessThan(testCase, abs(B_dn - Bc_crit), 2*step);
verifyLessThan(testCase, abs(B_dn - S.Bc_1z), 2*step);
% (iii) the ordered moment vanishes continuously toward the boundary
mo = S.m_1z(fm);
verifyTrue(testCase, all(diff(mo) < 0) && mo(end) < 0.5*mo(1));
% (iv) broadening/mesh refinement via a LOCAL BRANCH TRACKER on the inverse response
% (round-3 P1-6): the soft branch is the LOWEST-frequency LOCAL minimum of
% |1 - J0eff*chi~0(w)|, chi~0 = chi0cc_w./(1+Sigma_w) -- a branch identity, not a
% global argmin (which can hop between electronuclear satellites under refinement).
% The refined minimum must lie in the SAME basin and shift only within tolerance.
% NOTE (documented non-identity): at finite m the real-axis pipeline excludes the
% elastic sector, so the w -> 0 limit of this inverse response equals D_uni only in
% the m -> 0 boundary limit; D_ord/crit above remain the load-bearing closure gates --
% this subtest checks refinement stability of the finite-frequency pole only.
Blast = fields(find(fm, 1, 'last'));
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
ptj = invz_solve_point_ordered(ion, T, [Blast 0 0], Jnu, o);
w2 = (0.0025:0.0025:0.6).';
c1 = struct('Jsel', 6.4e-3, 'eta', 5e-3,   'Jxx0', ion.Jxx0, 'Jshape', 0, 'hyp', true, ...
            'transverse_mf', 'legacy_x', 'si', ptj.si);
c2 = c1;  c2.eta = 2.5e-3;
o1 = invz_chi_realaxis(ion, T, [Blast 0 0], ptj, w,  c1);
o2 = invz_chi_realaxis(ion, T, [Blast 0 0], ptj, w2, c2);
dinv1 = abs(1 - 6.4e-3 * o1.chi0cc_w(:) ./ (1 + o1.Sigma_w(:)));
dinv2 = abs(1 - 6.4e-3 * o2.chi0cc_w(:) ./ (1 + o2.Sigma_w(:)));
i1 = find(islocalmin(dinv1), 1, 'first');               % lowest-frequency LOCAL minimum
i2 = find(islocalmin(dinv2), 1, 'first');
verifyTrue(testCase, ~isempty(i1) && ~isempty(i2));
% same-basin check: the refined minimum must fall inside the coarse minimum's basin
% (between the coarse grid's neighboring local maxima around i1)
lm = find(islocalmax(dinv1));
blo = max([w(lm(lm < i1)); 0]);  bhi = min([w(lm(lm > i1)); w(end)]);
verifyTrue(testCase, w2(i2) > blo && w2(i2) < bhi, 'refined minimum left the basin');
verifyLessThan(testCase, abs(w(i1) - w2(i2)), max(4*(w(2)-w(1)), 0.1*max(w(i1), w2(i2))));
% (v) sweep-direction invariance: reversed field order reproduces labels and moments
S2 = invz_spectra_map(ion, T, fliplr(fields), w, ...
        struct('Jnu', Jnu, 'info', info, 'verbose', false));
verifyEqual(testCase, fliplr(S2.phase_1z), S.phase_1z);
verifyEqual(testCase, fliplr(S2.m_1z), S.m_1z, 'RelTol', 1e-6, 'AbsTol', 1e-9);
end
```
Tolerances are grid-limited by design; the knobs to tighten are `fields` density, `nH`, and the refinement pair in (iv). `Epeak` fits are deliberately NOT used as gates (they are argmax quantities — diagnosis Stage-3 warning); the spectra remain available for visual inspection.

- [ ] **Step 5: Run** the closure test once with `INVZ_SLOW=1` (must pass), the spectra-map test file, then the fast suite: expected **160 passed / 0 failed / 21 incomplete** (the reserved-field assert lives inside an existing test function so the passed count is unchanged; the new slow test registers Incomplete). Verify the exact count against the post-Task-4 baseline and record it.

- [ ] **Step 6: Commit**

```bash
git add invz_projected/invz_spectra_map.m invz_projected/invz_check_solve_opts.m invz_projected/invz_run_spectra.m invz_projected/README.html invzp_QCP_diagnosis.md invz_projected/tests/test_invz_qcp_closure.m invz_projected/tests/test_invz_spectra_map.m
git commit -m "feat(invz): stage-2 wired into the 1/z ordered leg -- pole-based FM/PM closure at Bc_1z (regression 2 un-waived via D_ord), ordered_mode reserved, 'bare' escape hatch (stage-2 task 5)"
```

---

### Task 6: Free-energy consistency validation (J 2.34) — EXECUTED, FAILED AS DESIGNED, SUPERSEDED

> **EXECUTION RECORD — do NOT execute these acceptance instructions.** This task ran (commit `bf9621a`): route B validated internally; the production route A diverged for the domain reasons diagnosed in §7a/§7b (numbers preserved in `task-6-report.md` and the diagnosis). **Task 6b is the only live acceptance contract** (second §7 review, P1). The text below is retained verbatim for provenance.

**Files:**
- Create: `invz_projected/invz_deltaF_ordered.m`
- Test: `invz_projected/tests/test_invz_deltaF.m` (one fast + one INVZ_SLOW localfunction)

**Interfaces:**
- Consumes: Task 3's `prof` (route A is a post-processing integral of it); PM solves (`invz_solve_point`) on a temperature grid for route B, whose δU uses `pt.alpha`, `pt.lambda(1)`, `pt.K`, `pt.Sigma`, `pt.tl`, plus `G0(iωn)` recomputed from `pt.si` on the SOLVER'S OWN Matsubara grid, reconstructed locally per temperature (round-2 P0-B — the test defines `wnk`, `wtsk`, `betak` itself by calling the same grid constructor `invz_solve_point` uses; the implementer greps the solver for its grid-building call — if it is a helper like `invz_matsubara(T, Ecut)`, reuse it; if inline, mirror the block verbatim — and uses `betak` from it, NOT a separate `kB` literal).
- Produces: `[dF, out] = invz_deltaF_ordered(ion, T, Bx, Jnu_flat, opts)` — route A; `out.dF_routeA` (per site, meV), `out.tail_est` (saturation-truncation ESTIMATE `|δh(end)|·(M0 − m(end))`, `M0 = ion.J` exactly — `max⟨Jz⟩ = J` — unless the repo exports a Jz operator constructor, in which case use `max(abs(eig(Jz)))` and state which; **P2-G contract: this is an estimate, not a proven bound, unless `|δh|` is verified monotone/plateaued beyond the last node — the slow test therefore ALSO requires route-A convergence under a larger `hmax_fac`**), `out.prof`.

- [ ] **Step 1: Fast test — bare limit vanishes**

```matlab
function test_bare_limit_deltaF_zero(testCase)
% Route A with force_bare: dh = h0 - hmf = 0 identically, so BOTH dF and the
% truncation tail estimate are exactly zero (P0-3: the round-1 draft asserted
% tail > 0, which contradicts dh == 0).
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[dF, out] = invz_deltaF_ordered(ion, T, [2.85 0 0], Jnu, ...
    struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'force_bare', true));
verifyEqual(testCase, dF, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, out.tail_est, 0, 'AbsTol', 1e-15);
end
```

- [ ] **Step 2: Implement route A**

```matlab
function [dF, out] = invz_deltaF_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_DELTAF_ORDERED Free-energy correction dF(m=0), route A (framework SS9.4, J 2.34):
%   dF = -int_0^{M0} (h0(m') - hmf(m')) d<Jz>'
% evaluated on the invz_hmf_ordered profile as -trapz(m, h0 - hgrid), with the profile
% extended toward saturation (hmax_fac default 4: fluctuations quench as the splitting
% grows, dh -> 0). out.tail_est = |dh(end)|*(M0 - m(end)) estimates the cut tail and must
% be reported alongside dF whenever dF is quoted. Validation-only -- no production path
% depends on this function.
if nargin < 5, opts = struct(); end
o = opts;  o.hmax_fac = getf(opts, 'hmax_fac', 4);  o.nH = getf(opts, 'nH', 65);
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, o);
% round-4 P1-C: integrate ONLY trusted profiles -- an 'unresolved' or 'node_failed'
% profile can carry a nonempty grid and finite arrays, but its contract forbids use.
trusted = strcmp(prof.status, 'ok') && ~isempty(prof.hgrid) && all(prof.node_conv);
if ~trusted
    dF = NaN;  out = struct('dF_routeA', NaN, 'tail_est', NaN, 'prof', prof);  return;
end
dh = prof.h0 - prof.hgrid;
dF = -trapz(prof.m, dh);
M0 = ion.J;                                          % max <Jz> = J exactly
out = struct('dF_routeA', dF, 'tail_est', abs(dh(end))*(M0 - prof.m(end)), 'prof', prof);
end
```
(`tail_est` is an ESTIMATE — P2-G: it bounds the cut tail only if `|δh|` does not grow past the last node; the slow test's `hmax_fac` doubling is the convergence check that makes quoting it honest.)

- [ ] **Step 3: Slow test — the two-route check (sign corrected, tails explicit — P0-3)**

Append to `test_invz_deltaF.m`:

```matlab
function test_two_routes_agree(testCase)
% Framework SS9.4: route A (field dependence of Sigma) vs route B (frequency structure,
% J 2.22 temperature-integrated) must agree -- the stringent global check.
%   dF(T) = + T * int_T^inf dU(T')/T'^2 dT'      (sign derived in plan SS7, P0-3)
% Tolerance 10% relative with BOTH truncation tails required small first; tighten grids
% before tolerances. FAILURE BEYOND TOLERANCE IS A BLOCKED ESCALATION (SS9).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (T-grid sweep).');
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[dFA, outA] = invz_deltaF_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfinite(dFA));
verifyLessThan(testCase, outA.tail_est, 0.05*abs(dFA));      % route-A tail must be small
% P2-G: tail_est is only honest if route A is CONVERGED in the saturation cutoff --
% doubling hmax_fac must not move dF beyond 5% relative.
[dFA2, ~] = invz_deltaF_ordered(ion, T, Bx, Jnu, setfield(o, 'hmax_fac', 8)); %#ok<SFLD>
verifyEqual(testCase, dFA, dFA2, 'RelTol', 0.05);
% Route B on a geometric T grid; dU per J 2.22 with G reconstructed via J 2.30, on the
% SOLVER'S OWN Matsubara grid rebuilt per temperature (round-2 P0-B, made concrete in
% round 3: invz_matsubara exists -- verify the solver calls it, then reuse).
Tg = T * (1.35.^(0:17));  dU = nan(size(Tg));        % 18 points: last 3 are the Tmax check
for k = 1:numel(Tg)
    pt = invz_solve_point(ion, Tg(k), Bx, Jnu, o);
    if ~pt.converged, continue; end
    [wnk, wtsk, betak] = invz_matsubara(Tg(k), getf(o, 'Ecut', 40));
    c0k = invz_chi0z(pt.si, Tg(k), 1i*wnk, struct('elastic', true));
    G0k = -real(squeeze(c0k(3,3,:)));
    Gk  = G0k ./ (1 + pt.Sigma + pt.K .* G0k);       % J 2.30
    tlk = pt.tl;
    dU(k) = 0.5*( pt.alpha*tlk.n01*tlk.Delta/(1 + pt.alpha) - tlk.M2*pt.lambda(1) ...
                  + real(sum(wtsk .* pt.K .* (Gk - G0k)))/betak );
end
fin = isfinite(dU);
verifyGreaterThanOrEqual(testCase, nnz(fin), 15);
Tf = Tg(fin);  Uf = dU(fin);
n15 = nnz(Tf <= T*1.35^14);                          % the shorter grid's reach
dFB15 = + T * trapz(Tf(1:n15), Uf(1:n15) ./ Tf(1:n15).^2);   % CORRECTED SIGN (P0-3)
dFB   = + T * trapz(Tf, Uf ./ Tf.^2);
% round-3 P2-7: the tail is an EMPIRICAL estimate, not a proven bound -- so require
% Tmax convergence explicitly (extending Tmax by 1.35^3 moves dFB < 5%), and report
% tail_est_B alongside.
verifyEqual(testCase, dFB15, dFB, 'RelTol', 0.05);   % Tmax-converged
tail_est_B = + T * abs(Uf(end)) / Tf(end);           % empirical remainder estimate
verifyLessThan(testCase, tail_est_B, 0.05*max(abs(dFA), abs(dFB)));
verifyEqual(testCase, dFA, dFB, 'RelTol', 0.10, ...
    sprintf('two-route dF mismatch (A=%.4g, B=%.4g, tailA~%.2g, tailB~%.2g) -- BLOCKED escalation', ...
            dFA, dFB, outA.tail_est, tail_est_B));
end
```
(Round-4 P2-E: the route-B code is fully concrete — `invz_matsubara` is the solver's own grid constructor; the implementer verifies that with one grep before relying on it. If the tail assertions fail, EXTEND the grids — `hmax_fac`/`nH` for route A, more/higher `Tg` points for route B — never weaken the 5% tail bounds.)

- [ ] **Step 4: Run** the fast δF test, the slow two-route test with `INVZ_SLOW=1` (MUST pass — this is the algebra's global verdict), then the fast suite: expected **161 passed / 0 failed / 22 incomplete** (160 + 1 fast; the slow test registers Incomplete).

- [ ] **Step 5: Commit**

```bash
git add invz_projected/invz_deltaF_ordered.m invz_projected/tests/test_invz_deltaF.m
git commit -m "feat(invz): free-energy consistency validation (J 2.34) -- corrected-sign temperature route, local Matsubara reconstruction, explicit truncation tails and hmax convergence, two-route agreement gate (stage-2 task 6)"
```

---

### Task 6b: Route-A resolution per §7b — closed-model J 2.34 validation + production relabel

**Files:**
- Modify: `invz_projected/invz_deltaF_ordered.m` (relabel to the partial-diagnostic contract)
- Modify: `invz_projected/invz_hmf_ordered.m` + `invz_projected/tests/test_invz_hmf_ordered.m` (the `hmax_abs` exact-override option and its test — third §7 review P1)
- Modify: `invz_projected/tests/test_invz_deltaF.m` (update the fast test; REPLACE the slow production two-route test with the FAST closed-model two-route test)
- Modify: `invzp_QCP_diagnosis.md` (scope-limitation note preserving the failed numbers)

**Interfaces:**
- Consumes: the production functions §7b lists; `invz_matsubara(T, Ecut)`; Task 1's fixture pattern.
- Produces: `invz_deltaF_ordered` returns `out.dF_partial` (REPLACING `out.dF_routeA`; the first output `dF` keeps its value but the docstring states the partial-diagnostic contract: finite-domain line integral over the hybrid profile, cutoff-dependent — `out.hmax_used` reported — NEVER `δF(m=0)`, never compared to route B as a gate, trust gate on `prof.status` unchanged); the closed-model test is self-contained in the test file.

- [ ] **Step 1: Relabel the production diagnostic (second-review P1 contract fixes included)**

In `invz_deltaF_ordered.m`: rename `out.dF_routeA` → `out.dF_partial`; RENAME `out.tail_est` → `out.endpoint_dh = abs(dh(end))` (a NON-CLAIMING endpoint diagnostic — the old `|δh(end)|·(ion.J − m(end))` product implied a sanctioned continuation to `ion.J`, which §7b rules out; do not describe any field as a bound or saturation remainder); add `out.hmax_used = max(prof.hgrid)` with `NaN` on empty/failed profiles (never `max` of an empty grid); add opts `hmax_abs` implemented IN `invz_hmf_ordered` (third §7 review P1 — the file is now in scope): a positive finite value is an EXACT OVERRIDE, `hmax = hmax_abs` (never a `min` cap, which would not give a common cutoff); nonpositive/nonfinite values error (`invz:hmfOpts`); the bare-order early-return check is retained separately from the ceiling selection; add a test in `test_invz_hmf_ordered.m` asserting two states with DIFFERENT bare `hz` values report identical `max(prof.hgrid)` under the same `hmax_abs`. Keep the profile trust/status gate. Rewrite the header: "PARTIAL HYBRID DIAGNOSTIC (plan §7b): finite-domain Legendre line integral over the H_MF profile. NOT δF(m=0) — the hybrid's saturation-normalized absolute free energy is outside its validated domain (route-A divergence recorded in the task-6 execution record); cross-state comparison is meaningful only at a COMMON cutoff (use opts.hmax_abs). The J 2.34 two-route identity is validated in the closed 2×2 model (test_invz_deltaF)." Update the fast bare test accordingly (`dF == 0`, `endpoint_dh == 0`).

- [ ] **Step 2: Write the closed-model two-route test (replaces the slow test)**

Replace `test_two_routes_agree` in `test_invz_deltaF.m` with (FAST — 2×2 model, no CF diagonalization):

```matlab
function tl = closed_tl(h, D0, M0op, beta)
% FROZEN closed 2x2 (SS7b): H = diag(+/-D0/2) - h*Jz2, Jz2 = [0 M0op; M0op 0] (FIXED
% operator, never re-selected from a CF solve). Purely off-diagonal Jz2: m(0) = 0 (true
% PM anchor) and M_sat = M0op ANALYTICALLY (the Jz2 eigenvalue).
Jz2 = [0, M0op; M0op, 0];
H = [D0/2, 0; 0, -D0/2] - h*Jz2;
[V, E] = eig((H+H')/2, 'vector');  [E, ix] = sort(real(E));  V = V(:, ix);
p = exp(-beta*(E - E(1)));  p = p/sum(p);
Mz = V'*Jz2*V;
tl = struct('m', real(Mz(1,1)), 'M2', abs(Mz(1,2))^2, 'n01', p(1) - p(2), ...
            'Delta', E(2) - E(1), 'g0', 0);
tl.g0 = 2*tl.n01/tl.Delta;
end

function [rk, Mk, S0k, l1g0, ok] = closed_node(h, T, Jnu, J0eff, D0, M0op)
% One closed-model node, COLD-STARTED (third-review P2: node independence explicit --
% no warm-start seed). The ordered Sigma<->EMT loop with TWO-LEVEL static weights: the
% exact SS3-check-(a) parametrization in which J 2.28 is verbatim. Mirrors the
% production round-5 final-tuple contract (third-review P1): every dynamic solve
% requires med.converged; the final static refresh keeps its K0 in K(1); lambdas and
% the ordered Sigma target are RECOMPUTED from that exported K and must close to the
% outer tolerance; the returned diagnostics come from the revalidated final tuple.
[wn, wts, beta] = invz_matsubara(T, 40);
tl = closed_tl(h, D0, M0op, beta);
h0el = beta*(1 - tl.n01^2);
g  = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);  G0(1) = -(tl.M2*tl.g0 + tl.m^2*h0el);      % elastic in the static slot
gi = -tl.M2*tl.g0;  ge = -tl.m^2*h0el;
Sigma = zeros(size(wn));  K = zeros(size(wn));  lam = [0; 0; 0];  K0s = 0;  ok = false;
for outer = 1:200
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    if ~med.converged, break; end
    K = med.K;
    [K0s, ~, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu, K0s, beta, ...
                                           J0eff, gi, ge, struct());
    K(1) = K0s;
    lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
    sg  = invz_sigma_ordered(tl, lam, K, g, beta);
    dS  = max(abs(sg.Sigma - Sigma));
    Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-8 && so.converged, ok = true; break; end
end
% Simultaneous FINAL tuple: refresh the static slot at the final Sigma, keep ITS K0,
% then REVALIDATE lambdas and the ordered Sigma target against the exported K.
[K0s, ~, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu, K0s, beta, ...
                                       J0eff, gi, ge, struct());
K(1) = K0s;
lam_chk = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg_chk  = invz_sigma_ordered(tl, lam_chk, K, g, beta);
ok = ok && so.converged && isfinite(so.resid) && so.resid < 1e-10 ...
        && max(abs(sg_chk.Sigma - Sigma)) < 1e-8;
rk = so.r;  Mk = tl.m*tl.n01;  S0k = Sigma(1);  l1g0 = abs(lam_chk(1)*tl.g0);
end

function [dFA, dgn] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax, nH, testCase)
% Route A on ONE grid, all gates inside -- called three times so the SS7b convergence
% sweeps are EXECUTED assertions, not comments (third-review P1).
hgrid = hmax * (1e-4^(1/(nH-1))).^((nH-1):-1:0);
[r0, ~, ~, ~, ok0] = closed_node(0, T, Jnu, J0eff, D0, M0op);
verifyTrue(testCase, ok0);
[rv, Mv, S0v, Lv] = deal(nan(1, nH));
for k = 1:nH
    [rv(k), Mv(k), S0v(k), Lv(k), okk] = closed_node(hgrid(k), T, Jnu, J0eff, D0, M0op);
    verifyTrue(testCase, okk, sprintf('closed node h = %.4g failed', hgrid(k)));
end
h0 = cumtrapz([0 hgrid], [r0 rv]);  h0 = h0(2:end);
dh = h0 - hgrid;
q = max(1e-4, 0.01*abs(r0 - 1));
verifyTrue(testCase, all(abs(rv(end-4:end) - 1) < q), 'fluctuations not quenched');
verifyTrue(testCase, all(abs(S0v(end-4:end)) < q), 'Sigma(0) not quenched on the tail');
verifyTrue(testCase, all(abs(Lv(end-4:end)) < q), 'lambda1*g0 not quenched on the tail');
verifyLessThan(testCase, M0op - Mv(end), 0.01*M0op);          % analytic saturation reached
dFA = -trapz([0 Mv], [0 dh]);                                  % one moment throughout
dgn = struct('tail_est_A', abs(dh(end))*(M0op - Mv(end)), 'dh_end', dh(end), 'r0', r0);
end

function test_two_routes_closed_model(testCase)
% J 2.34 two-route identity in the CLOSED 2x2 model (SS7b): one moment M = m*n01 serves
% J 2.31 (Gate 1), J 2.33, the work term, and J 2.34's measure -- the conjugacy chain the
% rejected SS7a draft broke. Saturation is VERIFIED against the ANALYTIC M_sat = M0op.
% FAST: seconds of 2x2 work. Disagreement beyond tolerance = BLOCKED escalation (SS9).
T = 0.31;  C = invz_const();  beta = 1/(C.kB*T);
D0 = 0.2;  M0op = 3.0;  J0eff = 6.4e-3;
Jnu = linspace(-4e-3, 4.0e-3, 24).';                  % ZERO-MEAN: enforces Jensen's
verifyLessThan(testCase, abs(mean(Jnu)), 1e-15);      % no-self-site identity J(ii) = 0
% (load-bearing, second SS7 review P0-1 -- a nonzero mean is NOT a valid Jensen lattice:
%  measured consequence was opposite-sign routes)
% ---- Route A with EXECUTED convergence sweeps (third-review P1) --------------------
hmax = 40*D0;  nH = 81;
[dFA,    dA] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax,   nH,     testCase);
[dFA_h2, ~ ] = closed_routeA(T, Jnu, J0eff, D0, M0op, 2*hmax, nH,     testCase);
[dFA_n2, ~ ] = closed_routeA(T, Jnu, J0eff, D0, M0op, hmax,   2*nH-1, testCase);
verifyEqual(testCase, dFA, dFA_h2, 'RelTol', 2e-2);   % hmax-doubling convergence
verifyEqual(testCase, dFA, dFA_n2, 'RelTol', 2e-2);   % grid-doubling convergence
% ---- Route B: temperature integral at h = 0 (PM, m = 0) ----------------------------
% DENSITY-CONVERGED grid (ratio sqrt(1.35)); EVERY node must converge -- silently
% deleting failed temperatures corrupts both convergence checks and can false-pass
% (third-review P1). The DENSE value enters the gate.
Tg = T * (1.35.^(0:0.5:17));  dU = nan(size(Tg));
for k = 1:numel(Tg)
    [wnk, wtsk, betak] = invz_matsubara(Tg(k), 40);
    tlk = closed_tl(0, D0, M0op, betak);
    gk  = real(invz_g(tlk, 1i*wnk));
    G0k = -(tlk.M2*gk);                               % m = 0: no elastic term
    Sk = zeros(size(wnk));  okT = false;
    for outer = 1:200                                 % PM loop: emt + invz_sigma (lam [1 2])
        medk = invz_emt_scalar(G0k, Sk, Jnu, struct());
        if ~medk.converged, break; end
        lamk = invz_lambdas(medk.K, gk, wtsk, betak, [1 2]);
        sgk  = invz_sigma(tlk, lamk, medk.K, gk, betak);
        dS = max(abs(sgk.Sigma - Sk));  Sk = Sk + 0.7*(sgk.Sigma - Sk);
        if dS < 1e-8, okT = true; break; end
    end
    verifyTrue(testCase, okT, sprintf('route-B node T = %.4g K did not converge', Tg(k)));
    % one SIMULTANEOUS final tuple, recomputed from the converged Sk (third-review P1)
    medk = invz_emt_scalar(G0k, Sk, Jnu, struct());
    verifyTrue(testCase, medk.converged);
    lamk = invz_lambdas(medk.K, gk, wtsk, betak, [1 2]);
    sgk  = invz_sigma(tlk, lamk, medk.K, gk, betak);
    Gk = G0k ./ (1 + Sk + medk.K .* G0k);             % J 2.30
    dU(k) = 0.5*( sgk.alpha*tlk.n01*tlk.Delta/(1 + sgk.alpha) - tlk.M2*lamk(1) ...
                  + real(sum(wtsk .* medk.K .* (Gk - G0k)))/betak );
end
verifyTrue(testCase, all(isfinite(dU)));              % NO skipped nodes in a blocking test
dFB   = + T * trapz(Tg, dU ./ Tg.^2);                 % the DENSE value: THE gate input
dFB_c = + T * trapz(Tg(1:2:end), dU(1:2:end) ./ Tg(1:2:end).^2);   % true 1.35-spaced subset
verifyEqual(testCase, dFB_c, dFB, 'RelTol', 0.05);    % temperature-DENSITY convergence
nTm = nnz(Tg <= T*1.35^14);                           % shorter-ceiling prefix (intact grid)
dFBt = + T * trapz(Tg(1:nTm), dU(1:nTm) ./ Tg(1:nTm).^2);
verifyEqual(testCase, dFBt, dFB, 'RelTol', 0.05);     % Tmax-extension convergence
tail_est_B = + T * abs(dU(end)) / Tg(end);
verifyLessThan(testCase, tail_est_B, 0.05*max(abs(dFA), abs(dFB)));
verifyLessThan(testCase, dA.tail_est_A, 0.05*max(abs(dFA), abs(dFB)));
% ---- THE GATE (PENDING -- SS7b P0-2) -----------------------------------------------
% The tolerance below is finalized ONLY after the measured 13-16% closed-model
% discrepancy is localized and resolved (Step 3 diagnostics -> Step 4 escalation).
% Do NOT commit this test while it fails; do NOT widen the tolerance by fiat.
verifyEqual(testCase, dFA, dFB, 'RelTol', 0.10, ...
    sprintf('closed-model two-route mismatch (A=%.4g, B=%.4g) -- BLOCKED escalation', dFA, dFB));
end
```
(`invz_sigma`'s exact signature/outputs: mirror `invz_solve_point`'s `local_sigma_loop` usage. The convergence sweeps are now EXECUTED calls, not comments — third-review P1.)

- [ ] **Step 3: Diagnosis scope note (wording per third-review P2 — no validation claim while blocked)**

In `invzp_QCP_diagnosis.md`, after the Stage-2 STATUS line, append: `The J 2.34 two-route validation has been MOVED to a closed 2x2 reference model exercising the shared kernels (plan SS7b) and is PENDING resolution of a recorded closed-model mismatch; the production hybrid's saturation-normalized absolute dF is a DOCUMENTED DOMAIN LIMITATION -- route A diverges once the full-manifold moment path enters CF-mixing scales (measured: dFA -0.0094/-0.0032/+0.49 at hmax_fac 4/8/300 vs route B -0.0376; task-6 execution record) -- and invz_deltaF_ordered is a partial, cutoff-dependent diagnostic only.` The affirmative "is validated" wording (and the matching `invz_deltaF_ordered` docstring phrase — say "designated validation route", not "has validated") is written ONLY in post-resolution Step 6.

- [ ] **Step 4: Localization diagnostics (second-review P0-2 — MANDATORY before any gate claim)**

Run the corrected harness and record `dFA`, `dFB`, and every convergence number. Then localize the expected residual mismatch (reviewer-measured 13–16% on this fixture):
1. **δU via J 2.21 vs J 2.22 (with the third-review P1 convergence-factor guard):** transcribe the framework's un-reduced internal-energy form (J 2.21; read the framework §8 around worktree eq 36 — do NOT guess it). J 2.21 contains the frequency-ODD `(Δ+ωn)G` contribution and the framework says it requires a convergence factor — a direct substitution into the repo's nonnegative doubled-weight Matsubara sum is INVALID (it can even destroy the bare-limit cancellation). Before using it diagnostically, implement ONE executable treatment: (a) reconstruct the signed grid and take the convergence-factor limit with its analytic high-frequency remainder, or (b) derive a symmetrized expression in which the convergence-factor piece is analytic. Then a ZERO-ORDER GATE: with `K = Σ = 0` the computed `δU_21` must vanish to numerical tolerance — only after that gate passes, compare J 2.21 and J 2.22 on the interacting final tuple. If the two δU's differ materially, the discrepancy sits in the ENERGY REDUCTION; if they agree, it sits in the FIELD route or the resummation order.
2. **Same-final-tuple audit:** recompute route-B's δU strictly from each temperature's final converged `(Sigma, K, lam)` — no in-loop values — and confirm bit-level agreement with the test's own computation.
3. **Coupling-scale exponent scan (quantity DEFINED — third-review P1):** repeat the two-route comparison at `Jnu·s` for `s = 1, 1/2, 1/4` and fit `log ε(s)` vs `log s` for the RELATIVE mismatch `ε(s) = |dFA(s) − dFB(s)| / max(|dFA(s)|, |dFB(s)|)` (the absolute difference has a different leading power and would misclassify a same-order defect as higher-order). Also record the individual scaling exponents of `dFA(s)` and `dFB(s)` to confirm the denominator's leading behavior. Interpretation: ε-exponent ≈ 0 → same-order transcription/order inconsistency; positive exponent → a defect vanishing toward weak coupling. Three points are DIAGNOSTIC, not a derivation — any new gate still requires the analytic resummation-order argument §7b demands.
Write all numbers to the report file.

- [x] **Step 5: BLOCKED escalation — COMPLETED AS DESIGNED.** Evidence delivered (task-6b-report, probes report); independently adjudicated (`invzp_task6b_blocker_review_Codex.md`, against the published Jensen paper): implementation NOT shown defective; the 13.65% is a probable same-retained-order static-elastic approximation residual; resolution = the three-outcome structure of §7b, proceeding on outcome 3's engineering path with the HoF3 discriminator as the follow-up scientific experiment.

- [ ] **Step 6 (ENGINEERING RESOLUTION — per the blocker review; the ONLY committable form):**

Amend the Step-2 harness test before committing: REMOVE the `verifyEqual(dFA, dFB, 'RelTol', 0.10)` equality gate entirely (no percentage equality — 10% or 15% — may be committed as a physics gate). In its place:
1. **Order-consistency gate:** run the harness at `s = 1, 1/2, 1/4` (`Jnu·s`); fit `p_A, p_B, p_Δ` (the |dFA|, |dFB|, |dFA−dFB| exponents) and `p_ε`; assert `1.9 ≤ p_A, p_B, p_Δ ≤ 2.1` and `|p_ε| ≤ 0.15`, commented as engineering test tolerances encoding the observed leading-order structure, NOT analytic bounds.
2. **Approximation-fingerprint regression (non-gating diagnostic):** pin `dFA` and `dFB` at T = 0.31 K, plus the T = 0.10 K cross-check pair, each within 2% of the recorded values (`−5.2317e-4 / −6.0585e-4` and `−5.1364e-4 / −6.0060e-4`), commented as "approximation regression — detects code drift; NOT a two-route validation; the ~13.7% residual is the documented same-retained-order static-elastic approximation fingerprint (Jensen PRB 49, 11833: published full-machinery check = 2–3% on the physical HoF3 lattice at low T; see §7b)".
3. If the s-scan makes the test slower than ~2 minutes, gate it `INVZ_SLOW` and note the count consequence.
Diagnosis wording (replaces Step 3's pending text): `The J 2.34 two-route comparison in the closed 2x2 model shows a ~13.7% same-retained-order static-elastic approximation residual (scale-free in the coupling; Jensen's own full-machinery published check achieved 2-3% on the physical HoF3 lattice at low T, with his stated high-T elastic caveat). The implemented J 2.26-2.29 chain is not shown defective; thermodynamic cross-route closure is carried as a KNOWN APPROXIMATION LIMITATION, scientifically unvalidated pending the published-benchmark (HoF3) discriminator -- a follow-up task. Exact limits, convergence, and leading-order scaling ARE enforced in test_invz_deltaF.` The `invz_deltaF_ordered` docstring says "designated diagnostic route" — never "validated".
Then: test file green; fast suite expected **163 passed / 0 failed / 21 incomplete** (old slow gate gone; the harness test +1 fast — or 162/0/22 if the s-scan forced an `INVZ_SLOW` gate; state which); rerun the two stage-2 slow tests (`test_invz_qcp_closure`, hmf near-boundary) with `INVZ_SLOW=1` to confirm NO remaining known-failing gates on the branch. Then:

```bash
git add invz_projected/invz_deltaF_ordered.m invz_projected/tests/test_invz_deltaF.m invzp_QCP_diagnosis.md
git commit -m "fix(invz): J 2.34 closed-model check lands as order-consistency gates + pinned approximation-fingerprint diagnostic (blocker-review resolution: same-retained-order static-elastic residual, implementation not shown defective; HoF3 discriminator deferred) -- old failing slow gate removed (stage-2 task 6b)"
```

## Cost and risk notes

- **Cost:** a jensen ordered solve ≈ `nH` (33) node solves + ≤ 12 refinement evaluations, each = one fixed-field single-ion + one warm-started Σ↔EMT loop with the scalar static closure (cheap: ≤ 200 scalar iterations over a 24-element mean) ≈ 15–30× a bare ordered solve. Production 301-point sweeps: the ordered side moves from seconds to tens of seconds per field; `parfor` over fields (already supported) keeps wall-clock in the tens-of-minutes range. The synthetic-coupling tests stay fast except the two INVZ_SLOW gates.
- **Near-boundary convergence:** node solves use order-free `hz_fixed` (no MF Picard on the ordering channel at all); the Σ loop's `max_outer`/`mix_outer` knobs apply; jensen acceptance is root-existence, so there is no `m_tol` masking band (P1-4).
- **Known residual approximations** (documented, not silent): the two-level ξ at finite m AND the transverse-feedback bucketing into `G0el` (§6, round-4 P2-F — both exact in the two load-bearing limits); the geometric grid's first panel (O(`hmax·hmin_frac`), 2nd-order, probed by Gate 8); saturation/temperature tails in Task 6b (verified/analytic in the closed model); **the production hybrid's saturation-normalized absolute δF is OUTSIDE its validated domain** (§7b — route-A divergence at CF-mixing scales, numbers preserved in the task-6 record; `invz_deltaF_ordered` is a partial cutoff-dependent diagnostic only).
- **What failure looks like:** identity-gate, closure-gate, or two-route failure = BLOCKED escalation with numbers (§9). Silent sign fixes are forbidden.

## For the independent reviewer — what to check

1. **Transcriptions** (§3, §4, §7) against `jensen_1z_framework.html` §§9.1–9.4 (working-tree eqs 39–46; J 2.26–2.34) — every formula in Tasks 1, 2, 3, 6 traces to one of them; the framework's own §9.2 box records a term-by-term check against Jensen's published paper. Task 1's gates run in the two-level parametrization, in which the hybrid form IS J 2.28 verbatim.
2. **The fixed-field contract** (§1, P0-1): `hz_fixed` semantics in `invz_single_ion` (held only without `order`), the node/final-state construction, and the `invz:hzFixed` asserts.
2a. **The two derivative conventions** (round-3 P0-1/P0-2, §3 bullet + §4a): the J 2.31 transcription gate runs on a CLOSED analytic 2×2 fixture (the repo's re-embedded doublet fails the identity by up to 31% — that failure is expected and correct); the production `G0bare` is the PATH derivative including transverse-MF feedback via the static-tensor chain rule, gated by finite-field FD against the actual node path. Check the chain-rule block algebra (three-way mode switch: `'none'` fb = 0, `'legacy_x'` 1×1, `'vector_ab'` 2×2 — round-4 P1-A), that the FD gate covers all three modes, and that `G0inel` stays fixed-Hamiltonian (consistent with the frequency-dependent EMT propagator).
3. **The boundary-preserving hybrid static propagator** (§3, round-2 P0-A — the load-bearing round-2 change): J 2.28's structure on the FULL-electronuclear static weights with only ξ two-level; the measured failure of the pure two-level closure (K0 23% off, crit −0.011 vs −0.062); the m→0 symmetry argument for the vanishing elastic weight; and Gate 6b — the h→0 loop must reproduce the ACTUAL `invz_solve_point` fixed point (`K(1)`, `Sigma0`, `crit`).
4. **The static-closure design** (§3, P0-2; Task 2): the elastic sector breaks the ordinary Dyson structure, so `K(0)` comes from the scalar fixed point enforcing ⟨G(q,0)⟩_q = G(0) with `K0 = ⟨J·Gq⟩/⟨Gq⟩`; the m→0 degeneracy with `invz_emt_scalar`'s direct solve is the exactness anchor (gate C1, both weight scales); the closed `K(0)` feeds back into λ and Σ each outer pass, and into `pt.G(1)` at the final state (P1-D).
5. **The onset-coincidence derivation** (§5, as corrected in round 2): the m→0 limit `r → 1+Σ(0)` via the K-cancellation, WITH the indirect `K(0)`-dependence of `Σ(0)` resolved by the hybrid (this was the round-1 flaw); the small-field slope `F(h)/h → crit`. Gates 6 and 6b test it numerically.
6. **The free-energy validation** (§7, §7b): the corrected sign `δF = +T∫_T^∞ δU/T′²dT′`; the closed 2×2 model's one-moment conjugacy chain (the rejected §7a draft's coordinate mixing and the re-embedded-doublet non-quench are the cautionary record); the verified-saturation criteria (absolute+relative, ≥5 tail nodes, primitives, analytic `M_sat`); and the production `dF_partial` scope limitation with the preserved failure numbers.
7. **The hybrid sector mapping** (§6, revised): full-sector static weights + two-level ξ and dynamics; boundary-exactness by construction; the finite-m choices are DOCUMENTED SCOPE LIMITATIONS — the closed-model J 2.34 test does not probe them (§7b, third review).
8. **Root robustness** (P1-4, P1-C, round-3 P0-3): geometric grid; the ADAPTIVE downward extension driven by the INDEPENDENT h=0 PM predictor (`slope_pred`, computed before any profile sampling — the round-2 predicate was self-contradictory and unreachable); the explicit `'unresolved'` status when the floor is hit (never a silent PM label, propagated to `converged = false` in the solver); direct-evaluation bisection; root-existence acceptance (no `m_tol` band); the grid-convergence gate (Gate 8); and the NON-CIRCULAR near-boundary test — boundary refined from the PM `crit = 0` side, coarse `hmin_frac` forcing `n_extend ≥ 1`.
9. **The pole-based closure regression** (P1-5, round-3 P1-6): `D_uni = 1 + (J0eff−K0)·Gstat` as the ordered pole observable, its m→0 degeneracy with `crit`; the refinement subtest is a LOCAL branch tracker (lowest local minimum + same-basin check), with the finite-m non-identity between the ω→0 inverse response and `D_uni` documented (the real-axis pipeline excludes the elastic sector).
10. **Plumbing and contracts** (P1-6, P1-F, P2-G, round-3 P1-4/P1-5/P2-7): full-opts threading into every node including `emt_static` and the bare bracket's MF knobs; the H_MF `max_outer` default ALIGNED with the solvers' 200; the `B(3)` and `forced_moment` guards; `ordered_mode` reserved in `invz_check_solve_opts`; the `pt.crit` (legacy diagnostic) vs `pt.D_uni` (ordered pole mass) contract; the DISCRIMINATIVE `pt.G(1) = Gstat` test (the ordinary-Dyson value must fail it); `endpoint_dh` (production — non-claiming) and `tail_est_A`/`tail_est_B` (closed model — estimates) with EXECUTED convergence sweeps and the `hmax_abs` exact-override common-cutoff option (§7b, Task 6b).
11. **Back-compat claims**: `'bare'` byte-identity in Task 4; the auto/overlay leg untouched in Task 5; the `'bare'` escape hatch semantics.
