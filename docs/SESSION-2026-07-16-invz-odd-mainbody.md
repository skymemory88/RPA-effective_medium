# Session note — 2026-07-16/17: ODD main body (off-diagonal dipolar coupling) — handoff

Branch: `invz-1z-lihof4`. This is the handoff for the **ODD main-body feature**
(`odd_implementation_plan.html`; plan `docs/superpowers/plans/2026-07-16-invz-odd-mainbody.md`)
— phases P0 → T1 → T2 → T3 → V4, Tasks 1–11. The off-diagonal dipolar (ODD)
mechanism is a **strictly additive, opt-in** extension of the scalar-cc 1/z
theory: with the flags off the module is bit-for-bit the published Rönnow-2007
route (fast baseline **143 / 0 / 19**, frozen at `5f4ff92`, re-verified every
task). The running measurement log with every headline number is
`docs/ODD-LOG.md` (§P0 … §V4). The physics rationale and the equation set
(E1)–(E5) live in README §1.9. This note is the module map + design-decision
index.

## Headline result

Zero-field ordering temperature with Tier-1 ODD, Richardson (12³, 24³):

| quantity | value | error bar (grid) |
|---|---|---|
| Tc(0), 1/z baseline (ODD off) | **1.74347 K** | ± 0.024 K |
| Tc(0), 1/z + Tier-1 ODD | **1.50937 K** | ± 0.023 K |
| ΔTc(0) | **0.2341 K** (13.4 %) | — |
| Σc(0): off → ODD | 0.29798 → **0.38880** (ΔΣc = +0.091) | — |
| d (E5 uniform reduction) at Tc | **0.483 μeV** | ± (dpRng) see §V4.2 |

Physics context (qualitative; models differ — never rescale): the no-ODD 1/z
baseline sits ≈ 0.21 K above experiment (1.743 K vs Tc_exp ≈ 1.53 K). Tier-1 ODD
delivers 0.234 K of suppression, landing at **1.509 K** — i.e. it closes ≈ 111 %
of that gap (a **21 mK overshoot**: Tc_ODD 1.509 K vs exp 1.53 K). This is *more*
suppression than DS2023's 3-state MF (≈ 5 %) because the 1/z fluctuation channel
amplifies the ODD coupling's q-structure; the comparison to DS2023 is
directional/qualitative only.

## Module map — new functions (all under `invz_projected/`)

| function | signature | one line |
|---|---|---|
| `invz_odd_blocks`   | `[Vca,Vcb,Vcc,info] = invz_odd_blocks(ion,qvec,opts)` | geometric ODD blocks `J^{ca},J^{cb} = -gfac·dip(3,1|2)` (E1a, dipole-only, not Hermitized) + the cc block `invz_jq_modes` eigendecomposes; cached under `odd1_`; `info.Jcc0/Jaa0`. Demag guard (`invz:oddDemag`, intrinsic-only). |
| `invz_chiperp`      | `[Xp,info] = invz_chiperp(ion,T,B,opts)` | 2×2 transverse (a,b) single-ion susceptibility (meV⁻¹), the symmetrized block of the full 136-state electronuclear χ₀; Van Vleck-dominated (χ_aa=χ_bb≈17.64 at 1.53 K, 0 T); evaluated ONCE per (T,Bx), never iterated inside the EMT loop. Accepts `z = iωn` for the retarded route. |
| `invz_odd_deltaJ`   | `[dJ,d,dinfo] = invz_odd_deltaJ(Vca,Vcb,Xp)` | ODD-mediated coupling δJ^cc(q) (E1 contraction, PSD by construction), E4 self-site subtraction (diagonal BZ-mean removed), E5 uniform reduction `d = mean_s mean_q dJpre(s,s,q)`. |
| `invz_odd_zero_field` | `[Tc,out] = invz_odd_zero_field(ion,opts)` | zero-field Tc(0) with/without ODD over grids (Richardson); modes `off\|full\|uniform_only\|qstruct_only`; `out.split` = the GOVERNING condition/Σ-space MF/fluctuation decomposition; `out.Sc_rich`, `out.d_at_Tc`, `out.nex`. The V4.1 zero-field endpoint source AND the driver's odd-aware Tc0 anchor. |
| `invz_odd_fieldvar` | `[C,info] = invz_odd_fieldvar(ion,pt,S,T,opts)` | Tier-2 internal transverse-field covariance C (E3, meV²) at a CONVERGED ODD point solve, from the equal-time cc structure factor of the converged EMT propagator; `info.heq_T ≈ 0.30 T`. |
| `invz_twolevel_avg` | `[tla,avg] = invz_twolevel_avg(ion,T,Bx,C,opts)` | Gauss–Hermite quenched-Gaussian dressing of the doublet over N(0,C) internal fields → effective (Δ,M²,n01) fit to the node-averaged response, + (opt-in `G0`) the disorder-averaged electronuclear Ḡ0. **The plan's flagged least rigorous step.** |

## Module map — modified library functions (all strictly additive; flag-off byte-identical)

- `invz_jq_modes` — `opts.odd = false | struct('Xp',[2×2])`; the odd branch rebuilds cc modes from `Vcc + δJ` and applies `info.Jcc0 = infoB.Jcc0 − d` (E5). Never reads/writes the `jq4_` cache on the odd path.
- `invz_solve_point` / `invz_solve_point_ordered` — `opts.odd` (+ `opts.odd_blocks`, `Jnu_flat=[]`); rebuild modes per q, apply `J0eff ← J0eff − d` exactly once; `pt.odd`. Plus `opts.odd_retarded(_exact)`, `opts.odd_rn_override`, `opts.odd_tier2` (+ `tol_tier2`, `max_tier2`); `pt.C/tier2_iters/tla/tier2_resid`.
- `invz_emt_scalar` — `Jnu_flat` may be `[nJ,nw]` (per-frequency mode spectra) for the retarded route; a vector input is textually the original path; identical columns are bitwise-equal.
- `invz_crit_at` — `invz:odd*` structural errors rethrow (not absorbed as an ordered-phase verdict). Zero threading changes (opts/Jf forward verbatim).
- `invz_critical` — untouched; the **para_edge fallback** is what lets B-cut finders survive with ODD on (the ordered side never re-converges).
- `invz_critical_T` — `invz:oddTc0` guard in `adaptive_anchor` (must not silently anchor at the no-ODD Tc0). **CANNOT bracket with ODD on** (no para_edge analog — the T-cut wall, §V4.1).
- `invz_critical_T0field(ion, Sc, J0eff)` — `Sc`/`J0eff` may each be numeric OR a `function_handle` of T (the ODD T-dependent Sc(T)/J0(T)); `Sc=0` gives the mean-field-RPA Tc(0).
- `invz_ion` — `ion.odd = 0` documented default (drivers read `ion.odd`; libraries read `opts.odd`).
- `invz_run_phase_diagram` — `ion.odd` branch (blocks built ONCE pre-parfor, P0.4) for the standard sweep; **and (V4.1, this task) an opt-in `overlay_quick` block** — the headline ODD overlay (four boundary curves + closed-form endpoints + experimental anchors + a Σ(0) panel), B-cuts only, `Data/Phase_ODD_overlay_quick.fig`.

## Flag surface (complete)

- **Drivers:** `ion.odd` (0). Intrinsic-only, requires `ion.demag = 0`. `overlay_quick` (script-local, V4.1 quick overlay).
- **Point solve:** `opts.odd` (false) + `opts.odd_blocks = struct('Vca','Vcb','Vcc','Jcc0'(unshifted))` + `Jnu_flat = []`.
- **Retardation (T2, static frozen as default):** `opts.odd_retarded` (requires odd), `opts.odd_retarded_exact` (full per-(q,n) eig cross-check, wins when both set), `opts.odd_rn_override` (scalar or [nwn,1] test hook).
- **Tier 2 (variable moments):** `opts.odd_tier2` (requires odd; errors `invz:oddArgs` if combined with `odd_retarded(_exact)`), `opts.tol_tier2` (1e-3 rel), `opts.max_tier2` (8).
- **`invz_chiperp` opts:** `transverse_mf` ('legacy_x'), `hyp` (true), `Jxx0`, plus `z = iωn` (Matsubara).
- **`invz_odd_fieldvar` opts:** `Ecut` (40, must match the solve's grid), `static_approx` (false; never silent — `info.static_approx`).
- **`invz_twolevel_avg` opts:** `avg` ('response' default | 'params'), `ngh` (7; machine-converged), `G0` (false; opt-in, 49 electronuclear diagonalizations), `wn`, `Ecut`.
- **Zero-field engine opts:** `grids` ({12,24}), `dpRng` (30), `cache` (true), `mode` ('full').

## Cache keys

ODD geometric blocks have their **own** namespace under `invz_projected/cache/`:
`odd1_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat`, schema tag `odd1`,
`pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]` (25 elems); one file
per grid stores `Vca,Vcb,Vcc,info,pkey,qvec`; the loader `isequal`-verifies BOTH
`pkey` and `qvec`. The `Vca/Vcb` blocks are dipole-only (J12-independent); only
the `Vcc` block moves with a J12 change. The existing `jq4_` caches are **never
touched** by any ODD path (verified by code inspection and by 0 new `jq4_*.mat`
after odd calls).

## Adjudication history (one paragraph)

The source plan's ΔTc-space decomposition was **ill-posed** (measured on
12³/24³ Richardson) and produced a 0.461 K "closure defect"; the controller
re-adjudicated it (2026-07-17, after Task 6's BLOCKED round) into a **sequential
condition/Σ-space factorial** in (J0-shift) × (Σ-source), neither leg of which
enters the invalid regime: legs (a) full 1.509 / (b) condition-level 1.615 /
(c) Σ-level 1.630 K, closure defect **+0.008 K** (3.4 % of the effect).
Diagnostics kept and REPORTED: `c1_literal ≡ a` (a closed-form theorem
validating the E4/E5 bookkeeping — a simultaneous +d shift of couplings AND J0
leaves R2007 criticality invariant), `c_factorial` 1.794, `b_naive` 1.283 (96
excluded modes at 24³ — exactly DS2023's naive-MF inconsistency, where Jcc0−d
drops below the peak finite-q Vcc mode). Three further adjudications: **M2
convexity** (Task 9) — the plan's "M²_eff ≤ M²₀ monotone in ‖C‖" is measurably
FALSE at 2 T (bare M2 is convex there, M2''(2T)=+0.751/T²), so the gate was moved
to χ0_2l suppression and M2_eff is REPORTED; the **T-cut wall** (T2.2) — with ODD
on the ordered side never converges to a PM fixed point, so `invz_critical_T`'s
converged-sign-change bracket is structurally unreachable and only the B-cut
finder's `para_edge` fallback survives (drives every ODD phase-boundary
measurement to B-cuts); and **PM fixtures** — convergence-gated tests were moved
to guaranteed-paramagnetic points (e.g. 1.80 K / 0.05 T) because the plan's
original near-boundary fixtures sit in the non-convergent band.

## T2.2 decision (retardation)

Measured at the 2 T proxy (16³/dpRng 30): despite substantial retardation in the
transverse channels (r(ω₁=0.828 meV) = **0.679**, monotone to r₄₉=0.006), the
boundary shift is **negligible** — |ΔTc| = **0.022 mK ≪ 10 mK**, ΔBc(1.2 K)=0 to
crossing tolerance. Criticality lives in the elastic n=0 sector where r₁=1
exactly and χ⊥(0) is exact, so Tc(0) is retardation-invariant by construction.
Per the plan's 10 mK rule the **static form is the default**; the retarded path
stays behind `opts.odd_retarded` (surrogate `χ⊥(iωn) ≈ r_n·χ⊥(0)`, exact in the
δJ matrix, first-order in the mode spectrum, surrogate-vs-exact error 2e-11 on
Σ0; exact per-(q,n) cross-check behind `opts.odd_retarded_exact`).

## T3.4 verdict (small-Bx / B = 0 Tier 2)

Zero-field / small-Bx Tier-2 results are **directional only** (sign and rough
magnitude), not quantitative: the dressed splitting at the 0.05 T proxy is
Δ_eff ≈ 11.5–11.8 μeV, only ≈ 3.4× the hyperfine scale A ≈ 3.36 μeV, and the
quenched-Gaussian single-pole recompression is the plan's flagged least-rigorous
step. Exact B = 0 still throws `invz:degenerateDoublet` from the bare Tier-1
two-level scaffold (deliberate — the Σ machinery needs a finite bare doublet):
**B = 0 is owned by the closed form** (`invz_odd_zero_field` /
`invz_critical_T0field`), with 0.05 T the practical small-field proxy. The
Tier-2 dressing machinery itself IS B=0-capable (degenerate-node limit), noted as
a Tier-3 design option, not implemented (YAGNI).

## Open items — Tier 3 (deferred; plan §7)

Tier 3 = full-tensor ODD mode-mixing RPA (ODD lines hybridizing the cc soft mode
with the transverse excitons), out of scope for the main body:

1. **Matrix EMT** — the effective-medium self-energy generalized from the scalar
   cc channel to the full sublattice/tensor propagator.
2. **Lost single-channel identities** — the scalar sum rule and PSD identities
   that the matrix generalization must recover.
3. **Three-level minimal sector** — the minimal CF sector (beyond the doublet)
   needed to capture the transverse-exciton hybridization.
4. **CF-state convergence** — how many CF levels the tensor route needs.
5. **Tier1+2 ↔ Tier3 cross-validation target** — the quantitative check that the
   Gaussian Tier-1+2 result is recovered as the Gaussian limit of the tensor RPA.

**Deferred full-tensor route (Appendix A):** the `invz_tensor/` scaffold and its
plan `docs/superpowers/plans/2026-07-16-invz-tensor-odd.md` — DEFERRED by the
user decision of record (2026-07-16); this main-body plan was the active stream.
No `invz_common/` refactor was done (single active branch).

## Production runs left to the user

- **Phase-diagram production sweep** — the V4.1 overlay ships only the `quick`
  mode (≈ 3 B-cut temperatures + closed-form endpoints, in-session). The dense
  boundary (fine T grid, full para-edge boundary, per config: off / +Tier-1 /
  +Tier-1+2) is hours-long and LEFT TO THE USER (repo precedent): set
  `ion.odd = 1` (and/or `opts.odd_tier2`), widen `Ts`, keep `overlay_quick`
  false, run the standard parfor path once per config.
- **Spectra with ODD** — real-axis χ(q,ν,ω) with the ODD δJ folded in is **out
  of scope** (the ODD δJ enters the ordering channel; a full spectral treatment
  is Tier 3 mode-mixing).

## Verification (V4.3 consolidation)

- Fast suite: **143 / 0 / 19** (matches the frozen `5f4ff92` baseline exactly;
  0 failed = flag-off parity, published benchmark digits unchanged).
- ODD fast-test additions total **22.35 s** across the 8 ODD test files
  (`test_invz_chiperp`, `test_invz_odd_blocks`, `test_invz_odd_fieldvar`,
  `test_invz_odd_physics`, `test_invz_odd_retarded`, `test_invz_odd_solve`,
  `test_invz_odd_tier2`, `test_invz_twolevel_avg`) — **under** the 30 s budget,
  so no INVZ_SLOW demotion was needed.
- Slow suite (`INVZ_SLOW=1`, everything): see `docs/ODD-LOG.md` §V4.3 for totals.
