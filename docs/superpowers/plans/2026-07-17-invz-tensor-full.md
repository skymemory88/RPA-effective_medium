# invz_tensor — Full-Tensor 1/z Expansion (A0→A4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **REVISED 2026-07-17 (v2), incorporating the accepted findings of `full_tensor_review_by_Codex.md`:**
> suite isolation (tensor core suite runs with `invz_projected` OFF the path; projected parity lives
> in an optional `tests/interop/` suite; architecture decision by user 2026-07-17: shared
> `invz_common` retained, NOT self-contained copies); A2 rewritten as a DIRECT per-frequency
> closure with the corrected positive scalar sign; the A1 zero-field call replaced by the small-B
> proxy + a tensor-owned extrapolator (the two-level helper throws `invz:degenerateDoublet` at
> B = 0); `Sigma_cc_equiv` sign corrected to `+V/G0`; one tensor-owned q-grid contract with
> explicit half-open vs legacy-inclusive conventions; invalid test fixtures repaired (single-q
> `invz_odd_deltaJ` contract, S4-on-arbitrary-q, AFM-exchange Lorentz assert, Task-8 plumbing);
> the three-state model given an explicit Hamiltonian/operator contract; the A3 performance gate
> moved INTO Task 11 with a required factorized contraction design; KMS log-scaled kernels and
> basis-content ladder rungs; the negative-frequency tensor transpose relation locked; a `lat`
> struct threading grid/weights/`Jaa0` provenance through every solver; λ-protocol for the
> Gaussian emergence test; tests/drivers never write to `docs/ODD-LOG.md`.

> **REVISED 2026-07-17 (v3), incorporating the second review pass of `full_tensor_review_by_Codex.md`:**
> (Blocking) the Task-11 O(N³) "factored" vertex is reframed as an UNPROVEN conjecture — the DENSE
> ordered path-sum is now the correctness reference AND the production engine; `'factored'` is opt-in,
> gated behind an explicit identity derived+oracle-verified in Task 10, and "factorization not
> established" is an explicit Task-11 stop trigger (dense-only, budget-refuse large rungs). (Blocking)
> the Task-12 three-state model gains an INDEPENDENT low-doublet tunnelling parameter `Delta1` (review
> option (a)): (Δ1, m0, ρ) solve the three targets (Δ, M2, χ⊥) — no longer overconstrained — and ρ→0
> now leaves a genuine two-level system (splitting from Δ1, not ρ), making the exact scalar-compat
> identity consistent. Also: Task-1 path audit widened repo-wide + `which -all` collision check; Task-6
> zero-field guard on the TRANSVERSE component (`abs(Bx)`), not `norm(B)`; Task-6 crit square root by
> Hermitian eigendecomposition (not `sqrtm`) with rank clipping; Task-6 odd-on/off crit DIRECTION
> demoted to a measurement (the PSD Schur gate stays hard in Task 4); Task-9 S4 fixture keeps the
> COMPLEX Hermitian DFT (no `real()`, which changes eigenvalues), scalar helper reads only the 4 cc
> branches; Task-10 KMS check compares the ENTIRE contracted value vs mpmath at large βΔ with
> complex-signed/log accumulation; Task-12 λ-emergence fits an asymptotic slope over ≥3 λ (ratio-at-1
> demoted to a report); Task-13 virtual deficit "diagnoses" (not "bounds"), rungs record ACTUAL
> multiplet-complete dimension; A3 true-zero-field Tc explicitly DEFERRED (small-Bx proxy for all rungs;
> only the projected closed form is truly B=0); Task-2 frozen baseline uses a committed generator +
> filtered dirty flag over physics-source paths.

> **SUPERSEDES** `2026-07-16-invz-tensor-odd.md` (the deferred A0+A1+Tier2 bridge plan). That plan
> was written before the ODD main body landed in `invz_projected/`; this plan absorbs its A0/A1
> layer, fixes its two physics defects (the Task-8 additive ΔTc "(a)/(b)/(c)" decomposition,
> proven ill-posed by the uniform-shift theorem, framework §11.3/§11.6; and the Task-7 K-bookkeeping
> formula, see Task 6 here), drops its Tier-2-additive re-implementation (scope box below), and
> adds the genuinely new content: A2 matrix effective medium and A3 tensor 1/z self-energy with
> the A4 state-space ladder. The old plan carries its SUPERSEDED banner.

**Goal:** Implement the full-tensor (3 Cartesian × 4 sublattice, multi-level) version of Jensen's 1/z effective-medium expansion for LiHoF4 in `invz_tensor/`, staged A0 (tensor-RPA parity layer) → A1 (projected-1/z bridge) → A2 (matrix effective medium) → A3 (genuine tensor 1/z self-energy from exact component-labelled four-point cumulants) → A4 (state-space ladder), such that A3 reduces to the projected ODD equations E1/E4/E5 (framework §11.8's designed cross-validation, SHARPENED in Task 12) in two sectors: the SCALAR sector EXACTLY (the A3 vertex → `invz_sigma` at ρ→0, ≤ 1e-8), and the ODD sector under MATCHED TRUNCATION (A3 with dominant-only dressing reproduces E1 up to the O(1/z²) resummation-scheme ambiguity of constraint 8 — measured collapse of the beyond-E1 excess) — while the FULL A3 additionally dresses the transverse spectator that Jensen's dominant-only rule leaves bare (a REPORTED beyond-E1 correction, not a defect). Note: the original "Gaussian truncation of A3 reproduces E1 to leading order in λ" premise was found to conflate the vertex emergence with the O(1/z²) resummation ambiguity; the matched-truncation formulation is the corrected, constraint-8-consistent statement. Compatibility claims come in two explicitly labelled categories: **exact scalar-compatibility identities** (synthetic/decoupled inputs, tolerances ≤ 1e-8) and **physical full-tensor no-ODD comparisons** (finite tolerances with named residual sources — cross-Cartesian χ0 elements, dominant-mask vs whole-cc division, grid convention).

**Architecture:** Three layers. (1) `invz_common/` receives the branch-agnostic single-ion engine moved out of `invz_projected/` (pure `git mv`, behavior-neutral; projected tests stay the suite of record for it) — one source of truth, per the user's 2026-07-17 architecture decision. (2) `invz_tensor/` holds `invzt_`-prefixed functions only: the q-grid builder, `[12,12,nq]` coupling tensor, page-wise tensor RPA, the A1 bridge solver, the A2 direct matrix-EMT closure, and the A3 vertex engine (stable KMS-scaled divided-difference kernels + component-labelled path sums, oracle-verified before MATLAB encoding). Its CORE test suite (`invz_tensor/tests`) runs with `invz_projected` absent from the path; live projected-parity tests live in `invz_tensor/tests/interop/` (a separate, optional suite) and frozen projected reference numbers live in committed fixtures. (3) Physics numbers are *reported, never tuned*; hard gates are exact identities, oracle rows, symmetry, convergence, reproducibility, and the small set of direction gates that carry a stated derivation.

**Tech Stack:** MATLAB R2025a (`pagemtimes`/`pagemldivide`), `functiontests(localfunctions)` suites, root lattice machinery (`MF_dipole.m`, `exchange.m`, `qVec_generator.m` stay at repo root), `python3` with numpy+mpmath (oracle generation only — verified working: `verify_three_level.py` 41/41 and `verify_v2_vertex.py` 42/42 pass on this machine, 2026-07-17; the MATLAB suites must NOT depend on python at run time — oracles are committed JSON fixtures).

## Physics base (what each stage encodes — sources pinned)

- **Scalar theory (implemented, verified):** `jensen_1z_framework.html` §§0–10 (= Jensen PRB 49, 11833 Secs. I–II; "J 2.x"), R2007 = Rønnow et al. PRB 75, 054426 (4-sublattice cc application, Σc closed form eq. 10). ODD within 1/z: framework §11 (E1/E4/E5, uniform-shift theorem §11.3, condition/Σ-space sequential decomposition §11.6, retardation null §11.4, quenched dressing §11.5).
- **A0/A1:** `odd_implementation_plan.html` Appendix A (A.1–A.2); Jensen's "dominant transition renormalized, weak transitions kept at RPA" rule (PRB 49 Sec. III: a transition contributes a rank-1 Cartesian dyadic M^α_pq M^β_qp, so all its Cartesian components carry the same 1+Σ factor).
- **A2:** matrix effective-medium closure (Appendix A.2; template J. Jensen, J. Phys. C 17, 5367 (1984); scalar limit = framework §5 eqs. (16)–(18) — a DIRECT closure, like `invz_emt_scalar`'s closed form).
- **A3:** `three_level_1z_extension.html` (exact multi-level first-order vertex: path-sum F_{nℓ} with Hermite divided-difference kernels 𝓘₂/𝓘₃; Γ⁽⁴⁾ = F − pairings; three-point L_n/T³_ℓ; sum-rule consistency; ε_el control; Blume–Capel, Δ→0, T=0 anchors — all verified by the two repo python scripts) + Yang–Wang PRB 10, 4714 (1974) & PRB 12, 1057 (1975) (SBO diagrammatics: transition-space Dyson form with irreducibility w.r.t. interaction lines; Λ/Δ structure; q=0-transfer deletion; kinematic zero-identities).
- **Targets/benchmarks:** DS2023 = Dollberg–Andresen–Schechter PRB 105, L180413 (2022) + arXiv:2308.10095; projected headline (docs/ODD-LOG.md): ΔTc(0) = 0.2341 K (1.743 → 1.509 K vs exp 1.53), d = 0.483 μeV, ΔΣc = +0.091, Tier1:Tier2 = 97.2%:2.8% at the 0.5 T proxy, retardation null 0.022 mK.

**Hard mathematical constraints from the derivation base (bind every A3 task):**
1. **No λ_p closure beyond two levels.** Shifted-resolvent kernels do not reduce to powers of channel propagators; `invz_lambdas`/`invz_sigma` algebra is exact ONLY in the two-level limit. The exact multi-level vertex needs the full frequency-resolved sum Σ_ℓ K(iωℓ)·Γ⁽⁴⁾(n,ℓ).
2. **Carry the additive correction 𝒱, never a ratio.** Multi-pole G₀ has real zeros between poles; Σ = 𝒱/G₀ has artificial poles there. Numerics stores 𝒱_{μν}(iωn); resummation uses the symmetric bracket solve of Task 12, never `inv(G0)` or componentwise division.
3. **Spectator normalization.** Embedding a two-level channel in a larger space requires normalized subspace probabilities p̃ = p/(p_a+p_b) plus the explicit deficit term ∝ (p_a+p_b)(1−p_a−p_b); raw populations inside Jensen's α/γ are wrong. Tested with a POPULATED spectator (Task 10), not only a far-removed one.
4. **KMS boundary layers — as algorithms, not warnings.** Exponentially small initial-state weights pair with exponentially growing time factors into O(1) endpoint layers. Kernels evaluate population-weighted products in COMPLEX-signed log-scaled/regrouped form — a bare real `exp(log p_r + β·x_acc)` cannot carry the sign/phase of the divided-difference cancellations, so accumulation uses a complex log (`exp(clog(p_r) + Σ complex exponents)`) or an explicit analytic KMS regrouping of conjugate boundary terms (Task 11); initial-state truncation only behind an explicit KMS-regrouped bound; large-βΔ endpoint layers are TESTED by comparing the ENTIRE contracted kernel/vertex value (not individual exponentials) against mpmath (Task 10), not only exact degeneracy.
5. **Divided differences are generic.** Return paths force repeated nodes for every ℓ; kernels use exprel-stable φ, Hermite limits 𝓘₂(x,0), 𝓘₃(x,y,−y), 𝓘₃(x,y,0), 𝓘₃(x,0,0), complex arguments. The path sum uses initial-state weights p_r — never population-difference ratios D_pq/(εq−εp) — so degenerate levels (the Bx=0 doublet) are regular by construction; assert no-NaN at exact degeneracy explicitly.
6. **Sum rule is a consistency TEST, not a closure.** The tensor equal-time rule reads (1/β)Σ_n χ_{μν}(iωn) = ½⟨{δJ_μ, δJ_ν}⟩; report residuals with a STATED Matsubara cutoff and tail treatment; never fix constants from it.
7. **Elastic-sector control.** Report ε_el = β|K_cc(0)|·c_d (c_d = elastic cc variance of the level scheme) wherever A3 numbers are quoted; Jensen's tanh form (J 2.29) is one partial resummation, not unique.
8. **Resummation ambiguity is O(1/z²).** Cartesian-space Dyson vs additive differ at O(1/z²); Task 12 measures the spread and reports it as the method error bar.
9. **Negative-frequency tensor relation (LOCKED).** From the Lehmann form, χ_{ρσ}(−iωn) = χ_{σρ}(iωn) — the TRANSPOSE, not componentwise evenness. The same relation holds for K and every channel-matrix object. Scalar/diagonal cases degenerate to evenness. Pinned by a field-on complex oracle case (Task 10) and used for all negative-ℓ reconstruction (Task 11).

## Global Constraints

- Repo root: `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion/`. The path contains spaces — ALWAYS quote it in shell commands.
- **Precondition gate (check before Task 1, again before every commit):** the working tree carries user WIP. If `git status` shows uncommitted modifications under `invz_projected/` or to root `.m` files at Task-1 time, STOP and report BLOCKED asking the user to commit/stash them first (the `git mv` refactor needs a clean projected tree). NEVER `git add -A` / `git add .` anywhere in this plan — stage explicitly named paths only.
- **Suite topology (three suites, run from repo root):**
  - CORE (after EVERY task from Task 3 on; must pass with `invz_projected` absent from the path):
    `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"`
    (`runtests` on a folder does not recurse into `interop/`.)
  - INTEROP (optional; run when a task adds/changes interop tests):
    `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests/interop'); disp(r); assertSuccess(r)"`
  - PROJECTED (after Task 1 and at Task 15 final verification only):
    `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"`
  - Slow gates: prefix `INVZ_SLOW=1` (tests use `assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only')`).
- Non-negotiables: (i) `invz_projected` results are bit-for-bit untouched — this plan only MOVES its files (Task 1) and never edits their logic; (ii) never edit published benchmark values or pinned anchors (`invz_projected/tests/invz_odd_anchors.m` is read-only; a mismatch is a finding to escalate); (iii) physics acceptance criteria marked *report* are measurements — record in `docs/ODD-LOG.md` via the explicit report steps, never tuned; (iv) the §0 stop-rule binds every task: when a measurement contradicts the physics summary (units, symmetry, χ⊥ band 16–17 meV⁻¹, d band 0.3–0.5 μeV), STOP and report BLOCKED.
- **Tests and drivers never write to `docs/ODD-LOG.md` or any tracked file.** Drivers return data structs; `fprintf` in tests is fine. Each task's "report" step is an explicit controller action (running a print/report helper and pasting into ODD-LOG).
- Naming: every new function in `invz_tensor/` uses the `invzt_` prefix; no file may shadow a name in `invz_projected/`, `invz_common/`, or repo root.
- **Tensor index convention (LOCKED):** composite index sublattice-major, Cartesian-minor: `i = 3*(s-1) + mu`, s = 1..4, mu = 1(a), 2(b), 3(c). Slices: cc `Jt(3:3:12, 3:3:12, iq)`, ca `Jt(3:3:12, 1:3:12, iq)`. Local response embeds as `kron(eye(4), chi0_3x3)`.
- **q-grid contract (LOCKED, Task 3):** ALL grids come from `invzt_qgrid(n, conv)` returning `g.qvec` (Γ-excluded), `g.w` (weights summing to 1 INCLUDING the dropped-Γ convention stated in `g.note`), `g.conv` ('halfopen' | 'legacy_inclusive'), `g.hash`. 'halfopen' = `ndgrid((0:n-1)/n)` minus Γ (the default for all new tensor work); 'legacy_inclusive' = exactly `qVec_generator(...,'grid',[n n n],'range',[-0.5 0.5])` + Γ filter (endpoint-inclusive, spacing 1/(n−1) — the projected production convention; used ONLY for parity with projected anchors, on BOTH sides). Grid conv + hash go into cache keys, `lat` provenance, and every reported anchor. Never mix conventions inside one comparison.
- **Lattice struct (LOCKED, Task 3):** solvers take `lat = struct('Jt', [12,12,nq], 'qvec', ..., 'w', ..., 'conv', ..., 'JtGamma', [12,12], 'info', info)` — never a bare array. Solvers thread `lat.info.Jaa0` into the single-ion transverse mean field (`Jxx0`) and `lat.info.Jcc0` where a uniform coupling is needed; `pt` retains `lat` provenance (minus the bulk `Jt` pages: keep `qvec` hash, conv, `JtGamma`, `info`).
- Units and signs (must match `invz_projected`): energies meV, T in K, B in Tesla, χ in meV⁻¹; `J^{mu nu}_{ss'}(q) = -C.gfac*dip(mu,nu,s,sp) + sign(ion.J12)*ex(mu,nu,s,sp)`; Lorentz cavity `lorz = 4*pi/(3*ion.Vc)*C.gfac` at Γ-equivalent q on Cartesian-diagonal entries only, all 16 sublattice pairs. χ = −G; criticality `J(0)·χ0(0) = 1+Σ(0)`, `crit > 0` in the PM phase. Matsubara: bosonic half grid with doubling weights via `invz_matsubara`; `wn(1) = 0` is the elastic slot.
- Demagnetization: tensor layer intrinsic-only; `invzt:demag` error if `ion.demag ~= 0`.
- Published anchors (NEVER rescale): `Jcc0_dipole = 6.821e-3` meV, `Jaa0_dipole = 3.912e-3` meV, `Jcc0 = 6.421e-3` meV (RelTol 0.03); no-ODD Σc Richardson(12³,24³) = 0.2980 (published 0.3004), Tc(0) = 1.743 K; ODD anchors from `invz_odd_anchors()` (interop suite only): χ⊥(1.53 K, 0 T) = 17.63784561529863 meV⁻¹, Tc_ODD(0) = 1.509 K, d(Tc) ≈ 0.483 μeV; DS2023 3-state params Δ = 11.5 K, ρ = 2.34.
- ODD pitfalls binding every task: (i) never grid-extrapolate the q→0 ODD blocks; (ii) ODD small-q decay assertions on HIGH-SYMMETRY on-axis rays only, decay LINEAR in q; (iii) never import DS2023's 0.805 rescaling or perturbative hyperfine; (iv) ODD blocks strictly c↔(a,b); (v) (gL·μB)² slips enter δJ SQUARED — block parity with `invz_odd_blocks` (interop) + the Γ anchors are the unit guards.
- **ΔTc decomposition language (LOCKED):** sequential condition/Σ-space legs (`uniform_only` / `qstruct_only` + closure defect) only; additive "(a)/(b)/(c)" language is FORBIDDEN (uniform-shift theorem).
- **Tier-2 scope box:** the tensor branch does NOT re-implement the Tier-2 additive layer. A1 cross-validates against projected *Tier-1-only* numbers; variable-moments physics re-enters ONLY at A3 via mixed cumulants, cross-validated against the projected Tier-2 share (≈2.8%) as a REPORT. Never stack an A3 result with projected Tier-1/Tier-2 corrections.
- Cache discipline: `invz_tensor/cache/` with `.gitignore` = `*`; keys via `invz_cache_key('jqt1', dpRng, [qhash; convcode], pkey)` with `pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]`; loader `isequal`-verifies stored `pkey`+`qvec` before trusting a hit. Never share files with `invz_projected/cache`. Libraries serial; lattice sums built before any `parfor`.
- Fast tests: n ≤ 8 per axis, `dpRng ≤ 20`; `'cache', true` allowed for shared grids; `'cache', false` where cache-independence is the point. Avoid the non-convergence island (T = 0.31 K, B ∈ [1, 2] T) for convergence-gated fixtures.
- MATLAB test boilerplate — CORE suite (`invz_tensor/tests/*.m`; NOTE: no `invz_projected` path):

```matlab
function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end
```

  INTEROP suite (`invz_tensor/tests/interop/*.m`) adds two lines:

```matlab
addpath(fullfile(here, '..', '..', '..', 'invz_projected'));         % parity targets
addpath(fullfile(here, '..', '..', '..', 'invz_projected', 'tests'));% invz_odd_anchors fixture
```

- Commit style: conventional prefix + scope (`refactor(invz):` Task 1; `feat(invzt):`/`test(invzt):`/`docs(invzt):` after); trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Reporting: each task that produces physics numbers has an explicit controller step appending a dated entry to `docs/ODD-LOG.md` (§A0…§A4): what was implemented/measured, suite status, headline numbers, grid conv + dpRng provenance.

## File structure (all new/moved files at a glance)

| Path | Responsibility |
|---|---|
| `invz_common/` (T1) | 16 moved engine files + `README.md` + new `invz_cache_key.m` |
| `invz_tensor/invzt_qgrid.m` (T3) | THE q-grid builder (halfopen / legacy_inclusive, weights, hash) |
| `invz_tensor/invzt_jq_tensor.m` (T3) | cached `[12,12,nq]` coupling tensor + `lat` assembly |
| `invz_tensor/invzt_chi_rpa.m`, `invzt_gcc_lattice.m` (T4) | page-wise 12×12 tensor RPA; site-diagonal cc average |
| `invz_tensor/invzt_chi0_split.m` (T5) | transition-mask χ0 = χ_dom + χ_rest (elastic convention stated) |
| `invz_tensor/invzt_solve_point.m` (T6, extended T9/T12/T13) | A1 bridge → A2 medium → A3 solver (mode-switched) |
| `invz_tensor/invzt_critical.m`, `invzt_tc_pm_extrap.m` (T7) | PM-side Bc bisection; handle-based small-B Tc extrapolator |
| `invz_tensor/invzt_chi_realaxis.m` (T8) | tensor real-axis spectra (A1 scalar-Σ continuation ONLY) |
| `invz_tensor/invzt_emt_matrix.m` (T9) | A2 DIRECT matrix effective-medium closure |
| `verify_tensor_vertex.py` + `invz_tensor/tests/fixtures/vertex_oracle.json` (T10) | component-labelled vertex oracle |
| `invz_tensor/invzt_kernels.m`, `invzt_vertex4.m`, `invzt_vertex3.m` (T11) | KMS-scaled kernels; factorized path-sum engine + PERF GATE |
| `invz_tensor/invzt_threestate.m`, `invzt_sigma_tensor.m` (T12) | explicit 3-state model contract; A3 Σ assembly |
| `invz_tensor/invzt_run_ladder.m`, `invzt_report_ladder.m` (T13) | A4 basis-defined ladder driver (data only) + report serializer |
| `invz_tensor/tests/*.m`, `tests/interop/*.m`, `tests/invzt_anchors.m` (each task) | core/interop suites; tensor anchors |
| `invz_tensor/tests/fixtures/projected_baseline.json` (T2) | frozen projected reference numbers (hash + dirty flag + grid conv) |
| `invz_tensor/README.md`, `docs/SESSION-*-invz-tensor-full.md` (T15) | docs/handoff |

---

### Task 1: `invz_common/` — factor the branch-shared engine out of `invz_projected/`

**Files:**
- Create: `invz_common/`, `invz_common/README.md`, `invz_common/invz_cache_key.m`
- Move (git mv from `invz_projected/`, exactly 16 files): `getf.m`, `invz_const.m`, `invz_ion.m`, `stevens_ops.m`, `invz_cfrot.m`, `invz_field_vec.m`, `invz_single_ion.m`, `invz_chi0z.m`, `invz_chiperp.m`, `invz_check_transverse_mf.m`, `invz_matsubara.m`, `invz_is_gamma_equiv.m`, `invz_twolevel.m`, `invz_g.m`, `invz_lambdas.m`, `invz_sigma.m`.
- Modify: every test file in `invz_projected/tests/*.m` — in `setupOnce`, add `addpath(fullfile(here, '..', '..', 'invz_common'));`.
- Test: the existing `invz_projected/tests` suite is the acceptance test.

**Interfaces:**
- Produces: the 16 functions, unchanged, resolvable from `<root>/invz_common`; signatures as today (`invz_single_ion(ion, T, B, opts)`, `invz_chi0z(si, T, z, opts)` `[3,3,nz]`, `invz_chiperp(ion, T, B, opts)`, `invz_twolevel(ion, T, Bx, opts)`, `invz_g(tl, z)`, `invz_matsubara(T, Ecut)`, `invz_lambdas(K, g, wts, beta, plist)`, `invz_sigma(tl, lam, K, g, beta)`).
- Produces: `key = invz_cache_key(prefix, dpRng, qmeta, pkey)` — ~15-line helper reproducing the projected weak-hash scheme (`sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'))`), `qmeta` = `[qvec(:); convcode]`. Header documents weak-hash + mandatory isequal-content-verification. Consumed by `invz_tensor` ONLY — never retrofitted into projected files.

- [ ] **Step 1: Precondition + REPO-WIDE blast-radius audit (v3: widened).** Moving 16 UNPREFIXED functions out of `invz_projected/` breaks any caller anywhere in the repo that added only `invz_projected` to the path. Audit the WHOLE tree — nested exploratory scripts, root drivers, `docs/` documented interactive commands — not just `invz_projected/*.m`:

```bash
git status --short
# every addpath anywhere (exclude .git and gitignored caches):
grep -rn --include=*.m "addpath" . | grep -v '/cache/'
# any doc that instructs an interactive addpath of invz_projected:
grep -rn "addpath" docs/ README* *.html 2>/dev/null
# mfilename-relative data/cache path construction in the 16 movers:
grep -rln "mfilename" invz_projected/getf.m invz_projected/invz_const.m invz_projected/invz_ion.m invz_projected/stevens_ops.m invz_projected/invz_cfrot.m invz_projected/invz_field_vec.m invz_projected/invz_single_ion.m invz_projected/invz_chi0z.m invz_projected/invz_chiperp.m invz_projected/invz_check_transverse_mf.m invz_projected/invz_matsubara.m invz_projected/invz_is_gamma_equiv.m invz_projected/invz_twolevel.m invz_projected/invz_g.m invz_projected/invz_lambdas.m invz_projected/invz_sigma.m
```

Then, in a MATLAB session with BOTH `invz_projected` and repo root on the path, run `which -all <name>` for each of the 16 basenames (e.g. `getf`, `invz_const`, `stevens_ops`, …) to detect any pre-existing shadow of a moved name elsewhere in the tree.

Expected: no uncommitted modifications under `invz_projected/` or to root `.m` files (else STOP — precondition gate); no moved file builds `mfilename`-relative data/cache paths; `which -all` shows exactly one definition per name (any collision is a STOP → escalate). Record the FULL caller list. Every caller that added only `invz_projected` MUST be updated in Step 4 (tests) or, if a non-test/non-driver supported entry point exists, gain a one-line path-bootstrap (`addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'))`). If any of these fails, STOP and report BLOCKED.

- [ ] **Step 2: Baseline.** Run the PROJECTED suite; save output to `.superpowers/sdd/baseline_fast_premove.log`. Expected: 0 failed.
- [ ] **Step 3: `git mv` the 16 files**; verify every move is tracked as a rename.
- [ ] **Step 4: Update every `setupOnce` in `invz_projected/tests/*.m`** (+ any driver addpath found in Step 1).
- [ ] **Step 5: Write `invz_common/invz_cache_key.m` and `invz_common/README.md`** (~15 lines: scope, tests-of-record remain `invz_projected/tests`, callers must have `invz_common` on the path, `invz_cache_key` is tensor-only).
- [ ] **Step 6: Run the PROJECTED suite.** Expected: counts identical to Step 2.
- [ ] **Step 7: Commit** (explicit paths) — `refactor(invz): factor branch-shared single-ion engine into invz_common/ (16 files, pure move)`.

---

### Task 2: A0 preflight — tensor anchors, frozen projected baseline, exploratory measurements

**Files:**
- Create: `invz_tensor/tests/invzt_anchors.m`, `invz_tensor/tests/fixtures/projected_baseline.json`, `invz_tensor/tests/exploratory/explore_tensor_blocks.m`, `invz_tensor/tests/exploratory/freeze_projected_baseline.m`, `invz_tensor/cache/.gitignore` (content: `*`)
- Modify: `docs/ODD-LOG.md` (controller appends §A0 preflight entry)
- Test: none yet.

**Interfaces:**
- Produces: `A = invzt_anchors()` — tensor-owned full-precision literals measured on THIS tree: `A.odd_onaxis_smallq` (`.q = [1e-1 1e-2 1e-3]` along [q 0 0], `.maxca` per q, dpRng 30 — LINEAR decay), `A.odd_tilted_limit` (q·[1 0 1]/√2 at q = 1e-3 — non-decaying, vs `4*pi*C.gfac/ion.Vc ≈ 3.7` μeV scale), `A.odd_generic_q_max` (q = [0.31 0.17 0.09]), `A.jcc0`/`A.jaa0` (Γ info at dpRng 30). Provenance comments (script, date, git hash, dpRng). CORE tests read these; interop tests may read `invz_odd_anchors()` directly.
- Produces: `projected_baseline.json` — frozen projected reference numbers for CORE-suite consumption WITHOUT the projected path: at minimum {git hash, FILTERED dirty flag, date, grid conv + n + dpRng} plus: `Jnu_flat` on the legacy_inclusive 8³/dpRng-15 grid, `Sigma0`/`crit` of `invz_solve_point` at (1.6 K, 0.5 T) on that grid with `Jxx0 = info.Jaa0` recorded, `Bc(1.2 K)` from `invz_critical`, `invz_odd_zero_field` 'full' + 'off' outputs on the legacy 8³/dpRng-15 grid (Tc, Sc_at_Tc, d_at_Tc). Generated by `freeze_projected_baseline.m` (dev-time script, projected on path; NEVER run by the suites).
- **Clean-state semantics (v3, LOCKED):** the generator and its output are themselves uncommitted while producing the JSON, so a raw repo dirty flag is permanently true and meaningless. Instead: (i) the generator `freeze_projected_baseline.m` is committed FIRST in a clean intermediate commit (Step 3a), and (ii) it records a FILTERED dirty flag = `git status --porcelain` restricted to PHYSICS-SOURCE paths only (`invz_projected/*.m`, `invz_common/*.m`, root `*.m`, `MF_dipole.m`/`exchange.m`/`qVec_generator.m`), EXCLUDING the fixture-generation paths (`invz_tensor/tests/fixtures/`, `invz_tensor/tests/exploratory/`) and `invz_tensor/cache/`. A nonempty filtered status means the physics inputs changed under the recorded git hash → STOP (the baseline would not be reproducible). The JSON stores both the raw hash and the filtered-clean boolean so a later reader can detect physics drift while tolerating the always-dirty scaffolding.

- [ ] **Step 1: Write + run `explore_tensor_blocks.m`** (uses `invz_odd_blocks` via an explicit addpath inside the script — dev-time only): on-axis rays, tilted ray, generic q at dpRng 30; smallest-shell 16³ values along a*/c* with dpRng 20/30/40 sensitivity.
- [ ] **Step 2: Write `invzt_anchors.m`** from the printed values (`format long g`).
- [ ] **Step 3a: Write `freeze_projected_baseline.m` and commit it FIRST** (clean intermediate commit `test(invzt): add projected-baseline generator (pre-freeze)`), so the git hash the JSON records corresponds to a tree whose physics sources + generator are committed.
- [ ] **Step 3b: Run `freeze_projected_baseline.m`**; it records `git rev-parse HEAD`, the FILTERED dirty flag (physics-source paths only, per the clean-state box), and the date inside the JSON; commit the JSON it writes. If the filtered flag is nonempty, STOP (physics inputs uncommitted).
- [ ] **Step 4: Controller report** — append §A0 to `docs/ODD-LOG.md`: anchors, dpRng sensitivity, baseline JSON provenance, suite baselines (projected fast + slow run once here).
- [ ] **Step 5: Commit** — `docs(invzt): A0 preflight — tensor anchors, frozen projected baseline, ODD-block exploration`.

---

### Task 3: `invzt_qgrid.m` + `invzt_jq_tensor.m` — grid contract and cached coupling tensor

**Files:**
- Create: `invz_tensor/invzt_qgrid.m`, `invz_tensor/invzt_jq_tensor.m`, `invz_tensor/tests/test_invzt_jq_tensor.m`, `invz_tensor/tests/interop/test_invzt_jq_parity.m`
- Reference (read, do not modify): `invz_projected/invz_jq_modes.m`, `invz_projected/invz_odd_blocks.m` (mirror geom-priming/per-q/cache skeleton)

**Interfaces:**
- Produces: `g = invzt_qgrid(n, conv)` per the LOCKED grid contract. For 'legacy_inclusive' it must reproduce `qVec_generator` + Γ-filter EXACTLY (verified in tests via `evalc`-wrapped comparison); weights uniform over kept points, `g.note` records the dropped-Γ and (for legacy) duplicated-face conventions.
- Produces: `[lat] = invzt_jq_tensor(ion, g, opts)` returning the LOCKED `lat` struct (`lat.Jt` `[12,12,nq]` Hermitian per page, `lat.qvec = g.qvec`, `lat.w = g.w`, `lat.conv = g.conv`, `lat.JtGamma` `[12,12]` built internally at q = [0 0 0], `lat.info` with `Jcc0/Jaa0/Jcc0_dipole/Jaa0_dipole/lorz/dpRng`). `opts`: `dpRng` 30, `cache` true, `parts` 'full' | 'dipole' (dipole-only never cached). Assembly and Γ-Lorentz placement exactly as the projected pair of functions; demag guard `invzt:demag` first. Also accepts a raw `[nq,3]` qvec in place of `g` for single-point/test use (then `lat.conv = 'explicit'`, `lat.w` uniform).

- [ ] **Step 1: Write the failing CORE tests** (`test_invzt_jq_tensor.m`, core boilerplate):

```matlab
function test_qgrid_conventions(testCase)
g1 = invzt_qgrid(8, 'halfopen');
verifyEqual(testCase, size(g1.qvec, 1), 8^3 - 1);          % Gamma dropped
verifyEqual(testCase, sum(g1.w), 1, 'AbsTol', 1e-14);
verifyLessThan(testCase, max(abs(g1.qvec(:))), 1 - 1/16);  % half-open: no +0.5 face duplicate range
g2 = evalc_call(@() invzt_qgrid(8, 'legacy_inclusive'));
[qv, ~, ~] = evalc_call(@() qVec_generator(invz_ion().a, 'mode', 'grid', ...
    'grid', [8 8 8], 'range', [-0.5 0.5], 'verbose', false));
qv = qv(any(abs(qv) > 1e-12, 2), :);
verifyEqual(testCase, g2.qvec, qv, 'AbsTol', 1e-15);       % EXACT legacy reproduction
verifyFalse(testCase, strcmp(g1.hash, g2.hash));
end

function test_shape_hermiticity_conjsym(testCase)
ion = invz_ion();
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09; -0.25 0 0; 0 0 -0.25; -0.31 -0.17 -0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
verifyEqual(testCase, size(lat.Jt), [12 12 6]);
for iq = 1:6
    verifyLessThan(testCase, norm(lat.Jt(:,:,iq) - lat.Jt(:,:,iq)', 'fro'), 1e-13);
end
for iq = 1:3
    verifyLessThan(testCase, norm(lat.Jt(:,:,iq+3) - conj(lat.Jt(:,:,iq)), 'fro'), 1e-12);
end
end

function test_gamma_lorentz_exact_placement(testCase)
% full - dipole = exchange + Lorentz. Subtract the INDEPENDENTLY built exchange
% (root function) and require the remainder == lorz on Cartesian diagonals of
% every pair, 0 elsewhere. (The old ">= lorz" form fails: J12 < 0 is AFM.)
ion = invz_ion();  C = invz_const();
lorz = 4*pi/(3*ion.Vc)*C.gfac;
latF = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false));
latD = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false, 'parts', 'dipole'));
[ex, ~] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
R = latF.Jt - latD.Jt;                                   % exchange + Lorentz
for s = 1:4, for sp = 1:4
    blk = R(3*(s-1)+(1:3), 3*(sp-1)+(1:3)) - sign(ion.J12)*ex(:,:,s,sp);
    verifyEqual(testCase, blk, lorz*eye(3), 'AbsTol', 1e-14);
end, end
end

function test_onaxis_smallq_odd_decay(testCase)
ion = invz_ion();  A = invzt_anchors();
q = [1e-1 0 0; 1e-2 0 0; 1e-3 0 0];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 30, 'cache', false));
m = arrayfun(@(iq) max(abs(lat.Jt(3:3:12, 1:3:12, iq)), [], 'all'), 1:3);
verifyEqual(testCase, m(:), A.odd_onaxis_smallq.maxca(:), 'RelTol', 1e-6);
verifyEqual(testCase, m(2)/m(1), 0.1, 'RelTol', 0.3);    % ~linear decay
end

function test_cache_roundtrip_selfverifying(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
l1 = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', true));
l2 = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', true));
verifyEqual(testCase, l2.Jt, l1.Jt);
ion2 = ion;  ion2.J12 = ion.J12 * 1.05;
l3 = invzt_jq_tensor(ion2, q, struct('dpRng', 10, 'cache', true));
verifyFalse(testCase, isequal(l3.Jt, l1.Jt));
cdir = fullfile(fileparts(mfilename('fullpath')), '..', 'cache');
verifyTrue(testCase, ~isempty(dir(fullfile(cdir, 'jqt*.mat'))));
end
```

(`evalc_call` is a small local helper wrapping `evalc` so printing functions stay quiet while returning their outputs; write it against MATLAB's `[T, out1, ...] = evalc('expr')` form.)

- [ ] **Step 2: Write the failing INTEROP tests** (`interop/test_invzt_jq_parity.m`, interop boilerplate):

```matlab
function test_block_parity_with_projected_odd_blocks(testCase)
% ca/cb/cc block equality with invz_odd_blocks — inherits its DS2023 geometry
% guards transitively. Exact same assembly path expected.
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
[Vca, Vcb, Vcc, infoS] = invz_odd_blocks(ion, q, struct('dpRng', 15, 'cache', false));
for iq = 1:3
    verifyEqual(testCase, lat.Jt(3:3:12, 1:3:12, iq), Vca(:,:,iq), 'AbsTol', 1e-15);
    verifyEqual(testCase, lat.Jt(3:3:12, 2:3:12, iq), Vcb(:,:,iq), 'AbsTol', 1e-15);
    verifyEqual(testCase, lat.Jt(3:3:12, 3:3:12, iq), Vcc(:,:,iq), 'AbsTol', 1e-14);
end
verifyEqual(testCase, lat.info.Jcc0, infoS.Jcc0, 'RelTol', 1e-12);
verifyEqual(testCase, lat.info.Jaa0, infoS.Jaa0, 'RelTol', 1e-12);
end

function test_cc_eigen_parity_with_jq_modes(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
[Jnu, ~] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
for iq = 1:2
    Jcc = (lat.Jt(3:3:12, 3:3:12, iq) + lat.Jt(3:3:12, 3:3:12, iq)')/2;
    verifyEqual(testCase, sort(real(eig(Jcc))).', Jnu(iq,:), 'AbsTol', 1e-12);
end
end
```

- [ ] **Step 3: Run both test files** — expected FAIL (functions undefined).
- [ ] **Step 4: Implement `invzt_qgrid` then `invzt_jq_tensor`** per the contracts (cache via `invz_cache_key('jqt1', dpRng, [qvec(:); convcode], pkey)`).
- [ ] **Step 5: CORE suite green (with `invz_projected` NOT on the path — verify by running the core suite from a clean MATLAB session); INTEROP suite green.**
- [ ] **Step 6: Commit** — `feat(invzt): q-grid contract + cached [12,12,nq] lattice struct with exact Lorentz/exchange placement guards`.

---

### Task 4: `invzt_chi_rpa.m` + `invzt_gcc_lattice.m` — tensor RPA with identity and Schur gates

**Files:**
- Create: `invz_tensor/invzt_chi_rpa.m`, `invz_tensor/invzt_gcc_lattice.m`, `invz_tensor/tests/test_invzt_chi_rpa.m`, `invz_tensor/tests/interop/test_invzt_rpa_parity.m`
- Reference (read, do not modify): `invz_projected/invz_chi_tensor_ref.m` (mirror its internal `si` options exactly in the interop parity test)

**Interfaces:**
- Produces: `X = invzt_chi_rpa(chi0, Jt)` (`chi0` `[3,3]` one frequency → `X` `[12,12,nq]`, `X(:,:,iq) = (eye(12) - X0*Jt(:,:,iq)) \ X0`, `X0 = kron(eye(4), chi0)`, page-wise). `[Gcc, diag4] = invzt_gcc_lattice(chi0, lat)` (`chi0` `[3,3,nz]` → weighted BZ + sublattice average of the site-diagonal cc entries using `lat.w`; `diag4` `[4,nz]`).

- [ ] **Step 1: Write the failing CORE tests**:

```matlab
function test_decoupled_cc_closes_to_scalar_branches(testCase)
% Cartesian-diagonal chi0 + Cartesian-diagonal-only Jt: the cc sector must equal
% scalar RPA over the cc-block eigenvalues COMPUTED LOCALLY from Jt (no
% projected dependency).
ion = invz_ion();  T = 1.6;
q = [0.25 0 0; 0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
c0 = invz_chi0z(si, T, 1i*0.05, struct('elastic', true));
c0d = diag(diag(c0(:,:,1)));
Jz = zero_offdiag_blocks(lat.Jt);
X = invzt_chi_rpa(c0d, Jz);
x0cc = real(c0d(3,3));
for iq = 1:2
    Jcc = (lat.Jt(3:3:12, 3:3:12, iq) + lat.Jt(3:3:12, 3:3:12, iq)')/2;
    Jnu = sort(real(eig(Jcc)));
    Xcc = X(3:3:12, 3:3:12, iq);
    got = sort(real(eig((Xcc + Xcc')/2)));
    verifyEqual(testCase, got, x0cc ./ (1 - Jnu * x0cc), 'RelTol', 1e-10);
end
end

function test_gcc_lattice_consistency(testCase)
% Brute-force equality of the weighted average (exact identity). S4 sublattice
% EQUALITY is NOT asserted here — 3 arbitrary q points are not a symmetry-
% complete set; the S4 check lives in the solver tests on full grids (report).
ion = invz_ion();  T = 1.6;
q = [0.25 0 0; 0.31 0.17 0.09; 0.1 0.2 0.3];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', false));
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
zn = 1i*[0.01 0.05 0.2];
c0 = invz_chi0z(si, T, zn, struct('elastic', true));
[Gcc, diag4] = invzt_gcc_lattice(c0, lat);
for k = 1:3
    X = invzt_chi_rpa(c0(:,:,k), lat.Jt);
    brute = zeros(4,1);
    for s = 1:4
        brute(s) = sum(lat.w(:).' .* squeeze(real(X(3*(s-1)+3, 3*(s-1)+3, :))).');
    end
    verifyEqual(testCase, diag4(:,k), brute, 'RelTol', 1e-12);
    verifyEqual(testCase, Gcc(k), mean(brute), 'RelTol', 1e-12);
end
end

function test_schur_complement_equals_E1_direct(testCase)
% THE convention gate: transverse-transverse blocks zeroed + synthetic diagonal
% chi0 -> Schur elimination of (a,b) equals scalar RPA with Jcc + dJpre, where
% dJpre is built DIRECTLY (xperp*(Vca*Vca' + Vcb*Vcb')). No invz_odd_deltaJ
% call: its caller contract requires a full uniform BZ mesh, which one generic
% q is not (v2 review finding 6a).
ion = invz_ion();
q = [0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
xperp = 11.05;  xcc = 40.0;
Jz = lat.Jt;
for s = 1:4, for sp = 1:4
    Jz(3*(s-1)+(1:2), 3*(sp-1)+(1:2), 1) = 0;
end, end
X0 = kron(eye(4), diag([xperp, xperp, xcc]));
X = (eye(12) - X0*Jz(:,:,1)) \ X0;
Xcc = X(3:3:12, 3:3:12);
Vca = Jz(3:3:12, 1:3:12, 1);  Vcb = Jz(3:3:12, 2:3:12, 1);
Jcc = Jz(3:3:12, 3:3:12, 1);
dJpre = xperp*(Vca*Vca') + xperp*(Vcb*Vcb');
Xcc_e1 = xcc * ((eye(4) - (Jcc + dJpre)*xcc) \ eye(4));
verifyEqual(testCase, Xcc, Xcc_e1, 'RelTol', 1e-10);
end

function test_full_schur_enhancement_reported(testCase)
% Exact Schur with real transverse-transverse blocks kept; E1 is the
% Jperp,perp -> 0 limit. REPORT the ratio (finite, < 1).
ion = invz_ion();
q = [0.31 0.17 0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
Xp = invz_chiperp(ion, 1.53, [0 0 0], struct());        % invz_common function
P = zeros(8);  Vc = zeros(4, 8);
pidx = @(s) 2*(s-1)+(1:2);
for s = 1:4, for sp = 1:4
    P(pidx(s), pidx(sp)) = lat.Jt(3*(s-1)+(1:2), 3*(sp-1)+(1:2), 1);
    Vc(s, pidx(sp)) = lat.Jt(3*(s-1)+3, 3*(sp-1)+(1:2), 1);
end, end
Jcc = lat.Jt(3:3:12, 3:3:12, 1);
S = Jcc + Vc * ((kron(eye(4), inv(Xp)) - P) \ Vc');
dJpre = Vc * kron(eye(4), Xp) * Vc';
r = norm(S - (Jcc + dJpre), 'fro') / max(norm(dJpre, 'fro'), 1e-30);
fprintf('full-Schur vs E1 correction ratio: %.3g (expect ~ Xp*Jaa0 ~ 0.05 scale)\n', r);
verifyTrue(testCase, isfinite(r) && r < 1);
end
```

- [ ] **Step 2: Write the failing INTEROP test** (`interop/test_invzt_rpa_parity.m`): the Σ=0 uniform-mode parity vs `invz_chi_tensor_ref` — build `Jd = kron(ones(4)/4, diag([ion.Jxx0, ion.Jxx0, ion.J0eff]))`, project the 12×12 RPA with `u = kron(ones(4,1)/2, eye(3))`, compare `imag(Xu(3,3))` against `R.chi_ten` (RelTol 1e-8). READ `invz_chi_tensor_ref.m` FIRST and mirror its si options exactly in a local wrapper that self-checks against `R.chi_sc`.
- [ ] **Step 3: Run to verify failure. Step 4: Implement both functions. Step 5: CORE + INTEROP green.** Controller logs the full-Schur ratio to ODD-LOG §A0.
- [ ] **Step 6: Commit** — `feat(invzt): page-wise 12x12 tensor RPA with direct Schur/E1 gate and weighted local averaging`.

---

### Task 5: `invzt_chi0_split.m` — transition-mask decomposition with a stated elastic convention

**Files:**
- Create: `invz_tensor/invzt_chi0_split.m`, `invz_tensor/tests/test_invzt_chi0_split.m`
- Reference (read, do not modify): `invz_common/invz_chi0z.m`

**Interfaces:**
- Produces: `[cdom, crest, mspec] = invzt_chi0_split(si, T, z, opts)` with `cdom + crest == invz_chi0z(si, T, z, opts_pass)` EXACTLY (crest = full − cdom by construction). Mask: transition (a,b) DOMINANT iff both `si.E(a) < Esplit` and `si.E(b) < Esplit` (`Esplit` default 0.4653 meV).
- **ELASTIC CONVENTION (LOCKED — v2 addition):** the z≈0 elastic term of `cdom` is the elastic expression of `invz_chi0z` restricted to dominant-group degenerate pairs, with the connected counterterm computed from DOMINANT-GROUP means: `−⟨Jμ⟩_dom⟨Jν⟩_dom`, `⟨J⟩_dom = Σ_{p∈dom} P_p M_pp`. The residual (full counterterm minus dominant counterterm, plus excited-pair elastic weight) lands in `crest` by construction. This is a convention, not physics; its size is reported as `mspec.elastic_conv_share` = |cdom_elastic convention-dependence| / |full elastic| (bounded by the excited-state population share, ~1e-3 at 1.6 K). Document in the header.
- `mspec`: `ndom` (16 with hyp, 2 without), `fdom_cc0`, `fdom_perp0`, `elastic_conv_share`.

- [ ] **Step 1: Write the failing tests** (core boilerplate):

```matlab
function test_split_exact_and_groupsize(testCase)
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
z = [0, 1i*0.05, 1i*0.5, 0.1 + 1i*5e-3];
full = invz_chi0z(si, T, z, struct('elastic', true));
[cdom, crest, mspec] = invzt_chi0_split(si, T, z, struct());
verifyEqual(testCase, cdom + crest, full, 'AbsTol', 1e-14*max(abs(full(:))));
verifyEqual(testCase, mspec.ndom, 16);
si0 = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', false));
[~, ~, m0] = invzt_chi0_split(si0, T, 0, struct());
verifyEqual(testCase, m0.ndom, 2);
end

function test_dominant_sector_is_longitudinal(testCase)
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
[~, crest, mspec] = invzt_chi0_split(si, T, 0, struct());
verifyGreaterThan(testCase, mspec.fdom_cc0, 0.90);
verifyLessThan(testCase, mspec.fdom_perp0, 0.10);
verifyGreaterThan(testCase, real(crest(1,1)), 0);
verifyLessThan(testCase, mspec.elastic_conv_share, 0.02);   % convention residual small
end

function test_esplit_insensitivity_band(testCase)
% v3 amend (Task-5 review): Esplit must sit ABOVE the full ground hyperfine
% manifold (0 -> ~0.13 meV) and BELOW the first excited CF level (~0.934 meV =
% 10.84 K). The original 0.1 was INSIDE the ground manifold, so it excluded the
% hyperfine levels in (0.1, 0.13) and genuinely changed the dominant group — a
% CORRECT split then differs at 0.1 vs 0.7. Both probes must lie in the
% (0.13, 0.93) insensitivity band; corrected 0.1 -> 0.2.
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
c1 = invzt_chi0_split(si, T, 1i*0.05, struct('Esplit', 0.2));
c2 = invzt_chi0_split(si, T, 1i*0.05, struct('Esplit', 0.7));
verifyEqual(testCase, c1, c2, 'AbsTol', 1e-15*max(abs(c1(:))));
end
```

- [ ] **Step 2: Run to verify failure. Step 3: Implement** (masked copy of the chi0z loop; header cites the convention box + Jensen's dominant-transition rule). **Step 4: CORE suite green.**
- [ ] **Step 5: Commit** — `feat(invzt): transition-mask chi0 decomposition with locked elastic-counterterm convention`.

---

### Task 6: `invzt_solve_point.m` — the A1 bridge point solver (lat-threaded)

**Files:**
- Create: `invz_tensor/invzt_solve_point.m`, `invz_tensor/tests/test_invzt_solve_point.m`, `invz_tensor/tests/interop/test_invzt_solve_parity.m`
- Reference (read, do not modify): `invz_projected/invz_solve_point.m` (loop skeleton, option names, pt contract), `jensen_1z_framework.html` §5, §7, §11.8 (the locked K relation)

**Interfaces:**
- Consumes: `lat` (T3), `invzt_gcc_lattice` (T4), `invzt_chi0_split` (T5), common engine. A1 needs NO `invz_chiperp`/`invz_odd_deltaJ` — the 12×12 RPA carries the transverse mediation automatically (and retarded).
- Produces: `pt = invzt_solve_point(ion, T, B, lat, opts)`. `opts`: `mode` 'a1' (T9/T12 add 'a2'/'a3'), `Ecut` 40, `hyp` true, `transverse_mf` 'legacy_x', `mix_outer` 0.7, `tol_outer` 1e-8, `max_outer` 60, `Sigma_seed`, `odd` (true; false zeroes Cartesian-off-diagonal entries of every block of a COPY of `lat.Jt`), `chi_rest` (true), `Esplit`, `chi0_diag` (false; TEST HOOK: zero cross-Cartesian elements of the local tensor before use — enables exact-identity gates).
- **Jaa0 threading (v2, LOCKED):** the single-ion solve receives `Jxx0 = lat.info.Jaa0` (not the `ion.Jxx0` fallback); `pt` records the value used. This makes interop parity apples-to-apples.
- **A1 K-BOOKKEEPING (LOCKED, framework §11.8):** `G0til = -(cdom_cc/(1+Sigma) + crest_cc)` site-local; `Gloc = -[weighted site-diagonal cc average]`; per frequency
  `K = 1./Gloc - 1./G0til`.
  Decoupled limit reduces EXACTLY to `invz_emt_scalar`'s relation (`K = 1./G - (1+Sigma)./G0`). Derivation in the header.
  ALGORITHM: (1) once: `[wn, wts, beta] = invz_matsubara(T, Ecut)`; `si = invz_single_ion(ion, T, B, struct('hyp', opts.hyp, 'transverse_mf', opts.transverse_mf, 'Jxx0', lat.info.Jaa0))`; `[cdom, crest] = invzt_chi0_split(si, T, 1i*wn, ...)`; `tl = invz_twolevel(ion, T, B, struct('Jxx0', lat.info.Jaa0, 'transverse_mf', opts.transverse_mf))` — NOTE: `invz_twolevel` throws `invz:degenerateDoublet` below 1e-4 meV splitting, so mode 'a1' REQUIRES a symmetry-breaking TRANSVERSE field; the solver errors `invzt:a1ZeroField` with a clear message when the transverse component is too small (v3: guard the transverse magnitude, NOT `norm(B)` — a purely longitudinal `[0 0 Bz]` has nonzero norm but does not split the doublet): for the supported `transverse_mf = 'legacy_x'` guard `abs(B(1)) < 1e-6` T; a future vector transverse mode guards `hypot(B(1), B(2)) < 1e-6` T (zero-field physics belongs to the closed form and to A3); `g = real(invz_g(tl, 1i*wn))`; (2) loop: `ctil(:,:,n) = cdom(:,:,n)/(1+Sigma(n)) + crest(:,:,n)` (after the `chi0_diag` hook if set) → `[Gcc_chi, diag4] = invzt_gcc_lattice(ctil, lat_eff)` → K per the boxed rule → `lam = invz_lambdas(K, g, wts, beta, [1 2])` → `sg = invz_sigma(tl, lam, K, g, beta)` → damped mix on Σ.
  Outputs: `Sigma0`, `Sigma`, `alpha`, `lambda`, `K`, `G` (= Gloc), `tl`, `si`, `chi0cc0`, `crit`, `sumrule_rel`, `converged`, `outer_iters`, `diag4_spread` (REPORT field — computed on full grids only), `mode`, `odd`, `lat` provenance (qvec hash, conv, `JtGamma`, `info`, `Jxx0_used`).
  CRIT (v3, Hermitian eigendecomposition — NO `sqrtm`): with `ctil0` the static renormalized 3×3 (elastic on) and `C12 = kron(eye(4), ctil0)` PSD, form the active-subspace PSD square root by Hermitian eigendecomposition, NOT `sqrtm` (which injects avoidable complex/non-Hermitian noise near criticality): `[U, D] = eig((C12+C12')/2)`; clip `d = real(diag(D))` with `d(d < rank_tol) = 0` (rank tol `1e-12`, clipped mass + active rank recorded in `pt`); `S = U*diag(sqrt(max(d,0)))*U'`. Then `M = eye(size(S,1)) - S*JtGamma*S`; `pt.crit = min(real(eig((M+M')/2)))`. Hermitian by construction, and shares the zeros of `I − C12·JtΓ` on the active subspace. (`lat.JtGamma`'s ODD blocks vanish by C2 symmetry — assert once `< 1e-10*Jcc0`.)
  SUM RULE: `sumrule_rel = |sum(wts.*G)/beta + si.JzJz_fluct| / si.JzJz_fluct` (report-quality).

- [ ] **Step 1: Write the failing CORE tests**:

```matlab
function test_zero_field_guard(testCase)
ion = invz_ion();
lat = invzt_jq_tensor(ion, [0.25 0 0], struct('dpRng', 10, 'cache', false));
verifyError(testCase, @() invzt_solve_point(ion, 1.6, [0 0 0], lat, struct()), ...
    'invzt:a1ZeroField');
end

function test_odd_on_converges_crit_direction_reported(testCase)
% v3 (review Other 8): the EXACT PSD Schur relation (dJpre = xperp*(Vca*Vca' +
% Vcb*Vcb') >= 0) is the HARD gate, and it lives at FIXED chi0 in Task 4
% (test_schur_complement_equals_E1_direct / test_full_schur_enhancement_reported).
% Here the FULL odd-on solve also moves the self-consistent medium and scalar
% Sigma, so monotonicity of the final pt.crit does NOT follow from the Gaussian
% Schur identity alone. Hard-assert convergence + sum rule; MEASURE the crit
% direction (expected p1.crit < p0.crit; a positive shift is a finding to
% investigate, not an automatic pass/fail here).
ion = invz_ion();  T = 1.6;  B = [0.1 0 0];
g = invzt_qgrid(8, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
p1 = invzt_solve_point(ion, T, B, lat, struct('odd', true));
p0 = invzt_solve_point(ion, T, B, lat, struct('odd', false));
verifyTrue(testCase, p1.converged && p0.converged);
verifyLessThan(testCase, p1.sumrule_rel, 0.10);
verifyTrue(testCase, isfinite(p1.crit) && isfinite(p0.crit));
fprintf('A1 odd on/off (crit direction REPORTED): crit %.5f / %.5f (d=%+.2e), Sigma0 %.4f / %.4f, diag4 spread %.2e\n', ...
    p1.crit, p0.crit, p1.crit - p0.crit, p1.Sigma0, p0.Sigma0, p1.diag4_spread);
end

function test_chi_rest_toggle_reported(testCase)
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 10, 'cache', true));
pa = invzt_solve_point(ion, T, B, lat, struct('chi_rest', true));
pb = invzt_solve_point(ion, T, B, lat, struct('chi_rest', false));
verifyTrue(testCase, pa.converged && pb.converged);
fprintf('chi_rest on/off: crit %.5f / %.5f, Sigma0 %.4f / %.4f\n', pa.crit, pb.crit, pa.Sigma0, pb.Sigma0);
end

function test_no_odd_vs_frozen_projected_baseline(testCase)
% PHYSICAL no-ODD comparison against the FROZEN baseline (core suite: no
% projected path). Category: physical full-tensor no-ODD mode; residual sources:
% cross-Cartesian chi0 elements, dominant-mask vs whole-cc division. Gate 2e-3.
fx = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'projected_baseline.json')));
ion = invz_ion();
g = invzt_qgrid(8, 'legacy_inclusive');                  % SAME convention as the baseline
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
pt = invzt_solve_point(ion, fx.point.T, [fx.point.Bx 0 0], lat, ...
    struct('odd', false, 'chi_rest', false));
verifyTrue(testCase, pt.converged);
verifyEqual(testCase, pt.Sigma0, fx.point.Sigma0, 'AbsTol', 2e-3);
verifyEqual(testCase, sign(pt.crit), sign(fx.point.crit));
fprintf('A1 no-ODD vs frozen projected: dSigma0 = %.3e\n', pt.Sigma0 - fx.point.Sigma0);
end
```

- [ ] **Step 2: Write the failing INTEROP test** (`interop/test_invzt_solve_parity.m`): live version of the frozen-baseline comparison (same grid convention both sides, `'Jxx0', info.Jaa0` passed to the projected leg) plus the odd-on comparison vs the projected Tier-1 solver (blocks via `invz_odd_blocks` on the SAME legacy grid; REPORT `dSigma0`, `dcrit`; assert finiteness).
- [ ] **Step 3: Run to verify failure. Step 4: Implement** per the boxed algorithm (if the frozen-baseline gate cannot be met, STOP and report the measured `dSigma0`). **Step 5: CORE + INTEROP green** (runtime: solver tests < 90 s warm-cache). Controller logs to ODD-LOG §A1.
- [ ] **Step 6: Commit** — `feat(invzt): A1 bridge point solver — lat-threaded, corrected K bookkeeping, zero-field guard, frozen-baseline gate`.

---

### Task 7: Critical finders — Bc bisection, tensor-owned Tc extrapolator, zero-field parity

**Files:**
- Create: `invz_tensor/invzt_critical.m`, `invz_tensor/invzt_tc_pm_extrap.m`, `invz_tensor/tests/test_invzt_critical.m`, `invz_tensor/tests/interop/test_invzt_critical_parity.m`
- Reference (read, do not modify): `invz_projected/invz_critical.m`, `invz_projected/invz_crit_at.m`, `invz_projected/invz_odd_zero_field.m`, `invz_projected/tests/invz_odd_tc_pm_extrap.m` (PATTERN only — its signature `(ion, B, Jf, o, Tg)` is hard-wired to `invz_crit_at` and CANNOT be reused here)

**Interfaces:**
- Produces: `[Bc, out] = invzt_critical(ion, T, lat, Brange, opts)` — PM-side bisection on `pt.crit` (warm-start `Sigma_seed` threaded; non-convergence = phase signal).
- Produces: `[tc, used] = invzt_tc_pm_extrap(critfun, Tg)` — tensor-owned, HANDLE-BASED: `critfun(T) -> [crit, ok]`; samples the fixed grid `Tg` (house step 1/30 K), keeps `ok` points, linearly extrapolates the two lowest converged PM points to crit = 0 (the projected fixture's MATH, re-owned with a handle interface). Errors `invzt:tcNoWindow` if < 2 converged points.
- **Zero-field policy (v3, LOCKED — resolves the ambiguity flagged in review Other 4):** EVERY mode-switched tensor Tc (a1/a2/a3) uses the SMALL-B PROXY `B = [0.05 0 0]` T (the projected T3.3/T3.4 proxy convention) — never B = 0 (Task 6 guard). Tensor A3 true-zero-field Tc is **DEFERRED**: only the projected closed-form chain (the parity test below) is truly at B = 0. This is honest about scope — although the A3 divided-difference kernels are degeneracy-regular by construction (constraint 5, initial-state weights `p_r`, no `D_pq/(ε_q−ε_p)`), so a future true-B=0 A3 path is feasible, it is not built or validated in this plan; the small-Bx proxy stands in for every tensor Tc quote.

- [ ] **Step 1: Write the failing CORE tests**:

```matlab
function test_tc_extrap_handle_math(testCase)
% Synthetic crit(T): linear through Tc = 1.5 with a non-converged hole.
critfun = @(T) deal(0.4*(T - 1.5), T > 1.52 || T > 1.6);
f = @(T) wrap2(critfun, T);   % local helper returning [c, ok]
Tg = 1.4:1/30:1.9;
tc = invzt_tc_pm_extrap(f, Tg);
verifyEqual(testCase, tc, 1.5, 'AbsTol', 1e-10);
end

function test_bc_finder_structure(testCase)
ion = invz_ion();  T = 1.6;
g = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 10, 'cache', true));
[Bc, out] = invzt_critical(ion, T, lat, [0.5 5], struct('odd', false, 'tol', 0.05));
verifyTrue(testCase, isfinite(Bc) && Bc > 0.5 && Bc < 5);
verifyTrue(testCase, out.iters(end).pt.converged);
end
```

(`wrap2` local helper converts the two-output anonymous into the `[c, ok]` contract; write it explicitly in the test file.)

- [ ] **Step 2: Write the failing INTEROP tests** (`interop/test_invzt_critical_parity.m`):

```matlab
function test_zero_field_closed_form_parity(testCase)
% Tensor lat slices through the projected closed-form chain vs invz_odd_zero_field
% — BOTH on the SAME legacy_inclusive 8^3/dpRng-15 grid (v2: grid contract).
% READ invz_odd_zero_field.m first; adapt out-field names if they differ.
ion = invz_ion();
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
Vca = lat.Jt(3:3:12, 1:3:12, :);  Vcb = lat.Jt(3:3:12, 2:3:12, :);  Vcc = lat.Jt(3:3:12, 3:3:12, :);
[TcP, outP] = invz_odd_zero_field(ion, struct('mode', 'full', 'grids', {{8}}, 'dpRng', 15));
Xp = invz_chiperp(ion, TcP(1), [0 0 0], struct());
[dJ, d] = invz_odd_deltaJ(Vca, Vcb, Xp);
Jnu = invz_odd_modes(Vcc, dJ);
verifyEqual(testCase, d, outP.d_at_Tc(1), 'RelTol', 1e-9);
Sc = invz_sigma_crit(lat.info.Jcc0 - d, Jnu(:));
verifyEqual(testCase, Sc, outP.Sc_at_Tc(1), 'RelTol', 1e-9);
end

function test_bc_parity_no_odd(testCase)
ion = invz_ion();  T = 1.2;
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
[Jnu, info] = invz_jq_modes(ion, g.qvec, struct('dpRng', 15, 'cache', true));
% v3 amend (Task-7 review): chi_rest=TRUE for apples-to-apples parity. The
% projected invz_solve_point keeps the FULL local cc chi0 (no dominant/rest
% split), so chi_rest=false drops the excited-state chi0 and injects a ~0.16 T
% Bc error at the high-field Bc(1.2 K)~2.7 T (negligible only at the low-field
% T6 baseline). The rest-dropped variant is a reported split diagnostic, not the
% parity configuration.
Bc_t = invzt_critical(ion, T, lat, [2 6], struct('odd', false, 'chi_rest', true));
Bc_s = invz_critical(ion, T, Jnu(:), struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0));
fprintf('Bc(1.2 K) no-ODD: tensor %.4f T, projected %.4f T\n', Bc_t, Bc_s);
verifyEqual(testCase, Bc_t, Bc_s, 'AbsTol', 0.05);
end

function test_a1_smallB_tc_odd_slow(testCase)
% HEADLINE (report): A1 Tc at the 0.05 T proxy with ODD on, production
% 16^3/dpRng 30, via invzt_tc_pm_extrap; compare projected closed form 1.509 K
% (grid-matched legacy convention). The gap = A1 enhancements (retardation
% ~null + transverse RPA + chi_rest).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
g = invzt_qgrid(16, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
critfun = @(T) crit_ok(invzt_solve_point(ion, T, [0.05 0 0], lat, struct('odd', true)));
Tc_a1 = invzt_tc_pm_extrap(critfun, 1.30:1/30:1.75);
[TcP, ~] = invz_odd_zero_field(ion, struct('mode', 'full', 'grids', {{16}}, 'dpRng', 30));
fprintf('Tc ODD: A1 tensor(0.05 T proxy) %.4f K vs projected closed form %.4f K\n', Tc_a1, TcP(1));
verifyTrue(testCase, isfinite(Tc_a1) && Tc_a1 > 1.3 && Tc_a1 < 1.75);
end
```

(`crit_ok` local helper: `function [c, ok] = crit_ok(pt); c = pt.crit; ok = pt.converged && isfinite(c); end`.)

- [ ] **Step 3: Run to verify failure. Step 4: Implement both functions. Step 5: CORE + INTEROP green; run the slow test.** Controller logs Bc parity + the A1 proxy-Tc vs closed form (and a `chi_rest=false` rerun for the split) to ODD-LOG §A1.
- [ ] **Step 6: Commit** — `feat(invzt): PM-side Bc finder + handle-based small-B Tc extrapolator; zero-field closed-form parity (grid-matched)`.

---

### Task 8: `invzt_chi_realaxis.m` — tensor real-axis spectra (A1 continuation ONLY)

**Files:**
- Create: `invz_tensor/invzt_chi_realaxis.m`, `invz_tensor/tests/test_invzt_chi_realaxis.m`, plus one interop test appended to `interop/test_invzt_rpa_parity.m`
- Reference: `invz_projected/invz_chi_realaxis.m` (Kw seeding pattern)

**Interfaces:**
- Produces: `out = invzt_chi_realaxis(ion, T, B, pt, w, opts)` — rebuilds `ctil(w)` on `z = w + 1i*eta` (elastic OFF) with the A1 scalar continuation `Sigma_w = alpha + pref*(lambda1 - (1-n01^2)*Kw).*g(z)` (`Kw = pt.K(1)` seed), then pages `invzt_chi_rpa`. Uses `pt.lat` provenance: `qsel` 'gamma_uniform' (default; builds the Task-4 `Jd` uniform page from `pt.lat.info`) | 'gamma' (uses `pt.lat.JtGamma`) | an explicit `[nq,3]` list (rebuilds `lat` via `invzt_jq_tensor` at `opts.dpRng`). Returns `out.chi_uniform` `[3,3,nw]`, `out.chi_cc_q` `[nq,nw]` (when q-resolved), `out.Sigma_w`. `opts`: `eta` 5e-3, `odd` (must equal `pt.odd` — error `invzt:oddMismatch`), `force_sigma0` (test hook).
- **SCOPE (v2, explicit):** this is the A1 scalar-Σ analytic continuation. It does NOT extend to A3: continuing the full `Vmat(iωn)` needs either direct real-axis kernel evaluation (the 𝓘 kernels accept complex frequency arguments) or a fitted continuation — a separate future work item, documented in the header and README.

- [ ] **Step 1: Write the failing CORE tests**:

```matlab
function test_sigma0_gamma_uniform_matches_bare_rpa(testCase)
% force_sigma0 + qsel 'gamma_uniform': the uniform-mode response must equal the
% single-site 3x3 RPA with J = diag(Jaa0, Jaa0, Jcc0) computed LOCALLY (exact
% identity; no projected dependency).
ion = invz_ion();  T = 1.6;  B = [2 0 0];
w = (0.05:0.05:0.7).';
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('odd', false));
out = invzt_chi_realaxis(ion, T, B, pt, w, struct('force_sigma0', true, 'odd', false));
si = pt.si;
c0 = invz_chi0z(si, T, w + 1i*5e-3, struct('elastic', false));
J3 = diag([lat.info.Jaa0, lat.info.Jaa0, lat.info.Jcc0]);
for k = 1:numel(w)
    ref = (eye(3) - c0(:,:,k)*J3) \ c0(:,:,k);
    verifyEqual(testCase, out.chi_uniform(:,:,k), ref, 'RelTol', 1e-8);
end
end

function test_odd_on_spectra_reported(testCase)
ion = invz_ion();  T = 1.55;  B = [0.5 0 0];
w = (0.02:0.005:0.5).';
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
p1 = invzt_solve_point(ion, T, B, lat, struct('odd', true));
p0 = invzt_solve_point(ion, T, B, lat, struct('odd', false));
o1 = invzt_chi_realaxis(ion, T, B, p1, w, struct('odd', true));
o0 = invzt_chi_realaxis(ion, T, B, p0, w, struct('odd', false));
[~, i1] = max(squeeze(imag(o1.chi_uniform(3,3,:))));
[~, i0] = max(squeeze(imag(o0.chi_uniform(3,3,:))));
fprintf('ODD spectra: peak %.4f (on) vs %.4f (off) meV, shift %+.1f ueV\n', ...
    w(i1), w(i0), 1e3*(w(i1) - w(i0)));
verifyTrue(testCase, isfinite(w(i1)) && isfinite(w(i0)));
end
```

- [ ] **Step 2: Interop test** (append): no-ODD PM point peak-parity vs `invz_chi_realaxis` (peak energy within 0.01 meV; PM-side field, e.g. (1.6 K, 2 T)).
- [ ] **Step 3: Run to verify failure. Step 4: Implement. Step 5: CORE + INTEROP green.** Controller logs the first ODD spectra numbers to ODD-LOG §A1.
- [ ] **Step 6: Commit** — `feat(invzt): tensor real-axis spectra (A1 continuation) with exact gamma-uniform identity gate`.

---

### Task 9: `invzt_emt_matrix.m` — the A2 DIRECT matrix effective-medium closure

**Files:**
- Create: `invz_tensor/invzt_emt_matrix.m`, `invz_tensor/tests/test_invzt_emt_matrix.m`
- Modify: `invz_tensor/invzt_solve_point.m` (add `mode` 'a2')
- Reference: `jensen_1z_framework.html` §5; `invz_projected/invz_emt_scalar.m` (the 1-channel closed form)

**Interfaces (v2 — REWRITTEN per review finding 2):**
- The closure is DIRECT, not iterative. Substituting the impurity relation `chi_imp = (ctil⁻¹ − K)⁻¹` into the lattice form cancels K exactly, so per frequency:
  1. `chi_lat(q) = (ctil⁻¹ − J(q))⁻¹` — the bare tensor RPA of `ctil` (page solves; NEVER a literal inverse of a possibly rank-deficient `ctil`: implement as `chi_lat(q) = ctil*(I − J(q)*ctil)⁻¹`-form page solves, identical algebra, regular for singular `ctil`),
  2. `chi_bar` = weighted BZ + sublattice average of the site-diagonal 3×3 blocks (per-sublattice blocks kept; `info.persub_spread` reported),
  3. `K = ctil⁻¹ − chi_bar⁻¹` — computed via solves on ONE COMMON ACTIVE SUBSPACE (v2): rank-reveal `ctil(iwn)` at `opts.rank_tol` (1e-12) with a FREQUENCY-CONSISTENT projector (the union of active channels over the frequency grid, so near-null rotations cannot produce discontinuous K); K on the null complement is set to 0 and the projector is recorded in `info`.
  Scalar reduction check (sign LOCKED, v2): with `ctil = −G0/D` and `chi_bar = −G`, `K = 1/G − D/G0` — exactly `invz_emt_scalar`'s `med.K`, POSITIVE relation, no sign flip.
- Produces: `[K, chi_bar, info] = invzt_emt_matrix(ctil, lat, opts)`; `K`, `chi_bar` `[3,3,nz]`. **Symmetry (v3, corrected):** since `chi0(iwn)` is NON-Hermitian off the static slot (the physical gyrotropic ~B term is imaginary-antisymmetric), `K`/`chi_bar` obey the LOCKED transpose relation (constraint 9) `X(-iwn) = X(iwn).'` and are Hermitian ONLY at the static slot (wn=0) — do NOT assert single-wn Hermiticity, and do NOT symmetrize the gyrotropic part. `info`: active-subspace (frequency-consistent) projector, `persub_spread`, and the transpose-relation + static-Hermiticity residuals.
- Produces (solve_point): `mode = 'a2'` — as 'a1' except K: run `invzt_emt_matrix` on `ctil` and take `K_scalar(n) = +K_cc(n)` (LOCKED POSITIVE — v2 sign fix) into `invz_lambdas`/`invz_sigma`; store `pt.Kmat` for A3. With the `chi0_diag` hook set, 'a2' must equal 'a1' EXACTLY (matrix inversion == componentwise on diagonal input) — the exact-identity gate; with physical `ctil` (cross-Cartesian elements, e.g. the symmetry-allowed yz), the difference IS the matrix-medium content (report).

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_scalar_limit_parity_with_emt_scalar(testCase)
% Nonsingular spectator channels (v2: avoids the null-channel inverse), zero
% a/b couplings -> K_aa = K_bb = 0 exactly, K_cc == med.K (POSITIVE sign).
% S4-symmetric synthetic cc blocks: circulant eigenvectors (commute with the
% sublattice cyclic shift), so the four sites stay equivalent. v3 (review Other
% 3): KEEP the COMPLEX Hermitian DFT result — `F` is complex, so `real(F*diag*F')`
% would change the eigenvalues away from `Jnu`. `F*diag(Jnu)*F'` is complex
% Hermitian, circulant, and has eigenvalues EXACTLY `Jnu` (F unitary). The full
% 12x12 page is then complex Hermitian, which the matrix EMT handles natively.
rng(7);
nq = 40;  nz = 5;
F = dftmtx(4)/2;                                    % unitary circulant eigenbasis
Jnu = 6e-3 * (rand(nq, 4) - 0.3);
Jt = zeros(12, 12, nq);
for iq = 1:nq
    Jt(3:3:12, 3:3:12, iq) = F*diag(Jnu(iq,:))*F';  % complex Hermitian, eig == Jnu
end
lat = struct('Jt', Jt, 'qvec', zeros(nq,3), 'w', ones(nq,1)/nq, 'conv', 'explicit', ...
    'JtGamma', zeros(12), 'info', struct());
g0 = -[0.8; 1.1; 1.6; 2.2; 3.0];
Sigma = 0.1*ones(nz, 1);
ctil = zeros(3, 3, nz);
ctil(1,1,:) = 1e-3;  ctil(2,2,:) = 1e-3;            % nonsingular spectators
ctil(3,3,:) = -g0 ./ (1 + Sigma);
[K, chi_bar] = invzt_emt_matrix(ctil, lat, struct());
med = invz_emt_scalar(g0, Sigma, sort(real(eig_all(Jt))), struct());   % v3: eig_all reads ONLY the 4x4 cc blocks (not the full 12x12 — that would inject 8 spurious zero transverse branches into invz_emt_scalar)
verifyEqual(testCase, squeeze(K(3,3,:)), med.K, 'RelTol', 1e-8);       % POSITIVE sign (v2)
verifyEqual(testCase, squeeze(-chi_bar(3,3,:)), med.G, 'RelTol', 1e-8);
verifyEqual(testCase, max(abs(squeeze(K(1,1,:)))), 0, 'AbsTol', 1e-12);% decoupled spectator
end

function test_direct_closure_and_transpose_symmetry(testCase)
% Defining identities on physical input: chi_lat is the bare RPA of ctil (K
% cancels), chi_imp(K) == chi_bar on the active subspace.
% v3 amend (Task-9 review, verified empirically): chi0(iwn) is NON-Hermitian at
% wn~=0 -- the gyrotropic (~B) part is a PHYSICAL imaginary-antisymmetric term
% (measured rel anti-Hermitian part ~6-17% at Bx=0.5 T). So chi_bar and K are
% non-Hermitian OFF the static slot; they obey the LOCKED transpose relation
% (constraint 9) X(-iwn) = X(iwn).', and are Hermitian ONLY at wn=0. The v2
% asserts (chi_bar == (bz+bz')/2; K == K') contradicted constraint 9 and are
% corrected below. The gyrotropic part is NOT symmetrized away. (Tolerances here
% are magnitude-scaled sketches -- the implementer adapts to the measured norms.)
ion = invz_ion();  T = 1.6;
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
si = invz_single_ion(ion, T, [0.5 0 0], struct('hyp', true));
% One call over +wn AND -wn so the frequency-consistent active-subspace projector
% is shared (the transpose relation is exact only under a shared projector):
zset = 1i*[0.01 0.1 -0.01 -0.1];
ctil = invz_chi0z(si, T, zset, struct('elastic', true));
[K, chi_bar, info] = invzt_emt_matrix(ctil, lat, struct());
for k = 1:2                                              % positive-wn slots
    X = invzt_chi_rpa(ctil(:,:,k), lat.Jt);
    acc = zeros(3, 3, size(X, 3));
    for s = 1:4, acc = acc + X(3*(s-1)+(1:3), 3*(s-1)+(1:3), :) / 4; end
    bz = sum(acc .* reshape(lat.w, 1, 1, []), 3);
    verifyEqual(testCase, chi_bar(:,:,k), bz, 'RelTol', 1e-10);   % RAW average (NOT Hermitized)
    P = info.projector;                                 % [3 x r] orthonormal columns
    R = ctil(:,:,k) \ chi_bar(:,:,k) - K(:,:,k)*chi_bar(:,:,k);
    verifyEqual(testCase, P'*R*P, eye(size(P, 2)), 'AbsTol', 1e-8);
    % transpose relation (constraint 9): slot 2+k is -wn of slot k
    verifyEqual(testCase, K(:,:,2+k),       K(:,:,k).',       'AbsTol', 1e-9);
    verifyEqual(testCase, chi_bar(:,:,2+k), chi_bar(:,:,k).', 'AbsTol', 1e-9);
end
% off the static slot K is genuinely non-Hermitian (asserting Hermiticity here is
% the v2 bug): confirm the anti-Hermitian part is O(1), not O(1e-12):
verifyGreaterThan(testCase, norm(K(:,:,2) - K(:,:,2)'), 1e-3*norm(K(:,:,2)));
% static slot IS Hermitian (real symmetric):
c0 = invz_chi0z(si, T, 0, struct('elastic', true));
[K0, cb0] = invzt_emt_matrix(c0, lat, struct());
verifyLessThan(testCase, norm(K0(:,:,1)  - K0(:,:,1)'),  1e-9*max(norm(K0(:,:,1)), 1e-6));
verifyLessThan(testCase, norm(cb0(:,:,1) - cb0(:,:,1)'), 1e-9*max(norm(cb0(:,:,1)), 1e-6));
end

function test_a2_map_equals_a1_on_diagonal(testCase)
% v3 amend (Task-9 review): the exact identity "matrix EMT == scalar EMT on
% DIAGONAL input" is a property of the K MAP, not the converged Sigma0. At this
% chi0_diag point the T6 self-consistent solve is MULTI-ROOT (both roots are
% valid a1 fixed points; Anderson amplifies a ~1e-14 map difference into a
% Delta~0.09 root split), so a fixed-point Sigma0 equality is fragile. Instead
% verify the a2-converged root is ALSO an a1 fixed point (shared fixed-point set
% == identical map on diagonal input), which is the real content of the identity.
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
p2 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a2', 'chi0_diag', true, 'odd', false));
% seed a1 at the a2 root: a1 must STAY there (the a2 root is an a1 fixed point):
pchk = invzt_solve_point(ion, T, B, lat, ...
    struct('mode', 'a1', 'chi0_diag', true, 'odd', false, 'Sigma_seed', p2.Sigma));
verifyEqual(testCase, pchk.Sigma0, p2.Sigma0, 'AbsTol', 1e-9);
% physical (non-diagonal) a2 - a1 difference is a REPORT (the matrix-medium content):
p3 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a2', 'odd', true));
p4 = invzt_solve_point(ion, T, B, lat, struct('mode', 'a1', 'odd', true));
fprintf('A2 vs A1 (physical, odd on): dSigma0 = %+.3e, dcrit = %+.3e\n', ...
    p3.Sigma0 - p4.Sigma0, p3.crit - p4.crit);
verifyTrue(testCase, p3.converged && p4.converged);
end
```

(Local helpers `eig_all` (flatten the eigenvalues of the `Jt(3:3:12,3:3:12,:)` CC BLOCKS over pages — the 4 cc branches per q, NOT the full 12x12) and `pageinv_solve` (solve-based inverse on a subspace) are ~5 lines each; write them in the test file. The impurity-relation check may be simplified to the equivalent solve-form assert — keep the identity, adjust the linear-algebra plumbing to the implementation's projector convention.)

- [ ] **Step 2: Run to verify failure. Step 3: Implement** (direct closure; header carries the K-cancellation derivation and the scalar sign check). **Step 4: CORE suite green.** Controller logs the A2-vs-A1 physical difference to ODD-LOG §A2.
- [ ] **Step 5: Commit** — `feat(invzt): A2 direct matrix effective-medium closure (K = ctil^-1 - chibar^-1, positive scalar sign, active-subspace solves)`.

---

### Task 10: `verify_tensor_vertex.py` — the component-labelled vertex oracle

**Files:**
- Create: `verify_tensor_vertex.py` (repo root), `invz_tensor/tests/fixtures/vertex_oracle.json`
- Reference (read): `verify_v2_vertex.py` (reuse its 𝓘₂/𝓘₃ mpmath kernels and GL quadrature verbatim)

**Interfaces:** definitions LOCKED as in v1 (component-labelled path sum `F_{μν;ρσ}(n,ℓ)` with (O, ξ) triples permuted into descending time order; pairing subtraction `Γ⁽⁴⁾ = F − β·C_μν(ζn)·C_ρσ(ζℓ) − β·δ_{n,−ℓ}·C_μρ(ζn)·C_σν(ζn) − β·δ_{n,ℓ}·C_μσ(ζn)·C_ρν(ζn)`, to be NUMERICALLY pinned; contraction `V_μν(n) = (1/2β)Σ_ℓ Σ_{ρσ} K_ρσ(ℓ)·Γ⁽⁴⁾`; three-point `T³`). v2 additions to the check list:

Checks (print PASS/FAIL, exit nonzero on FAIL):
  1. `F` vs brute-force GL quadrature over the six ordered simplices; random 3-level Hermitian operators WITH diagonals; (n,ℓ) ∈ {(0,0),(0,1),(1,0),(1,1),(1,−1),(2,1)}; tol 5e-8.
  2. `Γ⁽⁴⁾` vs brute-force connected part where the pairing transforms are computed INDEPENDENTLY by quadrature (v2: never subtract the same analytic expression from both sides); tol 5e-7. This pins the delta/index assignment — if the written form fails, correct it, record the proven form in the docstring, re-run.
  3. Scalar reduction == `verify_v2_vertex`'s `Γ⁽⁴⁾_{nℓ}` to 1e-12.
  4. Two-level Jensen limit: `V_n == G0_n·Σ_n^{Jensen}` (J 2.17–2.19 inline), tol 1e-9.
  5. Degeneracy regularity: E1 = E2 doublet with off-diagonal elements; finite + quadrature agreement; tol 5e-7.
  6. **Populated-spectator check (v2):** a 3-level system at β·Δ2 ≈ 1 (spectator population ~0.2): the two-level-channel content must match the spectator-normalized two-level formula (p̃ = p/q + deficit term, three-level doc Eq. 31/32 pattern) — the REAL test of constraint 3; far-removed spectators do not exercise it.
  7. **Negative-frequency transpose (v2):** for a field-on (complex-element) system, verify `C_ρσ(−ζℓ) = C_σρ(ζℓ)` and the contracted `V_μν(−n)` relation implied by it; store rows at negative n.
  8. Physical chain anchor (DS2023 3-state params) at 1.53 K: mixed-sector rows (cc;cc), (cc;aa), (ca;ac).
  9. Sum-rule consistency (doc Eq. 47 analogue, cc): with STATED cutoff `n ≤ n_max = 64`, the ωn⁻² analytic tail estimate added, and tolerance 1e-6 relative (v2: reproducibility spec) — report residual.
  10. **Factored-identity verification (v3 — GATES Task 11's `'factored'` path; review Blocking 1).** Derive the proposed chained-resolvent identity explicitly in the script's docstring and implement BOTH the dense ordered path-sum AND the factored `p·O·R(ξ)·O·R(ξ')·O·R(ξ'')·O` form in python. Compare them PER ORDERING and PER `(n,ℓ)` on the random-complex, degenerate, and populated-spectator fixtures; tol 1e-9. Print `FACTORIZATION ESTABLISHED` (write `oracle.factored_ok = true` to the JSON) only if every ordering/(n,ℓ) matches; otherwise print `FACTORIZATION NOT ESTABLISHED`, record which orderings failed, set `factored_ok = false`, and DO NOT fail the script (the dense route is still valid — this only tells Task 11 whether it may wire `'factored'`). The partial-fraction step must handle near-degenerate denominators via the cluster Hermite-limit fallback; test that fallback at `|denominator|` down to `1e-10`.
  11. **KMS whole-value check (v3 — constraint 4; review Other 7).** At βΔ ∈ {50, 100, 200}, compare the ENTIRE population-weighted contracted value (`V_μν(n)` and representative kernel products), NOT individual exponentials, against a full-precision mpmath evaluation; the accumulation carries complex sign/phase (complex-log or explicit KMS regrouping). tol 1e-8 relative. Store these large-βΔ rows in the fixture for the MATLAB kernel tests.
Fixture: JSON with `beta`, energies, operator matrices (re/im), K-lists, kernel rows (incl. repeated-node and large-βΔ KMS cases), conventions block (index order, frequency assignment, transpose relation), the `factored_ok` boolean + per-ordering factored-vs-dense residuals (check 10), and `{tags, n, l, value_re, value_im}` rows at 1e-12; provenance header (git hash, date).

- [ ] **Step 1: Write the script** (structure from `verify_v2_vertex.py`; checks 1–11, incl. the dense+factored implementations for check 10 and the large-βΔ whole-value rows for check 11).
- [ ] **Step 2: Run:** `python3 verify_tensor_vertex.py` → `ALL TENSOR-VERTEX CHECKS PASSED` (checks 1–9, 11 are hard; check 10 records `factored_ok` without failing the run), JSON written. Iterate check 2 until the pairing form is proven; docstring records it. Record the check-10 verdict — `FACTORIZATION ESTABLISHED` or `NOT` — it decides whether Task 11 may wire `'factored'`.
- [ ] **Step 3: Hygiene:** re-run the two existing verifiers (must still pass); sanity-scan the JSON (no NaN/Inf).
- [ ] **Step 4: Commit** — `feat(invzt): component-labelled four-point vertex oracle with populated-spectator and transpose-relation checks + fixture`.

---

### Task 11: `invzt_kernels.m` + `invzt_vertex4.m` + `invzt_vertex3.m` — engine + PERFORMANCE GATE

**Files:**
- Create: `invz_tensor/invzt_kernels.m`, `invz_tensor/invzt_vertex4.m`, `invz_tensor/invzt_vertex3.m`, `invz_tensor/tests/test_invzt_vertex.m`, `invz_tensor/tests/exploratory/perf_vertex_scaling.m`

**Interfaces:**
- `k = invzt_kernels()` — `k.phi/k.I2/k.I3` exprel-stable, complex-capable, Hermite repeated-node branches, AND (v2, constraint 4) **KMS-scaled evaluation**: kernel calls accept the initial-state log-weight and return the population-weighted product computed in log-regrouped form (`exp(log p_r + accumulated exponent)` internally, never `p_r * huge`), so large-βΔ endpoint layers evaluate without overflow/underflow. Tested at βΔ up to 200 against mpmath oracle rows (fixture includes such rows).
- `V = invzt_vertex4(es, ops, Kmat, wn, beta, opts)` — as v1 (external (μν), internal channel set contracted against `Kmat` on the fly; no initial-state truncation; matrix-element pruning `opts.mtol` with recorded pruned-weight bound), with (v2) negative-ℓ reconstruction via the LOCKED transpose relation `K_ρσ(−iωℓ) = K_σρ(iωℓ)` (constraint 9). **v3 (review Blocking 1): the DENSE path is primary — no factorization is assumed.** Two implementations behind `opts.impl`:
  - `'dense'` (DEFAULT) — the literal ordered path-sum: for each of the six time-orderings, sum over the quadruple of eigenstates weighted by the KMS-scaled 𝓘₃ kernel (constraint 4 accumulation). This is BOTH the correctness reference AND the production engine for every rung whose measured cost fits the budget (perf gate below / Task 13). Cost ~N⁴ per ordering; that is accepted, not optimized away by assumption.
  - `'factored'` (OPT-IN, UNPROVEN until derived) — the O(N³) chained-resolvent form `p·O·R(ξ)·O·R(ξ')·O·R(ξ'')·O` is a CONJECTURE, not an established identity: the sequential 𝓘₃ denominators `𝓘₃(E_r−E_s+ξ₁, E_s−E_t+ξ₂, E_t−E_u+ξ₃)` and the initial-state weight `p_r` live in transition/Liouville space and do NOT obviously reduce to ordinary N×N state-space resolvents. It becomes available ONLY after Task 10 check 10 prints `FACTORIZATION ESTABLISHED` (the identity derived in the verifier docstring and shown to reproduce the dense result per ordering / per (n,ℓ) on the random-complex, degenerate, and populated-spectator oracle systems to 1e-9); near-degenerate clusters (|denominator| < `opts.dd_tol`) route through a dense Hermite-limit fallback restricted to the cluster. Until that proof lands `'factored'` is DISABLED (selecting it errors `invzt:factoredUnproven`) and the O(N³) scaling is NOT assumed anywhere in the plan.
- `[dm, T3] = invzt_vertex3(...)` — three-point engine (sum-rule route), same kernel/KMS machinery.
- **PERFORMANCE GATE (v3 — blocking; DENSE is the costed path):** `perf_vertex_scaling.m` measures wall time of one full `Vmat` assembly (all 6 external pairs) with `'dense'` at N = 3, 6, 12, 24 (synthetic spectra), fits the DENSE scaling exponent (expected ~N⁴ per ordering), and projects the cost of one solver iteration at each ladder rung (three, e3, e6, e17, and the ×I8 products up to 136). If — and ONLY if — Task 10 established the `'factored'` identity, its scaling is measured and projected alongside. STOP RULE (two independent triggers): **(i) FACTORIZATION NOT ESTABLISHED** — record it, keep `'dense'` as the sole engine, and proceed (large rungs may then be budget-refused); **(ii) BUDGET** — if the projected time for ONE production A3 solve (≈ 30 outer iterations) at a given rung exceeds `budget_hours = 12` under the available engine, that rung is refused (Task 13) and the optimization backlog (proving/using the factored identity, transition/Liouville-space contraction, time-simplex matrix quadrature, symmetry blocking by electronuclear quantum numbers, ωℓ-grid compression, tail subtraction, cached transition sums) is documented as the next step. Reaching the budget limit is NOT a plan failure — it bounds which rungs Task 13 reports.

- [ ] **Step 1: Write the failing tests** (core): kernel rows vs oracle (incl. repeated-node and large-βΔ KMS rows, AbsTol 1e-11 relative, whole-value per constraint 4); `'dense'` `V` vs oracle on all fixture systems (random, degenerate, populated-spectator, physical chain; 1e-9) — the PRIMARY correctness gate; two-level Jensen bridge (as v1, code unchanged — `V == G0.*sg.Sigma`, RelTol 1e-6, using `invz_twolevel`/`invz_sigma` from `invz_common`); degenerate-doublet no-NaN (as v1, without the bogus centering line — Jc = diag(m, −m, 0) is already centered at equal populations); negative-n reconstruction vs stored negative-n oracle rows. **CONDITIONAL (only if the fixture's `factored_ok` is true):** `'factored'` `V` vs oracle + a dense-vs-factored agreement test (1e-9) on random/degenerate/populated systems. If `factored_ok` is false, instead add a test asserting `opts.impl='factored'` errors `invzt:factoredUnproven`.
- [ ] **Step 2: Run to verify failure. Step 3: Implement** kernels → the DENSE engine (production). Wire the `'factored'` path ONLY if Task 10's `factored_ok` is true; otherwise leave it erroring `invzt:factoredUnproven` and rely on `'dense'`. **Step 4: CORE suite green.**
- [ ] **Step 5: Run `perf_vertex_scaling.m`**; controller logs the DENSE scaling table + per-rung projections + gate verdict (and the factored projection iff established) to ODD-LOG §A3. If either STOP trigger fires, record it and continue with Task 12 (three-state, whose N ≤ 3 is cheap under dense regardless).
- [ ] **Step 6: Commit** — `feat(invzt): KMS-scaled kernels + oracle-verified factored vertex engine with performance gate`.

---

### Task 12: `invzt_threestate.m` + `invzt_sigma_tensor.m` + solver mode 'a3' — three-state A3 validation

**Files:**
- Create: `invz_tensor/invzt_threestate.m`, `invz_tensor/invzt_sigma_tensor.m`, `invz_tensor/tests/test_invzt_a3_threestate.m`
- Modify: `invz_tensor/invzt_solve_point.m` (add `mode` 'a3' and `opts.nlevels`)

**THREE-STATE MODEL CONTRACT (v3 — review Blocking 2, option (a): independent low-doublet tunnelling `Delta1`, `rho` reserved for spectator coupling; replaces the v2 overconstrained form):**
`ts = invzt_threestate(ion, T, B, opts)` builds, in the fixed basis {|1⟩,|2⟩,|3⟩}:
```
H3(B, hx) = diag(0, 0, Delta2) + (Delta1/2)*Sx12 - (gmuB*Bx + hx)*Ja0 - hz*Jc0
Sx12 = |1><2| + |2><1|                                  % DIRECT low-doublet tunnelling (v3): the independent splitting knob
Jc0  = diag(m0, -m0, 0)                                 % Ising doublet moment (cc channel)
Ja0(1,3)=Ja0(3,1)=Ja0(2,3)=Ja0(3,2)=rho/sqrt(2)        % SPECTATOR (transverse) coupling ONLY — does NOT split the doublet at first order
Jb0  = (rho/sqrt(2)) * (i*|1><3| - i*|3><1| - i*|2><3| + i*|3><2|)   % 90-deg partner, Hermitian
```
**Why the direct term (v3):** the v2 model made the doublet splitting a *second-order* effect through |3⟩ (∝ (gμB·Bx·ρ)²/Δ2), so a single `rho` had to set BOTH the splitting AND the transverse susceptibility while `m0` set the weight — three targets (Δ, M2, χ⊥) on two knobs (overconstrained), and `rho→0` erased the splitting mechanism entirely (a contradictory scalar-limit test). With the direct `Delta1` term the doublet block at `hz=0` is `[[0, Delta1/2],[Delta1/2, 0]]`: splitting `Delta1`, transition moment `|<+|Jc0|−>| = m0` (M2 = m0²), INDEPENDENT of `rho`. So `rho` is free to carry ONLY the spectator/transverse channel.
Parameters and the WELL-POSED match: `Delta2` = 0.9306 meV (`opts.far_excited` true → 40 meV); the constructor NUMERICALLY solves the three knobs `(Delta1, m0, rho)` against the three targets — `eig(H3)` doublet splitting = `tl.Delta`, `|Mz(1,2)|²` = `tl.M2` (both from `invz_twolevel(ion, T, Bx)`), and the computed transverse susceptibility χ⊥ = `opts.chiperp_target` (default 11.05 meV⁻¹). This is a 3×3, diagonally-dominant Newton solve (each knob mainly drives one target: `Delta1`→Δ, `m0`→M2, `rho`→χ⊥); seed `Delta1 = tl.Delta`, `m0 = sqrt(tl.M2)`, `rho = sqrt(chiperp_target*Delta2/2)` (the exact ρ=0 doublet + leading spectator estimate). The v2 closed-form splitting expression is SUPERSEDED by this numerical solve (no fragile `2*sqrt(...)` formula). Errors `invzt:threeStateMatch` if the residual exceeds `1e-10` after Newton (a bad toy must not silently proceed). `opts.chiperp_scale` rescales the matched `rho → chiperp_scale·rho` and the constructor RE-REFINES `(Delta1, m0)` at that rescaled `rho` so the doublet ALWAYS reproduces `(tl.Delta, tl.M2)` and only χ⊥ varies (∝ chiperp_scale²) — NOT a post-hoc multiply that would leave `(Delta1, m0)` carrying the old ρ²/Δ2 correction (documented: this changes both odd-on and odd-off physics — baselines recomputed per scale). At `chiperp_scale = 0` the re-refinement gives EXACTLY the two-level `(Delta1, m0) = (tl.Delta, √tl.M2)` with |3⟩ decoupled — the basis of the exact ρ→0 scalar gate below (exact at ANY Δ2, since disconnection not depopulation makes it exact). `hx` self-consistent transverse MF (legacy_x analogue) optional via `opts.transverse_mf`. Hyperfine: EXCLUDED at this rung (stated).
Returns a COMPLETE toy `si` struct — `E` (3×1, ground-shifted, from `eig(H3)`), `V` (eigenvectors), `P` (Boltzmann), `Mx/My/Mz` = eigenbasis transforms of `Ja0/Jb0/Jc0`, `Jexp`, `JzJz_fluct`, `hx/hz`, `mf_converged`, plus `match` (achieved Δ, M2, χ⊥ and residuals — REPORTED, since χ⊥ is now matched not imposed-by-luck) — every field the downstream chain (`invz_chi0z`, `invzt_chi0_split`, sum-rule reporting) reads. `nlevels='three'` in `invzt_solve_point` routes single-ion construction through `invzt_threestate` for ALL modes (a1/a2/a3), with `Esplit` selecting the doublet.

**A3 assembly (unchanged from v1 except the two v2 fixes):**
- `invzt_sigma_tensor` assembles `Vmat` `[3,3,nwn]` via `invzt_vertex4`. **`opts.dress` (v4, Task-12):** `'full'` (default — internal channels {a,b,c} against `pt.Kmat`, the genuine tensor A3) or `'dominant'` (E1's dominant-transition rule: the dominant cc vertex only, transverse a,b spectator BARE, spectator-normalized populations, E1 dominant/rest `ctil0` for crit). `'dominant'` is the matched-truncation used by the §11.8 ODD-emergence gate (it reduces to E1 up to the O(1/z²) resummation ambiguity); `'full'` adds the beyond-E1 transverse-spectator dressing.
- Resummation: `G0til = pagemtimes(G0, pagemldivide(G0 + Vmat, G0))` (the LOCKED symmetric bracket; active-subspace policy shared with A2); `chi_til = -G0til`; lattice + medium via `invzt_emt_matrix`; damped outer mix on `Vmat`.
- `pt.Sigma_cc_equiv = +Vmat(3,3,:)./G0_cc` where |G0_cc| > rank_tol (**v2 SIGN FIX**: 𝒱 = G0·Σ ⇒ Σ = 𝒱/G0; DIAGNOSTIC ONLY).
- `pt.resum_spread_crit` = crit(dyson) − crit(additive); monitors `eps_el`, `eps_cross`, sum-rule residual.
- crit via the Task-6 Hermitian similarity form.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_threestate_matches_twolevel_sector(testCase)
% Constructor contract: doublet sector reproduces invz_twolevel's (Delta, M2).
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
tl = invz_twolevel(ion, T, B(1), struct());
ts = invzt_threestate(ion, T, B, struct());
verifyEqual(testCase, ts.E(2) - ts.E(1), tl.Delta, 'RelTol', 1e-9);
verifyEqual(testCase, abs(ts.Mz(1,2))^2, tl.M2, 'RelTol', 1e-9);
verifyTrue(testCase, all(isfield(ts, {'P','Jexp','JzJz_fluct','Mx','My'})));
end

function test_a3_scalar_compat_exact_rho0(testCase)
% EXACT scalar-compatibility gate (category 1) — CONSISTENT under v3 option (a).
% chiperp_scale = 0 sets rho = 0, so Ja0 = Jb0 = 0 and |3> has ZERO matrix
% elements to the doublet in EVERY response (Mx = My = 0; Mz(*,3) = 0): the
% spectator is disconnected regardless of its thermal population, and the doublet
% KEEPS its splitting Delta1 = tl.Delta and moment m0 = sqrt(tl.M2) from the DIRECT
% tunnelling term. So the toy is EXACTLY the two-level system and A3 must reproduce
% the scalar two-level chain (invz_emt_scalar + invz_sigma on the same branch
% spectrum) to 1e-8. (This is the fix for the v2 contradiction, where rho -> 0
% also erased the splitting.) far_excited is kept as belt-and-suspenders but is no
% longer load-bearing: disconnection, not depopulation, makes the limit exact. The
% populated-spectator normalization (constraint 3) is the ORACLE's job (Task 10
% check 6), not an exact solver identity.
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('mode', 'a3', 'nlevels', 'three', ...
    'chiperp_scale', 0, 'far_excited', true, 'odd', false));
verifyTrue(testCase, pt.converged);
% scalar reference on the SAME cc branch spectrum and the SAME (Delta, M2):
tl = invz_twolevel(ion, T, B(1), struct());
[wn, wts, beta] = invz_matsubara(T, 40);
g = real(invz_g(tl, 1i*wn));
G0 = -tl.M2 * g;
Jnu = local_cc_branches(lat);                 % helper: sort(eig) of cc blocks, flattened
Sigma = zeros(size(wn));
for it = 1:60
    med = invz_emt_scalar(G0, Sigma, Jnu, struct());
    lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
    sg = invz_sigma(tl, lam, med.K, g, beta);
    dS = max(abs(sg.Sigma - Sigma));  Sigma = Sigma + 0.7*(sg.Sigma - Sigma);
    if dS < 1e-10, break, end
end
verifyEqual(testCase, real(pt.Sigma_cc_equiv(1)), Sigma(1), 'AbsTol', 1e-8);
end

function test_a3_emergence_matched_truncation(testCase)
% THE framework SS11.8 emergence gate (v4, Task-12 physics — SUPERSEDES the v3
% slope>=2.3 protocol, which was internally INCONSISTENT with LOCKED constraint 8:
% A3's whole-cc Dyson + matrix-EMT resummation and A1's dom/rest + K-bookkeeping
% differ at O(1/z^2), so their crit-shift lambda^2 COEFFICIENTS differ -> the
% crit-shift slope is ~2 inherently and >=2.3 is unreachable for ANY A3-vs-A1
% crit-shift). The emergence is validated in TWO sectors instead:
%   SCALAR (exact): A3 vertex -> invz_sigma at rho->0 (test_a3_scalar_compat_exact_rho0, 1e-11).
%   ODD (matched truncation, here): A3 with dress='dominant' (E1's dominant-transition
%   rule: cc vertex only, transverse spectator BARE) reduces to A1/E1 up to the
%   O(1/z^2) resummation ambiguity, while full-A3 (dress='full') ADDS the beyond-E1
%   transverse-spectator dressing. Stable PM anchor T=2.0 K, Bx=0.5 T (A1 crit>0,
%   single-root; Sigma_seed continuity across lambda).
ion = invz_ion();  T = 2.0;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
lams = [1.0 0.5 0.25];  rd = zeros(1,3);  rf = zeros(1,3);
o = @(m,dr,od,sd) struct('mode',m,'nlevels','three','dress',dr,'odd',od,'Sigma_seed',sd);
b1 = invzt_solve_point(ion, T, B, lat, o('a1','full',false,[]));
verifyGreaterThan(testCase, b1.crit, 0);          % stable single-root PM baseline (documented anchor)
bf = invzt_solve_point(ion, T, B, lat, o('a3','full',    false,[]));
bd = invzt_solve_point(ion, T, B, lat, o('a3','dominant',false,[]));
for i = 1:3
    latL = scale_odd_blocks(lat, lams(i));        % helper: scale ca/cb (+ ac/bc) ODD blocks by lambda
    a1 = invzt_solve_point(ion, T, B, latL, o('a1','full',    true, b1.Sigma));
    af = invzt_solve_point(ion, T, B, latL, o('a3','full',    true, bf.Sigma));
    ad = invzt_solve_point(ion, T, B, latL, o('a3','dominant',true, bd.Sigma));
    d1 = a1.crit - b1.crit;
    rf(i) = (af.crit - bf.crit)/d1;               % full-A3 / A1  (beyond-E1 dressing)
    rd(i) = (ad.crit - bd.crit)/d1;               % dominant-dress-A3 / A1  (matched to E1)
end
pf = polyfit(log(lams), log(abs(rf - 1) + eps), 1);
collapse_frac = 0.4;  band = 0.05;                % O(1/z^2) resummation window; measured |rd-1|~0.016, ratio~0.14
fprintf(['SS11.8 emergence (matched trunc): rd(1)=%.4f rf(1)=%.4f | collapse |rd-1|/|rf-1|=%.3f (<=%.2f) | ' ...
    'band |rd-1|=%.4f (<=%.2f) | full-A3 excess slope %.2f (O(1/z^2)) resum_spread %.2e  (rf/slope/resum are REPORTS)\n'], ...
    rd(1), rf(1), abs(rd(1)-1)/abs(rf(1)-1), collapse_frac, abs(rd(1)-1), band, pf(1), bf.resum_spread_crit);
verifyLessThanOrEqual(testCase, abs(rd(1)-1), collapse_frac*abs(rf(1)-1));  % dominant COLLAPSES the beyond-E1 excess (transverse dressing confirmed)
verifyLessThanOrEqual(testCase, abs(rd(1)-1), band);                       % dominant residual within the O(1/z^2) resummation band
end

function test_uniform_shift_response_reported(testCase)
% As v1 (report-only probe; the strict theorem lives at the closed-form Tc root).
ion = invz_ion();  T = 1.6;  B = [0.1 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
d = 5e-4;
latS = lat;
latS.Jt(3:3:12, 3:3:12, :) = latS.Jt(3:3:12, 3:3:12, :) + d*repmat(eye(4), 1, 1, size(lat.Jt, 3));
latS.JtGamma(3:3:12, 3:3:12) = latS.JtGamma(3:3:12, 3:3:12) + d*eye(4);
p0 = invzt_solve_point(ion, T, B, lat,  struct('mode', 'a3', 'nlevels', 'three', 'odd', false));
p1 = invzt_solve_point(ion, T, B, latS, struct('mode', 'a3', 'nlevels', 'three', 'odd', false));
r = (p1.crit - p0.crit) / d;
fprintf('uniform-shift crit response (A3, three-state): dcrit/d = %.3f\n', r);
verifyTrue(testCase, isfinite(r));  verifyLessThan(testCase, abs(r), 5);
end

function test_a3_monitors_reported(testCase)
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
pt = invzt_solve_point(ion, T, B, lat, struct('mode', 'a3', 'nlevels', 'three', 'odd', true));
fprintf('A3 monitors: sumrule_rel %.3f, eps_el %.3f, eps_cross %.4f, resum spread %.2e\n', ...
    pt.sumrule_rel, pt.eps_el, pt.eps_cross, pt.resum_spread_crit);
verifyTrue(testCase, pt.converged && isfinite(pt.eps_el));
verifyLessThan(testCase, pt.sumrule_rel, 0.25);
end
```

(Local helpers `local_cc_branches` and `scale_odd_blocks` ~6 lines each, in the test file.)

- [ ] **Step 2: Run to verify failure. Step 3: Implement** `invzt_threestate` (constructor solves the 3-param (Delta1, m0, ρ) → (Δ, M2, χ⊥) match, erroring `invzt:threeStateMatch` on non-convergence; header carries the contract box) → `invzt_sigma_tensor` → mode 'a3'. **Step 4: CORE suite green** (minutes; if the λ-emergence slope gate fails, STOP — derivation-level finding, report with the convergence evidence). Controller logs everything to ODD-LOG §A3.
- [ ] **Step 5: Commit** — `feat(invzt): A3 on the explicit three-state model — exact rho->0 scalar gate, lambda-protocol Gaussian cross-validation`.

---

### Task 13: A4 — basis-defined state-space ladder (data-only driver + separate report)

**Files:**
- Create: `invz_tensor/invzt_run_ladder.m` (returns data ONLY), `invz_tensor/invzt_report_ladder.m` (serializes a completed run for ODD-LOG pasting), `invz_tensor/invzt_rung_basis.m` (rung → basis projector + basis energies; used by the solver and directly testable), `invz_tensor/tests/test_invzt_a4_ladder.m`
- Modify: `invz_tensor/invzt_solve_point.m` (`nlevels` rungs)

**Interfaces:**
- **Rungs (v3 — basis-content-defined, LOCKED):** `'three'` (toy), `'e3'`, `'e6'`, `'e17'` (electronic CF subspaces: build the single-ion Hamiltonian IN the reduced basis — lowest 3/6/17 CF states of the zero-field CF Hamiltonian as basis, complete degeneracy multiplets ALWAYS included together — then diagonalize with field/MF terms), `'e3xI8'`, `'e6xI8'`, `'e17xI8'` (= full 136; tensor products with the complete nuclear space). NEVER a lowest-N-eigenstate cut of the 136-dim spectrum. **Multiplet-completion rule (v3, review Other 6):** if the nominal `eN` count would SPLIT a degeneracy multiplet, include the whole multiplet and record the ACTUAL basis dimension in `out.dim_actual` (the `eN` label is nominal, not an exact size). Each rung records: the basis projector, omitted THERMAL weight, and a VIRTUAL-completeness DIAGNOSTIC (v3): the static χ0 tensor of the rung vs the full-136 χ0 at the same (T,B) — the deficit DIAGNOSES which virtual intermediate states the rung is missing. It is a diagnostic ONLY, NOT a mathematical bound on the missing four-point vertex (no such inequality is derived; small thermal population does not bound virtual contributions).
- `out = invzt_run_ladder(ion, opts)` — scalar struct with array fields (`out.rungs` cellstr, `out.crit_shift_odd`, `out.converged`, `out.eps_el`, `out.eps_cross`, `out.sumrule_rel`, `out.chi0_virtual_deficit`, `out.npaths`, `out.t_sec` — all [nr,1]); validation point (1.6 K, [0.1 0 0], 6³/dpRng 10; production settings behind `opts.production`). Tc capability: `opts.tc` true adds the SMALL-Bx-PROXY Tc (0.05 T) per rung via `invzt_tc_pm_extrap` (Task 7) — per the v3 zero-field policy, tensor A3 true-B=0 Tc is DEFERRED, so EVERY rung's `opts.tc` output is the proxy, never a true-zero-field number (only the projected closed form is at B=0). NO file writes.
- Budget: consumes Task 11's measured scaling; refuses rungs whose projected cost exceeds `opts.budget_hours` (12) with a recorded verdict (`out.skipped_rungs`), instead of running them.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_ladder_fast_two_rungs(testCase)
ion = invz_ion();
out = invzt_run_ladder(ion, struct('rungs', {{'three', 'e3'}}));
verifyEqual(testCase, numel(out.rungs), 2);
verifyTrue(testCase, all(isfinite(out.crit_shift_odd)));
verifyTrue(testCase, all(out.converged));
verifyTrue(testCase, all(out.chi0_virtual_deficit >= 0));
fprintf('ladder fast: crit shifts %s\n', mat2str(out.crit_shift_odd(:).', 4));
end

function test_rung_basis_multiplet_complete(testCase)
% v3 (review Other 6): assert COMPLETE-MULTIPLET inclusion and record the ACTUAL
% dimension, rather than hard-coding the nominal label size. At Bx = 0 the CF
% ground doublet is degenerate; the e3 basis must hold BOTH doublet states (no
% split multiplet at the cut) and expose dim_actual + multiplet_complete. The
% nuclear product multiplies by the full I=7/2 space (8).
ion = invz_ion();
rb = invzt_rung_basis(ion, 'e3');
verifyTrue(testCase, isfield(rb, 'dim_actual') && rb.dim_actual == size(rb.projector, 2));
verifyGreaterThanOrEqual(testCase, rb.dim_actual, 3);                 % nominal 3; larger only if a multiplet completes
verifyLessThan(testCase, abs(rb.E_basis(2) - rb.E_basis(1)), 1e-6);   % ground doublet intact (not split by the cut)
verifyTrue(testCase, rb.multiplet_complete);                         % no partial multiplet at the top edge
rb8 = invzt_rung_basis(ion, 'e3xI8');
verifyEqual(testCase, rb8.dim_actual, 8*rb.dim_actual);              % complete nuclear product
end

function test_ladder_production_slow(testCase)
% HEADLINE (report): rung table at production settings + small-B-proxy Tc for
% the largest affordable rung. Cross-validation comparators (REPORT, never
% tune): projected Tier-1+2 = 1.509 K, DeltaTc = 0.2341 K; A3-beyond-Gaussian
% share vs the projected Tier-2 share (~2.8%).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
out = invzt_run_ladder(ion, struct('rungs', {{'three','e3','e6','e17','e3xI8'}}, ...
    'production', true, 'tc', true, 'budget_hours', 12));
verifyTrue(testCase, isstruct(out) && ~isempty(out.rungs));
end
```

- [ ] **Step 2: Run fast tests to verify failure. Step 3: Implement** rung-basis construction + driver + report serializer. **Step 4: CORE suite green; run the slow ladder** (`INVZ_SLOW=1`; honor the budget refusals). **Step 5: Controller runs `invzt_report_ladder(out)`** and pastes §A4 into ODD-LOG: per-rung table with rung-convergence error bars, virtual-deficit column, feasibility verdicts, the largest-rung Tc vs 1.509/1.743/exp 1.53 K, the beyond-Gaussian share vs 2.8%, all monitors.
- [ ] **Step 6: Commit** — `feat(invzt): A4 basis-defined ladder with virtual-completeness diagnostics, budget refusals, data-only driver`.

---

### Task 14: Reproducibility pinning + robustness spot-checks

**Files:**
- Modify: `invz_tensor/tests/invzt_anchors.m` (append measured headline values with ODD-LOG provenance), `invz_tensor/tests/interop/test_invzt_critical_parity.m` (add `test_headline_reproducibility_slow`)

- [ ] **Step 1:** Pin into `invzt_anchors.m`: the Task-7 A1 proxy-Tc, the Task-12 λ-emergence slope AND λ=1 A3/A1 ratio (report values), the Task-13 largest-rung Tc and crit shifts (full-precision literals).
- [ ] **Step 2:** Add `test_headline_reproducibility_slow` (INVZ_SLOW, interop): rerun the Task-7 A1 proxy-Tc measurement; assert reproduction of the pinned anchor to 1% (rerun-stability, never magnitude).
- [ ] **Step 3: Robustness spot-checks** (one cheap point each, controller-logged): A1 proxy-Tc at 12³ vs 16³ (+ Richardson note); ODD-block dpRng 20/30/40 smallest-shell; vertex `mtol` 1e-12 vs 1e-10 and `dd_tol` ×10 on one A3 point (pruning/cluster insensitivity).
- [ ] **Step 4: CORE + INTEROP green (incl. the new slow test). Commit** — `test(invzt): pin headline anchors, reproducibility gate, robustness error bars`.

---

### Task 15: Docs, supersession trail, handoff

**Files:**
- Create: `docs/SESSION-<date>-invz-tensor-full.md`
- Modify: `invz_tensor/README.md` (replace scaffold), `docs/ODD-LOG.md` (final summary)

- [ ] **Step 1: Rewrite `invz_tensor/README.md`:** module map; mode ladder (`a1`/`a2`/`a3` × `nlevels` rungs); flag surface (`odd`, `chi_rest`, `Esplit`, `chi0_diag`, `force_sigma0`, `resum`, `mtol`, `dd_tol`, `rank_tol`, `chiperp_scale`, `far_excited`, `impl`); the suite topology (core/interop/projected) and the grid contract; cache namespace; the Tier-2 scope box; the LOCKED conventions (index order, K bookkeeping, symmetric bracket, transpose relation, decomposition language, elastic convention, v3: Hermitian-eigendecomposition crit, three-state independent-tunnelling contract, DENSE-primary vertex with `'factored'` gated on the Task-10 proof); headline results with error bars; open items (ordered-phase tensor Σ via framework §9 + three-level Eq. 36; A3 real-axis continuation; A3 true-zero-field Tc (deferred — proxy only); the factored-vertex identity if `factored_ok` was false; skipped ladder rungs + the optimization backlog if the perf gate fired).
- [ ] **Step 2: Verify the supersession trail** (stamped 2026-07-17 during plan review — verify intact, re-add only if missing): SUPERSEDED banner atop `2026-07-16-invz-tensor-odd.md`; updated pointer + the two locked bookkeeping relations in `jensen_1z_framework.html` §11.8; updated pointer + staged Tier-1/Tier-2 wording in `invz_projected/README.html` §1.9.
- [ ] **Step 3: Write the session/handoff doc:** what exists, every function/flag/cache key, measured tables, deferred items with source pointers, production runs left to the user.
- [ ] **Step 4: Full verification:** CORE + INTEROP + PROJECTED suites green (fast; slow where defined); PROJECTED counts identical to the Task-1/2 frozen baseline; CORE verified once more with `invz_projected` off the path.
- [ ] **Step 5: Commit** — `docs(invzt): README, ODD-LOG summary, session handoff (A0–A4 v2)`.

---

## Execution notes for the controller

- Tasks strictly sequential; no parallel implementers.
- Model guidance: Tasks 1–2 mechanical; 3–5, 8 standard; 6–7, 9 physics-careful; **10–13 are the research core** (most capable tier, same for reviews); 14–15 standard.
- Before writing any test that calls a projected function (INTEROP only), READ its source for the true signature/out-field names (`invz_odd_zero_field`, `invz_chi_tensor_ref`, `invz_emt_scalar`) — the sketches mark every adaptation point.
- The CORE suite must be verified at least twice (Task 3, Task 15) from a session where `invz_projected` was never added to the path.
- MATLAB runtime discipline: warm `jqt` caches once; first 16³/dpRng-30 build takes minutes; libraries serial.
- The stop-rules are load-bearing: Task-1 precondition (dirty tree), §0 physics stop-rule, Task-6 frozen-baseline gate, Task-11 performance gate, Task-12 λ-gate, Task-13 budget refusals. When one fires, STOP and report — do not improvise.
