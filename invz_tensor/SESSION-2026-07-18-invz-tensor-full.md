# Session handoff — full-tensor 1/z (A0→A4) for LiHoF4

**Date:** 2026-07-18 · **Branch:** `invz-1z-lihof4` · **Plan:** `docs/superpowers/plans/2026-07-17-invz-tensor-full.md` (v4, committed `c011a4f`)
**Status:** All 15 plan tasks complete; final whole-branch review **Ready to merge: Yes** (0 Critical / 0 Important / 0 must-fix). Suites green.

This note is the human handoff. The module reference (module map, mode ladder, flag surface, LOCKED conventions, headline results, open items) lives in `invz_tensor/README.html`; the measured physics log is `docs/ODD-LOG.md` §A0–§A4. This note adds the *why*, the execution story, and what's left.

---

## What was built

A staged full-tensor (3 Cartesian × 4 sublattice, multi-level) generalization of Jensen's 1/z effective-medium expansion, in `invz_tensor/` (`invzt_` prefix), plus a shared `invz_common/` single-ion engine (pure `git mv` out of `invz_projected/`, zero logic change) and a first-principles Python vertex oracle:

- **A0** — cached `[12,12,nq]` coupling tensor + page-wise 12×12 tensor RPA (`invzt_qgrid`, `invzt_jq_tensor`, `invzt_chi_rpa`, `invzt_gcc_lattice`).
- **A1** — dominant-sector projected-1/z bridge point solver + critical finders + real-axis continuation (`invzt_chi0_split`, `invzt_solve_point` mode `'a1'`, `invzt_critical`, `invzt_tc_pm_extrap`, `invzt_chi_realaxis`).
- **A2** — DIRECT matrix effective-medium closure (`invzt_emt_matrix`, mode `'a2'`).
- **A3** — genuine tensor 1/z self-energy from the exact component-labelled four-point vertex (`invzt_kernels`, `invzt_vertex4`, `invzt_gamma4`, `invzt_vertex3`, `invzt_threestate`, `invzt_sigma_tensor`, mode `'a3'`); oracle `verify_tensor_vertex.py` → `tests/fixtures/vertex_oracle.json`.
- **A4** — basis-defined state-space ladder toward the full 136-dim electronuclear space (`invzt_rung_basis`, `invzt_run_ladder`, `invzt_report_ladder`).

**Suites (final):** CORE `runtests('invz_tensor/tests')` = 47 / 0 / 1 (the 1 Incomplete is `INVZ_SLOW`-gated); INTEROP = 8 / 0 / 2; PROJECTED `runtests('invz_projected/tests')` = **143 / 0 / 19, bit-identical to the pre-work baseline** (the projected results were never disturbed — the non-negotiable held across all 15 tasks).

---

## Key physics results

- **A1 proxy-Tc** (0.05 T proxy, 16³/dpRng-30) = **1.5599 K** vs the grid-matched projected closed form 1.5442 K → +0.016 K A1 enhancement (retardation ≈ null + transverse RPA + χ_rest). Zero-field parity `d`/`Sc` to ~1e-15 vs the projected chain.

- **χ0(iωn) is NON-Hermitian off the static slot.** The gyrotropic ∝B response is a physical imaginary-antisymmetric term (measured anti-Hermitian part 6–17 %, up to 9.75 % at Bx = 0.5 T). K / χ_bar / Vmat / Kmat obey the LOCKED transpose relation `X(−iωn)=X(iωn).'` and are Hermitian only at ωn=0; the gyrotropic part is preserved, never symmetrized. (This corrected an internal contradiction in the original plan, which had asserted Hermiticity — see "execution story".)

- **The O(N³) vertex factorization is DISPROVEN.** The Python oracle (`factored_ok=false`) showed the naive resolvent chain mismatches the dense path-sum by 1e15–1e17 on all 24/24 orderings (spurious poles at ζn=Ea, dropped 𝓘₃/KMS structure). The A3 engine is therefore **dense-only** (`opts.impl='factored'` errors `invzt:factoredUnproven`). Perf gate: `{three, e3, e6}` affordable under a 12 h budget; `{e17, all ×I8 up to 136}` budget-refused (dense-infeasible — see open items).

- **Framework §11.8 emergence — SHARPENED into a rigorous two-sector statement** (the plan's original "A3-Gaussian reproduces E1 to leading order in λ, slope ≥ 2.3" premise was found inconsistent with the plan's own constraint 8):
  - **Scalar sector, EXACT:** the A3 vertex reduces to the scalar `invz_sigma` at ρ→0 to **3.24e-11**.
  - **ODD sector, matched truncation:** with `dress='dominant'` (Jensen's dominant-only rule; the transverse spectator provably bare via the toy's reflection symmetry), A3 reproduces A1/E1 up to the O(1/z²) resummation-scheme ambiguity — the beyond-E1 excess **collapses** from full-A3 rf=1.1132 to dominant rf=1.0159 (86 % removed), with the 1.6 % residual matching the constraint-8 method error bar (`resum_spread_crit` = −0.039).
  - **`dress='full'` adds a genuine beyond-E1 correction:** A3 *dresses* the transverse spectator that E1's dominant-only rule leaves bare — A3 is more complete than E1, not wrong.

- **A4 basis climb — independent cross-validation of the tensor 1/z against the projected Tier-2.** As the basis grows toward the real CF content, the beyond-E1 transverse-dressing share **converges to the projected Tier-2 value**: full-A3/A1 ratio rf = three 1.1132 → e3 1.1140 → **e6 1.0263** (share 11.3 % → 2.6 %, landing on the projected Tier-2 ≈ 2.8 %); the virtual-completeness deficit (a diagnostic, not a bound) drops monotonically 0.041 → 0.002. The small toys *overestimate* the beyond-Gaussian content; e6's CF basis screens it to the projected value.

---

## Execution story (why the reviews mattered)

Built via subagent-driven-development (fresh implementer + independent reviewer per task). The per-task and final reviews caught — and the fixes/adjudications resolved — several *real* issues that a rubber-stamp pass would have shipped:

- **T5:** the plan's `Esplit=0.1` sat inside the ground hyperfine manifold → a correct split legitimately failed the test. Amended to 0.2 (band (0.13, 0.93) meV).
- **T7:** the Bc-parity test's `chi_rest=false` mismatched the projected full-χ0 solver (0.16 T high-field error). Amended to `true`.
- **T9:** the plan asserted K/χ_bar Hermitian — **verified false in MATLAB** (χ0 non-Hermitian off static, contradicting the plan's own constraint 9). Tests + interface corrected to the transpose relation.
- **T10:** the O(N³) factorization **disproven** by the oracle — vindicating the v3 decision to gate it behind proof rather than assume it (a v1 build on the false identity would have been wrong).
- **T12 (derivation-level, user-adjudicated):** the §11.8 λ-slope gate was internally inconsistent with constraint 8; reframed to the two-sector matched-truncation form above (user chose "verify then matched-truncation").
- **T13:** a hardcoded `eye(3)` in the A3 assembly (latent from T12, hidden because T12 only tested the N=3 toy) crashed every rung with dim ≠ 3 — the review caught it by *running e6*. Fixed to `eye(N)`, proven at e6, and guarded by a fast structural test.

A `/simplify` cleanup pass (4 parallel angle-agents) then removed genuine duplication (shared `invzt_active_projector` + `invzt_str`), page-vectorized `invzt_chi_rpa` (`pagemldivide`, machine-precision-identical), and hoisted loop-invariants — net −17 lines, CORE unchanged.

---

## Open items / future work

**Production runs left to the user** (the plan's design — the module is data-only and budget-honest):
- The full A4 production Tc scan at **e6** (~10 h per the T11 budget) and larger rungs. `invzt_run_ladder` runs them under `INVZ_SLOW`; the fast `{three, e3}` validation runs in-session.
- Rungs `{e17, all ×I8 up to 136}` are **budget-refused** under dense scaling (~196 h to ~5.9e5 h). Reaching 136 needs the optimization backlog below.

**The dense-only bottleneck (documented backlog):** the disproven naive factorization means the exact tensor vertex is O(N⁴)/ordering. A genuine speedup needs a *derived* transition/Liouville-space contraction (or time-simplex matrix quadrature, symmetry blocking by electronuclear quantum numbers, ωℓ-grid compression, cached transition sums) — proven in the oracle first, exactly as the dense path was.

**Deferred (documented in-code):** A3 real-axis continuation (the `invzt_chi_realaxis` A1 continuation does not extend to the full Vmat(iωn)); A3 true-zero-field Tc (the small-Bx proxy is used for every tensor Tc; only the projected closed form is truly at B=0).

**Fast-follow Minors (from the final review — non-blocking, ship-as-is):**
1. ~~`invzt_chi_realaxis` has no `pt.mode` guard~~ — **RESOLVED 2026-07-18** (drivers work, commit `6bfbd32`): runtime guard `invzt:realaxisMode` rejects any `pt.mode ~= 'a1'`, with a CORE regression test.
2. ~~`invzt_chi_realaxis` explicit-qvec `chi_cc_q` uses odd-on coupling regardless of `pt.odd`~~ — **RESOLVED 2026-07-18** (commit `1ae0c7f`): the odd=false Cartesian mask is now a shared helper (`invzt_odd_mask`) applied by both `invzt_solve_point` and the continuation (the mismatch was a 17.2% response error at the review probe point, not harmless at finite q); CORE regression test pins the masked response to 1e-12.
3. `invzt_emt_matrix.info.antisym_K` uses `pagetranspose` while its "anti-Hermitian" label implies `pagectranspose` (report-only, never asserted).
4. `invzt_run_ladder.rung_cost_hours` labels the T11 *scan-scale* projection as "one-solve cost" (the refusal decision is correct; the label is imprecise).

**Deferred simplifications (note-only, judged too risky/large for the cleanup pass):** factor the Anderson-mix scaffold into a shared `invzt_anderson_step` (convergence-critical — the outer solve is genuinely multi-root/seed-sensitive at some configs, so this needs its own careful validation); hoist the `contract_vertex` G4 re-slice (~7.2 GB of redundant strided copying over a 400-iter production solve) and the `gamma4` `pairings_grid` 9× redundancy (both mainly help the user's long production runs); generalize `invz_chi0z`/`invz_single_ion`/`invzt_kernels` to eliminate the `masked_chi0` / `local_rung_si` / gamma4-I2/I3 re-derivations (cross-module — touch shared `invz_common`, needs re-validation against the oracle + projected suites).

---

## Provenance & loose ends

- **Plan of record:** `docs/superpowers/plans/2026-07-17-invz-tensor-full.md` (committed `c011a4f`, with all v3/v4 amendments and the physics adjudications).
- **Superseded:** `docs/superpowers/plans/2026-07-16-invz-tensor-odd.md` carries a SUPERSEDED banner (this plan absorbed its A0/A1 layer and fixed its two physics defects).
- **Physics log:** `docs/ODD-LOG.md` §A0–§A4 (dated per-stage controller entries).
- **Framework docs (loose end — user-owned):** the `jensen_1z_framework.html` §11.8 pointer + two locked bookkeeping relations, and the `invz_projected/README.html` §1.9 caveat, are **intact in the working tree but left uncommitted** — they are interleaved with prior-session ODD-main-body framework edits (mixed work-streams), so committing them is the user's call, not the tensor branch's.
- **Commit trail (tensor feature):** `87e4cdf` (projected WIP precondition) → 15 task commits (T1 `8b0de0d` … T15 `7be2f2b`) + the ODD-LOG §A0–§A4 controller entries → plan `c011a4f` → `/simplify` `1f3b802` → README `80ae006`.
