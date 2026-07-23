# Task 2 report — causal-discriminator matrix, complete evidence + frozen-threshold classification

**Status: DONE — decision gate.** This report is analysis-only: it classifies the already-run
40-cell matrix under the frozen pre-registration and presents the classification. **It does not
select or start Task 3A, 3B, or 4.** No solver, checker, harness, or driver file was read for
anything other than understanding the data shape; none was modified.

---

## 0. Corrections applied (errata — added after independent review)

An independent review (`stage2_blocker_review_codex.md`, repo root) of this report identified six
findings, all applied at their cited locations below and re-verified directly against
`.superpowers/sdd/task2_matrix_results.mat` during this correction pass (MATLAB, read-only, no
cell re-run). **The operative verdict — LATTICE/MESH-UNRESOLVED (§7) — is UNCHANGED.**

The verdict does not depend on the questionable index-aligned, mixed-`dpRng` comparisons flagged by
C3 below. It is independently supported by the ONE pairwise comparison unaffected by any of the six
findings: **`baseline` (`unshifted`, `dpRng=30`) vs. `half_step`/`dpRng=30`**, whose H_MF `h`-grids
agree to machine precision at every field (re-confirmed this pass). §3 (new paragraph) gives the
exact, reproducible breakdown: on that pairwise comparison alone, restricted to nodes
checker-accepted on both sides, the frozen §D numeric test **fails at every common-accepted node,
at all 4 physical fields.**

Corrections applied (numbering matches the review):
- **C1** — §5's existence bar corrected: "MET at 1/4" → "MET at 0/4" under the literal frozen
  reading (the prior version miscounted `baseline` + `unshifted/dp40` — the SAME grid offset at two
  `dpRng` values — as "≥2 mesh offsets").
- **C2** — disclosed a second, previously-unflagged untested axis: no cell in this matrix varies BZ
  grid **size**; every rung holds the grid fixed at 16³, and `G3` downsampling is a subset of that
  same fixed grid, not an independent size refinement.
- **C3** — disclosed that the index-aligned §D ladder mixes a lattice change with a small H_MF
  evaluation-point change on the `dpRng=40` rungs, and added the exact-`h` pairwise breakdown
  (`baseline` vs. `half_step/dp30`) that the verdict actually rests on.
- **C4** — flagged that the field-level "existence" criterion (what it means for a field to "yield
  a state" when every cell's `hstar` is `NaN`) is not itself frozen by the prereg text, and marked
  §5's numbers as conditional on this report's own operationalization of it.
- **C5** — narrowed three overstated conclusions: (a) the §4 seed/continuation finding; (b) the §6
  distribution-shape/geometry finding, substantially rewritten after confirming from source that
  the ordered EMT map is permutation-invariant over the flat coupling multiset (so it cannot, by
  construction, distinguish "distribution shape" from "geometry"); (c) the F2581023 `G3` stride-8
  classification, corrected from "identical to G1" to "identical counts, two-node boundary shift"
  (re-verified against the results file, node ids 26/27).
- **C6** — corrected an incorrect floor/relative-tolerance explanation for the worst observed `K`
  discrepancy in §4 (re-verified against the results file: `9.405×10⁻¹⁰` is ≈14.6× the `K` absolute
  floor `1×10⁻⁸·|J0eff|` ≈ `6.424×10⁻¹¹`, not "orders of magnitude inside" it; it passes the
  combined tolerance through the `RelTol` term).

No code, test, driver, harness, or the frozen pre-registration was modified to produce these
corrections, and no cell of the matrix was re-run.

**Refreshed review (second pass) — two further wording corrections + one new technical finding.**
A second independent-review pass confirmed the six corrections above and the operative verdict, and
added:
- **W1 (§7)** — "residual-clean" vs. "lattice-converged" were conflated. The 18-node F3754215
  interval **is** residual-clean (its states are checker-accepted); it is separately **not**
  lattice-converged, and `nH` persistence is untested. Wording corrected in §7.
- **W2 (§6)** — "bitwise-identical under permutation" is too strong for floating-point `mean`
  reductions (reordering can change roundoff). Corrected to "mathematically permutation-invariant,
  up to floating-point reduction-order effects."
- **T1 (new, verified from source — see §10)** — the production BZ grid is the endpoint-inclusive
  `linspace(-0.5,0.5,N)` convention, which has **duplicate reciprocal-periodic boundary faces**: a
  nominal 16³ grid holds only 15³ = 3375 distinct periodic points (17.6 % redundant samples,
  numerically confirmed). Because the EMT map is a `mean` over the flat modes (the §6
  permutation-invariance), those duplicate faces are **over-weighted**, and a half-step offset moves
  *which* modes are duplicated — so the §3 offset comparison confounds sampling phase with
  duplicate-face reweighting and Γ/cardinality handling. This is a plausible mechanistic driver of
  the §3 sensitivity and **the first construction question to resolve** before extending the same
  legacy grid/cutoff ladder. It does not change the frozen-threshold verdict. My **immediate action
  plan** in response is recorded in §11.

---

## 1. Provenance

**What was run.** The full prereg §G experiment matrix: 40 cells (G1G2: 12, G3: 12, G4: 2, G5: 2,
G6: 12), enumerated by `invz_task2_matrix_enumerate.m` and executed by
`invz_task2_matrix_run.m`/`invz_task2_matrix.m`. Confirmed from `.superpowers/sdd/task2_matrix_run.log`:
`n_total_enumerated=40, n_requested=40, n_run=40, n_failed=0, n_checkpointed=40`. Results
checkpointed to `.superpowers/sdd/task2_matrix_results.mat` (git-ignored, 113 MB; `results`,
1×40 struct). This report loaded that file read-only (MATLAB foreground/blocking, `load(...,
'results')`, 3.7–3.8 s) and performed no solve of any kind.

**Frozen pre-registration.** `docs/invzp_task2_prereg.md`, committed `a70e7e9` ("FROZEN
pre-registration... set by the user 2026-07-23 and committed BEFORE any discriminator run").
Every threshold applied below (§A ±1×10⁻³ classification bands; §C `RelTol=1×10⁻⁶`,
`AbsTol_q=1×10⁻⁸·scale_q`; §D `AbsTol=1×10⁻⁶, RelTol=1×10⁻⁴`; §E ≥3-of-4 existence bar) is that
document's literal number, applied via the frozen tools it specifies
(`invz_task2_classify.m`, `invz_task2_agree.m`, `invz_task2_ladder_ok.m`, committed `8be2845`)
with **no re-tuning, no reinterpretation of a numeric value.**

**Harness/driver provenance.** Task 2a harness + classifier (`8be2845`, review-fixed `17a5181`);
Task 2b-driver enumeration/execution machinery (`2b23f53`, review-fixed `dce796d`, default-path
fix `22f7d1b`). Current branch `invzp-stage2c-diagnostic`, HEAD at analysis time
`22f7d1b77d435f574d09ae96160fb20d34698fdb` (base `main`@`0edc0ab`).

**`Bc_PM` and the four frozen physical fields** (Task 2a report, §E): `Bc_PM = 4.692769` T,
derived from `invz_critical` (PM-mass root `1 − J·chi0_PM(0) → 0`) on the real
`invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30))` lattice at `T=0.1` K. The four
physical fields used throughout this matrix:

| tag | field (T) | derivation |
|---|---|---|
| F1173192 | 1.173192 | 0.25·Bc_PM |
| F2581023 | 2.581023 | 0.55·Bc_PM |
| F3754215 | 3.754215 | 0.80·Bc_PM |
| F2850000 | 2.850000 | existing defect anchor |

**Data-integrity checks performed on all 1360 nodes (40 cells × 34 nodes) before any
classification claim below:**
- `replay_mismatch_any` is **false for every one of the 40 cells** (production-trace vs.
  harness-replay cross-check never disagreed).
- The §A accepted-only rule (`class=='unconverged'` iff `~accepted` or non-finite `D_uni`/`Dq_min`;
  `s == min(D_uni, min(Dq))` exactly whenever finite) holds with **0 violations across all 1360
  nodes** — the recorded `class` field is faithful to the frozen rule everywhere in this matrix, no
  exceptions found.
- **0 of 1360 nodes classify `marginal`** (`|s| ≤ 1×10⁻³`) anywhere in the entire 40-cell matrix —
  every accepted node's `s` sits clearly on one side of the ±1×10⁻³ band or the other; class counts
  matrix-wide: 302 stable, 0 marginal, 279 unstable, 779 unconverged (=1360).
- `node.id == array position j` for every node in every cell (position-based node matching, used
  throughout §D/§C/§E below, is therefore valid).
- The H_MF `h`-grid is **index-aligned across every cell at a given field**: same `n_nodes=34`,
  strictly non-decreasing `h` in visit order, and cross-cell `max|Δh|` at matched index is either
  exactly 0, at machine precision (~1×10⁻¹⁵–10⁻¹⁷, for the `half_step/dp30` rung), or ≤3.34×10⁻⁵
  (for `dp40`/synthetic-cardinality cells, where the coupling-dependent grid-construction constants
  shift very slightly) — never a reordering. Full table in §3 below.
- **⚠️ Confirmed evidence gap (stated up front, expanded in §7):** every one of the 40 cells has
  `n_nodes = 34` — the matrix varies seed, grid offset, and `dpRng`, but **never the H_MF grid
  density (`nH`)**.
- **⚠️ Second confirmed evidence gap, disclosed after independent review (§0/C2):** the matrix also
  never varies BZ grid **size** — every cell holds the grid fixed at `16×16×16`, and `G3`
  downsampling (§6) is a structured subset of that same fixed grid, not an independent size
  refinement. Prereg §D names "grid size, grid offset, `dpRng`" as the ladder's three axes; the
  ladder actually run (§3) covers only two of them. Detail and source-provenance citations in §3.

No production/harness/driver file was modified to produce this report. This report's own working
scripts (`.superpowers/sdd/task2b_part1.m`, `task2b_part2.m`, `task2b_debug_qbranch.m`,
`task2b_probe.m`) are read-only against the results file and are git-ignored scratch (under
`.superpowers/`), not part of this commit.

---

## 2. The complete 40-cell matrix (every cell)

`s` range is over **checker-accepted nodes only** (§A); `n/a (0 accepted)` means no node in that
cell was accepted. All 40 cells: `hmf_status` is never `'ok'` — **no cell's H_MF profile ever
reaches a converged root** (`hstar` is `NaN` everywhere); every swept cell's status is
`'node_failed'` (≥1 evaluated node failed to converge, per `invz_hmf_ordered.m`'s own "status must
be truthful on node failure" contract), and every isolated cell reports `'isolated'` (no profile,
by construction).

| # | cfg_id | group | variant | field_T (T) | hmf_status | n_nodes | n_accepted | stable | marginal | unstable | unconverged | s range (accepted) | replay_mismatch |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | g1_swept_F1173192 | G1G2 | swept | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3796, 0.9852] | 0 |
| 2 | g1_isolated_cold_F1173192 | G1G2 | isolated_cold | 1.173192 | isolated | 34 | 10 | 10 | 0 | 0 | 24 | [0.3796, 0.9852] | 0 |
| 3 | g1_isolated_seed2_F1173192 | G1G2 | isolated_seed2 | 1.173192 | isolated | 34 | 10 | 10 | 0 | 0 | 24 | [0.3796, 0.9852] | 0 |
| 4 | g1_swept_F2581023 | G1G2 | swept | 2.581023 | node_failed | 34 | 8 | 7 | 0 | 1 | 26 | [-0.4771, 0.8641] | 0 |
| 5 | g1_isolated_cold_F2581023 | G1G2 | isolated_cold | 2.581023 | isolated | 34 | 8 | 7 | 0 | 1 | 26 | [-0.4771, 0.8641] | 0 |
| 6 | g1_isolated_seed2_F2581023 | G1G2 | isolated_seed2 | 2.581023 | isolated | 34 | 8 | 7 | 0 | 1 | 26 | [-0.4771, 0.8641] | 0 |
| 7 | g1_swept_F3754215 | G1G2 | swept | 3.754215 | node_failed | 34 | 29 | 6 | 0 | 23 | 5 | [-0.4513, 0.6216] | 0 |
| 8 | g1_isolated_cold_F3754215 | G1G2 | isolated_cold | 3.754215 | isolated | 34 | 29 | 6 | 0 | 23 | 5 | [-0.4513, 0.6216] | 0 |
| 9 | g1_isolated_seed2_F3754215 | G1G2 | isolated_seed2 | 3.754215 | isolated | 34 | 29 | 6 | 0 | 23 | 5 | [-0.4513, 0.6216] | 0 |
| 10 | g1_swept_F2850000 | G1G2 | swept | 2.850000 | node_failed | 34 | 9 | 6 | 0 | 3 | 25 | [-1.019, 0.8207] | 0 |
| 11 | g1_isolated_cold_F2850000 | G1G2 | isolated_cold | 2.850000 | isolated | 34 | 9 | 6 | 0 | 3 | 25 | [-1.019, 0.8207] | 0 |
| 12 | g1_isolated_seed2_F2850000 | G1G2 | isolated_seed2 | 2.850000 | isolated | 34 | 9 | 6 | 0 | 3 | 25 | [-1.019, 0.8207] | 0 |
| 13 | g3_ds2_swept_F1173192 | G3 | downsample_stride2 | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3796, 0.9852] | 0 |
| 14 | g3_ds4_swept_F1173192 | G3 | downsample_stride4 | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3798, 0.9852] | 0 |
| 15 | g3_ds8_swept_F1173192 | G3 | downsample_stride8 | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3783, 0.9851] | 0 |
| 16 | g3_ds2_swept_F2581023 | G3 | downsample_stride2 | 2.581023 | node_failed | 34 | 8 | 7 | 0 | 1 | 26 | [-0.4771, 0.8641] | 0 |
| 17 | g3_ds4_swept_F2581023 | G3 | downsample_stride4 | 2.581023 | node_failed | 34 | 8 | 7 | 0 | 1 | 26 | [-0.4648, 0.8641] | 0 |
| 18 | g3_ds8_swept_F2581023 | G3 | downsample_stride8 | 2.581023 | node_failed | 34 | 8 | 7 | 0 | 1 | 26 | [-0.1912, 0.8641] | 0 |
| 19 | g3_ds2_swept_F3754215 | G3 | downsample_stride2 | 3.754215 | node_failed | 34 | 29 | 6 | 0 | 23 | 5 | [-0.4513, 0.6216] | 0 |
| 20 | g3_ds4_swept_F3754215 | G3 | downsample_stride4 | 3.754215 | node_failed | 34 | 31 | 6 | 0 | 25 | 3 | [-0.4533, 0.6216] | 0 |
| 21 | g3_ds8_swept_F3754215 | G3 | downsample_stride8 | 3.754215 | node_failed | 34 | 28 | 6 | 0 | 22 | 6 | [-0.4106, 0.6216] | 0 |
| 22 | g3_ds2_swept_F2850000 | G3 | downsample_stride2 | 2.850000 | node_failed | 34 | 9 | 6 | 0 | 3 | 25 | [-1.019, 0.8207] | 0 |
| 23 | g3_ds4_swept_F2850000 | G3 | downsample_stride4 | 2.850000 | node_failed | 34 | 9 | 6 | 0 | 3 | 25 | [-0.9844, 0.8207] | 0 |
| 24 | g3_ds8_swept_F2850000 | G3 | downsample_stride8 | 2.850000 | node_failed | 34 | 11 | 6 | 0 | 5 | 23 | [-1.552, 0.8206] | 0 |
| 25 | g4_cardsynth_swept_F2850000 | G4 | cardinality_synth | 2.850000 | node_failed | 34 | 33 | 8 | 0 | 25 | 1 | [-19.03, 0.751] | 0 |
| 26 | g4_cardsynth_swept_F1173192 | G4 | cardinality_synth | 1.173192 | node_failed | 34 | 18 | 17 | 0 | 1 | 16 | [-0.04593, 0.9606] | 0 |
| 27 | g5_histmatch_swept_F2850000 | G5 | histmatch_synth | 2.850000 | node_failed | 34 | 9 | 6 | 0 | 3 | 25 | [-1.019, 0.8207] | 0 |
| 28 | g5_histmatch_swept_F1173192 | G5 | histmatch_synth | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3796, 0.9852] | 0 |
| 29 | g6_offunshifted_dp40_swept_F1173192 | G6 | unshifted_dp40 | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3791, 0.9852] | 0 |
| 30 | g6_offhalfstep_dp30_swept_F1173192 | G6 | half_step_dp30 | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3771, 0.9851] | 0 |
| 31 | g6_offhalfstep_dp40_swept_F1173192 | G6 | half_step_dp40 | 1.173192 | node_failed | 34 | 10 | 10 | 0 | 0 | 24 | [0.3804, 0.9852] | 0 |
| 32 | g6_offunshifted_dp40_swept_F2581023 | G6 | unshifted_dp40 | 2.581023 | node_failed | 34 | 8 | 7 | 0 | 1 | 26 | [-0.4649, 0.864] | 0 |
| 33 | g6_offhalfstep_dp30_swept_F2581023 | G6 | half_step_dp30 | 2.581023 | node_failed | 34 | 8 | 7 | 0 | 1 | 26 | [-0.2099, 0.8634] | 0 |
| 34 | g6_offhalfstep_dp40_swept_F2581023 | G6 | half_step_dp40 | 2.581023 | node_failed | 34 | 9 | 7 | 0 | 2 | 25 | [-0.6192, 0.8641] | 0 |
| 35 | g6_offunshifted_dp40_swept_F3754215 | G6 | unshifted_dp40 | 3.754215 | node_failed | 34 | 29 | 6 | 0 | 23 | 5 | [-0.4476, 0.6213] | 0 |
| 36 | g6_offhalfstep_dp30_swept_F3754215 | G6 | half_step_dp30 | 3.754215 | node_failed | 34 | 32 | 6 | 0 | 26 | 2 | [-0.5803, 0.6194] | 0 |
| 37 | g6_offhalfstep_dp40_swept_F3754215 | G6 | half_step_dp40 | 3.754215 | node_failed | 34 | 31 | 6 | 0 | 25 | 3 | [-0.5607, 0.6213] | 0 |
| 38 | g6_offunshifted_dp40_swept_F2850000 | G6 | unshifted_dp40 | 2.850000 | node_failed | 34 | 9 | 6 | 0 | 3 | 25 | [-1.02, 0.8206] | 0 |
| 39 | g6_offhalfstep_dp30_swept_F2850000 | G6 | half_step_dp30 | 2.850000 | node_failed | 34 | 8 | 6 | 0 | 2 | 26 | [-0.2647, 0.8198] | 0 |
| 40 | g6_offhalfstep_dp40_swept_F2850000 | G6 | half_step_dp40 | 2.850000 | node_failed | 34 | 8 | 6 | 0 | 2 | 26 | [-0.2561, 0.8206] | 0 |

**Totals across the 4 physical G1 baseline cells (rows 1,4,7,10):** 136 nodes; 56 accepted (29
stable, 0 marginal, 27 unstable); 80 unconverged.

**Closest ±Dq q/branch at the most-marginal (min-`s`) checker-accepted node, per swept cell**
(G1/G2 use the harness's native `qbranch_pos/neg`; G3/G6 recover it via
`invz_ordered_trace_resolve(struct('is_synthetic',false,'Jnu_unflat',lp.Jnu_unflat,'nq',lp.nq),
idx_flat)` then `lp.qc(q_idx,:)` against `results(k).lattice_provenance`, per the brief's recipe;
G4/G5 have no physical `q` — fully synthetic, no BZ grid — and are marked accordingly). Format:
`(q_idx, branch_idx, J(q), [qx qy qz])`.

| cfg_id | node (min-s accepted) | h | s | D_uni | Dq_min | closest +Dq | closest −Dq |
|---|---|---|---|---|---|---|---|
| g1_swept_F1173192 | id=25 | 0.0061959 | 0.3796 | 0.3796 | 0.4274 | (1946,4,0.005985,[0.100 0.100 −0.033]) | n/a (no −Dq mode) |
| g1_swept_F2581023 | id=26 | 0.0065915 | −0.4771 | −0.4771 | −0.3491 | (1913,4,0.004764,[−0.033 0.033 −0.033]) | (1738,4,0.004848,[0.300 0.100 −0.100]) |
| g1_swept_F3754215 | id=25 | 0.0038270 | −0.4513 | −0.4513 | −0.2875 | (1702,4,0.005196,[0.167 −0.167 −0.100]) | (1877,4,0.005364,[−0.167 −0.233 −0.033]) |
| g1_swept_F2850000 | id=24 | 0.0040573 | −1.019 | −1.019 | −0.8411 | (1205,4,0.00387,[0.233 −0.233 −0.233]) | (1367,4,0.00394,[−0.167 −0.100 −0.167]) |
| g3_ds2_swept_F1173192 | id=25 | 0.0061959 | 0.3796 | 0.3796 | 0.4274 | (1076,4,0.005985,[−0.100 −0.100 0.033]) | n/a |
| g3_ds4_swept_F3754215 | id=25 | 0.0038270 | −0.4533 | −0.4533 | −0.2887 | (599,4,0.005196,[−0.167 0.167 0.100]) | (470,4,0.005364,[−0.167 −0.233 −0.033]) |
| g6_offhalfstep_dp30_swept_F3754215 | id=1 | 0 | −0.5803 | −0.5683 | −0.5803 | (2184,4,0.004837,[0.067 0.067 0.067]) | (1901,4,0.004889,[−0.067 0.333 0.000]) |
| g4_cardsynth_swept_F2850000 | id=21 | 0.0021212 | −19.03 | 5.266 | −19.03 | n/a (synthetic) | n/a (synthetic) |

(The full 40-row version of this table, including every G3/G6 rung, is preserved verbatim in the
working notes, `.superpowers/sdd/task2b_out/frag_qbranch.md`; the rows above are the ones cited by
number elsewhere in this report. Every G1/G3/G6 row resolved cleanly — no unresolved
`invz_ordered_trace_resolve` failures anywhere in the matrix.)

---

## 3. §D — lattice/mesh ladder (adjudicated FIRST, per the frozen table)

**Method, stated explicitly (an implementation choice within the frozen rule, not a
reinterpretation of it — flagged as such by Task 2a's own header: "the brief specifies the rule,
not the exact call shape... 2b's matrix driver will need to build results into this shape").**
`invz_task2_ladder_ok`'s own contract is generic: its `N`-axis is "one comparison point per
slot, in the same order across every rung." This report operationalizes that axis, **per physical
field**, as the field's own **34 index-aligned H_MF nodes** (validated below), rather than as the
4 physical fields themselves — because each field's own 34-node profile is exactly what "does
refining the mesh change the verdict at this h" asks. The four ladder rungs compared at each
field: **baseline** (`unshifted`, `dpRng=30` — the G1 swept cell), `unshifted/dpRng=40` (G6),
`half_step/dpRng=30` (G6), `half_step/dpRng=40` (G6).

**Node-index alignment, confirmed before use:** every cell at a given field has exactly 34 nodes,
`node.id == j` for all j, `h` strictly non-decreasing in `id` order, and cross-rung `h` at matched
index differs by at most ~2.3×10⁻⁵ (an order of magnitude smaller than the local node-to-node `h`
spacing everywhere checked) — confirmed no reordering, so position-`j` matching is valid. (E.g. at
F3754215: `g6_offunshifted_dp40` max|Δh|=2.3e-05, `g6_offhalfstep_dp30` max|Δh|=2.38e-15,
`g6_offhalfstep_dp40` max|Δh|=2.3e-05, vs. baseline.)

**Two axes this ladder does NOT test, disclosed after independent review (§0/C2, C3; see also
§1/§8):**

1. **BZ grid size.** Prereg §D names the ladder's axes as "grid size, grid offset, `dpRng`." Every
   rung compared above holds the BZ grid fixed at the same `16×16×16` mesh — confirmed directly
   from each cell's own stored construction provenance, e.g. `g6_offhalfstep_dp30`'s: `"half_step
   q-grid ([16 16 16], dpRng=30); 1 Gamma-equivalent point(s) excluded via
   invz_is_gamma_equiv."`. `G3`'s deterministic downsampling (§6) is built, per its own provenance
   string, as `"stride-N modular decimation ... of the [16 16 16]/dpRng=30 baseline q-grid"` — a
   **subset** of that same fixed grid, not an independently constructed coarser/finer BZ mesh — and
   `G3` is a density/distribution discriminator (§6), never passed through `invz_task2_ladder_ok`
   or any §D adjudication above. **No structured BZ grid-size refinement (e.g. `12³`/`16³`/`20³`)
   was run anywhere in this 40-cell matrix.** This is a second explicitly-untested axis alongside
   the H_MF-density (`nH`) gap already flagged in §1/§8 — the lattice ladder as actually run covers
   2 of the 3 axes prereg §D names, not all 3.
2. **A single explicit, shared `h_list`.** The comparison above matches rungs by **node index**
   (position `j`), not by a common, explicitly-shared H_MF evaluation-point list — a defensible
   operationalization of `invz_task2_ladder_ok`'s generic "one comparison point per slot" contract
   (Method, above), but with a real consequence: same-`dpRng` rungs share `h` at matched index to
   machine precision (confirmed above), while the `dpRng=40` rungs' own grid-construction constants
   shift `h` at matched index by up to ~2.3×10⁻⁵ relative to baseline. **Every `baseline`-vs-
   `dpRng=40` pairwise difference reported in this section therefore mixes a lattice-construction
   change with a (small) change in the physical H_MF point being evaluated — it is not a
   lattice-only comparison.**

**The §D verdict below does not depend on that mixing.** Restricting to the one rung-pair that
shares `h` EXACTLY — `baseline` (`unshifted/dp30`) vs. `half_step/dp30` — and to nodes
checker-accepted on **both** sides of that pair only (a strictly pairwise, 2-rung restriction,
distinct from the "ACCEPTED-IN-ALL-4-RUNGS" 4-rung set used in the table below), the frozen §D
numeric test still fails at **every** such node, at all 4 fields — re-derived this pass by calling
`invz_task2_ladder_ok` itself (unmodified) on this 2-rung restriction, one call per field, each
common-accepted node treated as one field-slot (the same axis-operationalization as the rest of
this section):

| field | max\|Δh\| (baseline vs. half_step/dp30) | common-accepted nodes | pass §D numeric test | worst margin (s/D_uni/Dq_min) | class disagreements |
|---|---|---|---|---|---|
| 1.173192 | 6.9×10⁻¹⁷ | 10 | **0** | 0.0503 | 0 |
| 2.581023 | 1.29×10⁻¹⁵ | 7 | **0** | 0.0851 | 0 |
| 3.754215 | 2.38×10⁻¹⁵ | 29 | **0** | 0.2747 | 0 |
| 2.850000 | 7.6×10⁻¹⁶ | 8 | **0** | 0.09937 | 0 |

This is the basis on which the §D/§F verdict is independent of the mismatched-evaluation-point
issue in (2) above: it uses only the one rung-pair with an exact, machine-precision-shared `h`, and
it fails identically to the full ladder. **Recommendation for any future lattice adjudication: use
one explicit, shared `h_list` across every rung** (not index-position matching), so no comparison
mixes a lattice change with an evaluation-point change — see §9's added Option E.

### §D summary table (per field)

| field | FULL 34-node: resolved | class_disagree | numeric_disagree | worst margins (s / D_uni / Dq_min) | ACCEPTED-IN-ALL-4-RUNGS: n_common | resolved | numeric_disagree | worst margins (restricted) |
|---|---|---|---|---|---|---|---|---|
| 1.173192 | **NO** | 0/34 | 34/34 | 16.05 / 7.941 / 16.05 | 10/34 | **NO** | 10/10 | 0.003289 / 0.001609 / 0.0503 |
| 2.581023 | **NO** | 2/34 | 34/34 | 180.9 / 62.11 / 180.9 | 7/34 | **NO** | 7/7 | 0.005844 / 0.0001893 / 0.0851 |
| 3.754215 | **NO** | 4/34 | 34/34 | 72.04 / 10.33 / 72.15 | 28/34 | **NO** | 28/28 | 0.1709 / 0.1589 / 0.2811 |
| 2.850000 | **NO** | 1/34 | 34/34 | 659.5 / 267.7 / 659.6 | 8/34 | **NO** | 8/8 | 0.03442 / 0.03442 / 0.09937 |

"Worst margin" = `max(|Δ| − tol)` over every rung-pair and every node at that field (positive ⇒
failed by that amount; `AbsTol=1×10⁻⁶, RelTol=1×10⁻⁴`, prereg §D, cited verbatim). The
**ACCEPTED-IN-ALL-4-RUNGS restricted view** — the interpretive breakdown described above,
restricting the same frozen test to only the nodes that are checker-accepted in every one of the 4
rungs (so the comparison is never contaminated by unconverged, path-dependent transient iterates)
— **still fails the numeric test at 100% of the common-accepted nodes, at every field.**

**Two clearly separate findings emerge, and both matter:**

1. **The discrete stability *classification* is remarkably robust across the ladder.** `class`
   disagreements are rare (0–4 of 34 nodes per field) and occur **only at accept/reject boundary
   nodes** — never inside the interior of an established stable or unstable run. E.g. at F2581023,
   the 2 disagreeing nodes are `j=26,27` (`h≈0.0066,0.0082`), exactly the transition zone between
   the low-h unconverged region and the stable tail. On the **accepted-in-all-4-rungs** restricted
   set, `class_disagree = 0` at **every field** — i.e., whenever all four rungs actually reach a
   checker-accepted state at the same node, they agree unanimously on stable/unstable.
2. **The frozen *numeric* convergence test (`s`, `D_uni`, `Dq_min` within `1×10⁻⁶+1×10⁻⁴·max`)
   fails everywhere, including on that same accepted-in-all-4-rungs set.** Concrete example: node
   `j=25` at F1173192 (h=0.0061959, deep in the interior of the reproducible 10-node stable tail):
   `s = 0.3796` (baseline) vs. `0.3791` (unshifted/dp40) vs. `0.3771` (half_step/dp30) vs.
   `0.3804` (half_step/dp40) — a ≈0.9% relative spread, about 90× the frozen 0.01% relative
   tolerance. At F3754215, node `j=1` (h=0, the start of the 18-node contiguous accepted-unstable
   run described in §6): `s = -0.4192` (baseline) vs. `-0.5803` (half_step/dp30) — a ≈38% relative
   shift.
3. **The `half_step` offset is a larger source of numeric sensitivity than the `dpRng 30→40`
   increase, at least in the checker-accepted stable region.** Isolating each axis (baseline vs.
   `unshifted/dp40` changes only `dpRng`; baseline vs. `half_step/dp30` changes only the offset)
   over every checker-accepted **stable** node at every field (the region where the comparison is
   physically meaningful — an unconverged node's transient `s` is path/history-dependent and
   isolating an axis there is not informative, consistent with §D's own framing): `half_step/dp30`
   is **consistently the largest deviation from baseline** of the three rungs (e.g. F1173192 node
   25: `Δs=0.0025`, vs. `unshifted/dp40`'s `Δs=0.0005` at the same node; F3754215 node 1:
   `Δs≈0.16`, vs. `unshifted/dp40`'s `Δs≈0.01`), and `unshifted/dp40` most closely tracks baseline
   in nearly every stable-region node checked (one minor exception: F2581023 node 28, where
   `half_step/dp40` is marginally closer, `Δs=0.00001` vs. `unshifted/dp40`'s `Δs=0.00005` — both
   far below the tolerance regardless). Consistent with this, `unshifted/dp40` is the *only* rung
   that ever numerically passed the pairwise §D-style test against baseline among the
   checker-accepted stable nodes checked in §5's existence-bar analysis (F1173192, nodes 30–34
   only; zero passes at the other 3 fields' stable nodes, and zero passes for `half_step/dp30` or
   `half_step/dp40` at any field, including at those same nodes 30–34).

**§D determination: LATTICE/MESH-UNRESOLVED at all 4 physical fields**, under the frozen rule
applied literally (both the full 34-node set and the accepted-only restricted subset fail the
numeric-convergence requirement). Per prereg §F, this is adjudicated first and is a **hard stop
that precedes 3A/3B — neither may be selected while it holds.**

### Full per-node ladder detail, all 4 fields (complete evidence, not excerpted)

<details><summary>F1173192 (1.173192 T) — full 34-node table</summary>

| node j | h(baseline) | class:g1 | class:u40 | class:h30 | class:h40 | class_agree | s:g1 | s:u40 | s:h30 | s:h40 | numeric_agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | unconverged | unconverged | unconverged | unconverged | YES | -0.9377 | -8.061 | -1.2 | -0.1164 | NO |
| 2 | 4.32368e-05 | unconverged | unconverged | unconverged | unconverged | YES | -1.828 | -2.479 | 0.6821 | -2.112 | NO |
| 3 | 5.36542e-05 | unconverged | unconverged | unconverged | unconverged | YES | -1.361 | -1.387 | -1.956 | -0.6237 | NO |
| 4 | 6.65815e-05 | unconverged | unconverged | unconverged | unconverged | YES | -1.249 | -1.898 | -8.372 | -0.3525 | NO |
| 5 | 8.26235e-05 | unconverged | unconverged | unconverged | unconverged | YES | -3.416 | -0.862 | -1.325 | -1.956 | NO |
| 6 | 0.000102531 | unconverged | unconverged | unconverged | unconverged | YES | -1.766 | -0.5569 | -1.892 | -2.636 | NO |
| 7 | 0.000127234 | unconverged | unconverged | unconverged | unconverged | YES | -1.346 | 0.3006 | -1.612 | -0.7761 | NO |
| 8 | 0.00015789 | unconverged | unconverged | unconverged | unconverged | YES | -3.845 | -0.871 | -2.979 | -2.476 | NO |
| 9 | 0.000195931 | unconverged | unconverged | unconverged | unconverged | YES | -9.488 | -4.094 | -1.201 | 0.1125 | NO |
| 10 | 0.000243138 | unconverged | unconverged | unconverged | unconverged | YES | -0.7937 | -0.4818 | -0.2036 | -7.187 | NO |
| 11 | 0.00030172 | unconverged | unconverged | unconverged | unconverged | YES | -0.8533 | -1.298 | -0.8152 | -0.5401 | NO |
| 12 | 0.000374415 | unconverged | unconverged | unconverged | unconverged | YES | -0.931 | -0.9917 | -0.5724 | -0.8318 | NO |
| 13 | 0.000464626 | unconverged | unconverged | unconverged | unconverged | YES | -0.8478 | -0.3593 | -0.1399 | -0.4968 | NO |
| 14 | 0.000576572 | unconverged | unconverged | unconverged | unconverged | YES | -1.229 | -0.7044 | -0.3396 | -2.362 | NO |
| 15 | 0.00071549 | unconverged | unconverged | unconverged | unconverged | YES | -0.584 | -0.4724 | -0.4823 | -0.3994 | NO |
| 16 | 0.000887879 | unconverged | unconverged | unconverged | unconverged | YES | -0.4635 | -1.534 | -0.05531 | -0.3456 | NO |
| 17 | 0.0011018 | unconverged | unconverged | unconverged | unconverged | YES | -0.2116 | 0.02927 | -0.1661 | -0.7606 | NO |
| 18 | 0.00136727 | unconverged | unconverged | unconverged | unconverged | YES | -0.5423 | -0.8481 | -0.3565 | -0.6483 | NO |
| 19 | 0.00169669 | unconverged | unconverged | unconverged | unconverged | YES | -15.3 | 0.7428 | -12.2 | -7.014 | NO |
| 20 | 0.00210549 | unconverged | unconverged | unconverged | unconverged | YES | -0.1018 | -0.2 | -0.296 | -0.2725 | NO |
| 21 | 0.00261279 | unconverged | unconverged | unconverged | unconverged | YES | -0.08753 | 0.00969 | -0.05188 | -0.1587 | NO |
| 22 | 0.0032423 | unconverged | unconverged | unconverged | unconverged | YES | -0.0341 | 0.904 | -0.8254 | -0.7217 | NO |
| 23 | 0.0040235 | unconverged | unconverged | unconverged | unconverged | YES | -0.1047 | -0.3389 | -2.179 | -0.07722 | NO |
| 24 | 0.00499291 | unconverged | unconverged | unconverged | unconverged | YES | 0.08914 | 0.3009 | 0.1447 | -1.548 | NO |
| 25 | 0.00619589 | stable | stable | stable | stable | YES | 0.3796 | 0.3791 | 0.3771 | 0.3804 | NO |
| 26 | 0.00768871 | stable | stable | stable | stable | YES | 0.6042 | 0.6039 | 0.603 | 0.605 | NO |
| 27 | 0.00954122 | stable | stable | stable | stable | YES | 0.7429 | 0.7428 | 0.7422 | 0.7434 | NO |
| 28 | 0.0118401 | stable | stable | stable | stable | YES | 0.8303 | 0.8302 | 0.8298 | 0.8306 | NO |
| 29 | 0.0146928 | stable | stable | stable | stable | YES | 0.8875 | 0.8875 | 0.8872 | 0.8877 | NO |
| 30 | 0.0182328 | stable | stable | stable | stable | YES | 0.9256 | 0.9255 | 0.9253 | 0.9257 | **YES** (u40 leg) |
| 31 | 0.0226258 | stable | stable | stable | stable | YES | 0.9508 | 0.9508 | 0.9506 | 0.9509 | **YES** (u40 leg) |
| 32 | 0.0280772 | stable | stable | stable | stable | YES | 0.9674 | 0.9674 | 0.9673 | 0.9675 | **YES** (u40 leg) |
| 33 | 0.0348421 | stable | stable | stable | stable | YES | 0.9782 | 0.9782 | 0.9781 | 0.9782 | **YES** (u40 leg) |
| 34 | 0.0432368 | stable | stable | stable | stable | YES | 0.9852 | 0.9852 | 0.9851 | 0.9852 | **YES** (u40 leg) |

(`numeric_agree` column shown here is the **4-rung-simultaneous** verdict per `invz_task2_ladder_ok`
(NO everywhere in the full table above by construction, since `h30`/`h40` never pass at this
field) — the "YES (u40 leg)" annotation on nodes 30–34 records the *pairwise*
baseline-vs-`unshifted/dp40` sub-result used in §5's existence-bar test, not a 4-way pass.)
</details>

<details><summary>F2581023 (2.581023 T) — full 34-node table</summary>

| node j | h(baseline) | class:g1 | class:u40 | class:h30 | class:h40 | class_agree | s:g1 | s:u40 | s:h30 | s:h40 | numeric_agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | unconverged | unconverged | unconverged | unconverged | YES | -23.08 | -10.17 | -9.943 | -3.365 | NO |
| 2 | 3.70666e-05 | unconverged | unconverged | unconverged | unconverged | YES | -0.1476 | -2.831 | -181 | -18.13 | NO |
| 3 | 4.59973e-05 | unconverged | unconverged | unconverged | unconverged | YES | -2.461 | -18.57 | -2.333 | -1.012 | NO |
| 4 | 5.70798e-05 | unconverged | unconverged | unconverged | unconverged | YES | -11.1 | -2.393 | -16.77 | 0.5454 | NO |
| 5 | 7.08325e-05 | unconverged | unconverged | unconverged | unconverged | YES | -0.2385 | -4.756 | -4.189 | -36.91 | NO |
| 6 | 8.78987e-05 | unconverged | unconverged | unconverged | unconverged | YES | -0.2774 | -1.815 | -0.5226 | -1.957 | NO |
| 7 | 0.000109077 | unconverged | unconverged | unconverged | unconverged | YES | -2.023 | -2.074 | -2.853 | -1.564 | NO |
| 8 | 0.000135358 | unconverged | unconverged | unconverged | unconverged | YES | -4.212 | -3.562 | -2.224 | -5.677 | NO |
| 9 | 0.00016797 | unconverged | unconverged | unconverged | unconverged | YES | -8.06 | -37.68 | -1.652 | -3.782 | NO |
| 10 | 0.000208441 | unconverged | unconverged | unconverged | unconverged | YES | -8.083 | -3.17 | -0.04559 | -0.6459 | NO |
| 11 | 0.000258662 | unconverged | unconverged | unconverged | unconverged | YES | -4.638 | -10.41 | -4.548 | -1.232 | NO |
| 12 | 0.000320983 | unconverged | unconverged | unconverged | unconverged | YES | -2.68 | -4.286 | -3.123 | -3.276 | NO |
| 13 | 0.00039832 | unconverged | unconverged | unconverged | unconverged | YES | -28.23 | -1.94 | -60.25 | -1.89 | NO |
| 14 | 0.000494291 | unconverged | unconverged | unconverged | unconverged | YES | -1.236 | -17.82 | -2.756 | -22.05 | NO |
| 15 | 0.000613384 | unconverged | unconverged | unconverged | unconverged | YES | -2.119 | -2.484 | -2.498 | -109.3 | NO |
| 16 | 0.000761172 | unconverged | unconverged | unconverged | unconverged | YES | -0.3372 | -3.606 | -3.307 | -1.102 | NO |
| 17 | 0.000944567 | unconverged | unconverged | unconverged | unconverged | YES | -1.811 | -1.632 | -0.9403 | -1.393 | NO |
| 18 | 0.00117215 | unconverged | unconverged | unconverged | unconverged | YES | -1.887 | -1.437 | -4.742 | -1.034 | NO |
| 19 | 0.00145456 | unconverged | unconverged | unconverged | unconverged | YES | -5.501 | -5.318 | -2.934 | -10.67 | NO |
| 20 | 0.00180502 | unconverged | unconverged | unconverged | unconverged | YES | -4.526 | -2.308 | -2.305 | -5 | NO |
| 21 | 0.00223992 | unconverged | unconverged | unconverged | unconverged | YES | -2.31 | -2.124 | -1.939 | -4.048 | NO |
| 22 | 0.0027796 | unconverged | unconverged | unconverged | unconverged | YES | -1.453 | -1.698 | -2.311 | -0.887 | NO |
| 23 | 0.00344931 | unconverged | unconverged | unconverged | unconverged | YES | -2.356 | -1.073 | -0.9813 | -2.377 | NO |
| 24 | 0.00428038 | unconverged | unconverged | unconverged | unconverged | YES | -1.132 | -0.9311 | -1.205 | -1.057 | NO |
| 25 | 0.00531169 | unconverged | unconverged | unconverged | unconverged | YES | -0.5663 | -0.7961 | -0.7243 | -0.6837 | NO |
| 26 | 0.00659147 | unstable | unstable | unconverged | unstable | **NO** | -0.4771 | -0.4649 | -0.5646 | -0.6192 | NO |
| 27 | 0.00817961 | unconverged | unconverged | unstable | unstable | **NO** | -0.1775 | -0.1438 | -0.2099 | -0.2018 | NO |
| 28 | 0.0101504 | stable | stable | stable | stable | YES | 0.06027 | 0.06022 | 0.05442 | 0.06026 | NO |
| 29 | 0.012596 | stable | stable | stable | stable | YES | 0.2585 | 0.2584 | 0.254 | 0.2583 | NO |
| 30 | 0.0156308 | stable | stable | stable | stable | YES | 0.4286 | 0.4285 | 0.4253 | 0.4285 | NO |
| 31 | 0.0193969 | stable | stable | stable | stable | YES | 0.5759 | 0.5758 | 0.5736 | 0.5758 | NO |
| 32 | 0.0240704 | stable | stable | stable | stable | YES | 0.6984 | 0.6983 | 0.6969 | 0.6984 | NO |
| 33 | 0.0298698 | stable | stable | stable | stable | YES | 0.7941 | 0.794 | 0.793 | 0.794 | NO |
| 34 | 0.0370666 | stable | stable | stable | stable | YES | 0.8641 | 0.864 | 0.8634 | 0.8641 | NO |

</details>

<details><summary>F3754215 (3.754215 T) — full 34-node table</summary>

| node j | h(baseline) | class:g1 | class:u40 | class:h30 | class:h40 | class_agree | s:g1 | s:u40 | s:h30 | s:h40 | numeric_agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | unstable | unstable | unstable | unstable | YES | -0.4192 | -0.4094 | -0.5803 | -0.5607 | NO |
| 2 | 2.67063e-05 | unstable | unstable | unstable | unstable | YES | -0.4192 | -0.4093 | -0.5803 | -0.5607 | NO |
| 3 | 3.31408e-05 | unstable | unstable | unstable | unstable | YES | -0.4191 | -0.4093 | -0.5802 | -0.5606 | NO |
| 4 | 4.11257e-05 | unstable | unstable | unstable | unstable | YES | -0.4191 | -0.4093 | -0.5802 | -0.5606 | NO |
| 5 | 5.10344e-05 | unstable | unstable | unstable | unstable | YES | -0.4191 | -0.4093 | -0.5801 | -0.5606 | NO |
| 6 | 6.33305e-05 | unstable | unstable | unstable | unstable | YES | -0.419 | -0.4092 | -0.5801 | -0.5605 | NO |
| 7 | 7.85893e-05 | unstable | unstable | unstable | unstable | YES | -0.419 | -0.4092 | -0.5799 | -0.5603 | NO |
| 8 | 9.75244e-05 | unstable | unstable | unstable | unstable | YES | -0.4189 | -0.4091 | -0.5797 | -0.5601 | NO |
| 9 | 0.000121022 | unstable | unstable | unstable | unstable | YES | -0.4187 | -0.4089 | -0.5794 | -0.5598 | NO |
| 10 | 0.00015018 | unstable | unstable | unstable | unstable | YES | -0.4184 | -0.4087 | -0.5789 | -0.5594 | NO |
| 11 | 0.000186364 | unstable | unstable | unstable | unstable | YES | -0.418 | -0.4084 | -0.5781 | -0.5587 | NO |
| 12 | 0.000231267 | unstable | unstable | unstable | unstable | YES | -0.4174 | -0.4078 | -0.5769 | -0.5576 | NO |
| 13 | 0.000286988 | unstable | unstable | unstable | unstable | YES | -0.4165 | -0.407 | -0.5751 | -0.5559 | NO |
| 14 | 0.000356134 | unstable | unstable | unstable | unstable | YES | -0.415 | -0.4057 | -0.5724 | -0.5533 | NO |
| 15 | 0.00044194 | unstable | unstable | unstable | unstable | YES | -0.4128 | -0.4037 | -0.5682 | -0.5493 | NO |
| 16 | 0.00054842 | unstable | unstable | unstable | unstable | YES | -0.4094 | -0.4006 | -0.5619 | -0.5433 | NO |
| 17 | 0.000680555 | unstable | unstable | unstable | unstable | YES | -0.4041 | -0.3959 | -0.5524 | -0.5342 | NO |
| 18 | 0.000844526 | unstable | unstable | unstable | unstable | YES | -0.3962 | -0.3888 | -0.5382 | -0.5204 | NO |
| 19 | 0.001048 | unconverged | unstable | unstable | unconverged | **NO** | -0.3845 | -0.3778 | -0.517 | -0.4601 | NO |
| 20 | 0.00130051 | unconverged | unconverged | unconverged | unconverged | YES | -72.48 | -0.4374 | -0.4935 | -0.4896 | NO |
| 21 | 0.00161385 | unconverged | unconverged | unstable | unstable | **NO** | -0.3666 | -0.4352 | -0.5024 | -0.4848 | NO |
| 22 | 0.00200269 | unconverged | unconverged | unconverged | unconverged | YES | -0.488 | -26.41 | -0.4629 | -0.4106 | NO |
| 23 | 0.00248521 | unstable | unconverged | unstable | unstable | **NO** | -0.3971 | -0.39 | -0.4724 | -0.456 | NO |
| 24 | 0.00308399 | unconverged | unconverged | unstable | unstable | **NO** | -0.2301 | -0.3397 | -0.3808 | -0.3692 | NO |
| 25 | 0.00382704 | unstable | unstable | unstable | unstable | YES | -0.4513 | -0.4476 | -0.3075 | -0.2983 | NO |
| 26 | 0.00474912 | unstable | unstable | unstable | unstable | YES | -0.2732 | -0.2736 | -0.2978 | -0.2874 | NO |
| 27 | 0.00589336 | unstable | unstable | unstable | unstable | YES | -0.1993 | -0.1989 | -0.2369 | -0.2276 | NO |
| 28 | 0.0073133 | unstable | unstable | unstable | unstable | YES | -0.1142 | -0.1136 | -0.1144 | -0.1065 | NO |
| 29 | 0.00907535 | stable | stable | stable | stable | YES | 0.01566 | 0.0157 | 0.0106 | 0.01618 | NO |
| 30 | 0.0112619 | stable | stable | stable | stable | YES | 0.1377 | 0.1376 | 0.1324 | 0.1376 | NO |
| 31 | 0.0139754 | stable | stable | stable | stable | YES | 0.2628 | 0.2627 | 0.2584 | 0.2626 | NO |
| 32 | 0.0173426 | stable | stable | stable | stable | YES | 0.3877 | 0.3875 | 0.3841 | 0.3874 | NO |
| 33 | 0.021521 | stable | stable | stable | stable | YES | 0.5086 | 0.5084 | 0.5058 | 0.5083 | NO |
| 34 | 0.0267063 | stable | stable | stable | stable | YES | 0.6216 | 0.6213 | 0.6194 | 0.6213 | NO |

</details>

<details><summary>F2850000 (2.850000 T) — full 34-node table</summary>

| node j | h(baseline) | class:g1 | class:u40 | class:h30 | class:h40 | class_agree | s:g1 | s:u40 | s:h30 | s:h40 | numeric_agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | unconverged | unconverged | unconverged | unconverged | YES | -0.03184 | -3.949 | -3 | -1.904 | NO |
| 2 | 3.51347e-05 | unconverged | unconverged | unconverged | unconverged | YES | -1.416 | -0.9661 | -5.546 | -660.5 | NO |
| 3 | 4.35999e-05 | unconverged | unconverged | unconverged | unconverged | YES | -0.9992 | -0.8202 | -0.4349 | -43.86 | NO |
| 4 | 5.41048e-05 | unconverged | unconverged | unconverged | unconverged | YES | -1.449 | -1.143 | -0.7315 | -1.043 | NO |
| 5 | 6.71407e-05 | unconverged | unconverged | unconverged | unconverged | YES | -0.03822 | -1.764 | -15.93 | -1.32 | NO |
| 6 | 8.33174e-05 | unconverged | unconverged | unconverged | unconverged | YES | -0.4756 | -15.14 | -79.11 | -1.005 | NO |
| 7 | 0.000103392 | unconverged | unconverged | unconverged | unconverged | YES | -0.9309 | -1.026 | -1.681 | -2.641 | NO |
| 8 | 0.000128303 | unconverged | unconverged | unconverged | unconverged | YES | -0.6388 | -6.19 | -0.5502 | -0.7233 | NO |
| 9 | 0.000159216 | unconverged | unconverged | unconverged | unconverged | YES | -9.923 | -1.362 | -7.148 | -0.3427 | NO |
| 10 | 0.000197577 | unconverged | unconverged | unconverged | unconverged | YES | -1.27 | -0.9631 | -5.743 | -266.2 | NO |
| 11 | 0.00024518 | unconverged | unconverged | unconverged | unconverged | YES | -0.5394 | -1.182 | -1.615 | -4.643 | NO |
| 12 | 0.000304254 | unconverged | unconverged | unconverged | unconverged | YES | -0.9748 | -10.99 | -2.676 | -1.563 | NO |
| 13 | 0.00037756 | unconverged | unconverged | unconverged | unconverged | YES | -1.246 | -0.5421 | -0.4055 | -1.359 | NO |
| 14 | 0.000468528 | unconverged | unconverged | unconverged | unconverged | YES | -1.963 | -0.1509 | -35.07 | -1.216 | NO |
| 15 | 0.000581414 | unconverged | unconverged | unconverged | unconverged | YES | -0.2877 | -1.297 | 0.05449 | -1.707 | NO |
| 16 | 0.000721499 | unconverged | unconverged | unconverged | unconverged | YES | -1.764 | -8.647 | -2.071 | -2.833 | NO |
| 17 | 0.000895335 | unconverged | unconverged | unconverged | unconverged | YES | -2.161 | -1.714 | -3.03 | -0.8 | NO |
| 18 | 0.00111106 | unconverged | unconverged | unconverged | unconverged | YES | -1.373 | -1.486 | -4.207 | -2.831 | NO |
| 19 | 0.00137875 | unconverged | unconverged | unconverged | unconverged | YES | -6.503 | -2.382 | -3.731 | -1.395 | NO |
| 20 | 0.00171094 | unconverged | unconverged | unconverged | unconverged | YES | -1.452 | -1.227 | -1.482 | -1.774 | NO |
| 21 | 0.00212317 | unconverged | unconverged | unconverged | unconverged | YES | -1.195 | -1.087 | -0.8873 | -1.84 | NO |
| 22 | 0.00263473 | unconverged | unconverged | unconverged | unconverged | YES | -3.092 | -3.009 | -1.412 | -1.276 | NO |
| 23 | 0.00326953 | unconverged | unconverged | unconverged | unconverged | YES | -12.41 | -1.039 | -1.749 | -1.596 | NO |
| 24 | 0.00405729 | unstable | unstable | unconverged | unconverged | **NO** | -1.019 | -1.02 | -0.769 | -1.386 | NO |
| 25 | 0.00503484 | unconverged | unconverged | unconverged | unconverged | YES | -1.219 | -0.7351 | -0.7608 | -0.8822 | NO |
| 26 | 0.00624792 | unconverged | unconverged | unconverged | unconverged | YES | -0.5542 | -0.4141 | -0.4159 | -0.4278 | NO |
| 27 | 0.00775328 | unstable | unstable | unstable | unstable | YES | -0.289 | -0.2905 | -0.2647 | -0.2561 | NO |
| 28 | 0.00962134 | unstable | unstable | unstable | unstable | YES | -0.02845 | -0.02848 | -0.03571 | -0.02934 | NO |
| 29 | 0.0119395 | stable | stable | stable | stable | YES | 0.1707 | 0.1706 | 0.1656 | 0.1705 | NO |
| 30 | 0.0148162 | stable | stable | stable | stable | YES | 0.342 | 0.3419 | 0.3381 | 0.3418 | NO |
| 31 | 0.0183859 | stable | stable | stable | stable | YES | 0.4939 | 0.4937 | 0.491 | 0.4937 | NO |
| 32 | 0.0228158 | stable | stable | stable | stable | YES | 0.6262 | 0.6261 | 0.6242 | 0.6261 | NO |
| 33 | 0.028313 | stable | stable | stable | stable | YES | 0.7358 | 0.7356 | 0.7344 | 0.7357 | NO |
| 34 | 0.0351347 | stable | stable | stable | stable | YES | 0.8207 | 0.8206 | 0.8198 | 0.8206 | NO |

</details>

---

## 4. §C — reproducibility (cold-vs-multistart, isolated-vs-swept)

Full exported state (`Sigma, K, lambda, K0, Gstat, D_uni, hstar, m_star, r_star`) compared via
`invz_task2_agree` (`RelTol=1×10⁻⁶`, `AbsTol_q=1×10⁻⁸·scale_q`), restricted to node-index pairs
**checker-accepted on both sides being compared** (testing reproducibility *of an accepted state*,
per prereg §C's own framing — an unconverged node has no fixed point to reproduce). `hstar`/`m_star`
/`r_star` are consistently `NaN` on both sides everywhere in this matrix (no cell ever reaches a
root — §1/§2), so they pass trivially via `invz_task2_agree`'s documented `isequaln` convention;
the comparison is therefore carried entirely by `Sigma/K/lambda/K0/Gstat/D_uni`.

| field | cold-vs-multistart: matched pairs | agree | disagree | isolated-vs-swept: matched pairs | agree | disagree |
|---|---|---|---|---|---|---|
| 1.173192 | 10 | **10** | 0 | 10 | **10** | 0 |
| 2.581023 | 8 | **8** | 0 | 8 | **8** | 0 |
| 3.754215 | 29 | **29** | 0 | 29 | **29** | 0 |
| 2.850000 | 9 | **9** | 0 | 9 | **9** | 0 |
| **total** | **56** | **56** | **0** | **56** | **56** | **0** |

**56/56 matched-accepted node-pairs agree on both reproducibility axes, at every field, with zero
disagreements.** The worst-case observed discrepancy across all 112 comparisons was of order
10⁻⁹–10⁻¹⁰ (e.g. F3754215 node 23, cold-vs-seed2, `K` differs by 9.405×10⁻¹⁰; F2850000 node 27,
isolated-vs-swept, `K` differs by 1.116×10⁻¹⁰). **Corrected (§0/C6):** these do NOT sit orders of
magnitude inside the `AbsTol_q` floor — re-derived directly from the results file this pass,
`J0eff = 0.0064244` (meV), so the `K`-quantity absolute floor `AbsTol_q = 1×10⁻⁸·|J0eff| ≈
6.424×10⁻¹¹`. The worst discrepancy, `9.405×10⁻¹⁰`, is **≈14.6× that absolute floor** (and the
second example, `1.116×10⁻¹⁰`, is ≈1.7× it) — both numerically EXCEED `AbsTol_q` alone. Both pass
the combined §C tolerance (`AbsTol_q + RelTol·max(|a|,|b|)`, `RelTol=1×10⁻⁶`) through the
**relative** term instead: at the worst node, `RelTol·max(|a|,|b|) ≈ 1.286×10⁻⁹`, giving a combined
bound `≈1.348×10⁻⁹ ≥ 9.405×10⁻¹⁰`. This remains ordinary `tol_outer=1×10⁻⁸`-level solver noise
relative to `K`'s own scale (a ~10⁻⁶–10⁻⁷ relative deviation), not a physically meaningful
difference — the correction is to the tolerance mechanism cited, not to the conclusion that these
are noise-level agreements. Most comparisons are far tighter still (many at 10⁻¹⁷–10⁻¹⁹, i.e.
floating-point noise).

Additionally, at **every one of the 4 fields**: **0 nodes** where `accepted(swept) ≠
accepted(isolated_cold)`, and **0 nodes** where `class(swept) ≠ class(isolated_cold)`. Isolated and
swept fail (and succeed) at **exactly the same nodes**, everywhere in this matrix.

**Conclusion (§C, narrowed after independent review, §0/C5a): a seed/state-management or
continuation defect is not observed for the two registered seed constructions and the tested
isolated-vs-swept comparison.** Full-state reproducibility across cold vs. multi-start seeds, and
across isolated vs. swept-continuation solves, is total (56/56, 0 disagreements) at every physical
field tested. This is the quantitative confirmation prereg §F itself anticipated ("Task 2 CONFIRMS
this quantitatively via C… rather than re-litigating it").

---

## 5. §E — existence bar

**Method** (an operationalization of "≥2 seeds and ≥2 mesh offsets (per C/D)", applying the C and D
tools exactly as defined, not a looser or tighter substitute — **leg (c) corrected after
independent review, §0/C1**): a node counts toward existence at a field iff (a) it is
checker-accepted **stable** in `g1_swept`, `isolated_cold`, **and** `isolated_seed2`
simultaneously; (b) the §C full-state test passes for **both** cold-vs-multistart and
isolated-vs-swept at that node (the "≥2 seeds" leg); **and** (c) **at least one rung using a grid
offset DIFFERENT from baseline's** — i.e. `half_step/dp30` or `half_step/dp40`, **NOT**
`unshifted/dp40`, which shares baseline's own `unshifted` offset and differs only in `dpRng` —
numerically agrees with the `g1` baseline at that node under the same §D-style pairwise test
(class-identical + `1×10⁻⁶+1×10⁻⁴·max` tolerance on `s`, `D_uni`, `Dq_min` — the "≥2 mesh offsets"
leg: baseline's `unshifted` offset + 1 agreeing `half_step` rung = 2 distinct offsets). **The
original version of this report incorrectly credited `unshifted/dp40` as satisfying leg (c);
corrected here.**

| field | existence-bar candidate nodes (leg c corrected) | MET |
|---|---|---|
| 1.173192 | none — `unshifted/dp40` agrees at nodes 30–34 (§3), but that is the SAME offset as baseline, not a second offset; neither `half_step` rung numerically agrees anywhere in this matrix (§3) | NO |
| 2.581023 | none | NO |
| 3.754215 | none | NO |
| 2.850000 | none | NO |

**§E result (corrected, §0/C1): MET at 0 of 4 fields — below the ≥3-of-4 bar.** Leg (b) (seeds)
passes at **every** stable node at **every** field (direct consequence of §4's 56/56 result), so
the entire bar outcome is carried by leg (c) (mesh offsets). Per §3, **no `half_step` rung ever
numerically agrees with baseline at any checker-accepted node, at any field, anywhere in this
matrix** — so leg (c), read literally (a genuinely different grid offset, not merely a different
`dpRng` at the same offset), never passes, and the strict existence bar is unmet at all 4 fields.
(The `unshifted/dp40` agreement at F1173192 nodes 30–34, §3/§5, remains true and is still the only
rung/field/node combination where any numeric ladder pass occurs against baseline anywhere in the
matrix — it demonstrates `dpRng`-axis stability specifically, not a second mesh offset, and so does
not satisfy leg (c) as frozen.)

**⚠️ Underspecification flagged, not resolved here (§0/C4):** prereg §E requires a field to "yield
a stable (class A) checker-accepted ordered state," but — since every one of the 40 cells' `hstar`
is `NaN` (§1/§2; no cell ever reaches a converged H_MF root) — the prereg text does not itself say
what counts as "the field yielding a state" at the node level. This report operationalizes it as
"the field has at least one stable, checker-accepted **internal** H_MF node" satisfying legs
(a)/(b)/(c) above. That is a defensible reading but not a frozen one. **The 0/4 (and 4/4
class-only, below) figures are conditional on this unfrozen operationalization** and should not be
treated as a frozen-criterion result until the user fixes either a specific endpoint condition or a
continuation-based existence criterion.

**Complementary, explicitly-labeled class-only reading (NOT the rule-satisfying operationalization
above, and NOT the verdict — offered only as context, since §3 showed class agreement is
separately excellent):** if leg (c) is read at the *class* level alone (ignoring §D's numeric
sub-test, and not requiring the agreeing rung to be a distinct offset from baseline), all of nodes
25–34 (10 nodes) qualify at F1173192, nodes 28–34 (7 nodes) at F2581023, nodes 29–34 (6 nodes) at
F3754215, and nodes 29–34 (6 nodes) at F2850000 — i.e. the bar would read **MET at 4/4 fields**
under this looser, non-frozen reading. This discrepancy (0/4 strict vs. 4/4 class-only) is not a
contradiction; it is the same §D numeric-non-convergence finding surfacing again downstream, and is
reported here for completeness — **not as an alternate verdict, and not as evidence that softens
§7's LATTICE/MESH-UNRESOLVED determination.**

---

## 6. Step 4 — density/distribution discriminators (G1 vs. G3 vs. G4 vs. G5)

**⚠️ Substantially rewritten after independent review (§0/C5b).** Before interpreting the table
below, a fact verified directly from source during this correction pass changes what this
comparison can and cannot show: **the ordered-phase EMT map is permutation-invariant over the flat
coupling array, by construction.** `invz_emt_scalar.m:43,48` computes `A = mean(Jf ./ (D +
Jf.*G0), 1)` and `invz_emt_static_ordered.m:47-50` computes `Gq = Gs./(1+(Jf-K0).*Gs)`, `Gbar =
mean(Gq)`, `K0 = mean(Jf.*Gq)/Gbar` — every operation is either elementwise in `Jf` or a `mean`
taken over `Jf`'s flat mode index; **neither function ever references `q`, branch identity, or any
other geometric label.** Consequently the calculation depends on the flat coupling array
`Jnu_flat` **only through its values as an unordered multiset** — any permutation of the same
values, or any other coupling array with the same multiset, produces the same `A`/`K`/`Gq`/`K0`,
**mathematically permutation-invariant, up to floating-point reduction-order effects** (the `mean`
reductions are not associative in IEEE arithmetic, so a reordering can perturb the result at the
roundoff level — the invariance is exact in real arithmetic, not bit-for-bit; refreshed-review
correction). **Geometry (which `q`, which branch, which BZ point a coupling came from) never enters
this map at all.**

| field | source | n_accepted | stable | marginal | unstable | unconverged | note |
|---|---|---|---|---|---|---|---|
| 1.173192 | g1_swept | 10 | 10 | 0 | 0 | 24 | G1 baseline (full 16384) |
| 1.173192 | g3_ds2 | 10 | 10 | 0 | 0 | 24 | density 1/2 (8192), real shape |
| 1.173192 | g3_ds4 | 10 | 10 | 0 | 0 | 24 | density 1/4 (4096), real shape |
| 1.173192 | g3_ds8 | 10 | 10 | 0 | 0 | 24 | density 1/8 (2048), real shape |
| 1.173192 | g4_cardsynth | **18** | **17** | 0 | **1** | **16** | cardinality 16384, non-real shape |
| 1.173192 | g5_histmatch | 10 | 10 | 0 | 0 | 24 | cardinality 16384, real shape |
| 2.581023 | g1_swept | 8 | 7 | 0 | 1 | 26 | G1 baseline |
| 2.581023 | g3_ds2/ds4 | 8 | 7 | 0 | 1 | 26 | node-by-node identical to G1 |
| 2.581023 | g3_ds8 | 8 | 7 | 0 | 1 | 26 | identical counts; two-node boundary shift — nodes 26/27 unstable/unconverged SWAP vs. G1 (verified this pass against the results file: G1 has node26=unstable/node27=unconverged, `g3_ds8` has node26=unconverged/node27=unstable; corrected from "identical", §0/C5c) |
| 3.754215 | g1_swept | 29 | 6 | 0 | 23 | 5 | G1 baseline |
| 3.754215 | g3_ds2 | 29 | 6 | 0 | 23 | 5 | density 1/2 — identical to G1 |
| 3.754215 | g3_ds4 | 31 | 6 | 0 | 25 | 3 | density 1/4 — 2-node boundary shift |
| 3.754215 | g3_ds8 | 28 | 6 | 0 | 22 | 6 | density 1/8 — 1-node boundary shift |
| 2.850000 | g1_swept | 9 | 6 | 0 | 3 | 25 | G1 baseline |
| 2.850000 | g3_ds2/ds4 | 9 | 6 | 0 | 3 | 25 | identical to G1 |
| 2.850000 | g3_ds8 | 11 | 6 | 0 | 5 | 23 | density 1/8 — 2-node boundary shift |
| 2.850000 | g4_cardsynth | **33** | **8** | 0 | **25** | **1** | cardinality 16384, non-real shape |
| 2.850000 | g5_histmatch | 9 | 6 | 0 | 3 | 25 | cardinality 16384, real shape |

**Density alone (G3, real distribution shape, `{1, 1/2, 1/4, 1/8}` cardinality via deterministic
modular decimation):** the classification pattern is essentially invariant down to 1/4 density at
every field, and remains qualitatively unchanged (same overall structure, ±1–2 nodes at the
accept/reject boundary) even at the coarsest 1/8 density (F2581023's `ds8` two-node boundary swap,
above, is exactly this kind of boundary-level effect). **Density is at most a secondary,
boundary-level effect** — consistent with the permutation-invariance fact above: a structured
modular subset leaves the multiset's bulk statistics (and hence `A`/`K0`) nearly unchanged until
the subset becomes coarse enough to perturb a handful of boundary nodes.

**Coupling-multiset sensitivity at fixed cardinality (G1 vs. G4 vs. G5, all exactly 16384
couplings):** **G5** (couplings histogram-matched to the real lattice's own `Jnu` distribution, at
real cardinality) **reproduces G1 exactly** at both tested fields — F1173192: `(10,0,0,24)` = G1's
own `(10,0,0,24)`; F2850000: `(6,0,3,25)` = G1's own `(6,0,3,25)`. **Given the permutation-invariance
fact above, this is EXPECTED BY CONSTRUCTION and is NOT independent evidence** (§0/C5b): `G5` is
explicitly built to match the real coupling multiset's histogram, so at real cardinality it
approximately reconstructs the same multiset the map is provably blind to reordering of —
reproducing `G1` is close to a tautology of the construction, not a discovery about what drives the
pattern. **G4** (the pinned 24-point synthetic fixture tiled to the same cardinality, a multiset
with a different shape/histogram) **diverges sharply from G1 at the same two fields**: F1173192
`(17,0,1,16)` vs. G1's `(10,0,0,24)` — 70% more stable nodes, an unstable node appears where G1 has
none; F2850000 `(8,0,25,1)` vs. G1's `(6,0,3,25)` — the unstable count jumps from 3 to **25** and
the unconverged count collapses from 25 to **1**. This divergence IS informative: it shows the
classification pattern is sensitive to the coupling multiset's composition (not just its
cardinality — G1/G4/G5 all hold cardinality fixed at 16384).

**What this can and cannot conclude (corrected, §0/C5b).** The data supports: the stability
classification is sensitive to the multiset of coupling values (G4 vs. G1/G5 at matched
cardinality). **The data does NOT support, and cannot support, a claim that this is "distribution
shape... not... geometry"** as if the two were competing, separable causes: **geometry never enters
the map** (per the source-verified permutation-invariance above), so no experiment in this matrix —
or any experiment on this map — can isolate "distribution shape" from "geometry" as distinct
drivers. "Coupling-value distribution shape" and "geometry" are not alternative explanations here;
geometry's *only* possible channel of influence on this calculation is by determining *which
multiset of values* gets sampled onto the flat array in the first place.

**This sharpens, rather than weakens, the §3/§D reading:** the grid-offset (`unshifted` vs.
`half_step`) and `dpRng` sensitivity documented in §3 must propagate entirely through the resulting
**changed sampled multiset** of `J(q)` at each rung (a different offset/cutoff samples different
`q`-points, hence a different multiset of coupling values) — there is no other channel available to
this map. **"Lattice/mesh-unresolved" (§7) can therefore be restated precisely: the sampled
coupling multiset is not converged under mesh refinement** (grid offset, and to a lesser extent
`dpRng`), which is exactly the quantity §3 measured failing to converge.

**Step-4 finding, corrected (§0/C5b): the data shows the classification pattern is sensitive to the
sampled coupling multiset's composition (G4-vs-G1/G5 at fixed cardinality; density/`G3` only a
secondary, boundary effect) — but, because the map is provably geometry-blind by construction, this
cannot be further attributed to "distribution shape" as opposed to "geometry"; the two are not
separable causes for this calculation.**

---

## 7. §F classification

Applying the frozen table (prereg §F) **in order**:

**1. LATTICE/MESH-UNRESOLVED — TRIGGERED.** Per §3: the stability class is highly (not perfectly)
reproducible across the grid-offset/`dpRng` ladder, but the frozen **numeric** convergence
requirement on `s`/`D_uni`/`Dq_min` (`|Δ| ≤ 1×10⁻⁶+1×10⁻⁴·max`) **fails at all 4 physical fields**,
including in the restricted view limited to nodes checker-accepted in every rung. Per prereg §F,
this verdict **is adjudicated first and is a hard stop preceding 3A/3B — neither may be selected
while it holds.** **This is the operative classification of this report.**

**The remaining rows of the table are reported below for completeness and to inform how the
lattice question should be resolved — they are explicitly NOT triggered/operative verdicts, since
(1) above already precludes selecting any of them:**

- **Task 3A's specific mechanism is absent from the evidence.** 3A requires a checker-accepted
  stable state at **the same node where the production (swept) iteration fails**, reached by
  isolated/cold/downsampled but not by the dense swept run. §4 found **zero** such nodes anywhere:
  isolated and swept accept/reject **exactly** the same nodes at every field (0/136 discrepancies),
  and full-state agreement is total (56/56). There is no node in this matrix where continuation
  "poisons" a state that isolation would otherwise reach — so even setting the lattice-unresolved
  hard stop aside, the data does not support a pure conditioning-repair (3A) story.
- **Task 3B's signature has real, but incomplete, qualitative evidence, and is blocked twice
  over.** A contiguous run of **18** checker-accepted nodes with `s<-1×10⁻³` exists at F3754215,
  `h∈[0, 0.000845]` (nodes `j=1`–`18`, baseline `g1_swept_F3754215`), and this run's **count and
  approximate location reproduce** across every lattice realization tested at that field (`G3`
  stride2/4/8: 18/20/18 nodes; `G6` unshifted-dp40/halfstep-dp30/halfstep-dp40: 19/19/18 nodes) and
  across both seed variants (§4: all 18 nodes are accepted+agreeing in cold-vs-multistart and
  isolated-vs-swept). Weaker, non-contiguous unstable touches (1–3 nodes) also appear at F2581023
  and F2850000. This is qualitatively the character of evidence 3B wants. **The interval IS
  residual-clean in the checker sense** — all 18 states are checker-accepted — but that is a
  distinct property from lattice convergence, and two separate obstacles block 3B (refreshed-review
  correction: do not conflate "residual-clean" with "lattice-converged"): (a) the interval's
  numerical values are **not lattice-converged** — this same interval **fails the strict §D numeric
  ladder test** that triggered LATTICE-MESH-UNRESOLVED above (e.g. node `j=1` at F3754215:
  `s=-0.4192` baseline vs. `-0.5803` half_step/dp30, a ≈38% relative shift); and (b) **HMF-grid
  (`nH`) persistence remains untested** — the refinement axis prereg §F explicitly requires
  ("persisting under HMF-grid refinement") was never run; every one of the 40 cells has
  `n_nodes=34` (§1). Concisely: *a residual-clean unstable interval exists on each tested
  realization, but its numerical values are not lattice-converged and HMF-grid persistence remains
  untested.* **3B is therefore not established, and is not even a clean, gap-flagged candidate while
  the lattice/mesh question (a) remains open** — both obstacles would need to clear before 3B could
  be assessed.
- **"3B then 3A" and "UNSUPPORTED"** are, likewise, not reachable determinations while (1) holds;
  neither is discussed further as an operative verdict here.

**§F verdict: LATTICE/MESH-UNRESOLVED at all 4 physical fields.** Per the frozen table, the lattice
construction/normalization must be resolved before any solver-conditioning (3A) or
theory/reformulation (3B) path is entertained.

---

## 8. What this report does NOT establish

- It does **not** establish 3A, 3B, "3B then 3A", or UNSUPPORTED — all are precluded by the §D
  hard stop reached in §7. If any of them is reached in a later pass, "UNSUPPORTED" (should it ever
  be reached) would mean *unsupported by the current derivation and completed search*, **not** a
  proof that no mathematical solution exists — that qualifier is not reached here, but is recorded
  per the brief's own requirement.
- It does **not** determine whether the measured `dpRng`/grid-offset numeric sensitivity (§3) is a
  genuine, expected property of a conditionally-convergent dipolar lattice sum at these cutoff
  radii/offsets, or a fixable artifact of the specific `dpRng∈{30,40}` / `{unshifted,half_step}`
  choices tested. Distinguishing these requires further lattice-convergence work (e.g. extending
  the `dpRng` ladder, or examining the sum's convergence order) that this report does not perform.
- It does **not** test H_MF-grid (`nH`) refinement **at all** — zero cells in this matrix vary
  `nH`; every cell has `n_nodes=34`. This is a completely unexercised axis, independent of the
  `dpRng`/offset axis that *was* tested (and failed) in §3.
- It does **not** test BZ grid-**size** refinement **at all** (disclosed after independent review,
  §0/C2) — every cell in this matrix (`G1` through `G6`) holds the BZ grid fixed at `16×16×16`,
  confirmed from each cell's own stored construction provenance (§3); `G3`'s downsampling is a
  structured subset of that same fixed grid, not an independent size refinement. This is a second
  completely unexercised axis alongside `nH` above — the lattice ladder actually run (§3) covers 2
  of the 3 axes prereg §D names ("grid size, grid offset, `dpRng`"), not all 3.
- It does **not** distinguish "coupling-value distribution shape" from "lattice geometry" as
  drivers of the classification pattern (§6, corrected §0/C5b) — the ordered EMT map is
  permutation-invariant over the flat coupling multiset by construction (verified from source,
  `invz_emt_scalar.m`/`invz_emt_static_ordered.m`), so geometry has no channel of influence on this
  calculation except by determining which multiset of values gets sampled. §6 supports multiset
  sensitivity; it cannot and does not establish "shape" as opposed to "geometry."
- It does **not** propose, evaluate, or apply any alternative numeric tolerance for §D/§E — the
  frozen `1×10⁻⁶+1×10⁻⁴·max` bound was applied exactly as written throughout. Whether that bound is
  the right one for a quantity as lattice-sensitive as this static pole observable is a question
  for the user, not something this report (or its author) may unilaterally re-open.
- It does **not** modify, repair, retune, or re-run any solver, checker, harness, or driver.
  Read-only analysis of the already-checkpointed 40-cell matrix only.
- It does **not** run the fast MATLAB test suite (`runtests('invz_projected/tests')`) — no
  production or test file was touched by this task, so there is nothing to regress; the suite was
  not re-run to avoid claiming an unnecessary verification step.

---

## 9. Recommended next options (for the USER to choose — not a decision taken here)

The frozen table requires the lattice/mesh question to be resolved before 3A/3B are entertained.
Options, with their evidentiary basis:

**Option A — resolve the lattice/mesh question directly (the frozen table's own prescribed next
step).** Extend the `dpRng` ladder further (e.g. `dpRng∈{50,60,...}`) to test whether
`s`/`D_uni`/`Dq_min` asymptote as the real-space dipole cutoff grows, and/or investigate the
`half_step` offset construction specifically — §3 found it, not the `dpRng` increase, is the
dominant source of the observed numeric sensitivity among the checker-accepted stable nodes
checked (the `unshifted/dp40` rung is the only one that ever numerically passed the §D bound
against baseline among those nodes, §3/§5). This directly targets what triggered §7's verdict.

**Option B — reconsider whether the current frozen §D numeric tolerance is calibrated for this
observable.** §D's `1×10⁻⁴` relative bound was set before any lattice-offset data existed; §3's
measurements (a static-pole scalar shifting by up to ~38% relative under a half-BZ-step offset,
even in the interior of an otherwise-robust, checker-accepted, class-reproducible stable/unstable
run) may indicate the bound is tighter than a finite-`nq`, conditionally-convergent dipolar sum can
be expected to satisfy under this offset construction — or may indicate a genuine defect in the
offset construction itself; the evidence here does not distinguish these. **This would be
re-opening a frozen threshold, which is explicitly the user's call alone**, not a change this
report or its author may make unilaterally.

**Option C — run the untested HMF-grid (`nH`) refinement axis.** Independent of the lattice-offset
question, this axis is required evidence for 3B under prereg §F and is currently completely
unexercised (0 of 40 cells). Worth running regardless of how Option A/B resolve, since it will be
needed either way if the lattice question clears.

**Option D — proceed informally/exploratory on the qualitative findings, accepting
lattice-unresolved as an open risk.** §4's total seed/continuation reproducibility (56/56) and
§3's robust *class*-level ladder agreement, together with §6's coupling-multiset sensitivity
finding (§0/C5b — **not** the "clean distribution-shape, not geometry" reading originally stated;
see §6's corrected framing) and the reproducible 18-node unstable interval at F3754215 (§7), are
genuinely informative and may be enough for the user to form an independent physics judgment.
**This is explicitly NOT what the frozen table recommends** ("No path is selected because it
produces fuller plots" — prereg §F; lattice-unresolved precedes 3A/3B by design) and is noted here
only because the user, not this report, owns the decision to follow the process as specified or to
override it with that judgment.

**Option E — run the independent reviewer's proposed exact-`h` lattice audit (added §0; an
elaboration of Option A, offered as a concrete design, not a decision taken here).** Re-run the
§D-style lattice ladder with the methodological fixes C2/C3 identified: (1) **one explicit, shared
`h_list`** evaluated identically on every rung (not index-position matching, so no comparison mixes
a lattice change with an evaluation-point change); (2) **actual BZ grid-size values** in the
ladder, e.g. `12³`, `16³`, `20³` (currently completely untested, §0/C2); (3) an **extended `dpRng`
ladder**, e.g. `dpRng ∈ {30, 40, 50, 60}`; (4) **both grid offsets** (`unshifted`, `half_step`) at
every grid-size/`dpRng` combination, with explicit, identical Γ handling; (5) comparisons
restricted to **checker-accepted states only**, under the frozen §A/§D tolerances, unchanged.
**Suggested pilot fields (before running all 4): `3.754215` T** (the reproducible 18-node unstable
interval and the largest observed cross-rung sensitivity, §3/§7) **and `1.173192` T** (the clean,
fully-accepted stable tail) — extending to all 4 physical fields only if the pilot does not
immediately expose a construction defect. If the observables converge under this design, the 3B
unstable-interval question (§7) can be reassessed (still pending the separate, untested `nH` axis,
Option C); if they do not converge, the `half_step`/Γ-normalization and real-space dipolar-sum
convergence should be investigated before any solver or theory path.

No path above has been selected. This report ends at the decision gate.

---

## 10. Main technical concern (refreshed review): the legacy BZ quadrature

**Verified from source, this pass:**
- `invz_bz_couplings.m:14` builds the coupling grid via
  `qVec_generator(ion.a, 'mode','grid', 'grid',grid, 'range',[-0.5 0.5], ...)` with the **default
  `endpoint=true`** — i.e. `linspace(-0.5, 0.5, N)` per axis, endpoint-inclusive.
- `qVec_generator.m:189–193` documents this convention's own hazard verbatim: for a
  reciprocal-periodic range like `[-0.5,0.5]` "the two boundary faces are duplicates."
- **Numerically confirmed:** a 16³ endpoint-inclusive grid contains **3375 = 15³ distinct periodic
  points**, i.e. **17.6 % redundant samples** (721 duplicate boundary-face points), computed by
  wrapping `mod(q+0.5,1)` and counting unique rows.
- A **half-open** alternative already exists in the same constructor (`endpoint=false` →
  `lo + (0:N-1)/N·(hi−lo)`, N distinct points per axis, no duplicate face), currently unused by the
  production coupling path.

**Why this matters mechanistically.** §6 established that the ordered EMT map depends on the flat
couplings only through a `mean` over modes (permutation-invariant). A `mean` over a grid with 721
duplicated boundary-face modes is therefore a **weight-distorted quadrature** — those faces are
over-counted. Shifting the grid by a half-step changes *which* modes land on the duplicated faces,
so the effective quadrature weights move between rungs; this feeds straight into `K0`/`Gstat`/`s`/
`D_uni`, which is exactly the offset sensitivity §3 measured. Consequently the current
`unshifted`-vs-`half_step` comparison is **not yet a clean BZ quadrature-convergence test**: it
mixes genuine sampling-phase change with duplicate-face reweighting, differing Γ/cardinality
handling, and (on the `dp40` rungs) evaluation-point drift. This is a plausible root cause of the
LATTICE/MESH-UNRESOLVED verdict — and it must be resolved **before** the frozen §D tolerance is
reopened (only after removing the duplicated-face/Γ/cardinality/evaluation-point confounders is a
tolerance judgement physically interpretable).

## 11. Immediate action plan (proposed — HARD STOP still active; not yet executed)

The hard stop remains in force; nothing below runs until the user approves it and settles the two
theory/threshold decisions flagged. The plan follows the refreshed review's recommended order, with
the cheapest diagnostic first.

**Phase 0 — report corrections (this commit).** The two wording fixes (W1/W2) and this
finding (T1/§10). Done.

**Phase 1 — coupling-only BZ-quadrature/Γ audit (NO ordered solves; cheap, read-only-ish).** Build a
pre-registered, tested coupling-grid audit that, for each convention/offset/`dpRng`, checks: point
uniqueness, cardinality, Γ count, reciprocal periodicity `J(q+G)=J(q)`, and coupling-multiset
statistics — with **no H_MF/ordered solve at all**. Adopt: (a) a **half-open** periodic grid
(`endpoint=false`, no duplicate ±0.5 face); (b) shifted coordinates **wrapped back into one BZ**;
(c) explicit **quadrature weights summing to 1** for every convention; (d) a single frozen **Γ
policy**. Test the baseline plus the eight `{0,½}³` subcell offsets (their union = a refined 2N
grid, separating true resolution convergence from single-offset dependence). Reference
implementation pattern for half-open grids + explicit conventions + Γ-exclusion provenance:
`invz_tensor/invzt_qgrid.m`. **This phase is additive diagnostic code; it does not touch the
production coupling path or the frozen prereg.**

**Phase 2 — freeze the corrected convention + weights**, then **recompute `Bc_PM`** under it and
**re-freeze the four physical fields** if the convention shifts them (they are derived from `Bc_PM`).

**Phase 3 — exact-`h` ordered lattice audit.** ONE explicit shared `h_list` for **every** rung (the
decisive fix for §3/C3); real BZ **grid sizes** {12³, 16³, 20³} (closes the §0/C2 gap);
`dpRng ∈ {30,40,50}`; both offsets (or the eight subcell offsets) with **identical** Γ handling;
compare **only checker-accepted** states under the frozen class + numeric tolerances. **Pilot
first**: F3754215 (interval start/end + one stable-tail node) and F1173192 (onset + interior);
expand to all four fields only if the pilot shows a coherent grid construction.

**Phase 4 — HMF-`nH` refinement** (the axis §F requires for 3B; independent of Phase 3).

**Phase 5 — reassess 3B** only after **both** lattice and `nH` convergence are established. If
`dpRng` convergence is still slow after the quadrature fix, the durable route is an Ewald /
convergence-accelerated dipolar sum with Lorentz + demagnetization separated analytically (real-
space cutoff growth alone may not converge a conditionally-convergent sum).

**Decisions required from the user before Phase 1 executes:**
1. **Go-ahead** to build the Phase-1 coupling-only audit (cheap; no ordered solves).
2. **The Γ policy** — a *theory* decision I will not make unilaterally: use the complete quadrature
   for the EMT average while treating the uniform Γ pole separately in the `D_uni`/`Dq` diagnostic,
   **or** drop Γ and renormalize weights consistently in every rung. (I recommend the former unless
   the derivation requires dropping Γ.)
3. **The §D tolerance stays frozen for now** — per the refreshed review, do *not* reopen it until
   Phase 1–3 remove the confounders. No action needed unless you disagree.

I recommend approving Phase 1 (and settling decision 2) as the next step: it is cheap, it directly
tests whether the duplicate-face/Γ construction is the culprit, and every later phase depends on its
outcome.

---

## Appendix: reproducibility

Every quantitative claim above is reproducible from `.superpowers/sdd/task2_matrix_results.mat`
(git-ignored, 113 MB) via the cell ids, node ids (`j`/`id`, identical), and `h` values cited
throughout. The frozen tools used — `invz_task2_classify.m`, `invz_task2_agree.m`,
`invz_task2_ladder_ok.m` — were called directly, unmodified; no threshold in them was edited for
this report.
