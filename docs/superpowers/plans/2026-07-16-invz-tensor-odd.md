# invz_tensor — Full-Tensor 1/z ODD Branch (Appendix A: A0+A1 + Tier 2) Implementation Plan

> **SUPERSEDED 2026-07-17** by `2026-07-17-invz-tensor-full.md` (A0–A4). That plan absorbs the
> A0/A1 layer (updated to consume the now-implemented projected ODD functions as parity targets),
> fixes two physics defects here — the Task-8 additive "(a)/(b)/(c)" ΔTc decomposition (ill-posed
> per the uniform-shift theorem, framework §11.3/§11.6) and the Task-7 K formula (must be defined
> against the renormalized propagator: K = 1/G_loc − 1/G̃₀, G̃₀ = −[χ_dom/(1+Σ)+χ_rest]_cc) —
> drops the Tier-2-additive re-implementation (Tier-2 physics re-enters via the A3 mixed
> cumulants), and adds A2 (matrix EMT) + A3 (tensor Σ) + A4 (state-space ladder). Do not execute
> this document.

> **STATUS: DEFERRED (user decision, 2026-07-16).** The active ODD stream is the
> main body of `odd_implementation_plan.html` executed against `invz_projected/`
> — see `2026-07-16-invz-odd-mainbody.md`. This Appendix-A plan is kept as the
> design record for the future `invz_tensor/` branch; its `invz_common/` refactor
> (Task 1) should be revisited when that branch starts. Do not execute now.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the off-diagonal dipolar (ODD) physics of `odd_implementation_plan.html` in a new `invz_tensor/` branch via the plan's Appendix A staged route — A0 (cached `[12,12,nq]` tensor builder + tensor-RPA parity layer) and A1 (projected-1/z bridge: dominant-sector Σc inside the full 12×12 RPA) — plus the Tier-2 internal-field renormalization kept additive per the Appendix scope box, with branch-shared functions factored into a new `invz_common/`.

**Architecture:** Three layers. (1) `invz_common/` receives the branch-agnostic single-ion + scalar-Σ engine moved out of `invz_projected/` (behavior-neutral: everything is resolved via the MATLAB path; existing tests keep living in `invz_projected/tests`). (2) `invz_tensor/` gets `invzt_`-prefixed functions only (both branches sit on the path simultaneously in parity tests — name collisions are forbidden): tensor builder, tensor RPA, E1/E4/E5 mediated-coupling diagnostics, transition-mask χ0 split, the A1 bridge point solver, critical finders, and Tier-2 field-variance dressing. (3) Physics measurements are *reported, never tuned*; the only pass/fail physics gates are ODD-off regressions against the published `invz_projected` benchmarks and internal identities (Parseval, Schur, PSD, sum rule). A2/A3 (matrix EMT, tensor cumulants) stay out of scope, documented as deferred.

**Tech Stack:** MATLAB R2025a (`pagemtimes`/`pagemldivide` available), `functiontests(localfunctions)` suites, existing root lattice machinery (`MF_dipole.m`, `exchange.m`, `qVec_generator.m` stay at repo root).

## Global Constraints

- Repo root: `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion/`. The path contains spaces — ALWAYS quote it in shell commands.
- Fast suites (run from repo root; EVERY task ends with BOTH green):
  `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"`
  and, from Task 3 on,
  `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_tensor/tests'); disp(r); assertSuccess(r)"`
  Slow gates: prefix `INVZ_SLOW=1` (tests use `assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), ...)`).
- Non-negotiables (ODD plan §0): (i) `invz_projected` results are bit-for-bit untouched — this feature only MOVES its files (Task 1) and never edits their logic; (ii) never edit published benchmark values in existing tests; (iii) physics acceptance criteria marked *report* are measurements — record the number, do not tune to a target.
- Naming: every new function in `invz_tensor/` uses the `invzt_` prefix. No file may shadow a name in `invz_projected/`, `invz_common/`, or repo root.
- Tensor index convention (LOCKED): the 12-dim composite index is sublattice-major, Cartesian-minor: `i = 3*(s-1) + mu`, s = 1..4 sublattices, mu = 1(a), 2(b), 3(c). The (s,s') block of `Jt` is the 3×3 Cartesian coupling `J^{mu nu}_{ss'}(q)`. Local response embeds as `kron(eye(4), chi0_3x3)`.
- Units and sign conventions (must match `invz_projected`): energies meV; `C = invz_const()`; ferromagnetic J > 0; `J^{mu nu}_{ss'}(q) = -C.gfac*dip(mu,nu,s,s') + sign(ion.J12)*ex(mu,nu,s,s')` with `dip` from `MF_dipole` (Å⁻³, gfac carries (gL·muB)²·mu0/4π = 0.08388 meV·Å³) and `ex` from `exchange`; Lorentz cavity `lorz = 4*pi/(3*ion.Vc)*C.gfac` added at Γ-equivalent q (via `invz_is_gamma_equiv`) on Cartesian-DIAGONAL entries only (`mu == nu`), all 16 sublattice pairs — exactly the tensor generalization of the scalar broadcast in `invz_jq_modes`. Susceptibility χ = −G, χ0 in meV⁻¹, `J*chi0` dimensionless.
- Demagnetization: the tensor layer is intrinsic-only. `invzt_jq_tensor` must `error('invzt:demag', ...)` if `ion.demag ~= 0` (mirrors `invz_chi_tensor_ref`).
- Published anchors (NEVER rescale): `info.Jcc0_dipole = 6.821e-3` meV, `info.Jaa0_dipole = 3.912e-3` meV, `info.Jcc0 = 6.421e-3` meV (RelTol 0.03, R2007 eq. 4); no-ODD 1/z: Σc(0) Richardson(12³,24³) = 0.2980 (published 0.3004, AbsTol 0.006), Tc(0) = 1.743 K vs published 1.74; DS2023 Suppl. Table I lattice sums (a = 5.175 Å): B_xz,xz = B_yz,yz = 36.73·a⁻⁶, B_zz,zz = 17.93·a⁻⁶; DS2023 3-state params Δ = 11.5 K, ρ = 2.34 (χ⊥(3-state) = 2ρ²/Δ ≈ 11.05 meV⁻¹); expected full-CF χ⊥(0) ≈ 16–17 meV⁻¹ at (1.53 K, 0 T); expected d ≈ 0.3–0.5 μeV (report-band, not a gate).
- ODD-plan pitfalls that BIND every task (plan §8 + A0): (i) never grid-extrapolate the q→0 ODD blocks — the only q=0 modification is the explicit constant −d on the uniform coupling (E5); (ii) ODD small-q decay assertions hold on HIGH-SYMMETRY on-axis rays only — the q→0 limit is direction-dependent (macroscopic term ∝ q̂_c q̂_α, element scale up to ≈ 1.8 μeV on tilted rays); (iii) never import DS2023's 0.805 longitudinal rescaling or their perturbative hyperfine treatment; (iv) the ODD blocks are strictly c↔(a,b): assert no aa/bb leakage into E1; (v) `(gL·muB)²` slips enter δJ SQUARED — the DS2023 geometry sums + the 6.821 μeV anchor are the mandatory unit guards.
- Anchors fixture: Task 2 creates `invz_tensor/tests/invzt_anchors.m` (a function returning a struct of controller-verified digits measured on THIS tree). Later tasks' tests reference `invzt_anchors()` fields instead of hard-coding unpinned digits. Never edit a pinned anchor to make a test pass — a mismatch is a finding to escalate.
- Cache discipline: `invz_tensor/cache/`, key format `jqt<schema>_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat`, `pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]` (trailing 1 = schema v1), file stores `pkey`+`qvec` and the loader `isequal`-verifies both before trusting a hit (cache-hardening precedent from the scalar branch). Never share cache files with `invz_projected/cache`. Library functions stay serial; any parallelism lives in drivers, and all lattice sums are computed before entering a `parfor` (workers do no disk I/O).
- Fast tests use small grids (n ≤ 8 per axis, dpRng ≤ 20, `'cache', false` unless the test IS the cache test); `qVec_generator`'s grid param is `'grid'` (NOT `'size'`) and it `fprintf`s diagnostics — wrap calls in `evalc` where output must be clean, or build uniform meshes directly with `ndgrid` in tests.
- Test fixture non-convergence island (scalar branch, inherited physics): avoid (T = 0.31 K, B ∈ [1, 2] T) for convergence-gated fixtures.
- MATLAB test boilerplate for `invz_tensor/tests` (every new test file):

```matlab
function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion + scalar-Sigma engine
addpath(fullfile(here, '..', '..', 'invz_projected'));  % parity targets (invz_jq_modes, invz_solve_point, ...)
end
```

- Commit style: conventional prefix + scope, e.g. `refactor(invz): ...` (Task 1), `feat(invzt): ...`, `test(invzt): ...`, `docs(invzt): ...`; every commit ends with the trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Reporting: each task that produces physics numbers appends a dated entry to `docs/ODD-LOG.md` (Task 2 creates it): what was implemented/measured, test status, headline numbers.

---

### Task 1: `invz_common/` — factor the branch-shared engine out of `invz_projected/`

**Files:**
- Create: `invz_common/` (git mv targets below), `invz_common/README.md`
- Move (git mv from `invz_projected/`): `getf.m`, `invz_const.m`, `invz_ion.m`, `stevens_ops.m`, `invz_cfrot.m`, `invz_field_vec.m`, `invz_single_ion.m`, `invz_chi0z.m`, `invz_check_transverse_mf.m`, `invz_matsubara.m`, `invz_is_gamma_equiv.m`, `invz_twolevel.m`, `invz_g.m`, `invz_lambdas.m`, `invz_sigma.m` — 15 files, exactly this list (they are what `invz_tensor` consumes; `invz_emt_scalar`, `invz_twolevel_ordered`, `invz_jq_modes` and everything else stays projected-specific).
- Modify: every test file in `invz_projected/tests/*.m` — in `setupOnce`, add `addpath(fullfile(here, '..', '..', 'invz_common'));` after the existing addpath lines. Also any file found by the greps below.
- Test: the existing `invz_projected/tests` suite is the acceptance test — no new tests.

**Interfaces:**
- Produces: the 15 functions above, unchanged in content (pure `git mv`), resolvable from the path `<root>/invz_common`. Later tasks call them exactly as `invz_projected` does today (`invz_single_ion(ion, T, B, opts)`, `invz_chi0z(si, T, z, opts)`, `invz_twolevel(ion, T, Bx, opts)`, `invz_g(tl, z)`, `invz_lambdas(K, g, wts, beta, [1 2])`, `invz_sigma(tl, lam, K, g, beta)`, `[wn, wts, beta] = invz_matsubara(T, Ecut)`).

- [ ] **Step 1: Audit the blast radius before moving anything.** From repo root run and record the output:

```bash
grep -rn "addpath" invz_projected/ --include="*.m" | grep -v tests/
grep -rln "mfilename" invz_projected/getf.m invz_projected/invz_const.m invz_projected/invz_ion.m invz_projected/stevens_ops.m invz_projected/invz_cfrot.m invz_projected/invz_field_vec.m invz_projected/invz_single_ion.m invz_projected/invz_chi0z.m invz_projected/invz_check_transverse_mf.m invz_projected/invz_matsubara.m invz_projected/invz_is_gamma_equiv.m invz_projected/invz_twolevel.m invz_projected/invz_g.m invz_projected/invz_lambdas.m invz_projected/invz_sigma.m
```

Expected: no moved file builds `mfilename`-relative paths (caching lives in `invz_jq_modes`/drivers, which stay). If a moved file DOES use `mfilename` for a data/cache path, STOP and report (BLOCKED) — do not improvise a fix.

- [ ] **Step 2: `git mv` the 15 files** into `invz_common/` (create the folder). One command per file or a loop; verify with `git status` that every move is tracked as a rename.

- [ ] **Step 3: Update every `setupOnce` in `invz_projected/tests/*.m`** to also add `invz_common` to the path (one added line per file, matching the boilerplate already there). Check drivers: `grep -rn "addpath" invz_projected/*.m` — if any driver addpaths explicitly (rather than assuming the caller's path), extend it the same way.

- [ ] **Step 4: Write `invz_common/README.md`** (~15 lines): what lives here (branch-agnostic single-ion + two-level/scalar-Σ engine shared by `invz_projected/` and `invz_tensor/`), the rule that its tests remain in `invz_projected/tests` (documented single suite of record for the common layer), and that callers must have `invz_common` on the path.

- [ ] **Step 5: Run the fast suite** (`invz_projected/tests`). Expected: pass/fail counts identical to the Task-0 baseline log (`.superpowers/sdd/baseline_fast_5f4ff92.log`) — 0 failed, same passed count, same incomplete count.

- [ ] **Step 6: Commit** — `refactor(invz): factor branch-shared single-ion + scalar-Sigma engine into invz_common/`.

---

### Task 2: P0 preflight — anchors fixture, exploratory ODD-block measurements, baseline freeze

**Files:**
- Create: `docs/ODD-LOG.md`, `invz_tensor/tests/invzt_anchors.m`, `invz_tensor/tests/exploratory/explore_odd_blocks.m`, `invz_tensor/tests/exploratory/explore_chiperp.m`
- Test: none (read-only phase — scripts + log; nothing enters the fast suite yet)

**Interfaces:**
- Produces: `A = invzt_anchors()` returning a struct with AT LEAST these fields, each measured on this tree and written as full-precision literals: `A.chiperp_1p53K_0T` (2×2, meV⁻¹, from `invz_chi0z` z=0 (a,b) block, `hyp=true`, `transverse_mf='legacy_x'`, elastic true), `A.chiperp_ab_asym_1p53K` (max abs antisymmetric part), `A.chiperp_0p31K_Bx` (struct: `Bx = [0 1 2 3 4 5 6]` T and the corresponding `chi_aa` values), `A.odd_generic_q_max` (max element magnitude in meV of the c-a block at q = [0.31 0.17 0.09], dpRng 30), `A.odd_onaxis_smallq` (struct with fields `.q = [1e-1 1e-2 1e-3]` (along [q 0 0]) and `.maxca` — max |J^{ca}| element in meV per q, dpRng 30), `A.jcc0` / `A.jaa0` (info values from `invz_jq_modes` at dpRng 30 for cross-reference).
- Produces: `docs/ODD-LOG.md` §P0 with: full-tensor availability statement (`MF_dipole` returns `[3,3,4,4]` — index layout `dip(mu,nu,s,s')`, units Å⁻³ pre-gfac, cc extraction at `invz_jq_modes` line ~66), the χ⊥ matrix and its interpretation vs the 16–17 meV⁻¹ DS2023 band, the on-axis decay + tilted-ray direction-dependence numbers with the `4*pi*gfac/ion.Vc` macroscopic scale for comparison, generic-q ODD magnitudes vs `Jcc0`, the cache-key + pre-parfor audit note, and the frozen baseline (suite outputs + git hash).

- [ ] **Step 1: Write `explore_chiperp.m`** (a plain script, `run`-able from repo root with the three addpaths of the boilerplate): build `ion = invz_ion()`; at (T, B) = (1.53 K, [0 0 0]) with `hyp=true` compute `si = invz_single_ion(ion, T, [0 0 0], struct('hyp', true))`, `c0 = invz_chi0z(si, T, 0, struct('elastic', true))`; print the (1:2,1:2) block, its symmetric/antisymmetric split, and the elastic-off variant (`'elastic', false`) for the elastic-share note. Sweep `Bx = 0:1:6` T at 0.31 K printing `chi_aa`, `chi_bb`. Expected: `chi_aa = chi_bb` at Bx = 0 to ~1e-10; magnitude in the 16–17 meV⁻¹ band (if it lands near 11, suspect truncation; if off by ×2-scale, suspect a convention slip — STOP and report per plan §0).

- [ ] **Step 2: Write `explore_odd_blocks.m`**: using `MF_dipole` directly (dpRng 30, `geom` reuse across q): (i) on-axis rays q = [q 0 0] and [0 0 q], q ∈ {1e-1, 1e-2, 1e-3} — print max |−gfac·dip(3,1,:,:)| and |−gfac·dip(3,2,:,:)| (meV), confirm decay toward 0 and record the rate; (ii) tilted ray q = q·[1 0 1]/√2, same q sequence — record the NON-decaying direction-dependent limit and compare its scale to `4*pi*C.gfac/ion.Vc` (≈ 3.7 μeV; element scale up to ≈ 1.8 μeV); (iii) generic q = [0.31 0.17 0.09] — record block element magnitudes vs `Jcc0 = 6.421e-3` meV; (iv) smallest-shell values on the standard 16³ grid (`q_min` along a* and c*) for the A0 dpRng-sensitivity note.

- [ ] **Step 3: Run both scripts**, then write `invzt_anchors.m` as a function returning the measured struct (full-precision `format long g` literals, each field commented with its provenance: script, date, git hash, dpRng).

- [ ] **Step 4: Baseline freeze.** Run the fast suite AND the slow suite (`INVZ_SLOW=1`) on this tree; append the result tables and the git hash to `docs/ODD-LOG.md` §P0 as the frozen baseline (Σc, Tc(0), Hc(0.31 K), soft-mode, timings — copy the console lines the slow tests print). Record the cache/parfor audit: `invz_jq_modes` key scheme (`jq4_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat`, pkey schema v4) and the statement that ODD/tensor caches get the separate `jqt` namespace in `invz_tensor/cache/`.

- [ ] **Step 5: Commit** — `docs(invzt): P0 preflight — anchors fixture, ODD-block exploration, frozen baseline (ODD-LOG)`.

---

### Task 3: `invzt_jq_tensor.m` — cached [12,12,nq] coupling tensor + DS2023 geometry guards

**Files:**
- Create: `invz_tensor/invzt_jq_tensor.m`, `invz_tensor/tests/test_invzt_jq_tensor.m`
- Reference (read, do not modify): `invz_projected/invz_jq_modes.m` (mirror its structure: geom priming call, per-q loop, Γ info block, cache pattern)

**Interfaces:**
- Consumes: `MF_dipole`, `exchange`, `invz_const`, `invz_is_gamma_equiv` (all on path).
- Produces: `[Jt, info] = invzt_jq_tensor(ion, qvec, opts)`. `Jt` `[12,12,nq]` complex Hermitian per page, sublattice-major index `i = 3*(s-1)+mu`. `opts`: `dpRng` (default 30), `cache` (default true), `parts` (default `'full'` = dipole + exchange + Γ-Lorentz; `'dipole'` = dipole only, no exchange, no Lorentz — for lattice-sum tests). `info`: `Jcc0`, `Jaa0`, `Jcc0_dipole`, `Jaa0_dipole` (uniform-mode projections, must equal `invz_jq_modes`'s to 1e-12 rel at demag = 0), `lorz`, `dpRng`, `geomD`, `geomX`. Later tasks slice blocks as `Jt(3*(s-1)+(1:3), 3*(sp-1)+(1:3), iq)`.

- [ ] **Step 1: Write the failing tests** — `test_invzt_jq_tensor.m` with the boilerplate `setupOnce` and these test functions (complete code; `ion = invz_ion()` in each):

```matlab
function test_shape_hermiticity_conjsym(testCase)
ion = invz_ion();
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09; -0.25 0 0; 0 0 -0.25; -0.31 -0.17 -0.09];
Jt = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
verifyEqual(testCase, size(Jt), [12 12 6]);
for iq = 1:6
    verifyLessThan(testCase, norm(Jt(:,:,iq) - Jt(:,:,iq)', 'fro'), 1e-14);
end
for iq = 1:3   % J(-q) = conj(J(q)): real-space couplings are real
    verifyLessThan(testCase, norm(Jt(:,:,iq+3) - conj(Jt(:,:,iq)), 'fro'), 1e-12);
end
end

function test_cc_block_parity_with_scalar_branch(testCase)
% The cc sub-matrix (mu = nu = 3) must reproduce invz_jq_modes' eigenvalues.
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];
[Jt, infoT] = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
[Jnu, infoS] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
cc = 3:3:12;
for iq = 1:3
    Jcc = (Jt(cc, cc, iq) + Jt(cc, cc, iq)')/2;
    verifyEqual(testCase, sort(real(eig(Jcc))).', Jnu(iq,:), 'AbsTol', 1e-12);
end
verifyEqual(testCase, infoT.Jcc0, infoS.Jcc0, 'RelTol', 1e-12);
verifyEqual(testCase, infoT.Jaa0, infoS.Jaa0, 'RelTol', 1e-12);
end

function test_gamma_lorentz_diagonal_only(testCase)
% At Gamma the Lorentz term sits on Cartesian-diagonal entries of every
% sublattice pair; off-diagonal (ODD and ab) entries carry NO shape term.
ion = invz_ion();
C = invz_const();
lorz = 4*pi/(3*ion.Vc)*C.gfac;
Jt  = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false));
Jtd = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false, 'parts', 'dipole'));
D = Jt - Jtd;   % full - dipole-only = exchange + Lorentz, both Cartesian-diagonal
for s = 1:4, for sp = 1:4
    blk = D(3*(s-1)+(1:3), 3*(sp-1)+(1:3));
    verifyLessThan(testCase, max(abs(blk(~eye(3)))), 1e-15);
    verifyGreaterThan(testCase, real(blk(3,3)), lorz - 1e-12);  % contains lorz (+ exchange)
end, end
end

function test_ds2023_geometry_sums(testCase)
% DS2023 Suppl. Table I (a = 5.175 Ang): pure-geometry real-space sums over the
% dpRng sphere, gfac-free — THE unit guard. Sum over all sites j (all 4
% sublattices) of the squared dipole tensor factors for a fixed central ion:
%   sum_j (3 x z / r^5)^2           = 36.73 / a^6   (B_xz,xz; = B_yz,yz)
%   sum_j (3 z z / r^5 - 1/r^3)^2   = 17.93 / a^6   (B_zz,zz)
ion = invz_ion();
[~, ~, geom] = MF_dipole([0 0 0], 30, ion.a, ion.tau);
a = 5.175;
[Sxz, Syz, Szz] = deal(0);
for sp = 1:4   % geom stores lower-triangular pairs; row 1 = central sublattice 1
    if sp >= 1
        nt = max(sp, 1); mt = min(sp, 1);
        Tf = geom.Tf{nt, mt};        % [nsites, 3, 3]: 3 r_n r_m / r^5 - delta_nm / r^3
        Sxz = Sxz + sum(Tf(:,3,1).^2);
        Syz = Syz + sum(Tf(:,3,2).^2);
        Szz = Szz + sum(Tf(:,3,3).^2);
    end
end
verifyEqual(testCase, Sxz, 36.73/a^6, 'RelTol', 0.01);
verifyEqual(testCase, Syz, 36.73/a^6, 'RelTol', 0.01);
verifyEqual(testCase, Szz, 17.93/a^6, 'RelTol', 0.01);
end

function test_onaxis_smallq_odd_decay(testCase)
% ODD blocks vanish on high-symmetry rays as q -> 0 (C2-about-c). ON-AXIS ONLY:
% tilted rays carry a direction-dependent macroscopic limit (plan A0 pitfall).
ion = invz_ion();
A = invzt_anchors();
q = [1e-1 0 0; 1e-2 0 0; 1e-3 0 0];
Jt = invzt_jq_tensor(ion, q, struct('dpRng', 30, 'cache', false));
ca = @(iq) Jt(3:3:12, 1:3:12, iq);   % c-row, a-column sublattice block
m = arrayfun(@(iq) max(abs(ca(iq)), [], 'all'), 1:3);
verifyEqual(testCase, m(:), A.odd_onaxis_smallq.maxca(:), 'RelTol', 1e-6);  % pinned P0 values
verifyLessThan(testCase, m(3), 1e-6 * 6.421e-3 + 5e-9);  % <= ~1e-6 * Jcc0 scale at q = 1e-3
end

function test_cache_roundtrip_selfverifying(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
cdir = fullfile(fileparts(mfilename('fullpath')), '..', 'cache');
J1 = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', true));
J2 = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', true));
verifyEqual(testCase, J2, J1);                       % bitwise round-trip
ion2 = ion;  ion2.J12 = ion.J12 * 1.05;              % key must miss on physics change
J3 = invzt_jq_tensor(ion2, q, struct('dpRng', 10, 'cache', true));
verifyFalse(testCase, isequal(J3, J1));
verifyTrue(testCase, ~isempty(dir(fullfile(cdir, 'jqt*.mat'))));
end
```

Note on `test_ds2023_geometry_sums`: `geom.Tf`/`geom.r` are stored lower-triangular (`{nt,mt}` for `mt <= nt`). The sums need ALL FOUR sublattice pairs seen from one central ion: pairs `{s,1}` for s ≥ 1 plus the self-sublattice `{1,1}`. The skeleton above is intentionally naive about that — as written it double-counts nothing but may MISS the `{s,1}` blocks for s > 1 if `geom` indexing differs; FIRST verify against the independently known cc anchor: `C.gfac * (Szz-type sum over the correct pair set)` must be consistent with the 6.821 μeV Γ benchmark's real-space content. Get the pair bookkeeping right (central ion on sublattice 1, sum over all sites of all 4 sublattices, self-site excluded — `MF_dipole` already excludes r² < 0.01), then assert the three DS2023 numbers. If after correct bookkeeping the sums are off by a uniform factor (2, 4, ...), STOP and report (BLOCKED) with the measured ratio — do not fold a fudge factor in.

- [ ] **Step 2: Run the test file to verify it fails** (function `invzt_jq_tensor` undefined).

- [ ] **Step 3: Implement `invzt_jq_tensor.m`** mirroring `invz_jq_modes.m`: demag guard first; geom priming call at q = [0 0 0] (capture `dip0` for the info block); per-q loop assembling all 9 Cartesian blocks (not just (3,3)): `blk = -C.gfac*dip(:,:,s,sp) + sign(ion.J12)*ex(:,:,s,sp)` for `parts='full'`, dipole-only for `'dipole'`; Γ-equivalent q get `+ lorz*eye(3)` per pair (only for `'full'`); Hermitize per page `Jt(:,:,iq) = (Jt(:,:,iq) + Jt(:,:,iq)')/2`; info block computed from `dip0` exactly as `invz_jq_modes` does (uniform-mode v = ones(4,1)/2 projections of the cc and aa blocks + `4*ion.J12`); cache per the Global Constraints key format (store `Jt`, `info`, `pkey`, `qvec`; loader isequal-checks `pkey` AND `qvec`). `parts='dipole'` results are NOT cached (they exist for tests).

- [ ] **Step 4: Run the tests to verify they pass**, then run BOTH fast suites.

- [ ] **Step 5 (INVZ_SLOW-gated, add to the same test file): Parseval cross-check** — BZ-average vs real-space sum with the SAME dpRng on both sides:

```matlab
function test_parseval_odd_vs_realspace_slow(testCase)
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();  C = invz_const();
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];   % uniform mesh incl. Gamma
Jt = invzt_jq_tensor(ion, qvec, struct('dpRng', 20, 'cache', true, 'parts', 'dipole'));
ca = Jt(3:3:12, 1:3:12, :);                                        % J^{ca} blocks, all q
lhs = mean(sum(abs(ca(1, :, :)).^2, 2), 3);                        % central sublattice row
[~, ~, geom] = MF_dipole([0 0 0], 20, ion.a, ion.tau);
rhs = 0;   % same pair bookkeeping as test_ds2023_geometry_sums, times gfac^2
% ... sum_j (C.gfac * Tf_xz)^2 over all sublattice pairs seen from sublattice 1
verifyEqual(testCase, lhs, rhs, 'RelTol', 0.01);
end
```

(The `rhs` assembly reuses the corrected pair bookkeeping from Step 1; the 1% tolerance absorbs superlattice folding at n = 8, which is r⁻⁶-suppressed.)

- [ ] **Step 6: Run the slow-gated test once** (`INVZ_SLOW=1`, this file only), record the Parseval residual in `docs/ODD-LOG.md`.

- [ ] **Step 7: Commit** — `feat(invzt): cached [12,12,nq] coupling tensor with DS2023 geometry + Parseval unit guards`.

---

### Task 4: `invzt_chi_rpa.m` + `invzt_gcc_lattice.m` — tensor RPA with parity gates

**Files:**
- Create: `invz_tensor/invzt_chi_rpa.m`, `invz_tensor/invzt_gcc_lattice.m`, `invz_tensor/tests/test_invzt_chi_rpa.m`
- Reference (read, do not modify): `invz_projected/invz_chi_tensor_ref.m` (the Σ=0 uniform-mode parity target), `invz_projected/invz_jq_modes.m`

**Interfaces:**
- Consumes: `Jt` `[12,12,nq]` from `invzt_jq_tensor` (Task 3); `chi0` `[3,3]` or `[3,3,nz]` from `invz_chi0z`.
- Produces:
  - `X = invzt_chi_rpa(chi0, Jt)`: `chi0` `[3,3]` (ONE frequency) → `X` `[12,12,nq]` with `X(:,:,iq) = (eye(12) - X0*Jt(:,:,iq)) \ X0`, `X0 = kron(eye(4), chi0)`. Implemented page-wise (`pagemtimes` + `pagemldivide`), no per-q MATLAB loop.
  - `[Gcc, diag4] = invzt_gcc_lattice(chi0, Jt)`: `chi0` `[3,3,nz]` → `Gcc` `[nz,1]` the site-diagonal cc lattice response averaged over q AND sublattices, `diag4` `[4,nz]` per-sublattice (for symmetry asserts). Sign: returns χ-convention (positive-definite-like at iωn); callers form `G = -Gcc`.

- [ ] **Step 1: Write the failing tests** (complete code):

```matlab
function test_uniform_mode_parity_with_tensor_ref(testCase)
% Sigma = 0 parity: at Gamma with the sublattice-uniform projection, the 12x12
% RPA must reproduce invz_chi_tensor_ref's single-site 3x3 with J = diag(Jaa0,Jaa0,Jcc0).
% Construct the parity limit EXACTLY: replace every 3x3 block of Jt(Gamma) by
% blockmean = (1/4) sum_{s,s'} block_{ss'} projected to diag(Jaa0,Jaa0,Jcc0)/4 —
% i.e. build Jd = kron(ones(4)/4, diag([Jaa0, Jaa0, Jcc0])) directly, so the
% uniform mode carries exactly the reference couplings.
ion = invz_ion();  T = 0.31;
[~, info] = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false));
si = invz_single_ion(ion, T, [4 0 0], struct('hyp', true, 'order', true, ...
    'J0z', ion.J0eff, 'Jxx0', ion.Jxx0));
w = (0.02:0.02:0.6).';  z = w + 1i*5e-3;
c0 = invz_chi0z(si, T, z, struct('elastic', false));
Jd = kron(ones(4)/4, diag([ion.Jxx0, ion.Jxx0, ion.J0eff]));
R = invz_chi_tensor_ref(ion, T, [4 0 0], w, struct('eta', 5e-3));
u = kron(ones(4,1)/2, eye(3));            % sublattice-uniform embedding, u'*u = I3
chi_ten = zeros(numel(w), 1);
for k = 1:numel(w)
    X = invzt_chi_rpa(c0(:,:,k), Jd);     % nq = 1 page
    Xu = u' * X * u;                      % uniform-mode 3x3
    chi_ten(k) = imag(Xu(3,3));
end
verifyEqual(testCase, chi_ten, R.chi_ten, 'RelTol', 1e-8, 'AbsTol', 1e-12);
end

function test_decoupled_cc_parity_with_scalar_branches(testCase)
% With Cartesian-off-diagonal chi0 entries zeroed AND ODD/ab blocks of Jt zeroed,
% the cc sector closes: q-resolved cc response == scalar RPA over the Jnu branches.
ion = invz_ion();  T = 1.6;
q = [0.25 0 0; 0.31 0.17 0.09];
[Jt, ~] = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
[Jnu, ~] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
c0 = invz_chi0z(si, T, 1i*0.05, struct('elastic', true));
c0d = diag(diag(c0(:,:,1)));
Jz = zero_offdiag_blocks(Jt);             % local helper: zero mu ~= nu entries of every block
X = invzt_chi_rpa(c0d, Jz);
x0cc = c0d(3,3);
for iq = 1:2
    Xcc = X(3:3:12, 3:3:12, iq);
    got = sort(real(eig((Xcc + Xcc')/2)));
    want = sort(x0cc ./ (1 - Jnu(iq,:).' * x0cc));   % scalar RPA per branch
    verifyEqual(testCase, got, want, 'RelTol', 1e-10);
end
end

function test_gcc_lattice_consistency(testCase)
% invzt_gcc_lattice must equal the brute-force per-q, per-frequency average.
ion = invz_ion();  T = 1.6;
q = [0.25 0 0; 0.31 0.17 0.09; 0.1 0.2 0.3];
Jt = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', false));
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
zn = 1i*[0.01 0.05 0.2];
c0 = invz_chi0z(si, T, zn, struct('elastic', true));
[Gcc, diag4] = invzt_gcc_lattice(c0, Jt);
for k = 1:3
    X = invzt_chi_rpa(c0(:,:,k), Jt);
    brute = zeros(4,1);
    for s = 1:4, brute(s) = mean(real(squeeze(X(3*(s-1)+3, 3*(s-1)+3, :)))); end
    verifyEqual(testCase, diag4(:,k), brute, 'RelTol', 1e-12);
    verifyEqual(testCase, Gcc(k), mean(brute), 'RelTol', 1e-12);
end
% per-sublattice equality (S4 site symmetry) — a physics assert, looser tol
verifyEqual(testCase, max(abs(diag4 - mean(diag4,1)), [], 'all')/max(abs(diag4(:))), 0, 'AbsTol', 1e-8);
end
```

(`zero_offdiag_blocks` is a ~6-line local function in the test file: loop s,s', keep `blk .* eye(3)`.)

- [ ] **Step 2: Run to verify failure**, **Step 3: implement both functions** (page-wise: `X0p = repmat(kron(eye(4), chi0), 1, 1, nq)` or implicit expansion; `A = eye(12) - pagemtimes(X0p, Jt)`; `X = pagemldivide(A, X0p)`; `invzt_gcc_lattice` loops frequencies calling the same page kernel, extracting `diag4`/`Gcc` without storing full X across frequencies), **Step 4: run tests + both fast suites**.

- [ ] **Step 5: Commit** — `feat(invzt): page-wise 12x12 tensor RPA with uniform-mode and decoupled-cc parity gates`.

---

### Task 5: `invzt_chiperp.m` + `invzt_odd_deltaJ.m` — E1/E4/E5 with the Schur-complement convention test

**Files:**
- Create: `invz_tensor/invzt_chiperp.m`, `invz_tensor/invzt_odd_deltaJ.m`, `invz_tensor/tests/test_invzt_odd_deltaJ.m`

**Interfaces:**
- Consumes: `Jt` from Task 3; `invz_single_ion`/`invz_chi0z` from common; `invzt_anchors()`.
- Produces:
  - `[Xp, info] = invzt_chiperp(ion, T, B, opts)`: `Xp` real symmetric 2×2 (meV⁻¹) — the symmetrized (a,b) block of the full electronuclear `invz_chi0z(si, T, 0, ...)` at the converged single-ion state; `opts`: `hyp` (true), `transverse_mf` ('legacy_x'), `elastic` (true). `info.asym` = max abs antisymmetric part (gyrotropic, discarded), `info.elastic_share`. DESIGN DECISION (record in header): χ⊥ is Van Vleck-dominated → computed once per (T,B), never inside a self-consistency loop. Must NOT route through `invz_twolevel` (regular across the CF gap; the `invz:degenerateDoublet` guard must not fire at Bx = 0).
  - `[dJ, d, info] = invzt_odd_deltaJ(Jt, Xp)`: `dJ` `[4,4,nq]` Hermitian — E1 contracted then E4 self-site-subtracted:
    E1: `dJpre(:,:,iq) = Vca*Xp(1,1)*Vca' + Vca*Xp(1,2)*Vcb' + Vcb*Xp(2,1)*Vca' + Vcb*Xp(2,2)*Vcb'` with `Vca(s,s') = Jt(3*(s-1)+3, 3*(s'-1)+1, iq)`, `Vcb(s,s') = Jt(3*(s-1)+3, 3*(s'-1)+2, iq)`;
    E4: `dJ(s,s,:) = dJpre(s,s,:) - mean_q(dJpre(s,s,:))` (diagonal only), off-diagonal passes through;
    E5: `d = mean over s of mean_q(dJpre(s,s,:))` — scalar, meV.
    `info`: `d_per_sublattice` [4,1] (assert equal to 1e-10 rel), `presub_min_eig` (assert ≥ −1e-12·max|dJpre|), `postsub_diag_bzavg` (assert 0 to machine), `dJ_max`.
    CALLER CONTRACT (document in header): the supplied `qvec` behind `Jt` must be a full uniform BZ mesh — E4/E5 are averages over the SUPPLIED grid. The −d uniform shift is applied by CALLERS to `Jcc0` (`Jcc0_odd = info.Jcc0 - d`), never by re-adding anything at q = 0 (plan §8 q=0 pitfall; cite E4/E5 + this plan in a comment to prevent future double counting — the subtracted on-site constant is a pure energy shift in the strict two-level limit; its physical residue is Tier 2's, no double counting).

- [ ] **Step 1: Write the failing tests** (complete code; the Schur test is THE decisive convention/double-counting gate):

```matlab
function test_chiperp_anchors_and_symmetry(testCase)
ion = invz_ion();  A = invzt_anchors();
[Xp, info] = invzt_chiperp(ion, 1.53, [0 0 0], struct());
verifyEqual(testCase, Xp, A.chiperp_1p53K_0T, 'RelTol', 1e-9);      % pinned P0 digits
verifyEqual(testCase, Xp(1,1), Xp(2,2), 'AbsTol', 1e-10*abs(Xp(1,1)));
verifyEqual(testCase, Xp, Xp.', 'AbsTol', 1e-15);
verifyLessThan(testCase, info.asym, 1e-8*abs(Xp(1,1)));
% smoothness along Bx at 0.31 K: no jumps > 25% between 1 T steps
x = arrayfun(@(b) getfield(invzt_chiperp(ion, 0.31, [b 0 0], struct()), {1,1}), 0:6); %#ok<GFLD>
verifyLessThan(testCase, max(abs(diff(x))./abs(x(1:end-1))), 0.25);
end

function test_deltaJ_structure(testCase)
% PSD before subtraction, zero diagonal BZ-average after, equal d on sublattices.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];
Jt = invzt_jq_tensor(ion, qvec, struct('dpRng', 15, 'cache', false));
Xp = invzt_chiperp(ion, 1.53, [0 0 0], struct());
[dJ, d, info] = invzt_odd_deltaJ(Jt, Xp);
verifyEqual(testCase, size(dJ), [4 4 n^3]);
verifyGreaterThan(testCase, info.presub_min_eig, -1e-12*max(1e-30, info.dJ_max));
verifyLessThan(testCase, abs(info.postsub_diag_bzavg), 1e-15*max(1e-30, info.dJ_max));
verifyEqual(testCase, info.d_per_sublattice, d*ones(4,1), 'RelTol', 1e-10);
verifyGreaterThan(testCase, d, 0);
% closure of E4/E5: BZ-avg of post-subtraction diagonal recovers -d ... it is 0
% by construction, so instead reconstruct: mean_q diag(dJpre) = d per sublattice.
verifyEqual(testCase, mean(real(dJ(1,1,:) + d), 3), d, 'RelTol', 1e-9);
% Xp = 0 -> dJ identically 0, d = 0
[dJ0, d0] = invzt_odd_deltaJ(Jt, zeros(2));
verifyEqual(testCase, max(abs(dJ0(:))), 0);
verifyEqual(testCase, d0, 0);
end

function test_schur_complement_equals_E1(testCase)
% THE convention gate (plan A0): eliminate the transverse block of the 12x12
% RPA by Schur complement. With the transverse-transverse lattice blocks
% REMOVED (set to zero), the effective cc coupling equals Jcc + deltaJ_pre(q)
% EXACTLY (E1 with the same Xp), at any q. Uses a synthetic diagonal chi0
% (chi_aa = chi_bb = xperp, chi_cc = xcc, cross terms 0) so the elimination is
% algebraically clean — this pins index conventions and factor-of-2/dagger slips.
ion = invz_ion();
q = [0.31 0.17 0.09];
Jt = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
xperp = 11.05;  xcc = 40.0;                    % meV^-1, 3-state-scale magnitudes
% zero the transverse-transverse blocks (aa, ab, ba, bb) of every pair,
% keep cc and the ODD ca/cb blocks:
Jz = Jt;
for s = 1:4, for sp = 1:4
    r = 3*(s-1); c = 3*(sp-1);
    Jz(r+(1:2), c+(1:2), 1) = 0;
end, end
X0 = kron(eye(4), diag([xperp, xperp, xcc]));
X = (eye(12) - X0*Jz(:,:,1)) \ X0;             % full tensor RPA at this q
% cc sublattice block of the tensor response:
Xcc = X(3:3:12, 3:3:12);
% E1 route: scalar-RPA with effective coupling Jcc + dJpre
Vca = Jz(3:3:12, 1:3:12, 1);  Vcb = Jz(3:3:12, 2:3:12, 1);
Jcc = Jz(3:3:12, 3:3:12, 1);
dJpre = xperp*(Vca*Vca') + xperp*(Vcb*Vcb');   % E1 with Xp = xperp*I2
Xcc_e1 = xcc * ((eye(4) - (Jcc + dJpre)*xcc) \ eye(4));
verifyEqual(testCase, Xcc, Xcc_e1, 'RelTol', 1e-10);
end

function test_full_tensor_vs_E1_correction_logged(testCase)
% With the REAL transverse-transverse blocks kept, Schur gives
% J_eff = Jcc + Jca*(Xp^-1 - Jaa_block)^-1*Jac + ... — E1 is its Jaa->0 limit.
% REPORT (not gate) the relative size of the difference at a generic q: this is
% the transverse-channel RPA enhancement that A1 keeps and Tier 1 drops.
ion = invz_ion();
q = [0.31 0.17 0.09];
Jt = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
Xp = invzt_chiperp(ion, 1.53, [0 0 0], struct());
[dJ, d] = invzt_odd_deltaJ(Jt, Xp);   %#ok<ASGLU>
% exact Schur at this q (2x2 transverse per sublattice pair kept):
% ... build P = transverse 8x8 block matrix, S = Jcc_block + Vc,perp*(inv(kron(I4,Xp)) - P)^-1*Vperp,c
% log norm(S - (Jcc + dJpre))/norm(dJpre) to docs/ODD-LOG.md via fprintf; assert only finiteness.
verifyTrue(testCase, true);
end
```

(The fourth test's body is a REPORT: implement the exact 8×8 Schur elimination as sketched, `fprintf` the ratio, assert it is finite and < 1 — the measured number goes to ODD-LOG. The E1-vs-Schur gap at χ⊥·Jaa0 ~ 0.05 scale is expected physics, not a bug.)

- [ ] **Step 2: Run to verify failure**, **Step 3: implement both functions** per the interface contract, **Step 4: run tests + both fast suites**. Record in `docs/ODD-LOG.md`: `d` at (1.53 K, 0 T) on the 6³/dpRng-15 test grid AND on a production 16³/dpRng-30 grid (expected band 0.3–0.5 μeV — report the number either way; an order-of-magnitude miss with Task 3's geometry sums green means the χ⊥ contraction is at fault — STOP and report), `info.asym`, `info.elastic_share`, and the Schur-vs-E1 correction ratio.

- [ ] **Step 5: Commit** — `feat(invzt): E1/E4/E5 mediated-coupling diagnostics with exact Schur-complement convention gate`.

---

### Task 6: `invzt_chi0_split.m` — transition-mask decomposition χ0 = χ_dom + χ_rest

**Files:**
- Create: `invz_tensor/invzt_chi0_split.m`, `invz_tensor/tests/test_invzt_chi0_split.m`
- Reference (read, do not modify): `invz_common/invz_chi0z.m` (the mask reuses its exact algebra)

**Interfaces:**
- Consumes: `si` from `invz_single_ion`.
- Produces: `[cdom, crest, mspec] = invzt_chi0_split(si, T, z, opts)`; `cdom + crest == invz_chi0z(si, T, z, opts_pass)` EXACTLY (enforced by construction: `cdom` = masked re-run of the chi0z algebra, `crest = chi_full - cdom`); `opts`: `Esplit` (default 0.4653 meV = half the 0.9306 meV CF gap), `elastic` (true; the elastic term belongs to `cdom` — both states of an elastic pair are in the ground group whenever the pair contributes), plus passthrough `degtol`/`ztol`. Mask: transition (a,b) is DOMINANT iff BOTH `E_a - E0 < Esplit` AND `E_b - E0 < Esplit`. `mspec`: `ndom` (count of ground-group states — 16 with hyp, 2 without), `fdom_cc0` (fraction `cdom_cc(z=0)/chi_full_cc(z=0)` — report), `fdom_perp0` (same for the aa element — expected SMALL: the transverse response is Van Vleck, i.e. lives in `crest`).

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_split_exact_and_groupsize(testCase)
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
z = [0, 1i*0.05, 1i*0.5, 0.1 + 1i*5e-3];
full = invz_chi0z(si, T, z, struct('elastic', true));
[cdom, crest, mspec] = invzt_chi0_split(si, T, z, struct());
verifyEqual(testCase, cdom + crest, full, 'AbsTol', 1e-14*max(abs(full(:))));
verifyEqual(testCase, mspec.ndom, 16);                      % 2 electronic x 8 nuclear
si0 = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', false));
[~, ~, m0] = invzt_chi0_split(si0, T, 0, struct());
verifyEqual(testCase, m0.ndom, 2);
end

function test_dominant_sector_is_longitudinal(testCase)
% cc weight lives in the ground group; transverse (Van Vleck) weight in crest.
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
[cdom, crest, mspec] = invzt_chi0_split(si, T, 0, struct());
verifyGreaterThan(testCase, mspec.fdom_cc0, 0.90);          % doublet dominates chi_cc
verifyLessThan(testCase, mspec.fdom_perp0, 0.10);           % chi_perp is inter-doublet
verifyGreaterThan(testCase, real(crest(1,1)), 0);           % Van Vleck part positive
end

function test_esplit_insensitivity_band(testCase)
% Any Esplit in (0.05, 0.8) meV selects the same 16-state group (the hyperfine
% ladder spans ~0.01 meV; the first CF level sits at 0.9306 meV).
ion = invz_ion();  T = 1.6;
si = invz_single_ion(ion, T, [0.1 0 0], struct('hyp', true));
[c1] = invzt_chi0_split(si, T, 1i*0.05, struct('Esplit', 0.1));
[c2] = invzt_chi0_split(si, T, 1i*0.05, struct('Esplit', 0.7));
verifyEqual(testCase, c1, c2, 'AbsTol', 1e-15*max(abs(c1(:))));
end
```

- [ ] **Step 2: Run to verify failure**, **Step 3: implement** (copy the `invz_chi0z` loop verbatim into a masked variant: weight `w(a,b)` zeroed for non-dominant pairs, elastic term restricted to dominant pairs; `crest = invz_chi0z(...) - cdom`; header documents that the split is exact by construction and cites A1), **Step 4: run tests + both fast suites**.

- [ ] **Step 5: Commit** — `feat(invzt): transition-mask chi0 decomposition for the A1 dominant-sector bridge`.

---

### Task 7: `invzt_solve_point.m` — the A1 bridge point solver

**Files:**
- Create: `invz_tensor/invzt_solve_point.m`, `invz_tensor/tests/test_invzt_solve_point.m`
- Reference (read, do not modify): `invz_projected/invz_solve_point.m` (loop skeleton, option names, pt-field contract), `jensen_1z_framework.html` §5 (K definition), §7 (sum rule)

**Interfaces:**
- Consumes: `Jt` + `info` (Task 3), `invzt_gcc_lattice` (Task 4), `invzt_chi0_split` (Task 6), common engine (`invz_matsubara`, `invz_single_ion`, `invz_chi0z`, `invz_twolevel`, `invz_g`, `invz_lambdas`, `invz_sigma`).
- Produces: `pt = invzt_solve_point(ion, T, B, Jt, opts)` with `opts`: `Ecut` 40, `hyp` true, `transverse_mf` 'legacy_x', `mix_outer` 0.7, `tol_outer` 1e-8, `max_outer` 60, `Sigma_seed`, `odd` (true; false ZEROES every Cartesian-off-diagonal entry of every 3×3 block of a COPY of `Jt` — `ca/cb/ab` — leaving `aa/bb/cc`), `chi_rest` (true; false drops `crest` — quantifies the weak-transition RPA effect), `JtGamma` (REQUIRED: `[12,12]` tensor at Γ from `invzt_jq_tensor(ion, [0 0 0], ...)` — BZ grids exclude Γ), `Esplit` passthrough.
  ALGORITHM (per iteration, all frequencies vectorized):
  1. once: `[wn, wts, beta] = invz_matsubara(T, Ecut)`; `si`; `[cdom, crest] = invzt_chi0_split(si, T, 1i*wn, ...)`; `tl = invz_twolevel(ion, T, B, ...)`; `g = real(invz_g(tl, 1i*wn))`; `G0d = -real(squeeze(cdom(3,3,:)))` (dominant-sector cc propagator, the object Σ renormalizes).
  2. loop: `ctil(:,:,n) = cdom(:,:,n)/(1 + Sigma(n)) + crest(:,:,n)` → `[Gcc, diag4] = invzt_gcc_lattice(ctil, Jt_eff)` → `Gloc = -Gcc` → K from the scalar-medium Dyson form used by the projected branch (framework §5): `K = (G0til - Gloc.*(1 + Sigma)) ./ (G0til .* Gloc)` where `G0til = G0d + Grest_cc` — CAREFUL: derive the bookkeeping so that in the decoupled no-ODD limit it reduces EXACTLY to the projected `invz_emt_scalar` relation `Gloc = G0/(1 + Sigma + K*G0)`; state the chosen G0til convention in the header and verify it in the regression test below → `lam = invz_lambdas(K, g, wts, beta, [1 2])` → `sg = invz_sigma(tl, lam, K, g, beta)` → mix, converge on `max|dSigma|`.
  3. outputs mirror `invz_solve_point`: `Sigma0`, `Sigma`, `alpha`, `lambda`, `K`, `G` (= Gloc), `tl`, `si`, `chi0cc0`, `crit`, `sumrule_rel`, `converged`, `outer_iters`, plus `diag4_spread` (max per-sublattice deviation) and `pt.odd`.
  CRIT (uniform-mode instability measure): `M = eye(12) - kron(eye(4), ctil0) * JtGamma_eff` with `ctil0` the static (z = 0, elastic true) renormalized tensor; `pt.crit = min(real(eig(M)))`; assert `max|imag(eig)| < 1e-8`. Sign contract: `crit > 0` in the PM phase, → 0 at the boundary (verified in the regression test against the projected `pt.crit` sign).
  SUM RULE: `sumrule_rel = |sum(wts.*G)/beta + si.JzJz_fluct| / si.JzJz_fluct` exactly as projected (equal-time cc sum rule on the LOCAL lattice propagator; report-quality diagnostic).

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_no_odd_regression_vs_projected(testCase)
% Decoupled limit: odd = false AND chi_rest = false AND Cartesian-diagonal chi0
% is NOT enforceable from outside — so the regression is on OBSERVABLES:
% Sigma0 and crit sign must track invz_solve_point closely (the residual
% difference = cross-Cartesian chi0 elements propagating in aa/bb channels,
% expected small in the PM phase at moderate Bx). Gate: |dSigma0| <= 2e-3 and
% same crit sign; report both numbers.
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];  % drop Gamma
[Jt, info] = invzt_jq_tensor(ion, qvec, struct('dpRng', 15, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false));
pt = invzt_solve_point(ion, T, B, Jt, struct('odd', false, 'chi_rest', false, 'JtGamma', JtG));
verifyTrue(testCase, pt.converged);
[Jnu, ~] = invz_jq_modes(ion, qvec, struct('dpRng', 15, 'cache', true));
ps = invz_solve_point(ion, T, B, Jnu(:), struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0));
verifyTrue(testCase, ps.converged);
verifyEqual(testCase, pt.Sigma0, ps.Sigma0, 'AbsTol', 2e-3);
verifyEqual(testCase, sign(pt.crit), sign(ps.crit));
fprintf('no-ODD bridge vs projected: dSigma0 = %.3e, crit_t = %.4f, crit_s = %.4f\n', ...
    pt.Sigma0 - ps.Sigma0, pt.crit, ps.crit);
end

function test_odd_on_converges_and_moves_crit(testCase)
% PM point with ODD on: converges, overhead sane, crit strictly SMALLER than
% odd-off at the same point (ODD adds fluctuation weight -> closer to critical).
ion = invz_ion();  T = 1.6;  B = [0.1 0 0];
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt = invzt_jq_tensor(ion, qvec, struct('dpRng', 15, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false));
t0 = tic;
p1 = invzt_solve_point(ion, T, B, Jt, struct('odd', true,  'JtGamma', JtG));
p0 = invzt_solve_point(ion, T, B, Jt, struct('odd', false, 'JtGamma', JtG));
verifyTrue(testCase, p1.converged && p0.converged);
verifyLessThan(testCase, p1.crit, p0.crit);
verifyLessThan(testCase, p1.sumrule_rel, 0.10);
verifyLessThan(testCase, p1.diag4_spread, 1e-6);
fprintf('odd on/off crit: %.5f / %.5f, Sigma0: %.4f / %.4f, t = %.1fs\n', ...
    p1.crit, p0.crit, p1.Sigma0, p0.Sigma0, toc(t0));
end

function test_chi_rest_toggle_reported(testCase)
% crest off vs on: report the Sigma0/crit shift (the weak-transition RPA
% content A1 adds over the strict two-level projection). Assert only finiteness
% and convergence.
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt = invzt_jq_tensor(ion, qvec, struct('dpRng', 10, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
pa = invzt_solve_point(ion, T, B, Jt, struct('chi_rest', true,  'JtGamma', JtG));
pb = invzt_solve_point(ion, T, B, Jt, struct('chi_rest', false, 'JtGamma', JtG));
verifyTrue(testCase, pa.converged && pb.converged);
fprintf('chi_rest on/off: crit %.5f / %.5f, Sigma0 %.4f / %.4f\n', pa.crit, pb.crit, pa.Sigma0, pb.Sigma0);
end
```

- [ ] **Step 2: Run to verify failure**, **Step 3: implement** per the algorithm block (derive the `G0til`/K bookkeeping FIRST on paper in the function header comment — the no-ODD reduction to `invz_emt_scalar`'s Dyson relation is the correctness anchor; if the regression gate cannot be met, STOP and report the measured `dSigma0` rather than loosening the gate), **Step 4: run tests + both fast suites** (watch runtime: the two solver tests should stay < 60 s combined at n = 8/dpRng 15 with warm cache; if not, drop to n = 6 and note it).

- [ ] **Step 5: Commit** — `feat(invzt): A1 dominant-sector bridge point solver with no-ODD regression gate`.

---

### Task 8: Critical finders + zero-field physics measurements (T1.5 analog)

**Files:**
- Create: `invz_tensor/invzt_critical.m`, `invz_tensor/invzt_critical_T0field.m`, `invz_tensor/tests/test_invzt_critical.m`
- Reference (read, do not modify): `invz_projected/invz_critical.m`, `invz_projected/invz_critical_T0field.m`, `invz_projected/tests/test_invz_sigma_crit.m` (the published Σc/Tc slow benchmarks — mirror their grids/tolerances)

**Interfaces:**
- Consumes: Tasks 3, 5, 7 outputs.
- Produces:
  - `[Tc, out] = invzt_critical_T0field(ion, opts)` — the Tier-1 CLOSED-FORM zero-field route: per grid n ∈ `opts.grids` (default {12, 24} + Richardson): build `Jt` (dpRng 30, cached) on the n³ mesh (Γ excluded exactly as the projected route does), `Xp(T) = invzt_chiperp(ion, T, [0 0 0])`, `[dJ, d] = invzt_odd_deltaJ(Jt, Xp)`, modes `Jnu_odd(q) = sort(eig(Jcc4(q) + dJ(q)))`, uniform `J0(T) = info.Jcc0 - d(T)`, `Sigma_c(0; T) = mean_{q,nu} Jnu_odd/(J0 - Jnu_odd)` (E2), then scalar root-find on `J0(T)*chi0cc(0; T) - 1 - Sigma_c(0; T) = 0` over T (χ⊥ nearly flat — secant from the no-ODD Tc, ≤ 6 iterations). `opts.odd = false` bypasses δJ/d entirely and MUST reproduce the projected route. `out`: per-grid Tc, Richardson value, `Sigma_c`, `d`, plus the three-way decomposition below.
  - DECOMPOSITION (T1.5, the key DS2023-mappable measurement) inside `out.split`: Tc computed three ways — (a) full: modes of `Jcc4 + dJ`, `J0 = Jcc0 - d`; (b) uniform-shift only: modes of `Jcc4`, `J0 = Jcc0 - d` (DS2023's MF mechanism); (c) q-structure only: modes of `Jcc4 + dJ + d*I4` per q, `J0 = Jcc0` (the fluctuation piece their MF misses). Report (a), (b), (c) and the closure defect `(a) - [(b) + (c) - baseline]`.
  - `[Bc, out] = invzt_critical(ion, T, Jt, JtGamma, Brange, opts)` — PM-side bisection on `pt.crit` from `invzt_solve_point` (Σ warm-start threaded between iterations like `invz_critical`; `Jt`/`JtGamma` precomputed by the caller ONCE outside any loop — P0.4 discipline).
- Tests: fast tests exercise structure on small grids; the published-benchmark reproduction and the ODD measurements are INVZ_SLOW-gated.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_t0field_no_odd_matches_projected_route(testCase)
% Fast structural gate on a small grid: odd = false must give the same Tc as
% the projected closed form evaluated on the same grid/dpRng (identical modes).
ion = invz_ion();
[TcT, outT] = invzt_critical_T0field(ion, struct('odd', false, 'grids', {{8}}, 'dpRng', 15));
TcS = invz_critical_T0field(ion, struct('grids', {{8}}, 'dpRng', 15));   % adapt to its true signature
verifyEqual(testCase, TcT, TcS, 'RelTol', 1e-9);
verifyEqual(testCase, outT.d, 0);
end

function test_t0field_odd_lowers_Tc(testCase)
ion = invz_ion();
[Tc0] = invzt_critical_T0field(ion, struct('odd', false, 'grids', {{8}}, 'dpRng', 15));
[Tc1, out1] = invzt_critical_T0field(ion, struct('odd', true, 'grids', {{8}}, 'dpRng', 15));
verifyLessThan(testCase, Tc1, Tc0);
verifyGreaterThan(testCase, out1.d, 0);
% split consistency: (b) and (c) each lower Tc; (a) is the largest suppression
verifyLessThan(testCase, out1.split.Tc_b, Tc0);
verifyLessThan(testCase, out1.split.Tc_c, Tc0);
verifyLessThan(testCase, Tc1, min(out1.split.Tc_b, out1.split.Tc_c) + 1e-6);
end

function test_published_benchmarks_odd_off_slow(testCase)
% ODD-off Richardson(12,24) must land on the published no-ODD anchors.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
[Tc, out] = invzt_critical_T0field(ion, struct('odd', false, 'grids', {{12, 24}}, 'dpRng', 30));
verifyEqual(testCase, Tc, 1.743, 'AbsTol', 0.006);
verifyEqual(testCase, out.Sigma_c_rich, 0.2980, 'AbsTol', 0.006);
end

function test_odd_measurements_slow(testCase)
% REPORT: the headline zero-field numbers. Assertions are reproducibility-only
% (rerun-stability), never magnitude — record every number in docs/ODD-LOG.md.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
[Tc1, out1] = invzt_critical_T0field(ion, struct('odd', true, 'grids', {{12, 24}}, 'dpRng', 30));
fprintf('ODD Tier1 closed form: Tc = %.4f K (dTc = %.4f), Sigma_c = %.4f, d = %.4g ueV\n', ...
    Tc1, 1.743 - Tc1, out1.Sigma_c_rich, out1.d*1e3);
fprintf('split: (a) %.4f  (b) %.4f  (c) %.4f K\n', Tc1, out1.split.Tc_b, out1.split.Tc_c);
verifyTrue(testCase, isfinite(Tc1) && Tc1 > 0.5 && Tc1 < 1.743);
end

function test_bc_finder_pm_side_slow(testCase)
% A1 boundary point: Bc(1.2 K) odd-off parity with the projected finder, then
% the odd-on shift (REPORT).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt = invzt_jq_tensor(ion, qvec, struct('dpRng', 30, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 30, 'cache', false));
[Bc0, ~] = invzt_critical(ion, 1.2, Jt, JtG, [2 6], struct('odd', false));
[Bc1, ~] = invzt_critical(ion, 1.2, Jt, JtG, [2 6], struct('odd', true));
fprintf('A1 Bc(1.2 K): odd off %.4f T, on %.4f T, shift %+.4f T\n', Bc0, Bc1, Bc1 - Bc0);
verifyTrue(testCase, isfinite(Bc0) && isfinite(Bc1));
verifyLessThan(testCase, Bc1, Bc0 + 1e-3);   % ODD cannot RAISE the PM boundary
end
```

(Adapt `invz_critical_T0field`'s actual signature when writing the parity test — read the source first; if its option names differ, mirror THEM, and note the mapping in the test comment.)

- [ ] **Step 2: Run fast tests to verify failure**, **Step 3: implement both finders**, **Step 4: fast suites green**, **Step 5: run the slow tests** (`INVZ_SLOW=1` on this file), **Step 6: append the full measurement block to `docs/ODD-LOG.md`**: ΔTc(0), ΔΣc(0), the (a)/(b)/(c) split + closure defect, d and its (T) dependence, the Tier-1-static vs A1-tensor comparison at Bx = 0.05 T proxy (run `invzt_solve_point` crit(T) root by hand for the log if cheap — otherwise defer the A1 zero-field number to Task 10's log entry and SAY SO), Bc(1.2 K) shift. State the T2-decision analog: A1 carries χ⊥(iωn) retardation automatically; the static-vs-retarded gap is |ΔTc(closed form) − ΔTc(A1)| — record it and whether it clears the 10 mK threshold of the original plan's T2.2.

- [ ] **Step 7: Commit** — `feat(invzt): critical finders + zero-field ODD measurements (Tier-1 closed form, MF/fluctuation split)`.

---

### Task 9: Tier 2 components — `invzt_odd_fieldvar.m` + `invzt_twolevel_avg.m`

**Files:**
- Create: `invz_tensor/invzt_odd_fieldvar.m`, `invz_tensor/invzt_twolevel_avg.m`, `invz_tensor/tests/test_invzt_tier2.m`
- Reference: `jensen_1z_framework.html` §7 (Matsubara tail handling), `invz_common/invz_twolevel.m`, `invz_common/invz_g.m`

**Interfaces:**
- Consumes: converged `pt` from `invzt_solve_point` (Task 7), `Jt`, common engine.
- Produces:
  - `[C, info] = invzt_odd_fieldvar(pt, ion, T, Jt, opts)` — E3: `C_ab = (1/(4*Nq)) * sum_q tr_subl[ Vac(q) * Scc(q) * Vcb(q) ]` with `Scc(q) [4,4] = (1/beta) * sum_n chi_cc_block(q, i wn)` from the converged bridge solution (rebuild `ctil` from `pt.Sigma`/`pt.si` exactly as the solver does; reuse `invz_matsubara` wts for the tail — the sum converges like ωn⁻²; handle the tail the same way the sum rule does). `Vac = Vca'` blocks from `Jt`. `C` real symmetric 2×2 PSD, meV². `info`: `heq_T` = `sqrt(C_aa)` converted to an equivalent transverse field in Tesla via `gL*muB` and `<Jx>`-scale (Dollberg Fig. 4 width h ≈ 0.4 T is the qualitative comparator — log only), `tail_share`.
  - `[tla, info] = invzt_twolevel_avg(ion, T, B, C, opts)` — Gauss–Hermite average (default `opts.ngh = 7` per axis; diagonalize `C = V*diag(s2)*V'` first, nodes in the rotated frame): at each node `(ha, hb)` the extra transverse field enters as a FIELD SHIFT `B_node = B + [ha, hb, 0]/(gL*muB)` with `gL*muB` from `invz_const` (`C.muB*5/4`) — this injects exactly `-ha*Ja - hb*Jb` through the existing Zeeman term, so NO single-ion code changes are needed; per node run `invz_twolevel` and `invz_g` on the node's params; average the RESPONSE `gbar(i wn) = sum_w g_node(i wn)` (response-averaging is the default; parameter-averaging sits behind `opts.avg = 'params'` for the one-shot comparison); then fit effective two-level params by matching (i) `gbar(0)`, (ii) `(1/beta)*sum_n gbar` (must equal 1 per node — assert it survives averaging to 1e-6), (iii) the leading ωn⁻² tail coefficient. Return `tla` with the SAME fields as `invz_twolevel` so `invz_g`/`invz_lambdas`/`invz_sigma` consume it unchanged; `info`: node params spread, `Delta_eff`, `M2_eff`, `n01_eff`, fit residuals. FLAG IN HEADER (verbatim requirement from the plan): this quenched-Gaussian dressing is the least rigorous step of the whole program.
  - Degenerate-doublet handling at B → 0 (T3.4): the off-origin GH nodes lift the degeneracy; the exact-origin node (only present for odd `ngh`) is evaluated in its `h → 0` limit — document which branch the code takes and confirm no `invz:degenerateDoublet` throw at `B = [0 0 0]`, `C > 0`.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_fieldvar_structure(testCase)
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt  = invzt_jq_tensor(ion, qvec, struct('dpRng', 10, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
pt = invzt_solve_point(ion, T, B, Jt, struct('JtGamma', JtG));
assumeTrue(testCase, pt.converged);
[C, info] = invzt_odd_fieldvar(pt, ion, T, Jt, struct());
verifyEqual(testCase, C, C.', 'AbsTol', 1e-15*max(abs(C(:))));
e = eig(C);  verifyGreaterThan(testCase, min(e), -1e-14*max(e));
verifyGreaterThan(testCase, max(e), 0);                    % nonzero variance
fprintf('field variance: C_aa = %.4g meV^2, heq = %.3f T (Dollberg ~0.4 T)\n', C(1,1), info.heq_T);
end

function test_fieldvar_aa_bb_smallB(testCase)
% C_aa = C_bb at (near-)zero transverse field.
ion = invz_ion();  T = 1.6;  B = [0.05 0 0];
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt  = invzt_jq_tensor(ion, qvec, struct('dpRng', 10, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
pt = invzt_solve_point(ion, T, B, Jt, struct('JtGamma', JtG));
assumeTrue(testCase, pt.converged);
C = invzt_odd_fieldvar(pt, ion, T, Jt, struct());
verifyEqual(testCase, C(1,1), C(2,2), 'RelTol', 5e-2);     % small-B, not exactly 0
end

function test_twolevel_avg_C0_bitwise(testCase)
ion = invz_ion();  T = 1.6;  B = [2 0 0];
tl  = invz_twolevel(ion, T, B(1), struct());
tla = invzt_twolevel_avg(ion, T, B, zeros(2), struct());
verifyEqual(testCase, tla, tl);                            % C = 0 short-circuits to invz_twolevel
end

function test_twolevel_avg_monotone_on_C_ray(testCase)
% Level repulsion: Delta_eff grows with ||C||; moment sum rule: M2_eff shrinks.
ion = invz_ion();  T = 1.6;  B = [2 0 0];
s = [1e-8 4e-8 1.6e-7];                                    % meV^2 scale ray
[D, M2] = deal(zeros(3,1));
for i = 1:3
    [~, info] = invzt_twolevel_avg(ion, T, B, s(i)*eye(2), struct());
    D(i) = info.Delta_eff;  M2(i) = info.M2_eff;
end
tl = invz_twolevel(ion, T, B(1), struct());
verifyTrue(testCase, all(diff(D) > 0) && D(1) >= tl.Delta - 1e-12);
verifyTrue(testCase, all(diff(M2) < 0) && M2(1) <= tl.M2 + 1e-12);
end

function test_twolevel_avg_zero_field_no_guard(testCase)
% T3.4: at B = 0 with C > 0 the GH nodes lift the doublet degeneracy; no throw.
ion = invz_ion();
tla = invzt_twolevel_avg(ion, 1.55, [0 0 0], 4e-8*eye(2), struct());
verifyTrue(testCase, isstruct(tla) && isfinite(tla.Delta) && tla.Delta > 0);
end

function test_avg_mode_comparison_reported(testCase)
% response-averaging (default) vs parameter-averaging: one-shot report.
ion = invz_ion();  T = 1.6;  B = [2 0 0];
[~, ia] = invzt_twolevel_avg(ion, T, B, 4e-8*eye(2), struct());
[~, ip] = invzt_twolevel_avg(ion, T, B, 4e-8*eye(2), struct('avg', 'params'));
fprintf('avg-mode Delta_eff: response %.6g vs params %.6g meV\n', ia.Delta_eff, ip.Delta_eff);
verifyTrue(testCase, isfinite(ia.Delta_eff) && isfinite(ip.Delta_eff));
end
```

(Check `invz_twolevel`'s actual field names — `tl.Delta`/`tl.M2`/`tl.n01` are assumed above; read the source and use ITS names everywhere, including in `tla`.)

- [ ] **Step 2: Run to verify failure**, **Step 3: implement both functions** (fieldvar: rebuild `ctil` from the pt struct, page-solve per frequency accumulating `Scc(q)` without storing `[12,12,nq,nz]`; twolevel_avg: C = 0 short-circuit FIRST — literal `tla = invz_twolevel(...); return`), **Step 4: run tests + both fast suites**. Log `C`, `heq_T`, the avg-mode comparison, and the C-ray monotonicity table to `docs/ODD-LOG.md`.

- [ ] **Step 5: Commit** — `feat(invzt): Tier-2 components — ODD field covariance (E3) and Gauss-Hermite-dressed doublet`.

---

### Task 10: Tier-2 outer self-consistency in `invzt_solve_point` + combined measurements

**Files:**
- Modify: `invz_tensor/invzt_solve_point.m` (add the Tier-2 outer loop behind `opts.odd_tier2`, default false)
- Create: `invz_tensor/tests/test_invzt_tier2_outer.m`

**Interfaces:**
- Consumes: Tasks 7 + 9.
- Produces: `opts.odd_tier2 = true` → after the Tier-1 bridge solve converges: `C = invzt_odd_fieldvar(...)` → `tla = invzt_twolevel_avg(ion, T, B, C, opts.tier2)` → re-solve the Σ loop with `tl = tla` (and `g` rebuilt from it) → iterate outer₂ until `max(abs(C - C_prev), [], 'all') < opts.tol_tier2 * max(abs(C(:)))` (default 1e-3), cap `opts.max_tier2 = 8`, mix 0.5 on C, non-convergence masked exactly like the EMT loop (pt.converged = false). New pt fields: `C`, `tier2_iters`, `tla`. `opts.odd_tier2 = false` path BYTE-IDENTICAL to Task 7 behavior (guard everything behind the flag).
- IR safety (T3.3, strong correctness check): C stays bounded as crit → 0 because the ODD blocks vanish at q = 0 — test by approaching the boundary at fixed T.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_flag_off_byte_identical(testCase)
ion = invz_ion();  T = 1.6;  B = [0.5 0 0];
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt  = invzt_jq_tensor(ion, qvec, struct('dpRng', 10, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
p1 = invzt_solve_point(ion, T, B, Jt, struct('JtGamma', JtG));
p2 = invzt_solve_point(ion, T, B, Jt, struct('JtGamma', JtG, 'odd_tier2', false));
p1 = rmfield_ifpresent(p1, {'C', 'tier2_iters', 'tla'});   % helper in test file
p2 = rmfield_ifpresent(p2, {'C', 'tier2_iters', 'tla'});
verifyEqual(testCase, p2, p1);
end

function test_tier2_converges_and_suppresses(testCase)
ion = invz_ion();  T = 1.55;  B = [0.05 0 0];
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt  = invzt_jq_tensor(ion, qvec, struct('dpRng', 10, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 10, 'cache', false));
p1 = invzt_solve_point(ion, T, B, Jt, struct('JtGamma', JtG, 'odd', true));
p2 = invzt_solve_point(ion, T, B, Jt, struct('JtGamma', JtG, 'odd', true, 'odd_tier2', true));
verifyTrue(testCase, p2.converged);
verifyLessThanOrEqual(testCase, p2.tier2_iters, 8);
verifyLessThan(testCase, p2.crit, p1.crit + 1e-9);          % variable moments push toward critical (or null)
fprintf('Tier2: iters = %d, crit %.5f -> %.5f, Delta_eff/Delta = %.4f\n', ...
    p2.tier2_iters, p1.crit, p2.crit, p2.tla.Delta / p1.tl.Delta);
end

function test_C_bounded_near_boundary_slow(testCase)
% IR safety: approach the boundary at fixed T; C must saturate, not diverge.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();  T = 1.2;
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Jt  = invzt_jq_tensor(ion, qvec, struct('dpRng', 30, 'cache', true));
JtG = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 30, 'cache', false));
[Bc, ~] = invzt_critical(ion, T, Jt, JtG, [2 6], struct('odd', true));
Cn = zeros(3,1);  db = [0.5 0.2 0.05];
for i = 1:3
    pt = invzt_solve_point(ion, T, [Bc + db(i), 0, 0], Jt, struct('JtGamma', JtG));
    assumeTrue(testCase, pt.converged);
    C = invzt_odd_fieldvar(pt, ion, T, Jt, struct());
    Cn(i) = norm(C);
end
verifyLessThan(testCase, Cn(3)/Cn(1), 20);                  % grows but bounded
fprintf('C-norm approaching Bc: %.3g / %.3g / %.3g (dB = 0.5/0.2/0.05 T)\n', Cn);
end
```

- [ ] **Step 2: Run to verify failure**, **Step 3: implement the outer loop**, **Step 4: run tests + both fast suites + the slow test on this file**.

- [ ] **Step 5: Measurements → `docs/ODD-LOG.md`** (report): at (1.55 K, Bx = 0.05 T proxy, 16³/dpRng 30): crit with (i) no ODD, (ii) Tier 1 (A1), (iii) Tier 1+2 — and the A1 zero-field Tc via crit(T) root-scan at 0.05 T for all three, giving the combined ΔTc(0) and the Tier1 : Tier2 split. Characterize Bx → 0 (T3.4): does Tier 2's internal field move zero-field results from "directional only" to quantitative, or does the hyperfine-scale caveat stand? Document the answer explicitly.

- [ ] **Step 6: Commit** — `feat(invzt): Tier-2 outer self-consistency with IR-safety gate + combined ODD measurements`.

---

### Task 11: Validation drivers, robustness sweeps, docs, handoff

**Files:**
- Create: `invz_tensor/invzt_run_phase_overlay.m`, `invz_tensor/README.md` (replace the scaffold), `docs/SESSION-2026-07-16-invz-tensor-odd.md`
- Modify: `docs/ODD-LOG.md` (final V4 table), `invz_tensor/tests/test_invzt_critical.m` (add the reproducibility slow test pinning the Task 8/10 logged numbers to 1%)

**Interfaces:**
- Consumes: everything above.
- Produces: `invzt_run_phase_overlay.m` — a driver (same style as `invz_run_phase_diagram`) that overlays: 1/z no-ODD baseline, 1/z+Tier1 (A1), 1/z+Tier1+2, MF, and the experimental points already used by the README §7 comparisons (Bitko 1996 / Babkevich 2016 values as in R2007 Fig. 1); second panel: Σ(0) along the boundary with/without ODD. Coarse preview mode (`opts.quick = true`: 12³ grid, 3 boundary points, minutes) runs in the session; the production sweep (16³/24³, full boundary, hours) is documented and LEFT TO THE USER (repo precedent).

- [ ] **Step 1: Write the driver** with `quick` and production modes; run `quick` once, save the figure to `invz_tensor/` (`.fig`), and record the quick-mode numbers in ODD-LOG.

- [ ] **Step 2: Robustness spot-checks** (each ONE cheap point, logged, not tests): closed-form ΔTc(0) at 12³ vs 16³ vs 24³ (+ Richardson) — note whether ODD convergence is slower than the baseline Σc's; ODD-block dpRng sensitivity 20/30/40 at the 16³ smallest shell (confirm the effective r⁻⁶ short-rangedness of the mediated coupling); GH nodes 5/7/9 on one Tier-2 point. Every headline number in the final table gets a numerical-error bar from these sweeps.

- [ ] **Step 3: Reproducibility slow test**: extend `test_invzt_critical.m` with `test_odd_headline_reproducibility_slow` asserting the Task 8 logged `Tc_ODD(0)` Richardson value and `d` reproduce to 1% (values read from `invzt_anchors()` — ADD them as anchors now, with provenance comments pointing at the ODD-LOG entries).

- [ ] **Step 4: Write `invz_tensor/README.md`**: module map (every `invzt_` function, one line each), flags (`odd`, `chi_rest`, `odd_tier2`, `parts`, `Esplit`, `ngh`), cache-key namespace, the three-layer scope statement (A0+A1+Tier2 implemented; A2/A3 tensor-EMT + tensor-cumulant route deferred — NEVER stack Tiers 1-2 with a future A2+A3, per the Appendix scope box; A0+A1 subsumes T2 retardation), the headline results table, and the "least rigorous step" flag on the Tier-2 dressing. Update `docs/ODD-LOG.md` with the final V4 summary: Tc(0) [1/z 1.74 K → 1/z+ODD x K vs exp 1.53 K], fraction of the 0.21 K gap closed, decomposition, error bars.

- [ ] **Step 5: Write the handoff doc** `docs/SESSION-2026-07-16-invz-tensor-odd.md`: what exists, every new function/flag/cache key, open Tier-3/A2/A3 items (with the framework §10 / J1984 / Part II pointers from the ODD plan §7), and the production runs left to the user.

- [ ] **Step 6: Full verification** — both fast suites AND both slow suites (`INVZ_SLOW=1`) green; the `invz_projected` numbers must match the Task 2 frozen baseline exactly.

- [ ] **Step 7: Commit** — `docs(invzt): phase-overlay driver, robustness error bars, README + handoff (V4)`.

---

## Execution notes for the controller

- Tasks are strictly sequential (each consumes the previous interfaces); no parallel implementers.
- Model guidance: Task 1 mechanical (cheap tier); Tasks 2, 3, 4, 6 standard; Tasks 5, 7, 8, 9 are the physics-critical core (most capable tier; same for their reviews); Tasks 10, 11 standard.
- MATLAB runtime discipline: warm the `jqt` caches once (the Task 7/8 tests do) and keep `'cache', true` in every slow test; expect the first 16³/dpRng-30 tensor build to take minutes (one-off).
- The plan's §0 stop-rule binds every task: when a P0/A0 finding contradicts the physics summary (units, symmetry, χ⊥ band, d band), STOP and report BLOCKED rather than improvising.
