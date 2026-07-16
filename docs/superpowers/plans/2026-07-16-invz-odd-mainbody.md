# ODD Main-Body (P0 → T1 → T2 → T3 → V4) in `invz_projected/` — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute the MAIN BODY of `odd_implementation_plan.html` (phases P0, T1, T2, T3, V4 — NOT Appendix A) against the existing scalar-cc module `invz_projected/`: off-diagonal dipolar (ODD) physics as a strictly additive, opt-in extension — Tier 1a (static ODD-mediated coupling δJ^cc via E1/E4/E5), Tier 1b (retardation sensitivity), Tier 2 (internal-transverse-field renormalization of the doublet), and validation/reporting.

**Architecture:** Four new functions (`invz_odd_blocks`, `invz_chiperp`, `invz_odd_deltaJ`, `invz_odd_fieldvar`) + one averaging wrapper (`invz_twolevel_avg`) + one zero-field measurement engine (`invz_odd_zero_field`), wired behind `opts.odd` / `opts.odd_retarded` / `opts.odd_tier2` into `invz_jq_modes`, `invz_solve_point(_ordered)`, the critical finders, and `invz_run_phase_diagram`. Flag off (default) = bit-for-bit current behavior — that is a hard gate on every task. Physics numbers are *reported, never tuned*; the only physics pass/fail gates are internal identities (PSD, Parseval, closure, sum rule) and flag-off regressions against published benchmarks.

**Tech Stack:** MATLAB R2025a, `functiontests(localfunctions)` suites, existing module + root lattice machinery (`MF_dipole.m`, `exchange.m`).

**User decision of record (2026-07-16):** the full-tensor Appendix A route (`invz_tensor/`, plan `2026-07-16-invz-tensor-odd.md`) is DEFERRED; this main-body plan is the active stream. `invz_tensor/` stays a scaffold; no `invz_common/` refactor now (single active branch).

## Global Constraints

- Repo root: `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion/`. The path contains spaces — ALWAYS quote it in shell commands.
- Fast suite (run from repo root; EVERY task ends with it green, counts vs baseline 109 passed / 0 failed / 12 incomplete at `5f4ff92` — new tests may raise the passed count, nothing may fail):
  `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "r = runtests('invz_projected/tests'); disp(r); assertSuccess(r)"`
  Slow gate: prefix `INVZ_SLOW=1`; slow tests use `assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), ...)`.
- Non-negotiables (ODD plan §0): (i) flag OFF (default) reproduces every existing README §7 benchmark bit-for-bit — the extension is strictly additive and opt-in; (ii) never edit published benchmark values in existing tests; (iii) criteria marked *report* are measurements — record the number, do not tune. When a P0 finding contradicts the ODD plan §1 (units, symmetry, χ⊥ band), STOP and report BLOCKED before improvising.
- State tracking (ODD plan §0): each task's commit ALSO flips the `[ ]` → `[x]` in the task headers of `odd_implementation_plan.html` for the task IDs it completes (P0.1, T1.1, …) and sets the matching `status` fields to `"done"` in the `#task-manifest` JSON at the bottom of that file. Task IDs go in commit messages.
- Reporting (ODD plan §0): every phase appends a dated entry to `docs/ODD-LOG.md` (Task 1 creates it): what was implemented, test status, headline numbers (ΔTc, Σc shifts, timings).
- Units & conventions (must match the codebase — ODD plan §1.3): energies meV; `C = invz_const()`; χ = −G, χ0 in meV⁻¹, crit dimensionless; ferromagnetic J > 0, δJ ⪰ 0 adds to J^cc(q); 4 sublattices, cc channel = 4×4 per q → `Jnu [nq,4]`; J(ii) = 0 enforced for δJ by per-sublattice q-average subtraction (E4); Lorentz/demag are Cartesian-diagonal → NO ODD-block shape terms, δJ is a pure lattice sum; `J^{cα}(q)_{ss'} = -C.gfac*dip(3,α,s,s')` with `dip` from `MF_dipole` (full `[3,3,4,4]`, Å⁻³; α: 1 = a, 2 = b), exchange is isotropic (Cartesian-diagonal) and contributes NOTHING off-diagonal.
- Published anchors (NEVER rescale): `info.Jcc0_dipole = 6.821e-3` meV, `info.Jaa0_dipole = 3.912e-3` meV, `info.Jcc0 = 6.421e-3` meV; Σc(0) Richardson(12³,24³) = 0.2980 (published 0.3004, AbsTol 0.006); Tc(0) = 1.743 K (published 1.74). DS2023 anchors (safe to assert: convention-free): lattice sums (Suppl. Table I, a = 5.175 Å) `B_xz,xz = B_yz,yz = 36.73/a^6`, `B_zz,zz = 17.93/a^6`; 3-state χ⊥ = 2ρ²/Δ ≈ 11.05 meV⁻¹ (Δ = 11.5 K, ρ = 2.34); full-CF expectation χ⊥(0) ≈ 16–17 meV⁻¹ at (1.53 K, 0 T); expected d ≈ 0.3–0.5 μeV (report-band). NEVER import DS2023's 0.805 longitudinal rescaling or their perturbative hyperfine treatment (plan §8).
- q = 0 pitfall (plan §8, binds Tasks 2, 4, 5): driver grids exclude q = 0; the ODD blocks' q→0 limit is direction-dependent (macroscopic term ∝ q̂_c q̂_α, scale `4*pi*C.gfac/ion.Vc` ≈ 3.7 μeV, elements up to ≈ 1.8 μeV on tilted rays) — decay assertions hold on HIGH-SYMMETRY on-axis rays ONLY; never special-case q = 0 by grid extrapolation: the ONLY q = 0 modification is the explicit constant −d applied to `info.Jcc0` / `J0eff` (E5), exactly once (post-subtraction δJ already carries −d on its diagonal; no other q = 0 handling).
- Physics sign contract (derived from `pt.crit = 1 + Sigma0 - J0eff*chi0cc0`): ODD suppresses ordering — at a fixed converged PM point, `crit` INCREASES with ODD on (J0eff loses d; Σ gains fluctuation weight); Tc(0) and Bc(T) DECREASE. Tests gate on these directions, never on magnitudes.
- Anchors fixture: Task 1 creates `invz_projected/tests/invz_odd_anchors.m` (a plain function returning a struct of controller-verified digits measured on THIS tree, full-precision literals with provenance comments). Later tests reference `invz_odd_anchors()` fields instead of hard-coding unpinned digits. Never edit a pinned anchor to make a test pass — a mismatch is a finding to escalate.
- Cache discipline: ODD geometric blocks get their OWN key namespace in `invz_projected/cache/`: `odd1_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat` (schema tag `odd1`; same `hash_vec` helper style as `invz_jq_modes`), `pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]`, file stores `pkey` + `qvec` and the loader `isequal`-verifies both (cache-hardening precedent). One cache file per grid holds ALL blocks (Vca, Vcb, Vcc). Existing `jq4_` caches are never touched. Library functions stay serial; drivers compute all lattice sums BEFORE entering `parfor` (workers do no disk I/O — P0.4).
- Fast tests use small grids (n ≤ 8 per axis, dpRng ≤ 15, `'cache', false` unless the test IS the cache test); fast-suite additions across ALL tasks total < 30 s (V4.3 budget). Uniform meshes in tests are built directly: `n = 6; [h, k, l] = ndgrid((0:n-1)/n); qvec = [h(:) k(:) l(:)];` (include or drop Γ per test — say which). `qVec_generator`'s grid param is `'grid'` (NOT `'size'`) and it `fprintf`s diagnostics — wrap in `evalc` where output must be clean.
- Known island: avoid convergence-gated fixtures at (T = 0.31 K, B ∈ [1, 2] T).
- Test boilerplate (every new test file in `invz_projected/tests/`):

```matlab
function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
end
```

- Commit style: conventional prefix + `(invz)` scope + task IDs, e.g. `feat(invz): T1.1 geometric ODD lattice sums`, ending with the trailer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.

---

### Task 1: P0 preflight — anchors fixture, exploratory measurements, baseline freeze (P0.1–P0.5)

**Files:**
- Create: `docs/ODD-LOG.md`, `invz_projected/tests/invz_odd_anchors.m`, `invz_projected/tests/exploratory/explore_chiperp.m`, `invz_projected/tests/exploratory/explore_odd_blocks.m`
- Modify: NOTHING in module code (read-only phase).

**Interfaces:**
- Produces: `A = invz_odd_anchors()` returning a struct with AT LEAST: `A.chiperp_1p53K_0T` (2×2 meV⁻¹: symmetrized (a,b)×(a,b) block of `invz_chi0z(si, T, 0, struct('elastic', true))` at T = 1.53 K, B = [0 0 0], `hyp = true`, default `transverse_mf`), `A.chiperp_asym_1p53K` (max abs antisymmetric part), `A.chiperp_elastic_share_1p53K`, `A.chiperp_0p31K_Bx` (struct `.Bx = 0:6`, `.chi_aa` values), `A.odd_onaxis_smallq` (struct `.q = [1e-1 1e-2 1e-3]` along [q 0 0], `.maxca` = max |J^{ca}| element in meV per q at dpRng 30), `A.odd_generic_q_max` (max |J^{ca}| element at q = [0.31 0.17 0.09], dpRng 30), `A.jcc0`/`A.jaa0` (info values from `invz_jq_modes`, dpRng 30).
- Produces: `docs/ODD-LOG.md` §P0 with the five numbered findings below + the frozen baseline.

- [ ] **Step 1 (P0.1): confirm the full dipole tensor upstream of the cc projection.** Read `MF_dipole.m` and `invz_jq_modes.m`; record in ODD-LOG: `MF_dipole` returns `[3,3,4,4]` (`dip(mu,nu,s,s')`, Å⁻³ pre-gfac; index layout stated), the cc extraction line (`invz_jq_modes.m` line ~66: `-squeeze(C.gfac*dip(3,3,:,:))`), how `info.Jcc0` is assembled (Γ dipole + Lorentz `4*pi/(3*ion.Vc)*C.gfac` + `4*ion.J12`), and run the grep from the plan (`grep -rn "Jxx0\|Jcc0\|MF_dipole" invz_projected/*.m`) summarizing where couplings flow.

- [ ] **Step 2 (P0.2): pin units + χ⊥.** Write and run `explore_chiperp.m` (plain script, addpaths as in the boilerplate): at (1.53 K, [0 0 0], hyp true) compute the (a,b) block of `invz_chi0z` at z = 0, print symmetric/antisymmetric split, elastic on/off variants; sweep Bx = 0:6 T at 0.31 K printing χ_aa, χ_bb. Verify dimensional closure: χ0 in meV⁻¹, Jnu/Jcc0 in meV (gfac carries (gL·μB)² — cross-check the 6.821 μeV anchor), so E1 needs NO extra g-factors. Expected χ_aa = χ_bb ≈ 16–17 meV⁻¹ (≈ 11 suggests truncation; ×2-scale off suggests convention slip → STOP/BLOCKED).

- [ ] **Step 3 (P0.3): ODD symmetry check.** Write and run `explore_odd_blocks.m` using `MF_dipole` directly (dpRng 30, geom reuse): (i) on-axis rays [q 0 0] and [0 0 q], q ∈ {1e-1, 1e-2, 1e-3}: max |−gfac·dip(3,1,:,:)| and (3,2) — confirm decay to 0, record the rate and residual scale vs `Jcc0`; (ii) tilted ray q·[1 0 1]/√2 — record the NON-decaying direction-dependent limit vs the `4*pi*C.gfac/ion.Vc` scale (this bounds what small-q assertions may claim); (iii) generic q = [0.31 0.17 0.09] magnitudes vs Jcc0; (iv) the smallest shells of the standard 16³ grid along a*/c*.

- [ ] **Step 4 (P0.4 + P0.5): cache audit + baseline freeze.** Record in ODD-LOG: the `invz_jq_modes` cache key scheme (`jq4_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat`, pkey schema v4) and the ODD namespace decision (`odd1_`, Global Constraints); the pre-parfor discipline point (drivers build all lattice sums before `parfor` — cite `invz_run_phase_diagram`'s structure). Run the fast suite AND the slow suite (`INVZ_SLOW=1`) on the clean tree; paste the totals + benchmark console lines + git hash into ODD-LOG as the frozen baseline. Commit nothing before this step's runs complete.

- [ ] **Step 5: write `invz_odd_anchors.m`** from the measured numbers (full-precision literals, provenance comments: script, date, git hash, dpRng).

- [ ] **Step 6: Commit** — `docs(invz): P0.1-P0.5 preflight — ODD-LOG, anchors fixture, exploratory scripts, frozen baseline` (flip P0.1–P0.5 checkboxes + manifest in `odd_implementation_plan.html` in this commit).

---

### Task 2: `invz_odd_blocks.m` — geometric ODD lattice sums, cached (T1.1)

**Files:**
- Create: `invz_projected/invz_odd_blocks.m`, `invz_projected/tests/test_invz_odd_blocks.m`
- Reference (read, do not modify): `invz_projected/invz_jq_modes.m` (mirror: geom priming, per-q loop, Γ handling, cache pattern, `hash_vec`)

**Interfaces:**
- Consumes: `MF_dipole`, `exchange`, `invz_const`, `invz_is_gamma_equiv`.
- Produces: `[Vca, Vcb, Vcc, info] = invz_odd_blocks(ion, qvec, opts)`.
  - `Vca, Vcb` `[4,4,nq]` complex: `Vca(s,s',iq) = -C.gfac*dip(3,1,s,s')` at q(iq,:) — DIPOLE-ONLY (exchange is Cartesian-diagonal; Lorentz/demag are diagonal → no ODD shape terms; do NOT Hermitize individual blocks — J^{ca}(q) is not Hermitian; the pair identity is `J^{ca}(-q) = conj(J^{ca}(q))` and `J^{ac}(q) = J^{ca}(q)'`).
  - `Vcc` `[4,4,nq]` Hermitian: the SAME assembly `invz_jq_modes` eigendecomposes — `-C.gfac*dip(3,3,:,:) + sign(ion.J12)*ex(3,3,:,:)` + `lorz` at Γ-equivalent q, Hermitized per q. (Controller-approved adaptation of the plan's 2-output spec: `invz_jq_modes` caches only eigenvalues, but the ODD path needs the cc MATRICES to add δJ before `eig`; caching Vcc alongside halves I/O per the plan's own one-file rule. Verified against `invz_jq_modes` by the parity test below.)
  - `info`: `dpRng`, `Jcc0`, `Jaa0`, `Jcc0_dipole`, `Jaa0_dipole` (same Γ info block as `invz_jq_modes`, computed from the priming call), `lorz`.
  - `opts`: `dpRng` (30), `cache` (true). Cache per Global Constraints (`odd1_` namespace, stores `Vca`,`Vcb`,`Vcc`,`info`,`pkey`,`qvec`, isequal-verified load). `ion.demag` must be 0 or absent-equivalent for the Vcc Γ info to match the intrinsic convention — if `ion.demag ~= 0`, `error('invz:oddDemag', ...)` (ODD layer is intrinsic-only; demag handling stays in `invz_jq_modes`).

- [ ] **Step 1: Write the failing tests** — `test_invz_odd_blocks.m` (boilerplate setupOnce; `ion = invz_ion()` in each):

```matlab
function test_shapes_and_conj_symmetry(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; -0.25 0 0; -0.31 -0.17 -0.09];
[Vca, Vcb, Vcc] = invz_odd_blocks(ion, q, struct('dpRng', 15, 'cache', false));
verifyEqual(testCase, size(Vca), [4 4 4]);  verifyEqual(testCase, size(Vcb), [4 4 4]);
for iq = 1:2   % real-space couplings real => J(-q) = conj(J(q))
    verifyLessThan(testCase, max(abs(Vca(:,:,iq+2) - conj(Vca(:,:,iq))), [], 'all'), 1e-12);
    verifyLessThan(testCase, max(abs(Vcb(:,:,iq+2) - conj(Vcb(:,:,iq))), [], 'all'), 1e-12);
    verifyLessThan(testCase, norm(Vcc(:,:,iq) - Vcc(:,:,iq)', 'fro'), 1e-14);   % cc Hermitian
end
end

function test_vcc_parity_with_jq_modes(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];      % includes a Gamma-equivalent point
[~, ~, Vcc, infoB] = invz_odd_blocks(ion, q, struct('dpRng', 15, 'cache', false));
[Jnu, infoS] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
for iq = 1:3
    verifyEqual(testCase, sort(real(eig(Vcc(:,:,iq)))).', Jnu(iq,:), 'AbsTol', 1e-12);
end
verifyEqual(testCase, infoB.Jcc0, infoS.Jcc0, 'RelTol', 1e-12);
verifyEqual(testCase, infoB.Jaa0, infoS.Jaa0, 'RelTol', 1e-12);
end

function test_ds2023_geometry_sums(testCase)
% DS2023 Suppl. Table I (a = 5.175 Ang): pure-geometry, gfac-free real-space
% sums — THE unit guard (gfac slips enter deltaJ SQUARED). Central ion on
% sublattice 1; geom stores lower-triangular pairs, so {s,1} exists for s=1..4
% and covers all four sublattice partners. Tf(:,n,m) = 3 r_n r_m/r^5 - d_nm/r^3.
ion = invz_ion();
[~, ~, geom] = MF_dipole([0 0 0], 30, ion.a, ion.tau);
a = 5.175;
[Sxz, Syz, Szz] = deal(0);
for s = 1:4
    Tf = geom.Tf{s, 1};
    Sxz = Sxz + sum(Tf(:,3,1).^2);
    Syz = Syz + sum(Tf(:,3,2).^2);
    Szz = Szz + sum(Tf(:,3,3).^2);
end
verifyEqual(testCase, Sxz, 36.73/a^6, 'RelTol', 0.01);
verifyEqual(testCase, Syz, 36.73/a^6, 'RelTol', 0.01);
verifyEqual(testCase, Szz, 17.93/a^6, 'RelTol', 0.01);
end

function test_onaxis_smallq_decay(testCase)
% C2-about-c kills the ODD blocks on-axis as q -> 0. ON-AXIS ONLY (tilted rays
% carry a direction-dependent macroscopic limit — plan SS8/P0.3).
% P0 AMENDMENT (ODD-LOG SSP0.3): element decay is LINEAR in q (sublattice phase
% factors; the macroscopic term vanishes on-axis), so the source plan's
% 1e-6*Jcc0 element gate at q = 1e-3 is unachievable — its own escape clause
% ("or the deviation is explained by grid geometry") applies. Gates: pinned
% values, linear decay structure, and E1-relevant smallness of the SQUARE
% (deltaJ ~ chi_perp*|Jca|^2 is what must vanish vs Jcc0).
ion = invz_ion();
A = invz_odd_anchors();
q = [1e-1 0 0; 1e-2 0 0; 1e-3 0 0];
Vca = invz_odd_blocks(ion, q, struct('dpRng', 30, 'cache', false));
m = arrayfun(@(iq) max(abs(Vca(:,:,iq)), [], 'all'), 1:3);
verifyEqual(testCase, m(:), A.odd_onaxis_smallq.maxca(:), 'RelTol', 1e-6);   % pinned P0 digits
verifyEqual(testCase, m(2)/m(1), 0.1, 'RelTol', 0.25);                       % ~linear decade steps
verifyEqual(testCase, m(3)/m(2), 0.1, 'RelTol', 0.25);
verifyLessThan(testCase, 18 * m(3)^2, 1e-5 * 6.421e-3);                      % chi_perp*|Jca|^2 << Jcc0
end

function test_cache_roundtrip_selfverifying(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
[V1a, V1b, V1c] = invz_odd_blocks(ion, q, struct('dpRng', 10, 'cache', true));
[V2a, V2b, V2c] = invz_odd_blocks(ion, q, struct('dpRng', 10, 'cache', true));
verifyEqual(testCase, {V2a, V2b, V2c}, {V1a, V1b, V1c});      % bitwise round-trip
ion2 = ion;  ion2.J12 = ion.J12 * 1.05;                        % physics change must miss
% AMENDED (Task 2 adjudication): Vca/Vcb are dipole-only, hence J12-INDEPENDENT
% by this plan's own interface — the miss is observable on the cc block, which
% carries exchange. A wrong cache HIT would leave V3c == V1c.
[~, ~, V3c] = invz_odd_blocks(ion2, q, struct('dpRng', 10, 'cache', true));
verifyFalse(testCase, isequal(V3c, V1c));
cdir = fullfile(fileparts(mfilename('fullpath')), '..', 'cache');
verifyTrue(testCase, ~isempty(dir(fullfile(cdir, 'odd1_*.mat'))));
end

function test_parseval_odd_vs_realspace_slow(testCase)
% T1.1(iv): BZ-average of the squared ca blocks == real-space squared sum
% (Parseval), SAME dpRng both sides; 1% tolerance absorbs superlattice folding
% at n = 8 (r^-6-suppressed).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();  C = invz_const();
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];   % uniform mesh INCLUDING Gamma
Vca = invz_odd_blocks(ion, qvec, struct('dpRng', 20, 'cache', true));
lhs = mean(sum(abs(Vca(1, :, :)).^2, 2), 3);                       % row s = 1, sum over s', avg q
[~, ~, geom] = MF_dipole([0 0 0], 20, ion.a, ion.tau);
rhs = 0;
for s = 1:4
    Tf = geom.Tf{s, 1};
    rhs = rhs + sum((C.gfac * Tf(:,3,1)).^2);
end
verifyEqual(testCase, lhs, rhs, 'RelTol', 0.01);
% independent conversion check via the cc channel (validates gfac placement):
[~, ~, Vcc] = invz_odd_blocks(ion, qvec, struct('dpRng', 20, 'cache', true));
% subtract exchange+Lorentz first: rebuild the dipole-only cc from MF_dipole per q
% is overkill — instead run the same Parseval on Vca vs B_xz,xz in absolute terms:
verifyEqual(testCase, lhs, C.gfac^2 * 36.73/5.175^6, 'RelTol', 0.015);
end
```

Note on `test_ds2023_geometry_sums`: verify the `{s,1}` cell indexing against `MF_dipole`'s storage before trusting it (`geom.Tf{nt,mt}` filled for `mt <= nt`; `{1,1}` self-sublattice with self-site excluded via the `r2 < 0.01` cut). If the sums land off the DS2023 values by a uniform integer-ish factor, STOP and report BLOCKED with the measured ratio — do not fold in a fudge factor.

- [ ] **Step 2: Run the test file — verify it fails** (undefined function).
- [ ] **Step 3: Implement `invz_odd_blocks.m`** mirroring `invz_jq_modes` (demag guard; geom priming at [0 0 0] capturing `dip0` for info; per-q loop filling the three block sets; Γ-equivalence via `invz_is_gamma_equiv` for Vcc's Lorentz; cache write/verify per Global Constraints).
- [ ] **Step 4: Run tests (fast) — pass; run the full fast suite — green.**
- [ ] **Step 5: Run the slow-gated Parseval test once** (`INVZ_SLOW=1`, this file only); record lhs/rhs and residual in ODD-LOG.
- [ ] **Step 6: Commit** — `feat(invz): T1.1 geometric ODD lattice sums (cached) with DS2023 + Parseval unit guards` (flip T1.1 checkbox + manifest).

---

### Task 3: `invz_chiperp.m` — static (and Matsubara-ready) transverse single-ion susceptibility (T1.2)

**Files:**
- Create: `invz_projected/invz_chiperp.m`, `invz_projected/tests/test_invz_chiperp.m`

**Interfaces:**
- Consumes: `invz_single_ion`, `invz_chi0z`, `invz_field_vec`, `getf`.
- Produces: `[Xp, info] = invz_chiperp(ion, T, B, opts)` — `Xp` `[2,2]` real symmetric (meV⁻¹): symmetrized (a,b)×(a,b) block of the full electronuclear `invz_chi0z(si, T, z, ...)`; `opts`: `hyp` (true), `transverse_mf` ('legacy_x'), `elastic` (true), `Jxx0` (ion.Jxx0, forwarded to the single-ion MF), `z` (default 0; a vector of complex frequencies returns `Xp [2,2,nz]` — symmetrize and real-part per slice; used by T2). `info`: `asym` (max abs antisymmetric part — gyrotropic, discarded; assert small downstream), `elastic_share`, `si` (the converged single-ion state, so callers can reuse it). DESIGN DECISION (record verbatim in header, plan T1.2): χ⊥ is Van Vleck-dominated and insensitive to the cc order parameter and K — computed ONCE per (T, Bx) BEFORE the EMT⇆Σ loop, never self-consistently inside it. Must NOT route through `invz_twolevel` (the transverse response crosses the ~10 K CF gap and is regular where the two-level object is not; the `invz:degenerateDoublet` guard must never fire here — this function works at Bx = 0).

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_anchors_and_symmetry(testCase)
ion = invz_ion();  A = invz_odd_anchors();
[Xp, info] = invz_chiperp(ion, 1.53, [0 0 0], struct());
verifyEqual(testCase, Xp, A.chiperp_1p53K_0T, 'RelTol', 1e-9);
verifyEqual(testCase, Xp(1,1), Xp(2,2), 'AbsTol', 1e-10*abs(Xp(1,1)));   % C4 at Bx = 0
verifyEqual(testCase, Xp, Xp.', 'AbsTol', 1e-15);
verifyLessThan(testCase, info.asym, 1e-8*abs(Xp(1,1)));
verifyGreaterThan(testCase, Xp(1,1), 10);   % the 16-17 meV^-1 band, floor guard
verifyLessThan(testCase, Xp(1,1), 25);
end

function test_zero_field_no_doublet_guard(testCase)
% Regular at Bx = 0 (never routes through invz_twolevel).
ion = invz_ion();
Xp = invz_chiperp(ion, 0.31, [0 0 0], struct());
verifyTrue(testCase, all(isfinite(Xp(:))));
end

function test_reproducible_along_Bx(testCase)
% P0 AMENDMENT (ODD-LOG SSP0.2): the 0.31 K chi_aa(Bx) curve has a PHYSICAL
% peak near 1 T (halves by 2 T; all points MF-converged), so a step-size
% "smoothness" gate is wrong physics. The anchor-pinned sweep IS the
% no-numerical-artifact guard.
ion = invz_ion();  A = invz_odd_anchors();
x = zeros(1, 7);
for i = 0:6, Xp = invz_chiperp(ion, 0.31, [i 0 0], struct()); x(i+1) = Xp(1,1); end
verifyEqual(testCase, x(:), A.chiperp_0p31K_Bx.chi_aa(:), 'RelTol', 1e-9);   % pinned P0 sweep
verifyTrue(testCase, all(isfinite(x)) && all(x > 0));
end

function test_matsubara_frequencies(testCase)
% z-vector form: real symmetric per slice, decaying with |w_n|, r_n in (0, 1].
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.53, 10);
[Xp, ~] = invz_chiperp(ion, 1.53, [0 0 0], struct('z', 1i*wn));
verifyEqual(testCase, size(Xp), [2 2 numel(wn)]);
r = squeeze(Xp(1,1,:)) / Xp(1,1,1);
verifyEqual(testCase, r(1), 1, 'AbsTol', 1e-14);
verifyTrue(testCase, all(diff(r) < 0) && all(r > 0));   % monotone decay along iw_n
% rough scale from the plan SS4: r(w1 = 0.828 meV) ~ eps1^2/(eps1^2 + w1^2) ~ 0.56
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement** (si once via `invz_single_ion(ion, T, invz_field_vec(B), struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf))`; `c0 = invz_chi0z(si, T, z, struct('elastic', elastic))`; slice (1:2,1:2,:); symmetrize; realness assert `imag < 1e-10*max` then `real()`). **Step 4: tests + full fast suite green.**
- [ ] **Step 5: Commit** — `feat(invz): T1.2 static transverse single-ion susceptibility invz_chiperp` (flip T1.2 + manifest).

---

### Task 4: `invz_odd_deltaJ.m` + `invz_jq_modes` ODD extension (T1.3)

**Files:**
- Create: `invz_projected/invz_odd_deltaJ.m`
- Modify: `invz_projected/invz_jq_modes.m` (opts.odd branch ONLY — the default path must stay byte-identical)
- Test: extend `invz_projected/tests/test_invz_odd_blocks.m` (plan T1.3 says extend this file)

**Interfaces:**
- Consumes: Tasks 2–3 outputs.
- Produces:
  - `[dJ, d, dinfo] = invz_odd_deltaJ(Vca, Vcb, Xp)` — E1 then E4/E5 over the SUPPLIED grid (caller contract in header: qvec behind the blocks must be a full uniform BZ mesh; driver grids excluding Γ are fine — δJ(Γ)→0 anyway, document the O(1/Nq) difference):
    E1 per q: `dJpre = Vca*Xp(1,1)*Vca' + Vca*Xp(1,2)*Vcb' + Vcb*Xp(2,1)*Vca' + Vcb*Xp(2,2)*Vcb'`, Hermitized.
    E4: `dJ = dJpre` with `dJ(s,s,:) = dJpre(s,s,:) - mean_q dJpre(s,s,:)` (diagonal only).
    E5: `d = mean_s mean_q dJpre(s,s,:)` (scalar, meV, > 0).
    `dinfo`: `d_per_sublattice` [4,1], `presub_min_eig` (min over q of min eig of dJpre), `postsub_diag_bzavg` (max |per-sublattice BZ-avg of diag(dJ)|), `dJ_max`. COMMENT TRAIL (verbatim requirement, plan T1.3/§8): the subtracted on-site constant multiplies (σᶻ)² = 1 in the strict two-level limit — a pure energy shift with no Tc content; its physical residue (internal field renormalizing Δ, M², n01) is exactly what Tier 2 owns — cite E4/E5 and this plan; no double counting.
  - `invz_jq_modes(ion, qvec, opts)` gains `opts.odd = false (default) | struct('Xp', [2 2] double)`:
    - `false`/absent: byte-identical current code path (the diff must not touch it).
    - struct: after resolving dpRng/cache, call `[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', dpRng, 'cache', useCache))`; `[dJ, d] = invz_odd_deltaJ(Vca, Vcb, opts.odd.Xp)`; per q `M = Vcc(:,:,iq) + dJ(:,:,iq)` (Hermitize), `Jnu(iq,:) = sort(real(eig(M)))`, `Juni(iq) = real(v'*M*v)`; `info.Jcc0 = infoB.Jcc0 - d` (THE single line carrying the MF-level DS2023 mechanism into MF, RPA denominator, and 1/z criticality — single-sourcing; bookkeeping rule: grid matrices carry post-subtraction δJ (already contains −d on the diagonal), Jcc0 carries the explicit −d, NO other q = 0 handling); `info.odd = struct('d', d, 'dJ_mean_diag', <[4,1] pre-sub diag means>, 'dJ_max', ..., 'uniform_residual', <max |dJpre| on the smallest-|q| shell of qvec>, 'Xp', opts.odd.Xp)`. The ODD path does NOT write to the `jq4_` cache (its Jnu depend on Xp; the geometric blocks are already cached under `odd1_`).

- [ ] **Step 1: Write the failing tests** (append to `test_invz_odd_blocks.m`):

```matlab
function test_deltaJ_identities(testCase)
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];  % driver-style, no Gamma
[Vca, Vcb] = invz_odd_blocks(ion, qvec, struct('dpRng', 15, 'cache', false));
Xp = invz_chiperp(ion, 1.53, [0 0 0], struct());
[dJ, d, dinfo] = invz_odd_deltaJ(Vca, Vcb, Xp);
verifyGreaterThan(testCase, d, 0);
verifyEqual(testCase, dinfo.d_per_sublattice, d*ones(4,1), 'RelTol', 1e-10);       % T1.3 assert
verifyGreaterThan(testCase, dinfo.presub_min_eig, -1e-12*max(1e-30, dinfo.dJ_max)); % PSD pre-sub
verifyLessThan(testCase, dinfo.postsub_diag_bzavg, 1e-15*max(1e-30, dinfo.dJ_max)); % E4 exact
for iq = [1, round(size(dJ,3)/2)]
    verifyLessThan(testCase, norm(dJ(:,:,iq) - dJ(:,:,iq)', 'fro'), 1e-14);         % Hermitian
    verifyLessThan(testCase, max(abs(imag(eig(dJ(:,:,iq))))), 1e-12);               % real eigs
end
% E4/E5 closure (T1.3 acceptance vi): BZ-avg of the pre-subtraction diagonal, as
% reconstructed from the assembled outputs dJ(s,s,:) + d, recovers d per sublattice.
verifyEqual(testCase, squeeze(mean(real(dJ(1,1,:)), 3)) + d, d, 'AbsTol', 1e-15*max(1e-30, d)); %#ok<NASGU> % diag avg is 0
% Xp = 0 -> dJ = 0, d = 0
[dJ0, d0] = invz_odd_deltaJ(Vca, Vcb, zeros(2));
verifyEqual(testCase, max(abs(dJ0(:))), 0);  verifyEqual(testCase, d0, 0);
end

function test_jq_modes_odd_off_bitwise(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09; 0 0 0];
[J1, i1, u1] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false));
[J2, i2, u2] = invz_jq_modes(ion, q, struct('dpRng', 15, 'cache', false, 'odd', false));
verifyTrue(testCase, isequaln({J1, i1, u1}, {J2, i2, u2}));                          % regression (i)
end

function test_jq_modes_odd_zero_Xp_equals_off(testCase)
ion = invz_ion();
n = 4;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[J1, i1] = invz_jq_modes(ion, qvec, struct('dpRng', 10, 'cache', false));
[J2, i2] = invz_jq_modes(ion, qvec, struct('dpRng', 10, 'cache', false, 'odd', struct('Xp', zeros(2))));
verifyEqual(testCase, J2, J1, 'AbsTol', 1e-14);                                      % regression (ii)
verifyEqual(testCase, i2.Jcc0, i1.Jcc0, 'AbsTol', 1e-18);
verifyEqual(testCase, i2.odd.d, 0);
end

function test_jq_modes_odd_on_structure(testCase)
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
Xp = invz_chiperp(ion, 1.53, [0 0 0], struct());
[J1, i1] = invz_jq_modes(ion, qvec, struct('dpRng', 15, 'cache', false));
[J2, i2] = invz_jq_modes(ion, qvec, struct('dpRng', 15, 'cache', false, 'odd', struct('Xp', Xp)));
verifyGreaterThan(testCase, i2.odd.d, 0);
verifyEqual(testCase, i2.Jcc0, i1.Jcc0 - i2.odd.d, 'AbsTol', 1e-18);                 % E5 shift
% ordering mode stays the uniform branch (T1.3 acceptance v): margin shrinks
% from both sides but must remain positive on this grid
verifyLessThan(testCase, max(J2(:)), i2.Jcc0);
% finite-q modes gain weight: mean top-branch coupling strictly increases
verifyGreaterThan(testCase, mean(J2(:,4)), mean(J1(:,4)));
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement** (`invz_odd_deltaJ` first; then the `opts.odd` branch in `invz_jq_modes` — keep the default path's code untouched line-for-line). **Step 4: tests + full fast suite green** (the bitwise regression test is the gate for non-negotiable (i)).
- [ ] **Step 5: Record in ODD-LOG**: d on the 6³/dpRng-15 test grid AND on the production 16³/dpRng-30 grid (expected band 0.3–0.5 μeV — report the number; an order-of-magnitude miss with Task 2's geometry sums green means the χ⊥ contraction is at fault → STOP/BLOCKED), `uniform_residual` vs the 1e-4·Jcc0 claim (on-axis shells only).
- [ ] **Step 6: Commit** — `feat(invz): T1.3 ODD-mediated coupling deltaJ (E1/E4/E5) behind invz_jq_modes opts.odd` (flip T1.3 + manifest).

---

### Task 5: Wire ODD into the point solve and critical finders (T1.4)

**Files:**
- Modify: `invz_projected/invz_solve_point.m`, `invz_projected/invz_solve_point_ordered.m`, `invz_projected/invz_crit_at.m` (read it first — thread opts through unchanged if it already passes opts to the solver), `invz_projected/invz_critical_T0field.m`, `invz_projected/invz_ion.m` (add `ion.odd = 0;` documented default, mirroring `ion.demag`), `invz_projected/invz_run_phase_diagram.m`
- Create: `invz_projected/tests/test_invz_odd_solve.m`

**Interfaces:**
- Produces:
  - `invz_solve_point(ion, T, Bx, Jnu_flat, opts)` gains `opts.odd` (false default) and `opts.odd_blocks`:
    - `opts.odd = true` REQUIRES `opts.odd_blocks = struct('Vca', ..., 'Vcb', ..., 'Vcc', ..., 'Jcc0', <unshifted info.Jcc0>)` (precomputed once by the caller from `invz_odd_blocks` — P0.4 discipline: no cache/disk reads inside solver loops) and REQUIRES `Jnu_flat == []` (assert `error('invz:oddArgs', ...)` otherwise — prevents a stale baseline Jnu silently overriding the rebuild).
    - When on: `Xp = invz_chiperp(ion, T, Bx, struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf))` (the SAME single-ion options as the rest of the solve — T1.2's "same converged single-ion state" requirement); `[dJ, d] = invz_odd_deltaJ(...)`; rebuild `Jnu_flat = sort-eig over nq` (4×4 eig per q, O(ms) at 16³); `J0eff_use = getf(opts, 'J0eff', ion.J0eff) - d` — CALLERS PASS THE UNSHIFTED `info.Jcc0` as `opts.J0eff`, exactly as today; the −d shift is applied HERE, once (assert-comment citing the T1.3 bookkeeping rule). `pt` gains `pt.odd = struct('d', d, 'Xp', Xp)`; `pt.crit` uses `J0eff_use`.
    - Flag off: byte-identical (all new code behind `if`).
  - `invz_solve_point_ordered`: same contract (its `J0eff` also seeds the single-ion `J0z` — the shifted value is the physical uniform coupling and is what both uses receive).
  - `invz_critical` / `invz_critical_T`: NO signature change — they already pass remaining opts through `invz_crit_at` to the solver; verify `invz_crit_at` forwards `opts` and `Jnu_flat = []` cleanly; with `opts.odd` on, callers pass `Jnu_flat = []`. `invz_critical_T`'s adaptive Tc0 anchor: with ODD on it must NOT silently use the no-ODD anchor — require `opts.Tc0` when `opts.odd` is set (`error('invz:oddTc0', ...)`), drivers compute it via `invz_odd_zero_field` (Task 6).
  - `invz_critical_T0field(ion, Sc, J0eff)`: generalized so `Sc` and `J0eff` may each be a numeric (current behavior, byte-identical) OR a `function_handle` of T — `f = @(T) J0(T)*static_chi_cc(ion, T) - (1 + Sc(T))` with numerics wrapped as constant handles internally. This is the plan's "replace direct inversion with a scalar root find on J(0)·χ0cc(0;T) − 1 − Σc(0;T) = 0"; χ⊥(T) is nearly flat so the existing bisection converges unchanged.
  - `invz_ion`: `ion.odd = 0;` + comment block (opt-in ODD switch, mirroring `ion.demag`; drivers read it, libraries read `opts.odd`).
  - `invz_run_phase_diagram`: when `ion.odd`, build `[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, ...)` ONCE at the same pre-parfor point where `invz_jq_modes` is called today, and thread `opts.odd = true`, `opts.odd_blocks = ...`, `Jnu_flat = []`, `opts.Tc0 = <odd-aware>` into the finder calls. Flag off: byte-identical driver behavior.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_solve_point_flag_off_bitwise(testCase)
% Default-opts call must take the identical code path (non-negotiable i).
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';               % synthetic fixture (house convention)
p1 = invz_solve_point(ion, 1.6, 0.5, Jnu, struct('J0eff', 6.4e-3));
p2 = invz_solve_point(ion, 1.6, 0.5, Jnu, struct('J0eff', 6.4e-3, 'odd', false));
verifyTrue(testCase, isequaln(p1, p2));
end

function test_solve_point_odd_args_guard(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [1e-3; 2e-3], ...
    struct('odd', true)), 'invz:oddArgs');          % blocks missing
S = struct('Vca', zeros(4,4,2), 'Vcb', zeros(4,4,2), 'Vcc', zeros(4,4,2), 'Jcc0', 6.4e-3);
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [1e-3; 2e-3], ...
    struct('odd', true, 'odd_blocks', S)), 'invz:oddArgs');   % Jnu_flat not []
end

function test_solve_point_odd_pm_point(testCase)
% Hard gate at a GUARANTEED-PM point (1.80 K > no-ODD Tc(0) = 1.743 K): the ODD
% plan's own T1.4 acceptance point (1.60 K, 0.1 T) sits on the ORDERED side of
% the no-ODD baseline, so it cannot be a convergence gate — it is measured and
% REPORTED below instead (whether it now converges is itself the ODD Tc-shift
% signal). Gates: convergence, crit INCREASES with ODD on (sign contract),
% overhead <= 20%.
ion = invz_ion();
n = 8;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 15, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
Jnu0 = zeros(size(Vcc,3), 4);
for iq = 1:size(Vcc,3), Jnu0(iq,:) = sort(real(eig(Vcc(:,:,iq)))).'; end
o0 = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0);
t0 = tic;  p0 = invz_solve_point(ion, 1.80, 0.1, Jnu0(:), o0);            t_off = toc(t0);
o1 = o0;  o1.odd = true;  o1.odd_blocks = S;
t1 = tic;  p1 = invz_solve_point(ion, 1.80, 0.1, [], o1);                 t_on = toc(t1);
verifyTrue(testCase, p0.converged && p1.converged);
verifyGreaterThan(testCase, p1.crit, p0.crit);                            % ODD suppresses ordering
verifyGreaterThan(testCase, p1.odd.d, 0);
verifyLessThan(testCase, p1.sumrule_rel, 0.10);
fprintf('odd overhead: %.2fs vs %.2fs (%.0f%%), crit %.4f -> %.4f, d = %.3g ueV\n', ...
    t_on, t_off, 100*(t_on/t_off - 1), p0.crit, p1.crit, p1.odd.d*1e3);
verifyLessThan(testCase, t_on, 1.2*t_off + 0.5);                          % +0.5 s absolute floor for timer noise
% REPORT the plan's (1.60 K, 0.1 T) point, both flags — no convergence gate:
r0 = invz_solve_point(ion, 1.60, 0.1, Jnu0(:), o0);
r1 = invz_solve_point(ion, 1.60, 0.1, [], o1);
fprintf('plan point (1.60 K, 0.1 T): off converged=%d crit=%.4g; odd converged=%d crit=%.4g\n', ...
    r0.converged, r0.crit, r1.converged, r1.crit);
end

function test_t0field_handles_backcompat(testCase)
% Numeric args byte-identical; constant handles reproduce numerics exactly.
ion = invz_ion();
Tc1 = invz_critical_T0field(ion, 0.30, 6.4e-3);
Tc2 = invz_critical_T0field(ion, @(T) 0.30, @(T) 6.4e-3);
verifyEqual(testCase, Tc2, Tc1);                    % same bisection path
end

function test_critical_T_requires_anchor_with_odd(testCase)
ion = invz_ion();
S = struct('Vca', zeros(4,4,2), 'Vcb', zeros(4,4,2), 'Vcc', zeros(4,4,2), 'Jcc0', 6.4e-3);
verifyError(testCase, @() invz_critical_T(ion, 2.0, [], ...
    struct('odd', true, 'odd_blocks', S)), 'invz:oddTc0');
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement** (read `invz_crit_at.m` first and state in the commit body how opts flow through it; keep every change behind the flag). **Step 4: tests + FULL fast suite green — this task touches the solver spine, so also re-run the two heaviest existing test files (`test_invz_solve_point.m`, `test_invz_critical.m`) explicitly and confirm identical outcomes.**
- [ ] **Step 5: Commit** — `feat(invz): T1.4 ODD wiring — solve_point/_ordered, critical finders, T0field handles, ion.odd driver switch` (flip T1.4 + manifest).

---

### Task 6: Zero-field engine + physics measurements + README (T1.5 + T1.6)

**Files:**
- Create: `invz_projected/invz_odd_zero_field.m`, `invz_projected/tests/test_invz_odd_physics.m`
- Modify: `invz_projected/README.html` (new §1.7 "ODD-mediated coupling (Tier 1)" + §8 scope row), `docs/ODD-LOG.md`

**Interfaces:**
- Consumes: Tasks 2–5.
- Produces: `[Tc, out] = invz_odd_zero_field(ion, opts)` — the T1.5 measurement engine AND the V4.1 zero-field endpoint source:
  - `opts`: `grids` (default `{12, 24}` → per-grid Tc + Richardson `Tc_rich = (n2^2*Tc2 - n1^2*Tc1)/(n2^2 - n1^2)`-style, mirroring the existing Σc benchmark's Richardson convention — READ `test_invz_sigma_crit.m` and reuse ITS formula), `dpRng` (30), `cache` (true), `mode` ∈ `'off' | 'full' | 'uniform_only' | 'qstruct_only'` (default 'full').
  - Per grid n: build a Γ-EXCLUDED uniform mesh exactly the way the existing Σc benchmark builds it (read `test_invz_sigma_crit.m` / `invz_bz_couplings.m` and reuse the same generator + Γ filter so ODD-off reproduces the published number bitwise); `[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(...)`;
    - `'off'`: `Sc(T) = invz_sigma_crit(infoB.Jcc0, Jnu0_flat)` (T-independent), `J0(T) = infoB.Jcc0` — must reproduce the published route;
    - `'full'` (T1.5 a): `Xp(T) = invz_chiperp(ion, T, [0 0 0])`; `[dJ, d] = invz_odd_deltaJ(Vca, Vcb, Xp)`; modes of `Vcc + dJ`; `J0(T) = infoB.Jcc0 - d(T)`; `Sc(T) = invz_sigma_crit(J0(T), Jnu_odd_flat(T))`;
    - `'uniform_only'` (T1.5 b, DS2023's MF mechanism): modes of `Vcc` (unmodified), `J0(T) = infoB.Jcc0 - d(T)`;
    - `'qstruct_only'` (T1.5 c, the fluctuation piece): modes of `Vcc + dJ + d*eye(4)` per q, `J0(T) = infoB.Jcc0`;
    then `Tc = invz_critical_T0field(ion, ScT_handle, J0T_handle)`.
  - `out`: per-grid `Tc`, `Tc_rich`, `Sc_at_Tc`, `Sc_rich`, `d_at_Tc`, `mode`, timings, and (mode 'full') `out.split` = the (a)/(b)/(c) trio + closure defect `(a) - [(b) + (c) - Tc_off]`.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_zero_field_off_matches_published_slow(testCase)
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
[Tc, out] = invz_odd_zero_field(ion, struct('mode', 'off'));
verifyEqual(testCase, Tc, 1.743, 'AbsTol', 0.006);          % published anchor
verifyEqual(testCase, out.Sc_rich, 0.2980, 'AbsTol', 0.006);
end

function test_zero_field_structure_fast(testCase)
% Small-grid structural gate: ODD lowers Tc; each split component lowers Tc;
% (a) is the largest suppression.
ion = invz_ion();
o = struct('grids', {{8}}, 'dpRng', 15);
Tc0 = invz_odd_zero_field(ion, setfield(o, 'mode', 'off'));   %#ok<SFLD>
[Tc1, out1] = invz_odd_zero_field(ion, setfield(o, 'mode', 'full')); %#ok<SFLD>
verifyLessThan(testCase, Tc1, Tc0);
verifyGreaterThan(testCase, out1.d_at_Tc, 0);
verifyLessThan(testCase, out1.split.Tc_b, Tc0);
verifyLessThan(testCase, out1.split.Tc_c, Tc0);
verifyLessThan(testCase, Tc1, min(out1.split.Tc_b, out1.split.Tc_c) + 1e-9);
end

function test_odd_headline_slow(testCase)
% T1.5 headline (REPORT + reproducibility): Richardson(12,24) ODD numbers.
% First run: print and RECORD in ODD-LOG + invz_odd_anchors (controller step);
% the anchor assert below activates once the anchor fields exist.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
[Tc1, out1] = invz_odd_zero_field(ion, struct('mode', 'full'));
fprintf('ODD zero field: Tc = %.4f K (dTc = %.4f), Sc = %.4f (dSc = %+.4f), d(Tc) = %.4g ueV\n', ...
    Tc1, 1.743 - Tc1, out1.Sc_rich, out1.Sc_rich - 0.2980, out1.d_at_Tc*1e3);
fprintf('split: (a) %.4f  (b) %.4f  (c) %.4f K, closure defect %.4g K\n', ...
    Tc1, out1.split.Tc_b, out1.split.Tc_c, out1.split.closure_defect);
A = invz_odd_anchors();
if isfield(A, 'odd_Tc_rich')                      % reproducibility gate (1%, plan T1.5)
    verifyEqual(testCase, Tc1, A.odd_Tc_rich, 'RelTol', 0.01);
    verifyEqual(testCase, out1.d_at_Tc, A.odd_d_at_Tc, 'RelTol', 0.01);
end
verifyTrue(testCase, isfinite(Tc1) && Tc1 > 0.8 && Tc1 < 1.743);
end

function test_boundary_shift_slow(testCase)
% T1.5 boundary table at the README grid: Bc(T) ODD on/off at T = 0.31, 0.8,
% 1.2, 1.5 K; REPORT the shifts + the expected attenuation toward high Bx
% (crossover ~3.5 T per DS2023) — gate only the sign (Bc_on <= Bc_off).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
Jnu0 = zeros(size(Vcc,3), 4);
for iq = 1:size(Vcc,3), Jnu0(iq,:) = sort(real(eig(Vcc(:,:,iq)))).'; end
Ts = [0.31 0.8 1.2 1.5];
for it = 1:numel(Ts)
    o0 = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0);
    b0 = invz_critical(ion, Ts(it), Jnu0(:), o0);
    o1 = o0;  o1.odd = true;  o1.odd_blocks = S;
    b1 = invz_critical(ion, Ts(it), [], o1);
    fprintf('Bc(%.2f K): off %.4f T, on %.4f T, shift %+.4f T\n', Ts(it), b0, b1, b1 - b0);
    verifyLessThanOrEqual(testCase, b1, b0 + 0.02);   % <= within crossing tol
end
end
```

- [ ] **Step 2: Run fast test — verify failure.** **Step 3: implement `invz_odd_zero_field.m`** (grid/Γ-filter/Richardson conventions copied from the existing Σc benchmark — name the source lines in a comment). **Step 4: fast suite green; run the slow tests** (`INVZ_SLOW=1` on this file; budget: Bc table ≈ 4 × 2 finder runs on warm caches — expect tens of minutes; run in background and continue only when green).
- [ ] **Step 5: Pin the measured headline numbers** as `A.odd_Tc_rich`, `A.odd_d_at_Tc` (+ `A.odd_Sc_rich`) in `invz_odd_anchors.m`; append the full T1.5 block to ODD-LOG: ΔTc(0), ΔΣc(0), (a)/(b)/(c) + closure defect, Bc-shift table + attenuation pattern vs the ~3.5 T crossover, d-band check (0.3–0.5 μeV), DS2023 qualitative anchors (their 5% 3-state MF reduction; MC Tc = 1.6295 K) with the models-differ caveat (plan §8).
- [ ] **Step 6 (T1.6): README** — new §1.7 with (E1)–(E5), the flag surface (`opts.odd`/`odd_blocks`, `ion.odd`), the χ⊥-held-fixed design decision, measured ΔTc; add `invz_odd_blocks`, `invz_chiperp`, `invz_odd_deltaJ`, `invz_odd_zero_field` to the function reference; §8 scope: move "ODD terms" to in-scope-(Gaussian) with the Tier-3 caveat; note Σc closed form is now T-dependent with ODD on.
- [ ] **Step 7: Commit** — `feat(invz): T1.5+T1.6 zero-field ODD measurements (MF/fluctuation split), boundary shifts, README SS1.7` (flip T1.5, T1.6 + manifest).

---

### Task 7: Tier 1b — retardation sensitivity (T2.1 + T2.2)

**Files:**
- Modify: `invz_projected/invz_emt_scalar.m` (Jnu_flat `[nJ,1]` → also `[nJ,nw]`), `invz_projected/invz_solve_point.m` (`opts.odd_retarded`)
- Create: `invz_projected/tests/test_invz_odd_retarded.m`

**Interfaces:**
- Produces:
  - `invz_emt_scalar(G0, Sigma, Jnu_flat, opts)`: `Jnu_flat` may be `[nJ, nw]` (per-frequency mode spectra). The block loop generalizes: for matrix input, column `n` pairs with frequency `n` (`A(idx) = mean(Jf(:,idx) ./ (D(idx).' + Jf(:,idx) .* G0(idx).'), 1).'`). `[nJ,1]` input takes the EXISTING code path byte-identical (regression-gated).
  - `invz_solve_point` `opts.odd_retarded = true` (requires `opts.odd`): after the static rebuild, evaluate `Xpn = invz_chiperp(ion, T, Bx, struct(..., 'z', 1i*wn))`; scalar surrogate `r_n = squeeze(Xpn(1,1,:)) / Xpn(1,1,1)` (assert `r_1 = 1`, `0 < r_n <= 1`, monotone; log `max |Xpn(1,2,:)|/Xpn(1,1,:)` — the χab smallness that justifies the scalar surrogate); per-frequency modes by FIRST-ORDER PERTURBATION around the static-ODD spectrum: `Jnu_n(q,nu) = Jnu_static(q,nu) + (r_n - 1) * w(q,nu)` with `w(q,nu) = u_nu' * dJ(q) * u_nu` (eigvec projections captured during the static rebuild — no per-n eig); assemble `Jnu_flat [nq*4, nwn]` and pass to the generalized EMT. `pt.odd.r_n` stored. Exact-diagonalization cross-check per plan: behind `opts.odd_retarded_exact = true`, do the full per-(q,n) eig — used ONCE by a test to log the surrogate+perturbation error, not in production.
- Acceptance (T2.1): with `r_n ≡ 1` forced, bitwise agreement with the static T1 output; sum-rule residual unchanged at the 1e-6 level. Decision (T2.2): report Tc(0.1 T) via `invz_critical_T` and Bc(1.2 K) via `invz_critical`, static vs retarded; if |ΔTc| < 10 mK freeze STATIC as default (retarded stays behind the flag), else retarded becomes the recommended default and T1.5 is re-measured — either way ODD-LOG records both and the README documents the decision. (Note: the zero-field closed form is intrinsically static — criticality is controlled by the elastic n = 0 sector where χ⊥(0) is exact, plan §4 mitigation — so the T2.2 Tc comparison runs at the 0.1 T proxy through the full solver.)

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_emt_matrix_column_constant_bitwise(testCase)
% [nJ,1] path untouched; [nJ,nw] with identical columns reproduces it exactly.
G0 = -linspace(30, 5, 12).';  Sigma = 0.05*ones(12,1);
Jf = linspace(-2e-3, 6e-3, 24).';
m1 = invz_emt_scalar(G0, Sigma, Jf, struct('debug', true));
m2 = invz_emt_scalar(G0, Sigma, repmat(Jf, 1, 12), struct('debug', true));
verifyEqual(testCase, m2.K, m1.K);           % bitwise
verifyEqual(testCase, m2.G, m1.G);
end

function test_retarded_rn_unity_bitwise(testCase)
% r_n forced to 1 (test hook opts.odd_rn_override) == static ODD solve, bitwise.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
p1 = invz_solve_point(ion, 1.8, 0.5, [], o);
o2 = o;  o2.odd_retarded = true;  o2.odd_rn_override = 1;    % force r_n = 1
p2 = invz_solve_point(ion, 1.8, 0.5, [], o2);
verifyEqual(testCase, p2.Sigma, p1.Sigma);   % bitwise (same EMT inputs)
verifyEqual(testCase, p2.crit,  p1.crit);
end

function test_retarded_physics_and_surrogate_error(testCase)
% Retarded weakens the n~=0 mediated coupling: Sigma0 decreases vs static.
% Also logs the scalar-surrogate + first-order-perturbation error vs exact eig.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
ps = invz_solve_point(ion, 1.8, 0.5, [], o);
o.odd_retarded = true;
pr = invz_solve_point(ion, 1.8, 0.5, [], o);
verifyTrue(testCase, ps.converged && pr.converged);
verifyLessThanOrEqual(testCase, pr.Sigma0, ps.Sigma0 + 1e-12);
verifyTrue(testCase, all(pr.odd.r_n > 0 & pr.odd.r_n <= 1 + 1e-12));
verifyLessThan(testCase, pr.sumrule_rel, ps.sumrule_rel + 1e-6);   % T2.1 acceptance
o.odd_retarded_exact = true;
pe = invz_solve_point(ion, 1.8, 0.5, [], o);
fprintf('retarded: Sigma0 static %.5f -> ret %.5f (exact %.5f); surrogate err %.2e\n', ...
    ps.Sigma0, pr.Sigma0, pe.Sigma0, abs(pr.Sigma0 - pe.Sigma0));
verifyEqual(testCase, pr.Sigma0, pe.Sigma0, 'AbsTol', 1e-3);       % surrogate sanity
end

function test_t2_decision_measurement_slow(testCase)
% T2.2: Tc(0.1 T) and Bc(1.2 K), static vs retarded. REPORT + the 10 mK rule.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
[TcA, outA] = invz_odd_zero_field(ion, struct('mode', 'full'));    % static anchor for Tc0
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S, 'Tc0', TcA);
tS = invz_critical_T(ion, 0.1, [], o);
o.odd_retarded = true;
tR = invz_critical_T(ion, 0.1, [], o);
bS = invz_critical(ion, 1.2, [], rmfield(o, 'odd_retarded'));
bR = invz_critical(ion, 1.2, [], o);
fprintf('T2.2: Tc(0.1T) static %.4f vs retarded %.4f (d = %+.1f mK); Bc(1.2K) %.4f vs %.4f T\n', ...
    tS, tR, 1e3*(tR - tS), bS, bR);
verifyTrue(testCase, all(isfinite([tS tR bS bR])));
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement** (EMT generalization FIRST with its bitwise gate, then the solver branch). **Step 4: fast suite green; run the slow decision test; apply the 10 mK decision rule and document it** (README §1.7 note + ODD-LOG; if retarded wins, re-run Task 6's headline slow test and re-pin anchors — say so explicitly in the log).
- [ ] **Step 5: Commit** — `feat(invz): T2.1+T2.2 retarded mediated coupling (scalar r_n surrogate) + static-vs-retarded decision` (flip T2.1, T2.2 + manifest).

---

### Task 8: `invz_odd_fieldvar.m` — Tier-2 field covariance (T3.1)

**Files:**
- Create: `invz_projected/invz_odd_fieldvar.m`, `invz_projected/tests/test_invz_odd_fieldvar.m`

**Interfaces:**
- Consumes: converged `pt` from `invz_solve_point` (needs `pt.G`, `pt.K`, `pt.Sigma`, `pt.odd`), blocks struct `S`, `invz_matsubara`.
- Produces: `[C, info] = invz_odd_fieldvar(ion, pt, S, T, opts)` — E3:
  `C_ab = (1/(4*Nq)) * sum_q tr[ V_ac(q) * Scc(q) * V_cb(q) ]` with `V_ac = Vca'` per q, `Scc(q) [4,4] = (1/beta) sum_n wts_n * chi_cc(q, i w_n)`;
  `chi_cc(q, i w_n) = -sum_nu u_nu(q) u_nu(q)' * Gq_nu(i w_n)`, `Gq_nu = G ./ (1 + (Jnu_odd(q,nu) - K).*G)` (the converged EMT lattice propagator, framework eq 7 — G = pt.G, K = pt.K), `u_nu`/`Jnu_odd` from `eig(Vcc + dJ)` recomputed inside (4×4 per q, cheap; `dJ` from `pt.odd.Xp` so the covariance is evaluated at the SAME converged state).
  Matsubara sum: reuse `[wn, wts, beta] = invz_matsubara(T, Ecut)` with the pt's Ecut (pass via opts, default 40) — same truncation the sum rule carries; `info.tail_share` = |last-frequency contribution|/|total| (report). `opts.static_approx = true` gates the cheap `Scc ≈ kB*T*chi_cc(q, 0)` variant — NEVER silent (plan T3.1): when used, `info.static_approx = true` and the caller logs the difference once (test below).
  `C` real symmetric 2×2 PSD (meV²); `info.heq_T = sqrt(max(eig(C))) / (ion.gL * C0.muB)` — the equivalent transverse field scale in Tesla for the qualitative Dollberg-width comparison (h ≈ 0.4 T, their Fig. 4), log-only.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_fieldvar_structure(testCase)
ion = invz_ion();  T = 1.8;   % PM-guaranteed (> no-ODD Tc(0) = 1.743 K)
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
pt = invz_solve_point(ion, T, 0.5, [], o);
assumeTrue(testCase, pt.converged);
[C, info] = invz_odd_fieldvar(ion, pt, S, T, struct());
verifyEqual(testCase, C, C.', 'AbsTol', 1e-15*max(abs(C(:)) + 1e-30));
e = eig(C);
verifyGreaterThan(testCase, min(e), -1e-14*max(e));       % PSD
verifyGreaterThan(testCase, max(e), 0);                    % nonzero variance
verifyLessThan(testCase, info.tail_share, 0.05);           % w_n^-2 tail under control
fprintf('T3.1: C_aa = %.4g meV^2, heq = %.3f T (Dollberg ~0.4 T), tail %.2e\n', ...
    C(1,1), info.heq_T, info.tail_share);
end

function test_fieldvar_aa_bb_at_zero_Bx(testCase)
% C_aa = C_bb at Bx -> 0 (C4). Small field + PM-guaranteed temperature.
ion = invz_ion();  T = 1.8;
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
pt = invz_solve_point(ion, T, 0.05, [], o);
assumeTrue(testCase, pt.converged);
C = invz_odd_fieldvar(ion, pt, S, T, struct());
verifyEqual(testCase, C(1,1), C(2,2), 'RelTol', 5e-2);
end

function test_fieldvar_static_approx_logged(testCase)
% The kB*T*chi(0) shortcut is gated and its error is measured once (plan T3.1).
ion = invz_ion();  T = 1.8;
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
pt = invz_solve_point(ion, T, 0.5, [], o);
assumeTrue(testCase, pt.converged);
[Cf, ~]  = invz_odd_fieldvar(ion, pt, S, T, struct());
[Cs, is] = invz_odd_fieldvar(ion, pt, S, T, struct('static_approx', true));
verifyTrue(testCase, is.static_approx);
fprintf('T3.1 static-approx rel diff: %.3f\n', norm(Cs - Cf)/norm(Cf));
verifyTrue(testCase, isfinite(norm(Cs - Cf)/norm(Cf)));
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement** (basis choice: contract in the MODE basis with eigenvectors — `tr[Vac*Scc*Vcb] = sum_nu (u' Vcb ... )` — or rotate to sublattice basis; pick whichever is cheaper given the 4×4 eig is recomputed anyway, and DOCUMENT the choice in the header per plan T3.1). **Step 4: tests + full fast suite green;** log C, heq_T, tail_share, static-approx delta to ODD-LOG.
- [ ] **Step 5: Commit** — `feat(invz): T3.1 ODD internal-field covariance (E3)` (flip T3.1 + manifest).

---

### Task 9: `invz_twolevel_avg.m` — Gaussian-dressed doublet (T3.2)

**Files:**
- Create: `invz_projected/invz_twolevel_avg.m`, `invz_projected/tests/test_invz_twolevel_avg.m`
- Reference: `invz_common` does not exist — `invz_projected/invz_twolevel.m`, `invz_projected/invz_g.m`, `jensen_1z_framework.html` §7 (two-level sum rule)

**Interfaces:**
- Consumes: `invz_twolevel` (fields `Delta`, `M2`, `m`, `n01`, `g0`, `transverse_mf`; guards `invz:degenerateDoublet`, `invz:orderedPhase`), `invz_g` (`g(z) = 2*n01*Delta ./ (Delta^2 - z.^2)`), `invz_single_ion`, `invz_chi0z`, `invz_const`.
- Produces: `[tla, avg] = invz_twolevel_avg(ion, T, Bx, C, opts)`:
  - `C = zeros(2)` (or all-zero) SHORT-CIRCUITS: `tla = invz_twolevel(ion, T, Bx, opts_pass)` — bitwise, no quadrature (`avg.nodes = 1`).
  - Else: 2-D Gauss–Hermite, `opts.ngh` (default 7) nodes per axis; diagonalize `C = V*diag(s2)*V'`; node fields `[ha; hb] = V*[sqrt(2*s2(1))*x_i; sqrt(2*s2(2))*x_j]` (meV, conjugate to J_a/J_b), weights `w_ij = (wi*wj)/pi` (normalized: sum = 1 — assert to 1e-12). Each node enters the single-ion problem as a FIELD SHIFT `B_node = invz_field_vec(Bx) + [ha, hb, 0]/(ion.gL*C0.muB)` — injects exactly `-ha*Ja - hb*Jb` through the existing Zeeman term; NO single-ion code changes.
  - DEFAULT `opts.avg = 'response'` (the plan's design decision): per node build `tl_node = invz_twolevel(ion, T, B_node, opts_pass)` and average `gbar(i w_n) = sum_nodes w * invz_g(tl_node, i w_n)` on the grid `opts.wn` (REQUIRED for response mode; from `invz_matsubara`); fit effective params by matching (i) `gbar(0)`, (ii) `(1/beta)*sum_n wts.*gbar` — equals 1 EXACTLY per node (two-level identity n01 = tanh(beta*Delta/2): verify analytically in a comment; assert the truncated-grid value is within 2e-2 of 1 per node and that averaging moves it by < 1e-6 — the plan's "check it survives averaging"), (iii) the leading `w_n^-2` tail coefficient `lim w^2*gbar = 2*n01*Delta`; solve `Delta_eff = sqrt(tail/gbar0)`, `n01_eff = gbar0*Delta_eff/2`, then `M2_eff` from node-averaged `M2*g` weight: `M2_eff = (sum_nodes w * tl.M2 * g_node(0)) / gbar0` (moment-weighted; document). `tla` carries the SAME field names as `invz_twolevel` (`Delta`→Delta_eff etc., plus `tla.gh = true`); `avg`: `Delta_eff`, `M2_eff`, `n01_eff`, `fit_resid` (max relative mismatch of the three conditions when reproduced by the fitted two-level form), `node_Delta` spread, `sumrule_avg`.
  - `opts.avg = 'params'`: average `(Delta, M2, n01)` directly (the inequivalent variant — kept for the one-shot comparison, plan T3.2).
  - `opts.wn` also triggers `avg.G0` = node-averaged FULL electronuclear cc propagator (`-real(squeeze(chi0(3,3,:)))` per node via `invz_single_ion(hyp=true)` + `invz_chi0z(elastic=true)`) and `avg.chi0cc0` — the disorder-averaged G0 the Tier-2 outer loop swaps in (controller-approved extension of the plan's T3.2 spec: the plan's "average ... χ0cc node-by-node" needs a home; this is it).
  - Degenerate-doublet handling (T3.4 groundwork): nodes where `invz_twolevel` would throw `invz:degenerateDoublet` (Δ_node < 1e-4 meV, e.g. the exact-origin node at Bx = 0) are evaluated in their h→0 LIMIT: catch the error, substitute the analytic degenerate-limit response `g_node(i w_n) = beta * (w_n == 0)`-type elastic weight — NO: instead use the limit of the two-level form: `g(0) = beta` (since `2*tanh(beta*D/2)/D -> beta` as D → 0) and `g(i w_n ~= 0) -> 2*D*n01/w_n^2 -> 0`; document this branch and count such nodes in `avg.n_degenerate`. (At Bx = 0 the OFF-origin nodes carry finite Δ_node ∝ Van-Vleck·h² — the internal field lifts the degeneracy node-wise, plan T3.4.)

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_C0_bitwise(testCase)
ion = invz_ion();
tl  = invz_twolevel(ion, 1.6, 2, struct());
tla = invz_twolevel_avg(ion, 1.6, 2, zeros(2), struct());
fn = fieldnames(tl);
for i = 1:numel(fn), verifyEqual(testCase, tla.(fn{i}), tl.(fn{i})); end   % bitwise per shared field
end

function test_gh_weights_and_monotonicity(testCase)
% Level repulsion: Delta_eff grows with ||C||; moments shrink. Realistic C scale:
% Dollberg h ~ 0.4 T <-> C ~ (0.4 * gL*muB)^2 ~ 8e-4 meV^2 — ray spans below that.
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
s = [1e-5 8e-5 6.4e-4];
[D, M2] = deal(zeros(3,1));
for i = 1:3
    [~, avg] = invz_twolevel_avg(ion, 1.6, 2, s(i)*eye(2), struct('wn', wn));
    D(i) = avg.Delta_eff;  M2(i) = avg.M2_eff;
    verifyLessThan(testCase, avg.fit_resid, 1e-6);
end
tl = invz_twolevel(ion, 1.6, 2, struct());
verifyTrue(testCase, all(diff(D) > 0));
verifyGreaterThanOrEqual(testCase, D(1), tl.Delta - 1e-12);   % Delta_eff >= Delta
verifyTrue(testCase, all(diff(M2) < 0));
verifyLessThanOrEqual(testCase, M2(1), tl.M2 + 1e-12);        % M2_eff <= M2
end

function test_sumrule_survives_averaging(testCase)
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
[~, avg] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn));
verifyEqual(testCase, avg.sumrule_avg, 1, 'AbsTol', 2e-2);    % truncation-level, not drift
end

function test_zero_field_origin_node_limit(testCase)
% T3.4: B = 0 with C > 0 must not throw invz:degenerateDoublet; origin node
% handled by its h -> 0 limit, off-origin nodes lift the degeneracy.
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.55, 40);
[tla, avg] = invz_twolevel_avg(ion, 1.55, 0, 2e-4*eye(2), struct('wn', wn));
verifyTrue(testCase, isfinite(tla.Delta) && tla.Delta > 0);
verifyGreaterThanOrEqual(testCase, avg.n_degenerate, 1);      % the origin node
fprintf('T3.4: B = 0 dressed doublet Delta_eff = %.4g meV (bare degenerate)\n', tla.Delta);
end

function test_avg_mode_comparison_reported(testCase)
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
[~, ir] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn));
[~, ip] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn, 'avg', 'params'));
fprintf('T3.2 avg-mode Delta_eff: response %.6g vs params %.6g meV\n', ir.Delta_eff, ip.Delta_eff);
verifyTrue(testCase, isfinite(ir.Delta_eff) && isfinite(ip.Delta_eff));
end

function test_gh_node_convergence(testCase)
% 5 vs 7 vs 9 nodes at a realistic C: Delta_eff differences shrink (V4.2 sweep seed).
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.6, 40);
D = zeros(3,1);  ns = [5 7 9];
for i = 1:3
    [~, avg] = invz_twolevel_avg(ion, 1.6, 2, 2e-4*eye(2), struct('wn', wn, 'ngh', ns(i)));
    D(i) = avg.Delta_eff;
end
verifyLessThan(testCase, abs(D(3) - D(2)), abs(D(2) - D(1)) + 1e-15);
fprintf('GH 5/7/9: Delta_eff = %.8g / %.8g / %.8g meV\n', D);
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement** (GH nodes/weights via a small local function — hardcode the standard recurrence (Golub–Welsch on the Hermite Jacobi matrix, `eig` of the symmetric tridiagonal), no toolbox dependency; FLAG IN HEADER verbatim: this quenched-Gaussian, static dressing is the least rigorous step of the whole plan — README must say so too (Task 10 wires that)). **Step 4: tests + full fast suite green;** ODD-LOG entry with the avg-mode and GH-convergence numbers.
- [ ] **Step 5: Commit** — `feat(invz): T3.2 Gauss-Hermite-dressed doublet invz_twolevel_avg (response-averaged, params variant gated)` (flip T3.2 + manifest).

---

### Task 10: Tier-2 outer self-consistency + small-Bx characterization (T3.3 + T3.4)

**Files:**
- Modify: `invz_projected/invz_solve_point.m` (`opts.odd_tier2` outer loop), `invz_projected/README.html` (§1.7 Tier-2 subsection + "least rigorous step" flag)
- Create: `invz_projected/tests/test_invz_odd_tier2.m`

**Interfaces:**
- Produces: `invz_solve_point` `opts.odd_tier2 = true` (requires `opts.odd`): after the Tier-1 solve converges — outer₂ loop (cap `opts.max_tier2` = 8, mix 0.5 on C, tol `opts.tol_tier2` = 1e-3 relative on `max|dC|`): `C = invz_odd_fieldvar(ion, pt, S, T, ...)` → `[tla, avg] = invz_twolevel_avg(ion, T, Bx, C_mixed, struct('wn', wn, ...))` → re-run the inner Σ loop with `tl = tla`, `g = invz_g(tla, 1i*wn)`, `G0 = avg.G0` (the disorder-averaged propagator) → recompute C; converge. χcc feedback is stabilizing (suppression reduces Scc — expect 2–4 iterations, plan T3.3); non-convergence masked exactly like the EMT loop (`pt.converged = false`). New pt fields: `pt.C`, `pt.tier2_iters`, `pt.tla`. Flag off: byte-identical to Task 7's output (guarded).
- IR safety (T3.3, strong correctness check): C remains bounded as crit → 0 at fixed T (the ODD blocks vanish at q = 0) — slow test approaches the boundary.
- T3.4 deliverable is CHARACTERIZATION (in ODD-LOG + README): the Bx → 0 behavior of the Tier-1+2 point (via the 0.05 T proxy per plan T3.3 and the `invz_twolevel_avg` B = 0 handling from Task 9), the `invz:degenerateDoublet` guard status (deliberate, not accidental), and whether zero-field results move from "directional only" to quantitative or the hyperfine-scale caveat stands.

- [ ] **Step 1: Write the failing tests**:

```matlab
function test_tier2_flag_off_byte_identical(testCase)
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
p1 = invz_solve_point(ion, 1.6, 0.5, [], o);
o2 = o;  o2.odd_tier2 = false;
p2 = invz_solve_point(ion, 1.6, 0.5, [], o2);
verifyTrue(testCase, isequaln(p1, p2));
end

function test_tier2_converges_and_suppresses(testCase)
% T3.3 convergence gate at a GUARANTEED-PM point (1.80 K, 0.05 T): the plan's
% (1.55 K, 0.05 T proxy) sits below the no-ODD Tc and is REPORTED separately
% below (whether it converges tells which side of the ODD-shifted boundary it
% lands on). Variable moments suppress ordering: crit increases vs Tier 1.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
p1 = invz_solve_point(ion, 1.80, 0.05, [], o);
o.odd_tier2 = true;
p2 = invz_solve_point(ion, 1.80, 0.05, [], o);
verifyTrue(testCase, p1.converged && p2.converged);
verifyLessThanOrEqual(testCase, p2.tier2_iters, 8);
verifyGreaterThan(testCase, p2.crit, p1.crit - 1e-12);     % suppression direction
verifyGreaterThanOrEqual(testCase, p2.tla.Delta, p1.tl.Delta - 1e-12);  % level repulsion
fprintf('T3.3: tier2 iters = %d, crit %.5f -> %.5f, Delta %.5g -> %.5g, C_aa = %.3g meV^2\n', ...
    p2.tier2_iters, p1.crit, p2.crit, p1.tl.Delta, p2.tla.Delta, p2.C(1,1));
% REPORT the plan's T3.3 point (1.55 K, 0.05 T), no gate:
r = invz_solve_point(ion, 1.55, 0.05, [], o);
fprintf('plan point (1.55 K, 0.05 T) Tier1+2: converged=%d, crit=%.4g\n', r.converged, r.crit);
end

function test_tier2_C_bounded_near_boundary_slow(testCase)
% T3.3 IR safety: C stays finite as crit -> 0 (ODD blocks vanish at q = 0).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();  T = 1.2;
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
Bc = invz_critical(ion, T, [], o);
Cn = zeros(3,1);  db = [0.5 0.2 0.05];
for i = 1:3
    pt = invz_solve_point(ion, T, Bc + db(i), [], o);
    assumeTrue(testCase, pt.converged);
    C = invz_odd_fieldvar(ion, pt, S, T, struct());
    Cn(i) = norm(C);
end
verifyLessThan(testCase, Cn(3)/Cn(1), 20);                  % grows but saturates
fprintf('T3.3 IR: |C| at Bc+0.5/0.2/0.05 T = %.3g / %.3g / %.3g meV^2\n', Cn);
end

function test_tier2_combined_measurement_slow(testCase)
% Combined dTc and the Tier1 : Tier2 split (REPORT). Tc via invz_critical_T at
% Bx = 0.5 T — NOT lower: below ~0.5 T the doublet is near-degenerate and the
% boundary suffers the known small-B non-convergence speckle (plan SS8 +
% invz_critical_T header); 0.5 T is the established reliable floor.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
Jnu0 = zeros(size(Vcc,3), 4);
for iq = 1:size(Vcc,3), Jnu0(iq,:) = sort(real(eig(Vcc(:,:,iq)))).'; end
Tc0f = invz_odd_zero_field(ion, struct('mode', 'off'));
o0 = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'Tc0', Tc0f);
t_off = invz_critical_T(ion, 0.5, Jnu0(:), o0);
o1 = o0;  o1.odd = true;  o1.odd_blocks = S;  o1.Tc0 = invz_odd_zero_field(ion, struct('mode', 'full'));
t_t1 = invz_critical_T(ion, 0.5, [], o1);
o2 = o1;  o2.odd_tier2 = true;
t_t12 = invz_critical_T(ion, 0.5, [], o2);
fprintf('T3.3 combined (0.5 T): Tc off %.4f, +T1 %.4f, +T1+T2 %.4f K; split %.1f%% : %.1f%%\n', ...
    t_off, t_t1, t_t12, 100*(t_off - t_t1)/max(t_off - t_t12, 1e-9), ...
    100*(t_t1 - t_t12)/max(t_off - t_t12, 1e-9));
verifyLessThanOrEqual(testCase, t_t1, t_off + 5e-3);
verifyLessThanOrEqual(testCase, t_t12, t_t1 + 5e-3);
end
```

- [ ] **Step 2: Run — verify failure.** **Step 3: implement the outer loop** (all behind the flag; the Σ inner loop is re-entered with swapped `tl`/`g`/`G0` — factor the inner loop into a local function if that keeps the flag-off path textually untouched). **Step 4: tests + full fast suite green; run the slow tests (background; tens of minutes on warm caches).**
- [ ] **Step 5: ODD-LOG + README**: combined ΔTc and Tier1 : Tier2 split; T3.4 characterization paragraph (Bx → 0 status, guard behavior deliberate, quantitative-vs-directional verdict); README §1.7 Tier-2 subsection with the "least rigorous step" flag and the response-vs-params inequivalence note.
- [ ] **Step 6: Commit** — `feat(invz): T3.3+T3.4 Tier-2 outer self-consistency, IR-safety gate, small-Bx characterization` (flip T3.3, T3.4 + manifest).

---

### Task 11: V4 — headline figure, robustness sweeps, consolidation, handoff (V4.1–V4.3)

**Files:**
- Modify: `invz_projected/invz_run_phase_diagram.m` (ODD overlay flag), `docs/ODD-LOG.md`, `invz_projected/tests/invz_odd_anchors.m` (final pins), `invz_projected/README.html` (§7 note that ODD benchmarks exist)
- Create: `docs/SESSION-2026-07-16-invz-odd-mainbody.md`

**Interfaces:**
- Produces:
  - `invz_run_phase_diagram` ODD mode (driven by `ion.odd` from Task 5 + a `quick` config): overlays 1/z baseline, 1/z+Tier1, 1/z+Tier1+2, MF, and the experimental points already used for README §7 (Bitko 1996 / Babkevich 2016 as in R2007 Fig. 1); zero-field endpoints from `invz_odd_zero_field` (closed-form route); second panel Σ(0) along the boundary with/without ODD. `quick` mode (coarse grid, ~3 boundary points + endpoints) runs in-session; the production sweep (hours) is documented and LEFT TO THE USER (repo precedent). Save the quick figure as `Data/Phase_ODD_overlay_quick.fig`.
  - Robustness sweeps (V4.2, each ONE logged point, some INVZ_SLOW): grid 12³/16³/24³ (+ Richardson) for ΔTc(0) — note whether ODD convergence in grid size is slower than the baseline Σc's (δJ shifts weight to finite q); dpRng 20/30/40 sensitivity of d and of the smallest-shell ODD blocks (confirm the effective r⁻⁶ short-rangedness after the χ⊥ contraction of two r⁻³ kernels); GH 5/7/9 on one Tier-2 point (Task 9's test already prints it — copy the numbers). Every headline number in the final table carries an error bar from these sweeps.
  - Consolidation (V4.3): confirm fast-suite ODD additions total < 30 s (time them; if over, demote the heaviest to INVZ_SLOW and say so); full fast + slow suites green; flag-off path re-verified against the Task 1 frozen baseline (same totals, same benchmark digits).
  - Handoff: `docs/SESSION-2026-07-16-invz-odd-mainbody.md` — module map of every new function/flag/cache key, the T2.2 decision, the T3.4 verdict, open Tier-3 items (plan §7: matrix EMT, lost single-channel identities, three-level minimal sector, CF-state convergence, Tier1+2 ↔ Tier3 cross-validation target) and the deferred Appendix-A/`invz_tensor` route pointer.
  - Final ODD-LOG entry: Tc(0) [1/z 1.74 K → 1/z+ODD x K vs exp 1.53 K], the fraction of the 0.21 K gap closed, ΔΣc, the (a)/(b)/(c) split, Bc-shift table with attenuation crossover, d ± error bar, timings.

- [ ] **Step 1: Extend the driver** (read its current structure first; keep flag-off byte-identical; pre-parfor block build per P0.4). Run `quick` mode; save the figure; log the numbers.
- [ ] **Step 2: Run the V4.2 sweeps** (background where slow); collect the error-bar table into ODD-LOG.
- [ ] **Step 3: Consolidation** — run both suites; verify baseline parity; adjust fast/slow gating if the 30 s budget is exceeded.
- [ ] **Step 4: Write the handoff doc + final ODD-LOG + README §7 note; pin any remaining anchors.**
- [ ] **Step 5: Commit** — `docs(invz): V4.1-V4.3 ODD overlay driver, robustness error bars, consolidation, handoff` (flip V4.1–V4.3 + manifest — the plan file should now show every main-body task `[x]`).

---

## Execution notes for the controller

- Tasks strictly sequential. Dispatch models: Task 1 standard; Tasks 2, 3 standard; Tasks 4, 5, 7, 8, 9, 10 physics-critical (most capable tier for implementation AND review); Tasks 6, 11 standard implementation but most-capable review of the physics numbers.
- Slow-test budget: Tasks 6, 7, 10 carry multi-ten-minute INVZ_SLOW runs (Bc tables) — run them in the background via the controller, not inside a blocking implementer turn, when practical.
- The ODD plan's §0 stop rule binds every task: a P0/T1 finding that contradicts §1 (units, χ⊥ band, d band, symmetry) is a BLOCKED report to the user, never an improvised fix.
- After Task 6 and after Task 7's decision rule, re-read the ODD-LOG numbers before dispatching the next task — they parameterize later briefs (anchors, expected scales).
