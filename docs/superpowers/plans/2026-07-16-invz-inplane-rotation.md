# In-Plane Field Rotation (phi_ab) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a C4-consistent in-plane field rotation (`phi_ab`) to the scalar-cc 1/z module via an opt-in vector transverse mean field `(hx, hy)`, with diagnostics, a Sigma=0 tensor reference over angle, and a pinned mapping to the experimentally calibrated external-stack convention (`ion.cfRot(Ho) = -11 deg`).

**Architecture:** Staged per the reviewed Codex plan (`IP-field-rotation_plan_byCodex.md`) with four controller amendments: (1) extend the existing `invz_chi_tensor_ref`/`invz_run_tensor_ref` machinery instead of rebuilding it; (2) library-level reject-with-error (never auto-switch) when a b-axis field component meets `legacy_x`; (3) the combined `theta_c + phi_ab` case is explicitly out of the validated scope; (4) scalar `Jaa0` on the existing `Jxx0` plumbing, no `Jperp0` matrix (YAGNI — tetragonal symmetry forces Jaa = Jbb).

**Tech Stack:** MATLAB R2025a, `functiontests(localfunctions)` test suites, existing `invz_projected/` module.

## Global Constraints

- Module root: `invz_projected/` under `/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion/`. The repo path contains spaces — ALWAYS quote it in shell commands.
- Test command (fast suite): `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests'); assertSuccess(results)"` run from the repo root. Baseline before this plan: 81 passed / 0 failed / 12 incomplete (incomplete = INVZ_SLOW-gated + src/ cross-check skips; they are expected). Every task ends with this suite green.
- Backward compatibility is bit-for-bit: `transverse_mf` defaults to `'legacy_x'` in EVERY function; a default-opts call must execute the identical code path (hy stays exactly 0.0; `H0 - 0*Jy` is bitwise `H0`).
- Mode vocabulary everywhere: `'legacy_x'` (current x-only MF) | `'none'` (hx = hy = 0, bare CF+Zeeman diagnostic) | `'vector_ab'` (self-consistent hx AND hy, both always active, never guarded on By). Invalid value → `error('invz:transverseMF', ...)`.
- Vector mode coupling: `hy_new = Jxx0*jy` using the SAME scalar `Jxx0` channel as hx (drivers pass demag-aware `info.Jaa0` through the existing `opts.Jxx0` plumbing; `ion.Jxx0` = fallback). No matrix `Jperp0`.
- `solve_opts` reserved fields stay exactly `J0eff/Jxx0/hyp` (error `invz:solveOpts`); `transverse_mf` is a LEGAL `solve_opts` field.
- Fail-loud guard (amendment 2): `invz_spectra_map` / `invz_spectra_qpath` MUST error (`invz:transverseMF`) when the field table has any nonzero b-component AND the mode is `'legacy_x'`. `'none'` and `'vector_ab'` pass.
- Fast tests use the synthetic fixture `Jnu = linspace(-2e-3, 6.0e-3, 24).'`, `info = struct('Jcc0', 6.4e-3)`. NEVER use fixture fields with B in [1, 2] T at T = 0.31 K (known non-convergence island). Slow/production-coupling tests gate on `assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), ...)`.
- Digit anchors (controller-verified on this codebase, 2026-07-16, `ion = invz_ion()`, T = 0.31 K, `hyp = false`):
  - legacy `[4 0 0]`: `E(2)-E(1) = 0.369235620278` meV, `Jexp(2) = -0.0689991117463`, `Jexp(1) = 3.5842413528471`
  - legacy `[0 4 0]`: `E(2)-E(1) = 0.352382559983` meV
  - bare (`Jxx0 = 0`) at 4 T: `Delta(15 deg) = 0.345693898281`, `Delta(60 deg) = 0.370662538778` meV
- Production context: the external MF/RPA stack (`Mean Field/LiReF4`) matches experiment with `ion.cfRot(Ho) = -11 deg` (in-plane crystal-field rotation, `cf.m` 'coefficient' method, 4*rot mixing of the m=4 Stevens pairs). Task 5 pins the equivalent `phi_ab` sign numerically. `phi_ab ~ 11 deg` is the production target angle and appears in every validation angle set.
- Physics facts usable in tests: C4 is EXACT for the crystal field (all non-axial terms are J±^4); `<Jz> = 0` exactly for in-plane fields on the para path (time-reversal x C2); B64s breaks only the sigma_v mirrors (extrema need not sit at 0/45 deg).
- Scope guard (amendment 3): nothing in this plan validates combined `theta_c ~= 0 AND phi_ab ~= 0`; the tilt bound theta_c <= 5 deg was measured under `legacy_x` and does not transfer. Docs must say so explicitly.
- Commit messages: conventional prefix + `(invz)` scope, e.g. `feat(invz): ...`, each ending with the Claude Code trailer.

---

### Task 1: `transverse_mf` modes in `invz_single_ion`

**Files:**
- Modify: `invz_projected/invz_single_ion.m`
- Test: `invz_projected/tests/test_invz_transverse_mf.m` (create)

**Interfaces:**
- Produces: `opts.transverse_mf` ('legacy_x'|'none'|'vector_ab', default 'legacy_x'); new return fields `si.hy` (meV, 0 unless vector mode finds `<Jy> ~= 0`) and `si.transverse_mf` (the resolved mode string). `si.F_mf` gains `0.5*hy*si.Jexp(2)`. All existing fields unchanged in meaning.

- [ ] **Step 1: Write the failing tests** — create `invz_projected/tests/test_invz_transverse_mf.m`:

```matlab
function tests = test_invz_transverse_mf
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_legacy_default_bitforbit(testCase)
% Omitted opt and explicit 'legacy_x' take the identical code path.
ion = invz_ion();
s1 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false));
s2 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'legacy_x'));
verifyTrue(testCase, isequaln(s1, s2));
% digit anchors (controller-verified 2026-07-16): guards against accidental rebaseline
verifyEqual(testCase, s1.E(2)-s1.E(1), 0.369235620278, 'AbsTol', 1e-11);
verifyEqual(testCase, s1.Jexp(2), -0.0689991117463, 'AbsTol', 1e-10);
verifyEqual(testCase, s1.hy, 0);
verifyEqual(testCase, s1.transverse_mf, 'legacy_x');
end

function test_invalid_mode_errors(testCase)
ion = invz_ion();
verifyError(testCase, ...
    @() invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'diagonal')), ...
    'invz:transverseMF');
end

function test_none_equals_bare(testCase)
% 'none' == the current Jxx0 = 0 bare CF+Zeeman calculation.
ion = invz_ion();
sn = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'none'));
sb = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'Jxx0', 0));
verifyEqual(testCase, sn.E,    sb.E,    'AbsTol', 1e-14);
verifyEqual(testCase, sn.Jexp, sb.Jexp, 'AbsTol', 1e-14);
verifyEqual(testCase, sn.hx, 0);  verifyEqual(testCase, sn.hy, 0);
end

function test_vector_c4_axes(testCase)
% C4 is exact for the CF: a/b-axis fields must be equivalent under vector MF.
% Rotation sense pinned: +90 deg about c maps x->y, so <J>([0 B 0]) = Rz(90)*<J>([B 0 0]).
ion = invz_ion();  o = struct('hyp', false, 'transverse_mf', 'vector_ab');
sx = invz_single_ion(ion, 0.31, [4 0 0], o);
sy = invz_single_ion(ion, 0.31, [0 4 0], o);
verifyEqual(testCase, sy.E, sx.E, 'AbsTol', 1e-10);
verifyEqual(testCase, sy.Jexp(2),  sx.Jexp(1), 'AbsTol', 1e-9);   % <Jy>' = <Jx>
verifyEqual(testCase, sy.Jexp(1), -sx.Jexp(2), 'AbsTol', 1e-9);   % <Jx>' = -<Jy>
end

function test_vector_hy_selfconsistent(testCase)
ion = invz_ion();
s = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'vector_ab'));
verifyTrue(testCase, s.mf_converged);
verifyEqual(testCase, s.hy, ion.Jxx0*s.Jexp(2), 'AbsTol', 1e-12);
verifyTrue(testCase, abs(s.hy) > 1e-6);   % B64s => <Jy> ~= 0 even for an x-axis field
end

function test_vector_inplane_Jz_zero(testCase)
% Theta x C2 protects <Jz> = 0 exactly for any in-plane field (para path).
ion = invz_ion();
s = invz_single_ion(ion, 0.31, 4*[cosd(30) sind(30) 0], ...
                    struct('hyp', false, 'transverse_mf', 'vector_ab'));
verifyEqual(testCase, s.Jexp(3), 0, 'AbsTol', 1e-9);
end

function test_vector_Fmf_c4_invariant(testCase)
% F_mf must include the 0.5*hy*<Jy> term: under C4 the (hx,hy) weight swaps
% between channels, so F_mf([4 0 0]) == F_mf([0 4 0]) only if hy is counted.
ion = invz_ion();  o = struct('hyp', false, 'transverse_mf', 'vector_ab');
sx = invz_single_ion(ion, 0.31, [4 0 0], o);
sy = invz_single_ion(ion, 0.31, [0 4 0], o);
verifyEqual(testCase, sy.F_mf, sx.F_mf, 'AbsTol', 1e-10);
% self-consistency identity: 0.5*(hx<Jx>+hy<Jy>) == (hx^2+hy^2)/(2*Jxx0)
lhs = 0.5*(sx.hx*sx.Jexp(1) + sx.hy*sx.Jexp(2));
verifyEqual(testCase, lhs, (sx.hx^2 + sx.hy^2)/(2*ion.Jxx0), 'AbsTol', 1e-12);
end

function test_vector_ordered_mode(testCase)
% hy iterates alongside hz in order mode.
ion = invz_ion();
s = invz_single_ion(ion, 0.31, [1 0.5 0], ...
    struct('hyp', false, 'order', true, 'transverse_mf', 'vector_ab'));
verifyTrue(testCase, s.mf_converged);
verifyEqual(testCase, s.hy, ion.Jxx0*s.Jexp(2), 'AbsTol', 1e-12);
end

function test_vector_hzfixed_mode(testCase)
% hy iterates alongside a held-fixed hz; F_mf stays NaN under hz_fixed.
ion = invz_ion();
s = invz_single_ion(ion, 0.31, [2 1 0], ...
    struct('hyp', false, 'hz_fixed', 0.01, 'transverse_mf', 'vector_ab'));
verifyTrue(testCase, s.mf_converged);
verifyEqual(testCase, s.hy, ion.Jxx0*s.Jexp(2), 'AbsTol', 1e-12);
verifyTrue(testCase, isnan(s.F_mf));
end
```

- [ ] **Step 2: Run to verify failure** — `"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests/test_invz_transverse_mf.m'); assertSuccess(results)"` from repo root. Expected: FAIL (`si.hy`/`si.transverse_mf` missing; no `invz:transverseMF` error).

- [ ] **Step 3: Implement.** In `invz_projected/invz_single_ion.m`:

(a) After the `Jxx0` line (currently line 19) insert:

```matlab
tmf = 'legacy_x'; if isfield(opts,'transverse_mf'), tmf = opts.transverse_mf; end
if ~(ischar(tmf) || isstring(tmf)) || ~any(strcmp(tmf, {'legacy_x','none','vector_ab'}))
    error('invz:transverseMF', ...
        'transverse_mf must be ''legacy_x'', ''none'' or ''vector_ab''.');
end
tmf = char(tmf);
vecmf  = strcmp(tmf, 'vector_ab');
nonemf = strcmp(tmf, 'none');
```

(b) Next to `hx = 0;` add `hy = 0;`.

(c) In the MF loop, change `H = H0 - hx*Jx - hz*Jz;` to `H = H0 - hx*Jx - hy*Jy - hz*Jz;` and replace the `hx_new = Jxx0*jx;` line with:

```matlab
    if nonemf
        hx_new = 0;  hy_new = 0;
    else
        hx_new = Jxx0*jx;
        hy_new = 0;
        if vecmf
            jy = real(diag(V'*Jy*V)).'*p;
            hy_new = Jxx0*jy;
        end
    end
```

(d) Residual and mixing include the y channel:

```matlab
    dmf = max([abs(hx_new - hx), abs(hy_new - hy), abs(hz_new - hz)]);
    if dmf < 1e-12
        hx = hx_new;  hy = hy_new;  hz = hz_new;
        converged = true;
        break;
    end
    hx = hx + mix*(hx_new - hx);
    hy = hy + mix*(hy_new - hy);
    hz = hz + mix*(hz_new - hz);
```

(e) Final recompute: `H = H0 - hx*Jx - hy*Jy - hz*Jz;`. After `si.hx = hx;` add `si.hy = hy;` and after the diagnostics block add `si.transverse_mf = tmf;`.

(f) F_mf: `si.F_mf = Fsite + 0.5*(hx*si.Jexp(1) + hy*si.Jexp(2) + hz*si.Jexp(3));` and update the comment above it to `equals Fsite + (hx^2+hy^2)/(2Jxx0) + hz^2/(2J0z)`.

(g) Header doc: document the three modes, `si.hy`, `si.transverse_mf`, and that `'none'` supersedes `Jxx0` (both channels forced to zero).

- [ ] **Step 4: Run the new test file (PASS), then the full fast suite (81+9 passed / 0 failed / 12 incomplete).**

- [ ] **Step 5: Commit** — `feat(invz): transverse_mf modes (legacy_x|none|vector_ab) in invz_single_ion`

---

### Task 2: Thread `transverse_mf` through the local-state chain

**Files:**
- Modify: `invz_projected/invz_twolevel.m`, `invz_projected/invz_twolevel_ordered.m`, `invz_projected/invz_solve_point.m`, `invz_projected/invz_solve_point_ordered.m`, `invz_projected/invz_solve_auto.m` (doc only), `invz_projected/invz_chi_realaxis.m`, `invz_projected/invz_chi_tensor_ref.m`
- Test: `invz_projected/tests/test_invz_transverse_mf_threading.m` (create)

**Interfaces:**
- Consumes: Task 1's `opts.transverse_mf` contract.
- Produces: every function above accepts `opts.transverse_mf` (default `'legacy_x'`) and forwards it to EVERY `invz_single_ion` call it makes, so the electronuclear state (chi0) and the electronic two-level state (Sigma) always share one MF model. `invz_twolevel`/`invz_twolevel_ordered` additionally record `tl.transverse_mf`. `invz_solve_auto` needs no code change (it forwards `opts` wholesale) — update its header doc only.

- [ ] **Step 1: Write the failing tests** — create `invz_projected/tests/test_invz_transverse_mf_threading.m`:

```matlab
function tests = test_invz_transverse_mf_threading
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'Jcc0', 6.4e-3);
end

function test_solve_point_legacy_default_identical(testCase)
fx = fixture();  ion = invz_ion();
p1 = invz_solve_point(ion, 0.31, [3 0 0], fx.Jnu, struct('J0eff', fx.Jcc0));
p2 = invz_solve_point(ion, 0.31, [3 0 0], fx.Jnu, ...
                      struct('J0eff', fx.Jcc0, 'transverse_mf', 'legacy_x'));
verifyTrue(testCase, isequaln(p1, p2));
end

function test_solve_point_vector_consistent(testCase)
% si and tl must carry the same MF model; hy is live in both.
fx = fixture();  ion = invz_ion();
pt = invz_solve_point(ion, 0.31, [3 1 0], fx.Jnu, ...
                      struct('J0eff', fx.Jcc0, 'transverse_mf', 'vector_ab'));
verifyEqual(testCase, pt.si.transverse_mf, 'vector_ab');
verifyEqual(testCase, pt.tl.transverse_mf, 'vector_ab');
verifyTrue(testCase, abs(pt.si.hy) > 1e-8);
end

function test_solve_auto_inplane_vector(testCase)
% In-plane vector field routes transversely (Bz = 0) and forwards the mode.
fx = fixture();  ion = invz_ion();
[pt, phase] = invz_solve_auto(ion, 0.31, [3 0.5 0], fx.Jnu, ...
                              struct('J0eff', fx.Jcc0, 'transverse_mf', 'vector_ab'));
verifyEqual(testCase, phase, 'para');
verifyEqual(testCase, pt.si.transverse_mf, 'vector_ab');
end

function test_solve_point_c4_axes(testCase)
% cc channel is C4-invariant: a/b-axis solves must agree in vector mode.
fx = fixture();  ion = invz_ion();
o = struct('J0eff', fx.Jcc0, 'transverse_mf', 'vector_ab');
px = invz_solve_point(ion, 0.31, [4 0 0], fx.Jnu, o);
py = invz_solve_point(ion, 0.31, [0 4 0], fx.Jnu, o);
verifyTrue(testCase, px.converged && py.converged);
verifyEqual(testCase, py.Sigma0, px.Sigma0, 'AbsTol', 1e-9);
verifyEqual(testCase, py.crit,   px.crit,   'AbsTol', 1e-9);
end

function test_chi_tensor_ref_vector_mode(testCase)
ion = invz_ion();  w = (0:0.01:0.6).';
o = struct('Jsel', 6.4e-3, 'Jaa0', 3.5e-3, 'eta', 0.02, 'transverse_mf', 'vector_ab');
R = invz_chi_tensor_ref(ion, 0.1, 3*[cosd(20) sind(20) 0], w, o);
verifyTrue(testCase, isfinite(R.eps_spec) && isfinite(R.Epeak_ten));
end

function test_twolevel_ordered_vector(testCase)
ion = invz_ion();
tl = invz_twolevel_ordered(ion, 0.31, [2 1 0], 0.01, ...
                           struct('transverse_mf', 'vector_ab'));
verifyEqual(testCase, tl.transverse_mf, 'vector_ab');
end
```

Note: check `invz_twolevel_ordered`'s actual signature (`(ion, T, Bx, hz, opts)` expected) and `invz_solve_point`'s opts names against the source before finalizing the tests; adjust the calls (NOT the assertions) if the argument order differs. If `invz_solve_point` takes couplings via a different opt than `J0eff`, mirror what `invz_projected/tests/test_invz_solve_auto.m` does for the fixture.

- [ ] **Step 2: Run to verify failure** (missing `tl.transverse_mf`, mode not forwarded).

- [ ] **Step 3: Implement.** Uniform pattern — in each function, next to the existing `Jxx0` opt parse, add:

```matlab
tmf = 'legacy_x';  if isfield(opts,'transverse_mf'), tmf = opts.transverse_mf; end
```

then append `'transverse_mf', tmf` to every `struct(...)` handed to `invz_single_ion` in that file. Specifics:

- `invz_twolevel.m`: forward into the line-8 call; add `tl.transverse_mf = tmf;` with the other tl fields.
- `invz_twolevel_ordered.m`: forward into the line-14 call (keeps `hz_fixed`); add `tl.transverse_mf = tmf;`.
- `invz_solve_point.m`: forward into both the `invz_single_ion` call and the `invz_twolevel` opts struct.
- `invz_solve_point_ordered.m`: add to the `siopts` struct construction; forward into the `invz_twolevel_ordered` opts.
- `invz_chi_realaxis.m`: parse `tmf` from its opts, forward into the paramagnet fallback `invz_single_ion` call (line ~38).
- `invz_chi_tensor_ref.m`: `tmf = getf(opts, 'transverse_mf', 'legacy_x');` forwarded into its order-mode `invz_single_ion` call.
- `invz_solve_auto.m`: header doc gains one line noting `opts.transverse_mf` passes through to both routes.

- [ ] **Step 4: Run the new test file (PASS) + full fast suite (0 failed).**
- [ ] **Step 5: Commit** — `feat(invz): thread transverse_mf through solvers, two-level, chi_realaxis, tensor ref`

---

### Task 3: `phi_ab` driver knob, spectra guard, metadata

**Files:**
- Modify: `invz_projected/invz_spectra_map.m`, `invz_projected/invz_spectra_qpath.m`, `invz_projected/invz_run_spectra.m`
- Test: `invz_projected/tests/test_invz_phi_spectra.m` (create)

**Interfaces:**
- Consumes: Task 2 (`transverse_mf` legal in `solve_opts`, forwarded by solvers).
- Produces: `invz_spectra_map`/`invz_spectra_qpath` error (`invz:transverseMF`) on any nonzero b-component under `'legacy_x'`; both return `S.transverse_mf`. Driver knobs `phi_ab` (deg) and `transverse_mf` with `dhat = [cosd(theta_c)*cosd(phi_ab), cosd(theta_c)*sind(phi_ab), sind(theta_c)]`.

- [ ] **Step 1: Write the failing tests** — create `invz_projected/tests/test_invz_phi_spectra.m`:

```matlab
function tests = test_invz_phi_spectra
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', ...
            'info', struct('Jcc0', 6.4e-3), 'verbose', false);
end

function test_map_rejects_legacy_with_by(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(20) sind(20) 0];
verifyError(testCase, @() invz_spectra_map(ion, 0.31, [3 5.5], w, fx), ...
            'invz:transverseMF');
end

function test_qpath_rejects_legacy_with_by(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
qp = [0 0 0; 0.5 0 0];
verifyError(testCase, ...
    @() invz_spectra_qpath(ion, 0.31, 3*[cosd(20) sind(20) 0], qp, w, struct()), ...
    'invz:transverseMF');
end

function test_map_vector_mode_accepted(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(20) sind(20) 0];
fx.solve_opts = struct('transverse_mf', 'vector_ab');
S = invz_spectra_map(ion, 0.31, [3 5.5], w, fx);
verifyEqual(testCase, S.transverse_mf, 'vector_ab');
verifyTrue(testCase, any(isfinite(S.chiz(:))));
end

function test_map_metadata_default_legacy(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
S = invz_spectra_map(ion, 0.31, [3 5.5], w, fixture());
verifyEqual(testCase, S.transverse_mf, 'legacy_x');
end

function test_map_c4_phi_plus_90(testCase)
% phi and phi+90 give the same cc spectrum in vector mode (C4 exact).
ion = invz_ion();  w = (0.05:0.05:0.6).';
fx = fixture();  fx.solve_opts = struct('transverse_mf', 'vector_ab');
fx.field_dir = [cosd(20) sind(20) 0];
S1 = invz_spectra_map(ion, 0.31, 3.0, w, fx);
fx.field_dir = [cosd(110) sind(110) 0];
S2 = invz_spectra_map(ion, 0.31, 3.0, w, fx);
verifyEqual(testCase, S2.chiz, S1.chiz, 'AbsTol', 1e-8);
end
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement.**

(a) `invz_spectra_map.m`: after `sopts` assembly and the `BvecM` dead-band block, add (reusing the file's `getf`):

```matlab
tmf = getf(sxtra, 'transverse_mf', 'legacy_x');
if strcmp(tmf, 'legacy_x') && any(abs(BvecM(:,2)) > 0)
    error('invz:transverseMF', ['field has a b-axis (y) component but transverse_mf is ' ...
        '''legacy_x'' (x-only mean field; C4-inconsistent, 17 ueV a/b asymmetry at 4 T). ' ...
        'Set opts.solve_opts.transverse_mf = ''vector_ab'' (or ''none'' for bare diagnostics).']);
end
```

and `S.transverse_mf = tmf;` next to `S.field_dir`. Document both in the header opts/returns lists.

(b) `invz_spectra_qpath.m`: same guard on `Bvec(2)` after its field normalization; add `S.transverse_mf = tmf;`. Header doc.

(c) `invz_run_spectra.m`: next to the `theta_c` knob add:

```matlab
phi_ab = 0.0;                        % deg -- IN-PLANE rotation of the swept field, a -> b.
                                     % phi_ab ~ 11 deg is the production experimental angle
                                     % (external stack ion.cfRot(Ho) = -11 deg; the exact
                                     % sign mapping is pinned in test_invz_cfrot_equiv).
                                     % Nonzero phi_ab REQUIRES transverse_mf = 'vector_ab'
                                     % below (the library errors otherwise, by design).
                                     % NOTE: vector_ab shifts even phi_ab = 0 results
                                     % slightly (~0.04 ueV at 4 T; grows at low field) --
                                     % never compare legacy_x and vector_ab runs as if only
                                     % the angle differed. Combined theta_c AND phi_ab is
                                     % NOT validated (tilt bound was measured under legacy_x).
transverse_mf = 'legacy_x';          % 'legacy_x' | 'none' | 'vector_ab'
```

update `dhat = [cosd(theta_c)*cosd(phi_ab), cosd(theta_c)*sind(phi_ab), sind(theta_c)];`, extend `tiltStr` to include `phi_ab` and (when not legacy) the MF mode, and merge `transverse_mf` into every `solve_opts`/opts struct handed to `invz_spectra_map`/`invz_spectra_qpath` in the driver body.

- [ ] **Step 4: Run the new test file (PASS) + full fast suite (0 failed).**
- [ ] **Step 5: Commit** — `feat(invz): phi_ab knob, legacy_x b-component guard, transverse_mf metadata`

---

### Task 4: In-plane diagnostics driver + tensor-reference angle scan

**Files:**
- Create: `invz_projected/invz_c4fit.m`, `invz_projected/invz_run_ip_scan.m`
- Modify: `invz_projected/invz_chi_tensor_ref.m` (add integrated-weight metrics)
- Test: `invz_projected/tests/test_invz_ip_scan.m` (create)

**Interfaces:**
- Consumes: Tasks 1–2.
- Produces: `[A, phi0, resid] = invz_c4fit(phi_deg, y)` — LSQ fit of `y = A(1) + A(2)*cosd(4*phi) + A(3)*sind(4*phi) + A(4)*cosd(8*phi)`, `phi0 = atan2d(A(3), A(2))/4` (principal-axis angle, deg), `resid` = max |fit - y|. `invz_chi_tensor_ref` gains `R.W_sc`, `R.W_ten` (trapz of the wmin-masked chi'' over w) and `R.eps_W = abs(W_sc - W_ten)/W_ten`. `invz_run_ip_scan` prints the IP0/IP3 tables.

- [ ] **Step 1: Write the failing tests** — create `invz_projected/tests/test_invz_ip_scan.m`:

```matlab
function tests = test_invz_ip_scan
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_c4fit_recovers_synthetic(testCase)
phi = (0:5:90).';
y = 1 + 0.1*cosd(4*(phi - 17));
[A, phi0, resid] = invz_c4fit(phi, y);
verifyEqual(testCase, A(1), 1.0, 'AbsTol', 1e-10);
verifyEqual(testCase, hypot(A(2), A(3)), 0.1, 'AbsTol', 1e-10);
verifyEqual(testCase, phi0, 17.0, 'AbsTol', 1e-8);
verifyLessThan(testCase, resid, 1e-10);
end

function test_bare_extrema_anchor(testCase)
% Controller-verified digits (2026-07-16): bare CF+Zeeman at 4 T, T = 0.31 K, hyp = false.
ion = invz_ion();
o = struct('hyp', false, 'transverse_mf', 'none');
s15 = invz_single_ion(ion, 0.31, 4*[cosd(15) sind(15) 0], o);
s60 = invz_single_ion(ion, 0.31, 4*[cosd(60) sind(60) 0], o);
verifyEqual(testCase, s15.E(2)-s15.E(1), 0.345693898281, 'AbsTol', 1e-10);
verifyEqual(testCase, s60.E(2)-s60.E(1), 0.370662538778, 'AbsTol', 1e-10);
end

function test_nonaxial_zero_kills_anisotropy(testCase)
% With all m=4 Stevens terms zeroed the in-plane rotation is exact to numerics.
ion = invz_ion();  ion.B44 = 0;  ion.B64c = 0;  ion.B64s = 0;
o = struct('hyp', false, 'transverse_mf', 'none');
d = zeros(7,1);  phis = 0:15:90;
for k = 1:7
    s = invz_single_ion(ion, 0.31, 4*[cosd(phis(k)) sind(phis(k)) 0], o);
    d(k) = s.E(2) - s.E(1);
end
verifyLessThan(testCase, max(d) - min(d), 1e-10);
end

function test_tensor_ref_weight_metrics(testCase)
ion = invz_ion();  w = (0:0.01:0.6).';
o = struct('Jsel', 6.4e-3, 'Jaa0', 3.5e-3, 'eta', 0.02, 'transverse_mf', 'vector_ab');
R = invz_chi_tensor_ref(ion, 0.1, 3*[cosd(20) sind(20) 0], w, o);
verifyTrue(testCase, isfinite(R.eps_W) && R.eps_W >= 0);
verifyTrue(testCase, R.W_ten > 0 && R.W_sc > 0);
end
```

- [ ] **Step 2: Run to verify failure.**

- [ ] **Step 3: Implement `invz_c4fit.m`:**

```matlab
function [A, phi0, resid] = invz_c4fit(phi_deg, y)
%INVZ_C4FIT Least-squares C4 harmonic fit y = A1 + A2*cos4phi + A3*sin4phi + A4*cos8phi.
% phi0 (deg) = principal-axis angle of the leading harmonic, atan2d(A3,A2)/4.
phi_deg = phi_deg(:);  y = y(:);
M = [ones(size(phi_deg)), cosd(4*phi_deg), sind(4*phi_deg), cosd(8*phi_deg)];
A = M \ y;
phi0 = atan2d(A(3), A(2))/4;
resid = max(abs(M*A - y));
end
```

- [ ] **Step 4: Implement the `invz_chi_tensor_ref.m` additions.** Where the existing wmin-masked `amp_sc`/`amp_ten` are computed (mask `mk = w >= peak_wmin`), add:

```matlab
R.W_sc  = trapz(w(mk), max(imag(chis_sc(mk)), 0));
R.W_ten = trapz(w(mk), max(imag(chis_ten(mk)), 0));
R.eps_W = abs(R.W_sc - R.W_ten)/max(R.W_ten, realmin);
```

using the file's actual variable names for the scalar/tensor cc spectra (read the file; keep naming consistent with `amp_sc`/`amp_ten`). Document the three new fields in the header. `invz_projected/tests/test_invz_tensor_ref.m` must stay green (its reproducibility constants are unaffected by additive fields).

- [ ] **Step 5: Implement `invz_run_ip_scan.m`** (measurement driver, mirrors `invz_run_tensor_ref.m` conventions incl. the addpath prologue and live-coupling block):

```matlab
%INVZ_RUN_IP_SCAN In-plane rotation diagnostics (IP0) + Sigma=0 tensor reference (IP3).
% Section A -- single-ion angular scans at T = 0.31 K, hyp = false, B = [2 4 6] T,
%   phi_ab = [0:5:90 union 11] deg, for the three transverse-MF models
%   ('none' = bare CF+Zeeman, 'legacy_x', 'vector_ab'):
%   Delta(phi) = E2-E1 and its C4 harmonic fit (invz_c4fit), span (max-min)/mean,
%   principal-axis angle phi0. The legacy_x rows exist to display the C4 violation,
%   never for production.
% Section B -- Sigma=0 scalar-vs-tensor cc comparison (invz_chi_tensor_ref,
%   transverse_mf = 'vector_ab') at T = 0.1 K, fields [2 4.95 6] T,
%   phi_ab = [0 5 11 15 30 45 60 75 90] deg, w = (0:0.005:0.6), eta = 0.02:
%   per row dE_peak, eps_amp, eps_W, Epeak_sc/ten; gate per the tilt criterion
%   (eps_amp <= 0.10 AND dE_peak <= max(0.02*Epeak_ten, eta)).
% Copy both printed tables into docs/SESSION-2026-07-16-inplane-rotation.md.
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
ion = invz_ion();

% ---- Section A: single-ion angular scans ----
T_A   = 0.31;  fieldsA = [2 4 6];  phis = unique([0:5:90, 11]);
modes = {'none', 'legacy_x', 'vector_ab'};
fprintf('%10s %6s %10s %10s %10s %8s\n', 'model', 'B(T)', 'span(%)', 'A4/A0(%)', 'A8/A0(%)', 'phi0');
for im = 1:numel(modes)
    for ib = 1:numel(fieldsA)
        B = fieldsA(ib);  d = zeros(numel(phis), 1);
        for k = 1:numel(phis)
            s = invz_single_ion(ion, T_A, B*[cosd(phis(k)) sind(phis(k)) 0], ...
                                struct('hyp', false, 'transverse_mf', modes{im}));
            d(k) = s.E(2) - s.E(1);
        end
        [A, phi0, ~] = invz_c4fit(phis, d);
        fprintf('%10s %6.1f %10.3f %10.3f %10.3f %8.2f\n', modes{im}, B, ...
            100*(max(d)-min(d))/mean(d), 100*hypot(A(2),A(3))/A(1), 100*abs(A(4))/A(1), phi0);
    end
end

% ---- Section B: Sigma=0 scalar-vs-tensor over angle (production couplings) ----
[qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
qc = qc(any(abs(qc) > 1e-12, 2), :);
[~, info] = invz_jq_modes(ion, qc, struct('dpRng', 30, 'cache', true));
Jaa0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end
T_B = 0.1;  fieldsB = [2 4.95 6];  phisB = [0 5 11 15 30 45 60 75 90];
w = (0:0.005:0.6).';  eta = 0.02;
ropts = struct('Jsel', info.Jcc0, 'Jaa0', Jaa0, 'eta', eta, 'transverse_mf', 'vector_ab');
fprintf('\n%8s %8s %12s %12s %12s %10s %10s %5s\n', 'phi', '|B| (T)', 'dE_peak', 'eps_amp', 'eps_W', 'Ep_sc', 'Ep_ten', 'ok');
supported = true(size(phisB));
for ib = 1:numel(fieldsB)
    for ia = 1:numel(phisB)
        ph = phisB(ia);
        R = invz_chi_tensor_ref(ion, T_B, fieldsB(ib)*[cosd(ph) sind(ph) 0], w, ropts);
        ok = R.eps_amp <= 0.10 && ...
             ( (isnan(R.dE_peak) && isnan(R.Epeak_sc) == isnan(R.Epeak_ten)) || ...
               R.dE_peak <= max(0.02*R.Epeak_ten, eta) );
        supported(ia) = supported(ia) && ok;
        fprintf('%8.1f %8.2f %12.4g %12.4g %12.4g %10.4f %10.4f %5d\n', ...
            ph, fieldsB(ib), R.dE_peak, R.eps_amp, R.eps_W, R.Epeak_sc, R.Epeak_ten, ok);
    end
end
fprintf('\nsupported in-plane angles (all fields): %s deg\n', mat2str(phisB(supported)));
```

- [ ] **Step 6: Run the new test file (PASS) + full fast suite (0 failed).** Do NOT run `invz_run_ip_scan` Section B in this task (Task 5 does the production run).
- [ ] **Step 7: Commit** — `feat(invz): invz_run_ip_scan IP0 diagnostics, invz_c4fit, tensor-ref weight metrics`

---

### Task 5: cfRot convention pin, production validation, docs

**Files:**
- Create: `invz_projected/tests/test_invz_cfrot_equiv.m`, `docs/SESSION-2026-07-16-inplane-rotation.md`
- Modify: `invz_projected/README.html` (knob + mode + supported-range documentation)

**Interfaces:**
- Consumes: everything above.
- Produces: the numerically pinned mapping `ion.cfRot = r` (external `cf.m` 'coefficient' convention) ⇔ `phi_ab = f(r)` with sign, asserted by test; the measured IP3 table; the supported-regime statement.

- [ ] **Step 1: Write the cfRot equivalence test** — create `invz_projected/tests/test_invz_cfrot_equiv.m`. Physics: the external stack's `cf.m` ('coefficient') rotates the m=4 coefficient pairs by `Br = [cos(4r) sin(4r); -sin(4r) cos(4r)]`. Rotating the CF coefficients by `r` must equal rotating the in-plane FIELD by `-r` or `+r` (one of the two — the test discovers the sign, then asserts it):

```matlab
function tests = test_invz_cfrot_equiv
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function E = spec_rotcf(ion, r_deg, B)
% Spectrum with the m=4 CF coefficient pairs rotated by r (cf.m 'coefficient' method).
% invz_ion has B44s = 0 implicitly; rotation generates it. Build Hcf manually.
o = stevens_ops(ion.J);
c4 = cosd(4*r_deg);  s4 = sind(4*r_deg);
B44c =  c4*ion.B44;                 B44s = -s4*ion.B44;
B64c =  c4*ion.B64c + s4*ion.B64s;  B64s = -s4*ion.B64c + c4*ion.B64s;
assert(isfield(o, 'O44s') && isfield(o, 'O64s'), 'stevens_ops must expose O44s/O64s');
Hcf = ion.B20*o.O20 + ion.B40*o.O40 + B44c*o.O44 + B44s*o.O44s ...
    + ion.B60*o.O60 + B64c*o.O64c + B64s*o.O64s;
C = invz_const();
H = Hcf - ion.gL*C.muB*(B(1)*o.Jx + B(2)*o.Jy + B(3)*o.Jz);
E = sort(real(eig((H + H')/2)));  E = E - E(1);
end

function E = spec_rotfield(ion, phi_deg, Bmag)
o = stevens_ops(ion.J);
Hcf = ion.B20*o.O20 + ion.B40*o.O40 + ion.B44*o.O44 ...
    + ion.B60*o.O60 + ion.B64c*o.O64c + ion.B64s*o.O64s;
C = invz_const();
B = Bmag*[cosd(phi_deg) sind(phi_deg) 0];
H = Hcf - ion.gL*C.muB*(B(1)*o.Jx + B(2)*o.Jy + B(3)*o.Jz);
E = sort(real(eig((H + H')/2)));  E = E - E(1);
end

function test_pin_cfrot_field_mapping(testCase)
% Discover then assert the sign: coefficient rotation by r == field rotation by s*r.
ion = invz_ion();  r = -11;  Bmag = 4;
Erot = spec_rotcf(ion, r, [Bmag 0 0]);
Ep = spec_rotfield(ion, +r, Bmag);   % candidate s = +1
Em = spec_rotfield(ion, -r, Bmag);   % candidate s = -1
dp = max(abs(Ep - Erot));  dm = max(abs(Em - Erot));
% exactly one candidate matches (B64s breaks the mirrors, so +11 and -11 differ)
verifyLessThan(testCase, min(dp, dm), 1e-10);
verifyGreaterThan(testCase, max(dp, dm), 1e-6);
% ASSERT THE DISCOVERED SIGN HERE: after the first run, replace the min() check with
% a hard assertion on the matching candidate and record the mapping in the comment:
%   external ion.cfRot = -11 deg  <=>  invz phi_ab = <SIGN>11 deg
end

function test_mapping_holds_off_axis(testCase)
% Same mapping at a second rotation angle and field, guards against 4r-aliasing.
ion = invz_ion();  r = 7;  Bmag = 6;
Erot = spec_rotcf(ion, r, [Bmag 0 0]);
dp = max(abs(spec_rotfield(ion, +r, Bmag) - Erot));
dm = max(abs(spec_rotfield(ion, -r, Bmag) - Erot));
verifyLessThan(testCase, min(dp, dm), 1e-10);
end
```

Implementation notes for the test: (i) check `stevens_ops` for the exact field names of the sine operators (`O44s`/`O64s`); if `O44s` is absent, add it to `stevens_ops.m` as `O44s = (Jp^4 - Jm^4)/(2i)` matching the cosine convention already there (additive change, existing outputs untouched); (ii) after the first passing run, HARDEN `test_pin_cfrot_field_mapping` to assert the specific matching sign and update the comment — the discovered mapping is a deliverable, not just a pass.

- [ ] **Step 2: Run the equivalence test; harden the sign assertion; re-run (PASS).**

- [ ] **Step 3: Production validation run.** Run `invz_run_ip_scan` end to end:

```bash
"/Applications/MATLAB_R2025a.app/bin/matlab" -batch "run('invz_projected/invz_run_ip_scan.m')"
```

Capture both printed tables verbatim. Expected shape (from the reviewed Codex diagnostics): Section A spans ~5-8% growing with field, extrema near 15/60 deg, legacy_x rows visibly C4-broken; Section B dE_peak of order a few ueV at 6 T (Codex measured max tensor-induced peak shift ~2.45 ueV at 6 T, inside the gate). If any Section B row FAILS the gate, do not massage thresholds — report BLOCKED with the numbers.

- [ ] **Step 4: Write `docs/SESSION-2026-07-16-inplane-rotation.md`:** feature summary (modes, knob, guard), both measured tables, the supported-regime statement (peak observables, angle x field grid, thresholds restated), the pinned cfRot⇔phi_ab mapping with sign and the sentence "phi_ab = <sign>11 deg reproduces the experimentally calibrated external-stack geometry (ion.cfRot(Ho) = -11 deg)", the vector-vs-legacy phi=0 shift caveat (never cross-compare MF models), and the explicit scope exclusion: combined theta_c + phi_ab is unvalidated; theta_c <= 5 deg was measured under legacy_x only.

- [ ] **Step 5: Update `invz_projected/README.html`:** add `phi_ab`/`transverse_mf` to the driver-knob documentation (Section 4), one paragraph in Section 1 (or 1.7) on the vector transverse MF and why legacy_x is C4-inconsistent for rotated fields (cite the 17 ueV a/b asymmetry), `invz_run_ip_scan` + `invz_c4fit` rows in the function reference (Section 6), and the supported-range + scope-exclusion statement in Section 8. Keep the existing HTML style (MathJax inline, same section markup).

- [ ] **Step 6: Full fast suite green; run INVZ_SLOW suite:**

```bash
INVZ_SLOW=1 "/Applications/MATLAB_R2025a.app/bin/matlab" -batch "results = runtests('invz_projected/tests'); assertSuccess(results)"
```

Expected: 0 failed (slow count grows only if slow tests were added — none planned).

- [ ] **Step 7: Commit** — `feat(invz): cfRot<->phi_ab mapping pinned, in-plane validation tables, docs`

---

## Self-Review Notes (controller)

- Spec coverage vs Codex plan + amendments: §5.1→Task 1, §5.2→Task 2, §5.3/5.4→Task 3, §6/IP0→Task 4 (amendment 1: extends `invz_chi_tensor_ref`), §7 acceptance tests distributed across task test files, IP3→Task 5. Descoped by design: `Jperp0` matrix (amendment 4), tensorized 1/z (explicitly deferred to `invz_tensor/`), combined-angle validation (amendment 3, documented not implemented).
- Type consistency: `transverse_mf` is a char row vector everywhere ('legacy_x' canonical); `si.hy` meV scalar; `tl.transverse_mf` char; `R.W_sc/W_ten/eps_W` doubles; `invz_c4fit` returns `A` 4x1, `phi0` scalar deg, `resid` scalar.
- Known unknowns delegated with instructions: exact `invz_twolevel_ordered`/`invz_solve_point` opts spelling (Task 2 note), `invz_chi_tensor_ref` internal spectrum variable names (Task 4 step 4), `stevens_ops` sine-operator field names (Task 5 note).
