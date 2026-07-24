# Ewald integration surface map — REVIEW-CORRECTED implementation contract

> **STATUS ADDENDUM (2026-07-24 — supersedes the pre-freeze status language below).** The
> preregistration `docs/invzp_ewald_prereg.md` was **FROZEN 2026-07-24** (with §12 Errata E1),
> incorporating this map's Blocker-5/5A/6 contracts (the wrapper API, the spectra-forwarder
> contract, and the `Jpath_base_cc`/`Jgamma_cc` q-path contract). The Step-4 primitive
> `invz_dipole_ewald.m` was subsequently **built, reviewed, and committed**
> (`086d102..fcb3031`, additive; production default still `bruteforce`), and the Step-4
> manifest/status closeout landed in `a58bdb4` (suite 270/0/23). Consequently, in the header
> below: **"No production `.m` file has been modified"** remains true only for the pre-existing
> production files this map discusses (`invz_jq_modes.m`, `invz_bz_couplings.m`,
> `invz_jq_path.m`, the spectra drivers) — `invz_dipole_ewald.m` itself now exists as a new,
> additive primitive file, not a modification of an existing one; **"`invz_dipole_ewald.m` and
> `docs/invzp_ewald_prereg.md` do not exist yet"** is a **historical pre-freeze snapshot and is
> no longer operative** — both now exist; and **"Numerical Γ and cutoff values cited here are
> design calibration, not formal post-preregistration gate results"** is likewise superseded —
> those values are now the FROZEN prereg §2/§3 values (`docs/invzp_ewald_prereg.md`), not mere
> calibration, though the caller-level Gate-C4/C6/C7 integration checks this map anticipates
> remain prospective until Step 5 lands them. Step 5 (opt-in `invz_jq_modes`/BZ/q-path/spectra
> integration per this map's Blocker 5/5A/6 contracts, production default unchanged) is the
> current authorized phase per `docs/invzp_ewald_design.md` §8. Current operative status and
> roadmap: `docs/invzp_ordered_1z_state.md`. Per the freeze rule, the original review text is
> retained unedited below this note rather than rewritten.

**Status: REVIEW-CORRECTED, READ-ONLY design closure.** This fixes the Gate-E wrapper contract,
`invz_jq_path` Γ/metadata contract, and the previously omitted user-facing spectra forwarders.
It proposes APIs/contracts only. **No production `.m` file has been modified.**
`invz_dipole_ewald.m` and `docs/invzp_ewald_prereg.md` do not exist yet. Numerical Γ and cutoff
values cited here are design calibration, not formal post-preregistration gate results.

---

## Blocker 5 — physics-anchor wrapper API (Gate E / design §5.3)

### 5.1 What each anchor computes today (exact call chains)

**`invz_projected/tests/test_invz_sigma_crit.m`** — two local test functions:

- `test_fcc_watson` (`test_invz_sigma_crit.m:12-18`): a **pure synthetic fcc lattice** check
  (`fcc_jq`, `:20-31`) — never touches `invz_ion`, `invz_jq_modes`, or any dipole primitive.
  Not a candidate for backend injection; out of scope for this blocker.
- `test_lihof4_sigma_crit` (`:33-50`) — **the Gate-E anchor**:
  `ion = invz_ion()` (`:39`) → loop `ns=[12 24]` (`:40-46`) → per grid size,
  `qVec_generator(ion.a,'mode','grid','grid',[ns(k) ns(k) ns(k)],'range',[-0.5 0.5])` via
  `evalc` (`:42`) → drop `abs(q)<=1e-12` rows (`:43`) → `[Jnu,info] =
  invz_jq_modes(ion,qvec,struct('dpRng',30,'cache',true))` (`:44`) → `S(k) =
  invz_sigma_crit(info.Jcc0, Jnu(:))` (`:45`). Richardson-extrapolates `Sc=2*S(2)-S(1)`
  (`:47`). `info` retained from the **last** loop iteration (24³) for the `Jcc0` assertion.
  Asserts (`:48-49`): `Sc` `AbsTol 0.006` vs `0.3004`; `info.Jcc0` `RelTol 0.03` vs `6.421e-3`
  meV.

**`invz_projected/tests/test_invz_critical.m`** — six local test functions; the design's Gate
E list (design §6 Gate E) names only **two**:

- `test_zero_field_tc` (`:24-46`, **Gate-E anchor**): its **own inline** `ns=[12 24]` Richardson
  loop (`:32-39`, structurally identical to sigma_crit's, **not** via the file's own
  `lihof4_couplings()` helper) → `Sc=2*S(2)-S(1)` (`:40`), `J0=inf_k.Jcc0` (`:41`, last-iteration
  info) → `Tc = invz_critical_T0field(ion,Sc,J0)` (`:42`), `TcMF =
  invz_critical_T0field(ion,0,J0)` (`:43`). Asserts: `Tc` `AbsTol 0.08` vs `1.74` K (`:44`);
  `TcMF > Tc` qualitatively (`:45`).
- `test_critical_field_at_310mK` (`:48-57`, **Gate-E anchor**): uses the file's
  `lihof4_couplings()` helper (`:12-22`: single **16³** grid, `dpRng=30`) → `bx =
  invz_critical(ion,0.31,Jf,struct('J0eff',J0))` (`:53`), asserted `4.0 < bx < 4.6` T (`:54`) →
  `pt = invz_solve_point(ion,0.31,bx,Jf,struct('hyp',true,'J0eff',J0))` (`:55`), asserted
  `pt.Sigma0` `AbsTol 0.02` vs `0.0932` (`:56`).
- `test_tc_small_field`, `test_tc_at_fixed_field_crossing`, `test_tc_boundary_is_smooth`
  (`:59-99`) — **not** in the design's Gate-E list; these are solver self-consistency/regression
  tests (smoothness, fixed-B/fixed-T crossing agreement), also built on `lihof4_couplings()`'s
  16³ grid. I treat these as **out of scope** for the Ewald candidate wrapper (the design names
  only zero-field `Tc` and 310 mK `Bc`); the same wrapper machinery would trivially cover them
  too if ever wanted, but nothing requires it.

**`invz_projected/tests/test_invz_ordered_phase.m`** — one local test function:

- `test_ordered_solve_and_soft_mode` (`:12-43`, **the Gate-E anchor**): `ion = invz_ion()`
  (`:18`) → single **16³** grid (`:19-21`) → `[Jnu,info] =
  invz_jq_modes(ion,qvec,struct('dpRng',30,'cache',true))` (`:22`) → `o =
  struct('hyp',true,'J0eff',info.Jcc0)` (`:24`) → `invz_solve_point_ordered(ion,T,B,Jnu(:),o)`
  at `B=2.0,4.0,6.0` (`:26,27,41`) → soft mode via `soft_peak` (`:45-49`) →
  `invz_chi_realaxis(ion,T,B,pt,w,struct('Jsel',Jsel))` with `Jsel=info.Jcc0` (`:46`). Asserts
  (`:28-39,42`) are **qualitative/ordinal**: `is_ordered && converged` at B=2,4; `crit>0` at
  both, `crit` decreasing toward `Bc`; `m0` decreasing toward `Bc`; soft-mode peak energy
  `E2>E4`; `E4>0.05` (gapped, never softens to zero); `E2<0.6`; paramagnet (`~is_ordered`) at
  B=6.

This test calls `invz_solve_point_ordered` without `opts.ordered_mode`; the solver default is
`'bare'`. It is therefore a useful legacy physical anchor but **not** a test of the Jensen/HMF
ordered branch that currently masks the FM side of `invz_spectra_map`. Gate E must retain this
legacy test and add a separate target-path Jensen acceptance wrapper whose exact fields and
continuity criteria are frozen in Step 3.

**Downstream consumers are opaque numeric sinks.** A repo-wide grep of every `MF_dipole(` and
`invz_jq_modes(` call site (see Bash transcript) shows `invz_sigma_crit.m`,
`invz_critical.m`, `invz_critical_T0field.m`, `invz_solve_point.m`,
`invz_solve_point_ordered.m`, and `invz_chi_realaxis.m` **never** call `MF_dipole` or
`invz_jq_modes` themselves — they only consume `Jnu(:)`/`info.Jcc0`/`opts.J0eff`/`opts.Jsel` as
plain numbers. None of them need to change for either backend; coupling injection happens
entirely upstream of these calls.

### 5.2 Key finding: none of the three anchors call `invz_bz_couplings.m`

All three inline their **own** `qVec_generator` + `invz_jq_modes` call rather than using the
shared driver entry point `invz_projected/invz_bz_couplings.m:1-19`. This matters because design
§2.3 (`docs/invzp_ewald_design.md:130`) says bz_couplings must forward the new options "so the
parameterized acceptance wrappers can actually inject the candidate" — but today there is no
route from bz_couplings to these frozen files at all.

I confirmed the two constructions are **value-for-value equivalent** at matching grid size:
`invz_bz_couplings.m:14` calls `qVec_generator(ion.a,'mode','grid','grid',grid,'range',[-0.5
0.5],'verbose',false)`; each anchor calls the same signature via `evalc(...)` purely to suppress
stdout (`'verbose'` isn't even passed, but `qVec_generator`'s default grid construction is
identical — the `evalc` wrapper changes nothing numeric). Both apply the identical Γ-row filter
`qc(any(abs(qc)>1e-12,2),:)` (`invz_bz_couplings.m:15`, and e.g. `test_invz_sigma_crit.m:43`).
So `invz_bz_couplings(ion, struct('grid',[N N N],'dpRng',30,'cache',true))` reproduces each
anchor's inline `qvec` bit-for-bit. **Consequence for the design:** the parameterized wrapper
can safely be built as new code that calls (an extended) `invz_bz_couplings.m` instead of
inlining `qVec_generator`+`invz_jq_modes` a second time — without changing what the legacy
files compute, because the two paths agree exactly at the shared defaults.

### 5.3 Existing prior art to reuse, not reinvent

The "half-open grid convention + Γ policy P-complete/P-drop" the task asks about **already has a
working implementation** — just not wired to bz_couplings or the anchors:

- `invz_projected/invz_phase1_offsets.m:1-25` — the eight `{0,1/2}³` offsets as a struct array
  (`.tag`, `.flags`).
- `invz_projected/invz_phase1_qgrid.m:1-148` —
  `g = invz_phase1_qgrid(ion, N, offsetFlags, convention, gammaPolicy)`, with **exactly** the
  literal strings `convention ∈ {'halfopen','legacy_inclusive'}` (`:91-93`) and `gammaPolicy ∈
  {'P_complete','P_drop'}` (`:94-96`, `:121-136`) the design's prose calls "P-complete"/"P-drop"
  (note: code uses an underscore, not a hyphen — I use the code's literal values throughout).
  Returns `g.qvec`, `g.w` (**uniform** weights either way: `1/nominal` for `P_complete`,
  `1/(nominal-n_gamma)` for `P_drop`, `:123,133`), `g.n_gamma`, `g.nominal`, full provenance.
- `invz_projected/invz_phase1_couplings.m:1-39` — `c = invz_phase1_couplings(ion, g, dpRng)`
  already bridges a `g` from `invz_phase1_qgrid` into `invz_jq_modes` (`:31`), exposing
  `c.Jnu_unflat`, `c.Jnu_flat`, `c.info`, `c.J0eff`, `c.Jcc0`, `c.maxJnu` — the exact
  "grid-object → couplings" bridge shape a wrapper needs, just missing dipole-backend
  threading (its `invz_jq_modes` call at `:31` hardcodes `struct('dpRng',dpRng,'cache',true)`,
  no `dipole`/`ewald` fields) and used today only by the Phase-1 structural driver/tests, not by
  Gate-E.

**Concrete finding on today's implicit convention:** `qVec_generator.m:23-25` documents
`'endpoint'` default `true` (inclusive linspace = `invz_phase1_qgrid`'s `'legacy_inclusive'`
label). Neither `invz_bz_couplings.m:14` nor any anchor's inline call ever passes `'endpoint'`,
so **today's production/anchor grid convention is `legacy_inclusive`, not `halfopen`** — the
halfopen construction Phase 1 recommends (design §0) is not exercised by any of Gate E's inputs
today.

### 5.4 Proposed API / contract

**(A) Extend `invz_projected/invz_bz_couplings.m`'s `opts`** (signature unchanged:
`[Jnu,info,Jaa0] = invz_bz_couplings(ion,opts)`; new fields only, all optional):

```matlab
% New opts fields (all optional; getf() default-accessor already used at :11-13, invz_common/getf.m:1-5)
opts.dipole         % absent => implicit legacy brute force; explicit 'bruteforce' | 'ewald'
opts.ewald          % absent by default; complete struct required only with explicit 'ewald'
opts.gridConvention % 'legacy_inclusive' (default) | 'halfopen'   -> invz_phase1_qgrid arg 4
opts.gridOffset     % [1x3] logical, default [0 0 0]              -> invz_phase1_qgrid arg 3 (offsetFlags)
opts.gammaPolicy    % 'P_drop' (default) | 'P_complete'           -> invz_phase1_qgrid arg 5
```

Behavior (proposed body shape, replacing `invz_bz_couplings.m:14-16`):

```matlab
useNewGrid = isfield(opts,'gridConvention') || isfield(opts,'gridOffset') || isfield(opts,'gammaPolicy');
if ~useNewGrid
    % UNCHANGED bit-identical path (today's :14-15) -- design Sec.7 regression mandate
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
else
    if any(grid ~= grid(1))
        error('invz:bzCouplingsAnisotropicHalfopen', ...   % see Risk 1, Sec.5.7
            'gridConvention/gridOffset/gammaPolicy require a cubic grid (invz_phase1_qgrid N is scalar); got %s.', mat2str(grid));
    end
    convention   = getf(opts, 'gridConvention', 'legacy_inclusive');
    gridOffset   = getf(opts, 'gridOffset',     [0 0 0]);
    gammaPolicy  = getf(opts, 'gammaPolicy',    'P_drop');
    g  = invz_phase1_qgrid(ion, grid(1), gridOffset, convention, gammaPolicy);
    qc = g.qvec;
end
jqOpts = struct('dpRng', dpRng, 'cache', cache);
% Preserve the exact old jq_modes call shape when neither field is present.
if isfield(opts, 'dipole'), jqOpts.dipole = opts.dipole; end
if isfield(opts, 'ewald'),  jqOpts.ewald  = opts.ewald;  end
[Jnu, info] = invz_jq_modes(ion, qc, jqOpts);
```

The `useNewGrid` guard is a **presence test** (`isfield`), not a default-value test, so an
explicit `opts.gammaPolicy='P_drop'` still routes through the new (equivalent but not
bit-identical, because `invz_phase1_qgrid` additionally wraps via `mod(qraw+0.5,1)-0.5` at
`:112`, which folds a `legacy_inclusive` grid's `+0.5` endpoint row onto `-0.5` — see Risk 2)
branch, while an absent field never does. `opts.dipole='ewald'` alone (no grid* fields) does
**not** by itself switch grid-construction branches under this proposal — it only changes the
dipole backend passed to `invz_jq_modes`; callers that want the recommended
`halfopen`+Ewald combination must request `gridConvention='halfopen'` explicitly. (An
alternative would auto-imply `halfopen` whenever `dipole='ewald'`; I did not choose that because
design §1.2/§6 treats primitive validity and BZ-quadrature convention as independently-gated
axes, and Gate E's own wrapper — part (C) below — pins the combination explicitly rather than
relying on an implicit default.)

Backend options are validated as a pair: any `opts.ewald` field without
`opts.dipole='ewald'`, `opts.dipole='bruteforce'` with an Ewald field, an incomplete Ewald-control
struct, and unknown backend/boundary labels all error. They are never ignored or silently
normalized into another request.

**(B) New, additive shared helper** (proposed new file — not created):
`invz_projected/tests/invz_anchor_couplings.m` (co-located with its only intended callers, the
new Gate-E wrapper files in part C; `invz_phase1_couplings.m` is the precedent for this kind of
thin, non-test "wiring" helper living next to production code, so `invz_projected/` is an
equally defensible location — a naming/placement detail, not a hard constraint):

```matlab
function [Jnu, info] = invz_anchor_couplings(ion, N, opts)
%INVZ_ANCHOR_COUPLINGS Gate-E parameterized coupling builder (candidate-backend injection point).
% Mirrors each frozen anchor exactly when no optional candidate fields are supplied.
if nargin < 3, opts = struct(); end
bzOpts = struct('grid', [N N N], ...
                'dpRng', getf(opts,'dpRng',30), ...
                'cache', getf(opts,'cache',true));
% Presence is semantic: do not add these fields at defaults, because doing so
% activates invz_bz_couplings' new wrapped-grid branch.
for name = {'dipole','ewald','gridConvention','gridOffset','gammaPolicy'}
    f = name{1};
    if isfield(opts,f), bzOpts.(f) = opts.(f); end
end
[Jnu, info, ~] = invz_bz_couplings(ion, bzOpts);
end
```

At its defaults this now genuinely reproduces each anchor's inline construction: none of the three
grid-policy fields is present, so `invz_bz_couplings` takes its unchanged legacy branch. Candidate
tests explicitly supply `dipole`/`ewald`/`gridConvention`/`gridOffset`/`gammaPolicy`. An additive
parity test must compare its default `Jnu` and all pre-existing `info` fields against each legacy
inline construction; a grid-builder unit test separately compares the generated q arrays. This is
a generalization of
a pattern the codebase already has locally: `test_invz_critical.m:12-22`'s own
`lihof4_couplings()` helper is a single-file, non-parameterized instance of exactly this
"thin coupling-builder used by several test functions in the same file" shape.

**(C) New, additive sibling test files** (proposed — not created; frozen legacy files
**untouched**, satisfying "candidate tests are additive" per design §7):

- `test_invz_sigma_crit_ewald.m` — mirrors `test_lihof4_sigma_crit`'s Richardson(12,24) body,
  replacing `test_invz_sigma_crit.m:42-44` with two calls to `invz_anchor_couplings`:

  ```matlab
  function test_lihof4_sigma_crit_ewald_pcomplete(testCase)
      assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
      assumeTrue(testCase, exist('invz_dipole_ewald','file') == 2, 'Ewald primitive not yet implemented');
      ion = invz_ion();
      S = zeros(1,2); ns = [12 24];
      eopts = struct('dipole','ewald','gridConvention','halfopen','gridOffset',[0 0 0], ...
                      'gammaPolicy','P_complete', ...
                      'ewald', frozen_ewald_opts());   % values generated from the frozen Step-3 prereg
      for k = 1:2
          [Jnu, info] = invz_anchor_couplings(ion, ns(k), eopts);
          S(k) = invz_sigma_crit(info.Jcc0, Jnu(:));
      end
      Sc = 2*S(2) - S(1);
      verifyEqual(testCase, Sc, 0.3004, 'AbsTol', 0.006);          % SAME anchor tolerance, unchanged
      verifyEqual(testCase, info.Jcc0, 6.421e-3, 'RelTol', 0.03);  % SAME anchor tolerance, unchanged
  end
  % test_lihof4_sigma_crit_ewald_pdrop: identical body, 'gammaPolicy','P_drop'
  ```

- `test_invz_critical_ewald.m` — mirrors **only** `test_zero_field_tc` and
  `test_critical_field_at_310mK` (the two Gate-E-named functions), each split into `_pcomplete`/
  `_pdrop` variants (4 new functions total), same assertions/tolerances as
  `test_invz_critical.m:44-45` and `:54,56`.
- `test_invz_ordered_phase_ewald.m` — mirrors `test_ordered_solve_and_soft_mode`, `_pcomplete`/
  `_pdrop` variants, same qualitative assertions as `test_invz_ordered_phase.m:28-39,42`.
- `test_invz_spectra_map_ewald_jensen.m` — new target-path test, `_pcomplete`/`_pdrop` variants,
  calling `invz_spectra_map(..., ordered_1z='jensen')`. Step 3 freezes its absolute field grid and
  acceptance rules. It must check accepted FM columns (`phase_1z=1`, finite `Sigma0`, `m_1z`,
  `D_ord`), stable PM columns, no unexplained masks in the declared window, and the frozen
  ordered-to-PM continuity/softening trend.

Every new test function opens with the legacy `assumeTrue(INVZ_SLOW)` gate **plus** an
additional `assumeTrue(exist('invz_dipole_ewald','file')==2, ...)` gate, so these files can be
added now (design §8 step 5) without failing before the primitive exists (steps 3-4) — inert
scaffolding, not a silent skip of a gate that should hard-fail once the primitive lands.

**What each wrapper overrides** (explicit, so the diff-from-legacy is auditable at a glance):
`dipole` (`'bruteforce'→'ewald'`), `ewald` (absent→frozen candidate `eopts`), `gridConvention`
(`'legacy_inclusive'→'halfopen'`), `gridOffset` (`[0 0 0]` production), and `gammaPolicy` (fixed
per `_pcomplete`/`_pdrop` variant).
**Unchanged:** grid size(s) `N` (same 12/24 or 16 as the legacy sibling), `dpRng` (still
forwarded, though physically inert once Ewald converges), all physics inputs (`T`, `B`, field
windows), and — critically — **all assertion tolerances**, reused verbatim from the legacy file
(§5.6).

### 5.5 How P-complete vs P-drop are both exercised

Each candidate anchor gets exactly **two** variants, not four: `gridConvention` is **pinned** to
`'halfopen'` (the corrected construction the Ewald escalation targets, design §0) ×
`dipole='ewald'`, crossed only with `gammaPolicy ∈ {'P_complete','P_drop'}`. Gate F requires both
to be run and applies its frozen eligibility/cross-policy decision tree; wrappers do not select a
winner. `legacy_inclusive × gammaPolicy` is
Gate D's/Phase-1's own territory, already covered by `test_invz_phase1_quadrature.m`; not
repeated here. This is safe with an **unweighted** `mean()` inside `invz_sigma_crit` (untouched)
because both `gammaPolicy` branches of `invz_phase1_qgrid.m` emit **uniform** weights over their
own row count (`:123,133`) — a plain mean over `Jnu(:)` is the correct quadrature either way. If
weighting is ever generalized to non-uniform grids, `invz_sigma_crit` would need a weighted-mean
overload; flagged as a future risk, not needed here.

`invz_sigma_crit` has an additional policy-sensitive operation: it removes entries satisfying
`J0-J<=1e-12` before taking its mean. On the unshifted half-open P-complete grid, the preregistered
contract is exactly one excluded uniform Γ eigenvalue (and the expected warning); P-drop removes
the Γ row upstream and must exclude zero entries internally. Candidate wrappers must assert these
counts rather than merely tolerate the warning. P-complete and P-drop remain distinct because
P-drop removes all four Γ branches, not only the singular uniform branch.

### 5.6 Each anchor's own documented tolerance (reused verbatim — no blanket bar)

| Anchor | Quantity | Tolerance | File:line |
|---|---|---|---|
| `test_lihof4_sigma_crit` | `Sc` | `AbsTol 0.006` vs `0.3004` | `test_invz_sigma_crit.m:48` |
| | `info.Jcc0` | `RelTol 0.03` vs `6.421e-3` meV | `test_invz_sigma_crit.m:49` |
| `test_zero_field_tc` | `Tc` | `AbsTol 0.08` vs `1.74` K | `test_invz_critical.m:44` |
| | `TcMF` | `TcMF > Tc` (qualitative) | `test_invz_critical.m:45` |
| `test_critical_field_at_310mK` | `bx` | `4.0 < bx < 4.6` T (open interval) | `test_invz_critical.m:54` |
| | `pt.Sigma0` | `AbsTol 0.02` vs `0.0932` | `test_invz_critical.m:56` |
| `test_ordered_solve_and_soft_mode` | `is_ordered`, `converged`, `crit>0`, `crit`/`m0` ordering, `E2>E4`, `0.05<E4`, `E2<0.6`, high-B paramagnet | qualitative/ordinal only | `test_invz_ordered_phase.m:28-39,42` |

### 5.7 Risks

1. **Cubic-grid constraint.** `invz_phase1_qgrid.m:84-86` requires `N` to be a **scalar**
   (cubic grid only); `invz_bz_couplings`'s `opts.grid` is `[1x3]` and some callers use
   anisotropic grids. All three Gate-E anchors use cubic grids (12/16/24), so this doesn't block
   Gate E itself, but the proposed bz_couplings branch (5.4-A) must explicitly **error** (not
   silently truncate) if a non-cubic grid is combined with a non-default `gridConvention`/
   `gammaPolicy`/`gridOffset`, until `invz_phase1_qgrid` is extended to accept `[1x3] N`.
   (Sketched as the `any(grid~=grid(1))` guard above.)
2. **Not bit-identical once the new branch is taken.** `invz_phase1_qgrid.m:112` wraps every
   point via `mod(qraw+0.5,1)-0.5`; for `legacy_inclusive`'s inclusive `+0.5` endpoint row this
   maps it onto `-0.5`, **duplicating** the `-0.5` row instead of keeping a numerically distinct
   `+0.5` row — physically equivalent (periodic image) but not bit-identical to today's plain
   `qVec_generator` output. This is fine because the new branch only fires on an explicit,
   non-default opt-in (never for existing callers), but it means `opts.gridConvention=
   'legacy_inclusive'` requested *through the new branch* is not itself a useful "sanity" path —
   it's really the `halfopen` branch that matters; flagged so no one mistakes the new branch's
   `legacy_inclusive` output for a bit-identical replica of the old default path.
3. **`info.Jcc0` is Γ-only and grid-set-independent.** Both anchors that assert on `info.Jcc0`
   get it from a Γ-point-only priming computation inside `invz_jq_modes` (`invz_jq_modes.m:78,
   93-99`) that does not depend on which non-Γ rows are in `qvec` — so `gammaPolicy` cannot
   silently change the meaning of the `Jcc0` assertion; only the `Sc`/soft-mode quantities that
   integrate over the full `Jnu(:)` multiset are `gammaPolicy`-sensitive. No action needed, noted
   for the reviewer.
4. **ODD+Ewald guard is unreachable through this surface.** Design (`invzp_ewald_design.md:44-
   51`) mandates `opts.dipole='ewald'` + active `opts.odd` must error. `invz_bz_couplings.m`
   never forwards `opts.odd` today (only `dpRng`/`cache`, `:16`), so this combination cannot be
   reached through bz_couplings or the anchors — explicitly out of scope here, not a gap in this
   proposal.
5. **Two "Γ removals" must stay separate (design §4.1).** `gammaPolicy` (`P_complete`/`P_drop`)
   is the **BZ quadrature row-inclusion** policy; it is unrelated to `opts.ewald.boundary=
   'conducting_k0_omitted'` (the lattice sum's exact `k=0` omission plus no primitive surface
   term). Never conflate the two under one flag.
6. **Candidate controls are calibrated but not frozen.** The design recommends
   `alpha0=sqrt(pi)/Vc^(1/3)`, the `C_r={4.5,5.0,5.5}` and `C_g={9,10,11}` ladders, and the
   `conducting_k0_omitted` boundary. Step 3 must freeze exact comparison axes, norms, samples, and
   resource caps. The wrapper API deliberately treats `opts.ewald` as an opaque provenance-bearing
   struct until then.

---

## Blocker 5A — user-facing spectra option/provenance flow

The original mapping stopped at `invz_bz_couplings` and `invz_jq_path`, but the target driver strips
unknown options before either call:

- `invz_spectra_map.m` rebuilds a struct containing only `grid` and `dpRng`;
- `invz_spectra_qpath.m` independently builds one BZ struct and one q-path struct, neither carrying
  a backend; and
- `invz_run_spectra.m` exposes no backend, Ewald, grid-convention, offset, or Γ-policy knob.

That omission would make Ewald unreachable from the requested user-facing calculation. Worse, a
caller supplying precomputed Ewald BZ `Jnu/info` to `invz_spectra_qpath` would still obtain a
brute-force `invz_jq_path` dispersion. The corrected contract is:

1. `invz_spectra_map` and `invz_spectra_qpath` accept optional
   `opts.dipole`, `opts.ewald`, `opts.gridConvention`, `opts.gridOffset`, and
   `opts.gammaPolicy`.
2. Each copies those fields by **presence**, not by synthesizing defaults, into
   `invz_bz_couplings`. This preserves the old path when all fields are absent.
3. `invz_spectra_qpath` also forwards the same resolved `dipole` and `ewald` values to
   `invz_jq_path`. Grid convention/offset/Γ policy affect only the BZ medium, not the 1-D path.
4. `invz_bz_couplings`/`invz_jq_modes` return lossless provenance:
   `info.dipole.backend`, `info.dipole.ewald`, and `info.grid` (convention, offset, Γ policy,
   dimensions, exact q hash/schema). Precomputed `Jnu/info` must be internally complete.
5. If precomputed `Jnu/info` are supplied with an explicit conflicting request, error. If no
   explicit backend is requested, use the backend and controls recorded in `info` for the q-path
   calculation. Never silently mix backends.
6. `invz_run_spectra.m` exposes user knobs such as:

   ```matlab
   dipoleBackend = 'bruteforce';  % remains default until Gate F adoption
   ewaldOpts      = struct();      % populated only for the Ewald candidate
   gridConvention = '';            % empty => do not add a grid-policy field
   gridOffset     = [];            % empty => unchanged legacy grid path
   gammaPolicy    = '';            % empty => unchanged legacy Gamma filtering
   ```

   Before the production flip, the driver should add no grid-policy fields on its unchanged
   brute-force default path; candidate runs construct and pass them explicitly. After a successful
   Gate-F flip, these displayed defaults change to Ewald/half-open/frozen controls, while an
   explicit `dipoleBackend='bruteforce'` remains a supported diagnostic.

Required tests call both field-map and q-path routes with cold/warm caches, explicit Ewald,
explicit brute force, and precomputed couplings. They assert identical BZ/path provenance and the
documented mismatch errors.

---

## Blocker 6 — `invz_jq_path` Ewald Γ/metadata contract (design §2.3, Gate C)

### 6.1 Current mechanics (file:line walkthrough)

`invz_projected/invz_jq_path.m` (82 lines total):

- Public options (`:36-39`): `opts.dpRng` (30), `opts.cache` (true), `opts.snapfac` (2.5). No
  `dipole`/`ewald` option exists today.
- One coupling call for the **entire** path (`:43`): `[Jnu, info, Juni] =
  invz_jq_modes(ion, qpath, struct('dpRng', dpRng, 'cache', useCache))` — this line already
  determines the backend for every row once `invz_jq_modes` gains backend dispatch; it simply
  isn't forwarding `dipole`/`ewald` yet.
- Trust radius (`:47`): `ksnap = snapfac * 2*pi / (dpRng * min(vecnorm(ion.a,2,2)))` — **derived
  from `dpRng`**, i.e. from the brute-force real-space cutoff. This is the exact "`dpRng`-
  dependent Γ snap logic" the design (`invzp_ewald_design.md:132-133`) names for replacement.
- Per-point loop (`:50-77`): for each Γ-equivalent path point (`G=round(q)`,
  `invz_is_gamma_equiv`, `:51-52`), compute the Cartesian offset `k=(q-G)*Brec` (`:53`); **skip
  the whole snap** if `norm(k)>=ksnap` (`:54`, i.e. only points strictly inside the trust radius
  are touched); if `norm(k)<1e-12` (exactly at `G`), derive a **local path direction** from the
  adjacent path point (`:55-62`); else use the actual offset direction. `kz2=(k(3)/norm(k))^2`
  (`:63`).
- **The hazard** — lazy `Greg` build, once per call (`:64-71`):
  ```matlab
  dip0 = MF_dipole([0 0 0], dpRng, ion.a, ion.tau, info.geomD);
  ex0  = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau, info.geomX);
  Greg = -squeeze(C.gfac*dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));
  ```
  This is `invz_jq_path.m`'s **only** direct `MF_dipole` call, and it consumes `info.geomD`
  (produced by whatever backend `invz_jq_modes` used at `:43`).
- Directional broadcast (`:72-76`): `Jm = Greg + C.gfac*(4*pi/ion.Vc)*(1/3 - kz2)`, Hermitized,
  eigendecomposed into `Jnu(iq,:)`/`Juni(iq)`, `snapped(iq)=true`. This is **exactly** the
  derivation note's "existing convention" form (`invzp_ewald_derivation.md:73`):
  `J_cc(q→0,q̂)[existing] = J_reg′,cc + gfac(4π/Vc)(1/3−q̂_c²)` — confirming `Greg` here **is**
  `J_reg′,cc` (no Lorentz folded in; the `1/3` term supplies it isotropically at `kz2=1/3`,
  matching `invz_jq_modes`'s own `+lorz` at exact Γ, `invz_jq_modes.m:86,93`).
- Outputs (`:79-81`): `P.Jnu`, `P.Juni`, `P.snapped`, `P.ksnap`, `P.s`, `P.s_cart`.

**Callers/coverage:** `invz_projected/invz_spectra_qpath.m:91` (production driver) and
`invz_projected/tests/test_invz_spectra_qpath.m` (its only test coverage — no dedicated
`test_invz_jq_path.m` exists). Neither reads `P.ksnap`'s numeric value; only `.snapped`
(boolean) is asserted (`test_invz_spectra_qpath.m:31,152`).

### 6.2 The hazard, concretized

`MF_dipole(q,N,a,tau,geom)` (`MF_dipole.m:2,16-18`) trusts a 5th-argument `geom` **blindly**
whenever `nargin>=5 && ~isempty(geom)` — it indexes `geom.r{nt,mt}`, `geom.Tf{nt,mt}`, `geom.b`,
`geom.ntau` (`MF_dipole.m:21,26-27`) with **no schema/fingerprint check at all**. Design §2.1
(`invzp_ewald_design.md:95-96`) requires `invz_dipole_ewald`'s own `geom` to carry "validated
lattice metadata... a schema version" and to error on an incompatible `geom` — but that
self-check lives inside `invz_dipole_ewald`, not inside `MF_dipole`. If an Ewald-backend
`invz_jq_modes` call ever populates `info.geomD` with an Ewald-shaped struct (same field name,
different shape/schema), `invz_jq_path.m:68`'s `MF_dipole(...,info.geomD)` call would try to
dereference Ewald-geom fields as MF_dipole's own — most likely an "unrecognized field" crash
with **no diagnostic pointing at the real cause** (backend mismatch), and triggered **only** for
q-paths that happen to pass near/through a Γ-equivalent point — a latent, intermittent, confusing
failure exactly matching the design's phrase "ambiguous-geometry hazard"
(`invzp_ewald_design.md:137-138`).

### 6.3 Selected contract

**Foundational assumption** (owned by `invz_jq_modes`'s own backend-dispatch work, design §2.3
bullet 1 — not re-litigated here, but both this blocker and Blocker 5 need a consistent shape
from it): `info` carries an explicit, self-describing backend tag, e.g.
`info.dipole = struct('backend', 'bruteforce'|'ewald', ...)`. I use `info.dipole.backend` below;
any single unambiguous field works as long as Blockers 5 and 6 agree on its name.

**Selected contract — eliminate `invz_jq_path`'s `MF_dipole` call entirely:**

`invz_jq_modes`, for **both** backends, additionally exports a backend-agnostic
`info.Jpath_base_cc` `[4x4]` — exactly the quantity `invz_jq_path` currently rebuilds locally as
`Greg` (§6.1: `-gfac*dip0(3,3,:,:) + sign(J12)*ex0(3,3,:,:)`, the "existing convention"
`J_reg′,cc`, no Lorentz folded in):

- For **bruteforce**, this is a trivial hoist: `invz_jq_modes.m:78,93` already computes the full
  `dip0` `[3,3,4,4]` matrix and its `Jcc0d` slice in the priming call; today only the **scalar**
  uniform-mode projection `info.Jcc0_dipole` (`:97`) is exported. Add one line exporting the
  **full** `4x4` cc-matrix alongside it (same inputs already in scope, negligible cost).
- For **Ewald**, raw `-gfac*dip_reg_cc+Jex_cc` already contains one Lorentz broadcast relative to
  the legacy `Greg`. Therefore the export is explicitly

  ```text
  info.Jpath_base_cc =
      -gfac*dip_reg_cc(0) + Jex_cc(0) - lorz*ones(4).
  ```

  This subtraction is metadata normalization for the legacy q-path formula, not a primitive
  surface/Lorentz correction. The exact-Γ production matrix is exported separately as
  `info.Jgamma_cc = -gfac*dip_reg_cc(0)+Jex_cc(0)` for Ewald and
  `info.Jgamma_cc = Greg+lorz*ones(4)` for brute force.

`invz_jq_path` then replaces `:64-71` with simply:

```matlab
if isempty(Greg), Greg = info.Jpath_base_cc; end
```

— no primitive call, no `geom` field of any kind consumed, for **either** backend. The hazard is
**deleted**, not guarded: `invz_jq_path` would no longer read `info.geomD`/`info.geomX` at all.
The downstream broadcast formula (`:72`) is **unchanged** and stays backend-agnostic, because
`info.Jpath_base_cc` is defined to already be in the existing-convention normalization for both
backends by construction. Regression impact for bruteforce is bit-identical: `MF_dipole`/
`exchange` are deterministic pure functions of `(q=0, dpRng, a, tau[, geom])`, so "freshly
recomputed via a 2nd `MF_dipole` call" and "hoisted from the priming call already inside
`invz_jq_modes`" produce numerically identical values — this should be pinned by a bit-identity
regression test (extending the pattern already used by `test_invz_dipole_geometry_reuse.m`)
comparing `invz_jq_path`'s old locally-rebuilt `Greg` against the new `info.Jpath_base_cc` before
retiring the old code path.

No backend-dispatched geometry alternate remains in the approved design. Centralizing both
`Jpath_base_cc` and `Jgamma_cc` in `invz_jq_modes` removes the ambiguous-geometry hazard and
single-sources the Lorentz normalization.

### 6.4 `opts`/behavior changes to `invz_jq_path` itself

- New optional `opts.dipole` (absent means the current implicit brute-force default;
  explicit `'bruteforce'` or `'ewald'`) and conditionally present `opts.ewald` (a complete struct
  required with explicit `'ewald'`). Presence is preserved when forwarding to the single
  `invz_jq_modes(...)` call at `:43`; no empty Ewald field is synthesized.
- **`ksnap` policy becomes backend-conditional** (Gate C, `invzp_ewald_design.md:330-332`:
  *"remove the current brute-force `dpRng` trust-radius behavior from the Ewald q-path branch.
  Nonzero q must be evaluated directly"*):
  - `bruteforce` (default): **unchanged** — `ksnap` formula at `:47` and the
    `0 < norm(k) < ksnap` snap band exactly as today.
  - `ewald`: no trust radius applies. The trigger that currently gates entry into the whole
    directional-limit block (`:54`, `if norm(k) >= ksnap, continue`) must become
    backend-conditional: under `ewald`, only `norm(k) < 1e-12` (machine-exact Γ) should ever
    proceed; every genuinely nonzero offset, however small, is evaluated directly by the raw
    `invz_jq_modes` call already made at `:43`. **Implementation note:** simply setting
    `ksnap=0` and reusing the existing `>=` comparison does **not** work — `norm(k)>=0` is true
    even at `norm(k)==0`, which would skip the exact-Γ directional-prescription branch too. The
    **trigger condition itself** must branch on backend, not just its numeric threshold:
    ```matlab
    if strcmp(backend, 'ewald'), trigger = norm(k) < 1e-12;
    else,                         trigger = norm(k) < ksnap;   end
    if ~trigger, continue; end
    ```
  - `opts.snapfac` becomes a documented no-op under `ewald` (not silently misleading);
    `opts.dpRng` likewise stops feeding `ksnap` and is not forwarded to the Ewald branch of
    `invz_jq_modes`. It remains accepted only so common driver option structs need not be split;
    `info.dpRng`/`P.ksnap` are `NaN` under Ewald and the Ewald controls are reported instead.
  - `P.ksnap` (reported output) is `NaN` under `ewald` rather than `0`, so a reader
    never mistakes "no trust radius" for "a physically-tuned radius of zero." `bruteforce`
    unchanged.
- `P.Jnu`/`P.Juni`/`P.snapped`/`P.s`/`P.s_cart`: unchanged shape/meaning. `P.snapped`'s
  semantics under `ewald` narrow correctly to "this row was exactly Γ" (vs. "within the `dpRng`
  trust radius" for bruteforce) — no consumer break: `test_invz_spectra_qpath.m:150-152,31`
  only checks `.snapped` at path points that are exact-integer Γ-equivalent endpoints either way.
- `opts.odd` interaction: `invz_jq_path` does not read or forward `opts.odd` today, and this
  proposal does not add that; the Ewald+ODD unsupported-combination guard belongs where
  `opts.odd` is actually dispatched (`invz_jq_modes.m:35`, design §1.1) — no new guard needed
  in this file.

### 6.5 What `info` must identify (summary)

- `info.dipole.backend` — string, `'bruteforce'`\|`'ewald'` — consumed by `invz_jq_path` **only**
  to pick the `ksnap` trigger policy (§6.4); **not** consumed for the `Greg`/broadcast math
  (backend-agnostic given `info.Jpath_base_cc`).
- `info.Jpath_base_cc` `[4x4]` — **new**, backend-agnostic legacy-normalized q-path base, with the
  explicit Ewald `-lorz*ones(4)` normalization in §6.3.
- `info.Jgamma_cc` `[4x4]` — **new**, exact-Γ isotropic production matrix; never substituted for
  `Jpath_base_cc` without changing the analytic q-path formula.
- `info.dipole.ewald` and `info.dipole.q_reduction` — exact controls/provenance needed by caches,
  precomputed-input validation, and the user-facing q-path.

### 6.6 Files / functions / lines touched by the selected contract

- `invz_projected/invz_jq_path.m` — replace `:47` (`ksnap` formula → backend-conditional trigger,
  §6.4), `:54` (trigger condition), `:64-71` (`Greg` build → `info.Jpath_base_cc` read),
  `:36-39` (new `opts.dipole`/`opts.ewald`), `:43` (forward them
  into the `invz_jq_modes` call), `:79` (`P.ksnap` → `NaN` under ewald).
- `invz_projected/invz_jq_modes.m` — backend dispatch/cache validation plus
  `info.dipole.backend`, exact Ewald/q-reduction provenance, `info.Jpath_base_cc`, and
  `info.Jgamma_cc` alongside its existing Γ scalars.
- `MF_dipole.m` — **untouched**; the selected contract stops calling it from `invz_jq_path`.
- No change to `exchange(...)`/`info.geomX` handling — `exchange` is not part of the dipolar
  Ewald primitive (design §1.1 scope is dipole-only); it is either left exactly as-is or folded
  into the same metadata construction (`Jpath_base_cc` and `Jgamma_cc` both include the
  `sign(J12)*ex0(3,3,:,:)` term). No new exchange backend is introduced.

### 6.7 Risks

1. **Hermitization order.** `invz_jq_path` Hermitizes only the **sum** `Jm` after adding the
   broadcast term (`:73`, `Jm=(Jm+Jm')/2`), never `Greg` alone. `info.Jpath_base_cc` must preserve
   this exactly (not pre-Hermitized a second time inside `invz_jq_modes` in a way that changes
   the final bit pattern vs. today's construction order) — pin with a bit-identity regression
   test, not an assumption.
2. **Two Γ-removal concepts must stay separate here too (design §4.1).** `invz_jq_path`'s
   directional guard is a **path-local, single-point** δ-neighborhood prescription (1-D
   approach direction to one q-path point at Γ); Blocker 5's `gammaPolicy` is a **3-D BZ-volume
   quadrature** row-inclusion policy over a discrete grid. No shared flag: `opts.gammaPolicy`
   must never be threaded into `invz_jq_path`, and `invz_jq_path`'s Γ-direction logic must never
   be reused to justify a BZ quadrature choice, or vice versa.
3. **Gate-C arithmetic is fixed but still prospective.** Calibration selected “Ewald adds zero
   extra Lorentz at Γ”; the metadata formulas in §6.3 are therefore no longer placeholders.
   Formal Gate C must reproduce the result after the preregistration freeze or stop.
4. **Dispatch/cache/ODD behavior is part of the implementation contract.** `invz_jq_modes` must
   validate `opts.dipole`/`opts.ewald`, reject active ODD+Ewald, use lossless backend-separated
   cache keys and exact-payload validation, canonicalize q before reciprocal enumeration, restore
   gauge covariance, and export the provenance fields listed in §6.5. These are not deferred to an
   unspecified later design.

---

## Summary / remaining Step-3 details

- The Gate-C convention is no longer unresolved: calibration selects
  `conducting_k0_omitted`, no extra Γ Lorentz block, and the exact `Jpath_base_cc`/
  `Jgamma_cc` formulas in §6.3. Formal acceptance is still prospective and must rerun after the
  freeze.
- Backend dispatch, cache separation, ODD rejection, canonical q reduction/gauge restoration,
  metadata, and all three user-facing spectra forwarders are part of the closed implementation
  contract rather than unspecified dependencies.
- `invz_dipole_ewald.m` does not yet exist, so the integration claims remain static analysis.
  That is appropriate before preregistration and implementation.
- `invz_phase1_qgrid.m` remains cubic-only. Candidate Gate-D/E configurations are cubic; any
  non-cubic request through the new branch must error until separately supported.
- Step 3 still must freeze exact q/K/ray/basis samples, norms and absolute floors, α/cutoff
  comparison axes, hard resource caps, the full N×offset×policy run matrix, Jensen field/continuity
  criteria, expected singular-Γ exclusions, retained oracle artifacts, and the exact `Bc_PM`
  re-freeze outputs. Those are deliberately prospective preregistration details, not open API or
  mathematical design choices.
