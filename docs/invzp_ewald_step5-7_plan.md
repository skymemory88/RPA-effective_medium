# Ewald Step-5 Integration Implementation Plan

**Status:** reviewed and corrected 2026-07-24.

**Goal:** Integrate the additive Ewald primitive through the projected production path while the
production default remains `bruteforce`. Step 5 completes the remaining integration-level Gate-C
checks (full-Cartesian C4, demag C6, and cache/provenance C7) and preserves the legacy numerical path
bit-for-bit.

**Authority order:** `docs/invzp_ewald_prereg.md` (FROZEN, including §12 Errata E1) >
`docs/invzp_ewald_design.md` > `docs/invzp_ewald_integration_map.md` > this plan.

**Confirmed phase boundary:** Step 5 may create and parity-test
`invz_projected/tests/invz_anchor_couplings.m`, but it must not run or create the Gate-E Ewald
physics-anchor wrappers (`test_invz_sigma_crit_ewald.m`, Tc/Bc/Jensen tests). Gate D remains Step 6;
Gate E and any default flip remain Step 7.

## Entry state and non-negotiable constraints

- The Step-5 base is `a58bdb4`, which commits the first manifest/design closeout over `5701065`.
  The last independently verified pre-closeout projected-suite result is
  **269 passed / 0 failed / 23 incomplete**; Task 1 re-establishes the post-closeout count before
  production integration.
- The worktree is dirty. In particular, `invz_projected/invz_run_spectra.m` contains the user's
  unrelated parameter edits. Preserve those edits exactly. Never stage that whole file while they
  remain mixed with Step-5 work.
- Commit `a58bdb4` contains the first Step-4 closeout edits in `invz_dipole_ewald.m`,
  `invz_projected/tests/test_invz_dipole_ewald.m`, and `docs/invzp_ewald_design.md`. Task 1 audits
  that committed inventory against the stricter literal preregistration rule and adds only a
  follow-up if the source-to-manifest audit proves something is still omitted.
- No default change, no ODD/tensor-path adoption, no Jensen-field inspection, and no claim that the
  user-facing masking symptom is fixed.
- Do not modify `MF_dipole.m`, `exchange.m`, `invz_common/*`, `invz_tensor/*`, the ODD solver
  internals, or the frozen preregistration.
- After every production task, run the portable legacy regression plus the task's focused tests
  before committing. Run the full projected suite at the major integration checkpoints and closure.
- Stage only explicit paths/hunks. Do not add a `Co-Authored-By` trailer.

## Frozen behavior contracts

### Legacy regression

For both an absent `opts.dipole` and explicit `opts.dipole='bruteforce'`:

- `Jnu` and `Juni` from `invz_jq_modes` are `isequaln` to the frozen legacy reference.
- Every legacy `info` field exists and is individually `isequaln` to the reference. New fields are
  removed before comparing the legacy field set; the test must fail if any legacy field is missing.
- New additive fields (`info.dipole`, `info.Jpath_base_cc`, `info.Jgamma_cc`, and BZ-layer
  `info.grid`) are tested separately.
- `invz_jq_path`'s pre-existing output fields are `isequaln` to the frozen legacy reference in the
  brute-force branch. Its new provenance field, if added, is tested separately.

Use a retained, test-only MATLAB reference copied from the pre-Step-5 algorithm and evaluated at
runtime. Do **not** commit a binary `.mat` numerical baseline: eigensolver/BLAS bit patterns can be
platform-dependent, whereas a runtime reference tests exact equivalence on the active platform.

### Γ metadata

Let `lorz = gfac*4*pi/(3*Vc)` and let `Jex_cc(0)` include `sign(J12)`:

```text
Jgamma_cc:
  bruteforce = -gfac*dip_sphere_cc(0) + Jex_cc(0) + lorz*ones(4)
  ewald      = -gfac*dip_reg_cc(0)    + Jex_cc(0)

Jpath_base_cc:
  bruteforce = -gfac*dip_sphere_cc(0) + Jex_cc(0)
  ewald      = -gfac*dip_reg_cc(0)    + Jex_cc(0) - lorz*ones(4)

J(q→0,qhat) =
  Jpath_base_cc + gfac*(4*pi/Vc)*(1/3 - qhat_c^2)*ones(4)
```

Preserve the existing arithmetic and Hermitization order. In particular, `invz_jq_path` Hermitizes
the final reconstructed matrix, not `Jpath_base_cc` on its own.

### Options and Ewald q-path

- `opts.dipole` is absent or exactly `'bruteforce'|'ewald'`.
- `opts.ewald` is present only for Ewald and has exactly
  `{alpha,r_cut,g_cut,boundary}`. The only accepted boundary is
  `'conducting_k0_omitted'`.
- Top-level `opts.alpha` remains the demagnetization aspect ratio and is never confused with
  `opts.ewald.alpha`.
- Active ODD (`opts.odd` neither absent nor exactly `false`) plus Ewald errors. ODD plus legacy
  brute force is unchanged.
- Under Ewald, `dpRng` does not affect calculation, provenance, or cache identity;
  `info.dpRng=NaN`. A common caller may still contain `dpRng`, but `invz_jq_path` does not forward it
  into the Ewald `invz_jq_modes` call.
- Under Ewald, every genuinely nonzero q is evaluated directly. Only
  `norm(k)<1e-12` enters the exact-Γ local-direction prescription, and `P.ksnap=NaN`. The trigger
  branches on backend; setting `ksnap=0` is not sufficient.
- `gammaPolicy` is BZ quadrature provenance only and is never forwarded to `invz_jq_path`.

### Cache and provenance

Use a new cache schema/prefix; do not extend an old `jq4_` file in place:

```text
jq5_bruteforce_...
jq5_ewald_...
schema = invz_jq_modes/v5
```

The literal backend in the filename guarantees backend-separated candidate keys. A filename digest
may be compact, but a cache hit is accepted only after exact `isequaln` validation of a structured
payload containing:

- exact `qvec`, lattice, basis, `Vc`, exchange `J12`, `gfac`;
- demag and top-level aspect ratio;
- backend and exact Ewald controls (or an empty canonical brute-force value);
- brute-force `dpRng` or the Ewald `NaN` sentinel;
- exact BZ cache context (route, grid dimensions, convention, offset, Γ policy, q digest) or a
  canonical direct-call context;
- schema version and required output fields/shapes.

The BZ layer passes a private, validated `jqOpts.cacheContext` into `invz_jq_modes`. This is required:
`qvec` alone does not identify convention/offset/Γ-policy provenance, and the frozen preregistration
requires those fields to be validated. `info.grid` remains owned by `invz_bz_couplings`, not
`invz_jq_modes`.

`info.grid.qhash` is a SHA-256 (or equivalently collision-resistant) digest over a canonical
byte-level serialization of q-array class, shape, and exact IEEE-754 contents. Do not reuse the
existing weak single-precision weighted-sum `hash_vec` as provenance. Cache safety still rests on
exact stored-`qvec` validation, not on digest collision resistance.

## Files in scope

Production modifications:

- `invz_projected/invz_jq_modes.m`
- `invz_projected/invz_bz_couplings.m`
- `invz_projected/invz_jq_path.m`
- `invz_projected/invz_spectra_map.m`
- `invz_projected/invz_spectra_qpath.m`
- `invz_projected/invz_run_spectra.m`

Step-4 closeout modifications:

- `invz_dipole_ewald.m`
- `invz_projected/tests/test_invz_dipole_ewald.m`
- dated status addenda only in `docs/invzp_ewald_design.md` and
  `docs/invzp_ewald_integration_map.md`

New test/support files:

- `invz_projected/tests/invz_legacy_coupling_reference.m`
- `invz_projected/tests/test_invz_ewald_baselines.m`
- `invz_projected/tests/test_invz_jq_modes_ewald.m`
- `invz_projected/tests/test_invz_bz_couplings_ewald.m`
- `invz_projected/tests/invz_anchor_couplings.m`
- `invz_projected/tests/test_invz_jq_path_ewald.m`
- `invz_projected/tests/test_invz_spectra_forward_ewald.m`
- `invz_projected/tests/test_invz_ewald_integration.m`
- optionally `invz_projected/tests/test_invz_run_spectra_wiring.m` for the interactive script's
  static wiring contract

---

## Task 1 — Finish Step-4 closeout and freeze a portable legacy reference

**Purpose:** Close the strict resource-accounting/document-status issues before the first production
integration edit.

- [ ] Review commit `a58bdb4`'s manifest patch against the actual live arrays in
  `local_build_geom`, `local_gab`, `local_boxmin_dist`, and `local_assemble`. The committed expected-
  name test is evidence, not proof of completeness, because it can omit the same row as the
  implementation.
- [ ] Make the preflight inventory cover every planned retained/temporary size-dependent array,
  including simultaneously live `ndgrid` outputs plus concatenated meshes, real-space
  `local_gab` work arrays, reciprocal retain masks and retained `Ghkl/Gcart`, all box-minimizer
  arrays, and the complex `w.*kK`/matrix-product work in q assembly. Include fixed-size numeric
  scratch too, or document and test an authority-supported reason for excluding it. Keep the
  estimate conservative, cube-bound based, pre-allocation, and exactly
  `1.25*sum([manifest.bytes])`; do not change any cap.
- [ ] Strengthen `test_manifest_names_are_complete` so the expected inventory is derived from a
  maintained source-to-manifest table, checks representative class/complexity/shape bounds, and
  fails on a missing or duplicate row. Retain the existing cap tests.
- [ ] Verify the committed dated status addendum in `docs/invzp_ewald_design.md`. Add a matching
  dated status addendum to `docs/invzp_ewald_integration_map.md`, whose historical body still says
  the primitive/preregistration do not exist. Preserve the historical text below each addendum.
- [ ] Create `invz_legacy_coupling_reference.m`, a cache-free test-only copy of the pre-Step-5
  brute-force `invz_jq_modes` and `invz_jq_path` arithmetic. Keep call order and Hermitization
  order exact. It may call unchanged `MF_dipole`, `exchange`, and other production primitives, but
  it must not call the production `invz_jq_modes` or `invz_jq_path` being modified.
- [ ] In `test_invz_ewald_baselines.m`, prove the reference equals current production at a small
  deterministic q array containing ordinary, near-Γ, exact-Γ, and Γ-equivalent points. Compare
  `Jnu`, `Juni`, every `info` field, and every pre-existing `P` field with `isequaln`.
- [ ] Run the three Step-4 Ewald test files and the full projected suite. Commit only the reviewed
  closeout/docs/reference paths.

**Commit:** `fix(invz): close Ewald manifest and freeze portable legacy integration reference`

---

## Task 2 — Atomic `invz_jq_modes` dispatch, Γ metadata, and v5 cache

Backend dispatch, metadata construction, and cache replacement land in one commit. No intermediate
revision may evaluate Ewald through `jq4_`.

- [ ] Add failing tests for absent-backend and explicit-brute reference parity, including all three
  numerical outputs and every legacy `info` field.
- [ ] Add failing validation tests for unknown/non-scalar backend, `opts.ewald` without Ewald,
  Ewald controls that are missing or extra, invalid boundary, brute force plus `opts.ewald`, and
  active ODD plus Ewald. Pin stable namespaced error identifiers.
- [ ] Resolve the backend before the ODD diversion. Permit `opts.odd=false` with Ewald; reject every
  active ODD form with Ewald. Leave the existing ODD/brute path and its caches untouched.
- [ ] Preserve the brute-force code path's operation order. Hoist `Jpath_base_cc` from the existing
  priming values without pre-Hermitizing it, construct `Jgamma_cc` by the frozen formula, and append
  additive metadata only after the legacy fields have been produced.
- [ ] For Ewald, build/reuse one Ewald geometry across q, keep `exchange` unchanged, and add no
  extra Lorentz block at Γ. Construct the existing `Jcc0_dipole`, `Jaa0_dipole`, `Jcc0`, `Jaa0`,
  and `Jshape_cc` with the same observable semantics as the brute branch. In particular, Ewald's
  regular tensor already contains the isotropic term; do not add the brute-force `+lorz` to either
  cc or aa.
- [ ] Export
  `info.dipole = struct('backend',...,'ewald',...,'q_reduction',...,'primitive_schema',...)`.
  Brute force uses a canonical empty Ewald value and a documented legacy q convention; Ewald uses
  the exact validated controls and the primitive fingerprint's q-reduction/schema values.
- [ ] Implement `jq5_<backend>_...` with structured `cacheMeta` and exact payload validation per the
  global cache contract. Missing, legacy, malformed, or mismatched payloads are misses and are
  recomputed. Absent backend and explicit brute share the same canonical cache identity. Cold/warm
  calls for each backend must be `isequaln` across `Jnu`, whole `info`, and `Juni`.
- [ ] Under Ewald set `info.dpRng=NaN`, ignore `dpRng` in the key, and never expose Ewald geometry
  through `info.geomD`. `info.geomX` may remain the exchange geometry. The brute-force legacy
  `geomD/geomX` fields remain bit-identical.
- [ ] Run the portable legacy regression, the new jq-mode tests, all existing jq/ODD/cache tests,
  and the full projected suite.

**Commit:** `feat(invz): add atomic Ewald jq dispatch metadata and v5 cache`

---

## Task 3 — Integration Gate-C4 full Cartesian and Gate-C6 demag

### Gate-C4

`invz_jq_modes` exposes only the cc channel, so an all-nine-component test cannot claim to extract a
full Cartesian tensor from that API. Use the primitive plus caller-level coupling normalization for
all components, then cross-check the cc slice through `invz_jq_modes`.

- [ ] On the five frozen rays and `s∈{±1e-3,±1e-4}`, use the exact finite-s isolated projector
  `P0_ab=(4*pi/Vc)*qhat_a*qhat_b*exp(-|q|^2/(4*alpha^2))` for all nine components and all
  sublattice pairs. Apply the caller's `-gfac` conversion so the coupling sign is also tested.
- [ ] For every `(a,b)`, form the Ewald caller-level base
  `Jbase_ab=-gfac*dip_reg_ab(0)+Jex_ab(0)-gfac*(4*pi/(3*Vc))*delta_ab*ones(4)` and verify the
  complete limit reconstruction
  `Jbase_ab+gfac*(4*pi/Vc)*(delta_ab/3-qhat_a*qhat_b)*ones(4)`.
- [ ] Separately verify the cc member is exactly the `info.Jpath_base_cc` convention and its
  reconstructed limit agrees with the production q-path formula. Do not compare a finite-s
  projector to the no-Gaussian small-q limit at `M_id`; prereg §12 E1 explicitly forbids that
  mistake.

### Gate-C6

- [ ] For each backend, test exactly three caller cases: off (`demag=0`), sphere (`demag=1`,
  top-level aspect `alpha=1`), and c-axis needle (`demag=1`, top-level aspect `alpha=0`).
- [ ] Within each backend, require `Jnu`, `Jcc0`, `Sigma_c`, and
  `invz_critical_T0field(...)` to be identical across the three shapes. Test `Jshape_cc` and the
  demag-aware `Jaa0` against the existing analytic `ellipsoid_demagn` formulas.
- [ ] Test backend agreement only where a frozen numerical tolerance exists; do not invent a new
  brute/Ewald absolute-agreement tolerance for Gate C6. Confirm that the primitive options contain
  no demag/surface control.
- [ ] Run and commit the tests. No production edit is expected in this task.

**Commit:** `test(invz): complete integration Gate-C4 and Gate-C6`

---

## Task 4 — `invz_bz_couplings` grid/backend threading and exact grid provenance

- [ ] Add failing tests proving the no-grid-field route reproduces the portable legacy reference
  exactly and keeps `info.grid` absent.
- [ ] Add presence-routing tests for each of `gridConvention`, `gridOffset`, and `gammaPolicy`.
  Any one activates `invz_phase1_qgrid`; Ewald alone does not. Explicit
  `gridConvention='legacy_inclusive'` still activates the wrapped new route and is not claimed
  bit-identical to the absent-field route.
- [ ] Preserve the old `qVec_generator` and Γ-drop statements exactly in the absent-grid route.
  Forward `dipole` and `ewald` by presence through a `jqOpts` struct.
- [ ] Reject non-cubic grids on the new route with
  `invz:bzCouplingsAnisotropicHalfopen`. Let `invz_phase1_qgrid` retain authority over
  convention/offset/policy value validation.
- [ ] Construct complete `info.grid` only on the new route. At minimum include schema, convention,
  logical `1x3` offset, Γ policy, requested dimensions, nominal rows, retained rows, `n_gamma`, and
  the exact-byte q digest.
- [ ] Pass a canonical private `cacheContext` on both BZ routes:
  `kind='legacy_bz'` plus exact grid dimensions/q digest and explicit absent-policy sentinels for
  the legacy route; `kind='phase1_qgrid'` plus all resolved grid provenance for the new route.
  `invz_jq_modes` validates and stores this context but does not export `info.grid`.
- [ ] Test `P_complete` and `P_drop` row counts/provenance, exact cache-context separation, and
  absent/explicit-brute legacy parity. Run the portable regression and full suite.

**Commit:** `feat(invz): thread Ewald BZ options and exact grid provenance`

---

## Task 5 — Cheap `invz_anchor_couplings` helper only

- [ ] Create the thin helper exactly as a presence-preserving wrapper over
  `invz_bz_couplings(ion,...)`. Defaults are `dpRng=30`, `cache=true`, and no
  backend/grid-policy fields.
- [ ] On a cheap grid (`N=4` or `N=6`, `cache=false`), compare default helper `Jnu(:)` and every
  legacy `info` field bit-for-bit against the actual inline legacy
  `qVec_generator` + Γ-drop + `invz_jq_modes` construction.
- [ ] Add a cheap injection test with Ewald + half-open + `[0 0 0]` +
  `P_complete`, verifying complete backend/grid provenance. A second mechanical call may verify
  `P_drop`; do not calculate `Sigma_c`, Tc, Bc, Jensen fields, or any anchor inequality.
- [ ] Run and commit.

**Commit:** `feat(invz): add parity-tested Gate-E coupling injection helper`

---

## Task 6 — `invz_jq_path` metadata base and Ewald Γ behavior

- [ ] Add failing brute-force tests comparing every pre-existing `P` field to the portable old-path
  reference, including paths with ordinary points, inside/outside the old snap band, exact Γ, and
  Γ-equivalent endpoints.
- [ ] Add an explicit regression proving the old locally rebuilt `Greg` is `isequaln` to the new
  brute-force `info.Jpath_base_cc` before deleting the local `MF_dipole`/`exchange` reconstruction.
- [ ] Forward `dipole`/`ewald` by presence. Preserve the old jq call shape for absent backend;
  explicit brute adds only `dipole`. For Ewald, forward `cache`, `dipole`, and `ewald`, but not
  `dpRng`.
- [ ] Replace the local Γ base construction with `info.Jpath_base_cc`; remove every read of
  `info.geomD`, `info.geomX`, `MF_dipole`, and `exchange` from this file.
- [ ] Branch the trigger: brute force retains the existing trust radius exactly; Ewald triggers
  only for `norm(k)<1e-12`. Set `P.ksnap=NaN` under Ewald and retain `snapfac` as a documented no-op.
- [ ] Export additive `P.dipole=info.dipole` so downstream integration tests and users can verify
  the backend actually used. Preserve all existing `P` fields and meanings.
- [ ] Test a tiny nonzero Ewald q below the former brute trust radius against direct
  `invz_jq_modes` output (`P.snapped=false`), and exact Γ with local direction
  (`P.snapped=true`). Run the portable regression and full suite.

**Commit:** `feat(invz): use jq metadata for backend-correct Gamma paths`

---

## Task 7 — Spectra forwarding and strict precomputed resolution

Modify both spectra drivers symmetrically and document the new optional
`dipole/ewald/gridConvention/gridOffset/gammaPolicy/cache` fields.

- [ ] Reject a partial precomputed pair: exactly one of `opts.Jnu` and `opts.info` is an error.
- [ ] When couplings are computed, copy backend/grid/cache fields by presence into
  `invz_bz_couplings`. In q-path spectra, resolve one backend and pass it to `invz_jq_path`;
  never pass grid convention/offset/Γ policy to the path.
- [ ] When precomputed couplings have no explicit backend/grid request:
  - complete `info.dipole` provenance is inherited by the q-path;
  - provenance-less legacy synthetic inputs remain accepted and resolve to brute force;
  - existing tests whose `info` contains only fields such as `Jcc0` remain unchanged.
- [ ] When any explicit backend request is present, require complete `info.dipole` provenance for
  precomputed inputs and compare backend plus Ewald controls exactly. Missing provenance errors;
  disagreement errors with `invz:spectraBackendConflict`.
- [ ] When any explicit grid-policy request is present with precomputed inputs, require complete
  `info.grid`, resolve omitted grid-policy members using the same BZ defaults, and compare exact
  convention/offset/Γ policy/dimensions. Missing provenance and conflicts use stable, distinct
  namespaced errors.
- [ ] Validate explicit options even when precomputed inputs bypass `invz_jq_modes`; invalid
  Ewald controls must not escape checking simply because no lattice sum is run.
- [ ] In q-path spectra, require `isequaln(P.dipole,info.dipole)` whenever complete BZ provenance
  exists, and expose `S.path_dipole=P.dipole`. A mismatch is an error, never a warning.
- [ ] Add cheap tests for computed, precomputed-complete, and legacy-provenance-less routes; matching
  and conflicting explicit requests; both Γ policies; and no Jensen/symptom assertions.
- [ ] Run the portable regression, existing spectra tests, focused new tests, and full suite.

**Commit:** `feat(invz): enforce spectra backend provenance and precomputed conflicts`

---

## Task 8 — `invz_run_spectra` user knobs without absorbing user edits

- [ ] Record the pre-task diff of `invz_projected/invz_run_spectra.m` and treat it as protected.
- [ ] Add editable knobs:

  ```matlab
  dipoleBackend = 'bruteforce';
  ewaldOpts      = struct();
  gridConvention = '';
  gridOffset     = [];
  gammaPolicy    = '';
  ```

- [ ] Build one coupling-option struct and merge it into every map/q-path call. Always pass the
  explicit backend knob. Pass `ewald` only when Ewald is selected. Add each grid-policy field only
  when its knob is nonempty; the default script therefore activates no new grid route.
- [ ] Cover scalar q-path, multi-field q-path loop, and field-map calls. Do not change the user's
  current T/field/frequency/eta/angle/MF/plot settings.
- [ ] Because the script is interactive and its defaults are expensive, do not use a full default
  script run as the wiring test. Use `checkcode`, a small static source-contract test if needed, and
  the callable map/q-path integration tests from Tasks 7 and 9.
- [ ] Stage only the integration lines using an index-only patch or reviewed interactive hunk
  selection. Inspect `git diff --cached` and prove none of the user's pre-existing parameter edits
  is staged. If clean separation is impossible, leave this task uncommitted and request a user
  decision; never co-mingle the edits silently.

**Commit:** `feat(invz): expose opt-in Ewald spectra driver controls`

---

## Task 9 — Gate-C7 end-to-end cache/provenance matrix

Use small grids, one stable high-field point, a minimal frequency vector, and a short q-path. This
is a mechanical integration test, not a physics anchor.

- [ ] Exercise field-map and q-path routes for explicit brute force and explicit Ewald, with both
  `P_complete` and `P_drop` on the Ewald BZ route.
- [ ] For each backend, exercise cold then warm cache calls and require `isequaln` on all numerical
  outputs and whole metadata. Clean up only the exact test-created `jq5_*` files; never delete the
  cache directory or unrelated caches.
- [ ] Prove brute/Ewald filenames have distinct literal backend prefixes and that changing each of
  q, lattice, basis, exchange, demag/aspect, Ewald control, grid convention, offset, Γ policy, and
  schema prevents a false cache hit. A deliberately colliding compact digest, if practical, must
  still be rejected by exact `cacheMeta`.
- [ ] Exercise computed couplings, complete precomputed couplings, and legacy provenance-less
  synthetic precomputed inputs under their permitted rules. Pin partial-pair, missing-provenance,
  backend-conflict, grid-conflict, and mixed BZ/path errors.
- [ ] Re-run end-to-end absent-backend and explicit-brute legacy regression. Compare new provenance
  separately from all pre-existing result fields.
- [ ] Assert `S.info.dipole` and `S.path_dipole` agree for q-path runs, and that `S.info.grid`
  reports the exact selected BZ policy. Do not inspect Jensen-specific result fields and do not
  claim the masking symptom is resolved.
- [ ] Run all new tests individually and the full projected suite; report pass/fail/incomplete and
  commit.

**Commit:** `test(invz): close Gate-C7 cache provenance integration matrix`

---

## Task 10 — Whole-branch closure

- [ ] Run the Step-4 Ewald tests, every new/modified Step-5 test file independently, and
  `runtests('invz_projected/tests')`. Record pass/fail/incomplete counts and any assumption-filtered
  tests.
- [ ] Review `git diff <recorded-step5-base>..HEAD --stat` and the complete diff. Confirm only scoped
  production/docs/tests changed, the frozen preregistration is untouched, and no default flipped.
- [ ] Confirm absent/explicit brute are bit-identical to the runtime legacy reference; Ewald cold
  and warm are identical; cache contexts are exact; all five forwarding surfaces are covered; and
  full-Cartesian C4 plus C6/C7 pass.
- [ ] Compare the protected pre-task `invz_run_spectra.m` diff with the final worktree/index and
  confirm every user edit remains intact and unstaged unless the user explicitly chose otherwise.
- [ ] Perform a final whole-range review against the frozen preregistration, design, and integration
  map. Fix all correctness/contract findings before declaring Step 5 complete.
- [ ] Hand off Step 6 (Gate D, both Γ policies). Default remains `bruteforce`; Gate-E anchors and
  symptom/adoption claims remain prohibited until Step 7.

## Acceptance summary

Step 5 is complete only when:

1. The opt-in Ewald backend reaches jq modes, BZ couplings, q paths, both spectra drivers, and the
   interactive run script.
2. Legacy absent/explicit-brute behavior is runtime-reference bit-identical for all pre-existing
   outputs/fields.
3. Cache schema v5 is backend-separated and validates exact physics, backend, grid, and schema
   payloads.
4. Full-Cartesian Gate-C4, all three demag cases in Gate-C6, and the complete Gate-C7 matrix pass.
5. The helper for Step-7 anchors exists and is parity-tested, but no Gate-D/E physics conclusion has
   been taken out of order.
