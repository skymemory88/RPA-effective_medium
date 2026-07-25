function tests = test_invz_ewald_integration
%TEST_INVZ_EWALD_INTEGRATION Step-5 Task 9: Gate-C7 end-to-end cache/provenance matrix
% (docs/superpowers/plans/2026-07-24-ewald-step5-integration.md Task 9; authority
% docs/invzp_ewald_prereg.md Sec.5 Gate-C check 7). TEST-ONLY: no production edit.
%
% This is a MECHANICAL integration test exercising the already-committed opt-in Ewald
% plumbing (Tasks 1-8) through the top-level spectra drivers invz_spectra_map/
% invz_spectra_qpath -- it proves cache/provenance wiring, not physics. Every call in this
% file uses ONE fixed, deliberately HIGH field (T = 0.31 K, B = 8 T, well above both the
% bare-MF Bc (~4.3 T at this T) and the 1/z theory's own Bc_1z (~3 T)), confirmed by direct
% probe (not just literature inference) to give a clean, stable, non-masked PM point
% (phase/phase_1z = 2) on every small grid/backend/policy combination used below. This
% file NEVER inspects S.phase_1z/S.m_1z/S.D_ord or any other Jensen-specific field, and
% makes NO claim about the (separately tracked, out-of-Step-5-scope) low-field ordered-leg
% masking symptom -- it only needs a field high enough to stay comfortably clear of it.
%
% Covers (checklist bullets from the plan's Task 9):
%   1. Field-map (invz_spectra_map) and q-path (invz_spectra_qpath) routes, explicit
%      'bruteforce' and explicit 'ewald', with both P_complete and P_drop on the Ewald
%      half-open BZ grid route.
%   2. Cold-then-warm cache calls per backend/route: whole-struct isequaln on every
%      numerical output AND metadata field, with cleanup restricted to the exact new
%      jq5_* files each test creates (never the cache directory, never unrelated caches).
%   3. Backend-prefix separation (jq5_bruteforce_ vs jq5_ewald_) plus a discriminator
%      matrix: q, lattice, basis, exchange J12, demag/aspect, an Ewald control, grid
%      convention, offset, Gamma policy, and schema each independently prevent a false
%      cache hit -- including one deliberately hash-colliding q-vector pair (same
%      jq5_bruteforce_* filename by construction) and one hand-corrupted on-disk schema
%      tag, both proven rejected by the exact cacheMeta isequaln check, never a stale
%      return.
%   4. Computed, complete-precomputed, and legacy provenance-less synthetic precomputed
%      couplings under their permitted rules; the partial-pair, missing-provenance,
%      backend-conflict, grid-conflict, and mixed BZ/path-dipole-mismatch error paths.
%   5. Absent-opts.dipole vs explicit 'bruteforce' end-to-end legacy regression (whole-
%      struct parity), with the additive provenance fields (info.dipole/info.grid/
%      S.path_dipole) separately, explicitly inspected -- not just swept into the blanket
%      equality.
%   6. isequaln(S.info.dipole, S.path_dipole) for q-path runs, and S.info.grid reporting
%      the exact requested BZ policy (convention/offset/gammaPolicy) on the new route.
%
% Every cache-writing test snapshots invz_projected/cache/jq5_*.mat before it runs and
% deletes only the file names that are new afterward (helper snapshot_cache/
% cleanup_new_cache_files below) -- the SAME "capture the dir listing before/after"
% discipline the task brief requires. Distinct small dpRng "cache tags" are used per
% bruteforce-backend test (matching the test_invz_jq_modes_ewald.m/
% test_invz_bz_couplings_ewald.m precedent) purely to keep this file's own cache
% footprint from colliding with itself across test functions; dpRng is NOT a free
% parameter of the physics (it is MF_dipole's real-space cutoff) and is unused/ignored
% by the Ewald backend.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                        % invz_projected
addpath(fullfile(here, '..', '..'));                  % repo root
addpath(fullfile(here, '..', '..', 'invz_common'));   % invz_ion, invz_const, getf

ion = invz_ion();
testCase.TestData.ion      = ion;
testCase.TestData.T        = 0.31;             % K
testCase.TestData.B        = 8;                % T -- one stable HIGH field, well above Bc/Bc_1z
testCase.TestData.w        = (0.05:0.05:0.2).';        % minimal frequency vector (meV)
testCase.TestData.qpath    = [1 0 0; 2 0 0];            % short q-path
testCase.TestData.grid     = [4 4 4];                   % small BZ grid
testCase.TestData.cacheDir = fullfile(here, '..', 'cache');
testCase.TestData.eoDefault = default_eopts(ion.a);
end

% =====================================================================
% shared helpers (names deliberately free of "test" so functiontests skips them)
% =====================================================================
function a0 = alpha0_of(a)
a0 = sqrt(pi)/abs(det(a))^(1/3);
end

function eo = mk_eopts(alpha, r_cut, g_cut, boundary)
eo = struct('alpha', alpha, 'r_cut', r_cut, 'g_cut', g_cut, 'boundary', boundary);
end

function eo = default_eopts(a)
a0 = alpha0_of(a);
eo = mk_eopts(a0, 5.5/a0, 11*a0, 'conducting_k0_omitted');   % frozen production defaults (prereg Sec.2)
end

function [Jnu, info] = legacy_precomputed_fixture()
% The pre-Task-7 synthetic fixture already used throughout test_invz_spectra_map.m /
% test_invz_spectra_qpath.m / test_invz_spectra_forward_ewald.m: provenance-less (no
% .dipole, no .grid) -- the backward-compatibility case that MUST keep working unchanged.
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
end

function names = snapshot_cache(cacheDir)
d = dir(fullfile(cacheDir, 'jq5_*.mat'));
names = sort({d.name});
end

function cleanup_new_cache_files(cacheDir, before)
% Deletes ONLY the jq5_*.mat files that appeared since `before` was captured -- never
% the cache directory itself, never a jq3_/jq4_/odd1_ file, never another test's
% pre-existing entry (task brief: "clean up ONLY the exact test-created jq5_* files").
after = snapshot_cache(cacheDir);
newOnes = setdiff(after, before);
for i = 1:numel(newOnes)
    delete(fullfile(cacheDir, newOnes{i}));
end
end

function tok = stash_matching_caches(cacheDir, pattern)
% Move any pre-existing cache files matching `pattern` aside so a test can create a known-fresh
% file at that exact name WITHOUT destroying unrelated (e.g. production) caches; restore with
% restore_stashed_caches(tok). Respects the brief: never permanently deletes an unrelated cache
% (review finding 2).
d = dir(fullfile(cacheDir, pattern));
tok = struct('cacheDir', cacheDir, 'orig', {{}}, 'bak', {{}});
for i = 1:numel(d)
    orig = fullfile(cacheDir, d(i).name);
    bak  = [orig '.mask9stash'];
    movefile(orig, bak);
    tok.orig{end+1} = orig; tok.bak{end+1} = bak;
end
end

function restore_stashed_caches(tok)
for i = 1:numel(tok.orig)
    if exist(tok.bak{i}, 'file'), movefile(tok.bak{i}, tok.orig{i}); end
end
end

function h = local_hash_vec_replica(v)
% Byte-for-byte replica of invz_jq_modes.m's PRIVATE hash_vec(v) local function, used
% ONLY to self-verify a deliberately constructed q-vector cache-filename collision
% below (test_discriminator_q_collision_rejected_by_exact_cachemeta) -- never used to
% compute an actual cache key; production's own hash_vec is exercised unmodified
% through the real invz_spectra_qpath call.
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end

% =====================================================================
% SECTION 1 -- routes x backends x Gamma policies, cold then warm (bullets 1-2)
% =====================================================================
function test_map_bruteforce_computed_cold_warm(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B;
w = testCase.TestData.w; grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
dpRng = 12;
opts = struct('grid', grid, 'dpRng', dpRng, 'cache', true, 'verbose', false, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
Scold = invz_spectra_map(ion, T, B, w, opts);
afterCold = cache_count(cacheDir, 'jq5_bruteforce_*.mat');
Swarm = invz_spectra_map(ion, T, B, w, opts);
afterWarm = cache_count(cacheDir, 'jq5_bruteforce_*.mat');

verifyTrue(testCase, isequaln(Scold, Swarm), 'bruteforce map: cold vs warm S struct differ.');
verifyEqual(testCase, afterWarm, afterCold, 'warm call must not write an additional cache file.');
verifyEqual(testCase, Scold.info.dipole.backend, 'bruteforce');
verifyFalse(testCase, isfield(Scold.info, 'grid'));

cleanup_new_cache_files(cacheDir, before);
end

function n = cache_count(cacheDir, pattern)
d = dir(fullfile(cacheDir, pattern));
n = numel(d);
end

function test_map_ewald_pcomplete_cold_warm(testCase)
ewald_pgamma_cold_warm_map(testCase, 'P_complete');
end

function test_map_ewald_pdrop_cold_warm(testCase)
ewald_pgamma_cold_warm_map(testCase, 'P_drop');
end

function ewald_pgamma_cold_warm_map(testCase, gammaPolicy)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B;
w = testCase.TestData.w; grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
eo = testCase.TestData.eoDefault;
opts = struct('grid', grid, 'dpRng', 6, 'cache', true, 'verbose', false, ...
    'dipole', 'ewald', 'ewald', eo, 'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', gammaPolicy);

before = snapshot_cache(cacheDir);
Scold = invz_spectra_map(ion, T, B, w, opts);
Swarm = invz_spectra_map(ion, T, B, w, opts);

verifyTrue(testCase, isequaln(Scold, Swarm), ...
    sprintf('ewald(%s) map: cold vs warm S struct differ.', gammaPolicy));
verifyEqual(testCase, Scold.info.dipole.backend, 'ewald');
verifyTrue(testCase, isnan(Scold.info.dpRng));
verifyEqual(testCase, Scold.info.grid.convention, 'halfopen');
verifyEqual(testCase, Scold.info.grid.gammaPolicy, gammaPolicy);

cleanup_new_cache_files(cacheDir, before);
end

function test_qpath_bruteforce_computed_cold_warm(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B;
w = testCase.TestData.w; qpath = testCase.TestData.qpath; grid = testCase.TestData.grid;
cacheDir = testCase.TestData.cacheDir;
opts = struct('grid', grid, 'dpRng', 14, 'cache', true, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
Scold = invz_spectra_qpath(ion, T, B, qpath, w, opts);
Swarm = invz_spectra_qpath(ion, T, B, qpath, w, opts);

verifyTrue(testCase, isequaln(Scold, Swarm), 'bruteforce qpath: cold vs warm S struct differ.');
verifyEqual(testCase, Scold.path_dipole.backend, 'bruteforce');

cleanup_new_cache_files(cacheDir, before);
end

function test_qpath_ewald_pcomplete_cold_warm(testCase)
ewald_pgamma_cold_warm_qpath(testCase, 'P_complete');
end

function test_qpath_ewald_pdrop_cold_warm(testCase)
ewald_pgamma_cold_warm_qpath(testCase, 'P_drop');
end

function ewald_pgamma_cold_warm_qpath(testCase, gammaPolicy)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B;
w = testCase.TestData.w; qpath = testCase.TestData.qpath; grid = testCase.TestData.grid;
cacheDir = testCase.TestData.cacheDir; eo = testCase.TestData.eoDefault;
opts = struct('grid', grid, 'dpRng', 6, 'cache', true, 'dipole', 'ewald', 'ewald', eo, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', gammaPolicy);

before = snapshot_cache(cacheDir);
Scold = invz_spectra_qpath(ion, T, B, qpath, w, opts);
Swarm = invz_spectra_qpath(ion, T, B, qpath, w, opts);

verifyTrue(testCase, isequaln(Scold, Swarm), ...
    sprintf('ewald(%s) qpath: cold vs warm S struct differ.', gammaPolicy));
verifyEqual(testCase, Scold.path_dipole.backend, 'ewald');
verifyEqual(testCase, Scold.info.grid.gammaPolicy, gammaPolicy);
verifyTrue(testCase, isequaln(Scold.path_dipole, Scold.info.dipole));

cleanup_new_cache_files(cacheDir, before);
end

function test_backend_cache_filenames_have_distinct_prefixes(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B;
w = testCase.TestData.w; grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
eo = testCase.TestData.eoDefault;

before = snapshot_cache(cacheDir);
invz_spectra_map(ion, T, B, w, struct('grid', grid, 'dpRng', 15, 'cache', true, ...
    'verbose', false, 'dipole', 'bruteforce'));
invz_spectra_map(ion, T, B, w, struct('grid', grid, 'dpRng', 6, 'cache', true, ...
    'verbose', false, 'dipole', 'ewald', 'ewald', eo));
after = snapshot_cache(cacheDir);
newFiles = setdiff(after, before);

isBrute = ~cellfun('isempty', regexp(newFiles, '^jq5_bruteforce_', 'once'));
isEwald = ~cellfun('isempty', regexp(newFiles, '^jq5_ewald_', 'once'));
verifyGreaterThanOrEqual(testCase, sum(isBrute), 1, 'expected at least one new jq5_bruteforce_* file.');
verifyGreaterThanOrEqual(testCase, sum(isEwald), 1, 'expected at least one new jq5_ewald_* file.');
verifyTrue(testCase, all(isBrute | isEwald), ...
    'every new cache file must carry the literal jq5_bruteforce_ or jq5_ewald_ backend prefix.');

cleanup_new_cache_files(cacheDir, before);
end

function test_warm_call_returns_ondisk_payload_proving_a_read(testCase)
% Review finding 1: every cold->warm test above asserts only isequaln(Scold,Swarm) and an
% unchanged cache-file count -- BOTH pass identically whether the warm call genuinely read
% the cache OR silently recomputed (a cache=true recompute overwrites the same deterministic
% filename, so the count is unchanged either way). This is the INVERSE of
% test_discriminator_schema_corruption_rejected below: it corrupts the on-disk PAYLOAD while
% leaving cacheMeta completely untouched, so the warm call's own cache-hit check ACCEPTS the
% entry -- the only way S2.info.Jcc0 can come back as the poisoned sentinel is if the warm
% call actually read it off disk (a silent recompute would return the true, freshly computed
% Jcc0 instead). info.Jcc0 flows unmodified from invz_jq_modes' cached `info` all the way
% through invz_bz_couplings/invz_spectra_map to S.info.Jcc0 (invz_jq_modes.m: on a cache hit
% `info = S.info` verbatim; invz_spectra_map.m: `Jcc0 = info.Jcc0; ... S.info = info;`) -- the
% same passthrough test_discriminator_schema_corruption_rejected already relies on. The
% pre-existing slate at this dpRng is STASHED (never permanently deleted -- finding 2) and
% restored by onCleanup even if an assertion below throws.
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
dpRng = 27;                       % unused by any sibling test/pre-existing cache in this file
% A same-order-of-magnitude sentinel (true Jcc0 here is ~6.4e-3 meV, matching the codebase's
% own legacy_precomputed_fixture) -- NOT the schema test's -999/-777-style astronomically
% off-scale sentinel: this call goes through the FULL physics solve (the cache hit is
% genuinely ACCEPTED, unlike the schema test where it is rejected before ever reaching the
% solver), and an off-scale J0eff sends invz_single_ion's mean-field iteration into
% non-convergence (a printed warning, breaking the pristine-output requirement) without
% making the read-vs-recompute proof any stronger -- any value != trueJcc0 proves the read.
SENTINEL = -0.00777;
pattern = sprintf('jq5_bruteforce_%d_*.mat', dpRng);
stashTok = stash_matching_caches(cacheDir, pattern);
cleaner  = onCleanup(@() restore_stashed_caches(stashTok));   % restores even if an assertion throws

opts = struct('grid', grid, 'dpRng', dpRng, 'cache', true, 'verbose', false, 'dipole', 'bruteforce');

% CONTROL -- independent, cache-bypassing ground truth for this exact configuration.
Sfresh = invz_spectra_map(ion, T, B, w, setfield(opts, 'cache', false)); %#ok<SFLD>
trueJcc0 = Sfresh.info.Jcc0;
verifyTrue(testCase, isfinite(trueJcc0), 'control: fresh recompute must give a finite Jcc0.');

before = snapshot_cache(cacheDir);
S1 = invz_spectra_map(ion, T, B, w, opts);          % cold call: populates the cache

d = dir(fullfile(cacheDir, pattern));
verifyEqual(testCase, numel(d), 1, 'expected exactly one new cache file from the cold map call.');
fpath = fullfile(cacheDir, d(1).name);

Sload = load(fpath);
verifyTrue(testCase, isfield(Sload, 'cacheMeta') && isfield(Sload.cacheMeta, 'schema'));
verifyEqual(testCase, Sload.cacheMeta.schema, 'invz_jq_modes/v5');
Sload.info.Jcc0 = SENTINEL;               % poison ONLY the payload
Sload.Jnu(:) = SENTINEL;
save(fpath, '-struct', 'Sload');          % cacheMeta left completely untouched -> hit accepted

S2 = invz_spectra_map(ion, T, B, w, opts);   % warm call, identical opts

verifyEqual(testCase, S2.info.Jcc0, SENTINEL, ...
    ['genuine-read proof: the warm call did not return the poisoned on-disk payload -- ' ...
     'it recomputed instead of reading the cache.']);
verifyNotEqual(testCase, S2.info.Jcc0, trueJcc0);

cleanup_new_cache_files(cacheDir, before);
end

% =====================================================================
% SECTION 2 -- discriminator matrix: q, lattice, basis, exchange, demag/aspect, Ewald
% control, grid convention, offset, Gamma policy, schema each prevent a false cache hit
% (bullet 3)
% =====================================================================
function test_discriminator_q_collision_rejected_by_exact_cachemeta(testCase)
% Two DIFFERENT single-row q-paths whose hash_vec digest collides (identical weighted
% sum under weights 1:3), so BOTH resolve to the IDENTICAL jq5_bruteforce_* filename
% despite genuinely different q content -- a real, non-contrived "deliberately
% colliding compact filename digest" (task brief). Proves the frozen contract that
% filename-hash collision resistance is NOT the safety net; exact cacheMeta.qvec
% isequaln validation is.
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
cacheDir = testCase.TestData.cacheDir;
dpRng = 16;

q1 = [0.5 0.25 0.125];
q2 = [0.75 0.125 0.125];   % Dq_x=+0.25, Dq_y=-0.125: 1*Dq_x + 2*Dq_y = 0 -> identical weighted sum
verifyEqual(testCase, local_hash_vec_replica(q1(:)), local_hash_vec_replica(q2(:)), ...
    'test construction error: q1/q2 must collide under hash_vec for this test to be meaningful.');
verifyNotEqual(testCase, q1, q2);

opts = struct('grid', testCase.TestData.grid, 'dpRng', dpRng, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
S1 = invz_spectra_qpath(ion, T, B, q1, w, opts);
S2 = invz_spectra_qpath(ion, T, B, q2, w, opts);   % must NOT silently reuse q1's cached file/answer

Pref2 = invz_jq_path(ion, q2, struct('dpRng', dpRng, 'cache', false, 'dipole', 'bruteforce'));
verifyEqual(testCase, S2.Jq(1), Pref2.Juni(1), 'AbsTol', 0, ...
    ['q discriminator: the colliding-filename q2 call did not reproduce an independent fresh ' ...
     '(cache=false) q2 computation -- possible false cache hit on q1''s stale entry.']);
verifyNotEqual(testCase, S2.Jq(1), S1.Jq(1));

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_lattice_prevents_false_hit(testCase)
T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
ionBase = testCase.TestData.ion;
ionMod = ionBase;  ionMod.a = ionBase.a * (1 + 3e-3);  ionMod.Vc = abs(det(ionMod.a));
dpRng = 18;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ionBase, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, true));
Smod_fresh  = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'lattice discriminator: cached modified-lattice call did not match an independent fresh recomputation.');
verifyGreaterThan(testCase, abs(Smod_cached.info.Jcc0 - Sbase.info.Jcc0), 1e-8, ...
    'lattice discriminator: modified-lattice Jcc0 is suspiciously identical to the base lattice''s.');

cleanup_new_cache_files(cacheDir, before);
end

function o2 = setcache(o, tf)
o2 = o;  o2.cache = tf;
end

function test_discriminator_basis_prevents_false_hit(testCase)
T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
ionBase = testCase.TestData.ion;
ionMod = ionBase;  ionMod.tau(2, 1) = ionMod.tau(2, 1) + 0.01;
dpRng = 20;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ionBase, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, true));
Smod_fresh  = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'basis discriminator: cached modified-basis call did not match an independent fresh recomputation.');
verifyGreaterThan(testCase, abs(Smod_cached.info.Jcc0 - Sbase.info.Jcc0), 1e-8, ...
    'basis discriminator: modified-basis Jcc0 is suspiciously identical to the base basis''s.');

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_exchange_J12_prevents_false_hit(testCase)
T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
ionBase = testCase.TestData.ion;
ionMod = ionBase;  ionMod.J12 = ionBase.J12 * 3;
dpRng = 21;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ionBase, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, true));
Smod_fresh  = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'J12 discriminator: cached modified-J12 call did not match an independent fresh recomputation.');
% Jcc0 = Jcc0_dipole + 4*J12, and Jcc0_dipole does not depend on J12 -- an exact relation.
verifyEqual(testCase, Smod_cached.info.Jcc0 - Sbase.info.Jcc0, 4*(ionMod.J12 - ionBase.J12), 'AbsTol', 1e-13);

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_demag_aspect_prevents_false_hit(testCase)
T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
ionBase = testCase.TestData.ion;   % demag = 0 (off)
ionMod = ionBase;  ionMod.demag = 1;  ionMod.alpha = 0.5;
dpRng = 22;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce');

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ionBase, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, true));
Smod_fresh  = invz_spectra_map(ionMod, T, B, w, setcache(optsBase, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'demag/aspect discriminator: cached modified-shape call did not match an independent fresh recomputation.');
% R2007 demag-invariance of the criticality coupling: Jcc0 unchanged...
verifyEqual(testCase, Smod_cached.info.Jcc0, Sbase.info.Jcc0, 'AbsTol', 1e-12);
% ...while Jshape_cc/Jaa0 DO change -- confirms demag was genuinely applied, not silently ignored.
verifyNotEqual(testCase, Smod_cached.info.Jshape_cc, Sbase.info.Jshape_cc);
verifyNotEqual(testCase, Smod_cached.info.Jaa0, Sbase.info.Jaa0);

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_ewald_control_prevents_false_hit(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
eoBase = testCase.TestData.eoDefault;
eoMod = eoBase;  eoMod.r_cut = eoBase.r_cut * 1.2;
optsBase = struct('grid', grid, 'dpRng', 6, 'verbose', false, 'dipole', 'ewald', 'ewald', eoBase);
optsMod  = struct('grid', grid, 'dpRng', 6, 'verbose', false, 'dipole', 'ewald', 'ewald', eoMod);

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ion, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ion, T, B, w, setcache(optsMod, true));
Smod_fresh  = invz_spectra_map(ion, T, B, w, setcache(optsMod, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'Ewald-control discriminator: cached modified-r_cut call did not match an independent fresh recomputation.');
verifyEqual(testCase, Smod_cached.info.dipole.ewald.r_cut, eoMod.r_cut);
verifyNotEqual(testCase, Smod_cached.info.dipole.ewald.r_cut, Sbase.info.dipole.ewald.r_cut);

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_grid_convention_prevents_false_hit(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
dpRng = 23;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce', ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_drop');
optsMod = optsBase;  optsMod.gridConvention = 'legacy_inclusive';

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ion, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ion, T, B, w, setcache(optsMod, true));
Smod_fresh  = invz_spectra_map(ion, T, B, w, setcache(optsMod, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'grid-convention discriminator: cached modified-convention call did not match an independent fresh recomputation.');
verifyEqual(testCase, Smod_cached.info.grid.convention, 'legacy_inclusive');
verifyNotEqual(testCase, Smod_cached.info.grid.convention, Sbase.info.grid.convention);

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_grid_offset_prevents_false_hit(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
dpRng = 24;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce', ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_drop');
optsMod = optsBase;  optsMod.gridOffset = [1 0 0];

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ion, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ion, T, B, w, setcache(optsMod, true));
Smod_fresh  = invz_spectra_map(ion, T, B, w, setcache(optsMod, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'grid-offset discriminator: cached modified-offset call did not match an independent fresh recomputation.');
verifyEqual(testCase, Smod_cached.info.grid.offset, logical([1 0 0]));
verifyNotEqual(testCase, Smod_cached.info.grid.offset, Sbase.info.grid.offset);

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_gamma_policy_prevents_false_hit(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
dpRng = 25;
optsBase = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce', ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_complete');
optsMod = optsBase;  optsMod.gammaPolicy = 'P_drop';

before = snapshot_cache(cacheDir);
Sbase = invz_spectra_map(ion, T, B, w, setcache(optsBase, true));
Smod_cached = invz_spectra_map(ion, T, B, w, setcache(optsMod, true));
Smod_fresh  = invz_spectra_map(ion, T, B, w, setcache(optsMod, false));

verifyTrue(testCase, isequaln(Smod_cached.info, Smod_fresh.info), ...
    'Gamma-policy discriminator: cached modified-policy call did not match an independent fresh recomputation.');
verifyEqual(testCase, Smod_cached.info.grid.gammaPolicy, 'P_drop');
verifyEqual(testCase, Sbase.info.grid.retained, Sbase.info.grid.nominal);                    % P_complete
verifyEqual(testCase, Smod_cached.info.grid.retained, Smod_cached.info.grid.nominal - Smod_cached.info.grid.n_gamma); % P_drop

cleanup_new_cache_files(cacheDir, before);
end

function test_discriminator_schema_corruption_rejected(testCase)
% Hand-corrupt an on-disk jq5_bruteforce_* payload's cacheMeta.schema tag (simulating a
% stale/legacy schema version) plus an obviously-wrong Jnu/info sentinel, at the EXACT
% filename the next identical call will look up -- confirms the exact-cacheMeta check
% rejects it (recompute, not a stale/sentinel return), mirroring the
% test_invz_jq_modes_ewald.m v5-cache precedent but reached through the spectra driver.
% The pre-existing slate at this dpRng is STASHED (never permanently deleted -- review
% finding 2: the cache dir is shared with production/sibling tests, e.g. a future v5
% production cache at dpRng ~= 26) and restored by onCleanup even if an assertion throws.
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid; cacheDir = testCase.TestData.cacheDir;
dpRng = 26;
pattern = sprintf('jq5_bruteforce_%d_*.mat', dpRng);
stashTok = stash_matching_caches(cacheDir, pattern);
cleaner  = onCleanup(@() restore_stashed_caches(stashTok));   % restores even if an assertion throws

opts = struct('grid', grid, 'dpRng', dpRng, 'verbose', false, 'dipole', 'bruteforce', 'cache', true);
before = snapshot_cache(cacheDir);
S1 = invz_spectra_map(ion, T, B, w, opts);

d = dir(fullfile(cacheDir, pattern));
verifyEqual(testCase, numel(d), 1, 'expected exactly one new cache file from the cold map call.');
fpath = fullfile(cacheDir, d(1).name);

Sload = load(fpath);
verifyTrue(testCase, isfield(Sload, 'cacheMeta') && isfield(Sload.cacheMeta, 'schema'));
verifyEqual(testCase, Sload.cacheMeta.schema, 'invz_jq_modes/v5');
Sload.cacheMeta.schema = 'invz_jq_modes/v4';    % corrupt: simulate a stale schema tag
Sload.info.Jcc0 = -999;                          % obviously-wrong sentinel
Sload.Jnu = zeros(size(Sload.Jnu)) - 999;
save(fpath, '-struct', 'Sload');

S2 = invz_spectra_map(ion, T, B, w, opts);   % same opts: must detect the schema mismatch, recompute
verifyTrue(testCase, isequaln(S2.info, S1.info), ...
    'schema discriminator: a corrupted on-disk schema tag was NOT rejected -- returned a stale/sentinel value.');
verifyNotEqual(testCase, S2.info.Jcc0, -999);

cleanup_new_cache_files(cacheDir, before);
end

% =====================================================================
% SECTION 3 -- computed / complete-precomputed / legacy provenance-less couplings +
% pinned error paths (bullet 4)
% =====================================================================
function test_map_precomputed_complete_matches_computed(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid;
[Jnu, info] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', 6, 'cache', false, 'dipole', 'bruteforce'));
Sprecomp = invz_spectra_map(ion, T, B, w, struct('Jnu', Jnu, 'info', info, 'verbose', false));
Scomputed = invz_spectra_map(ion, T, B, w, struct('grid', grid, 'dpRng', 6, 'cache', false, ...
    'verbose', false, 'dipole', 'bruteforce'));

verifyEqual(testCase, Sprecomp.info, info);
verifyTrue(testCase, isequaln(Sprecomp.chiz, Scomputed.chiz));
verifyTrue(testCase, isequaln(Sprecomp.chirpa, Scomputed.chirpa));
end

function test_map_legacy_provenanceless_precomputed_accepted(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
[Jnu, info] = legacy_precomputed_fixture();
S = invz_spectra_map(ion, T, B, w, struct('Jnu', Jnu, 'info', info, 'verbose', false));
verifyEqual(testCase, S.info, info);
verifySize(testCase, S.chiz, [numel(w) 1]);
end

function test_qpath_legacy_provenanceless_precomputed_accepted(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath;
[Jnu, info] = legacy_precomputed_fixture();
S = invz_spectra_qpath(ion, T, B, qpath, w, struct('Jnu', Jnu, 'info', info));
verifyEqual(testCase, S.info, info);
verifyEqual(testCase, S.path_dipole.backend, 'bruteforce');
verifySize(testCase, S.chiz, [numel(w) size(qpath, 1)]);
end

function test_map_partial_precomputed_pair_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
[Jnu, info] = legacy_precomputed_fixture();
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, struct('Jnu', Jnu, 'verbose', false)), ...
    'invz:spectraPrecomputedPartial');
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, struct('info', info, 'verbose', false)), ...
    'invz:spectraPrecomputedPartial');
end

function test_qpath_partial_precomputed_pair_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath;
[Jnu, info] = legacy_precomputed_fixture();
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, struct('Jnu', Jnu)), ...
    'invz:spectraPrecomputedPartial');
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, struct('info', info)), ...
    'invz:spectraPrecomputedPartial');
end

function test_map_explicit_request_missing_provenance_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
[Jnu, info] = legacy_precomputed_fixture();
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'dipole', 'bruteforce')), ...
    'invz:spectraBackendProvenanceMissing');
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'gammaPolicy', 'P_complete')), ...
    'invz:spectraGridProvenanceMissing');
end

function test_qpath_explicit_request_missing_provenance_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath;
[Jnu, info] = legacy_precomputed_fixture();
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'dipole', 'bruteforce')), ...
    'invz:spectraBackendProvenanceMissing');
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'gridOffset', [0 0 0])), ...
    'invz:spectraGridProvenanceMissing');
end

function test_map_conflicting_backend_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
eo = testCase.TestData.eoDefault;
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false));   % bruteforce
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'dipole', 'ewald', 'ewald', eo)), ...
    'invz:spectraBackendConflict');

[JnuE, infoE] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
eo2 = eo;  eo2.r_cut = eo.r_cut * 1.1;
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, ...
    struct('Jnu', JnuE, 'info', infoE, 'verbose', false, 'dipole', 'ewald', 'ewald', eo2)), ...
    'invz:spectraBackendConflict');
end

function test_qpath_conflicting_backend_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath; eo = testCase.TestData.eoDefault;
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false));   % bruteforce
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'dipole', 'ewald', 'ewald', eo)), ...
    'invz:spectraBackendConflict');
end

function test_map_conflicting_grid_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen'));
verifyError(testCase, @() invz_spectra_map(ion, T, B, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'grid', [4 4 4], ...
           'gridConvention', 'legacy_inclusive')), ...
    'invz:spectraGridConflict');
end

function test_qpath_conflicting_grid_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath;
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen'));
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'grid', [6 6 6], 'gridConvention', 'halfopen')), ...
    'invz:spectraGridConflict');
end

function test_qpath_mixed_bz_path_backend_mismatch_errors(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
info = struct('Jcc0', 6.4e-3, 'dipole', struct('backend', 'bruteforce', ...
    'ewald', struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', ''), ...
    'q_reduction', 'STALE-CORRUPTED-PROVENANCE-NOT-THE-REAL-STRING', ...
    'primitive_schema', 'MF_dipole+exchange (legacy, unversioned)'));
verifyError(testCase, @() invz_spectra_qpath(ion, T, B, qpath, w, struct('Jnu', Jnu, 'info', info)), ...
    'invz:spectraPathDipoleMismatch');
end

function test_qpath_ewald_provenance_inherited(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath; eo = testCase.TestData.eoDefault;
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
S = invz_spectra_qpath(ion, T, B, qpath, w, struct('Jnu', Jnu, 'info', info));
verifyEqual(testCase, S.path_dipole.backend, 'ewald');
verifyTrue(testCase, isequaln(S.path_dipole, S.info.dipole));
end

% =====================================================================
% SECTION 4 -- absent-opts.dipole vs explicit 'bruteforce' end-to-end legacy regression;
% new provenance fields compared SEPARATELY from pre-existing result fields (bullet 5)
% =====================================================================
function test_map_absent_vs_explicit_bruteforce_legacy_regression(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid;
optsAbsent   = struct('grid', grid, 'dpRng', 6, 'cache', false, 'verbose', false);
optsExplicit = struct('grid', grid, 'dpRng', 6, 'cache', false, 'verbose', false, 'dipole', 'bruteforce');

Sabsent   = invz_spectra_map(ion, T, B, w, optsAbsent);
Sexplicit = invz_spectra_map(ion, T, B, w, optsExplicit);

% Pre-existing legacy result fields (everything except the additive S.info.dipole/
% S.info.grid provenance carried inside S.info) must be bit-identical between the
% absent-backend and explicit-bruteforce calls.
legacyFieldsOfS = setdiff(fieldnames(Sabsent), {'info'});
for i = 1:numel(legacyFieldsOfS)
    f = legacyFieldsOfS{i};
    verifyTrue(testCase, isequaln(Sabsent.(f), Sexplicit.(f)), ...
        sprintf('map legacy regression: S.%s differs between absent-backend and explicit-bruteforce.', f));
end
verifyFalse(testCase, isfield(Sabsent.info, 'grid'));
verifyFalse(testCase, isfield(Sexplicit.info, 'grid'));

% NEW provenance field, compared SEPARATELY (both content and cross-call equality):
verifyEqual(testCase, Sabsent.info.dipole.backend, 'bruteforce');
verifyEqual(testCase, Sexplicit.info.dipole.backend, 'bruteforce');
verifyTrue(testCase, isequaln(Sabsent.info.dipole, Sexplicit.info.dipole));

% Whole-struct confirmation (both legacy and additive fields together): with backend
% resolving identically either way, the two calls should in fact be fully isequaln.
verifyTrue(testCase, isequaln(Sabsent, Sexplicit));
end

function test_qpath_absent_vs_explicit_bruteforce_legacy_regression(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath; grid = testCase.TestData.grid;
optsAbsent   = struct('grid', grid, 'dpRng', 6, 'cache', false);
optsExplicit = struct('grid', grid, 'dpRng', 6, 'cache', false, 'dipole', 'bruteforce');

Sabsent   = invz_spectra_qpath(ion, T, B, qpath, w, optsAbsent);
Sexplicit = invz_spectra_qpath(ion, T, B, qpath, w, optsExplicit);

legacyFieldsOfS = setdiff(fieldnames(Sabsent), {'info', 'path_dipole'});
for i = 1:numel(legacyFieldsOfS)
    f = legacyFieldsOfS{i};
    verifyTrue(testCase, isequaln(Sabsent.(f), Sexplicit.(f)), ...
        sprintf('qpath legacy regression: S.%s differs between absent-backend and explicit-bruteforce.', f));
end
verifyFalse(testCase, isfield(Sabsent.info, 'grid'));

% NEW provenance fields, compared SEPARATELY:
verifyEqual(testCase, Sabsent.path_dipole.backend, 'bruteforce');
verifyEqual(testCase, Sexplicit.path_dipole.backend, 'bruteforce');
verifyTrue(testCase, isequaln(Sabsent.path_dipole, Sexplicit.path_dipole));
verifyTrue(testCase, isequaln(Sabsent.info.dipole, Sexplicit.info.dipole));

verifyTrue(testCase, isequaln(Sabsent, Sexplicit));
end

% =====================================================================
% SECTION 5 -- S.info.dipole == S.path_dipole for q-path runs; S.info.grid reports the
% exact selected BZ policy (bullet 6)
% =====================================================================
function test_qpath_path_dipole_matches_info_dipole_bruteforce(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath; grid = testCase.TestData.grid;
S = invz_spectra_qpath(ion, T, B, qpath, w, struct('grid', grid, 'dpRng', 6, 'cache', false, ...
    'dipole', 'bruteforce'));
verifyTrue(testCase, isequaln(S.info.dipole, S.path_dipole));
verifyEqual(testCase, S.path_dipole.backend, 'bruteforce');
end

function test_qpath_path_dipole_matches_info_dipole_ewald(testCase)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
qpath = testCase.TestData.qpath; grid = testCase.TestData.grid; eo = testCase.TestData.eoDefault;
S = invz_spectra_qpath(ion, T, B, qpath, w, struct('grid', grid, 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
verifyTrue(testCase, isequaln(S.info.dipole, S.path_dipole));
verifyEqual(testCase, S.path_dipole.backend, 'ewald');
end

function test_map_info_grid_reports_exact_policy_pcomplete(testCase)
info_grid_reports_exact_policy(testCase, 'P_complete');
end

function test_map_info_grid_reports_exact_policy_pdrop(testCase)
info_grid_reports_exact_policy(testCase, 'P_drop');
end

function info_grid_reports_exact_policy(testCase, gammaPolicy)
ion = testCase.TestData.ion; T = testCase.TestData.T; B = testCase.TestData.B; w = testCase.TestData.w;
grid = testCase.TestData.grid;
S = invz_spectra_map(ion, T, B, w, struct('grid', grid, 'dpRng', 6, 'cache', false, ...
    'verbose', false, 'dipole', 'bruteforce', 'gridConvention', 'halfopen', ...
    'gridOffset', [0 0 0], 'gammaPolicy', gammaPolicy));
verifyEqual(testCase, S.info.grid.convention, 'halfopen');
verifyEqual(testCase, S.info.grid.offset, logical([0 0 0]));
verifyEqual(testCase, S.info.grid.gammaPolicy, gammaPolicy);
verifyEqual(testCase, S.info.grid.requested, grid);
if strcmp(gammaPolicy, 'P_complete')
    verifyEqual(testCase, S.info.grid.retained, S.info.grid.nominal);
else
    verifyEqual(testCase, S.info.grid.retained, S.info.grid.nominal - S.info.grid.n_gamma);
end
end
