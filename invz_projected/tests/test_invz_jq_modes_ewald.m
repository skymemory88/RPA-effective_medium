function tests = test_invz_jq_modes_ewald
%TEST_INVZ_JQ_MODES_EWALD Step-5 Task 2: invz_jq_modes backend dispatch, Gamma
% metadata (info.dipole / info.Jpath_base_cc / info.Jgamma_cc), and the v5
% cache contract (docs/superpowers/plans/2026-07-24-ewald-step5-integration.md
% Task 2; authority docs/invzp_ewald_prereg.md Sec.5 Gate-C check 7,
% docs/invzp_ewald_design.md Sec.2.3/4.2, docs/invzp_ewald_integration_map.md
% Sec.6.3).
%
% Legacy-field parity methodology follows the PLAN's "Frozen behavior
% contracts / Legacy regression" contract literally: "Every legacy info field
% exists and is individually isequaln to the reference. New fields are
% removed before comparing the legacy field set." This is a STRIP-then-compare
% (subset) check, deliberately different from test_invz_ewald_baselines.m's
% raw whole-struct field-SET equality helper (that frozen file is untouched
% here -- see this task's report for the interaction the new additive fields
% have with that helper).
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                        % invz_projected
addpath(fullfile(here, '..', '..'));                  % repo root: invz_dipole_ewald, MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));   % invz_ion, invz_const
testCase.TestData.cacheDir = fullfile(here, '..', 'cache');
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

function fn = legacy_info_fields()
% The 8 pre-existing invz_jq_modes info fields (pre-Step-5; matches
% invz_legacy_coupling_reference.m's local_jq_modes_ref output exactly).
c = {'Jcc0_dipole','Jaa0_dipole','Jcc0','Jaa0','Jshape_cc','dpRng','geomD','geomX'};
fn = sort(c(:));
end

function s2 = strip_fields(s, names)
s2 = s;
for i = 1:numel(names)
    if isfield(s2, names{i}), s2 = rmfield(s2, names{i}); end
end
end

function n = cache_file_count(cacheDir, pattern)
d = dir(fullfile(cacheDir, pattern));
n = numel(d);
end

function local_clear_matching(cacheDir, pattern)
% Explicit, non-wildcard-ambiguous delete of every file matching pattern
% (cache files are disposable/gitignored; used only to guarantee a clean
% slate for tests that need to observe a genuinely fresh cold write).
d = dir(fullfile(cacheDir, pattern));
for i = 1:numel(d)
    delete(fullfile(cacheDir, d(i).name));
end
end

% =====================================================================
% Validation: backend / opts.ewald cross-field / boundary / ODD
% (stable invz:jqModes* identifiers)
% =====================================================================
function test_unknown_backend_errors(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole','bogus','cache',false)), ...
    'invz:jqModesBackend');
end

function test_nonscalar_backend_errors(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole',{{'bruteforce'}},'cache',false)), ...
    'invz:jqModesBackend');
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole',char('bruteforce','ewald'),'cache',false)), ...
    'invz:jqModesBackend');
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole',1,'cache',false)), ...
    'invz:jqModesBackend');
end

function test_ewald_opts_without_ewald_backend_errors_absent_dipole(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('ewald',eo,'cache',false)), ...
    'invz:jqModesEwaldOptsUnexpected');
end

function test_bruteforce_with_ewald_opts_errors(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole','bruteforce','ewald',eo,'cache',false)), ...
    'invz:jqModesEwaldOptsUnexpected');
end

function test_ewald_missing_opts_entirely_errors(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','cache',false)), ...
    'invz:jqModesEwaldOptsFields');
end

function test_ewald_missing_one_control_errors(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
eo = rmfield(eo, 'g_cut');
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','ewald',eo,'cache',false)), ...
    'invz:jqModesEwaldOptsFields');
end

function test_ewald_extra_control_errors(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
eo.extra_field = 1;
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','ewald',eo,'cache',false)), ...
    'invz:jqModesEwaldOptsFields');
end

function test_invalid_boundary_errors(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
eo.boundary = 'open_surface';
verifyError(testCase, @() invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','ewald',eo,'cache',false)), ...
    'invz:jqModesEwaldBoundary');
end

function test_active_odd_plus_ewald_errors(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
odd = struct('Xp', eye(2));
verifyError(testCase, ...
    @() invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','ewald',eo,'odd',odd,'cache',false)), ...
    'invz:jqModesOddEwald');
end

function test_odd_false_plus_ewald_does_not_error(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[Jnu, info] = invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','ewald',eo,'odd',false,'cache',false));
verifyEqual(testCase, info.dipole.backend, 'ewald');
verifySize(testCase, Jnu, [1 4]);
end

function test_odd_brute_force_still_unaffected(testCase)
% ODD + (absent/explicit) bruteforce must remain completely unchanged: no
% error, dispatch reaches jq_modes_odd exactly as before Task 2.
ion = invz_ion();
odd = struct('Xp', zeros(2));
[Jnu, info] = invz_jq_modes(ion, [0.1 0 0], struct('odd', odd, 'dpRng', 10, 'cache', false)); %#ok<ASGLU>
verifyTrue(testCase, isfield(info, 'odd'));   % jq_modes_odd's own info contract, untouched
end

% =====================================================================
% Brute-force path: legacy-field parity (strip-then-compare) + additive
% Gamma metadata, constructed in the frozen operation order
% =====================================================================
function test_bruteforce_legacy_fields_match_reference_nocache(testCase)
ion = invz_ion();
qvec = [0.310 -0.180 0.240; 0 0 0; 2.000 0 0; 1.985 0.010 -0.005];
ref = invz_legacy_coupling_reference();
opts = struct('dpRng', 30, 'cache', false);

[Jnu_ref, info_ref, Juni_ref]    = ref.jq_modes(ion, qvec, opts);
[Jnu_prod, info_prod, Juni_prod] = invz_jq_modes(ion, qvec, opts);

verifyEqual(testCase, Jnu_prod, Jnu_ref);
verifyEqual(testCase, Juni_prod, Juni_ref);

legacyFn = legacy_info_fields();
info_prod_legacy = strip_fields(info_prod, {'dipole','Jpath_base_cc','Jgamma_cc'});
verifyEqual(testCase, sort(fieldnames(info_prod_legacy)), legacyFn);
verifyEqual(testCase, sort(fieldnames(info_ref)), legacyFn);
for i = 1:numel(legacyFn)
    f = legacyFn{i};
    verifyTrue(testCase, isequaln(info_prod_legacy.(f), info_ref.(f)), ...
        sprintf('info.%s differs from the frozen legacy reference.', f));
end
verifyTrue(testCase, isfield(info_prod, 'dipole'));
verifyTrue(testCase, isfield(info_prod, 'Jpath_base_cc'));
verifyTrue(testCase, isfield(info_prod, 'Jgamma_cc'));
end

function test_bruteforce_legacy_fields_match_reference_cached(testCase)
ion = invz_ion();
qvec = [0.310 -0.180 0.240; 0 0 0; 2.000 0 0];
ref = invz_legacy_coupling_reference();
opts = struct('dpRng', 30, 'cache', true);

[Jnu_cold, info_cold] = invz_jq_modes(ion, qvec, opts);
[Jnu_warm, info_warm] = invz_jq_modes(ion, qvec, opts);
verifyTrue(testCase, isequaln(Jnu_cold, Jnu_warm));
verifyTrue(testCase, isequaln(info_cold, info_warm));

[Jnu_ref, ~] = ref.jq_modes(ion, qvec, opts);
verifyEqual(testCase, Jnu_warm, Jnu_ref);

legacyFn = legacy_info_fields();
info_warm_legacy = strip_fields(info_warm, {'dipole','Jpath_base_cc','Jgamma_cc'});
verifyEqual(testCase, sort(fieldnames(info_warm_legacy)), legacyFn);
end

function test_absent_and_explicit_bruteforce_share_cache_identity(testCase)
ion = invz_ion();
qvec = [0.05 0.11 -0.02; 0.33 0 0];
cacheDir = testCase.TestData.cacheDir;
opts1 = struct('dpRng', 9, 'cache', true);                          % absent opts.dipole
[Jnu1, info1] = invz_jq_modes(ion, qvec, opts1);
afterFirst = cache_file_count(cacheDir, 'jq5_bruteforce_*.mat');

opts2 = struct('dpRng', 9, 'cache', true, 'dipole', 'bruteforce');  % explicit
[Jnu2, info2] = invz_jq_modes(ion, qvec, opts2);
afterSecond = cache_file_count(cacheDir, 'jq5_bruteforce_*.mat');

verifyEqual(testCase, afterSecond, afterFirst, ...
    ['explicit dipole=bruteforce with identical opts must hit the SAME cache file as ' ...
     'absent-backend (shared canonical identity), not write a new one.']);
verifyEqual(testCase, Jnu2, Jnu1);
verifyEqual(testCase, info2, info1);
end

function test_bruteforce_gamma_metadata_formula(testCase)
ion = invz_ion();
C = invz_const();
opts = struct('dpRng', 30, 'cache', false);
[~, info] = invz_jq_modes(ion, [0.1 0.2 -0.05], opts);

[dip0, ~, ~] = MF_dipole([0 0 0], 30, ion.a, ion.tau);
ex0 = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
lorz = 4*pi/(3*ion.Vc)*C.gfac;

expectedBase  = -C.gfac*squeeze(dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));
expectedGamma = expectedBase + lorz*ones(4);

verifyEqual(testCase, info.Jpath_base_cc, expectedBase, 'AbsTol', 1e-14);
verifyEqual(testCase, info.Jgamma_cc, expectedGamma, 'AbsTol', 1e-14);
verifyEqual(testCase, info.Jgamma_cc, info.Jpath_base_cc + lorz*ones(4), 'AbsTol', 0);
end

function test_info_dipole_shape_bruteforce(testCase)
ion = invz_ion();
[~, info] = invz_jq_modes(ion, [0 0 0], struct('cache', false));
verifyEqual(testCase, sort(fieldnames(info.dipole)), sort({'backend';'ewald';'q_reduction';'primitive_schema'}));
verifyEqual(testCase, info.dipole.backend, 'bruteforce');
verifyEqual(testCase, info.dipole.ewald, struct('alpha',[],'r_cut',[],'g_cut',[],'boundary',''));
verifyTrue(testCase, ischar(info.dipole.q_reduction) && ~isempty(info.dipole.q_reduction));
verifyTrue(testCase, ischar(info.dipole.primitive_schema));
end

% =====================================================================
% Ewald path: additive Gamma metadata, no extra Lorentz at Gamma,
% info.geomD never exposed, info.dpRng == NaN
% =====================================================================
function test_ewald_gamma_metadata_formula_and_zero_extra_lorentz(testCase)
ion = invz_ion();
C = invz_const();
eo = default_eopts(ion.a);
opts = struct('dipole', 'ewald', 'ewald', eo, 'cache', false);
[~, info] = invz_jq_modes(ion, [0.1 0.2 -0.05], opts);

dip0 = invz_dipole_ewald([0 0 0], ion.a, ion.tau, eo);
ex0  = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
lorz = 4*pi/(3*ion.Vc)*C.gfac;

Jex0cc = sign(ion.J12)*squeeze(ex0(3,3,:,:));
JgammaExpected = -C.gfac*squeeze(dip0(3,3,:,:)) + Jex0cc;      % "adds 0" extra Lorentz at Gamma
JbaseExpected  = JgammaExpected - lorz*ones(4);

verifyEqual(testCase, info.Jgamma_cc, JgammaExpected, 'AbsTol', 1e-12);
verifyEqual(testCase, info.Jpath_base_cc, JbaseExpected, 'AbsTol', 1e-12);
verifyEqual(testCase, info.Jgamma_cc, info.Jpath_base_cc + lorz*ones(4), 'AbsTol', 0);

v = ones(4,1)/2;
Jcc0d_expected = -squeeze(C.gfac*dip0(3,3,:,:));               % NO +lorz under Ewald
Jcc0d_expected = (Jcc0d_expected + Jcc0d_expected')/2;
verifyEqual(testCase, info.Jcc0_dipole, real(v.'*Jcc0d_expected*v), 'AbsTol', 1e-10);

Jaa0d_expected = -squeeze(C.gfac*dip0(1,1,:,:));                % demag=0 -> dm_aa=0, NO +lorz
Jaa0d_expected = (Jaa0d_expected + Jaa0d_expected')/2;
verifyEqual(testCase, info.Jaa0_dipole, real(v.'*Jaa0d_expected*v), 'AbsTol', 1e-10);
end

function test_ewald_geomD_absent_dprng_nan_geomX_present(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[~, info] = invz_jq_modes(ion, [0.02 -0.01 0.03], struct('dipole','ewald','ewald',eo,'cache',false));
verifyFalse(testCase, isfield(info, 'geomD'));
verifyTrue(testCase, isfield(info, 'geomX'));
verifyTrue(testCase, isnan(info.dpRng));
end

function test_info_dipole_shape_ewald(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[~, info] = invz_jq_modes(ion, [0 0 0], struct('dipole','ewald','ewald',eo,'cache',false));
verifyEqual(testCase, sort(fieldnames(info.dipole)), sort({'backend';'ewald';'q_reduction';'primitive_schema'}));
verifyEqual(testCase, info.dipole.backend, 'ewald');
verifyEqual(testCase, info.dipole.ewald, eo);
verifyTrue(testCase, ischar(info.dipole.q_reduction) && ~isempty(info.dipole.q_reduction));

[~, ~, geomCheck] = invz_dipole_ewald([0.01 0 0], ion.a, ion.tau, eo);
verifyEqual(testCase, info.dipole.q_reduction, geomCheck.fingerprint.qconv);
verifyEqual(testCase, info.dipole.primitive_schema, geomCheck.fingerprint.schema);
verifyEqual(testCase, info.dipole.primitive_schema, 'invz_dipole_ewald/v1');
end

function test_ewald_modes_real_and_finite(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09];
[Jnu, ~] = invz_jq_modes(ion, q, struct('dipole','ewald','ewald',eo,'cache',false));
verifySize(testCase, Jnu, [3 4]);
verifyTrue(testCase, all(isfinite(Jnu(:))));
verifyLessThan(testCase, max(abs(imag(Jnu(:)))), 1e-9);
end

function test_ewald_demag_shape_term_matches_bruteforce_semantics(testCase)
ion0 = invz_ion();
eo = default_eopts(ion0.a);
opts0 = struct('dipole', 'ewald', 'ewald', eo, 'cache', false);
qs = [0 0 0; 0.31 0.05 -0.12];
[J0nu, i0] = invz_jq_modes(ion0, qs, opts0);
verifyEqual(testCase, i0.Jshape_cc, 0);

C = invz_const();  lorz4 = 4*(4*pi/(3*ion0.Vc)*C.gfac);
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
[JSnu, iS] = invz_jq_modes(ionS, qs, opts0);
verifyEqual(testCase, JSnu, J0nu, 'AbsTol', 1e-12);          % branch spectrum demag-invariant (Ewald too)
verifyEqual(testCase, iS.Jcc0, i0.Jcc0, 'AbsTol', 1e-12);    % criticality coupling demag-invariant
verifyEqual(testCase, iS.Jshape_cc, lorz4, 'RelTol', 1e-12);
verifyEqual(testCase, iS.Jaa0_dipole, i0.Jaa0_dipole - lorz4, 'RelTol', 1e-9);
end

% =====================================================================
% v5 cache: schema, cold/warm, backend-separated filenames, safety
% (mismatched/legacy/malformed payloads are misses, not stale hits)
% =====================================================================
function test_v5_cache_schema_and_cold_warm_bruteforce(testCase)
ion = invz_ion();
qvec = [0.271 0.014 -0.333];
opts = struct('dpRng', 8, 'cache', true);
[Jnu_cold, info_cold, Juni_cold] = invz_jq_modes(ion, qvec, opts);
[Jnu_warm, info_warm, Juni_warm] = invz_jq_modes(ion, qvec, opts);
verifyTrue(testCase, isequaln(Jnu_cold, Jnu_warm));
verifyTrue(testCase, isequaln(info_cold, info_warm));
verifyTrue(testCase, isequaln(Juni_cold, Juni_warm));
verifyEqual(testCase, info_cold.dipole.backend, 'bruteforce');
end

function test_v5_cache_schema_and_cold_warm_ewald(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
qvec = [0.271 0.014 -0.333];
opts = struct('dipole', 'ewald', 'ewald', eo, 'cache', true);
[Jnu_cold, info_cold, Juni_cold] = invz_jq_modes(ion, qvec, opts);
[Jnu_warm, info_warm, Juni_warm] = invz_jq_modes(ion, qvec, opts);
verifyTrue(testCase, isequaln(Jnu_cold, Jnu_warm));
verifyTrue(testCase, isequaln(info_cold, info_warm));
verifyTrue(testCase, isequaln(Juni_cold, Juni_warm));
verifyEqual(testCase, info_cold.dipole.backend, 'ewald');
verifyTrue(testCase, isnan(info_cold.dpRng));
verifyFalse(testCase, isfield(info_cold, 'geomD'));
end

function test_v5_cache_backend_filenames_separated(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
qvec = [0.061 -0.129 0.404];
cacheDir = testCase.TestData.cacheDir;
invz_jq_modes(ion, qvec, struct('dpRng', 7, 'cache', true));
invz_jq_modes(ion, qvec, struct('dipole', 'ewald', 'ewald', eo, 'cache', true));
verifyGreaterThanOrEqual(testCase, cache_file_count(cacheDir, 'jq5_bruteforce_*.mat'), 1);
verifyGreaterThanOrEqual(testCase, cache_file_count(cacheDir, 'jq5_ewald_*.mat'), 1);
end

function test_v5_cache_rejects_mismatched_payload(testCase)
ion = invz_ion();
qvec = [0.09 0.02 -0.14];
cacheDir = testCase.TestData.cacheDir;
dpRng = 11;   % small (dpRng feeds MF_dipole's real-space meshgrid(-N:N,...) range --
              % NOT a free cache tag) but distinctive vs. other tests, so the pre-clear below is surgical
opts = struct('dpRng', dpRng, 'cache', true);
pattern = sprintf('jq5_bruteforce_%d_*.mat', dpRng);

% Guarantee a clean slate: this test needs a genuinely FRESH cold write (so it
% can then corrupt that specific file), independent of any leftover file a
% PRIOR, separate invocation of this same test may have left on disk (the
% cache directory is real, persistent disk state, not reset between runs).
local_clear_matching(cacheDir, pattern);

[Jnu1, ~] = invz_jq_modes(ion, qvec, opts);
d = dir(fullfile(cacheDir, pattern));
verifyEqual(testCase, numel(d), 1, 'expected exactly one new cache file from the cold call.');
fpath = fullfile(cacheDir, d(1).name);

S = load(fpath);
verifyTrue(testCase, isfield(S, 'cacheMeta'));
S.cacheMeta.Vc = S.cacheMeta.Vc + 1;      % corrupt the stored payload identity
S.Jnu = zeros(size(S.Jnu)) - 999;         % obviously-wrong sentinel data
save(fpath, '-struct', 'S');

[Jnu2, ~] = invz_jq_modes(ion, qvec, opts);   % same opts as the original cold call
verifyEqual(testCase, Jnu2, Jnu1, 'AbsTol', 0);   % must recompute, not trust the mismatched cache
verifyFalse(testCase, any(Jnu2(:) == -999));
end

function test_v5_cache_missing_cachemeta_is_miss(testCase)
ion = invz_ion();
qvec = [0.222 -0.011 0.093];
cacheDir = testCase.TestData.cacheDir;
dpRng = 13;   % small (see test_v5_cache_rejects_mismatched_payload) but distinctive vs. other tests
opts = struct('dpRng', dpRng, 'cache', true);
pattern = sprintf('jq5_bruteforce_%d_*.mat', dpRng);

local_clear_matching(cacheDir, pattern);   % see test_v5_cache_rejects_mismatched_payload

[Jnu1, ~] = invz_jq_modes(ion, qvec, opts);
d = dir(fullfile(cacheDir, pattern));
verifyEqual(testCase, numel(d), 1);
fpath = fullfile(cacheDir, d(1).name);

% Simulate a legacy/malformed payload (no cacheMeta field) at this v5 filename.
Jnu = zeros(size(Jnu1)) - 555; info = struct('bogus', 1); Juni = zeros(size(Jnu1,1),1); pkey = 0; %#ok<NASGU>
save(fpath, 'Jnu', 'info', 'Juni', 'pkey');

[Jnu2, info2] = invz_jq_modes(ion, qvec, opts);
verifyEqual(testCase, Jnu2, Jnu1, 'AbsTol', 0);
verifyTrue(testCase, isfield(info2, 'Jcc0'));
end
