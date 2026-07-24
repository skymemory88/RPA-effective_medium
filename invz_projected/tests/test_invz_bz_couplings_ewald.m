function tests = test_invz_bz_couplings_ewald
%TEST_INVZ_BZ_COUPLINGS_EWALD Step-5 Task 4: invz_bz_couplings grid/backend threading and
% exact grid provenance (docs/superpowers/plans/2026-07-24-ewald-step5-integration.md Task 4;
% authority docs/invzp_ewald_integration_map.md Sec.5.4, docs/invzp_phase1_quadrature_prereg.md).
%
% Covers: absent-grid-field legacy parity (info.grid absent, bit-identical to an inline
% pre-Task-4 reference reconstruction) and absent==explicit-bruteforce parity; presence
% routing for gridConvention/gridOffset/gammaPolicy (any one activates invz_phase1_qgrid;
% Ewald alone does not; explicit gridConvention='legacy_inclusive' activates the new route
% but is NOT claimed bit-identical); non-cubic grid rejection
% (invz:bzCouplingsAnisotropicHalfopen); complete info.grid field presence/types;
% P_complete vs P_drop row-count/provenance divergence; and cacheContext separation
% (legacy_bz vs phase1_qgrid vs a bare direct invz_jq_modes call's canonical direct_call
% sentinel) with an explicit no-false-cross-hit stress test.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                        % invz_projected: invz_bz_couplings, invz_jq_modes, invz_phase1_qgrid
addpath(fullfile(here, '..', '..'));                  % repo root: qVec_generator, invz_dipole_ewald, MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));   % invz_ion, invz_const, getf
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

function local_clear_matching(cacheDir, pattern)
% Explicit, non-wildcard-ambiguous delete of every file matching pattern (cache files are
% disposable/gitignored; used only to guarantee a clean slate for a test that needs to
% observe a genuinely fresh cold write).
d = dir(fullfile(cacheDir, pattern));
for i = 1:numel(d)
    delete(fullfile(cacheDir, d(i).name));
end
end

function qc = legacy_qc(ion, grid)
% Independent, byte-for-byte reproduction of invz_bz_couplings.m's pre-Task-4
% (and post-Task-4 absent-grid-route) qVec_generator + Gamma-drop construction --
% the "pre-edit baseline" the production edit must keep reproducing exactly.
[qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
qc = qc(any(abs(qc) > 1e-12, 2), :);
end

% =====================================================================
% Absent-grid-field route: pre-edit baseline parity + info.grid absent
% =====================================================================
function test_absent_route_matches_inline_legacy_baseline(testCase)
ion = invz_ion();
grid = [4 4 4]; dpRng = 6;
qcRef = legacy_qc(ion, grid);
[JnuRef, infoRef] = invz_jq_modes(ion, qcRef, struct('dpRng', dpRng, 'cache', false));
JnuRef = JnuRef(:);
Jaa0Ref = ion.Jxx0; if isfield(infoRef, 'Jaa0'), Jaa0Ref = infoRef.Jaa0; end

[JnuProd, infoProd, Jaa0Prod] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng, 'cache', false));

verifyEqual(testCase, JnuProd, JnuRef);
verifyEqual(testCase, Jaa0Prod, Jaa0Ref);
verifyEqual(testCase, infoProd, infoRef);
verifyFalse(testCase, isfield(infoProd, 'grid'), 'the absent-grid-field route must keep info.grid ABSENT.');
end

function test_absent_dipole_matches_explicit_bruteforce_legacy_parity(testCase)
ion = invz_ion();
grid = [4 4 4]; dpRng = 6;
[Jnu1, info1, Jaa01] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng, 'cache', false));
[Jnu2, info2, Jaa02] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng, 'cache', false, 'dipole', 'bruteforce'));
verifyEqual(testCase, Jnu2, Jnu1);
verifyEqual(testCase, info2, info1);
verifyEqual(testCase, Jaa02, Jaa01);
verifyFalse(testCase, isfield(info1, 'grid'));
verifyFalse(testCase, isfield(info2, 'grid'));
end

% =====================================================================
% Presence routing: gridConvention / gridOffset / gammaPolicy each alone
% activate invz_phase1_qgrid; Ewald alone does not; legacy_inclusive still
% activates the new (non-bit-identical) route
% =====================================================================
function test_gridConvention_alone_activates_new_route(testCase)
ion = invz_ion();
[~, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen'));
verifyTrue(testCase, isfield(info, 'grid'));
verifyEqual(testCase, info.grid.convention, 'halfopen');
verifyEqual(testCase, info.grid.offset, logical([0 0 0]));    % default gridOffset
verifyEqual(testCase, info.grid.gammaPolicy, 'P_drop');       % default gammaPolicy
end

function test_gridOffset_alone_activates_new_route(testCase)
ion = invz_ion();
[~, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridOffset', [1 0 0]));
verifyTrue(testCase, isfield(info, 'grid'));
verifyEqual(testCase, info.grid.offset, logical([1 0 0]));
verifyEqual(testCase, info.grid.convention, 'legacy_inclusive');   % default gridConvention
verifyEqual(testCase, info.grid.gammaPolicy, 'P_drop');            % default gammaPolicy
end

function test_gammaPolicy_alone_activates_new_route(testCase)
ion = invz_ion();
[~, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gammaPolicy', 'P_complete'));
verifyTrue(testCase, isfield(info, 'grid'));
verifyEqual(testCase, info.grid.gammaPolicy, 'P_complete');
verifyEqual(testCase, info.grid.convention, 'legacy_inclusive');   % default gridConvention
verifyEqual(testCase, info.grid.offset, logical([0 0 0]));         % default gridOffset
end

function test_ewald_alone_does_not_activate_new_route(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
verifyFalse(testCase, isfield(info, 'grid'), ...
    'opts.dipole=''ewald'' alone (no grid-policy field) must NOT activate the new grid route.');
verifyEqual(testCase, info.dipole.backend, 'ewald');
verifyTrue(testCase, all(isfinite(Jnu(:))));
end

function test_ewald_plus_new_grid_route_can_combine(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo, 'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_complete'));
verifyTrue(testCase, isfield(info, 'grid'));
verifyEqual(testCase, info.grid.convention, 'halfopen');
verifyEqual(testCase, info.grid.gammaPolicy, 'P_complete');
verifyEqual(testCase, info.dipole.backend, 'ewald');
verifyTrue(testCase, all(isfinite(Jnu(:))));
end

function test_gridConvention_legacy_inclusive_activates_new_route_not_bitidentical(testCase)
ion = invz_ion();
grid = [4 4 4]; dpRng = 6;
[JnuAbsent, infoAbsent] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng, 'cache', false));
verifyFalse(testCase, isfield(infoAbsent, 'grid'));

[JnuNew, infoNew] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng, 'cache', false, ...
    'gridConvention', 'legacy_inclusive'));
verifyTrue(testCase, isfield(infoNew, 'grid'));
verifyEqual(testCase, infoNew.grid.convention, 'legacy_inclusive');

% NOT claimed bit-identical: invz_phase1_qgrid wraps every point via
% mod(q+0.5,1)-0.5 (its header, and docs/invzp_ewald_integration_map.md Sec.5.7
% Risk 2), folding the endpoint-inclusive convention's +0.5 face onto -0.5 --
% a genuinely different retained q-set than the raw, unwrapped qVec_generator
% + abs(q)>1e-12 Gamma-drop the absent route uses.
notBitIdentical = ~isequal(size(JnuAbsent), size(JnuNew)) || ~isequal(JnuAbsent, JnuNew);
verifyTrue(testCase, notBitIdentical, ...
    'explicit gridConvention=''legacy_inclusive'' through the new route must not be bit-identical to the absent route.');
end

% =====================================================================
% Non-cubic grid rejection on the new route
% =====================================================================
function test_noncubic_grid_rejected_on_new_route(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_bz_couplings(ion, struct('grid', [4 4 6], 'dpRng', 6, ...
    'cache', false, 'gammaPolicy', 'P_drop')), 'invz:bzCouplingsAnisotropicHalfopen');
verifyError(testCase, @() invz_bz_couplings(ion, struct('grid', [4 6 4], 'dpRng', 6, ...
    'cache', false, 'gridConvention', 'halfopen')), 'invz:bzCouplingsAnisotropicHalfopen');
end

function test_noncubic_grid_absent_route_still_unaffected(testCase)
% Sanity: the SAME anisotropic grid is perfectly legal on the absent-grid-field
% (legacy) route -- only the new grid-policy route is cubic-only.
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 6], 'dpRng', 6, 'cache', false));
verifyFalse(testCase, isfield(info, 'grid'));
verifyTrue(testCase, all(isfinite(Jnu(:))));
end

% =====================================================================
% Complete info.grid field presence/types (new route only)
% =====================================================================
function test_info_grid_complete_fields_and_types(testCase)
ion = invz_ion();
N = 4; dpRng = 6;
[~, info] = invz_bz_couplings(ion, struct('grid', [N N N], 'dpRng', dpRng, 'cache', false, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_complete'));
verifyTrue(testCase, isfield(info, 'grid'));
gr = info.grid;
expected = sort({'schema', 'convention', 'offset', 'gammaPolicy', 'requested', 'nominal', ...
                 'retained', 'n_gamma', 'qhash'});
verifyEqual(testCase, sort(fieldnames(gr)), expected(:));

verifyTrue(testCase, (ischar(gr.schema) || isstring(gr.schema)) && ~isempty(gr.schema));
verifyEqual(testCase, gr.convention, 'halfopen');
verifyTrue(testCase, islogical(gr.offset) && isequal(size(gr.offset), [1 3]));
verifyEqual(testCase, gr.offset, logical([0 0 0]));
verifyEqual(testCase, gr.gammaPolicy, 'P_complete');
verifyEqual(testCase, gr.requested, [N N N]);
verifyEqual(testCase, gr.nominal, N^3);
verifyEqual(testCase, gr.retained, N^3);            % P_complete: nothing dropped
verifyTrue(testCase, isnumeric(gr.n_gamma) && isscalar(gr.n_gamma) && gr.n_gamma >= 1);
verifyTrue(testCase, ischar(gr.qhash) && numel(gr.qhash) == 64, 'SHA-256 must be a 64-hex-char digest.');
verifyTrue(testCase, all(ismember(lower(gr.qhash), '0123456789abcdef')));
verifyEqual(testCase, info.dpRng, dpRng);   % dpRng still forwarded through the new route
end

% =====================================================================
% P_complete vs P_drop row-count/provenance divergence
% =====================================================================
function test_pcomplete_vs_pdrop_row_counts_and_provenance(testCase)
ion = invz_ion();
N = 4; grid = [N N N]; dpRng = 6;
common = struct('grid', grid, 'dpRng', dpRng, 'cache', false, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0]);
optsComplete = common; optsComplete.gammaPolicy = 'P_complete';
optsDrop     = common; optsDrop.gammaPolicy     = 'P_drop';

[~, infoComplete] = invz_bz_couplings(ion, optsComplete);
[~, infoDrop]     = invz_bz_couplings(ion, optsDrop);

verifyEqual(testCase, infoComplete.grid.nominal, N^3);
verifyEqual(testCase, infoDrop.grid.nominal, N^3);
verifyEqual(testCase, infoComplete.grid.n_gamma, infoDrop.grid.n_gamma);
verifyGreaterThanOrEqual(testCase, infoDrop.grid.n_gamma, 1, ...
    'the [0 0 0]-offset halfopen grid must contain the exact Gamma row.');
verifyEqual(testCase, infoComplete.grid.retained, N^3);
verifyEqual(testCase, infoDrop.grid.retained, N^3 - infoDrop.grid.n_gamma);
verifyNotEqual(testCase, infoComplete.grid.qhash, infoDrop.grid.qhash, ...
    'different retained q row sets must hash differently.');
end

% =====================================================================
% cacheContext separation: legacy_bz / phase1_qgrid / direct_call must never
% cross-hit, even when the underlying qvec/physical payload coincide exactly.
% =====================================================================
function test_cachecontext_kind_prevents_false_cross_hit(testCase)
ion = invz_ion();
grid = [4 4 4]; dpRng = 17;   % distinctive dpRng: isolates this test's own cache file(s)
cacheDir = testCase.TestData.cacheDir;
pattern = sprintf('jq5_bruteforce_%d_*.mat', dpRng);
local_clear_matching(cacheDir, pattern);

% direct_call route: writes a file whose (backend,dpRng,qvec,pkey) digest is
% EXACTLY what invz_bz_couplings' absent-grid route independently computes
% below (identical qVec_generator + Gamma-drop construction, identical ion/dpRng).
qcDirect = legacy_qc(ion, grid);
[JnuTrue, ~] = invz_jq_modes(ion, qcDirect, struct('dpRng', dpRng, 'cache', true));

d = dir(fullfile(cacheDir, pattern));
verifyEqual(testCase, numel(d), 1, 'expected exactly one direct_call cache file.');
fpath = fullfile(cacheDir, d(1).name);
S = load(fpath);
verifyEqual(testCase, S.cacheMeta.cacheContext.kind, 'direct_call');

% Corrupt the stored payload with an obviously-wrong sentinel while leaving
% cacheMeta.cacheContext at its direct_call identity: if invz_bz_couplings'
% own, DIFFERENT cacheContext.kind='legacy_bz' were ever ignored, the call
% below would silently return this -999 sentinel instead of recomputing.
S.Jnu = zeros(size(S.Jnu)) - 999;
save(fpath, '-struct', 'S');

[JnuBz, infoBz] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', dpRng, 'cache', true));
verifyEqual(testCase, JnuBz, JnuTrue(:), 'AbsTol', 0);
verifyFalse(testCase, any(JnuBz(:) == -999), ...
    'legacy_bz cacheContext must never false-hit a direct_call-identity cache entry.');
verifyFalse(testCase, isfield(infoBz, 'grid'));

S2 = load(fpath);
verifyEqual(testCase, S2.cacheMeta.cacheContext.kind, 'legacy_bz', ...
    'the recompute must overwrite the file with invz_bz_couplings'' OWN cacheContext identity.');
end

function test_cachecontext_phase1_qgrid_cold_warm_and_distinct_identity(testCase)
ion = invz_ion();
grid = [4 4 4]; dpRng = 19;   % distinctive dpRng: isolates this test's own cache file(s)
cacheDir = testCase.TestData.cacheDir;
pattern = sprintf('jq5_bruteforce_%d_*.mat', dpRng);
local_clear_matching(cacheDir, pattern);

optsNew = struct('grid', grid, 'dpRng', dpRng, 'cache', true, 'gridConvention', 'halfopen', ...
    'gridOffset', [0 0 0], 'gammaPolicy', 'P_drop');
[JnuCold, infoCold] = invz_bz_couplings(ion, optsNew);
[JnuWarm, infoWarm] = invz_bz_couplings(ion, optsNew);
verifyEqual(testCase, JnuWarm, JnuCold);
verifyEqual(testCase, infoWarm, infoCold);

d = dir(fullfile(cacheDir, pattern));
verifyGreaterThanOrEqual(testCase, numel(d), 1);
found = false;
for i = 1:numel(d)
    S = load(fullfile(cacheDir, d(i).name));
    if isfield(S, 'cacheMeta') && isfield(S.cacheMeta, 'cacheContext') ...
            && isfield(S.cacheMeta.cacheContext, 'kind') ...
            && strcmp(S.cacheMeta.cacheContext.kind, 'phase1_qgrid')
        found = true;
        verifyEqual(testCase, S.cacheMeta.cacheContext.convention, 'halfopen');
        verifyEqual(testCase, S.cacheMeta.cacheContext.gammaPolicy, 'P_drop');
        verifyEqual(testCase, S.cacheMeta.cacheContext.offset, logical([0 0 0]));
    end
end
verifyTrue(testCase, found, 'expected at least one on-disk cacheContext.kind=phase1_qgrid entry.');
end
