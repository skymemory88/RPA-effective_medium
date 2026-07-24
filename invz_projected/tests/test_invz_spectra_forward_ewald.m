function tests = test_invz_spectra_forward_ewald
%TEST_INVZ_SPECTRA_FORWARD_EWALD Step-5 Task 7: invz_spectra_map/invz_spectra_qpath backend/
% grid forwarding and strict precomputed-coupling provenance/conflict rules
% (docs/superpowers/plans/2026-07-24-ewald-step5-integration.md Task 7; authority
% docs/invzp_ewald_integration_map.md Sec.5A).
%
% Covers, for BOTH drivers unless noted: partial-precomputed-pair rejection; legacy
% provenance-less precomputed acceptance (backward compat); invalid opts.dipole/opts.ewald
% caught even when precomputed inputs bypass invz_jq_modes entirely; an explicit backend/grid
% request against precomputed inputs lacking that provenance; matching explicit requests (no
% error); conflicting explicit backend (invz:spectraBackendConflict) and grid policy (its own
% distinct id); the computed route forwarding backend+grid to invz_bz_couplings; both Gamma
% policies; and (q-path only) backend inheritance from precomputed provenance, the
% isequaln(P.dipole, info.dipole) enforcement, and S.path_dipole exposure.
%
% Deliberately CHEAP and MECHANICAL: tiny grids/dpRng, short w/qpath arrays, a single
% comfortably-paramagnetic field (5.5 T at 0.31 K, matching the existing spectra test
% fixtures). No Jensen/ordered-phase/masking-symptom assertion anywhere in this file.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
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

function w = cheap_w()
w = (0.05:0.05:0.2).';
end

function [Jnu, info] = legacy_precomputed()
% The pre-Task-7 synthetic fixture already used throughout test_invz_spectra_map.m /
% test_invz_spectra_qpath.m: provenance-less (no .dipole, no .grid) -- the backward-
% compatibility case that MUST keep working unchanged.
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
end

% =====================================================================
% Partial precomputed pair (opts.Jnu XOR opts.info) -- BOTH drivers
% =====================================================================
function test_map_partial_precomputed_pair_errors(testCase)
ion = invz_ion();  w = cheap_w();
[Jnu, info] = legacy_precomputed();
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, struct('Jnu', Jnu, 'verbose', false)), ...
    'invz:spectraPrecomputedPartial');
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, struct('info', info, 'verbose', false)), ...
    'invz:spectraPrecomputedPartial');
end

function test_qpath_partial_precomputed_pair_errors(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
[Jnu, info] = legacy_precomputed();
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('Jnu', Jnu)), ...
    'invz:spectraPrecomputedPartial');
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('info', info)), ...
    'invz:spectraPrecomputedPartial');
end

% =====================================================================
% Legacy provenance-less precomputed input: accepted unchanged, resolves to brute force
% (the backward-compatibility case exercised by every pre-existing spectra test)
% =====================================================================
function test_map_legacy_provenanceless_precomputed_accepted(testCase)
ion = invz_ion();  w = cheap_w();
[Jnu, info] = legacy_precomputed();
S = invz_spectra_map(ion, 0.31, 5.5, w, struct('Jnu', Jnu, 'info', info, 'verbose', false));
verifyEqual(testCase, S.info, info);
verifySize(testCase, S.chiz, [numel(w) 1]);
end

function test_qpath_legacy_provenanceless_precomputed_accepted(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
[Jnu, info] = legacy_precomputed();
S = invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('Jnu', Jnu, 'info', info));
verifyEqual(testCase, S.info, info);
verifyEqual(testCase, S.path_dipole.backend, 'bruteforce');
verifySize(testCase, S.chiz, [numel(w) 2]);
end

% =====================================================================
% Invalid opts.dipole/opts.ewald must be rejected even though a precomputed opts.Jnu/opts.info
% bypasses invz_jq_modes entirely (no lattice sum ever runs in either assertion below)
% =====================================================================
function test_map_invalid_ewald_controls_rejected_even_precomputed(testCase)
ion = invz_ion();  w = cheap_w();
[Jnu, info] = legacy_precomputed();
eo = default_eopts(ion.a);
eoBad = eo;  eoBad.boundary = 'open_surface';
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'dipole', 'ewald', 'ewald', eoBad)), ...
    'invz:jqModesEwaldBoundary');
eoMissing = rmfield(eo, 'g_cut');
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'dipole', 'ewald', 'ewald', eoMissing)), ...
    'invz:jqModesEwaldOptsFields');
end

function test_qpath_invalid_ewald_controls_rejected_even_precomputed(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
[Jnu, info] = legacy_precomputed();
eo = default_eopts(ion.a);
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'ewald', eo)), ...          % ewald opts without dipole='ewald'
    'invz:jqModesEwaldOptsUnexpected');
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'dipole', 'bogus')), ...
    'invz:jqModesBackend');
end

% =====================================================================
% Explicit backend/grid request against precomputed inputs that lack that provenance
% =====================================================================
function test_map_explicit_request_missing_provenance_errors(testCase)
ion = invz_ion();  w = cheap_w();
[Jnu, info] = legacy_precomputed();   % no .dipole, no .grid
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'dipole', 'bruteforce')), ...
    'invz:spectraBackendProvenanceMissing');
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'gammaPolicy', 'P_complete')), ...
    'invz:spectraGridProvenanceMissing');
end

function test_qpath_explicit_request_missing_provenance_errors(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
[Jnu, info] = legacy_precomputed();
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'dipole', 'bruteforce')), ...
    'invz:spectraBackendProvenanceMissing');
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'gridOffset', [0 0 0])), ...
    'invz:spectraGridProvenanceMissing');
end

% =====================================================================
% Matching explicit backend+grid request against complete precomputed provenance: no error
% (also exercises the precomputed-complete route itself)
% =====================================================================
function test_map_matching_explicit_request_no_conflict(testCase)
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_complete'));
w = cheap_w();
opts = struct('Jnu', Jnu, 'info', info, 'verbose', false, 'grid', [4 4 4], ...
    'dipole', 'bruteforce', 'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_complete');
S = invz_spectra_map(ion, 0.31, 5.5, w, opts);
verifyEqual(testCase, S.info, info);
end

function test_qpath_matching_explicit_request_no_conflict_path_dipole_matches(testCase)
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_complete'));
w = cheap_w();  qpath = [1 0 0; 2 0 0];
opts = struct('Jnu', Jnu, 'info', info, 'grid', [4 4 4], 'dipole', 'bruteforce', ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_complete');
S = invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, opts);
verifyEqual(testCase, S.info, info);
verifyEqual(testCase, S.path_dipole.backend, 'bruteforce');
verifyTrue(testCase, isequaln(S.path_dipole, S.info.dipole));
end

% =====================================================================
% Conflicting explicit backend -> invz:spectraBackendConflict (backend-string conflict, and
% same-backend-different-Ewald-controls conflict)
% =====================================================================
function test_map_conflicting_explicit_backend_errors(testCase)
ion = invz_ion();  w = cheap_w();
eo = default_eopts(ion.a);
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false));  % bruteforce
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'dipole', 'ewald', 'ewald', eo)), ...
    'invz:spectraBackendConflict');

[JnuE, infoE] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
eo2 = eo;  eo2.r_cut = eo.r_cut * 1.1;   % same backend string, different Ewald controls
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', JnuE, 'info', infoE, 'verbose', false, 'dipole', 'ewald', 'ewald', eo2)), ...
    'invz:spectraBackendConflict');
end

function test_qpath_conflicting_explicit_backend_errors(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
eo = default_eopts(ion.a);
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false));  % bruteforce
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'dipole', 'ewald', 'ewald', eo)), ...
    'invz:spectraBackendConflict');

[JnuE, infoE] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
eo2 = eo;  eo2.alpha = eo.alpha * 1.1;
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', JnuE, 'info', infoE, 'dipole', 'ewald', 'ewald', eo2)), ...
    'invz:spectraBackendConflict');
end

% =====================================================================
% Conflicting grid policy -> its own distinct id (invz:spectraGridConflict)
% =====================================================================
function test_map_conflicting_grid_policy_errors(testCase)
ion = invz_ion();  w = cheap_w();
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen'));
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'grid', [4 4 4], ...
           'gridConvention', 'legacy_inclusive')), ...
    'invz:spectraGridConflict');
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, ...
    struct('Jnu', Jnu, 'info', info, 'verbose', false, 'grid', [6 6 6], ...
           'gridConvention', 'halfopen')), ...
    'invz:spectraGridConflict');
end

function test_qpath_conflicting_grid_policy_errors(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'gridConvention', 'halfopen'));
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'grid', [4 4 4], 'gridConvention', 'legacy_inclusive')), ...
    'invz:spectraGridConflict');
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, ...
    struct('Jnu', Jnu, 'info', info, 'grid', [6 6 6], 'gridConvention', 'halfopen')), ...
    'invz:spectraGridConflict');
end

% =====================================================================
% Computed route: backend + grid policy actually reach invz_bz_couplings
% =====================================================================
function test_map_computed_route_forwards_backend_and_grid(testCase)
ion = invz_ion();  w = cheap_w();
opts = struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, 'verbose', false, ...
    'dipole', 'bruteforce', 'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_complete');
S = invz_spectra_map(ion, 0.31, 5.5, w, opts);
verifyEqual(testCase, S.info.dipole.backend, 'bruteforce');
verifyEqual(testCase, S.info.grid.convention, 'halfopen');
verifyEqual(testCase, S.info.grid.gammaPolicy, 'P_complete');
verifyEqual(testCase, S.info.grid.offset, logical([0 0 0]));
verifyEqual(testCase, S.info.grid.requested, [4 4 4]);
end

function test_qpath_computed_route_forwards_backend_and_grid_path_dipole_matches(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
opts = struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'bruteforce', 'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_drop');
S = invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, opts);
verifyEqual(testCase, S.info.dipole.backend, 'bruteforce');
verifyEqual(testCase, S.info.grid.convention, 'halfopen');
verifyEqual(testCase, S.info.grid.gammaPolicy, 'P_drop');
verifyEqual(testCase, S.path_dipole.backend, 'bruteforce');
verifyTrue(testCase, isequaln(S.path_dipole, S.info.dipole));
end

% =====================================================================
% Both Gamma policies forwarded on the computed route
% =====================================================================
function test_map_both_gamma_policies_forwarded(testCase)
ion = invz_ion();  w = cheap_w();
base = struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, 'verbose', false, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0]);
oC = base;  oC.gammaPolicy = 'P_complete';
oD = base;  oD.gammaPolicy = 'P_drop';
Sc = invz_spectra_map(ion, 0.31, 5.5, w, oC);
Sd = invz_spectra_map(ion, 0.31, 5.5, w, oD);
verifyEqual(testCase, Sc.info.grid.gammaPolicy, 'P_complete');
verifyEqual(testCase, Sd.info.grid.gammaPolicy, 'P_drop');
verifyEqual(testCase, Sc.info.grid.retained, Sc.info.grid.nominal);
verifyEqual(testCase, Sd.info.grid.retained, Sd.info.grid.nominal - Sd.info.grid.n_gamma);
end

% =====================================================================
% q-path: complete Ewald provenance on a precomputed input, with NO explicit request, is
% INHERITED (not silently dropped to brute force) -- the "precomputed-complete route" bullet
% =====================================================================
function test_qpath_precomputed_ewald_provenance_inherited(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [4 4 4], 'dpRng', 6, 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
w = cheap_w();  qpath = [1 0 0; 2 0 0];
S = invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('Jnu', Jnu, 'info', info));
verifyEqual(testCase, S.path_dipole.backend, 'ewald');
verifyTrue(testCase, isequaln(S.path_dipole, S.info.dipole));
end

% =====================================================================
% q-path: a stale/corrupted (but structurally complete) info.dipole must be caught, never
% silently trusted -- invz:spectraPathDipoleMismatch
% =====================================================================
function test_qpath_path_dipole_mismatch_errors(testCase)
ion = invz_ion();  w = cheap_w();  qpath = [1 0 0; 2 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
info = struct('Jcc0', 6.4e-3, 'dipole', struct('backend', 'bruteforce', ...
    'ewald', struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', ''), ...
    'q_reduction', 'STALE-CORRUPTED-PROVENANCE-NOT-THE-REAL-STRING', ...
    'primitive_schema', 'MF_dipole+exchange (legacy, unversioned)'));
verifyError(testCase, @() invz_spectra_qpath(ion, 0.31, 5.5, qpath, w, struct('Jnu', Jnu, 'info', info)), ...
    'invz:spectraPathDipoleMismatch');
end
