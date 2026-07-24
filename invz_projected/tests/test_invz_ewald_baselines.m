function tests = test_invz_ewald_baselines
% Step-5 Task 1 (docs/superpowers/plans/2026-07-24-ewald-step5-integration.md):
% proves invz_legacy_coupling_reference.m -- a cache-free, test-only,
% independently-maintained copy of the pre-Step-5 brute-force invz_jq_modes
% non-ODD branch and invz_jq_path arithmetic -- is numerically IDENTICAL to
% the current, unmodified production invz_jq_modes/invz_jq_path at a small
% deterministic q array covering ordinary, near-Gamma, exact-Gamma, and
% Gamma-equivalent (nonzero-integer) points.
%
% This is the regression spine every later Step-5 task (which DOES modify
% invz_jq_modes.m/invz_jq_path.m) must continue to satisfy for an absent or
% explicit opts.dipole='bruteforce' request: Jnu, Juni, every info field,
% and every pre-existing P field must remain isequaln to this frozen
% reference. Task 1 itself does not modify invz_jq_modes.m/invz_jq_path.m,
% so this test currently compares the reference against itself in all but
% implementation (both are, today, the identical pre-Step-5 arithmetic);
% its value is as the frozen baseline later tasks regress against.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));            % invz_projected: invz_jq_modes, invz_jq_path
addpath(fullfile(here, '..', '..'));      % repo root: MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));  % invz_ion, invz_const
end

% =====================================================================
% shared deterministic q array (also reused as a q-path for invz_jq_path)
% (name deliberately free of "test" so functiontests skips it)
% =====================================================================
function q = qtest_array()
% Row-by-row intended category (LiHoF4 production ion = invz_ion()):
%   1  ordinary            generic point, far from any Gamma-equivalent K
%   2  near-Gamma          small nonzero offset from the origin
%   3  exact-Gamma         q = [0 0 0] exactly
%   4  Gamma-equivalent    q = [2 0 0] exactly: invz_is_gamma_equiv([2 0 0],tau)
%                          is true for the LiHoF4 basis (tau(:,1)=[0 0 .5 .5],
%                          so sum(exp(2i*pi*2*tau(:,1)))/4 == 1 exactly, i.e.
%                          [2 0 0] is a genuine structure-factor-1 Bragg point,
%                          unlike e.g. [1 0 0] whose structure factor is 0)
%   5  near a Gamma-equiv. small nonzero offset from the [2 0 0] point above
%   6  ordinary            a second generic point, well outside the snap band
q = [ ...
    0.310  -0.180   0.240
    0.020   0.010  -0.015
    0       0       0
    2.000   0       0
    1.985   0.010  -0.005
   -0.450   0.330   0.150];
end

% =====================================================================
% invz_jq_modes: Jnu, Juni, every info field (prod vs the frozen reference)
% =====================================================================
function test_jq_modes_reference_matches_production_nocache(testCase)
ion = invz_ion();
qvec = qtest_array();
ref = invz_legacy_coupling_reference();
opts = struct('dpRng', 30, 'cache', false);

[Jnu_ref, info_ref, Juni_ref]   = ref.jq_modes(ion, qvec, opts);
[Jnu_prod, info_prod, Juni_prod] = invz_jq_modes(ion, qvec, opts);

verifyTrue(testCase, isequaln(Jnu_ref, Jnu_prod),   'Jnu (cache=false) differs between reference and production.');
verifyTrue(testCase, isequaln(Juni_ref, Juni_prod), 'Juni (cache=false) differs between reference and production.');
verify_struct_fields_isequaln(testCase, info_ref, info_prod, 'info', 'cache=false');
end

function test_jq_modes_reference_matches_production_cached_cold_and_warm(testCase)
ion = invz_ion();
qvec = qtest_array();
ref = invz_legacy_coupling_reference();
opts = struct('dpRng', 30, 'cache', true);

% cold (first call may write a jq4_* cache file; harmless -- excluded from
% git by invz_projected/cache/.gitignore's '*' rule).
[Jnu_prod_cold, info_prod_cold, Juni_prod_cold] = invz_jq_modes(ion, qvec, opts);
% warm (second call should hit the just-written cache).
[Jnu_prod_warm, info_prod_warm, Juni_prod_warm] = invz_jq_modes(ion, qvec, opts);
verifyTrue(testCase, isequaln(Jnu_prod_cold, Jnu_prod_warm),   'production cold vs warm Jnu differ.');
verifyTrue(testCase, isequaln(Juni_prod_cold, Juni_prod_warm), 'production cold vs warm Juni differ.');
verify_struct_fields_isequaln(testCase, info_prod_cold, info_prod_warm, 'info', 'production cold vs warm');

% cache-free reference vs the (now warm) cached production call.
[Jnu_ref, info_ref, Juni_ref] = ref.jq_modes(ion, qvec, opts);
verifyTrue(testCase, isequaln(Jnu_ref, Jnu_prod_warm),   'Jnu (cache=true) differs between reference and warm production.');
verifyTrue(testCase, isequaln(Juni_ref, Juni_prod_warm), 'Juni (cache=true) differs between reference and warm production.');
verify_struct_fields_isequaln(testCase, info_ref, info_prod_warm, 'info', 'cache=true (reference vs warm production)');
end

% =====================================================================
% invz_jq_path: every pre-existing P field (prod vs the frozen reference)
% =====================================================================
function test_jq_path_reference_matches_production(testCase)
ion = invz_ion();
qpath = qtest_array();
ref = invz_legacy_coupling_reference();
opts = struct('dpRng', 30, 'cache', false);

P_ref  = ref.jq_path(ion, qpath, opts);
P_prod = invz_jq_path(ion, qpath, opts);

% Exact-Gamma (row 3) and exact Gamma-equivalent (row 4) points are
% MATHEMATICALLY GUARANTEED to fall inside the (always-positive) snap
% radius -- their offset from the nearest integer point is exactly zero --
% independent of any hand-computed boundary estimate for the other rows.
% This is a non-tautological confirmation that the array actually exercises
% the exact-Gamma[-equivalent] directional-limit branch, not merely that
% ref==prod on whatever branch happened to run.
verifyTrue(testCase, P_prod.snapped(3), 'row 3 (exact Gamma) must be snapped in production.');
verifyTrue(testCase, P_prod.snapped(4), 'row 4 (exact Gamma-equivalent [2 0 0]) must be snapped in production.');
verifyTrue(testCase, any(~P_prod.snapped), 'expected at least one ordinary (non-snapped) row in the array.');

verify_struct_fields_isequaln(testCase, P_ref, P_prod, 'P', 'jq_path');
end

% =====================================================================
% shared helper (name deliberately free of "test" so functiontests skips it)
% =====================================================================
function verify_struct_fields_isequaln(testCase, s_ref, s_prod, label, ctx)
% Field-SET equality first (fails loudly if the reference or production
% side is missing/gains a field), then isequaln on every individual field
% (fails loudly, and per-field, on any numeric/content divergence).
fn_ref  = sort(fieldnames(s_ref));
fn_prod = sort(fieldnames(s_prod));
verifyEqual(testCase, fn_ref, fn_prod, sprintf( ...
    '%s field-name set differs between reference and production (%s).', label, ctx));
for i = 1:numel(fn_prod)
    f = fn_prod{i};
    verifyTrue(testCase, isequaln(s_ref.(f), s_prod.(f)), sprintf( ...
        '%s.%s differs between reference and production (%s).', label, f, ctx));
end
end
