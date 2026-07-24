function tests = test_invz_jq_path_ewald
%TEST_INVZ_JQ_PATH_EWALD Step-5 Task 6: invz_jq_path Gamma-metadata base and
% backend-conditional Gamma behavior (docs/superpowers/plans/2026-07-24-
% ewald-step5-integration.md Task 6; docs/invzp_ewald_integration_map.md Sec.6).
%
% Covers:
%   1. Brute-force P-field parity vs the frozen invz_legacy_coupling_reference
%      oracle, on a q-path explicitly covering an ordinary point, a point
%      INSIDE the old snap (trust-radius) band, a point OUTSIDE the old snap
%      band, exact Gamma, a Gamma-equivalent [2 0 0] endpoint, and a second
%      ordinary point -- non-vacuously confirmed (each row's snap/no-snap
%      outcome is checked against the actual trust-radius arithmetic, not
%      merely inferred from reference==production agreement).
%   2. The pre-Task-6 locally-rebuilt Greg (MF_dipole/exchange 5-arg reuse
%      form, from an absent-backend info that still carries geomD/geomX) is
%      isequaln to invz_jq_modes' brute-force info.Jpath_base_cc -- the
%      regression that justifies deleting invz_jq_path's own MF_dipole/
%      exchange reconstruction (this does not call invz_jq_path at all).
%   3. Forwarding by presence: an explicit opts.dipole='bruteforce' request
%      reproduces the absent-backend P bit-for-bit (every field, including
%      the additive P.dipole provenance).
%   4. A tiny nonzero Ewald q, well inside the FORMER brute-force trust
%      radius, is evaluated DIRECTLY (P.snapped=false), matching a direct
%      invz_jq_modes Ewald call bit-for-bit -- proving Ewald never uses the
%      dpRng-derived snap band.
%   5. Exact Gamma under Ewald (single-point path, in-plane default
%      direction) is still replaced by the directional limit (P.snapped=
%      true), P.ksnap is NaN, and the numeric value matches a hand-built
%      reconstruction from a live invz_jq_modes Ewald call.
%   6. The same, but approached along c* (kz2=1, genuinely anisotropic, not
%      the degenerate in-plane case) -- confirms the LOCAL PATH DIRECTION is
%      actually used under Ewald, not a hardcoded isotropic value.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                        % invz_projected: invz_jq_path, invz_jq_modes
addpath(fullfile(here, '..', '..'));                  % repo root: MF_dipole, exchange, invz_dipole_ewald
addpath(fullfile(here, '..', '..', 'invz_common'));   % invz_ion, invz_const

fx = invz_ewald_fixtures();
% Frozen production-default Ewald controls (same C_r=5.5, C_g=11 combo as
% test_invz_ewald_gateC_integration.m's mk_eopts).
testCase.TestData.eo = struct('alpha', fx.alpha0, 'r_cut', 5.5/fx.alpha0, ...
    'g_cut', 11*fx.alpha0, 'boundary', 'conducting_k0_omitted');
end

% =====================================================================
% shared fixtures/helpers (names free of "test" so functiontests skips them)
% =====================================================================
function q = bruteforce_parity_qpath()
% Row-by-row intended category (LiHoF4 production ion = invz_ion(); default
% dpRng=30, snapfac=2.5 -> ksnap = 2.5*2*pi/(30*5.175) = 0.10118 Ang^-1):
%   1  ordinary              generic point, far from any Gamma-equivalent K
%   2  INSIDE old snap band  offset 0.02 r.l.u. along a* -> |k|=0.0243 < ksnap
%   3  OUTSIDE old snap band offset 0.20 r.l.u. along a* -> |k|=0.2428 > ksnap
%      (still rounds to Gamma, so the guard's Gamma-equivalence gate fires
%      but the trust-radius gate then declines to replace -- this exercises
%      the "outside the band" branch, not merely "not near any integer")
%   4  exact Gamma           q = [0 0 0] exactly
%   5  Gamma-equivalent      q = [2 0 0] exactly (structure factor 1 for the
%                            LiHoF4 basis, tau(:,1) = [0 0 .5 .5])
%   6  ordinary               a second generic point, well outside the band
q = [ ...
    0.310  -0.180   0.240
    0.020   0        0
    0.200   0        0
    0       0        0
    2.000   0        0
   -0.450   0.330    0.150];
end

function verify_P_fields_isequaln(testCase, P_ref, P_prod, ctx)
% Every field of P_ref must be present in P_prod and isequaln to it. P_prod
% may additively carry fields P_ref lacks (e.g. P.dipole vs the frozen
% oracle) -- that is not a regression, mirroring test_invz_ewald_baselines.m's
% verify_struct_fields_isequaln contract.
fn_ref  = sort(fieldnames(P_ref));
missing = setdiff(fn_ref, fieldnames(P_prod));
verifyEmpty(testCase, missing, sprintf( ...
    '%s: production is MISSING legacy P field(s) {%s}.', ctx, strjoin(missing(:).', ', ')));
for i = 1:numel(fn_ref)
    f = fn_ref{i};
    if ~isfield(P_prod, f), continue; end   % already reported as missing above
    verifyTrue(testCase, isequaln(P_ref.(f), P_prod.(f)), ...
        sprintf('P.%s differs between reference and production (%s).', f, ctx));
end
end

% =====================================================================
% 1. brute-force P-field parity vs the frozen ref.jq_path oracle
% =====================================================================
function test_bruteforce_P_fields_match_frozen_reference(testCase)
ion = invz_ion();
qpath = bruteforce_parity_qpath();
ref = invz_legacy_coupling_reference();
opts = struct('dpRng', 30, 'cache', false);

P_ref  = ref.jq_path(ion, qpath, opts);
P_prod = invz_jq_path(ion, qpath, opts);

% Non-vacuous confirmation that the array actually exercises every branch it
% claims to, independent of any ref==prod agreement on whatever branch
% happened to run.
Brec  = 2*pi*inv(ion.a).';
ksnap = 2.5*2*pi/(30*min(vecnorm(ion.a, 2, 2)));
k2 = norm((qpath(2,:) - round(qpath(2,:)))*Brec);
k3 = norm((qpath(3,:) - round(qpath(3,:)))*Brec);
verifyTrue(testCase, k2 < ksnap,  'row 2 fixture is not actually inside the old snap band.');
verifyTrue(testCase, k3 >= ksnap, 'row 3 fixture is not actually outside the old snap band.');

verifyFalse(testCase, P_prod.snapped(1), 'row 1 (ordinary) must not be snapped.');
verifyTrue(testCase,  P_prod.snapped(2), 'row 2 (inside old snap band) must be snapped in production.');
verifyFalse(testCase, P_prod.snapped(3), 'row 3 (outside old snap band) must NOT be snapped in production.');
verifyTrue(testCase,  P_prod.snapped(4), 'row 4 (exact Gamma) must be snapped in production.');
verifyTrue(testCase,  P_prod.snapped(5), 'row 5 (Gamma-equivalent [2 0 0]) must be snapped in production.');
verifyFalse(testCase, P_prod.snapped(6), 'row 6 (ordinary) must not be snapped.');

verifyEqual(testCase, P_prod.dipole.backend, 'bruteforce', ...
    'P.dipole.backend must read ''bruteforce'' for an absent opts.dipole request.');

verify_P_fields_isequaln(testCase, P_ref, P_prod, 'bruteforce parity');
end

% =====================================================================
% 2. old locally-rebuilt Greg == invz_jq_modes' brute-force info.Jpath_base_cc
%    (does not call invz_jq_path; justifies deleting the local reconstruction)
% =====================================================================
function test_old_local_Greg_equals_Jpath_base_cc(testCase)
ion = invz_ion();
C = invz_const();
q = bruteforce_parity_qpath();
% Absent backend -> bruteforce; info still carries geomD/geomX (bruteforce-only
% fields), exactly as invz_jq_path used to consume before Task 6.
[~, info] = invz_jq_modes(ion, q, struct('dpRng', 30, 'cache', false));

dip0 = MF_dipole([0 0 0], 30, ion.a, ion.tau, info.geomD);
ex0  = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau, info.geomX);
Greg_old = -squeeze(C.gfac*dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));

verifyTrue(testCase, isequaln(Greg_old, info.Jpath_base_cc), ...
    'old locally-rebuilt Greg is not isequaln to info.Jpath_base_cc -- deletion is not justified.');
end

% =====================================================================
% 3. forwarding by presence: explicit 'bruteforce' == absent backend
% =====================================================================
function test_explicit_bruteforce_matches_absent_backend(testCase)
ion = invz_ion();
qpath = bruteforce_parity_qpath();
P_absent   = invz_jq_path(ion, qpath, struct('dpRng', 30, 'cache', false));
P_explicit = invz_jq_path(ion, qpath, struct('dpRng', 30, 'cache', false, 'dipole', 'bruteforce'));

verifyEqual(testCase, P_explicit.dipole.backend, 'bruteforce');
verify_P_fields_isequaln(testCase, P_absent, P_explicit, 'absent vs explicit bruteforce');
end

% =====================================================================
% 4. Ewald: tiny nonzero q inside the FORMER brute trust radius is evaluated
%    directly, never snapped
% =====================================================================
function test_ewald_tiny_q_below_old_trust_radius_not_snapped(testCase)
ion = invz_ion();
eo  = testCase.TestData.eo;
qtiny = [1e-4 0 0];    % |k| ~ 1.2e-4 Ang^-1, well inside the old ksnap ~ 0.101

Brec = 2*pi*inv(ion.a).';
knorm    = norm(qtiny*Brec);
oldKsnap = 2.5*2*pi/(30*min(vecnorm(ion.a, 2, 2)));
verifyTrue(testCase, knorm < oldKsnap, ...
    'qtiny fixture is not actually inside the FORMER brute trust radius.');

ewaldOpts = struct('dipole', 'ewald', 'ewald', eo, 'cache', false);
P = invz_jq_path(ion, qtiny, ewaldOpts);
verifyFalse(testCase, P.snapped(1), 'a tiny nonzero Ewald q must NOT be snapped (evaluated directly).');
verifyEqual(testCase, P.dipole.backend, 'ewald');
verifyTrue(testCase, isnan(P.ksnap), 'P.ksnap must be NaN under the Ewald backend.');

[JnuDirect, ~, JuniDirect] = invz_jq_modes(ion, qtiny, struct('dipole', 'ewald', 'ewald', eo, 'cache', false));
verifyTrue(testCase, isequaln(P.Jnu(1,:), JnuDirect(1,:)), ...
    'un-snapped Ewald row must match a direct invz_jq_modes evaluation bit-for-bit.');
verifyTrue(testCase, isequaln(P.Juni(1), JuniDirect(1)), ...
    'un-snapped Ewald Juni must match a direct invz_jq_modes evaluation bit-for-bit.');
end

% =====================================================================
% 5. Ewald: exact Gamma (single-point path, in-plane default) is snapped;
%    P.ksnap is NaN
% =====================================================================
function test_ewald_exact_gamma_single_point_snapped_ksnap_nan(testCase)
ion = invz_ion();
C   = invz_const();
eo  = testCase.TestData.eo;
ewaldOpts = struct('dipole', 'ewald', 'ewald', eo, 'cache', false);

P = invz_jq_path(ion, [0 0 0], ewaldOpts);
verifyTrue(testCase, P.snapped(1), 'exact Gamma under Ewald must still be snapped (local-direction limit).');
verifyTrue(testCase, isnan(P.ksnap), 'P.ksnap must be NaN under the Ewald backend.');
verifyEqual(testCase, P.dipole.backend, 'ewald');

% Cross-check the numeric value against the hand-built single-point-path
% in-plane default direction limit (kz2=0): Jm = Jpath_base_cc + gfac*(4pi/Vc)*(1/3).
[~, infoDirect] = invz_jq_modes(ion, [0 0 0], struct('dipole', 'ewald', 'ewald', eo, 'cache', false));
v  = ones(4,1)/2;
Jm = infoDirect.Jpath_base_cc + C.gfac*(4*pi/ion.Vc)*(1/3 - 0);
Jm = (Jm + Jm')/2;
expectedJnu  = sort(real(eig(Jm))).';
expectedJuni = real(v.'*Jm*v);
verifyTrue(testCase, isequaln(P.Jnu(1,:), expectedJnu), ...
    'exact-Gamma Ewald P.Jnu does not match the hand-built in-plane directional limit.');
verifyTrue(testCase, isequaln(P.Juni(1), expectedJuni), ...
    'exact-Gamma Ewald P.Juni does not match the hand-built in-plane directional limit.');
end

% =====================================================================
% 6. Ewald: exact Gamma approached along c* (kz2=1, anisotropic) -- confirms
%    the LOCAL PATH DIRECTION is genuinely used, not a hardcoded isotropic value
% =====================================================================
function test_ewald_exact_gamma_along_c_axis_direction_snapped(testCase)
ion = invz_ion();
C   = invz_const();
eo  = testCase.TestData.eo;
ewaldOpts = struct('dipole', 'ewald', 'ewald', eo, 'cache', false);

qpath = [0 0 -0.05; 0 0 0];    % approach Gamma purely along c* (out-of-plane)
P = invz_jq_path(ion, qpath, ewaldOpts);

verifyFalse(testCase, P.snapped(1), 'row 1 (finite offset along c*) must not be snapped (Ewald has no trust radius).');
verifyTrue(testCase, P.snapped(2), 'exact Gamma (row 2) must be snapped under Ewald.');
verifyTrue(testCase, isnan(P.ksnap), 'P.ksnap must be NaN under the Ewald backend.');

[~, infoDirect] = invz_jq_modes(ion, [0 0 0], struct('dipole', 'ewald', 'ewald', eo, 'cache', false));
v = ones(4,1)/2;
JmC = infoDirect.Jpath_base_cc + C.gfac*(4*pi/ion.Vc)*(1/3 - 1);   % kz2=1: pure c-axis approach
JmC = (JmC + JmC')/2;
expectedJnu  = sort(real(eig(JmC))).';
expectedJuni = real(v.'*JmC*v);
verifyTrue(testCase, isequaln(P.Jnu(2,:), expectedJnu), ...
    'exact-Gamma Ewald P.Jnu (c-axis approach) does not match the hand-built directional limit.');
verifyTrue(testCase, isequaln(P.Juni(2), expectedJuni), ...
    'exact-Gamma Ewald P.Juni (c-axis approach) does not match the hand-built directional limit.');

% The c-axis (kz2=1) and in-plane (kz2=0, test 5) directional limits must
% genuinely differ -- otherwise both tests could coincidentally be checking
% the same hardcoded isotropic value instead of the real direction-dependent
% formula actually being exercised.
JmInPlane = infoDirect.Jpath_base_cc + C.gfac*(4*pi/ion.Vc)*(1/3 - 0);
JmInPlane = (JmInPlane + JmInPlane')/2;
inPlaneJuni = real(v.'*JmInPlane*v);
verifyNotEqual(testCase, expectedJuni, inPlaneJuni, ...
    'the c-axis and in-plane directional limits coincide -- fixture does not distinguish direction-dependence.');
end
