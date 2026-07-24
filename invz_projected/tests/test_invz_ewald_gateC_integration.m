function tests = test_invz_ewald_gateC_integration
%TEST_INVZ_EWALD_GATEC_INTEGRATION Step-5 Task 3: integration-level Gate-C4
% (full Cartesian) and Gate-C6 (demag), per
% docs/superpowers/plans/2026-07-24-ewald-step5-integration.md Task 3 and
% docs/invzp_ewald_prereg.md (FROZEN) sec 5 Gate-C checks 4/6 +
% sec 12 Errata E1. TEST-ONLY: this file introduces no production edit.
%
% GATE-C4 (full Cartesian). invz_jq_modes exposes only the cc coupling
% channel, so the full nine-component Cartesian tensor cannot be extracted
% from that API. All nine components are instead built here directly from
% the primitive invz_dipole_ewald plus the SAME caller-level normalization
% invz_jq_modes.m uses (-gfac, +Jex(0), -lorz*delta_ab); only the cc slice is
% cross-checked against a live invz_jq_modes call (info.Jpath_base_cc /
% info.Jgamma_cc), per the task brief.
%
%   1. Isolated-projector-with-coupling-sign check (prereg sec 5 item 3,
%      generalized to all nine components and scaled by -gfac): the
%      primitive's isolated G=[0,0,0] reciprocal contribution, extracted
%      from PUBLIC geom data exactly as
%      test_invz_dipole_ewald_gammaC.m's C3 check does (an independent
%      from-scratch reimplementation -- MATLAB local functions are
%      file-scoped, so that file's helper is unreachable from here), equals
%      the closed-form P0_ab(q)=(4*pi/Vc)*qhat_a*qhat_b*exp(-|q|^2/(4*alpha^2))
%      at M_id, BOTH sides scaled by -gfac. Both sides carry the SAME
%      exp(-|q|^2/(4*alpha^2)) factor, so this is a valid finite-s-vs-finite-s
%      M_id comparison -- it does not touch E1.
%   2. Full-Cartesian base-reconstruction algebra: Jbase_ab is built once
%      from dip_reg_ab(0)/Jex_ab(0) (both evaluated ONLY at exact Gamma, no
%      finite-s primitive value involved), and the reconstruction formula
%      Jbase_ab + gfac*(4*pi/Vc)*(delta_ab/3-qhat_a*qhat_b)*ones(ntau) is
%      verified (M_id, pure algebra) to reduce to the clean closed form
%      -gfac*dip_reg_ab(0)+Jex_ab(0)-gfac*(4*pi/Vc)*qhat_a*qhat_b*ones(ntau)
%      -- i.e. the delta_ab/3 Lorentz-cancellation is exact -- for all nine
%      components, all sublattice pairs, on the five frozen rays at
%      s in {+-1e-3,+-1e-4}.
%   3. cc cross-check: the cc slice of Jbase (built independently from the
%      primitive) is verified to equal info.Jpath_base_cc from a live
%      invz_jq_modes(...,'dipole','ewald',...) call, and the cc slice of the
%      full reconstruction is verified to equal the production q-path
%      formula info.Jpath_base_cc + gfac*(4*pi/Vc)*(1/3-qhat_c^2)*ones(ntau).
%
%   E1 COMPLIANCE (prereg sec 12): checks 2/3 NEVER compare a finite-s
%   primitive/projector value against the no-Gaussian small-q limit formula
%   at M_id -- Jbase_ab/info.Jpath_base_cc/info.Jgamma_cc are all evaluated
%   ONLY at exact Gamma (s=0, no qhat/exp ambiguity), and the reconstruction
%   formula is verified purely algebraically against its own closed-form
%   reduction, never against a directly-measured finite-s dip(q) value. This
%   sidesteps exactly the mistake E1 warns about (a genuine, non-vanishing
%   O(s^2) gap from the missing exp(-|q|^2/(4*alpha^2)) factor would appear
%   at M_id scale if the no-exp reconstruction were compared to a real
%   finite-s primitive value -- confirmed by hand derivation during test
%   design, not merely assumed).
%
% GATE-C6 (demag). For EACH backend (bruteforce, ewald), exactly three
% caller cases (off: demag=0; sphere: demag=1,alpha=1; needle: demag=1,
% alpha=0) are exercised on the SAME q array. Within each backend, Jnu,
% Jcc0, Sigma_c (invz_sigma_crit), and invz_critical_T0field(...) are
% required to be bit-identical across the three shapes -- demag/alpha never
% enter those computations (traced through invz_jq_modes.m: dm_cc/dm_aa feed
% ONLY Jshape_cc and Jaa0, never Jcc0_dipole/Jcc0/the per-q Jnu loop), so
% plain verifyEqual (no tolerance) is correct and mirrors the existing
% test_invz_demag_invariance.m precedent exactly. Jshape_cc and the
% demag-aware Jaa0 are independently checked against the analytic
% ellipsoid_demagn(...) formula (M_id). No brute/Ewald numerical-agreement
% tolerance is invented anywhere in this file: the one brute-vs-Ewald
% comparison present (Jshape_cc) is an EXACT algebraic identity by
% construction (computed once in invz_jq_modes.m's shared pre-dispatch code,
% strictly before the backend branch, reading no dip0/dip quantity at all),
% not a physics/backend-agreement claim. The Ewald primitive's eopts field
% set is confirmed to carry no demag/surface control.
%
% Authority: docs/invzp_ewald_prereg.md (FROZEN) sec 3 (M_id) + sec 5 (Gate C
% items 3/4/6) + sec 12 Errata E1; docs/invzp_ewald_integration_map.md sec
% 6.3 (Jpath_base_cc/Jgamma_cc caller normalization);
% docs/superpowers/plans/2026-07-24-ewald-step5-integration.md Task 2 ("Gamma
% metadata" formulas, already committed) and Task 3 (this file's checklist).
% Reused frozen helpers: invz_ewald_metrics.m (M_id), invz_ewald_fixtures.m
% (alpha0). invz_jq_path.m is NOT wired for Ewald yet (Task 6, later in the
% Step-5 plan) so it is never called here; the "production q-path formula"
% is verified directly against invz_jq_modes' exported Jpath_base_cc/
% Jgamma_cc fields plus hand-built primitive quantities, per the task
% brief's explicit instruction to cross-check the cc slice only through
% invz_jq_modes.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                        % invz_projected
addpath(fullfile(here, '..', '..'));                  % repo root: invz_dipole_ewald, exchange, ellipsoid_demagn
addpath(fullfile(here, '..', '..', 'invz_common'));   % invz_ion, invz_const

ion = invz_ion();
C   = invz_const();
fx  = invz_ewald_fixtures();
eo  = mk_eopts(fx.alpha0, 5.5/fx.alpha0, 11*fx.alpha0, 'conducting_k0_omitted');  % frozen production defaults

[rays, raylabels] = frozen_rays();
smags = frozen_smags();
[Q, rayOf, sOf] = ray_mag_grid(rays, smags);
Qall = [0 0 0; Q];                                    % row 1 = exact Gamma, rows 2:21 = the 20 ray/mag points

[dipAll, ~, geom0] = invz_dipole_ewald(Qall, ion.a, ion.tau, eo);
ex0 = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);

testCase.TestData.ion       = ion;
testCase.TestData.C         = C;
testCase.TestData.eo        = eo;
testCase.TestData.raylabels = raylabels;
testCase.TestData.Q         = Q;
testCase.TestData.rayOf     = rayOf;
testCase.TestData.sOf       = sOf;
testCase.TestData.dip0      = dipAll(:,:,:,:,1);      % exact-Gamma primitive tensor, [3,3,ntau,ntau]
testCase.TestData.geom0     = geom0;                  % q-independent geometry (valid at any q, same lattice/eopts)
testCase.TestData.ex0       = ex0;                    % exact-Gamma exchange tensor, [3,3,ntau,ntau]
end

% =====================================================================
% shared helpers (names deliberately free of "test" so functiontests skips
% them -- test_invz_dipole_ewald.m / test_invz_dipole_ewald_gammaC.m precedent)
% =====================================================================
function eo = mk_eopts(alpha, r_cut, g_cut, boundary)
eo = struct('alpha', alpha, 'r_cut', r_cut, 'g_cut', g_cut, 'boundary', boundary);
end

function [rays, raylabels] = frozen_rays()
% Five frozen Gamma-approach rays in reduced (Miller/hkl) coordinates:
% a*, b*, c*, [1 1 1], [2 1 -1] (prereg sec 5 item 3).
rays = [1 0 0; 0 1 0; 0 0 1; 1 1 1; 2 1 -1];
raylabels = {'a*', 'b*', 'c*', '[1 1 1]', '[2 1 -1]'};
end

function smags = frozen_smags()
smags = [1e-3, -1e-3, 1e-4, -1e-4];   % prereg sec 5 item 3
end

function [Q, rayOf, sOf] = ray_mag_grid(rays, smags)
% Q = [nrays*nsmags,3]; rayOf/sOf are parallel index/value arrays.
nr = size(rays, 1); ns = numel(smags);
Q = zeros(nr*ns, 3); rayOf = zeros(nr*ns, 1); sOf = zeros(nr*ns, 1);
row = 0;
for r = 1:nr
    for si = 1:ns
        row = row + 1;
        Q(row, :) = smags(si)*rays(r, :);
        rayOf(row) = r;
        sOf(row) = smags(si);
    end
end
end

function Dab = delta3(ntau)
% [3,3,ntau,ntau] real array, Dab(a,b,n,m) = delta_ab (independent of n,m).
Dab = reshape(eye(3), 3, 3, 1, 1) .* ones(1, 1, ntau, ntau);
end

function [qhat, qcart] = q_to_qhat(qrow, B)
% Canonical q-domain reduction (prereg sec 1) then Cartesian direction.
% qrow must be nonzero after reduction (qhat undefined at Gamma) -- never
% called at Gamma in this file.
K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*B;
qn = norm(qcart);
assert(qn > 0, 'q_to_qhat: qcart must be nonzero (qhat undefined at Gamma).');
qhat = qcart/qn;
end

function P0full = projector_replicated(qcart, Vc, alpha, ntau)
% Closed-form isolated reciprocal G=0 summand (prereg sec 5 item 3),
% P0_ab(q) = (4*pi/Vc)*qhat_a*qhat_b*exp(-|q|^2/(4*alpha^2)), replicated
% over every ordered sublattice pair (the G=0 plane wave carries no phase).
qn = norm(qcart);
qhat = qcart/qn;
P0 = (4*pi/Vc) * (qhat.' * qhat) * exp(-(qn^2)/(4*alpha^2));
P0full = complex(zeros(3, 3, ntau, ntau));
for n = 1:ntau
    for m = 1:ntau
        P0full(:, :, n, m) = P0;
    end
end
end

function dG0 = isolated_g0(geom, qrow, alpha, g_cut)
% Isolated G=[0,0,0] reciprocal candidate's contribution alone, from PUBLIC
% geom fields. Independent, trimmed, from-scratch reimplementation of the
% dG0 branch of test_invz_dipole_ewald_gammaC.m's reconstruct_parts (MATLAB
% local functions are file-scoped, so that helper is unreachable from here);
% this version omits the unused dR/dG-total/dS branches.
ntau = geom.ntau; B = geom.B; Vc = geom.Vc; taucart = geom.taucart; tau = geom.tau;
Gcart = geom.Gcart; Ghkl = geom.Ghkl;
K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*B;
k = Gcart + qcart; kk = sum(k.^2, 2); keep = (kk <= g_cut^2) & (kk > 0);
isG0 = all(Ghkl == 0, 2); g0k = isG0(keep);
dG0 = complex(zeros(3, 3, ntau, ntau));
if ~any(g0k), return; end
ksel = k(keep, :); kk2 = kk(keep); kernel = (4*pi/Vc)*exp(-kk2/(4*alpha^2))./kk2;
Gk = Gcart(keep, :);
ker0 = kernel(g0k); k0 = ksel(g0k, :);
for n = 1:ntau
    for m = 1:ntau
        d = taucart(m, :) - taucart(n, :);
        gph = exp(-1i*2*pi*(K*(tau(m, :) - tau(n, :)).'));
        ph0 = exp(1i*(Gk(g0k, :)*d.'));
        b0 = zeros(3, 3);
        for aa = 1:3
            for bb = 1:3
                b0(aa, bb) = sum(ker0.*k0(:, aa).*k0(:, bb).*ph0);
            end
        end
        dG0(:, :, n, m) = gph*b0;
    end
end
end

function Jbase = build_jbase_all(ion, C, dip_reg0, ex0)
% Jbase_ab = -gfac*dip_reg_ab(0) + Jex_ab(0) - gfac*(4*pi/(3*Vc))*delta_ab*ones(ntau),
% for ALL nine Cartesian (a,b) pairs and all sublattice pairs at once
% (task brief bullet 2), exactly mirroring invz_jq_modes.m's own Ewald-branch
% cc-only construction generalized to the full tensor. "Jex_ab(0)" includes
% sign(J12), matching the plan's "Gamma metadata" convention.
ntau = size(ion.tau, 1);
lorz = C.gfac*4*pi/(3*ion.Vc);
Dab = delta3(ntau);
Jbase = -C.gfac*dip_reg0 + sign(ion.J12)*ex0 - lorz*Dab;
end

% =====================================================================
% Gate-C4 (full Cartesian)
% =====================================================================
function test_gateC4_isolated_projector_with_gfac_all_nine_components(testCase)
ion = testCase.TestData.ion; C = testCase.TestData.C;
geom = testCase.TestData.geom0; eo = testCase.TestData.eo;
Q = testCase.TestData.Q; rayOf = testCase.TestData.rayOf; sOf = testCase.TestData.sOf;
raylabels = testCase.TestData.raylabels;
M = invz_ewald_metrics();
ntau = size(ion.tau, 1);

worst_margin = -inf;
for i = 1:size(Q, 1)
    qrow = Q(i, :);
    lbl = sprintf('ray %s, s=%.4g', raylabels{rayOf(i)}, sOf(i));

    dG0 = isolated_g0(geom, qrow, eo.alpha, eo.g_cut);
    [~, qcart] = q_to_qhat(qrow, geom.B);
    P0full = projector_replicated(qcart, geom.Vc, eo.alpha, ntau);

    Acoupling = -C.gfac*dG0;
    Bcoupling = -C.gfac*P0full;
    mres = M.mid(Acoupling, Bcoupling);
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C4: -gfac*(primitive isolated G=0 contribution) != -gfac*P0_ab(q) at M_id, ' ...
         'all nine components/sublattice pairs (%s, worst_margin=%.3e).'], lbl, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf(['Gate-C4 isolated-projector coupling-sign check: worst M_id margin over all 20 ' ...
    'ray/magnitude points (nine components x all sublattice pairs) = %.3e.\n'], worst_margin);
end

function test_gateC4_full_cartesian_base_reconstruction_algebra(testCase)
ion = testCase.TestData.ion; C = testCase.TestData.C;
geom = testCase.TestData.geom0;
Q = testCase.TestData.Q; rayOf = testCase.TestData.rayOf; sOf = testCase.TestData.sOf;
raylabels = testCase.TestData.raylabels;
dip_reg0 = testCase.TestData.dip0;
ex0 = testCase.TestData.ex0;
M = invz_ewald_metrics();
ntau = size(ion.tau, 1);

Jbase = build_jbase_all(ion, C, dip_reg0, ex0);
Dab = delta3(ntau);

worst_margin = -inf;
for i = 1:size(Q, 1)
    qrow = Q(i, :);
    lbl = sprintf('ray %s, s=%.4g', raylabels{rayOf(i)}, sOf(i));
    [qhat, ~] = q_to_qhat(qrow, geom.B);
    outer = reshape(qhat.'*qhat, 3, 3, 1, 1) .* ones(1, 1, ntau, ntau);

    % "Verify the complete limit reconstruction" (brief bullet 2): the
    % delta_ab/3 Lorentz term inside Jbase must cancel EXACTLY against the
    % delta_ab/3 term in the reconstruction formula, leaving the clean
    % closed form below. This is pure algebra on quantities evaluated ONLY
    % at exact Gamma (dip_reg0/ex0) plus a q-direction-only outer product
    % (qhat is s-independent) -- no finite-s primitive/projector value is
    % ever compared to the no-exp formula here (E1-safe by construction).
    % NB (information content): this reconstruction is an algebraic IDENTITY in
    % qhat -- it holds for every ray/s -- so the 20-point loop confirms the
    % reconstruction formula is implemented without a typo, but is NOT 20
    % independent finite-q validations. The independent finite-q/finite-s content
    % lives in the isolated-P0-projector tests, which E1 keeps separate.
    Jrecon         = Jbase + C.gfac*(4*pi/ion.Vc)*(Dab/3 - outer);
    JreconExpected = -C.gfac*dip_reg0 + sign(ion.J12)*ex0 - C.gfac*(4*pi/ion.Vc)*outer;

    mres = M.mid(Jrecon, JreconExpected);
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C4: Jbase_ab + gfac*(4pi/Vc)*(delta_ab/3-qhat_a*qhat_b) does not reduce to the ' ...
         'clean closed form -gfac*dip_reg_ab(0)+Jex_ab(0)-gfac*(4pi/Vc)*qhat_a*qhat_b at M_id ' ...
         '(delta_ab/3 cancellation), all nine components/sublattice pairs (%s, worst_margin=%.3e).'], ...
        lbl, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf(['Gate-C4 full-Cartesian base-reconstruction algebra: worst M_id margin over all 20 ' ...
    'ray/magnitude points (nine components x all sublattice pairs) = %.3e.\n'], worst_margin);
end

function test_gateC4_cc_matches_jq_modes_and_qpath_formula(testCase)
ion = testCase.TestData.ion; C = testCase.TestData.C;
geom = testCase.TestData.geom0; eo = testCase.TestData.eo;
dip_reg0 = testCase.TestData.dip0; ex0 = testCase.TestData.ex0;
Q = testCase.TestData.Q; rayOf = testCase.TestData.rayOf; sOf = testCase.TestData.sOf;
raylabels = testCase.TestData.raylabels;
M = invz_ewald_metrics();
ntau = size(ion.tau, 1);
lorz = C.gfac*4*pi/(3*ion.Vc);

Jbase = build_jbase_all(ion, C, dip_reg0, ex0);
Jbase_cc_hand = squeeze(Jbase(3, 3, :, :));

% Live invz_jq_modes call, exact Gamma only -- no finite-s value involved.
[~, info] = invz_jq_modes(ion, [0 0 0], struct('dipole', 'ewald', 'ewald', eo, 'cache', false));

mres0 = M.mid(Jbase_cc_hand, info.Jpath_base_cc);
verifyTrue(testCase, mres0.pass, sprintf(...
    ['Gate-C4: hand-built Jbase_cc (from the primitive directly) != info.Jpath_base_cc at M_id ' ...
     '(worst_margin=%.3e).'], mres0.worst_margin));

Jgamma_cc_hand = Jbase_cc_hand + lorz*ones(ntau);
mres1 = M.mid(Jgamma_cc_hand, info.Jgamma_cc);
verifyTrue(testCase, mres1.pass, sprintf(...
    'Gate-C4: hand-built Jgamma_cc != info.Jgamma_cc at M_id (worst_margin=%.3e).', mres1.worst_margin));

Dab = delta3(ntau);
% The genuine cross-checks are mres0/mres1 above: hand-built Jbase_cc/Jgamma_cc
% (from the primitive directly) vs the live info.Jpath_base_cc/info.Jgamma_cc.
% Given mres0 (info.Jpath_base_cc == Jbase_cc), the loop below is an algebraic
% IDENTITY in qhat -- it re-confirms the cc slice of the nine-component
% reconstruction equals the production q-path formula for every qhat (a
% formula-typo catcher), NOT 20 independent finite-q validations.
worst_margin = -inf;
for i = 1:size(Q, 1)
    qrow = Q(i, :);
    lbl = sprintf('ray %s, s=%.4g', raylabels{rayOf(i)}, sOf(i));
    [qhat, ~] = q_to_qhat(qrow, geom.B);
    kz2 = qhat(3)^2;

    % Production q-path formula, built purely from info.Jpath_base_cc (a
    % live invz_jq_modes output) -- never from a finite-s dip(q) value.
    JreconViaInfo = info.Jpath_base_cc + C.gfac*(4*pi/ion.Vc)*(1/3 - kz2)*ones(ntau);

    % The SAME point's cc slice from the full nine-component reconstruction
    % (test_gateC4_full_cartesian_base_reconstruction_algebra's formula),
    % rebuilt independently here for this test's self-containment.
    outer = reshape(qhat.'*qhat, 3, 3, 1, 1) .* ones(1, 1, ntau, ntau);
    Jrecon_full = Jbase + C.gfac*(4*pi/ion.Vc)*(Dab/3 - outer);
    JreconFull_cc = squeeze(Jrecon_full(3, 3, :, :));

    mres = M.mid(JreconViaInfo, JreconFull_cc);
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C4: cc slice of the full nine-component reconstruction != production q-path formula ' ...
         'info.Jpath_base_cc + gfac*(4pi/Vc)*(1/3-qhat_c^2)*ones(ntau) at M_id (%s, worst_margin=%.3e).'], ...
        lbl, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf(['Gate-C4 cc cross-check: Jbase_cc/Jgamma_cc vs info M_id margins = %.3e / %.3e; worst ' ...
    'reconstructed-limit margin over 20 ray/magnitude points = %.3e.\n'], ...
    mres0.worst_margin, mres1.worst_margin, worst_margin);
end

% =====================================================================
% Gate-C6 (demag)
% =====================================================================
function q = c6_qvec()
% Small deterministic, non-Gamma q array. Not a physically-meaningful BZ
% quadrature: Gate-C6 tests SHAPE-INVARIANCE identities, not a converged
% Sigma_c/Tc value, so grid quality/size is irrelevant (mirrors
% test_invz_demag_invariance.m's own reasoning for its small_grid()).
q = [ 0.13  0.27 -0.09
     -0.31  0.05  0.22
      0.44 -0.18  0.11
      0.02  0.36 -0.27
     -0.19 -0.24  0.38];
end

function ion = ion_shape(kind)
ion = invz_ion();
switch kind
    case 'off'
        ion.demag = 0;
    case 'sphere'
        ion.demag = 1; ion.alpha = 1;
    case 'needle'
        ion.demag = 1; ion.alpha = 0;
    otherwise
        error('ion_shape: unknown kind ''%s''.', kind);
end
end

function opts = bf_c6_opts()
opts = struct('dpRng', 8, 'cache', false);
end

function opts = ew_c6_opts(eo)
opts = struct('dipole', 'ewald', 'ewald', eo, 'cache', false);
end

function test_gateC6_bruteforce_three_shapes_identical_core_quantities(testCase)
c6_three_shapes_identical(testCase, bf_c6_opts(), 'bruteforce');
end

function test_gateC6_ewald_three_shapes_identical_core_quantities(testCase)
eo = testCase.TestData.eo;
c6_three_shapes_identical(testCase, ew_c6_opts(eo), 'ewald');
end

function c6_three_shapes_identical(testCase, opts, backendLabel)
% Within ONE backend, demag/alpha (off/sphere/needle) must never move Jnu,
% Jcc0, Sigma_c, or the zero-field Tc -- traced through invz_jq_modes.m,
% dm_cc/dm_aa feed ONLY Jshape_cc/Jaa0, never the per-q Jnu loop or
% Jcc0_dipole/Jcc0. verifyEqual with no tolerance mirrors
% test_invz_demag_invariance.m's test_tc0_and_sigma_crit_demag_invariant
% precedent exactly (that test already passes for bruteforce+sphere; this
% extends the SAME claim to both backends and the needle shape).
q = c6_qvec();
kinds = {'off', 'sphere', 'needle'};
Jnu_all = cell(1, 3); Jcc0_all = zeros(1, 3); Sc_all = zeros(1, 3); Tc_all = zeros(1, 3);
for k = 1:3
    ion = ion_shape(kinds{k});
    [Jnu, info] = invz_jq_modes(ion, q, opts);
    Jnu_all{k} = Jnu;
    Jcc0_all(k) = info.Jcc0;
    Sc_all(k) = invz_sigma_crit(info.Jcc0, Jnu(:));
    Tc_all(k) = invz_critical_T0field(ion, Sc_all(k), info.Jcc0);
end
for k = 2:3
    verifyEqual(testCase, Jnu_all{k}, Jnu_all{1}, sprintf(...
        'Gate-C6 (%s): Jnu differs for shape ''%s'' vs ''off''.', backendLabel, kinds{k}));
    verifyEqual(testCase, Jcc0_all(k), Jcc0_all(1), sprintf(...
        'Gate-C6 (%s): Jcc0 differs for shape ''%s'' vs ''off'' (R2007 demag-invariant criticality coupling).', ...
        backendLabel, kinds{k}));
    verifyEqual(testCase, Sc_all(k), Sc_all(1), sprintf(...
        'Gate-C6 (%s): Sigma_c differs for shape ''%s'' vs ''off''.', backendLabel, kinds{k}));
    verifyEqual(testCase, Tc_all(k), Tc_all(1), sprintf(...
        'Gate-C6 (%s): invz_critical_T0field(...) differs for shape ''%s'' vs ''off'' (Tc(B=0) must be demag-invariant).', ...
        backendLabel, kinds{k}));
end
fprintf('Gate-C6 (%s): Jnu/Jcc0/Sigma_c/Tc(B=0) identical across off/sphere/needle -- Sc=%.6g, Tc=%.6g K.\n', ...
    backendLabel, Sc_all(1), Tc_all(1));
end

function test_gateC6_bruteforce_jshape_and_jaa0_match_ellipsoid_demagn(testCase)
c6_ellipsoid_demagn_check(testCase, bf_c6_opts(), 'bruteforce');
end

function test_gateC6_ewald_jshape_and_jaa0_match_ellipsoid_demagn(testCase)
eo = testCase.TestData.eo;
c6_ellipsoid_demagn_check(testCase, ew_c6_opts(eo), 'ewald');
end

function c6_ellipsoid_demagn_check(testCase, opts, backendLabel)
C = testCase.TestData.C;
M = invz_ewald_metrics();
q = c6_qvec();

ionOff = ion_shape('off');
[~, infoOff] = invz_jq_modes(ionOff, q, opts);
verifyEqual(testCase, infoOff.Jshape_cc, 0, sprintf(...
    'Gate-C6 (%s): off-shape (demag=0) Jshape_cc must be exactly 0.', backendLabel));

shapes = {'sphere', 'needle'};
alphas = [1 0];
for k = 1:2
    ionS = ion_shape(shapes{k});
    [~, infoS] = invz_jq_modes(ionS, q, opts);
    Nd = ellipsoid_demagn(alphas(k));

    expected_Jshape_cc = 4*C.gfac*(4*pi/ionS.Vc)*ionS.demag*Nd(3, 3);
    mres = M.mid(infoS.Jshape_cc, expected_Jshape_cc);
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C6 (%s): Jshape_cc(%s) != 4*gfac*(4*pi/Vc)*demag*Nd(3,3) [ellipsoid_demagn] at M_id ' ...
         '(worst_margin=%.3e).'], backendLabel, shapes{k}, mres.worst_margin));

    expected_Jaa0 = infoOff.Jaa0 - 4*C.gfac*(4*pi/ionS.Vc)*ionS.demag*Nd(1, 1);
    mres2 = M.mid(infoS.Jaa0, expected_Jaa0);
    verifyTrue(testCase, mres2.pass, sprintf(...
        ['Gate-C6 (%s): Jaa0(%s) != Jaa0(off) - 4*gfac*(4*pi/Vc)*demag*Nd(1,1) [ellipsoid_demagn] at ' ...
         'M_id (worst_margin=%.3e).'], backendLabel, shapes{k}, mres2.worst_margin));

    fprintf(['Gate-C6 (%s) ellipsoid_demagn cross-check (%s): Jshape_cc M_id margin=%.3e, Jaa0 M_id ' ...
        'margin=%.3e.\n'], backendLabel, shapes{k}, mres.worst_margin, mres2.worst_margin);
end
end

function test_gateC6_jshape_cc_backend_agnostic_by_construction(testCase)
% info.Jshape_cc = 4*gfac*(4*pi/Vc)*demag*Nd(3,3) is computed ONCE in
% invz_jq_modes.m's shared pre-dispatch code, strictly BEFORE the
% bruteforce/ewald branch split, and never reads dip0/dip or any
% backend-specific quantity -- an EXACT algebraic identity by construction.
% This does not invent a new brute/Ewald absolute-agreement tolerance for
% Gate C6 (the task's explicit constraint): it pins an already-shared,
% backend-independent code path, not a physics comparison between two
% independently-computed dipolar sums.
eo = testCase.TestData.eo;
q = c6_qvec();

ionS = ion_shape('sphere');
[~, infoBFs] = invz_jq_modes(ionS, q, bf_c6_opts());
[~, infoEWs] = invz_jq_modes(ionS, q, ew_c6_opts(eo));
verifyEqual(testCase, infoEWs.Jshape_cc, infoBFs.Jshape_cc, ...
    'Gate-C6: Jshape_cc must be backend-independent by construction (sphere).');

ionN = ion_shape('needle');
[~, infoBFn] = invz_jq_modes(ionN, q, bf_c6_opts());
[~, infoEWn] = invz_jq_modes(ionN, q, ew_c6_opts(eo));
verifyEqual(testCase, infoEWn.Jshape_cc, infoBFn.Jshape_cc, ...
    'Gate-C6: Jshape_cc must be backend-independent by construction (needle).');
end

function test_gateC6_ewald_eopts_has_no_demag_or_surface_control(testCase)
ion = testCase.TestData.ion;
eo = testCase.TestData.eo;
verifyEqual(testCase, sort(fieldnames(eo)), sort({'alpha'; 'r_cut'; 'g_cut'; 'boundary'}));

% The primitive's own opts validation actively REJECTS an attempted
% demag/surface control field, at both the primitive and the jq_modes layer.
eoBad = eo; eoBad.demag = 1;
verifyError(testCase, @() invz_dipole_ewald([0 0 0], ion.a, ion.tau, eoBad), 'invz:ewaldArgs');

eoBad2 = eo; eoBad2.surface = true;
verifyError(testCase, @() invz_dipole_ewald([0 0 0], ion.a, ion.tau, eoBad2), 'invz:ewaldArgs');

verifyError(testCase, ...
    @() invz_jq_modes(ion, [0 0 0], struct('dipole', 'ewald', 'ewald', eoBad, 'cache', false)), ...
    'invz:jqModesEwaldOptsFields');
end
