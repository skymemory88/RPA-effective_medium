function tests = test_invz_dipole_ewald_ref
% Gate-B (Step-4 Task 8): retained acceptance test for the INDEPENDENT
% scalar-Coulomb Ewald oracle invz_scalar_ewald_ref.m, docs/invzp_ewald_prereg.md
% sec 4 (FROZEN). The oracle is a genuinely separate code path (a differentiated
% scalar lattice potential, phase on R only) from the tensor primitive
% invz_dipole_ewald.m; this file compares the two at M_FD -- the deliberately
% looser, dimensioned (Angstrom^-3) FD-limited tolerance -- and NEVER reports
% an M_FD pass as an M_T pass. The primitive's own tighter alpha/cutoff
% self-convergence is independently covered by Gate-A #9
% (test_gateA9_alpha_bracket_independence / test_gateA9_cutoff_ladders_separate_axes
% / test_gateA9_default_vs_joint_refinement in test_invz_dipole_ewald.m) and is
% neither re-derived nor re-labelled here.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));            % invz_projected
addpath(fullfile(here, '..', '..'));      % repo root: invz_dipole_ewald, MF_dipole
addpath(fullfile(here, '..', '..', 'invz_common'));  % invz_ion, invz_const
end

% =====================================================================
% shared helpers (names deliberately free of "test" so functiontests skips them)
% =====================================================================
function eo = mk_eopts(alpha, r_cut, g_cut, boundary)
eo = struct('alpha', alpha, 'r_cut', r_cut, 'g_cut', g_cut, 'boundary', boundary);
end

function eo = eopts_alpha_matched(alpha)
% alpha-matched generous cutoffs (prereg sec 2 / sec 4 Gate-B): r_cut=6.5/alpha,
% g_cut=13*alpha -- mirrors test_invz_dipole_ewald.m's own gateA9 helper of the
% same name (that file's local functions are not reachable from here, so this
% is an independent one-line reimplementation of the same frozen rule).
eo = mk_eopts(alpha, 6.5/alpha, 13*alpha, 'conducting_k0_omitted');
end

% =====================================================================
% Task 8 -- Gate-B: oracle struct schema (cheap, non-comparative red/green
% anchor -- this is the test that fails first, with invz_scalar_ewald_ref
% absent, per the TDD brief)
% =====================================================================
function test_gateB_oracle_struct_schema(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
fx = invz_ewald_fixtures();
q = fx.q_int(1,:);
alpha = fx.alpha0;
eo = eopts_alpha_matched(alpha);

oracle = invz_scalar_ewald_ref(q, a, tau, eo.alpha, eo.r_cut, eo.g_cut);

% ---- exact field set (brief: "Match the exact field set the plan's Task 8
% lists") ---------------------------------------------------------------
want = {'adjacent_residual','h','raw','richardson_coarse','richardson_fine','self_analytic'};
have = fieldnames(oracle);
verifyTrue(testCase, all(ismember(want, have)), ...
    'Gate-B: invz_scalar_ewald_ref struct is missing one or more required fields.');

% ---- shapes --------------------------------------------------------------
verifyEqual(testCase, oracle.h, [4e-3 2e-3 1e-3]);
verifyEqual(testCase, size(oracle.raw),               [3 3 ntau ntau 3]);
verifyEqual(testCase, size(oracle.richardson_coarse),  [3 3 ntau ntau]);
verifyEqual(testCase, size(oracle.richardson_fine),    [3 3 ntau ntau]);
verifyEqual(testCase, size(oracle.self_analytic),      [3 3 ntau ntau]);
verifyEqual(testCase, size(oracle.adjacent_residual),  [3 3 ntau ntau]);

% ---- finiteness of every exposed raw/extrapolated/residual array ---------
verifyTrue(testCase, all(isfinite(oracle.h(:))), 'oracle.h not finite.');
verifyTrue(testCase, all(isfinite(oracle.raw(:))), 'oracle.raw not finite.');
verifyTrue(testCase, all(isfinite(oracle.richardson_coarse(:))), 'oracle.richardson_coarse not finite.');
verifyTrue(testCase, all(isfinite(oracle.richardson_fine(:))), 'oracle.richardson_fine not finite.');
verifyTrue(testCase, all(isfinite(oracle.self_analytic(:))), 'oracle.self_analytic not finite.');
verifyTrue(testCase, all(isfinite(oracle.adjacent_residual(:))), 'oracle.adjacent_residual not finite.');

% ---- adjacent_residual is exactly richardson_fine - richardson_coarse ----
verifyEqual(testCase, oracle.adjacent_residual, ...
    oracle.richardson_fine - oracle.richardson_coarse, 'AbsTol', 0);

% ---- self_analytic is EXACTLY -delta_nm delta_ab 4*alpha^3/(3*sqrt(pi)) --
% (closed form, no FD involved -- exact-identity tolerance, not M_FD).
M = invz_ewald_metrics();
selfval = 4*eo.alpha^3/(3*sqrt(pi));
expected_self = complex(zeros(3,3,ntau,ntau));
for n = 1:ntau
    expected_self(:,:,n,n) = -selfval*eye(3);
end
mres = M.mid(oracle.self_analytic, expected_self);
verifyTrue(testCase, mres.pass, sprintf( ...
    'Gate-B: oracle.self_analytic is not exactly -4*alpha^3/(3*sqrt(pi))*eye(3) on-diagonal / 0 off-diagonal (M_id worst_margin=%.3e).', ...
    mres.worst_margin));
end

% =====================================================================
% Task 8 -- Gate-B: retained physics acceptance (prereg sec 4, FROZEN)
% =====================================================================
function test_gateB_scalar_oracle_agreement(testCase)
% For all three fx.q_int, all five frozen alpha multipliers (with the
% alpha-matched generous cutoffs r_cut=6.5/alpha, g_cut=13*alpha), every
% ordered sublattice pair, and all nine Cartesian components:
%   (1) the two Richardson estimates (R12 vs R23) agree at M_FD;
%   (2) off-diagonal (n~=m) blocks of the finer Richardson estimate (R23)
%       agree with the PRIMITIVE at M_FD;
%   (3) on n=m, primitive minus the regularized non-self oracle (R23) equals
%       the analytic self tensor at M_FD.
% M_FD's tolerance (AbsTol_FD=2e-8 Angstrom^-3 + RelTol_FD=1e-7*max(|A|,|B|))
% is applied COMPONENTWISE with no shared cross-component scale (unlike
% M_T's single T_scale), so restricting a comparison to a masked subset of
% components (e.g. only the off-diagonal or only the same-site blocks) is
% exactly equivalent to padding the excluded components with matching zeros
% and comparing the complete tensor -- both give the identical worst_margin.
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau,1);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();

verifyEqual(testCase, size(fx.q_int), [3 3]);
verifyGreaterThan(testCase, fx.alpha0, 0);

mult = [0.6 0.8 1.0 1.2 1.5];

sameMask = false(3,3,ntau,ntau);
for n = 1:ntau
    sameMask(:,:,n,n) = true;
end
offMask = ~sameMask;

worst_rich_margin = -inf; worst_rich_ratio = 0;
worst_off_margin  = -inf; worst_off_ratio  = 0;
worst_self_margin = -inf; worst_self_ratio = 0;

for qi = 1:size(fx.q_int,1)
    q = fx.q_int(qi,:);
    for ai = 1:numel(mult)
        alpha = mult(ai)*fx.alpha0;
        eo = eopts_alpha_matched(alpha);

        [dip, ~, ~] = invz_dipole_ewald(q, a, tau, eo);
        oracle = invz_scalar_ewald_ref(q, a, tau, eo.alpha, eo.r_cut, eo.g_cut);

        % ---- every exposed raw/extrapolated/residual array is finite -----
        verifyTrue(testCase, all(isfinite(oracle.h(:))) && all(isfinite(oracle.raw(:))) ...
            && all(isfinite(oracle.richardson_coarse(:))) && all(isfinite(oracle.richardson_fine(:))) ...
            && all(isfinite(oracle.self_analytic(:))) && all(isfinite(oracle.adjacent_residual(:))), ...
            sprintf('Gate-B: a non-finite value in the oracle struct at q_int(%d,:), alpha=%.6g (x%.2g).', ...
            qi, alpha, mult(ai)));

        % ---- (1) R12 vs R23, every ordered pair, all nine components -----
        mrich = M.mfd(oracle.richardson_coarse, oracle.richardson_fine);
        verifyTrue(testCase, mrich.pass, sprintf( ...
            'Gate-B: R12/R23 Richardson estimates disagree at M_FD, q_int(%d,:) alpha=%.6g (x%.2g), worst_margin=%.3e.', ...
            qi, alpha, mult(ai), mrich.worst_margin));
        worst_rich_margin = max(worst_rich_margin, mrich.worst_margin);
        worst_rich_ratio  = max(worst_rich_ratio,  mrich.worst_ratio);

        % ---- (2) off-diagonal (n~=m): finer Richardson vs primitive -------
        moff = M.mfd(dip(offMask), oracle.richardson_fine(offMask));
        verifyTrue(testCase, moff.pass, sprintf( ...
            'Gate-B: off-diagonal R23 vs primitive fails M_FD, q_int(%d,:) alpha=%.6g (x%.2g), worst_margin=%.3e.', ...
            qi, alpha, mult(ai), moff.worst_margin));
        worst_off_margin = max(worst_off_margin, moff.worst_margin);
        worst_off_ratio  = max(worst_off_ratio,  moff.worst_ratio);

        % ---- (3) n=m: primitive - regularized non-self oracle == self_analytic
        actual_self_diag   = dip(sameMask) - oracle.richardson_fine(sameMask);
        expected_self_diag = oracle.self_analytic(sameMask);
        mself = M.mfd(actual_self_diag, expected_self_diag);
        verifyTrue(testCase, mself.pass, sprintf( ...
            'Gate-B: self-term residual (primitive - regularized non-self oracle) fails M_FD vs the analytic self tensor, q_int(%d,:) alpha=%.6g (x%.2g), worst_margin=%.3e.', ...
            qi, alpha, mult(ai), mself.worst_margin));
        worst_self_margin = max(worst_self_margin, mself.worst_margin);
        worst_self_ratio  = max(worst_self_ratio,  mself.worst_ratio);

        fprintf(['Gate-B (q_int(%d,:), alpha=%.6g x%.2g): R12/R23 M_FD margin=%.3e, ' ...
            'off-diag-vs-primitive M_FD margin=%.3e, self-residual M_FD margin=%.3e.\n'], ...
            qi, alpha, mult(ai), mrich.worst_margin, moff.worst_margin, mself.worst_margin);
    end
end

fprintf(['Gate-B scalar-oracle agreement (pre-freeze calibration worst off-diag 9.2e-9, ' ...
    'self 4.1e-9 Angstrom^-3 -- reproduced as a PASS, not asserted bit-exact): worst R12/R23 ' ...
    'M_FD margin=%.3e (ratio %.3e), worst off-diagonal M_FD margin=%.3e (ratio %.3e), worst ' ...
    'self-residual M_FD margin=%.3e (ratio %.3e).\n'], ...
    worst_rich_margin, worst_rich_ratio, worst_off_margin, worst_off_ratio, ...
    worst_self_margin, worst_self_ratio);

% Gate-B is M_FD-only by definition (prereg sec 4). The primitive's own
% tighter alpha/cutoff self-convergence at M_T is independently retained and
% asserted by Gate-A #9 in test_invz_dipole_ewald.m; this file makes no M_T
% claim anywhere, so an M_FD pass here is never reported or usable as an M_T
% pass.
end
