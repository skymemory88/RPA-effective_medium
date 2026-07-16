function tests = test_invz_demag_invariance
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function q = small_grid()
% 6^3 grid, Gamma excluded: the invariances below are exact, so grid quality is irrelevant.
cmd = "qVec_generator(invz_ion().a, 'mode', 'grid', 'grid', [6 6 6], 'range', [-0.5 0.5])";
[~, q] = evalc(cmd);
q = q(any(abs(q) > 1e-12, 2), :);
end

function test_tc0_and_sigma_crit_demag_invariant(testCase)
% R 2007: the demagnetizing field cancels from the critical condition. Sigma_c and the
% zero-field Tc read ONLY Jnu/Jcc0, so with the shape knob on (sphere) they must be
% BIT-IDENTICAL to the intrinsic values. (At B = 0 the transverse channel is inert:
% <Jx> = 0, so the demag-aware info.Jaa0 cannot enter either.)
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
q = small_grid();
[J0nu, i0] = invz_jq_modes(ion0, q, struct('dpRng', 5, 'cache', false));
[JSnu, iS] = invz_jq_modes(ionS, q, struct('dpRng', 5, 'cache', false));
verifyEqual(testCase, JSnu, J0nu);
verifyEqual(testCase, iS.Jcc0, i0.Jcc0);
Sc0 = invz_sigma_crit(i0.Jcc0, J0nu(:));
ScS = invz_sigma_crit(iS.Jcc0, JSnu(:));
verifyEqual(testCase, ScS, Sc0);
verifyEqual(testCase, invz_critical_T0field(ionS, ScS, iS.Jcc0), ...
                      invz_critical_T0field(ion0, Sc0, i0.Jcc0));
end

function test_crit_demag_enters_only_via_transverse_channel(testCase)
% With the transverse coupling PINNED to the intrinsic value (opts.Jxx0), the full 1/z
% criticality at finite field must be demag-invariant -- the ordering channel carries no
% shape dependence. Unpinned (Jxx0 = demag-aware info.Jaa0), crit legitimately moves:
% that is the internal-vs-applied transverse field correction (canonical bullet 3). SLOW.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion0 = invz_ion();
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
q = small_grid();
[J0nu, i0] = invz_jq_modes(ion0, q, struct('dpRng', 5, 'cache', false));
[JSnu, iS] = invz_jq_modes(ionS, q, struct('dpRng', 5, 'cache', false));
T = 0.31;  B = 5.0;
p0 = invz_solve_point(ion0, T, B, J0nu(:), struct('J0eff', i0.Jcc0, 'Jxx0', ion0.Jxx0));
pS = invz_solve_point(ionS, T, B, JSnu(:), struct('J0eff', iS.Jcc0, 'Jxx0', ion0.Jxx0));
verifyEqual(testCase, pS.crit, p0.crit);          % pinned: bit-identical
pU = invz_solve_point(ionS, T, B, JSnu(:), struct('J0eff', iS.Jcc0, 'Jxx0', iS.Jaa0));
verifyGreaterThan(testCase, abs(pU.crit - p0.crit), 1e-9);   % unpinned: physical shift
end
