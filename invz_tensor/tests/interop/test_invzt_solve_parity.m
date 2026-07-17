function tests = test_invzt_solve_parity
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));                                 % invz_tensor
addpath(fullfile(here, '..', '..', '..'));                           % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', '..', 'invz_common'));            % shared single-ion engine
addpath(fullfile(here, '..', '..', '..', 'invz_projected'));         % parity targets
addpath(fullfile(here, '..', '..', '..', 'invz_projected', 'tests'));% invz_odd_anchors fixture
end

function test_no_odd_parity_live(testCase)
% Live version of the frozen-baseline comparison: no-ODD both sides, SAME grid
% convention (legacy_inclusive), and 'Jxx0' = info.Jaa0 passed to the projected
% leg so the single-ion transverse mean field is apples-to-apples. REPORT
% dSigma0/dcrit; assert finiteness + convergence.
ion = invz_ion();  T = 1.6;  Bx = 0.5;
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
[Jnu, info] = invz_jq_modes(ion, g.qvec, struct('dpRng', 15, 'cache', true));
pj = invz_solve_point(ion, T, [Bx 0 0], Jnu(:), ...
    struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0));            % projected leg (Jxx0 = info.Jaa0)
pt = invzt_solve_point(ion, T, [Bx 0 0], lat, ...
    struct('odd', false, 'chi_rest', false));                 % tensor leg
verifyTrue(testCase, pt.converged && pj.converged);
verifyTrue(testCase, isfinite(pt.Sigma0) && isfinite(pt.crit));
fprintf(['INTEROP no-ODD live: dSigma0 = %.3e, dcrit = %.3e ' ...
    '(Sigma0 tensor %.5f / proj %.5f; crit tensor %.5f / proj %.5f)\n'], ...
    pt.Sigma0 - pj.Sigma0, pt.crit - pj.crit, pt.Sigma0, pj.Sigma0, pt.crit, pj.crit);
end

function test_odd_on_parity_live(testCase)
% Odd-on comparison vs the projected Tier-1 solver: geometric blocks via
% invz_odd_blocks on the SAME legacy grid, projected solve rebuilds cc modes
% from Vcc + deltaJ (E1/E4/E5). REPORT dSigma0/dcrit; assert finiteness (the
% tensor A1 12x12 RPA carries the transverse mediation, the projected leg uses
% the scalar-Sigma deltaJ closure -- a physical, not exact, comparison).
ion = invz_ion();  T = 1.6;  Bx = 0.5;
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
[Vca, Vcb, Vcc, infoS] = invz_odd_blocks(ion, g.qvec, struct('dpRng', 15, 'cache', true));
ob = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoS.Jcc0);
pj = invz_solve_point(ion, T, [Bx 0 0], [], ...
    struct('odd', true, 'odd_blocks', ob, 'J0eff', infoS.Jcc0, 'Jxx0', infoS.Jaa0));
pt = invzt_solve_point(ion, T, [Bx 0 0], lat, struct('odd', true));
verifyTrue(testCase, pt.converged && pj.converged);
verifyTrue(testCase, isfinite(pt.Sigma0) && isfinite(pt.crit) ...
    && isfinite(pj.Sigma0) && isfinite(pj.crit));
fprintf(['INTEROP odd-on live: dSigma0 = %.3e, dcrit = %.3e ' ...
    '(Sigma0 tensor %.5f / proj %.5f; crit tensor %.5f / proj %.5f)\n'], ...
    pt.Sigma0 - pj.Sigma0, pt.crit - pj.crit, pt.Sigma0, pj.Sigma0, pt.crit, pj.crit);
end
