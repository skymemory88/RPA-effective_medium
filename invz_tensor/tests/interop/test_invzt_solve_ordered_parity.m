function tests = test_invzt_solve_ordered_parity
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', '..'));
addpath(fullfile(here, '..', '..', '..', 'invz_common'));
addpath(fullfile(here, '..', '..', '..', 'invz_projected'));
end

function test_ordered_no_odd_parity_live(testCase)
% Ordered-phase analog of the PM no-ODD live parity: both legs on the SAME grid
% and couplings (J0eff = info.Jcc0, Jxx0 = info.Jaa0 -> identical single-ion MF
% fixed point, hence identical m0 by construction). Both legs now renormalize the
% WHOLE local cc susceptibility in one piece (2026-07-20 amendment -- the tensor
% leg has no chi_rest knob anymore), so the tensor medium differs from the scalar
% one only by the remaining named residual: cross-Cartesian chi0 elements. Gate
% Sigma0/alpha_m at the PM frozen-baseline scale and m0 tight.
ion = invz_ion();  T = 0.1;  Bx = 3.0;
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
[Jnu, info] = invz_jq_modes(ion, g.qvec, struct('dpRng', 15, 'cache', true));
pj = invz_solve_point_ordered(ion, T, [Bx 0 0], Jnu(:), ...
    struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0));
pt = invzt_solve_point_ordered(ion, T, [Bx 0 0], lat, struct('odd', false));
verifyTrue(testCase, pj.is_ordered && pj.converged);
verifyTrue(testCase, pt.is_ordered && pt.converged);
verifyLessThan(testCase, abs(pt.m0 - pj.m0), 1e-6);
verifyLessThan(testCase, abs(pt.Sigma0 - pj.Sigma0), 5e-3);
verifyLessThan(testCase, abs(pt.alpha_m - pj.alpha_m), 5e-3);
fprintf(['INTEROP ordered no-ODD live: dm0 = %.3e, dSigma0 = %.3e, dalpha_m = %.3e, ' ...
    'dcrit = %.3e (m0 %.4f; Sigma0 tensor %.5f / proj %.5f)\n'], ...
    pt.m0 - pj.m0, pt.Sigma0 - pj.Sigma0, pt.alpha_m - pj.alpha_m, ...
    pt.crit - pj.crit, pt.m0, pt.Sigma0, pj.Sigma0);
end
