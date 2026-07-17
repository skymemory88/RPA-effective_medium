function tests = test_invzt_critical_parity
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

function test_zero_field_closed_form_parity(testCase)
% Tensor lat slices through the projected closed-form chain vs invz_odd_zero_field
% — BOTH on the SAME legacy_inclusive 8^3/dpRng-15 grid (v2: grid contract).
% READ invz_odd_zero_field.m first; adapt out-field names if they differ.
ion = invz_ion();
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
Vca = lat.Jt(3:3:12, 1:3:12, :);  Vcb = lat.Jt(3:3:12, 2:3:12, :);  Vcc = lat.Jt(3:3:12, 3:3:12, :);
[TcP, outP] = invz_odd_zero_field(ion, struct('mode', 'full', 'grids', {{8}}, 'dpRng', 15));
Xp = invz_chiperp(ion, TcP(1), [0 0 0], struct());
[dJ, d] = invz_odd_deltaJ(Vca, Vcb, Xp);
Jnu = invz_odd_modes(Vcc, dJ);
verifyEqual(testCase, d, outP.d_at_Tc(1), 'RelTol', 1e-9);
Sc = invz_sigma_crit(lat.info.Jcc0 - d, Jnu(:));
verifyEqual(testCase, Sc, outP.Sc_at_Tc(1), 'RelTol', 1e-9);
end

function test_bc_parity_no_odd(testCase)
% No-ODD Bc parity: tensor A1 vs the projected scalar solver, SAME legacy 8^3 grid.
% DEVIATION (task-7 report Deviation 1, flagged -- not a weakened assertion; the
% AbsTol stays 0.05 and PASSES at ~0.002 T once the comparison is apples-to-apples):
% the parity leg uses chi_rest = TRUE, NOT the brief-sketch's chi_rest = false. The
% projected invz_solve_point keeps the FULL local cc chi0 (invz_chi0z(si,T,1i*wn)
% (3,3), NO dominant/rest split -- invz_projected/invz_solve_point.m line 229), so
% the apples-to-apples tensor leg must ALSO keep crest. chi_rest = false drops the
% excited-state cc chi0: ~4e-5 negligible at the low-field T6 baseline (0.5 T,
% test_no_odd_parity_live) but growing with field to a 0.16 T Bc error at this
% HIGH-field Bc(1.2 K) ~ 2.7 T (probe: chi_rest=true -> tensor 2.789 T vs projected
% 2.787 T; chi_rest=false -> 2.627 T). The chi_rest=false Bc is reported below as
% the "split" diagnostic (brief Step 5), not asserted.
ion = invz_ion();  T = 1.2;
g = invzt_qgrid(8, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 15, 'cache', true));
[Jnu, info] = invz_jq_modes(ion, g.qvec, struct('dpRng', 15, 'cache', true));
Bc_t = invzt_critical(ion, T, lat, [2 6], struct('odd', false, 'chi_rest', true));
Bc_s = invz_critical(ion, T, Jnu(:), struct('J0eff', info.Jcc0, 'Jxx0', info.Jaa0));
Bc_t_split = invzt_critical(ion, T, lat, [2 6], struct('odd', false, 'chi_rest', false));
fprintf('Bc(1.2 K) no-ODD: tensor %.4f T, projected %.4f T (chi_rest=false split diag: %.4f T)\n', ...
    Bc_t, Bc_s, Bc_t_split);
verifyEqual(testCase, Bc_t, Bc_s, 'AbsTol', 0.05);
end

function test_a1_smallB_tc_odd_slow(testCase)
% HEADLINE (report): A1 Tc at the 0.05 T proxy with ODD on, production
% 16^3/dpRng 30, via invzt_tc_pm_extrap; compare projected closed form 1.509 K
% (grid-matched legacy convention). The gap = A1 enhancements (retardation
% ~null + transverse RPA + chi_rest).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
g = invzt_qgrid(16, 'legacy_inclusive');
lat = invzt_jq_tensor(ion, g, struct('dpRng', 30, 'cache', true));
critfun = @(T) crit_ok(invzt_solve_point(ion, T, [0.05 0 0], lat, struct('odd', true)));
Tc_a1 = invzt_tc_pm_extrap(critfun, 1.30:1/30:1.75);
[TcP, ~] = invz_odd_zero_field(ion, struct('mode', 'full', 'grids', {{16}}, 'dpRng', 30));
fprintf('Tc ODD: A1 tensor(0.05 T proxy) %.4f K vs projected closed form %.4f K\n', Tc_a1, TcP(1));
verifyTrue(testCase, isfinite(Tc_a1) && Tc_a1 > 1.3 && Tc_a1 < 1.75);
end

% ------------------------------------------------------------------------------
function [c, ok] = crit_ok(pt)
% [crit, ok] contract for a single invzt_solve_point result: a sample votes only
% when it converged to a finite crit (the classifier precondition; the two lowest
% converged PM points, crit > 0, are what invzt_tc_pm_extrap extrapolates).
c = pt.crit;
ok = pt.converged && isfinite(c);
end
