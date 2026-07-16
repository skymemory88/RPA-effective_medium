function tests = test_invz_field_angle_slow
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', ...
            'info', struct('Jcc0', 6.4e-3), 'verbose', false);
end

function test_crossover_continuity(testCase)
% Spec test 5: under a 0.5-degree tilt the former sharp transition is a rounded
% crossover -- no censored/masked points across the old Bc, sum rule intact.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'slow test: set INVZ_SLOW=1');
ion = invz_ion();  T = 0.31;  w = (0:0.02:0.5).';
fields = 4.6:0.05:5.3;
fx = fixture();  fx.field_dir = [cosd(0.5) 0 sind(0.5)];
S = invz_spectra_map(ion, T, fields, w, fx);
verifyTrue(testCase, all(S.phase == 1));                     % moment-form everywhere
verifyTrue(testCase, all(isfinite(S.Epeak)));                % no censored crossover gap
% sum rule at every field (pt-level; same fixture and route as the map)
sopts = struct('hyp', true, 'J0eff', fx.info.Jcc0, 'Jxx0', ion.Jxx0);
for k = 1:numel(fields)
    [pt, phase] = invz_solve_auto(ion, T, fields(k)*fx.field_dir, fx.Jnu, sopts);
    verifyEqual(testCase, phase, 1, sprintf('field %.2f', fields(k)));
    verifyLessThan(testCase, pt.sumrule_rel, 5e-2, sprintf('field %.2f', fields(k)));
end
end

function test_theta_to_zero_continuity(testCase)
% Spec test 6 (amended 2026-07-16, recorded justification -- see spec):
%   B = 6 T (PM at zero tilt): chi_cc even in Bz -> flat 1e-6 bound (measured
%     4.2e-7 at 1e-3 deg). Exercises forced moment-form -> strict-PM reduction.
%   B = 3 T (spontaneous FM at zero tilt; 2 T sits in a fixture-specific
%     non-convergence island B in [1,2] T): the single-domain response is
%     LINEAR in Bz (aligned branch breaks Z2; soft-mode coefficient ~ delta/eta,
%     measured 7.7e-3 at 1e-3 deg) -- a flat 1e-6 bound is the wrong physics.
%     Assert continuity as linear scaling: no O(1) jump at the spontaneous ->
%     forced routing boundary.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'slow test: set INVZ_SLOW=1');
ion = invz_ion();  T = 0.31;  w = (0:0.02:0.5).';
% --- PM side: B = 6 T, flat bound ---
S0 = invz_spectra_map(ion, T, 6.0, w, fixture());
ft = fixture();  ft.field_dir = [cosd(1e-3) 0 sind(1e-3)];
St = invz_spectra_map(ion, T, 6.0, w, ft);
a = S0.chiz(:);  b = St.chiz(:);
verifyTrue(testCase, all(isfinite(a)) && all(isfinite(b)), 'B = 6');
verifyLessThan(testCase, max(abs(b - a)) / max(abs(a)), 1e-6, 'B = 6');
% --- FM side: B = 3 T, linear-scaling assertion ---
S0 = invz_spectra_map(ion, T, 3.0, w, fixture());
angs = [1e-3 1e-4];  d = nan(1, 2);
for k = 1:2
    ft = fixture();  ft.field_dir = [cosd(angs(k)) 0 sind(angs(k))];
    St = invz_spectra_map(ion, T, 3.0, w, ft);
    a = S0.chiz(:);  b = St.chiz(:);
    verifyTrue(testCase, all(isfinite(a)) && all(isfinite(b)), sprintf('B = 3, %g deg', angs(k)));
    d(k) = max(abs(b - a)) / max(abs(a));
end
verifyLessThan(testCase, d(1), 3e-2);
r = d(1) / d(2);                                   % pure linear response -> 10
verifyGreaterThan(testCase, r, 6);
verifyLessThan(testCase, r, 15);
end
