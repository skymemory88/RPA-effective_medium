function tests = test_invz_field_angle_spectra
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function fx = fixture()
% Synthetic-coupling fixture (matches test_invz_spectra_map: fast, no lattice sum).
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', ...
            'info', struct('Jcc0', 6.4e-3), 'verbose', false);
end

function test_default_dir_equals_explicit_x(testCase)
% Spec test 2 (map form): default field_dir is exactly [1 0 0].
ion = invz_ion();  T = 0.31;  fields = [2.0 5.5];  w = (0.02:0.02:0.6).';
S1 = invz_spectra_map(ion, T, fields, w, fixture());
fx = fixture();  fx.field_dir = [1 0 0];
S2 = invz_spectra_map(ion, T, fields, w, fx);
verifyEqual(testCase, S2.chiz, S1.chiz);
verifyEqual(testCase, S2.chirpa, S1.chirpa);
verifyEqual(testCase, S1.field_dir, [1 0 0]);
verifyEqual(testCase, S1.Bvec, [fields(:) zeros(2, 2)]);
end

function test_api_validation(testCase)
ion = invz_ion();  w = (0.1:0.1:0.5).';
fx = fixture();  fx.solve_opts = struct('Jxx0', 1e-3);            % reserved field
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, fx), 'invz:solveOpts');
fx = fixture();  fx.field_dir = [0 0 0];
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 5.5, w, fx), 'invz:fieldDir');
verifyError(testCase, @() invz_spectra_map(ion, 0.31, -1.0, w, fixture()), 'invz:fields');
end

function test_deadband_zeroes_tiny_tilt(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [1 0 1e-12];       % Bz = 5.5e-12 T << bz_tol
S = invz_spectra_map(ion, 0.31, 5.5, w, fx);
verifyEqual(testCase, S.Bvec(3), 0);               % dead band applied BEFORE the solve
verifyEqual(testCase, S.phase, 2);                 % genuinely transverse: strict PM
end

function test_tilted_sweep_is_moment_form(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(1) 0 sind(1)];
S = invz_spectra_map(ion, 0.31, [2.0 5.5], w, fx);
verifyEqual(testCase, S.phase, [1 1]);             % field-induced moment at BOTH fields
verifyTrue(testCase, all(isfinite(S.chiz(:))));
verifyEqual(testCase, S.Bvec, [2.0; 5.5] * [cosd(1) 0 sind(1)], 'RelTol', 1e-12);
end

function test_failed_longitudinal_masks_not_crashes(testCase)
% Spec test 8 (review finding 5): a non-converged longitudinal point must yield a
% masked 1/z column plus an RPA-only overlay -- and must NOT abort the parfor by
% falling through to the strict-paramagnet invz_twolevel gate.
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(1) 0 sind(1)];
fx.solve_opts = struct('max_outer', 1);            % cripple the outer EMT loop
S = invz_spectra_map(ion, 0.31, 2.0, w, fx);
verifyEqual(testCase, S.phase, 0);
verifyTrue(testCase, all(~isfinite(S.chiz(:))));   % 1/z column masked
verifyTrue(testCase, any(isfinite(S.chirpa(:))));  % RPA overlay from the failed pto
end
