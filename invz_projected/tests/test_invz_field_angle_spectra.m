function tests = test_invz_field_angle_spectra
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
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

function test_qpath_scalar_vs_vector_and_metadata(testCase)
ion = invz_ion();  w = (0.1:0.1:0.5).';
qp  = [1 0 0; 1.5 0 0; 2 0 0];
fx  = fixture();  fx.dpRng = 8;                    % small real-space cutoff: fast path sums
S1 = invz_spectra_qpath(ion, 0.31, 5.5, qp, w, fx);
S2 = invz_spectra_qpath(ion, 0.31, [5.5 0 0], qp, w, fx);
verifyEqual(testCase, S2.chiz, S1.chiz);
verifyEqual(testCase, S2.Epeak, S1.Epeak);
verifyEqual(testCase, S1.Bvec, [5.5 0 0]);
verifyEqual(testCase, S1.Bmag, 5.5);
end

function test_qpath_error_message_wellformed_for_vector(testCase)
% Review finding 7: a 3-vector B through a scalar %.3f recycles the format string;
% the message must instead carry one mat2str-rendered vector.
ion = invz_ion();  w = (0.1:0.1:0.5).';
qp  = [1 0 0; 2 0 0];
% (1.9 K, 0.04 T): PM two-level raises invz:degenerateDoublet -> phase 0 -> noSolution,
% thrown BEFORE any path lattice sum (cheap).
try
    invz_spectra_qpath(ion, 1.9, [0.04 0 0], qp, w, fixture());
    verifyTrue(testCase, false, 'expected invz:noSolution');
catch err
    verifyEqual(testCase, err.identifier, 'invz:noSolution');
    verifyTrue(testCase, contains(err.message, mat2str([0.04 0 0], 4)));
end
end

function test_plot_map_axis_label(testCase)
% Review finding 7: the sweep axis is a magnitude; 'B_x (T)' is wrong under a tilt.
S = struct('fields', [1 2], 'w', (0:0.1:1).');
f = figure('Visible', 'off');
restore = onCleanup(@() close(f));
ax = axes(f);
invz_plot_spectra_map(ax, S, rand(11, 2), 'ttl');
verifyEqual(testCase, ax.XLabel.String, '|B| (T)');
end

function test_pm_bz_mirror_symmetry(testCase)
% Spec test 4: chi''_cc is even in Bz (Z2). Metric per second-review refinement 2.
ion = invz_ion();  w = (0:0.02:0.5).';
fp = fixture();  fp.field_dir = [cosd(1) 0 +sind(1)];
fm = fixture();  fm.field_dir = [cosd(1) 0 -sind(1)];
Sp = invz_spectra_map(ion, 0.31, 5.0, w, fp);
Sm = invz_spectra_map(ion, 0.31, 5.0, w, fm);
verifyEqual(testCase, Sp.phase, 1);
verifyEqual(testCase, Sm.phase, 1);
a = Sp.chiz(:);  b = Sm.chiz(:);
verifyTrue(testCase, all(isfinite(a)) && all(isfinite(b)));
verifyLessThan(testCase, max(abs(a - b)) / max(max(abs(a)), 1e-12), 1e-8);
end
