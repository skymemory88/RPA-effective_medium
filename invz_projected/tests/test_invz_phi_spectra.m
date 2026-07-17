function tests = test_invz_phi_spectra
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function fx = fixture()
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', ...
            'info', struct('Jcc0', 6.4e-3), 'verbose', false);
end

function test_map_rejects_legacy_with_by(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(20) sind(20) 0];
verifyError(testCase, @() invz_spectra_map(ion, 0.31, [3 5.5], w, fx), ...
            'invz:transverseMF');
end

function test_qpath_rejects_legacy_with_by(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
qp = [0 0 0; 0.5 0 0];
verifyError(testCase, ...
    @() invz_spectra_qpath(ion, 0.31, 3*[cosd(20) sind(20) 0], qp, w, struct()), ...
    'invz:transverseMF');
end

function test_map_vector_mode_accepted(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
fx = fixture();  fx.field_dir = [cosd(20) sind(20) 0];
fx.solve_opts = struct('transverse_mf', 'vector_ab');
S = invz_spectra_map(ion, 0.31, [3 5.5], w, fx);
verifyEqual(testCase, S.transverse_mf, 'vector_ab');
verifyTrue(testCase, any(isfinite(S.chiz(:))));
end

function test_map_metadata_default_legacy(testCase)
ion = invz_ion();  w = (0.02:0.02:0.6).';
S = invz_spectra_map(ion, 0.31, [3 5.5], w, fixture());
verifyEqual(testCase, S.transverse_mf, 'legacy_x');
end

function test_map_c4_phi_plus_90(testCase)
% phi and phi+90 give the same cc spectrum in vector mode (C4 exact).
ion = invz_ion();  w = (0.05:0.05:0.6).';
fx = fixture();  fx.solve_opts = struct('transverse_mf', 'vector_ab');
fx.field_dir = [cosd(20) sind(20) 0];
S1 = invz_spectra_map(ion, 0.31, 3.0, w, fx);
fx.field_dir = [cosd(110) sind(110) 0];
S2 = invz_spectra_map(ion, 0.31, 3.0, w, fx);
verifyEqual(testCase, S2.chiz, S1.chiz, 'AbsTol', 1e-8);
end

function test_qpath_c4_phi_plus_90(testCase)
% C1 regression: the 1/z q-path chi0 must follow transverse_mf (a legacy_x
% rebuild at a rotated field breaks C4 between phi and phi+90).
ion = invz_ion();  T = 0.31;  Bmag = 5.5;  w = (0.05:0.05:0.6).';
fx = fixture();  fx.dpRng = 10;
fx.solve_opts = struct('transverse_mf', 'vector_ab');
qp = [1 0 0; 2 0 0];
S1 = invz_spectra_qpath(ion, T, Bmag*[cosd(20)  sind(20)  0], qp, w, fx);
S2 = invz_spectra_qpath(ion, T, Bmag*[cosd(110) sind(110) 0], qp, w, fx);
verifyEqual(testCase, S1.phase, 2);   % strict-PM: phase==2 branch under test
verifyEqual(testCase, S2.chiz, S1.chiz, 'AbsTol', 1e-8);
end
