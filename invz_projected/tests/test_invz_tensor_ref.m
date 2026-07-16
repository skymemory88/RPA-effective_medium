function tests = test_invz_tensor_ref
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_cross_channel_hierarchy_at_zero_tilt(testCase)
% Measured physics (Task-9 blocked round, amended 2026-07-16): at theta = 0,
% B || a, the yz cross channel of chi0 is SYMMETRY-ALLOWED (B64s) and large,
% while xz is strongly suppressed (measured yz/zz = 0.183, xz/zz = 2.8e-3 at
% 6 T, hyp = false). The scalar-vs-tensor baseline discrepancy at zero tilt is
% therefore real and finite -- it is REPORTED as a baseline, never gated as a
% tilt error.
ion = invz_ion();
si = invz_single_ion(ion, 0.31, [6 0 0], struct('hyp', false));
z = (0:0.01:0.5).' + 1i*5e-3;
c0 = invz_chi0z(si, 0.31, z, struct('elastic', false));
xz = max(abs(squeeze(c0(1,3,:))));  yz = max(abs(squeeze(c0(2,3,:))));
zz = max(abs(squeeze(c0(3,3,:))));
verifyGreaterThan(testCase, yz, 10*xz);          % yz-dominated cross-channel hierarchy
verifyLessThan(testCase, xz, 0.02*zz);           % xz strongly suppressed at theta = 0
w = (0:0.01:0.5).';
R = invz_chi_tensor_ref(ion, 0.31, [6 0 0], w, struct('hyp', false, 'eta', 0.02));
verifyGreaterThan(testCase, R.eps_spec, 0);      % finite yz-driven baseline ...
verifyLessThan(testCase, R.eps_spec, 1);         % ... bounded (sanity)
end

function test_tilt_metric_wellformed(testCase)
% eps_tilt (invz_tilt_err): error in the tilt-induced change, baseline-
% differenced so the theta-independent yz discrepancy drops out. Reported as a
% diagnostic (round-2 resolution, 7f9f16b): every L2 lineshape metric,
% including eps_tilt, is dominated by the zero-tilt yz peak-offset artifact for
% these sharp lines, so the GATED metric is the peak-observable eps_amp
% (asserted finite/nonnegative below as a cheap sanity check of the new field).
ion = invz_ion();
w = (0:0.01:0.5).';
o = struct('hyp', false, 'eta', 0.02);
R0 = invz_chi_tensor_ref(ion, 0.31, [6 0 0], w, o);
R1 = invz_chi_tensor_ref(ion, 0.31, 6*[cosd(1) 0 sind(1)], w, o);
e1 = invz_tilt_err(R1, R0);
verifyGreaterThan(testCase, e1, 0);
verifyLessThan(testCase, e1, 2);
verifyGreaterThanOrEqual(testCase, R1.eps_amp, 0);
verifyTrue(testCase, isfinite(R1.eps_amp));
end

function test_demag_guard(testCase)
ion = invz_ion();  ion.demag = 1;
verifyError(testCase, ...
    @() invz_chi_tensor_ref(ion, 0.31, [6 0 0], (0:0.1:0.4).'), 'invz:tensorRef');
end

function test_reproducibility_of_logged_table(testCase)
% Slow: re-measures three logged (angle, field) points and asserts the eps_amp
% (GATED metric) and eps_tilt (diagnostic) values match
% docs/SESSION-2026-07-16-field-angle.md to 1% (reproducibility, NOT a size
% target).
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'slow test: set INVZ_SLOW=1');
ion = invz_ion();
w = (0:0.005:0.6).';
o = struct('eta', 0.02);      % default couplings -- must match the logged convention
pts           = {0.5, 6;  2, 6;  5, 4.95};        % {theta_deg, B} spot checks
expected_amp  = [0.006983726953; 0.01134694549; 0.05398975621];   % measured eps_amp (GATED metric)
expected_tilt = [0.1139439749;   0.08803939767;  0.07649402387];  % measured eps_tilt (diagnostic)
for k = 1:size(pts, 1)
    a = pts{k, 1};  B = pts{k, 2};
    R0 = invz_chi_tensor_ref(ion, 0.1, [B 0 0], w, o);
    R  = invz_chi_tensor_ref(ion, 0.1, B*[cosd(a) 0 sind(a)], w, o);
    verifyEqual(testCase, R.eps_amp, expected_amp(k), 'RelTol', 0.01);
    verifyEqual(testCase, invz_tilt_err(R, R0), expected_tilt(k), 'RelTol', 0.01);
end
end
