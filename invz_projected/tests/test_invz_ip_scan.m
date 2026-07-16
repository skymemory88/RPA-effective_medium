function tests = test_invz_ip_scan
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_c4fit_recovers_synthetic(testCase)
phi = (0:5:90).';
y = 1 + 0.1*cosd(4*(phi - 17));
[A, phi0, resid] = invz_c4fit(phi, y);
verifyEqual(testCase, A(1), 1.0, 'AbsTol', 1e-10);
verifyEqual(testCase, hypot(A(2), A(3)), 0.1, 'AbsTol', 1e-10);
verifyEqual(testCase, phi0, 17.0, 'AbsTol', 1e-8);
verifyLessThan(testCase, resid, 1e-10);
end

function test_bare_extrema_anchor(testCase)
% Controller-verified digits (2026-07-16): bare CF+Zeeman at 4 T, T = 0.31 K, hyp = false.
ion = invz_ion();
o = struct('hyp', false, 'transverse_mf', 'none');
s15 = invz_single_ion(ion, 0.31, 4*[cosd(15) sind(15) 0], o);
s60 = invz_single_ion(ion, 0.31, 4*[cosd(60) sind(60) 0], o);
verifyEqual(testCase, s15.E(2)-s15.E(1), 0.345693898281, 'AbsTol', 1e-10);
verifyEqual(testCase, s60.E(2)-s60.E(1), 0.370662538778, 'AbsTol', 1e-10);
end

function test_nonaxial_zero_kills_anisotropy(testCase)
% With all m=4 Stevens terms zeroed the in-plane rotation is exact to numerics.
ion = invz_ion();  ion.B44 = 0;  ion.B64c = 0;  ion.B64s = 0;
o = struct('hyp', false, 'transverse_mf', 'none');
d = zeros(7,1);  phis = 0:15:90;
for k = 1:7
    s = invz_single_ion(ion, 0.31, 4*[cosd(phis(k)) sind(phis(k)) 0], o);
    d(k) = s.E(2) - s.E(1);
end
verifyLessThan(testCase, max(d) - min(d), 1e-10);
end

function test_tensor_ref_weight_metrics(testCase)
ion = invz_ion();  w = (0:0.01:0.6).';
o = struct('Jsel', 6.4e-3, 'Jaa0', 3.5e-3, 'eta', 0.02, 'transverse_mf', 'vector_ab');
R = invz_chi_tensor_ref(ion, 0.1, 3*[cosd(20) sind(20) 0], w, o);
verifyTrue(testCase, isfinite(R.eps_W) && R.eps_W >= 0);
verifyTrue(testCase, R.W_ten > 0 && R.W_sc > 0);
end
