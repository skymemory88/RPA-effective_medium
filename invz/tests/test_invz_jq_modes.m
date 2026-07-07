function tests = test_invz_jq_modes
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));   % repo root: MF_dipole, exchange
end

function test_gamma_point_constants(testCase)
% R 2007 eq (4): J_D*D_aa(0) = 3.912 ueV, J_D*D_cc(0) = 6.821 ueV;
% uniform cc mode with exchange: J(0) = 6.821 + 4*(-0.1) = 6.421 ueV.
ion = invz_ion();
[~, info] = invz_jq_modes(ion, [0 0 0], struct('dpRng', 30, 'cache', false));
verifyEqual(testCase, info.Jcc0_dipole*1e3, 6.821e-3*1e3, 'RelTol', 0.03);  % meV -> ueV displayed
verifyEqual(testCase, info.Jaa0_dipole*1e3, 3.912e-3*1e3, 'RelTol', 0.03);
verifyEqual(testCase, info.Jcc0*1e3, 6.421e-3*1e3, 'RelTol', 0.03);
end

function test_modes_real_and_bounded(testCase)
ion = invz_ion();
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09];
[Jnu, info] = invz_jq_modes(ion, q, struct('dpRng', 20, 'cache', false));
verifyEqual(testCase, size(Jnu), [3 4]);
verifyLessThan(testCase, max(abs(imag(Jnu(:)))), 1e-12);
verifyLessThan(testCase, max(Jnu(:)), info.Jcc0 + 1e-4);   % Γ uniform mode is the max coupling
end

function test_cache_roundtrip(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
opts = struct('dpRng', 10, 'cache', true);
cacheDir = fullfile(fileparts(mfilename('fullpath')), '..', 'cache');
[J1, i1] = invz_jq_modes(ion, q, opts);   % cold: computes and saves
[J2, i2] = invz_jq_modes(ion, q, opts);   % warm: must load identical values
verifyEqual(testCase, J2, J1, 'AbsTol', 0);
verifyEqual(testCase, i2.Jcc0, i1.Jcc0, 'AbsTol', 0);
% different physics params must MISS the cache (different key), not reuse it
ion2 = ion;  ion2.J12 = -0.2e-3;
[~, i3] = invz_jq_modes(ion2, q, opts);
verifyEqual(testCase, i3.Jcc0 - i1.Jcc0, 4*(ion2.J12 - ion.J12), 'RelTol', 1e-9);
% 5% J12 retune collides in the filename hash (single-precision sum) — the
% content verification must force a recompute, not return stale values:
ion3 = ion;  ion3.J12 = -0.105e-3;
[~, i4] = invz_jq_modes(ion3, q, opts);
verifyEqual(testCase, i4.Jcc0 - i1.Jcc0, 4*(ion3.J12 - ion.J12), 'RelTol', 1e-9);
end
