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
verifyEqual(testCase, info.Jaa0*1e3, 3.512e-3*1e3, 'RelTol', 0.03);  % live transverse J(0) matches R2007
end

function test_modes_real_and_bounded(testCase)
ion = invz_ion();
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09];
[Jnu, info] = invz_jq_modes(ion, q, struct('dpRng', 20, 'cache', false));
verifyEqual(testCase, size(Jnu), [3 4]);
verifyLessThan(testCase, max(abs(imag(Jnu(:)))), 1e-12);
verifyLessThan(testCase, max(Jnu(:)), info.Jcc0 + 1e-4);   % Γ uniform mode is the max coupling
end

function test_demag_shape_term(testCase)
% Ordering-channel couplings are demag-INVARIANT: per R2007 the demagnetization
% field cancels from the critical condition (ordering occurs at q -> 0+, not the
% strict uniform mode), so Jcc0 and the branch spectrum Jnu never see the sample
% shape. The shape enters only through (a) info.Jshape_cc, the strict-uniform
% OBSERVABLE correction, and (b) info.Jaa0, the transverse (non-critical) channel.
ion0 = invz_ion();                                  % demag = 0 default
qs = [0 0 0; 1 0 0; 0.25 0 0.1];                    % Gamma, non-Gamma zone point, generic q
[J0nu, i0] = invz_jq_modes(ion0, qs, struct('dpRng',30,'cache',false));
verifyEqual(testCase, i0.Jshape_cc, 0);             % off: no shape correction
verifyEqual(testCase, i0.Jaa0, i0.Jaa0_dipole + 4*ion0.J12);

C = invz_const();  lorz4 = 4 * (4*pi/(3*ion0.Vc)*C.gfac);  % uniform-mode Lorentz share

% sphere (Nz = Nx = 1/3): ordering channel must NOT move; observable and
% transverse channels must. For a sphere dm_cc = dm_aa = lorz, so
% Jshape_cc = 4*lorz and the transverse dipole share drops by the same amount.
ionS = invz_ion();  ionS.demag = 1;  ionS.alpha = 1;
[JSnu, iS] = invz_jq_modes(ionS, qs, struct('dpRng',30,'cache',false));
verifyEqual(testCase, JSnu, J0nu);                  % branch spectrum bit-identical
verifyEqual(testCase, iS.Jcc0, i0.Jcc0);            % criticality coupling bit-identical
verifyEqual(testCase, iS.Jcc0_dipole, i0.Jcc0_dipole);
verifyEqual(testCase, iS.Jshape_cc, lorz4, 'RelTol', 1e-12);
verifyEqual(testCase, iS.Jaa0_dipole, i0.Jaa0_dipole - lorz4, 'RelTol', 1e-9);
verifyEqual(testCase, iS.Jaa0, iS.Jaa0_dipole + 4*ionS.J12);

% c-axis needle (Nz = 0, Nx = 1/2): no cc shape term; transverse drops by
% 4*dm_aa = 4*gfac*(4pi/Vc)*(1/2) = 1.5*lorz4.
ionN = invz_ion();  ionN.demag = 1;  ionN.alpha = 0;
[~, iN] = invz_jq_modes(ionN, [0 0 0], struct('dpRng',30,'cache',false));
verifyEqual(testCase, iN.Jshape_cc, 0, 'AbsTol', 1e-15);
verifyEqual(testCase, iN.Jcc0, i0.Jcc0);
verifyEqual(testCase, iN.Jaa0_dipole, i0.Jaa0_dipole - 1.5*lorz4, 'RelTol', 1e-9);

% cache must distinguish demag settings via the SHAPE fields (Jcc0 no longer differs)
opts = struct('dpRng',10,'cache',true);
[~, c0] = invz_jq_modes(invz_ion(), [0 0 0], opts);
ionS2 = invz_ion();  ionS2.demag = 1;  ionS2.alpha = 1;
[~, cS] = invz_jq_modes(ionS2, [0 0 0], opts);
verifyGreaterThan(testCase, abs(cS.Jshape_cc - c0.Jshape_cc), 1e-6);
verifyEqual(testCase, cS.Jcc0, c0.Jcc0);
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
