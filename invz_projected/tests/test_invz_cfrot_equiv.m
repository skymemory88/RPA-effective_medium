function tests = test_invz_cfrot_equiv
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function E = diag_spectrum(Hcf, B, ion, o)
% Diagonalize Hcf + Zeeman, sort ascending, shift to the ground state at 0.
% Shared tail for spec_rotcf/spec_rotfield -- these two are cross-checked bit-for-bit,
% so a future change to the diagonalization convention must not drift between them.
C = invz_const();
H = Hcf - ion.gL*C.muB*(B(1)*o.Jx + B(2)*o.Jy + B(3)*o.Jz);
E = sort(real(eig((H + H')/2)));  E = E - E(1);
end

function E = spec_rotcf(ion, r_deg, B)
% Spectrum with the m=4 CF coefficient pairs rotated by r (cf.m 'coefficient' method).
% invz_ion has B44s = 0 implicitly; rotation generates it. Build Hcf manually.
o = stevens_ops(ion.J);
c4 = cosd(4*r_deg);  s4 = sind(4*r_deg);
B44c =  c4*ion.B44;                 B44s = -s4*ion.B44;
B64c =  c4*ion.B64c + s4*ion.B64s;  B64s = -s4*ion.B64c + c4*ion.B64s;
assert(isfield(o, 'O44s') && isfield(o, 'O64s'), 'stevens_ops must expose O44s/O64s');
Hcf = ion.B20*o.O20 + ion.B40*o.O40 + B44c*o.O44 + B44s*o.O44s ...
    + ion.B60*o.O60 + B64c*o.O64c + B64s*o.O64s;
E = diag_spectrum(Hcf, B, ion, o);
end

function E = spec_rotfield(ion, phi_deg, Bmag)
o = stevens_ops(ion.J);
Hcf = ion.B20*o.O20 + ion.B40*o.O40 + ion.B44*o.O44 ...
    + ion.B60*o.O60 + ion.B64c*o.O64c + ion.B64s*o.O64s;
B = Bmag*[cosd(phi_deg) sind(phi_deg) 0];
E = diag_spectrum(Hcf, B, ion, o);
end

function test_pin_cfrot_field_mapping(testCase)
% Discovered and pinned (2026-07-16): coefficient rotation by r == field rotation by +r
% (SAME sign -- not the naive "opposite sign" guess).
%   external ion.cfRot = -11 deg  <=>  invz phi_ab = -11 deg
ion = invz_ion();  r = -11;  Bmag = 4;
Erot = spec_rotcf(ion, r, [Bmag 0 0]);
Ep = spec_rotfield(ion, +r, Bmag);   % s = +1 candidate: matches (hardened, dp ~ 1.8e-14)
Em = spec_rotfield(ion, -r, Bmag);   % s = -1 candidate: must NOT match (dm ~ 0.435)
dp = max(abs(Ep - Erot));  dm = max(abs(Em - Erot));
verifyLessThan(testCase, dp, 1e-10);
verifyGreaterThan(testCase, dm, 1e-6);
end

function test_mapping_holds_off_axis(testCase)
% Same +r mapping at a second rotation angle and field, guards against 4r-aliasing.
ion = invz_ion();  r = 7;  Bmag = 6;
Erot = spec_rotcf(ion, r, [Bmag 0 0]);
dp = max(abs(spec_rotfield(ion, +r, Bmag) - Erot));
dm = max(abs(spec_rotfield(ion, -r, Bmag) - Erot));
verifyLessThan(testCase, dp, 1e-10);
verifyGreaterThan(testCase, dm, 1e-6);
end

function test_pipeline_cfrot_equals_field_rot(testCase)
% End-to-end covariance: rotated CF (invz_cfrot) + field along a must equal the
% unrotated CF + rotated field (phi_ab route) through the FULL 1/z spectra
% pipeline. Requires vector_ab: the cc channel and hyperfine are rotation-
% invariant about c and the vector transverse MF is in-plane isotropic, so the
% CF is the only in-plane-anisotropic ingredient. Catches any hard-wired x-axis
% assumption anywhere in the chain (legacy_x fails this by construction).
ion = invz_ion();  r = -11;  w = (0.05:0.05:0.6).';
fx = struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', 6.4e-3), ...
            'verbose', false, 'solve_opts', struct('transverse_mf', 'vector_ab'));
fxA = fx;  fxA.field_dir = [cosd(r) sind(r) 0];          % rotated field, plain CF
SA = invz_spectra_map(ion, 0.31, [3 5.5], w, fxA);
fxB = fx;  fxB.field_dir = [1 0 0];                      % rotated CF, field along a
SB = invz_spectra_map(invz_cfrot(ion, r), 0.31, [3 5.5], w, fxB);
verifyEqual(testCase, SB.chiz,   SA.chiz,   'AbsTol', 1e-8);
verifyEqual(testCase, SB.chirpa, SA.chirpa, 'AbsTol', 1e-8);
end
