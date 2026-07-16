function tests = test_invz_transverse_mf
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_legacy_default_bitforbit(testCase)
% Omitted opt and explicit 'legacy_x' take the identical code path.
ion = invz_ion();
s1 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false));
s2 = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'legacy_x'));
verifyTrue(testCase, isequaln(s1, s2));
% digit anchors (controller-verified 2026-07-16): guards against accidental rebaseline
verifyEqual(testCase, s1.E(2)-s1.E(1), 0.369235620278, 'AbsTol', 1e-11);
verifyEqual(testCase, s1.Jexp(2), -0.0689991117463, 'AbsTol', 1e-10);
verifyEqual(testCase, s1.hy, 0);
verifyEqual(testCase, s1.transverse_mf, 'legacy_x');
end

function test_invalid_mode_errors(testCase)
ion = invz_ion();
verifyError(testCase, ...
    @() invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'diagonal')), ...
    'invz:transverseMF');
end

function test_none_equals_bare(testCase)
% 'none' == the current Jxx0 = 0 bare CF+Zeeman calculation.
ion = invz_ion();
sn = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'none'));
sb = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'Jxx0', 0));
verifyEqual(testCase, sn.E,    sb.E,    'AbsTol', 1e-14);
verifyEqual(testCase, sn.Jexp, sb.Jexp, 'AbsTol', 1e-14);
verifyEqual(testCase, sn.hx, 0);  verifyEqual(testCase, sn.hy, 0);
end

function test_vector_c4_axes(testCase)
% C4 is exact for the CF: a/b-axis fields must be equivalent under vector MF.
% Rotation sense pinned: +90 deg about c maps x->y, so <J>([0 B 0]) = Rz(90)*<J>([B 0 0]).
ion = invz_ion();  o = struct('hyp', false, 'transverse_mf', 'vector_ab');
sx = invz_single_ion(ion, 0.31, [4 0 0], o);
sy = invz_single_ion(ion, 0.31, [0 4 0], o);
verifyEqual(testCase, sy.E, sx.E, 'AbsTol', 1e-10);
verifyEqual(testCase, sy.Jexp(2),  sx.Jexp(1), 'AbsTol', 1e-9);   % <Jy>' = <Jx>
verifyEqual(testCase, sy.Jexp(1), -sx.Jexp(2), 'AbsTol', 1e-9);   % <Jx>' = -<Jy>
end

function test_vector_hy_selfconsistent(testCase)
ion = invz_ion();
s = invz_single_ion(ion, 0.31, [4 0 0], struct('hyp', false, 'transverse_mf', 'vector_ab'));
verifyTrue(testCase, s.mf_converged);
verifyEqual(testCase, s.hy, ion.Jxx0*s.Jexp(2), 'AbsTol', 1e-12);
verifyTrue(testCase, abs(s.hy) > 1e-6);   % B64s => <Jy> ~= 0 even for an x-axis field
end

function test_vector_inplane_Jz_zero(testCase)
% Theta x C2 protects <Jz> = 0 exactly for any in-plane field (para path).
ion = invz_ion();
s = invz_single_ion(ion, 0.31, 4*[cosd(30) sind(30) 0], ...
                    struct('hyp', false, 'transverse_mf', 'vector_ab'));
verifyEqual(testCase, s.Jexp(3), 0, 'AbsTol', 1e-9);
end

function test_vector_Fmf_c4_invariant(testCase)
% F_mf must include the 0.5*hy*<Jy> term: under C4 the (hx,hy) weight swaps
% between channels, so F_mf([4 0 0]) == F_mf([0 4 0]) only if hy is counted.
ion = invz_ion();  o = struct('hyp', false, 'transverse_mf', 'vector_ab');
sx = invz_single_ion(ion, 0.31, [4 0 0], o);
sy = invz_single_ion(ion, 0.31, [0 4 0], o);
verifyEqual(testCase, sy.F_mf, sx.F_mf, 'AbsTol', 1e-10);
% self-consistency identity: 0.5*(hx<Jx>+hy<Jy>) == (hx^2+hy^2)/(2*Jxx0)
lhs = 0.5*(sx.hx*sx.Jexp(1) + sx.hy*sx.Jexp(2));
verifyEqual(testCase, lhs, (sx.hx^2 + sx.hy^2)/(2*ion.Jxx0), 'AbsTol', 1e-12);
end

function test_vector_ordered_mode(testCase)
% hy iterates alongside hz in order mode.
ion = invz_ion();
s = invz_single_ion(ion, 0.31, [1 0.5 0], ...
    struct('hyp', false, 'order', true, 'transverse_mf', 'vector_ab'));
verifyTrue(testCase, s.mf_converged);
verifyEqual(testCase, s.hy, ion.Jxx0*s.Jexp(2), 'AbsTol', 1e-12);
end

function test_vector_hzfixed_mode(testCase)
% hy iterates alongside a held-fixed hz; F_mf stays NaN under hz_fixed.
ion = invz_ion();
s = invz_single_ion(ion, 0.31, [2 1 0], ...
    struct('hyp', false, 'hz_fixed', 0.01, 'transverse_mf', 'vector_ab'));
verifyTrue(testCase, s.mf_converged);
verifyEqual(testCase, s.hy, ion.Jxx0*s.Jexp(2), 'AbsTol', 1e-12);
verifyTrue(testCase, isnan(s.F_mf));
end
