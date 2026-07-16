function tests = test_invz_stevens
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_constants(testCase)
C = invz_const();
verifyEqual(testCase, C.kB,  0.08617333, 'RelTol', 1e-6);
verifyEqual(testCase, C.muB, 0.05788382, 'RelTol', 1e-6);
verifyEqual(testCase, C.gfac, 0.08388, 'RelTol', 2e-3);   % mu0/4pi*(gL*muB)^2 in meV*Ang^3
end

function test_angular_momentum_algebra(testCase)
ops = stevens_ops(8);
verifyEqual(testCase, size(ops.Jz), [17 17]);
comm = ops.Jx*ops.Jy - ops.Jy*ops.Jx;
verifyLessThan(testCase, max(abs(comm - 1i*ops.Jz), [], 'all'), 1e-12);
Jsq = ops.Jx^2 + ops.Jy^2 + ops.Jz^2;
verifyLessThan(testCase, max(abs(Jsq - 8*9*eye(17)), [], 'all'), 1e-10);
end

function test_stevens_J2_explicit(testCase)
% For J=2 (X=6): O20 diagonal = 3m^2-6 for m=2..-2; O44(1,5)=12 (=(1/2)*24)
ops = stevens_ops(2);
verifyEqual(testCase, diag(ops.O20), [6;-3;-6;-3;6], 'AbsTol', 1e-12);
verifyEqual(testCase, real(ops.O44(1,5)), 12, 'AbsTol', 1e-10);
for f = {'O20','O40','O44','O60','O64c','O64s'}
    M = ops.(f{1});
    verifyLessThan(testCase, max(abs(M - M'), [], 'all'), 1e-10);  % Hermitian
end
end
