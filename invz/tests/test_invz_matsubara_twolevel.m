function tests = test_invz_matsubara_twolevel
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_grid(testCase)
[wn, wts, beta] = invz_matsubara(2.0, 30);
verifyEqual(testCase, wn(1), 0);
verifyEqual(testCase, wn(2), 2*pi/beta, 'RelTol', 1e-12);
verifyEqual(testCase, wts(1), 1);  verifyEqual(testCase, wts(2), 2);
verifyGreaterThan(testCase, wn(end), 30);
end

function test_g_sum_identities(testCase)
% (1/beta)*sum_n g = 1  and  (1/beta)*sum_n g^2 = (1/2)[g(0)+beta(1-n01^2)]  (HTML eqs 26, Sec 7)
tl.Delta = 0.4;  T = 1.3;  C = invz_const();  beta = 1/(C.kB*T);
tl.n01 = tanh(beta*tl.Delta/2);  tl.M2 = 1;  tl.m = 0;  tl.g0 = 2*tl.n01/tl.Delta;
[wn, wts, ~] = invz_matsubara(T, 2000*tl.Delta);   % tail ~3e-4 < RelTol 1e-3 (brief's 400*Delta left a 1.5e-3 tail)
g = real(invz_g(tl, 1i*wn));
s1 = sum(wts.*g)/beta;
s2 = sum(wts.*g.^2)/beta;
verifyEqual(testCase, s1, 1, 'RelTol', 1e-3);
verifyEqual(testCase, s2, 0.5*(tl.g0 + beta*(1-tl.n01^2)), 'RelTol', 1e-6);
end

function test_twolevel_from_ion(testCase)
ion = invz_ion();  C = invz_const();
T = 0.31;  Bx = 4.0;
tl = invz_twolevel(ion, T, Bx);
si = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false));
verifyEqual(testCase, tl.Delta, si.E(2)-si.E(1), 'AbsTol', 1e-12);
verifyEqual(testCase, tl.M2, abs(si.Mz(1,2))^2, 'RelTol', 1e-12);
verifyLessThan(testCase, abs(tl.m), 1e-6);
verifyEqual(testCase, tl.n01, tanh(tl.Delta/(2*C.kB*T)), 'RelTol', 1e-12);
verifyGreaterThan(testCase, tl.M2, 1);       % sizeable Ising matrix element
end
