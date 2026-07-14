function tests = test_invz_sigma_ordered
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
end

function [tl, g, K, beta, lam] = fixture(m)
% Synthetic two-level + moments; m is the ordered diagonal moment.
C = invz_const();  T = 1.0;  beta = 1/(C.kB*T);
tl.Delta = 0.3;  tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.M2 = 25;  tl.m = m;  tl.g0 = 2*tl.n01/tl.Delta;
[wn, wts, ~] = invz_matsubara(T, 800*tl.Delta);
g = real(invz_g(tl, 1i*wn));
K = 2e-3 ./ (1 + (wn/0.5).^2);
lam = invz_lambdas(K, g, wts, beta, [1 2 3]);
end

function test_reduces_to_paramagnet_at_m0(testCase)
% At m = 0 the ordered self-energy must equal the paramagnetic invz_sigma exactly.
[tl, g, K, beta, lam] = fixture(0);
so = invz_sigma_ordered(tl, lam, K, g, beta);
sp = invz_sigma(tl, lam(1:2), K, g, beta);
verifyEqual(testCase, so.alpha_m, 0, 'AbsTol', 1e-14);
verifyEqual(testCase, so.Sigma, sp.Sigma, 'RelTol', 1e-12, 'AbsTol', 1e-14);
end

function test_alpha_m_scales_as_m_squared(testCase)
% alpha_m is linear in m^2 (its overall prefactor): doubling m quadruples alpha_m.
[tl1, g1, K1, b1, l1] = fixture(1.0);
[tl2, g2, K2, b2, l2] = fixture(2.0);
s1 = invz_sigma_ordered(tl1, l1, K1, g1, b1);
s2 = invz_sigma_ordered(tl2, l2, K2, g2, b2);
verifyGreaterThan(testCase, abs(s1.alpha_m), 0);
verifyEqual(testCase, s2.alpha_m, 4*s1.alpha_m, 'RelTol', 1e-12);
end

function test_ordered_sigma_differs_from_para(testCase)
% With m ~= 0 the ordered self-energy must actually differ from the paramagnetic one.
[tl, g, K, beta, lam] = fixture(1.5);
so = invz_sigma_ordered(tl, lam, K, g, beta);
sp = invz_sigma(tl, lam(1:2), K, g, beta);
verifyGreaterThan(testCase, max(abs(so.Sigma - sp.Sigma)), 1e-6);
end
