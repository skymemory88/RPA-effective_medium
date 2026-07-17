function tests = test_invz_sigma
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [tl, wn, wts, beta, g, K] = fixture(Ecut_factor)
C = invz_const();  T = 1.0;
tl.Delta = 0.3;  tl.n01 = tanh(tl.Delta/(2*C.kB*T));
tl.M2 = 30;  tl.m = 0;  tl.g0 = 2*tl.n01/tl.Delta;
[wn, wts, beta] = invz_matsubara(T, Ecut_factor*tl.Delta);
g = real(invz_g(tl, 1i*wn));
K = 2e-3 ./ (1 + (wn/0.5).^2);          % synthetic even K decaying like wn^-2
end

function test_lambda_analytic(testCase)
% With K = const = c: lambda_1 = c*(1/beta)*sum g = c; lambda_2 = c*(1/2)[g0+beta(1-n01^2)].
% The (1/beta)*sum g = 1 sum rule is only exact for the infinite Matsubara sum, so
% Ecut_factor=800 (not 400) is needed to bring the truncation residual under the 1e-3 RelTol.
[tl, ~, wts, beta, g, ~] = fixture(800);
c = 3.7e-3;
lam = invz_lambdas(c*ones(size(g)), g, wts, beta, [1 2]);
verifyEqual(testCase, lam(1), c, 'RelTol', 1e-3);
verifyEqual(testCase, lam(2), c*0.5*(tl.g0 + beta*(1-tl.n01^2)), 'RelTol', 1e-6);
end

function test_sumrule_cancellation(testCase)
% THE identity that fixes alpha (HTML eq 24): (1/beta) sum_n G0.*(K.*G0 + Sigma) = 0.
[tl, ~, wts, beta, g, K] = fixture(400);
G0 = -tl.M2*g;
lam = invz_lambdas(K, g, wts, beta, [1 2]);
sig = invz_sigma(tl, lam, K, g, beta);
lhs   = sum(wts .* G0 .* (K.*G0 + sig.Sigma)) / beta;
scale = sum(wts .* abs(G0 .* K .* G0)) / beta;
verifyLessThan(testCase, abs(lhs)/scale, 1e-3);
% convergence order: residual falls when the grid doubles
[tl2, ~, wts2, beta2, g2, K2] = fixture(800);
G02 = -tl2.M2*g2;
lam2 = invz_lambdas(K2, g2, wts2, beta2, [1 2]);
sig2 = invz_sigma(tl2, lam2, K2, g2, beta2);
lhs2 = sum(wts2 .* G02 .* (K2.*G02 + sig2.Sigma)) / beta2;
verifyLessThan(testCase, abs(lhs2), abs(lhs));
end

function test_sigma_zero_when_K_zero(testCase)
[tl, ~, wts, beta, g, ~] = fixture(400);
lam = invz_lambdas(zeros(size(g)), g, wts, beta, [1 2]);
sig = invz_sigma(tl, lam, zeros(size(g)), g, beta);
verifyEqual(testCase, sig.alpha, 0, 'AbsTol', 1e-15);
verifyLessThan(testCase, max(abs(sig.Sigma)), 1e-15);
end
