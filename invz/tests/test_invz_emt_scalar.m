function tests = test_invz_emt_scalar
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function [G0, Jf, wn, wts, beta] = synthetic_case()
% Two-level G0 on a simple-cubic nn lattice, subcritical.
C = invz_const();  T = 1.0;
tl.Delta = 0.3;  tl.n01 = tanh(tl.Delta/(2*C.kB*T));  tl.M2 = 30;
[wn, wts, beta] = invz_matsubara(T, 60*tl.Delta);
G0 = -tl.M2*real(invz_g(tl, 1i*wn));
n = 8;  f = 2*pi*(0:n-1)/n;
[qx,qy,qz] = ndgrid(f,f,f);
J0 = 0.8/(tl.M2*2*tl.n01/tl.Delta);                 % J(0)*chi0(0) = 0.8
Jf = (J0/6)*2*(cos(qx(:))+cos(qy(:))+cos(qz(:)));
end

function test_rpa_recovery_at_sigma_zero(testCase)
% HTML Sec 5: with Sigma=0 the effective-medium G(q) equals the RPA G0/(1+J*G0) exactly.
[G0, Jf, wn, ~, ~] = synthetic_case();
med = invz_emt_scalar(G0, zeros(size(G0)), Jf, struct());
verifyTrue(testCase, med.converged);
Grpa_avg = mean(G0.'./(1 + Jf.*G0.'), 1).';        % (1/N)sum_q RPA
verifyEqual(testCase, med.G, Grpa_avg, 'RelTol', 1e-8);
verifyLessThan(testCase, med.closure, 1e-8);
% K decays like wn^-2: check last point much smaller than K(0)
verifyLessThan(testCase, abs(med.K(end)), 0.05*max(abs(med.K(1)), eps));
end

function test_dyson_shift_with_constant_sigma(testCase)
% A constant Sigma acts exactly like G0 -> G0/(1+Sigma) (HTML eq 30).
[G0, Jf, ~, ~, ~] = synthetic_case();
s0 = 0.2;
medS = invz_emt_scalar(G0, s0*ones(size(G0)), Jf, struct());
medR = invz_emt_scalar(G0/(1+s0), zeros(size(G0)), Jf, struct());
verifyEqual(testCase, medS.G, medR.G, 'RelTol', 1e-8);
verifyEqual(testCase, medS.K, medR.K, 'RelTol', 1e-6);
end
