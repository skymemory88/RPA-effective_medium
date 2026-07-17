function tests = test_invz_emt_scalar
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
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
med = invz_emt_scalar(G0, zeros(size(G0)), Jf, struct('debug', true));
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

function res = emt_residual(G0, Sigma, Jf, med)
% Global-relative residual of the returned (G,K) against R eq 8, recomputed
% independently of how the solver found them. Normalising by max|K| (rather
% than per-entry with a floor) keeps the metric meaningful where K decays to
% ~1e-7 at high Matsubara frequency.
G0 = G0(:); Sigma = Sigma(:); Jf = Jf(:);
G  = G0 ./ (1 + Sigma + med.K.*G0);
Gq = G.' ./ (1 + (Jf - med.K.').*G.');
Kn = (Jf.'*Gq).'/numel(Jf) ./ G;
res = max(abs(Kn - med.K)) / max(abs(med.K));
end

function test_solution_is_exact_fixed_point(testCase)
% The returned (G,K) must satisfy the EMT self-consistency to machine
% precision, not merely to the iterative residual tolerance. This is the
% property the closed-form solve guarantees (finding #1); an iteration
% stopped at tol=1e-10 cannot meet a 1e-12 residual check.
[G0, Jf, ~, ~, ~] = synthetic_case();
Sigma = zeros(size(G0));
med = invz_emt_scalar(G0, Sigma, Jf, struct('debug', true));
verifyTrue(testCase, med.converged);
verifyLessThan(testCase, emt_residual(G0, Sigma, Jf, med), 1e-12);
verifyLessThan(testCase, med.closure, 1e-12);          % mean_q Gq == G exactly
% Also exact with a non-zero, frequency-dependent Sigma.
Sigma2 = 0.15 ./ (1:numel(G0)).';
med2 = invz_emt_scalar(G0, Sigma2, Jf, struct('debug', true));
verifyLessThan(testCase, emt_residual(G0, Sigma2, Jf, med2), 1e-12);
verifyLessThan(testCase, med2.closure, 1e-12);
end

function Kref = emt_reference_K(G0, Sigma, Jf)
% Independent brute-force fixed-point solve of R eqs 7-9 for K, run to a far
% tighter tolerance than production, used to pin the closed-form equivalence.
G0 = G0(:); Sigma = Sigma(:); Jf = Jf(:);  nJ = numel(Jf);
K = zeros(size(G0));
for it = 1:5000
    G  = G0 ./ (1 + Sigma + K.*G0);
    Gq = G.' ./ (1 + (Jf - K.').*G.');
    Kn = (Jf.'*Gq).'/nJ ./ G;
    if max(abs(Kn - K)) < 1e-14, K = Kn; break; end
    K = K + 0.5*(Kn - K);
end
Kref = K;
end

function test_matches_reference_iteration(testCase)
% The solver must land on the same fixed point an independent, tightly
% converged iteration finds, across sub- and super-critical couplings. This
% is the equivalence the closed-form replacement (finding #1) must preserve.
[G0, Jf, ~, ~, ~] = synthetic_case();
Sigma = zeros(size(G0));
for f = [0.5, 0.95, 1.6]
    med = invz_emt_scalar(G0, Sigma, f*Jf, struct());
    verifyTrue(testCase, med.converged);
    Kref = emt_reference_K(G0, Sigma, f*Jf);
    verifyEqual(testCase, med.K, Kref, 'RelTol', 1e-7, ...
        sprintf('K mismatch at coupling factor %.2f', f));
end
end

function test_nonfinite_denominator_not_converged(testCase)
% A genuinely singular medium (1 + J*G0 = 0 at some frequency, here G0=-1/J)
% makes K non-finite. Both the iterative and closed-form solves must report
% this as NOT converged so invz_solve_point can treat it as no solution.
Jf = 0.5;                           % single q-point
G0 = [-1/Jf; -0.4; -0.25];          % first entry is singular, rest are fine
med = invz_emt_scalar(G0, zeros(size(G0)), Jf, struct());
verifyFalse(testCase, med.converged);
verifyFalse(testCase, all(isfinite(med.K)));
end
