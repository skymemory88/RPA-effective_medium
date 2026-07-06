function tests = test_invz_solve_point
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function Jf = toy_couplings(J0)
% Cheap stand-in lattice (sc-like) with uniform-mode coupling J0, so the test
% does not depend on slow dipole sums.
n = 10;  f = 2*pi*(0:n-1)/n;
[qx,qy,qz] = ndgrid(f,f,f);
Jf = (J0/6)*2*(cos(qx(:))+cos(qy(:))+cos(qz(:)));
Jf = Jf(2:end);                                   % drop uniform mode (measure zero)
end

function test_converges_and_is_physical(testCase)
ion = invz_ion();
Jf = toy_couplings(ion.J0eff);
pt = invz_solve_point(ion, 1.0, 4.0, Jf, struct('hyp', true));
verifyTrue(testCase, pt.converged);
verifyGreaterThan(testCase, pt.Sigma0, 0);        % fluctuations suppress the response
verifyLessThan(testCase, pt.Sigma0, 1.0);
verifyGreaterThan(testCase, pt.crit, 0);          % paramagnetic at 1 K, 4 T
verifyLessThan(testCase, pt.sumrule_rel, 0.10);   % Jensen: resummed G obeys sum rule approximately
end

function test_approach_to_criticality(testCase)
ion = invz_ion();
Jf = toy_couplings(ion.J0eff);
s = zeros(1,3);  c = zeros(1,3);  Ts = [2.0 1.0 0.5];
for k = 1:3
    pt = invz_solve_point(ion, Ts(k), 4.0, Jf, struct('hyp', true));
    s(k) = pt.Sigma0;
    c(k) = pt.crit;
end
% Sigma0 rises from its high-T decay and saturates to a finite zero-point value
% near T->0 (it may peak in between — monotonicity is NOT guaranteed; the brief's
% original assertion was wrong physics, see task report):
verifyGreaterThan(testCase, s(2), s(1));
verifyGreaterThan(testCase, s(3), s(1));
% The unambiguous "approach to criticality" statement is that the critical
% denominator decreases monotonically on cooling:
verifyGreaterThan(testCase, c(1), c(2));
verifyGreaterThan(testCase, c(2), c(3));
verifyGreaterThan(testCase, c(3), 0);
end

function test_interaction_off_gives_zero_sigma(testCase)
ion = invz_ion();
pt = invz_solve_point(ion, 1.0, 4.0, 1e-12*ones(100,1), struct('hyp', false));
verifyLessThan(testCase, abs(pt.Sigma0), 1e-8);
end
