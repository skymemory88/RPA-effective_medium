function tests = test_invz_chi0z
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
end

function test_matsubara_equal_time_sum_rule(testCase)
% (1/beta)*sum_n chi_cc(iwn) = <Jz~^2>  (HTML eq 31 with chi = -G), electronuclear ion.
ion = invz_ion(); C = invz_const();
T = 4.0; beta = 1/(C.kB*T);
si = invz_single_ion(ion, T, [2 0 0], struct('hyp', true));
Ecut = 150;                                    % meV; CF spectrum tops out near 40 meV
nmax = ceil(Ecut*beta/(2*pi));
wn  = 2*pi*(0:nmax).'/beta;  wts = [1; 2*ones(nmax,1)];
chi = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
s   = real(squeeze(chi(3,3,:)));
lhs = sum(wts.*s)/beta;
verifyEqual(testCase, lhs, si.JzJz_fluct, 'RelTol', 0.03);
end

function test_even_real_on_matsubara_axis(testCase)
ion = invz_ion();
si = invz_single_ion(ion, 1.0, [3 0 0], struct('hyp', false));
wn = [0.3 0.9];
cp = invz_chi0z(si, 1.0,  1i*wn, struct());
cm = invz_chi0z(si, 1.0, -1i*wn, struct());
verifyLessThan(testCase, max(abs(imag(cp(3,3,:)))), 1e-12);
verifyEqual(testCase, cp(3,3,:), cm(3,3,:), 'AbsTol', 1e-12);
end

function test_static_positive_and_hyp_enhancement(testCase)
% chi_J(0)/chi(0) ≈ 0.77 at 0.31 K near the critical field (R 2007, Sec IIC).
ion = invz_ion();
T = 0.31; Bx = 4.24;
sJ = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false));
sF = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));
cJ = real(invz_chi0z(sJ, T, 0, struct()));  cJ = cJ(3,3);
cF = real(invz_chi0z(sF, T, 0, struct()));  cF = cF(3,3);
verifyGreaterThan(testCase, cJ, 0);
r = cJ/cF;
verifyGreaterThan(testCase, r, 0.68);  verifyLessThan(testCase, r, 0.86);
end
