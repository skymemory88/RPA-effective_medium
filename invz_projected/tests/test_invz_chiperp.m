function tests = test_invz_chiperp
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_anchors_and_symmetry(testCase)
ion = invz_ion();  A = invz_odd_anchors();
[Xp, info] = invz_chiperp(ion, 1.53, [0 0 0], struct());
% AbsTol floor for the nominally-zero off-diagonals (machine noise on a nominal zero can
% never satisfy RelTol); the ~17.6 diagonals stay governed by RelTol (1e-9*17.6 >> 1e-12).
verifyEqual(testCase, Xp, A.chiperp_1p53K_0T, 'RelTol', 1e-9, 'AbsTol', 1e-12);
verifyEqual(testCase, Xp(1,1), Xp(2,2), 'AbsTol', 1e-10*abs(Xp(1,1)));   % C4 at Bx = 0
verifyEqual(testCase, Xp, Xp.', 'AbsTol', 1e-15);
verifyLessThan(testCase, info.asym, 1e-8*abs(Xp(1,1)));
verifyGreaterThan(testCase, Xp(1,1), 10);   % the 16-17 meV^-1 band, floor guard
verifyLessThan(testCase, Xp(1,1), 25);
end

function test_zero_field_no_doublet_guard(testCase)
% Regular at Bx = 0 (never routes through invz_twolevel).
ion = invz_ion();
Xp = invz_chiperp(ion, 0.31, [0 0 0], struct());
verifyTrue(testCase, all(isfinite(Xp(:))));
end

function test_reproducible_along_Bx(testCase)
% P0 AMENDMENT (ODD-LOG SSP0.2): the 0.31 K chi_aa(Bx) curve has a PHYSICAL
% peak near 1 T (halves by 2 T; all points MF-converged), so a step-size
% "smoothness" gate is wrong physics. The anchor-pinned sweep IS the
% no-numerical-artifact guard.
ion = invz_ion();  A = invz_odd_anchors();
x = zeros(1, 7);
for i = 0:6, Xp = invz_chiperp(ion, 0.31, [i 0 0], struct()); x(i+1) = Xp(1,1); end
verifyEqual(testCase, x(:), A.chiperp_0p31K_Bx.chi_aa(:), 'RelTol', 1e-9);   % pinned P0 sweep
verifyTrue(testCase, all(isfinite(x)) && all(x > 0));
end

function test_matsubara_frequencies(testCase)
% z-vector form: real symmetric per slice, decaying with |w_n|, r_n in (0, 1].
ion = invz_ion();
[wn, ~, ~] = invz_matsubara(1.53, 10);
[Xp, ~] = invz_chiperp(ion, 1.53, [0 0 0], struct('z', 1i*wn));
verifyEqual(testCase, size(Xp), [2 2 numel(wn)]);
r = squeeze(Xp(1,1,:)) / Xp(1,1,1);
verifyEqual(testCase, r(1), 1, 'AbsTol', 1e-14);
verifyTrue(testCase, all(diff(r) < 0) && all(r > 0));   % monotone decay along iw_n
% rough scale from the plan SS4: r(w1 = 0.828 meV) ~ eps1^2/(eps1^2 + w1^2) ~ 0.56
end
