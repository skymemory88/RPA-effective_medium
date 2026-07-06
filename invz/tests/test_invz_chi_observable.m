function tests = test_invz_chi_observable
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function test_rpa_limit_and_positivity(testCase)
% With Sigma forced to zero the output must equal scalar RPA built directly from chi0.
ion = invz_ion();
T = 1.0;  Bx = 4.0;  w = (0.02:0.02:0.8).';  eta = 5e-3;
Jsel = [5e-3; 6.4e-3];                       % two representative couplings (meV)
pt = struct('alpha', 0, 'lambda', [0;0], 'tl', invz_twolevel(ion, T, Bx));
out = invz_chi_realaxis(ion, T, Bx, pt, w, struct('Jsel', Jsel, 'eta', eta, 'npass', 1));
si  = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', true));
c0  = squeeze(invz_chi0z(si, T, w + 1i*eta, struct('elastic', false)));
c0cc = squeeze(c0(3,3,:));
rpa = c0cc ./ (1 - Jsel(2)*c0cc);
verifyEqual(testCase, out.chi_cc_q(2,:).', rpa, 'RelTol', 1e-6);
verifyGreaterThan(testCase, min(imag(out.chi_cc_q(2,:))), -1e-10);  % chi'' >= 0 for w > 0
end

function test_soft_mode_near_criticality(testCase)
% R 2007 Fig 2: at T=0.31 K, H~Hc the lowest mode bottoms at ~0.19 meV (calc), never zero. SLOW.
%
% Controller adaptation: the brief calls invz_critical(ion, 0.31, ...) to locate
% Hc before solving; that bisection re-solves invz_solve_point O(log2((7-2)/0.02))
% ~ 8 times at full 16^3-grid EMT cost and takes ~10-20 minutes. To keep this test
% bounded while still self-contained, we use the fixed field Bx = 4.3 T instead
% (published Hc(0.31 K) ~ 42.4-43 kOe = 4.24-4.3 T, R 2007), which is within the
% bisection tolerance of the true critical field. The assert's peak-position band
% [0.10, 0.28] meV absorbs the resulting small mistuning.
assumeTrue(testCase, strcmp(getenv('INVZ_SLOW'), '1'), 'Set INVZ_SLOW=1 for slow tests');
ion = invz_ion();
cmd = "qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5])";
[~, qvec] = evalc(cmd);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
T = 0.31;
bx = 4.3;   % controller-approved substitute for invz_critical(ion, T, Jnu(:), ...); see note above
pt = invz_solve_point(ion, T, bx, Jnu(:), struct('hyp', true, 'J0eff', info.Jcc0));
w  = (0.01:0.005:0.6).';
out = invz_chi_realaxis(ion, T, bx, pt, w, struct('Jsel', info.Jcc0, 'eta', 5e-3));
[~, ipk] = max(imag(out.chi_cc_q(1,:)));
Epk = w(ipk);
verifyGreaterThan(testCase, Epk, 0.10);  verifyLessThan(testCase, Epk, 0.28);
end
