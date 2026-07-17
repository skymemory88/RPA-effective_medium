function tests = test_invz_chi_observable
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
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
% R 2007 Fig 2: at T=0.31 K near Hc the lowest mode bottoms at ~0.19 meV (calc),
% never zero. Uses a fixed field Bx = 4.3 T (published Hc(0.31 K) ~ 4.24-4.3 T)
% instead of bisecting invz_critical (too slow at full grid cost); the wide
% [0.10, 0.28] meV band absorbs the resulting mistuning. SLOW.
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

function test_chi0cc_passthrough_is_consumed(testCase)
% invz_chi_realaxis reuses a precomputed bare chi0cc(w) when given, so one field
% point need not rebuild the single ion / chi0 for both its RPA and 1/z calls
% (finding #2). A correct precomputed value reproduces the internal result exactly;
% a perturbed one changes the output, proving the option is actually consumed.
ion = invz_ion();
T = 1.0;  Bx = 3.0;  w = (0.05:0.02:0.7).';  eta = 5e-3;
pt = struct('alpha', 0.1, 'lambda', [0.2; 0.1], 'tl', invz_twolevel(ion, T, Bx), 'K', []);
copts = struct('hyp', false, 'eta', eta, 'npass', 1, 'Jsel', [5e-3; 6.4e-3]);
o_ref = invz_chi_realaxis(ion, T, Bx, pt, w, copts);
si  = invz_single_ion(ion, T, [Bx 0 0], struct('hyp', false, 'Jxx0', ion.Jxx0));
c0  = invz_chi0z(si, T, w + 1i*eta, struct('elastic', false));
chi0cc = squeeze(c0(3,3,:));
copts2 = copts;  copts2.chi0cc_w = chi0cc;
o_new = invz_chi_realaxis(ion, T, Bx, pt, w, copts2);
verifyEqual(testCase, o_new.chi_cc_q, o_ref.chi_cc_q, 'AbsTol', 1e-12);
copts3 = copts;  copts3.chi0cc_w = 1.5*chi0cc;
o_pert = invz_chi_realaxis(ion, T, Bx, pt, w, copts3);
verifyGreaterThan(testCase, max(abs(o_pert.chi_cc_q(:) - o_ref.chi_cc_q(:))), 1e-6);
end

function test_demag_observable_rescale(testCase)
% The measured strict-uniform response chi_meas = chi_int/(1 + Jshape*chi_int) must
% equal the OLD-convention shifted pole chit/(1 - (Jcc0 - Jshape)*chit) exactly
% (algebraic identity), and Jshape = 0 must be a byte-identical no-op.
ion = invz_ion();
T = 0.31;  Bx = 5.0;  w = (0:0.01:0.4).';
tl0 = invz_twolevel(ion, T, Bx);
pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
Jcc0 = ion.J0eff;  Jsh = 1e-3;
oi = invz_chi_realaxis(ion, T, Bx, pt0, w, struct('Jsel', Jcc0, 'npass', 1));
om = invz_chi_realaxis(ion, T, Bx, pt0, w, struct('Jsel', Jcc0, 'npass', 1, 'Jshape', Jsh));
chit  = oi.chi0cc_w(:).' ./ (1 + oi.Sigma_w(:).');
expct = chit ./ (1 - (Jcc0 - Jsh)*chit);
verifyEqual(testCase, om.chi_cc_q(1,:), expct, 'RelTol', 1e-10);
verifyEqual(testCase, oi.Jshape, 0);
verifyEqual(testCase, om.Jshape, Jsh);
end
