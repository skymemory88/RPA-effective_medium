%INVZ_RUN_SPECTRA chi''_cc(omega) at the uniform mode vs field at T = 0.31 K (cf. R 2007 Fig 2,
% Kovacevic 2016 Fig 3d). RPA vs 1/z overlay.
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
T = 0.31;  w = (0.01:0.005:0.7).';
fields = [3.6 4.2 4.8 5.4 6.0];
figure; hold on;
for B = fields
    pt = invz_solve_point(ion, T, B, Jnu(:), struct('hyp', true, 'J0eff', info.Jcc0));
    out  = invz_chi_realaxis(ion, T, B, pt, w, struct('Jsel', info.Jcc0));
    pt0  = struct('alpha', 0, 'lambda', [0;0], 'tl', pt.tl, 'K', 0*pt.K);
    out0 = invz_chi_realaxis(ion, T, B, pt0, w, struct('Jsel', info.Jcc0, 'npass', 1));
    plot(w, imag(out.chi_cc_q(1,:)), '-', 'DisplayName', sprintf('1/z, %.1f T', B));
    plot(w, imag(out0.chi_cc_q(1,:)), '--', 'DisplayName', sprintf('RPA, %.1f T', B));
    fprintf('B = %.1f T : Sigma(0) = %.4f, 1+Sigma(0) = %.4f\n', B, pt.Sigma0, 1+pt.Sigma0);
end
xlabel('\omega (meV)'); ylabel('\chi''''_{cc}'); legend show;
