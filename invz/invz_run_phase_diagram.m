%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
%
% NOTE: this is an interactive driver, not a test. Each invz_critical(...) call
% bisects over full invz_solve_point EMT solves at the 16^3-grid resolution;
% a full run over all seven temperatures below takes roughly 1-2 hours on a
% single core. qVec_generator prints its lattice diagnostics to stdout on each
% call below, which is expected/acceptable here (unlike in the test suite,
% where it is wrapped in evalc to keep test output clean).
addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Ts = [0.15 0.31 0.6 0.9 1.2 1.4 1.6];
Bc = nan(size(Ts));
for k = 1:numel(Ts)
    fprintf('[%d/%d] T = %.2f K : solving for Bc ...\n', k, numel(Ts), Ts(k));
    tic;
    try
        % Wide bracket: high-T points near Tc need a wider window to bracket
        % the boundary (e.g. Hc(1.6 K) ~ 1.5-2 T, close to the default lower
        % edge), while low-T points need the upper edge to reach ~4-5 T.
        Bc(k) = invz_critical(ion, Ts(k), Jnu(:), struct('J0eff', info.Jcc0, 'window', [0.5 7]));
    catch err
        warning('T=%.2f K: %s', Ts(k), err.message);
    end
    fprintf('T = %.2f K  ->  Bc = %.2f T   (%.1f s)\n', Ts(k), Bc(k), toc);
end
Tc0 = invz_critical_T0field(ion, invz_sigma_crit(info.Jcc0, Jnu(:)), info.Jcc0);
fprintf('Zero-field Tc (1/z) = %.3f K  [target 1.74 K]\n', Tc0);
figure; plot(Ts, Bc*10, 'o-'); hold on; plot(Tc0, 0, 'ks');
xlabel('T (K)'); ylabel('B_c (kOe)'); title('LiHoF_4 1/z phase boundary (paramagnetic side)');
