%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
%
% Two-regime search:
%   low-T:  bisect the critical field Bc(T) at each fixed T in Ts
%           (invz_critical, vertical cuts);
%   high-T: bisect the critical temperature Tc(B) at each fixed B in Bs
%           (invz_critical_T, horizontal cuts, window [Tlo Tmax] below).
% Near the classical critical point (B -> 0, T -> the zero-field Tc) the
% boundary is nearly parallel to the field axis, so a vertical cut crosses it
% at a glancing angle: brackets fail and tiny T errors give huge Bc errors.
% A horizontal cut crosses it transversally and is well-conditioned there.
% That is why the default Ts stops at 1.6 K -- vertical cuts above ~1.7 K
% have no bracket; the Bs list owns that part of the boundary.
%
% Bracket geometry: at T = Tlo a point is ordered exactly when B < Bc(Tlo)
% (~2.8 T at 1.0 K), so every Bs entry below that brackets cleanly in
% [Tlo, Tmax]; an entry above it fails its bracket assert to NaN without
% affecting other jobs.
%
% Bs has a 0.5 T floor: at 0.2-0.3 T the paramagnetic solve develops
% non-convergence patches near the boundary; the ordered/para classifier
% reads them as ordered and biases Tc upward by ~0.04-0.05 K (and at
% Bx ~ 0 invz_twolevel raises invz:degenerateDoublet outright). The
% 0 < B < 0.5 T boundary segment spans only ~4 mK below the zero-field Tc
% and is represented by the closed-form Tc0 endpoint on the plot.
%
% Parallelism: all nT+nB boundary points are INDEPENDENT 1-D bisections, so a
% single flat parfor covers both regimes -- near-linear speedup up to the job
% count. The bisection *inside* one point cannot be parallelized (each bracket
% halving needs the sign of the previous EMT solve). parfor degrades to a
% serial loop automatically when the Parallel Computing Toolbox is absent
% (nWorkers = 0); with it, a local pool is auto-created on first use
% (~10-30 s once). The Jq lattice sum is computed ONCE up front, so workers
% do no disk I/O and never touch the invz/cache -- no contention. Progress
% lines interleave across workers when running in parallel; that is expected.

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Jf = Jnu(:);
J0 = info.Jcc0;   % scalar hoist: avoids broadcasting the whole info struct to workers

% ---- knobs --------------------------------------------------------------
Ts = [0.05 0.2 0.3 0.4 0.6 0.8 1.0 1.2 1.4 1.6];   % low-T regime: Bc(T) points
Bs = [0.5 0.75 1.0 1.25 1.5];   % high-T regime: Tc(B) points (0.5 T floor, header)
useParallel = true;             % false -> force a serial run
% -------------------------------------------------------------------------
% Tc(B) bisection window: Tlo must be ordered for every Bs entry
% (max(Bs) < Bc(Tlo) ~ 2.8 T at 1.0 K); Tmax must be paramagnetic
% (anything above the zero-field boundary, ~1.78 K on this grid).
Tlo = 1.0;  Tmax = 2.0;

nT = numel(Ts);  nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];           % one independent bisection per job
kind = [ones(1,nT) 2*ones(1,nB)];   % 1 = Bc(T) at fixed T, 2 = Tc(B) at fixed B

nWorkers = 0;                       % 0 = serial (also the no-toolbox fallback)
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

out = nan(1, nT + nB);
parfor (k = 1:nT+nB, nWorkers)
    t0 = tic;
    v = jobs(k);  val = NaN;
    if kind(k) == 1
        try
            % Wide bracket: the low-T points need the upper edge to reach
            % ~4-5 T, while the 1.6 K point sits close to the lower edge.
            val = invz_critical(ion, v, Jf, struct('J0eff', J0, 'window', [0.5 7]));
        catch err
            fprintf('  T = %.2f K: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    else
        try
            val = invz_critical_T(ion, v, Jf, struct('J0eff', J0, 'window', [Tlo Tmax]));
        catch err
            fprintf('  B = %.2f T: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] B = %.2f T  ->  Tc = %.3f K   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    end
    out(k) = val;
end
Bc  = out(1:nT);                    % low-T branch:  Bc at each Ts
TcB = out(nT+1:end);                % high-T branch: Tc at each Bs

Tc0 = invz_critical_T0field(ion, invz_sigma_crit(J0, Jf), J0);
fprintf('Zero-field Tc (1/z) = %.3f K  [target 1.74 K]\n', Tc0);

% Merged boundary, T-sorted, finite points only -- workspace export for
% downstream use ('boundary' would shadow the built-in of that name).
% Columns [T(K) B(Tesla)] (the plot converts to kOe). Near the regime join
% (~1.5-1.6 K) both branches contribute a point, so the curve is not
% strictly single-valued in T there.
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);

figure; hold on;
plot(Ts, Bc*10, 'o-');
plot(TcB, Bs*10, 's-');
plot(Tc0, 0, 'ks');
xlabel('T (K)'); ylabel('B_c (kOe)');
title('LiHoF_4 1/z phase boundary (paramagnetic side)');
legend({'B_c(T) bisection', 'T_c(B) bisection', 'closed-form T_c(B=0)'}, 'Location', 'southwest');
