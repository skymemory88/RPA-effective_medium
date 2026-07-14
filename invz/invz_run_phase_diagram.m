%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
%
% Two-regime search:
%   low-T:  find the critical field Bc(T) at each fixed T in Ts
%           (invz_critical, vertical cuts);
%   high-T: find the critical temperature Tc(B) at each fixed B in Bs
%           (invz_critical_T, horizontal cuts, self-adapting window).
% Near the classical critical point (B -> 0, T -> the zero-field Tc) the
% boundary is nearly parallel to the field axis, so a vertical cut crosses it
% at a glancing angle: it becomes ill-conditioned and tiny T errors give large
% Bc errors. A horizontal cut crosses it transversally and is well-conditioned
% there. That is why the low-T Ts list is best kept below ~1.6 K and the Bs
% list owns the part of the boundary just under Tc0.
%
% Robustness (both cuts): within a narrow band around the boundary the
% paramagnetic self-consistency suffers critical slowing down (the outer loop
% does not reach tolerance; the soft-mode Dyson denominator blows up to NaN).
% Both invz_critical_T (samples crit on a T grid) and invz_critical (scans the
% field down from the paramagnet) classify from CONVERGED points ONLY and
% interpolate the crossing, so those failures no longer masquerade as "ordered"
% and scatter the boundary. (The previous bisection counted them as ordered and
% produced a rugged boundary, some Tc(B) points even above Tc0.) invz_critical_T
% additionally raises invz:multipleCrossings on more than one converged crossing
% -- a candidate hyperfine re-entrant nose worth inspecting rather than hiding.
%
% Small B: below ~0.5 T the doublet is near-degenerate (invz_twolevel raises
% invz:degenerateDoublet as Bx -> 0) so few points converge; the 0 < B < 0.5 T
% segment spans only ~4 mK below Tc0 and is best read from the closed-form Tc0
% endpoint on the plot.
%
% Parallelism: all nT+nB boundary points are INDEPENDENT 1-D root finds, so a
% single flat parfor covers both regimes -- near-linear speedup up to the job
% count. Each point is solved on its own worker; the per-point root find (a
% serial scan/grid of EMT solves plus interpolation) is left serial within a
% point, since the job-level parfor already saturates the pool. parfor degrades to a
% serial loop automatically when the Parallel Computing Toolbox is absent
% (nWorkers = 0); with it, a local pool is auto-created on first use
% (~10-30 s once). The Jq lattice sum is computed ONCE up front, so workers
% do no disk I/O and never touch the invz/cache -- no contention. Progress
% lines interleave across workers when running in parallel; that is expected.

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
% ion.demag = 1;  ion.alpha = 0.25;  % OPTIONAL sample-shape (demagnetization) knob; default off.
%   - info.Jcc0, Jnu, and the ordering-channel contribution to criticality are demag-
%     INVARIANT (R 2007: the demagnetizing field cancels from the critical condition,
%     since ordering occurs at q -> 0+, not strict q = 0).
%   - Tc(B = 0) is EXACTLY demag-invariant: the transverse moment vanishes there.
%   - Bc(T) vs APPLIED transverse field CAN shift through the demag-aware transverse
%     coupling info.Jaa0 (hoisted into Jxx0 below): internal-vs-applied field relation.
%   demag = 0 (default) is the intrinsic / internal-field boundary matching the R 2007
%   benchmark.
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5]);
qvec = qvec(any(abs(qvec) > 1e-12, 2), :);
[Jnu, info] = invz_jq_modes(ion, qvec, struct('dpRng', 30, 'cache', true));
Jf = Jnu(:);
J0 = info.Jcc0;   % scalar hoist: avoids broadcasting the whole info struct to workers
Jxx0 = ion.Jxx0;  if isfield(info, 'Jaa0'), Jxx0 = info.Jaa0; end   % live transverse J(0)
% (demag-aware; at demag = 0 it differs from the hardcoded ion.Jxx0 by <0.1% -- the live
% dipole sum supersedes the pasted constant). Tc0 below needs no Jxx0: at B = 0, <Jx> = 0.
% Zero-field Tc, computed ONCE up front (invz_sigma_crit warns once here rather
% than in every worker): it anchors the Tc(B) adaptive window and is the B=0
% endpoint on the plot.
Tc0 = invz_critical_T0field(ion, invz_sigma_crit(J0, Jf), J0);

% ---- knobs --------------------------------------------------------------
% -------------------------------------------------------------------------
% Tc(B) search window: invz_critical_T now self-adapts per field -- it anchors
% the window top at Tc0+0.05 K, spans 0.5 K down, samples crit on a grid,
% classifies from CONVERGED points only, and interpolates the highest-T
% crossing (see invz_critical_T header). No fixed [Tlo Tmax] bracket is needed
% and the old non-convergence ruggedness is gone. To force an explicit window
% instead, pass 'window',[Tlo Tmax] in the opts struct below.
Ts = linspace(0.05, 1.85, 26);   % low-T regime: Bc(T) points
Bs = [];   % high-T regime: Tc(B) points

% Ts = [];
% Bs = linspace(0.05,5,26);

useParallel = true;             % false -> force a serial run

nT = numel(Ts);  nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];           % one independent 1-D root find per job
kind = [ones(1,nT) 2*ones(1,nB)];   % 1 = Bc(T) at fixed T, 2 = Tc(B) at fixed B

nWorkers = 0;                       % 0 = serial (also the no-toolbox fallback)
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end

out = nan(1, nT + nB);
parfor (k = 1:nT+nB, nWorkers)
    t0 = tic;
    v = jobs(k);  val = NaN;
    if kind(k) == 1
        try
            % Field window [0.1 6]: the top (6 T) is paramagnetic for every Ts
            % (Bc(T=0) ~ 5 T); invz_critical scans DOWN from there to the
            % converged ordered/paramagnet crossing (see its header).
            val = invz_critical(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'window', [0.1 6]));
        catch err
            fprintf('  T = %.2f K: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    else
        try
            val = invz_critical_T(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'Tc0', Tc0));
        catch err
            fprintf('  B = %.2f T: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] B = %.2f T  ->  Tc = %.3f K   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    end
    out(k) = val;
end
Bc  = out(1:nT);                    % low-T branch:  Bc at each Ts
TcB = out(nT+1:end);                % high-T branch: Tc at each Bs

fprintf('Zero-field Tc (1/z) = %.3f K  [target 1.74 K]\n', Tc0);   % computed up front

% Merged boundary, T-sorted, finite points only -- workspace export for
% downstream use ('boundary' would shadow the built-in of that name).
% Columns [T(K) B(Tesla)]; the plot uses tesla directly. Near the regime join
% (~1.5-1.6 K) both branches contribute a point, so the curve is not
% strictly single-valued in T there.
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);

figure; hold on;
plot(Ts, Bc, 'o-');
plot(TcB, Bs, 's-');
plot(Tc0, 0, 'ks');
xlabel('T (K)'); ylabel('B_c (T)');
title('LiHoF_4 1/z phase boundary (paramagnetic side)');
legend({'B_c(T): fixed-T field cut', 'T_c(B): fixed-B temperature cut', 'closed-form T_c(B=0)'}, 'Location', 'southwest');
