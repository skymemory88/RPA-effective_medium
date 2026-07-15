%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
%
% Two-regime search:
%   Bc search:  find the critical field Bc(T) at each fixed T in Ts
%           (invz_critical, vertical cuts);
%   Tc search: find the critical temperature Tc(B) at each fixed B in Bs
%           (invz_critical_T, horizontal cuts, self-adapting window).

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
[qvec, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
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
% crossing (see invz_critical_T header). 

Ts = linspace(0.05, 1.95, 28);   % low-T regime: Bc(T) points
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
