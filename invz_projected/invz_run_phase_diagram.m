%INVZ_RUN_PHASE_DIAGRAM Reproduce R 2007 Fig 1 (paramagnetic-side boundary).
% Two-regime search: Bc(T) at each fixed T in Ts (invz_critical, vertical field cuts), and
% Tc(B) at each fixed B in Bs (invz_critical_T, horizontal temperature cuts, self-adapting window).

addpath(fileparts(mfilename('fullpath')));  addpath(fullfile(fileparts(mfilename('fullpath')),'..'));
ion = invz_ion();
% ion.demag = 1;  ion.alpha = 0.25;  % OPTIONAL sample-shape (demagnetization) knob; default off.
%   Ordering-channel criticality (info.Jcc0/Jnu) and Tc(B=0) are demag-INVARIANT (R2007); Bc(T)
%   vs APPLIED field can shift via the demag-aware info.Jaa0 (hoisted into Jxx0 below). See
%   invz_ion.m for the full explanation. demag = 0 (default) matches the R2007 benchmark.
% ion.odd = 1;  % OPTIONAL ODD (off-diagonal dipolar) Tier-1 knob (odd_implementation_plan.html,
%   T1.4): geometric blocks are built ONCE pre-parfor below (P0.4: workers do no disk I/O); the
%   point solvers rebuild the cc modes with deltaJ(T,B) and apply the E5 -d to J0eff internally
%   (callers keep passing the UNSHIFTED J0); the Tc(B) anchor Tc0 becomes odd-aware via
%   T-dependent Sc(T)/J0(T) handles. Intrinsic-only (requires demag = 0). Default off = the
%   published benchmark run, byte-identical.
[Jf, info, Jxx0] = invz_bz_couplings(ion);   % shared BZ-grid coupling branches (Jaa0-aware)
J0 = info.Jcc0;   % scalar hoist: avoids broadcasting the whole info struct to workers
% Jxx0 is demag-aware; at demag = 0 it differs from the hardcoded ion.Jxx0 by <0.1% -- the live
% dipole sum supersedes the pasted constant. Tc0 below needs no Jxx0: at B = 0, <Jx> = 0.
% Zero-field Tc, computed ONCE up front (invz_sigma_crit warns once here rather
% than in every worker): it anchors the Tc(B) adaptive window and is the B=0
% endpoint on the plot.
Sodd = struct();                            % ODD blocks (empty when ion.odd is off; parfor-safe)
if ion.odd
    % ---- ODD (T1.4): geometric blocks ONCE, pre-parfor, on the SAME 16^3 Gamma-less
    % grid invz_bz_couplings just used (its defaults: grid [16 16 16], dpRng 30, cache).
    [qodd, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', [16 16 16], 'range', [-0.5 0.5], 'verbose', false);
    qodd = qodd(any(abs(qodd) > 1e-12, 2), :);
    [VcaO, VcbO, VccO, infoO] = invz_odd_blocks(ion, qodd, struct('dpRng', 30, 'cache', true));
    Sodd = struct('Vca', VcaO, 'Vcb', VcbO, 'Vcc', VccO, 'Jcc0', infoO.Jcc0);   % Jcc0 UNSHIFTED
    % Odd-aware zero-field anchor -- T1.4 SEAM: T-dependent Sc(T)/J0(T) handles through the
    % generalized invz_critical_T0field (same algebra as invz_odd_zero_field mode 'full';
    % replace this block with a call there once T1.5 lands). invz_critical_T refuses to
    % anchor adaptively with opts.odd on (invz:oddTc0), so Tc0 MUST be odd-aware here.
    J0T = @(T) J0 - odd_d_at(ion, Sodd, Jxx0, T);
    ScT = @(T) odd_Sc_at(ion, Sodd, J0, Jxx0, T);
    ScT(2.0);   % probe: surfaces the known 16^3 invz:sigmaCritExcluded warning ONCE (ODD-LOG T1.3:
                % a few near-Gamma transverse-shell modes exceed J0 already without ODD)
    wsOdd = warning('off', 'invz:sigmaCritExcluded');   % silence the same known warning inside the bisection
    Tc0 = invz_critical_T0field(ion, ScT, J0T);
    warning(wsOdd);
else
Tc0 = invz_critical_T0field(ion, invz_sigma_crit(J0, Jf), J0);
end

% ---- knobs ----------------------------------------------------------------
% Tc(B) search window is not set here: invz_critical_T self-adapts per field
% (anchors at Tc0+0.05 K, spans 0.5 K down; see its header for details).

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
            if ion.odd
                % ODD (T1.4): Jnu_flat = [] (modes are rebuilt in the solver from the
                % blocks + deltaJ(T,B)); J0 stays UNSHIFTED (solver applies the E5 -d).
                val = invz_critical(ion, v, [], struct('J0eff', J0, 'Jxx0', Jxx0, ...
                    'window', [0.1 6], 'odd', true, 'odd_blocks', Sodd));
            else
            % Field window [0.1 6]: the top (6 T) is paramagnetic for every Ts
            % (Bc(T=0) ~ 5 T); invz_critical scans DOWN from there to the
            % converged ordered/paramagnet crossing (see its header).
            val = invz_critical(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'window', [0.1 6]));
            end
        catch err
            fprintf('  T = %.2f K: FAILED (%s)\n', v, err.message);
        end
        fprintf('  [%2d/%2d] T = %.2f K  ->  Bc = %.3f T   (%.0f s)\n', k, nT+nB, v, val, toc(t0));
    else
        try
            if ion.odd
                % ODD (T1.4): Tc0 here is the odd-aware anchor computed above
                % (invz_critical_T errors invz:oddTc0 without it).
                val = invz_critical_T(ion, v, [], struct('J0eff', J0, 'Jxx0', Jxx0, ...
                    'Tc0', Tc0, 'odd', true, 'odd_blocks', Sodd));
            else
            val = invz_critical_T(ion, v, Jf, struct('J0eff', J0, 'Jxx0', Jxx0, 'Tc0', Tc0));
            end
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

% ---------------------------------------------------------------------------
% ODD anchor helpers (ion.odd = 1 only; T1.4 seam, superseded by
% invz_odd_zero_field once T1.5 lands). Script-local functions: only reached
% via the handles built in the ion.odd branch above.
function d = odd_d_at(ion, S, Jxx0, T)
%ODD_D_AT E5 uniform reduction d(T) at B = 0 (chi_perp-mediated).
Xp = invz_chiperp(ion, T, [0 0 0], struct('Jxx0', Jxx0));
[~, d] = invz_odd_deltaJ(S.Vca, S.Vcb, Xp);
end

function Sc = odd_Sc_at(ion, S, J0, Jxx0, T)
%ODD_SC_AT Critical self-energy Sc(T) on the ODD-rebuilt modes with the shifted J(0)
% (invz_odd_zero_field mode 'full' algebra: modes of Vcc + deltaJ, J0(T) = J0 - d(T)).
Xp = invz_chiperp(ion, T, [0 0 0], struct('Jxx0', Jxx0));
[dJ, d] = invz_odd_deltaJ(S.Vca, S.Vcb, Xp);
nq = size(S.Vcc, 3);
Jnu = zeros(nq, 4);
for iq = 1:nq
    M = S.Vcc(:,:,iq) + dJ(:,:,iq);
    M = (M + M')/2;                            % both terms Hermitian; cleans rounding only
    Jnu(iq,:) = sort(real(eig(M))).';
end
Sc = invz_sigma_crit(J0 - d, Jnu(:));
end
