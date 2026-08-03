%INVZ_RUN_PHASE_DIAGRAM Projected 1/z phase boundary.
% Combines Bc(T) field cuts with Tc(B) temperature cuts. ODD is selected by
% one boolean knob.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs ------------------------------------------------------------------
Ts = linspace(0.05, 1.45, 20);    % fixed-temperature Bc cuts (K)
Bs = linspace(0.05, 5.0, 26);     % fixed-field Tc cuts (T)
useODD = true;
useParallel = true;
gridN = 16;                         % debug resolution; converge reported boundaries separately
dipoleBackend = 'ewald';           % non-ODD: 'ewald' | 'bruteforce'
dpRng = 30;                        % brute-force cutoff and retained ODD geometry cutoff
% -----------------------------------------------------------------------------

couplingOpts = struct('grid', [gridN gridN gridN], 'dpRng', dpRng, ...
    'cache', false, 'dipole', dipoleBackend, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], ...
    'gammaPolicy', 'P_drop');
if strcmp(dipoleBackend, 'ewald')
    couplingOpts.ewald = struct('alpha', 0.3, 'r_cut', 16, 'g_cut', 3, ...
        'boundary', 'conducting_k0_omitted');
end

if useODD
    g = invz_phase1_qgrid(ion, gridN, [0 0 0], 'halfopen', 'P_drop');
    [Vca, Vcb, Vcc, oddInfo] = invz_odd_blocks( ...
        ion, g.qvec, struct('dpRng', dpRng, 'cache', false));
    blocks = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, ...
        'Jcc0', oddInfo.Jcc0, 'Jaa0', oddInfo.Jaa0);
    solveOpts = struct('J0eff', oddInfo.Jcc0, 'Jxx0', oddInfo.Jaa0, ...
        'odd', true, 'odd_blocks', blocks);
    Jsolve = [];

    % The zero-field engine supplies the ODD anchor.
    Tc0 = invz_odd_zero_field(ion, struct( ...
        'mode', 'full', 'grids', {{gridN}}, 'blocks', blocks));
else
    [Jsolve, info, Jxx0] = invz_bz_couplings(ion, couplingOpts);
    solveOpts = struct('J0eff', info.Jcc0, 'Jxx0', Jxx0);
    J0 = info.Jcc0;
    Sc0 = invz_sigma_crit(J0, Jsolve);
    Tc0 = invz_critical_T0field(ion, Sc0, J0);
end

nT = numel(Ts);
nB = numel(Bs);
jobs = [Ts(:).' Bs(:).'];
kind = [ones(1,nT) 2*ones(1,nB)];   % 1=Bc(T), 2=Tc(B)

nWorkers = 0;
if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end
out = nan(size(jobs));
parfor (k = 1:numel(jobs), nWorkers)
    value = jobs(k);
    try
        if kind(k) == 1
            o = solveOpts;
            o.window = [0.1 6];
            out(k) = invz_critical(ion, value, Jsolve, o);
            fprintf('[%2d/%2d] T=%.3f K -> Bc=%.4f T\n', ...
                k, numel(jobs), value, out(k));
        else
            o = solveOpts;
            o.Tc0 = Tc0;
            out(k) = invz_critical_T(ion, value, Jsolve, o);
            fprintf('[%2d/%2d] B=%.3f T -> Tc=%.4f K\n', ...
                k, numel(jobs), value, out(k));
        end
    catch err
        fprintf('[%2d/%2d] %.3f failed: %s\n', ...
            k, numel(jobs), value, err.message);
    end
end

Bc = out(1:nT);
TcB = out(nT+1:end);
phase_boundary = sortrows([Ts(:) Bc(:); TcB(:) Bs(:)], 1);
phase_boundary = phase_boundary(all(isfinite(phase_boundary), 2), :);

fprintf('Zero-field Tc = %.6f K (ODD %s)\n', Tc0, onoff(useODD));
figure;
hold on;
plot(Ts, Bc, 'o-', 'DisplayName', 'B_c(T)');
plot(TcB, Bs, 's-', 'DisplayName', 'T_c(B)');
plot(Tc0, 0, 'ks', 'DisplayName', 'T_c(0)');
xlabel('T (K)');
ylabel('B_c (T)');
title('1/z phase boundary');
legend('Location', 'southwest');

function s = onoff(tf)
if tf, s = 'on'; else, s = 'off'; end
end
