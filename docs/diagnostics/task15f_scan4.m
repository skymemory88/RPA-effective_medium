% Task-15 final fix pass -- I2 hunt round 4: 'solver_failed' at B = 0 (uncounted) TOGETHER
% with an ordered/pm bracket, i.e. a finite Bc_1z -> the pre-fix sweep must be fully SILENT.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @(s) struct('Jnu', s*linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', s*6.4e-3), ...
                  'verbose', false);
j   = @(c) strjoin(c, ' | ');

for s = [0.30 0.35 0.40 0.45 0.50 0.55]
  for T = [0.31 0.10 0.05]
    o = syn(s);  o.static_medium = 'strict_1z_dyson_ref';
    F = [0 0.25 0.50 1.00 1.50 2.85 5.50];
    lastwarn('', 'invz:none');
    ws = warning('off', 'invz:spectraMapMasked');
    try
        S = invz_spectra_map(ion, T, F, w, o);
    catch ME
        warning(ws);  fprintf('s=%g T=%g THREW %s\n', s, T, ME.identifier);  continue;
    end
    warning(ws);
    [~, wid] = lastwarn();
    fprintf('s=%.2f T=%.2f Bc=%-8.4g %-12s WARN=%-24s reason = %s\n', s, T, S.Bc_1z, ...
            S.Bc_1z_status, wid, j(S.phase_1z_reason));
    if isfinite(S.Bc_1z) && any(strcmp(S.phase_1z_reason, 'solver_failed'))
        fprintf('   *** JACKPOT (silent uncounted mask + finite Bc_1z) ***\n');
    end
  end
end
fprintf('DONE_SCAN4\n');
