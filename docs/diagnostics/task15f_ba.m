% Task-15 final fix pass -- BEFORE/AFTER on one process each.
% Usage: SIDE = 'head' shadows 73dc55c's invz_spectra_map.m; SIDE = 'work' uses the worktree.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
SIDE = getenv('SIDE');
if strcmp(SIDE, 'head'), addpath(fullfile(ROOT, '.superpowers', 'sdd', 'shadow_head'), '-begin'); end
fprintf('##### SIDE = %s (which invz_spectra_map: %s) #####\n', SIDE, which('invz_spectra_map'));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @(s) struct('Jnu', s*linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', s*6.4e-3), ...
                  'verbose', false);
j   = @(c) strjoin(c, ' | ');

% (1) the weak-coupling 'solver_failed' sweep -- the undercount case
o = syn(0.3);  o.static_medium = 'strict_1z_dyson_ref';
lastwarn('', 'invz:none');
S = invz_spectra_map(ion, 0.31, [0 2.85 3.30 5.50], w, o);
[wm, wi] = lastwarn();
fprintf('BA1 reason  = %s\n', j(S.phase_1z_reason));
fprintf('BA1 n_msk=%d Bc_1z=%g (%s)\n', nnz(S.phase_1z == 0), S.Bc_1z, S.Bc_1z_status);
fprintf('BA1 WARNID  = %s\n', wi);
fprintf('BA1 WARNMSG = %s\n', wm);

% (1b) THE SILENT CASE: only-uncounted mask + a finite Bc_1z bracket
o = syn(0.40);  o.static_medium = 'strict_1z_dyson_ref';
Fs = [0 1.5 2.0 5.50];
lastwarn('', 'invz:none');
S = invz_spectra_map(ion, 0.31, Fs, w, o);
[wm, wi] = lastwarn();
fprintf('\nBA2 reason  = %s\n', j(S.phase_1z_reason));
fprintf('BA2 n_msk=%d Bc_1z=%g (%s)\n', nnz(S.phase_1z == 0), S.Bc_1z, S.Bc_1z_status);
fprintf('BA2 WARNID  = %s\n', wi);
fprintf('BA2 WARNMSG = %s\n', wm);

% (2) hunt: uncounted mask + FINITE Bc_1z == fully silent before the fix
fprintf('\n### hunt: uncounted masked reason together with a finite Bc_1z ###\n');
UNC = {'unstable_endpoint','solver_failed','not_attempted_longitudinal', ...
       'bare_not_ordered','response_failed'};
for s = [0.40 0.45 0.50 0.55]
    o = syn(s);  o.static_medium = 'strict_1z_dyson_ref';
    F = [0 0.5 1.0 1.5 2.0 2.85 5.50];
    lastwarn('', 'invz:none');
    ws = warning('off', 'invz:spectraMapMasked');
    try
        S = invz_spectra_map(ion, 0.31, F, w, o);
    catch ME
        warning(ws);  fprintf('s=%.2f THREW %s\n', s, ME.identifier);  continue;
    end
    warning(ws);
    [~, wid] = lastwarn();
    nmsk = nnz(S.phase_1z == 0);
    hit  = isfinite(S.Bc_1z) && any(ismember(S.phase_1z_reason, UNC));
    fprintf('s=%.2f Bc=%-8.4g %-12s nmsk=%d WARN=%-24s HIT=%d  %s\n', s, S.Bc_1z, ...
            S.Bc_1z_status, nmsk, wid, hit, j(S.phase_1z_reason));
end
fprintf('DONE_BA_%s\n', SIDE);
