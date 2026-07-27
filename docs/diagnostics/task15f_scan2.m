% Task-15 final fix pass -- I2 hunt round 2: an UNCOUNTED masked reason.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @() struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', 6.4e-3), ...
                 'verbose', false);
j = @(c) strjoin(c, ' | ');

%% ---- complexity narrowing, explained ------------------------------------------------------
fprintf('########## complexity narrowing through the driver prologue ##########\n');
fq = complex([3 5], [0 0]);
fprintf('isreal(fq)=%d  isreal(fq(:).'')=%d  isreal(fq(:)*[1 0 0])=%d\n', ...
        isreal(fq), isreal(fq(:).'), isreal(fq(:)*[1 0 0]));
fn = [3+1i 5+2i];
fprintf('isreal(fn)=%d  isreal(fn(:).'')=%d  isreal(fn(:)*[1 0 0])=%d\n', ...
        isreal(fn), isreal(fn(:).'), isreal(fn(:)*[1 0 0]));

%% ---- I2a: tilt sweep straddling the bare boundary -> not_attempted_longitudinal -----------
fprintf('\n########## I2a: tilt, mixed bare_escape_hatch + not_attempted_longitudinal ##########\n');
o = syn();  o.static_medium = 'strict_1z_dyson_ref';  o.field_dir = [cosd(3) 0 sind(3)];
lastwarn('', 'invz:none');
ws = warning('off', 'invz:spectraMapMasked');
S = invz_spectra_map(ion, 0.31, [2.85 5.50], w, o);
warning(ws);
[wm, wi] = lastwarn();
fprintf('I2a reason = %s\n', j(S.phase_1z_reason));
fprintf('I2a ph1z   = [%s]   n_msk = %d\n', num2str(S.phase_1z), nnz(S.phase_1z == 0));
fprintf('I2a Bc_1z  = %g (%s)\n', S.Bc_1z, S.Bc_1z_status);
fprintf('I2a WARNID %s\nI2a WARNMSG %s\n', wi, wm);

%% ---- I2b: tilt sweep, ALL columns above the bare boundary -> no bare escape hatch --------
fprintf('\n########## I2b: tilt, PM-side only -> masked with NO bare escape hatch ##########\n');
o = syn();  o.static_medium = 'strict_1z_dyson_ref';  o.field_dir = [cosd(3) 0 sind(3)];
lastwarn('', 'invz:none');
ws = warning('off', 'invz:spectraMapMasked');
S = invz_spectra_map(ion, 0.31, [5.50 8.00], w, o);
warning(ws);
[wm, wi] = lastwarn();
fprintf('I2b reason = %s\n', j(S.phase_1z_reason));
fprintf('I2b ph1z   = [%s]   n_msk = %d   n_bare = %d\n', num2str(S.phase_1z), ...
        nnz(S.phase_1z == 0), nnz(strcmp(S.phase_1z_reason, 'bare_escape_hatch')));
fprintf('I2b Bc_1z  = %g (%s)\n', S.Bc_1z, S.Bc_1z_status);
fprintf('I2b WARNID %s\nI2b WARNMSG %s\n', wi, wm);

%% ---- I2c: real 4^3 couplings under strict, hunting solver_failed/unstable_endpoint -------
fprintf('\n########## I2c: real 4^3 couplings, strict ##########\n');
Fs = [0 1.0 2.0 3.0 3.5 4.0 5.0 8.0];
o = struct('grid', [4 4 4], 'dpRng', 12, 'cache', false, 'verbose', false, ...
           'static_medium', 'strict_1z_dyson_ref', 'ordered_1z', 'jensen');
lastwarn('', 'invz:none');
ws = warning('off', 'invz:spectraMapMasked');
tic;  S = invz_spectra_map(ion, 0.31, Fs, w, o);  el = toc;
warning(ws);
[wm, wi] = lastwarn();
fprintf('I2c elapsed %.1f s   Bc_1z=%g (%s)\n', el, S.Bc_1z, S.Bc_1z_status);
for k = 1:numel(Fs)
    fprintf('   B=%5.2f ph1z=%d reason=%-28s diag=%-22s stab=%d crit=%.4g\n', Fs(k), ...
        S.phase_1z(k), S.phase_1z_reason{k}, S.ordered_diag_reason{k}, ...
        S.stability_1z(k), S.crit_pm(k));
end
fprintf('I2c WARNID %s\nI2c WARNMSG %s\n', wi, wm);
fprintf('DONE_SCAN2\n');
