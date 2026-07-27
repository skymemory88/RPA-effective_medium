% Task-15 final fix pass -- scan for (a) an UNCOUNTED masked reason with a finite Bc_1z (I2),
% (b) whether genuinely-complex `fields` really reach the line-348 screen (M1).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @() struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', 6.4e-3), ...
                 'verbose', false);

%% ---- M1b: complex variants, where exactly do they die? -----------------------------------
fprintf('########## M1b: complex `fields` variants ##########\n');
cv = {'complex_nonzero_imag', [3+1i 5+2i], ...
      'complex_zero_imag',    complex([3 5], [0 0]), ...
      'char',                 char([3 5]), ...
      'int32',                int32([3 5])};
for k = 1:2:numel(cv)
    nm = cv{k};  fv = cv{k+1};
    fprintf('M1b %-22s isnumeric=%d isreal=%d\n', nm, isnumeric(fv), isreal(fv));
    try
        Sx = invz_spectra_map(ion, 0.31, fv, w, syn());
        fprintf('M1b %-22s RETURNED  Bc_1z=%g status=%s\n', nm, Sx.Bc_1z, Sx.Bc_1z_status);
    catch ME
        fprintf('M1b %-22s THREW %s\n', nm, ME.identifier);
    end
end
% and what the products look like, to explain WHY
fq = complex([3 5], [0 0]);
Bq = fq(:) * [1 0 0];
fprintf('M1b product isreal(fields(:)*fdir) = %d\n', isreal(Bq));

%% ---- I2: hunt for unstable_endpoint / solver_failed / bare_not_ordered -------------------
fprintf('\n########## I2 scan: reasons over a fine strict sweep ##########\n');
Ts = [0.31 0.10 1.00];
Fs = [0.10 0.50 1.00 1.50 2.00 2.50 2.85 3.00 3.10 3.30 4.00 5.50 8.00];
for t = Ts
    o = syn();  o.static_medium = 'strict_1z_dyson_ref';
    ws = warning('off', 'invz:spectraMapMasked');
    S = invz_spectra_map(ion, t, Fs, w, o);
    warning(ws);
    fprintf('T=%.2f  Bc_1z=%g (%s)\n', t, S.Bc_1z, S.Bc_1z_status);
    for k = 1:numel(Fs)
        fprintf('   B=%5.2f  ph1z=%d  reason=%-28s diag=%-22s pmstat=%-18s crit=%.4g\n', ...
            Fs(k), S.phase_1z(k), S.phase_1z_reason{k}, S.ordered_diag_reason{k}, ...
            S.pm_probe_status{k}, S.crit_pm(k));
    end
end
fprintf('DONE_SCAN\n');
