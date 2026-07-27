ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();  T = 0.31;  w = (0.02:0.02:0.6).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
base = struct('Jnu', Jnu, 'info', info, 'verbose', false);

fprintf('=== resummed (default), fields [2.85 3.30 5.50] ===\n');
S = invz_spectra_map(ion, T, [2.85 3.30 5.50], w, base);
fprintf('static_medium   = %s\n', S.static_medium);
fprintf('phase           = %s\n', mat2str(S.phase));
fprintf('phase_1z        = %s\n', mat2str(S.phase_1z));
fprintf('reason          = %s\n', strjoin(S.phase_1z_reason, ' | '));
fprintf('sm_used         = %s\n', strjoin(S.static_medium_used, ' | '));
fprintf('pm_status       = %s\n', strjoin(S.pm_probe_status, ' | '));
fprintf('pm_err          = [%s]\n', strjoin(S.pm_probe_error_id, ' | '));
fprintf('stability_1z    = %s\n', mat2str(S.stability_1z));
fprintf('Bc_1z           = %g   Bc_auto = %g\n', S.Bc_1z, S.Bc_auto);
fprintf('Bc_1z_interval  = %s   status = %s\n', mat2str(S.Bc_1z_interval), S.Bc_1z_status);

for scheme = {'strict_1z_dyson_ref', 'strict_1z_bare_ref'}
    fprintf('\n=== %s ===\n', scheme{1});
    o = base;  o.static_medium = scheme{1};
    try
        S2 = invz_spectra_map(ion, T, [2.85 3.30 5.50], w, o);
        fprintf('static_medium   = %s\n', S2.static_medium);
        fprintf('phase_1z        = %s\n', mat2str(S2.phase_1z));
        fprintf('reason          = %s\n', strjoin(S2.phase_1z_reason, ' | '));
        fprintf('sm_used         = %s\n', strjoin(S2.static_medium_used, ' | '));
        fprintf('pm_status       = %s\n', strjoin(S2.pm_probe_status, ' | '));
        fprintf('crit_pm         = %s\n', mat2str(S2.crit_pm, 6));
        fprintf('stability_1z    = %s\n', mat2str(S2.stability_1z));
        fprintf('Bc_1z           = %g  interval = %s  status = %s\n', ...
                S2.Bc_1z, mat2str(S2.Bc_1z_interval), S2.Bc_1z_status);
        fprintf('any finite chiz = %d ; any finite chirpa = %d\n', ...
                any(isfinite(S2.chiz(:))), any(isfinite(S2.chirpa(:))));
    catch err
        fprintf('STRICT MAP ERROR %s : %s\n', err.identifier, err.message);
        for s = 1:min(4, numel(err.stack))
            fprintf('   at %s:%d\n', err.stack(s).name, err.stack(s).line);
        end
    end
end

fprintf('\n=== bare escape hatch under strict ===\n');
o = base;  o.static_medium = 'strict_1z_dyson_ref';  o.ordered_1z = 'bare';
S3 = invz_spectra_map(ion, T, [2.85 3.30 5.50], w, o);
fprintf('phase_1z = %s\nreason   = %s\nsm_used  = %s\n', mat2str(S3.phase_1z), ...
        strjoin(S3.phase_1z_reason, ' | '), strjoin(S3.static_medium_used, ' | '));

fprintf('\n=== solve_opts reservation ===\n');
for f = {'static_medium', 'ref_margin'}
    try
        o = base;  o.solve_opts = struct(f{1}, 'resummed');
        invz_spectra_map(ion, T, 5.5, w, o);
        fprintf('%-14s NO ERROR (BAD)\n', f{1});
    catch err
        fprintf('%-14s -> %s\n', f{1}, err.identifier);
    end
end
