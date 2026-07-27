function diag6()
% User-facing symptom: what invz_spectra_map reports per field, resummed vs strict.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag6.log'); if exist(LOG,'file'), delete(LOG); end
ion = invz_ion();  T = 0.31;
fields = [1 2 3 3.6 3.8 4.0 4.2 4.6];
w = (0.005:0.02:0.7).';
for scheme = {'resummed','strict_1z_dyson_ref'}
    o = struct('parallel', false, 'verbose', false, 'static_medium', scheme{1});
    try
        S = invz_spectra_map(ion, T, fields, w, o);
        say(LOG, sprintf('===== invz_spectra_map, static_medium = %s =====', scheme{1}));
        say(LOG, '   B   phase phase_1z  reason                 ordered_diag        chiz_finite  m_1z    D_ord');
        for k = 1:numel(fields)
            say(LOG, sprintf('%5.2f %5d %7d  %-22s %-20s %7d %8.4g %9.4g', fields(k), S.phase(k), S.phase_1z(k), ...
                S.phase_1z_reason{k}, S.ordered_diag_reason{k}, nnz(isfinite(S.chiz(:,k))), S.m_1z(k), S.D_ord(k)));
        end
        say(LOG, sprintf('  Bc_1z = %s ; Bc_auto = %s', num2str(S.Bc_1z), num2str(S.Bc_auto)));
    catch ME
        say(LOG, sprintf('spectra_map(%s) THREW %s : %s', scheme{1}, ME.identifier, ME.message));
    end
end
say(LOG,'done');
end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
