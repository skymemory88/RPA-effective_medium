% Pre-fix measurement for the Task-15 review findings F1 (empty sweep) and F2 (counters).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
info = struct('Jcc0', 6.4e-3);
Jnu  = linspace(-2e-3, 6.0e-3, 24).';
syn  = struct('Jnu', Jnu, 'info', info, 'verbose', false);

fprintf('########## F1: empty field sweep, DEFAULT resummed ##########\n');
try
    S = invz_spectra_map(ion, 0.31, [], w, syn);
    fprintf('F1 RETURNED  size(chiz)=[%d %d]  Bc_1z=%g  status=%s\n', ...
            size(S.chiz,1), size(S.chiz,2), S.Bc_1z, S.Bc_1z_status);
catch e
    fprintf('F1 THREW %s : %s\n', e.identifier, e.message);
end

fprintf('########## F1b: empty field sweep, STRICT ##########\n');
o = syn;  o.static_medium = 'strict_1z_dyson_ref';
try
    S = invz_spectra_map(ion, 0.31, [], w, o);
    fprintf('F1b RETURNED status=%s\n', S.Bc_1z_status);
catch e
    fprintf('F1b THREW %s : %s\n', e.identifier, e.message);
end

fprintf('########## F2: B = 0 column, real 4^3 couplings ##########\n');
ropts = struct('grid',[4 4 4],'dpRng',12,'cache',false,'verbose',false,'ordered_1z','jensen');
Sr = invz_spectra_map(ion, 0.31, [0 8], w, ropts);
fprintf('resummed:  reason = %s\n', strjoin(Sr.phase_1z_reason, ' | '));
os = ropts;  os.static_medium = 'strict_1z_dyson_ref';
lastwarn('','invz:none');
Ss = invz_spectra_map(ion, 0.31, [0 8], w, os);
[wmsg, wid] = lastwarn();
fprintf('strict:    reason = %s\n', strjoin(Ss.phase_1z_reason, ' | '));
fprintf('strict:    pm_eid = %s | %s\n', Ss.pm_probe_error_id{1}, Ss.pm_probe_error_id{2});
fprintf('strict:    phase_1z = [%s]  Sigma0 = [%s]\n', num2str(Ss.phase_1z), num2str(Ss.Sigma0));
fprintf('strict:    WARNID %s\n', wid);
fprintf('strict:    WARNMSG %s\n', wmsg);
fprintf('resummed:  Sigma0 = [%s]  Bc_1z = %g\n', num2str(Sr.Sigma0), Sr.Bc_1z);

fprintf('########## F7: strict + ordered_1z = bare ##########\n');
ob = ropts;  ob.static_medium = 'strict_1z_dyson_ref';  ob.ordered_1z = 'bare';
lastwarn('','invz:none');
Sb = invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, ob);
[wmsg2, wid2] = lastwarn();
fprintf('bare strict: reason = %s\n', strjoin(Sb.phase_1z_reason, ' | '));
fprintf('bare strict: Bc_1z = %g  status = %s  WARNID = %s\n', Sb.Bc_1z, Sb.Bc_1z_status, wid2);
obr = ob;  obr.static_medium = 'resummed';
Sbr = invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, obr);
fprintf('bare resummed: Bc_1z = %g\n', Sbr.Bc_1z);
fprintf('PREDONE\n');
