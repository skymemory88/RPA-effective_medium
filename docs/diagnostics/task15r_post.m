% Post-fix measurement for the Task-15 review findings F1, F2, F3, F7.
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

fprintf('########## F1: empty field sweep ##########\n');
S = invz_spectra_map(ion, 0.31, [], w, syn);
fprintf('F1 resummed RETURNED size(chiz)=[%d %d] Bc_1z=%g status=%s nreason=%d\n', ...
        size(S.chiz,1), size(S.chiz,2), S.Bc_1z, S.Bc_1z_status, numel(S.phase_1z_reason));
o = syn;  o.static_medium = 'strict_1z_dyson_ref';
Ss = invz_spectra_map(ion, 0.31, [], w, o);
fprintf('F1 strict   RETURNED size(chiz)=[%d %d] Bc_1z=%g status=%s\n', ...
        size(Ss.chiz,1), size(Ss.chiz,2), Ss.Bc_1z, Ss.Bc_1z_status);

fprintf('########## F2: B = 0 column, real 4^3 couplings ##########\n');
ropts = struct('grid',[4 4 4],'dpRng',12,'cache',false,'verbose',false,'ordered_1z','jensen');
Sr = invz_spectra_map(ion, 0.31, [0 8], w, ropts);
fprintf('resummed: reason      = %s\n', strjoin(Sr.phase_1z_reason, ' | '));
fprintf('resummed: diag_reason = %s\n', strjoin(Sr.ordered_diag_reason, ' | '));
os = ropts;  os.static_medium = 'strict_1z_dyson_ref';
lastwarn('','invz:none');
Sx = invz_spectra_map(ion, 0.31, [0 8], w, os);
[wmsg, wid] = lastwarn();
fprintf('strict:   reason      = %s\n', strjoin(Sx.phase_1z_reason, ' | '));
fprintf('strict:   diag_reason = %s\n', strjoin(Sx.ordered_diag_reason, ' | '));
fprintf('strict:   pm_eid      = <%s> | <%s>\n', Sx.pm_probe_error_id{1}, Sx.pm_probe_error_id{2});
fprintf('strict:   resp_eid    = <%s> | <%s>\n', Sx.response_error_id{1}, Sx.response_error_id{2});
fprintf('strict:   phase_1z    = [%s]\n', num2str(Sx.phase_1z));
fprintf('strict:   WARNID  %s\n', wid);
fprintf('strict:   WARNMSG %s\n', wmsg);

fprintf('########## F7: strict + ordered_1z = bare, synthetic ##########\n');
ob = syn;  ob.ordered_1z = 'bare';  ob.static_medium = 'strict_1z_dyson_ref';
lastwarn('','invz:none');
Sb = invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, ob);
[wmsg2, wid2] = lastwarn();
fprintf('bare strict:   reason = %s\n', strjoin(Sb.phase_1z_reason, ' | '));
fprintf('bare strict:   Bc_1z = %g  status = %s\n', Sb.Bc_1z, Sb.Bc_1z_status);
fprintf('bare strict:   WARNID  %s\n', wid2);
fprintf('bare strict:   WARNMSG %s\n', wmsg2);
obr = ob;  obr.static_medium = 'resummed';
Sbr = invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, obr);
fprintf('bare resummed: Bc_1z = %g (finite -> the asymmetry the warning must surface)\n', Sbr.Bc_1z);

fprintf('########## F3: response boundaries absorb nothing ##########\n');
for sc = {'resummed','strict_1z_dyson_ref'}
    oc = syn;  oc.static_medium = sc{1};
    Sc = invz_spectra_map(ion, 0.31, [2.85 5.50], w, oc);
    fprintf('%-20s resp_eid all empty = %d  any response_failed = %d  reason = %s\n', ...
        sc{1}, all(cellfun(@isempty, Sc.response_error_id)), ...
        any(strcmp(Sc.phase_1z_reason,'response_failed')), strjoin(Sc.phase_1z_reason,' | '));
end
fprintf('POSTDONE\n');
