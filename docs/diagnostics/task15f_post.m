% Task-15 final fix pass -- POST-fix measurements (I1, I2, M1, M5).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @(s) struct('Jnu', s*linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', s*6.4e-3), ...
                  'verbose', false);
j = @(c) strjoin(c, ' | ');
b = @(v) strjoin(arrayfun(@(x) sprintf('%d', x), v, 'uni', 0), ' ');

%% ---------------- I1: the corrected sentence, measured ------------------------------------
fprintf('########## I1: stability_1z on ''pm'' columns ##########\n');
schemes = {'resummed', 'strict_1z_dyson_ref', 'strict_1z_bare_ref'};
allpm = true;  anyOrderedRanOnPm = false;
for a = 1:numel(schemes)
    for o1z = {'jensen', 'bare'}
        o = syn(1);  o.static_medium = schemes{a};  o.ordered_1z = o1z{1};
        ws = warning('off', 'invz:spectraMapMasked');
        S = invz_spectra_map(ion, 0.31, [2.85 3.30 5.50], w, o);
        warning(ws);
        ispm = strcmp(S.phase_1z_reason, 'pm');
        allpm = allpm && all(S.stability_1z(ispm));
        anyOrderedRanOnPm = anyOrderedRanOnPm || ...
            any(~strcmp(S.ordered_diag_reason(ispm), 'not_attempted'));
        fprintf('%-20s o1z=%-7s reason=%-38s stab=[%s] diag=%s\n', schemes{a}, o1z{1}, ...
                j(S.phase_1z_reason), b(S.stability_1z), j(S.ordered_diag_reason));
    end
end
fprintf('I1 VERDICT  every ''pm'' column stability_1z==true : %d\n', allpm);
fprintf('I1 VERDICT  ordered leg ever ran on a ''pm'' column : %d  (0 => not pt.stable_1z)\n', ...
        anyOrderedRanOnPm);

%% ---------------- I2: the previously-undercounted masked column ---------------------------
fprintf('\n########## I2: masked ''solver_failed'' column, weak synthetic couplings ##########\n');
o = syn(0.3);  o.static_medium = 'strict_1z_dyson_ref';
F = [0 2.85 3.30 5.50];
lastwarn('', 'invz:none');
S = invz_spectra_map(ion, 0.31, F, w, o);
[wm, wi] = lastwarn();
fprintf('I2 reason  = %s\n', j(S.phase_1z_reason));
fprintf('I2 n_msk   = %d   Bc_1z = %g (%s)\n', nnz(S.phase_1z == 0), S.Bc_1z, S.Bc_1z_status);
fprintf('I2 WARNID  = %s\nI2 WARNMSG = %s\n', wi, wm);

fprintf('\n########## I2b: real 4^3 sweep, headline vs breakdown ##########\n');
o = struct('grid', [4 4 4], 'dpRng', 12, 'cache', false, 'verbose', false, ...
           'static_medium', 'strict_1z_dyson_ref');
lastwarn('', 'invz:none');
S = invz_spectra_map(ion, 0.31, [0 1 3 5], w, o);
[wm, wi] = lastwarn();
fprintf('I2b reason  = %s\n', j(S.phase_1z_reason));
fprintf('I2b diag    = %s\n', j(S.ordered_diag_reason));
fprintf('I2b WARNMSG = %s\n', wm);

%% ---------------- M1: char / complex on the DEFAULT path ----------------------------------
fprintf('\n########## M1: char / complex `fields` ##########\n');
probes = {'char', char([3 5]), 'complex_zero_imag', complex([3 5], [0 0]), ...
          'complex_nonzero_imag', [3+1i 5+2i]};
for k = 1:2:numel(probes)
    nm = probes{k};  fv = probes{k+1};
    for sch = {'resummed', 'strict_1z_dyson_ref'}
        o = syn(1);  o.static_medium = sch{1};
        try
            ws = warning('off', 'invz:spectraMapMasked');
            Sx = invz_spectra_map(ion, 0.31, fv, w, o);
            warning(ws);
            fprintf('M1 %-22s %-20s RETURNED Bc_1z=%-8g status=%s\n', nm, sch{1}, ...
                    Sx.Bc_1z, Sx.Bc_1z_status);
        catch ME
            warning(ws);
            fprintf('M1 %-22s %-20s THREW %s\n', nm, sch{1}, ME.identifier);
        end
    end
end

%% ---------------- M5: the ungated fixture --------------------------------------------------
fprintf('\n########## M5: synthetic fixture [0 5.50] under strict ##########\n');
o = syn(1);  o.ordered_1z = 'jensen';  o.static_medium = 'strict_1z_dyson_ref';
lastwarn('', 'invz:none');
tic;  S = invz_spectra_map(ion, 0.31, [0 5.50], w, o);  el = toc;
[wm, wi] = lastwarn();
fprintf('M5 elapsed %.2f s\n', el);
fprintf('M5 phase_1z=[%s] reason=%s  diag=%s  pm_eid=<%s>\n', b(S.phase_1z), ...
        j(S.phase_1z_reason), j(S.ordered_diag_reason), j(S.pm_probe_error_id));
fprintf('M5 WARNID %s\nM5 WARNMSG %s\n', wi, wm);
fprintf('M5 contains ''1 degenerate-doublet'' = %d ; ''0 unknown-PM-probe'' = %d\n', ...
        contains(wm, '1 degenerate-doublet'), contains(wm, '0 unknown-PM-probe'));
fprintf('DONE_POST\n');
