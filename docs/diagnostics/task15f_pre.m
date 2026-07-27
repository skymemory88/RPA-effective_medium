% Task-15 final fix pass -- PRE-fix measurements (I1, I2, M1, M5).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

ion = invz_ion();
w   = (0.02:0.04:0.42).';
syn = @() struct('Jnu', linspace(-2e-3, 6.0e-3, 24).', 'info', struct('Jcc0', 6.4e-3), ...
                 'verbose', false);
T   = 0.31;

j = @(c) strjoin(c, ' | ');
b = @(v) strjoin(arrayfun(@(x) sprintf('%d', x), v, 'uni', 0), ' ');

%% ---------------- I1: stability_1z on 'pm' columns ----------------------------------------
fprintf('########## I1: stability_1z per column, all schemes x ordered_1z ##########\n');
schemes = {'resummed', 'strict_1z_dyson_ref', 'strict_1z_bare_ref'};
o1zs    = {'jensen', 'bare'};
F       = [2.85 3.30 5.50];
for a = 1:numel(schemes)
    for c = 1:numel(o1zs)
        o = syn();  o.static_medium = schemes{a};  o.ordered_1z = o1zs{c};
        ws = warning('off', 'invz:spectraMapMasked');
        S = invz_spectra_map(ion, T, F, w, o);
        warning(ws);
        fprintf('%-20s o1z=%-7s reason = %-45s stab = [%s]  crit_pm = [%s]\n', ...
            schemes{a}, o1zs{c}, j(S.phase_1z_reason), b(S.stability_1z), ...
            strjoin(arrayfun(@(x) sprintf('%.4g', x), S.crit_pm, 'uni', 0), ' '));
        ispm = strcmp(S.phase_1z_reason, 'pm');
        fprintf('%-20s o1z=%-7s   PM columns: n=%d  all stab==1 ? %d\n', ...
            schemes{a}, o1zs{c}, nnz(ispm), all(S.stability_1z(ispm)));
    end
end

%% ---------------- I1b: does the ordered leg ever run on a 'pm' column? ---------------------
fprintf('\n########## I1b: ordered_diag_reason on the same columns ##########\n');
for a = 1:numel(schemes)
    o = syn();  o.static_medium = schemes{a};
    ws = warning('off', 'invz:spectraMapMasked');
    S = invz_spectra_map(ion, T, F, w, o);
    warning(ws);
    fprintf('%-20s reason=%-32s diag=%-40s stab=[%s]\n', schemes{a}, ...
        j(S.phase_1z_reason), j(S.ordered_diag_reason), b(S.stability_1z));
end

%% ---------------- I2: which masked reasons appear, and is the sweep silent? ----------------
fprintf('\n########## I2: masked reasons vs. the warning trigger ##########\n');
cases = { ...
  'tiltB',   struct('field_dir', [cosd(3) 0 sind(3)], 'fields', [2.85 3.30 5.50]), ...
  'wide',    struct('field_dir', [1 0 0],             'fields', [0.5 1.5 2.85 3.30 4.0 5.50 8.0]) };
for k = 1:2:numel(cases)
    nm = cases{k};  cf = cases{k+1};
    o = syn();  o.static_medium = 'strict_1z_dyson_ref';  o.field_dir = cf.field_dir;
    lastwarn('', 'invz:none');
    ws = warning('off', 'invz:spectraMapMasked');
    S = invz_spectra_map(ion, T, cf.fields, w, o);
    warning(ws);
    [wm, wi] = lastwarn();
    fprintf('%-8s reason = %s\n', nm, j(S.phase_1z_reason));
    fprintf('%-8s Bc_1z = %g (%s)  WARNID %s\n', nm, S.Bc_1z, S.Bc_1z_status, wi);
    if ~isempty(wm), fprintf('%-8s WARNMSG %s\n', nm, wm); end
end

%% ---------------- M1: char / complex fields on the DEFAULT path ---------------------------
fprintf('\n########## M1: char and complex `fields` on the default resummed path ##########\n');
probes = {'char', char([3 5]), 'complex', [3+1i 5+2i]};
for k = 1:2:numel(probes)
    nm = probes{k};  fv = probes{k+1};
    try
        Sx = invz_spectra_map(ion, T, fv, w, syn());
        fprintf('M1 %-8s RETURNED  Bc_1z=%g status=%s nB=%d\n', nm, Sx.Bc_1z, ...
                Sx.Bc_1z_status, numel(Sx.phase_1z));
    catch ME
        fprintf('M1 %-8s THREW %s : %s\n', nm, ME.identifier, ME.message);
    end
end

%% ---------------- M5: the synthetic fixture the gate should use --------------------------
fprintf('\n########## M5: synthetic fixture [0 5.50] under strict ##########\n');
o = syn();  o.ordered_1z = 'jensen';  o.static_medium = 'strict_1z_dyson_ref';
lastwarn('', 'invz:none');
tic;  S = invz_spectra_map(ion, T, [0 5.50], w, o);  el = toc;
[wm, wi] = lastwarn();
fprintf('M5 elapsed  = %.2f s\n', el);
fprintf('M5 phase_1z = [%s]\n', b(S.phase_1z));
fprintf('M5 reason   = %s\n', j(S.phase_1z_reason));
fprintf('M5 diag     = %s\n', j(S.ordered_diag_reason));
fprintf('M5 pm_eid   = <%s>\n', j(S.pm_probe_error_id));
fprintf('M5 WARNID   = %s\n', wi);
fprintf('M5 WARNMSG  = %s\n', wm);
fprintf('M5 stab     = [%s]\n', b(S.stability_1z));
fprintf('DONE_PRE\n');
