SCR = ['/private/tmp/claude-503/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-' ...
       'Programming-scripts-Matlab-Simulation-invZ-expansion/' ...
       'e72c63b7-dc11-4663-9d4b-867feace2d26/scratchpad'];
A = load(fullfile(SCR, 'task14_g9_ref.mat'));  A = A.out;
B = load(fullfile(SCR, 'task14_g9_new.mat'));  B = B.out;
fx = fieldnames(A);
nbad = 0;
for i = 1:numel(fx)
    a = A.(fx{i});  b = B.(fx{i});
    fa = fieldnames(a);  fb = fieldnames(b);
    miss = setdiff(fa, fb);
    added = setdiff(fb, fa);
    for k = 1:numel(fa)
        if ~isequaln(a.(fa{k}), b.(fa{k}))
            nbad = nbad + 1;
            fprintf('MISMATCH %s.%s\n', fx{i}, fa{k});
        end
    end
    fprintf('%-16s shared=%d bitwise_equal=%s  dropped={%s}  added={%s}\n', fx{i}, ...
        numel(fa), string(isempty(setdiff(fa,fb)) && all(cellfun(@(f) isequaln(a.(f), b.(f)), fa))), ...
        strjoin(miss', ','), strjoin(added', ','));
end
fprintf('G9COMPARE mismatches=%d\n', nbad);
% headline numbers, printed at full precision so the report can quote them
fprintf('F1 crit  ref=%.17g new=%.17g\n', A.F1_pm_std.crit, B.F1_pm_std.crit);
fprintf('F1 Sigma0 ref=%.17g new=%.17g\n', A.F1_pm_std.Sigma0, B.F1_pm_std.Sigma0);
fprintf('F5 hmf   ref=%.17g new=%.17g\n', A.F5_ord_jensen.hmf, B.F5_ord_jensen.hmf);
fprintf('F5 D_uni ref=%.17g new=%.17g\n', A.F5_ord_jensen.D_uni, B.F5_ord_jensen.D_uni);
fprintf('F5 conv  ref=%d new=%d\n', A.F5_ord_jensen.converged, B.F5_ord_jensen.converged);
fprintf('NEWFIELDS F5: %s\n', strjoin(setdiff(fieldnames(B.F5_ord_jensen), ...
                                              fieldnames(A.F5_ord_jensen))', ', '));
fprintf('NEWFIELDS F6(early): %s\n', strjoin(setdiff(fieldnames(B.F6_ord_early), ...
                                              fieldnames(A.F6_ord_early))', ', '));
fprintf('F5 stable_1z=%d crit_1z=%.17g Dq_min=%.17g path_omit_max=%.17g omit_max=%.17g\n', ...
    B.F5_ord_jensen.stable_1z, B.F5_ord_jensen.crit_1z, B.F5_ord_jensen.Dq_min, ...
    B.F5_ord_jensen.path_omit_max, B.F5_ord_jensen.omit_max);
fprintf('F5 static_medium=%s medium_status=%s denom=%.17g margin=%.17g\n', ...
    B.F5_ord_jensen.static_medium, B.F5_ord_jensen.medium_status, ...
    B.F5_ord_jensen.medium_denom, B.F5_ord_jensen.medium_margin);
