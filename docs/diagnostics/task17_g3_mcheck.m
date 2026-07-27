% Task 17 characterisation item (1): controller-notes.md §7(1) -- G3 pins G0el0=0, not m=0.
% Rerun the G3 loop with tl.m = 0.3 (everything else unchanged) and report whether out.r is
% still 1.25. Read-only measurement; does not modify the committed test.
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(fullfile(REPO, 'invz_projected'));
addpath(REPO);
addpath(fullfile(REPO, 'invz_common'));

beta = 1/(0.0862*0.31);
fprintf('--- baseline (tl.m = 0, as committed) ---\n');
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0 = 2*tl.n01/tl.Delta;
for K0 = [0, 1e-3, 5e-3, 0.05]
    [~, out] = invz_gstat_ordered(tl, [0.01; 0.02], K0, 0.25, beta, -300, 0);
    fprintf('  K0=%-8g  out.r=%.15g  (== 1.25 ? %d)\n', K0, out.r, out.r == 1.25);
end

fprintf('\n--- tl.m = 0.3, everything else unchanged (G0el0 argument still 0) ---\n');
tl2 = tl;  tl2.m = 0.3;
for K0 = [0, 1e-3, 5e-3, 0.05]
    [~, out2] = invz_gstat_ordered(tl2, [0.01; 0.02], K0, 0.25, beta, -300, 0);
    fprintf('  K0=%-8g  out.r=%.15g  (== 1.25 ? %d)  (RelTol 1e-12 match: %d)\n', ...
        K0, out2.r, out2.r == 1.25, abs(out2.r - 1.25) <= 1e-12*max(1,abs(1.25)));
end

fprintf('\n--- for contrast: tl.m = 0.3 WITH a nonzero G0el0 argument (0.4 instead of 0), to show m DOES matter once G0el0 != 0 ---\n');
for K0 = [0, 1e-3, 5e-3, 0.05]
    [~, out3] = invz_gstat_ordered(tl2, [0.01; 0.02], K0, 0.25, beta, -300, 0.4);
    fprintf('  K0=%-8g  out.r=%.15g  (== 1.25 ? %d)\n', K0, out3.r, out3.r == 1.25);
end
fprintf('\n### done ###\n');
