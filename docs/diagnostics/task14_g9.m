% Task 14 G9 probe: run the SAME fixtures against the HEAD (e58ec17) point solvers and against
% the task-14 ones, and compare every shared member BITWISE (isequaln, no tolerance).
% MODE = 'ref' shadows invz_projected with the extracted e58ec17 files; MODE = 'new' does not.
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
SCR  = ['/private/tmp/claude-503/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-' ...
        'Programming-scripts-Matlab-Simulation-invZ-expansion/' ...
        'e72c63b7-dc11-4663-9d4b-867feace2d26/scratchpad'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT, 'invz_common'));
addpath(genpath(fullfile(ROOT, 'invz_projected')));
if strcmp(MODE, 'ref')
    addpath(fullfile(SCR, 'task14_legacy'));      % prepended => shadows invz_projected
end
fprintf('WHICH invz_solve_point -> %s\n', which('invz_solve_point'));
fprintf('WHICH invz_solve_point_ordered -> %s\n', which('invz_solve_point_ordered'));

ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o   = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);

out = struct();
out.F1_pm_std      = invz_solve_point(ion, 0.31, [2.85 0 0], Jnu, o);
out.F2_pm_bitfix   = invz_solve_point(ion, 1.6, 0.5, Jnu, struct('J0eff', 6.4e-3));
out.F3_pm_ordside  = invz_solve_point(ion, 0.31, [0.5 0 0], Jnu, o);
out.F4_ord_bare    = invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
oj = o;  oj.ordered_mode = 'jensen';
out.F5_ord_jensen  = invz_solve_point_ordered(ion, 0.31, [2.85 0 0], Jnu, oj);
out.F6_ord_early   = invz_solve_point_ordered(ion, 2.0, [2.85 0 0], Jnu, o);

sf = fullfile(SCR, ['task14_g9_' MODE '.mat']);
save(sf, 'out');
fprintf('G9 %s saved -> %s\n', MODE, sf);
