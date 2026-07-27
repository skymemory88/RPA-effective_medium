ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
files = {'invz_projected/invz_spectra_map.m', ...
         'invz_projected/tests/test_invz_spectra_map_phase_reasons.m'};
for k = 1:numel(files)
    m = checkcode(files{k});
    fprintf('=== WORKTREE %s : %d messages ===\n', files{k}, numel(m));
    for i = 1:numel(m), fprintf('  L %d: %s\n', m(i).line, m(i).message); end
end
d = fullfile(tempdir, 't15head');
if ~exist(d, 'dir'), mkdir(d); end
system(sprintf('git show 73dc55c:invz_projected/invz_spectra_map.m > "%s"', ...
       fullfile(d, 'invz_spectra_map.m')));
system(sprintf('git show 73dc55c:invz_projected/tests/test_invz_spectra_map_phase_reasons.m > "%s"', ...
       fullfile(d, 'test_invz_spectra_map_phase_reasons.m')));
hf = {'invz_spectra_map.m', 'test_invz_spectra_map_phase_reasons.m'};
for k = 1:numel(hf)
    m = checkcode(fullfile(d, hf{k}));
    fprintf('=== 73dc55c %s : %d messages ===\n', hf{k}, numel(m));
    for i = 1:numel(m), fprintf('  L %d: %s\n', m(i).line, m(i).message); end
end
