% Task 17 Step 3 runner: INVZ_SLOW=1 run of the two new strict-medium gate test files.
% Absolute paths throughout because `run()` executes this script with the script's own
% directory context, not the repo root (standing rule, task-17-controller-notes.md §9).
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(genpath(fullfile(REPO, 'invz_projected')));
addpath(fullfile(REPO, 'invz_common'));

testfiles = {fullfile(REPO, 'invz_projected', 'tests', 'test_invz_strict_identities.m'), ...
             fullfile(REPO, 'invz_projected', 'tests', 'test_invz_strict_contracts.m')};

fprintf('INVZ_SLOW env = [%s]\n', getenv('INVZ_SLOW'));
fprintf('=== runtests start: %s ===\n', datestr(now));
t0 = tic;
r = runtests(testfiles);
wall_total = toc(t0);
fprintf('=== runtests end ===\n\n');

disp(table(r));

fprintf('\n--- Per-test details ---\n');
for i = 1:numel(r)
    fprintf('%-95s Passed=%d Failed=%d Incomplete=%d Duration=%.4f s\n', ...
        r(i).Name, r(i).Passed, r(i).Failed, r(i).Incomplete, r(i).Duration);
end

names = {r.Name};
prefixes = cell(size(names));
for i = 1:numel(names)
    parts = strsplit(names{i}, '/');
    prefixes{i} = parts{1};
end
ufiles = unique(prefixes, 'stable');
fprintf('\n--- Per-file summary (grouped by TestResult Name prefix) ---\n');
for k = 1:numel(ufiles)
    mask = strcmp(prefixes, ufiles{k});
    fprintf('%s: n=%d Passed=%d Failed=%d Incomplete=%d SumOfTestDurations=%.4f s\n', ...
        ufiles{k}, sum(mask), sum([r(mask).Passed]), sum([r(mask).Failed]), ...
        sum([r(mask).Incomplete]), sum([r(mask).Duration]));
end

fprintf('\nTOTAL wall-clock for combined runtests() call (tic/toc): %.4f s\n', wall_total);
fprintf('Aggregate: Passed=%d Failed=%d Incomplete=%d Total=%d\n', ...
    sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]), numel(r));
