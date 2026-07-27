ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

r = runtests('invz_projected/tests');
fprintf('SUITE P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), ...
        nnz([r.Incomplete]), numel(r));
f = find([r.Failed]);
if ~isempty(f)
    for k = 1:numel(f), fprintf('FAILEDNAME %s\n', r(f(k)).Name); end
else
    fprintf('FAILEDNONE\n');
end
% record the incomplete SET by name, so "24" can be compared as identity, not just count
inc = find([r.Incomplete]);
fprintf('INCOMPLETECOUNT %d\n', numel(inc));
for k = 1:numel(inc), fprintf('INC %s\n', r(inc(k)).Name); end
