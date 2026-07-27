ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(genpath(fullfile(ROOT, 'invz_projected')));
addpath(fullfile(ROOT, 'invz_common'));

r = runtests('invz_projected/tests/test_invz_twolevel_ordered_domain.m');
disp(table(r));
fprintf('GREEN P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), ...
        nnz([r.Incomplete]), numel(r));
for k = 1:numel(r)
    fprintf('NAME %s PASSED=%d\n', r(k).Name, r(k).Passed);
end
