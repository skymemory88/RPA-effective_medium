ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
r = runtests('invz_projected/tests/test_invz_emt_scalar_strict.m');
disp(table(r));
fprintf('P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r));
