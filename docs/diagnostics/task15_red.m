ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));

r = runtests('invz_projected/tests/test_invz_spectra_map_phase_reasons.m');
fprintf('RED P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), ...
        nnz([r.Incomplete]), numel(r));
for k = 1:numel(r)
    st = 'PASS'; if r(k).Failed, st = 'FAIL'; elseif r(k).Incomplete, st = 'INC'; end
    fprintf('RESULT %-6s %s\n', st, r(k).Name);
end
