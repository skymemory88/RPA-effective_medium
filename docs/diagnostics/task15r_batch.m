ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(ROOT); addpath(fullfile(ROOT,'invz_common'));
addpath(genpath(fullfile(ROOT,'invz_projected')));
r = runtests({'invz_projected/tests/test_invz_spectra_map_phase_reasons.m', ...
              'invz_projected/tests/test_invz_spectra_map.m', ...
              'invz_projected/tests/test_invz_ewald_integration.m'});
fprintf('BATCH P=%d F=%d I=%d T=%d\n', nnz([r.Passed]), nnz([r.Failed]), ...
        nnz([r.Incomplete]), numel(r));
for k = 1:numel(r)
    if r(k).Failed,      fprintf('RESULT FAIL %s\n', r(k).Name);
    elseif r(k).Incomplete, fprintf('RESULT INC  %s\n', r(k).Name);
    end
end
fprintf('BATCHDONE\n');
