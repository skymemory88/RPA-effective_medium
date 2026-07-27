REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};
targets = {'g3_ds2_swept_F1173192','g3_ds4_swept_F1173192','g3_ds8_swept_F1173192'};
for t = 1:numel(targets)
    k = find(strcmp(cfg_ids, targets{t}), 1);
    lp = results(k).lattice_provenance;
    fprintf('--- %s ---  nq=%d  numel(keep_idx)=%d\n', targets{t}, lp.nq, numel(lp.keep_idx));
    fprintf('  construction: %s\n', lp.construction);
end
fprintf('\nVERIFY_C2C_DONE\n');
