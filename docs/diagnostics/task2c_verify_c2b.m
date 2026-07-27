REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

k1 = find(strcmp(cfg_ids, 'g1_swept_F1173192'), 1);
fprintf('g1 cfg.couplings.grid = '); disp(results(k1).cfg.couplings.grid);
fprintf('g1 cfg.couplings.dpRng = %g\n', results(k1).cfg.couplings.dpRng);

targets = {'g6_offunshifted_dp40_swept_F1173192','g6_offhalfstep_dp30_swept_F1173192', ...
    'g6_offhalfstep_dp40_swept_F1173192','g3_ds8_swept_F1173192','g3_ds2_swept_F1173192'};
for t = 1:numel(targets)
    k = find(strcmp(cfg_ids, targets{t}), 1);
    lp = results(k).lattice_provenance;
    fprintf('--- %s ---\n', targets{t});
    fprintf('  nq=%d  n_gamma_dropped=%s\n', lp.nq, mat2str(lp.n_gamma_dropped));
    if isfield(lp, 'construction')
        c = lp.construction;
        if ischar(c) || isstring(c)
            fprintf('  construction: %s\n', c);
        elseif isstruct(c)
            fprintf('  construction fields: '); disp(fieldnames(c)');
            disp(c);
        end
    end
    if isfield(lp, 'keep_idx'), fprintf('  numel(keep_idx)=%d (stride subset size)\n', numel(lp.keep_idx)); end
end
fprintf('\nVERIFY_C2B_DONE\n');
