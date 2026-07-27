REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

targets = {'g1_swept_F1173192','g6_offunshifted_dp40_swept_F1173192','g6_offhalfstep_dp30_swept_F1173192', ...
    'g6_offhalfstep_dp40_swept_F1173192','g3_ds8_swept_F1173192'};
for t = 1:numel(targets)
    k = find(strcmp(cfg_ids, targets{t}), 1);
    if isempty(k), fprintf('NOT FOUND %s\n', targets{t}); continue; end
    cf = results(k).cfg;
    lp = results(k).lattice_provenance;
    fprintf('--- %s ---\n', targets{t});
    fprintf('  fieldnames(cfg): '); disp(fieldnames(cf)');
    if isstruct(cf.couplings)
        fprintf('  fieldnames(cfg.couplings): '); disp(fieldnames(cf.couplings)');
    end
    if isstruct(lp) && ~isempty(fieldnames(lp))
        fprintf('  fieldnames(lattice_provenance): '); disp(fieldnames(lp)');
        if isfield(lp,'nq'), fprintf('  lp.nq = %d\n', lp.nq); end
        if isfield(lp,'grid'), disp(lp.grid); end
    else
        fprintf('  lattice_provenance: empty/none\n');
    end
end
fprintf('\nVERIFY_C2_DONE\n');
