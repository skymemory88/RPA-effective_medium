% task2c_verify_c5c.m -- read-only verification of C5(c): at 2.581023 T, does G3 stride-8
% have identical node-by-node classification vs G1 baseline, or do labels swap at nodes 26/27?
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};
targets = {'g1_swept_F2581023','g3_ds2_swept_F2581023','g3_ds4_swept_F2581023','g3_ds8_swept_F2581023'};
for t = 1:numel(targets)
    k = find(strcmp(cfg_ids, targets{t}), 1);
    if isempty(k), fprintf('NOT FOUND: %s\n', targets{t}); continue; end
    o = results(k).out;
    fprintf('--- %s (cell %d) ---\n', targets{t}, k);
    for j = 24:29
        n = o.nodes(j);
        fprintf('  node id=%d h=%.8g accepted=%d class=%-12s s=%.6g D_uni=%.6g Dq_min=%.6g term=%s\n', ...
            n.id, n.h, n.accepted, n.class, n.s, n.D_uni, n.Dq_min, n.term_reason);
    end
end
fprintf('\nVERIFY_DONE\n');
