REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
k = find(strcmp({results.cfg_id}, 'g3_ds2_swept_F1173192'), 1);
o = results(k).out; lp = results(k).lattice_provenance;
disp('lp fields:'); disp(fieldnames(lp));
disp('size(lp.Jnu_unflat):'); disp(size(lp.Jnu_unflat));
disp('lp.nq:'); disp(lp.nq);
disp('size(lp.qc):'); disp(size(lp.qc));
n = o.nodes(25);
fprintf('idx_pos_flat=%g idx_neg_flat=%g\n', n.idx_pos_flat, n.idx_neg_flat);
fprintf('qbranch_pos.q_idx=%g\n', n.qbranch_pos.q_idx);
meta = struct('is_synthetic', false, 'Jnu_unflat', lp.Jnu_unflat, 'nq', lp.nq);
try
    [qi,bi,Jq] = invz_ordered_trace_resolve(meta, n.idx_pos_flat);
    fprintf('SUCCESS: qi=%d bi=%d Jq=%g\n', qi, bi, Jq);
catch ME
    fprintf('CAUGHT: %s : %s\n', ME.identifier, ME.message);
end
grp = results(k).group;
fprintf('grp=%s class test: %d\n', grp, (isnan(n.qbranch_pos.q_idx) && strcmp(grp,'G3')) || (isnan(n.qbranch_pos.q_idx) && strcmp(grp,'G6')));
