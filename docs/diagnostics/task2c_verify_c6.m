% task2c_verify_c6.m -- read-only verification of C6: is the worst K discrepancy (9.405e-10 at
% F3754215 node 23, cold-vs-seed2) inside or outside the K absolute floor (1e-8*|J0eff|)?
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

kcold = find(strcmp(cfg_ids, 'g1_isolated_cold_F3754215'), 1);
kseed2 = find(strcmp(cfg_ids, 'g1_isolated_seed2_F3754215'), 1);
fprintf('cold cell=%d, seed2 cell=%d\n', kcold, kseed2);

ocold = results(kcold).out;
oseed2 = results(kseed2).out;

% try to locate J0eff wherever it lives
fprintf('--- fieldnames(results) ---\n'); disp(fieldnames(results));
fprintf('--- fieldnames(out) cold ---\n'); disp(fieldnames(ocold));
if isfield(ocold, 'extra'), fprintf('--- fieldnames(out.extra) cold ---\n'); disp(fieldnames(ocold.extra)); end
fprintf('--- fieldnames(nodes) cold ---\n'); disp(fieldnames(ocold.nodes(1)));

ncold = ocold.nodes(23);
nseed2 = oseed2.nodes(23);
fprintf('node 23: cold.accepted=%d seed2.accepted=%d cold.h=%.8g seed2.h=%.8g\n', ...
    ncold.accepted, nseed2.accepted, ncold.h, nseed2.h);

if isfield(ncold, 'K')
    Kc = ncold.K; Ks = nseed2.K;
    dK = abs(Kc(:) - Ks(:));
    [mx, ix] = max(dK);
    fprintf('max|K diff| over K vector = %.6e at index %d (Kc=%.10g, Ks=%.10g)\n', mx, ix, Kc(ix), Ks(ix));
else
    fprintf('node has no .K field directly -- checking node fieldnames again:\n');
    disp(fieldnames(ncold));
end

% try to find J0eff somewhere
for src = {'results(kcold)', 'results(kcold).out', 'results(kcold).out.summary'}
end
if isfield(results(kcold), 'J0eff'), fprintf('results.J0eff = %.10g\n', results(kcold).J0eff); end
if isfield(ocold, 'J0eff'), fprintf('out.J0eff = %.10g\n', ocold.J0eff); end
if isfield(ocold.summary, 'J0eff'), fprintf('out.summary.J0eff = %.10g\n', ocold.summary.J0eff); end
if isfield(ncold, 'J0eff'), fprintf('node.J0eff = %.10g\n', ncold.J0eff); end

fprintf('\nVERIFY_C6_DONE\n');
