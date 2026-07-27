REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

kswept = find(strcmp(cfg_ids, 'g1_swept_F2850000'), 1);
kiso   = find(strcmp(cfg_ids, 'g1_isolated_cold_F2850000'), 1);
fprintf('swept cell=%d, isolated_cold cell=%d\n', kswept, kiso);

oswept = results(kswept).out;
oiso   = results(kiso).out;
nswept = oswept.nodes(27);
niso   = oiso.nodes(27);
fprintf('node 27: swept.accepted=%d iso.accepted=%d swept.h=%.8g iso.h=%.8g\n', ...
    nswept.accepted, niso.accepted, nswept.h, niso.h);

Kw = nswept.state.K; Ki = niso.state.K;
dK = abs(Kw(:) - Ki(:));
[mx, ix] = max(dK);
fprintf('max|K diff| = %.6e at index %d (Kswept=%.10g, Kiso=%.10g)\n', mx, ix, Kw(ix), Ki(ix));

J0eff = oswept.meta.J0eff;
floorAbs = 1e-8 * abs(J0eff);
relTerm = 1e-6 * max(abs(Kw(ix)), abs(Ki(ix)));
fprintf('J0eff=%.10g  AbsTol_q floor=%.6e  ratio(maxdK/floor)=%.4f\n', J0eff, floorAbs, mx/floorAbs);
fprintf('RelTol*max(|a|,|b|)=%.6e  combined bound=%.6e  passes=%d\n', relTerm, floorAbs+relTerm, mx <= floorAbs+relTerm);

fprintf('\nVERIFY_C6C_DONE\n');
