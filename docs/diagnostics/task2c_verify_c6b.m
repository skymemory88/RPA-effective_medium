REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

kcold = find(strcmp(cfg_ids, 'g1_isolated_cold_F3754215'), 1);
kseed2 = find(strcmp(cfg_ids, 'g1_isolated_seed2_F3754215'), 1);

ocold = results(kcold).out;
oseed2 = results(kseed2).out;
ncold = ocold.nodes(23);
nseed2 = oseed2.nodes(23);

fprintf('--- fieldnames(node.state) ---\n'); disp(fieldnames(ncold.state));
fprintf('--- fieldnames(results(kcold).cfg) ---\n'); disp(fieldnames(results(kcold).cfg));

Kc = ncold.state.K; Ks = nseed2.state.K;
dK = abs(Kc(:) - Ks(:));
[mx, ix] = max(dK);
fprintf('numel(K)=%d\n', numel(Kc));
fprintf('max|K diff| = %.6e at index %d (Kc=%.10g, Ks=%.10g)\n', mx, ix, Kc(ix), Ks(ix));

cf = results(kcold).cfg;
if isfield(cf, 'J0eff'), fprintf('cfg.J0eff = %.10g\n', cf.J0eff); end
% search recursively for J0eff anywhere in cfg/out
function search_j0eff(s, path)
    if isstruct(s)
        fn = fieldnames(s);
        for i = 1:numel(fn)
            f = fn{i};
            if numel(s) > 1
                v = s(1).(f);
            else
                v = s.(f);
            end
            newpath = [path '.' f];
            if strcmpi(f, 'J0eff') && isnumeric(v)
                fprintf('FOUND %s = %.10g\n', newpath, v);
            elseif isstruct(v)
                search_j0eff(v, newpath);
            end
        end
    end
end
search_j0eff(cf, 'cfg');
search_j0eff(ocold.meta, 'out.meta');
search_j0eff(ocold.summary, 'out.summary');
search_j0eff(ncold, 'node');
search_j0eff(ncold.state, 'node.state');

fprintf('\nVERIFY_C6B_DONE\n');
