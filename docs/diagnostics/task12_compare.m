root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
A = load(fullfile(root,'.superpowers','sdd','task12_cap_pre.mat'));   pre  = A.C;
B = load(fullfile(root,'.superpowers','sdd','task12_cap_post.mat'));  post = B.C;
% prof fields that existed BEFORE task 12 must be BITWISE unchanged on every case.
old = {'hgrid','r','h0','m','Sigma0','K0','D_uni','G0bare','Gstat','node_conv','F', ...
       'slope0','Sigma0_pm0','K0_pm0','J0eff','n_extend','hmin_initial','status', ...
       'redensified','m_star','D_uni_star','r_star','Gstat_star'};
cases = {'A','B','C','D'};
nbad = 0;
for ci = 1:numel(cases)
    c = cases{ci};
    if ~isequaln(pre.(c).h, post.(c).h)
        fprintf('CMP %s h  MISMATCH pre=%.17g post=%.17g\n', c, pre.(c).h, post.(c).h); nbad=nbad+1;
    end
    for k = 1:numel(old)
        f = old{k};
        if ~isequaln(pre.(c).prof.(f), post.(c).prof.(f))
            fprintf('CMP %s prof.%s MISMATCH\n', c, f);  nbad = nbad + 1;
        end
    end
    fprintf('CMP %s : h + %d pre-existing prof fields checked\n', c, numel(old));
end
% trace-node invariants (report item 2) + full pre-existing trace-node field parity
tf = {'n_nodes','n_iters','node_h','node_phase','node_term','node_ok','node_id','node_Duni', ...
      'node_resid','node_seedkind','node_seedfrom','node_outer','node_hitmax','node_dSbreak'};
for c = {'C','D'}
    for k = 1:numel(tf)
        if ~isequaln(pre.(c{1}).(tf{k}), post.(c{1}).(tf{k}))
            fprintf('CMP %s trace.%s MISMATCH\n', c{1}, tf{k});  nbad = nbad + 1;
        end
    end
    fprintf('CMP %s trace: numel(nodes) pre=%d post=%d | numel(iters) pre=%d post=%d\n', ...
        c{1}, pre.(c{1}).n_nodes, post.(c{1}).n_nodes, pre.(c{1}).n_iters, post.(c{1}).n_iters);
    if ~isequaln(pre.(c{1}).node_K0, post.(c{1}).node_K0)
        fprintf('CMP %s trace.node_K0 differs (expected ONLY on fbare nodes)\n', c{1});
    end
end
fprintf('CMP TOTAL MISMATCHES = %d\n', nbad);
