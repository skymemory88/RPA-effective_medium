% task2c_verify_hardstop.m -- read-only verification of the review's central claim: the
% exact-h pairwise comparison baseline(unshifted/dp30) vs half_step/dp30 fails the frozen
% §D numeric test (AbsTol=1e-6, RelTol=1e-4 on s, D_uni, Dq_min) at every node that is
% checker-accepted on BOTH sides, at all 4 physical fields.
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

fields = [1.173192, 2.581023, 3.754215, 2.850000];
ftag = {'F1173192','F2581023','F3754215','F2850000'};
AbsTol = 1e-6; RelTol = 1e-4;

for f = 1:numel(fields)
    kbase = find(strcmp(cfg_ids, sprintf('g1_swept_%s', ftag{f})), 1);
    kh30  = find(strcmp(cfg_ids, sprintf('g6_offhalfstep_dp30_swept_%s', ftag{f})), 1);
    ob = results(kbase).out; oh = results(kh30).out;
    nb = ob.nodes; nh = oh.nodes;
    assert(numel(nb) == numel(nh));
    maxdh = max(abs([nb.h] - [nh.h]));
    n_common_accepted = 0; n_pass = 0; n_class_disagree = 0;
    worst_margin = -Inf;
    for j = 1:numel(nb)
        if nb(j).accepted && nh(j).accepted
            n_common_accepted = n_common_accepted + 1;
            if ~strcmp(nb(j).class, nh(j).class), n_class_disagree = n_class_disagree + 1; end
            vals_a = [nb(j).s, nb(j).D_uni, nb(j).Dq_min];
            vals_b = [nh(j).s, nh(j).D_uni, nh(j).Dq_min];
            d = abs(vals_a - vals_b);
            tol = AbsTol + RelTol .* max(abs(vals_a), abs(vals_b));
            margin = max(d - tol); % >0 means fails by that amount
            worst_margin = max(worst_margin, margin);
            if all(d <= tol), n_pass = n_pass + 1; end
        end
    end
    fprintf('field=%s  max|Δh|(baseline vs h30)=%.3e  n_common_accepted=%d  n_pass=%d  n_class_disagree=%d  worst_margin=%.4g\n', ...
        ftag{f}, maxdh, n_common_accepted, n_pass, n_class_disagree, worst_margin);
end
fprintf('\nVERIFY_HARDSTOP_DONE\n');
