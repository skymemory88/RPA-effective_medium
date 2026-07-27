% task2c_verify_hardstop2.m -- re-verify the C3 pairwise breakdown by calling the ACTUAL frozen
% invz_task2_ladder_ok.m tool (not a hand-rolled reimplementation), on the 2-rung
% (baseline vs half_step/dp30) restriction, treating each common-accepted node as one "field"
% slot -- exactly the report's own §3 operationalization, applied to the pairwise subset.
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(fullfile(REPO, 'invz_projected'));
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
cfg_ids = {results.cfg_id};

fields = [1.173192, 2.581023, 3.754215, 2.850000];
ftag = {'F1173192','F2581023','F3754215','F2850000'};

for f = 1:numel(fields)
    kbase = find(strcmp(cfg_ids, sprintf('g1_swept_%s', ftag{f})), 1);
    kh30  = find(strcmp(cfg_ids, sprintf('g6_offhalfstep_dp30_swept_%s', ftag{f})), 1);
    nb = results(kbase).out.nodes; nh = results(kh30).out.nodes;
    common = find([nb.accepted] & [nh.accepted]);

    rung1 = struct('class', {{nb(common).class}}, 's', [nb(common).s], ...
        'D_uni', [nb(common).D_uni], 'Dq_min', [nb(common).Dq_min], 'label', 'baseline');
    rung2 = struct('class', {{nh(common).class}}, 's', [nh(common).s], ...
        'D_uni', [nh(common).D_uni], 'Dq_min', [nh(common).Dq_min], 'label', 'half_step_dp30');
    ladder = [rung1, rung2];

    [resolved, detail] = invz_task2_ladder_ok(ladder);
    n_field_pass = sum(~detail.class_disagree & ~detail.numeric_disagree);
    fprintf(['field=%s  n_common_accepted=%d  resolved(ALL)=%d  n_pass(per-node)=%d  ' ...
        'n_class_disagree=%d  worst.s=%.4g  worst.D_uni=%.4g  worst.Dq_min=%.4g\n'], ...
        ftag{f}, numel(common), resolved, n_field_pass, sum(detail.class_disagree), ...
        detail.worst.s, detail.worst.D_uni, detail.worst.Dq_min);
end
fprintf('\nVERIFY_HARDSTOP2_DONE\n');
