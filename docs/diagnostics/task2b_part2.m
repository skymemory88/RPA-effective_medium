% task2b_part2.m -- Sec. D ladder (invz_task2_ladder_ok), Sec. C agreement (invz_task2_agree),
% Sec. E existence bar, Step-4 density/distribution table. Read-only against
% task2_matrix_results.mat. Uses the FROZEN tools exactly as provided (no reimplementation).
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
outdir = fullfile(REPO, '.superpowers', 'sdd', 'task2b_out');
if ~exist(outdir, 'dir'), mkdir(outdir); end

tic;
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
fprintf('loaded %d cells in %.1f s\n', numel(results), toc);
cfg_ids = {results.cfg_id};

FIELDS = [1.173192, 2.581023, 3.754215, 2.850000];
FTAGS  = {'F1173192', 'F2581023', 'F3754215', 'F2850000'};

%% =================================================================================================
%% Sec. D -- ladder (grid offset x dpRng), per field, 4 rungs: g1(unshifted,dp30 baseline),
%% g6 unshifted/dp40, g6 halfstep/dp30, g6 halfstep/dp40. N-axis = node index (1..34) --
%% documented methodology note (see report body): invz_task2_ladder_ok's N-axis is generically
%% "comparison points sharing identical order across rungs"; here operationalized as the 34
%% index-aligned H_MF nodes WITHIN each field (h-grid confirmed index-aligned across all cells,
%% Part 1 frag_hgrid_check.md), rather than across the 4 physical fields, since each field's own
%% 34-node profile is exactly what Sec. D's mesh-ladder question concerns.
fid = fopen(fullfile(outdir, 'frag_ladder.md'), 'w');
ladder_verdicts = struct('field', {}, 'resolved_full', {}, 'resolved_acceptedonly', {}, ...
    'n_common_accepted', {}, 'class_disagree_idx', {}, 'numeric_disagree_idx', {});

for iF = 1:numel(FIELDS)
    F = FIELDS(iF); tag = FTAGS{iF};
    id_g1  = sprintf('g1_swept_%s', tag);
    id_u40 = sprintf('g6_offunshifted_dp40_swept_%s', tag);
    id_h30 = sprintf('g6_offhalfstep_dp30_swept_%s', tag);
    id_h40 = sprintf('g6_offhalfstep_dp40_swept_%s', tag);
    k_g1  = idxc(results, id_g1);
    k_u40 = idxc(results, id_u40);
    k_h30 = idxc(results, id_h30);
    k_h40 = idxc(results, id_h40);

    rung(1) = mk_rung(results(k_g1).out.nodes,  'unshifted_dp30(baseline)');
    rung(2) = mk_rung(results(k_u40).out.nodes, 'unshifted_dp40');
    rung(3) = mk_rung(results(k_h30).out.nodes, 'halfstep_dp30');
    rung(4) = mk_rung(results(k_h40).out.nodes, 'halfstep_dp40');

    [resolved_full, detail_full] = invz_task2_ladder_ok(rung);

    % accepted-in-ALL-4-rungs restricted view (interpretive breakdown; see report methodology note)
    accmask = true(1, 34);
    for r = 1:4
        accmask = accmask & [rung(r).accepted];
    end
    n_common_acc = nnz(accmask);
    if n_common_acc > 0
        rung_acc = rung;
        for r = 1:4
            rung_acc(r).class  = rung(r).class(accmask);
            rung_acc(r).s      = rung(r).s(accmask);
            rung_acc(r).D_uni  = rung(r).D_uni(accmask);
            rung_acc(r).Dq_min = rung(r).Dq_min(accmask);
        end
        [resolved_acc, detail_acc] = invz_task2_ladder_ok(rung_acc);
    else
        resolved_acc = NaN; detail_acc = struct();
    end

    fprintf(fid, '\n### field = %.6f T (%s)\n\n', F, tag);
    fprintf(fid, '- FULL 34-node ladder (all nodes, accepted or not): resolved=%d ', resolved_full);
    fprintf(fid, '(class_disagree at %d/34 nodes, numeric_disagree at %d/34 nodes)\n', ...
        nnz(detail_full.class_disagree), nnz(detail_full.numeric_disagree));
    fprintf(fid, '  worst margins (positive = failed by this much): s=%.4g, D_uni=%.4g, Dq_min=%.4g\n', ...
        detail_full.worst.s, detail_full.worst.D_uni, detail_full.worst.Dq_min);
    if n_common_acc > 0
        fprintf(fid, '- ACCEPTED-IN-ALL-4-RUNGS restricted view: n_common_accepted=%d/34, resolved=%d ', ...
            n_common_acc, resolved_acc);
        fprintf(fid, '(class_disagree at %d/%d, numeric_disagree at %d/%d)\n', ...
            nnz(detail_acc.class_disagree), n_common_acc, nnz(detail_acc.numeric_disagree), n_common_acc);
        fprintf(fid, '  worst margins: s=%.4g, D_uni=%.4g, Dq_min=%.4g\n', ...
            detail_acc.worst.s, detail_acc.worst.D_uni, detail_acc.worst.Dq_min);
    else
        fprintf(fid, '- ACCEPTED-IN-ALL-4-RUNGS restricted view: n_common_accepted=0 (no comparison possible)\n');
    end

    fprintf(fid, '\n| node j | h(baseline) | class:g1 | class:u40 | class:h30 | class:h40 | class_agree | s:g1 | s:u40 | s:h30 | s:h40 | numeric_agree |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|\n');
    hh = [results(k_g1).out.nodes.h];
    for j = 1:34
        fprintf(fid, '| %d | %.6g | %s | %s | %s | %s | %s | %.4g | %.4g | %.4g | %.4g | %s |\n', ...
            j, hh(j), rung(1).class{j}, rung(2).class{j}, rung(3).class{j}, rung(4).class{j}, ...
            tf2str(~detail_full.class_disagree(j)), rung(1).s(j), rung(2).s(j), rung(3).s(j), rung(4).s(j), ...
            tf2str(~detail_full.numeric_disagree(j)));
    end

    ladder_verdicts(iF) = struct('field', F, 'resolved_full', resolved_full, ...
        'resolved_acceptedonly', resolved_acc, 'n_common_accepted', n_common_acc, ...
        'class_disagree_idx', find(detail_full.class_disagree), ...
        'numeric_disagree_idx', find(detail_full.numeric_disagree));
end
fclose(fid);
fprintf('wrote frag_ladder.md\n');

%% =================================================================================================
%% Sec. C -- cold-vs-multistart, isolated-vs-swept, full-state agreement (invz_task2_agree),
%% restricted to nodes ACCEPTED on BOTH sides of each comparison (testing "reproducibility of
%% an accepted state", per prereg C's own framing -- an unconverged node has no state to
%% reproduce).
fid = fopen(fullfile(outdir, 'frag_agree.md'), 'w');
agree_summary = struct('field', {}, 'cvm_n_pairs', {}, 'cvm_n_pass', {}, 'cvm_n_fail', {}, ...
    'ivs_n_pairs', {}, 'ivs_n_pass', {}, 'ivs_n_fail', {});
for iF = 1:numel(FIELDS)
    F = FIELDS(iF); tag = FTAGS{iF};
    id_g1   = sprintf('g1_swept_%s', tag);
    id_cold = sprintf('g1_isolated_cold_%s', tag);
    id_sd2  = sprintf('g1_isolated_seed2_%s', tag);
    k_g1   = idxc(results, id_g1);
    k_cold = idxc(results, id_cold);
    k_sd2  = idxc(results, id_sd2);
    J0eff  = results(k_g1).out.meta.J0eff;   % same lattice across all three -> single J0eff
    aopts  = struct('J0eff', J0eff);

    n_g1   = results(k_g1).out.nodes;
    n_cold = results(k_cold).out.nodes;
    n_sd2  = results(k_sd2).out.nodes;

    fprintf(fid, '\n### field = %.6f T (%s), J0eff=%.6g\n\n', F, tag, J0eff);

    fprintf(fid, '**Cold vs multi-start** (isolated_cold vs isolated_seed2), nodes accepted on BOTH sides:\n\n');
    fprintf(fid, '| node j | h | accepted(cold) | accepted(seed2) | comparable | agree | worst_q | worst_diff | worst_AbsTol+RelTol\\*max |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
    n_pairs = 0; n_pass = 0; n_fail = 0;
    for j = 1:34
        accA = n_cold(j).accepted; accB = n_sd2(j).accepted;
        comparable = accA && accB;
        row_ok = ''; wq = ''; wdiff = NaN; wtol = NaN;
        if comparable
            n_pairs = n_pairs + 1;
            [ok, detail] = invz_task2_agree(n_cold(j).state, n_sd2(j).state, aopts);
            if ok, n_pass = n_pass + 1; else, n_fail = n_fail + 1; end
            row_ok = tf2str(ok);
            wq = detail.worst;
            if ~isempty(wq)
                wdiff = detail.(wq).diff;
                wtol  = detail.(wq).AbsTol + detail.(wq).RelTol * max(1, abs(wdiff));  % display only
            end
            fprintf(fid, '| %d | %.6g | %d | %d | 1 | %s | %s | %.4g | (Abs=%.3g,Rel=%.1e) |\n', ...
                j, n_cold(j).h, accA, accB, row_ok, wq, wdiff, detail.(safe_field(wq)).AbsTol, ...
                detail.(safe_field(wq)).RelTol);
        else
            fprintf(fid, '| %d | %.6g | %d | %d | 0 | n/a | n/a | n/a | n/a |\n', j, n_cold(j).h, accA, accB);
        end
    end
    fprintf(fid, '\n**Cold-vs-multistart summary @ %s**: matched-accepted pairs=%d, agree=%d, disagree=%d\n', tag, n_pairs, n_pass, n_fail);

    fprintf(fid, '\n**Isolated vs swept** (isolated_cold vs g1_swept), nodes accepted on BOTH sides:\n\n');
    fprintf(fid, '| node j | h | accepted(swept) | accepted(iso_cold) | comparable | agree | worst_q | worst_diff |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
    n_pairs2 = 0; n_pass2 = 0; n_fail2 = 0;
    for j = 1:34
        accA = n_g1(j).accepted; accB = n_cold(j).accepted;
        comparable = accA && accB;
        if comparable
            n_pairs2 = n_pairs2 + 1;
            [ok, detail] = invz_task2_agree(n_g1(j).state, n_cold(j).state, aopts);
            if ok, n_pass2 = n_pass2 + 1; else, n_fail2 = n_fail2 + 1; end
            wq = detail.worst;
            wdiff = NaN; if ~isempty(wq), wdiff = detail.(wq).diff; end
            fprintf(fid, '| %d | %.6g | %d | %d | 1 | %s | %s | %.4g |\n', ...
                j, n_g1(j).h, accA, accB, tf2str(ok), wq, wdiff);
        else
            fprintf(fid, '| %d | %.6g | %d | %d | 0 | n/a | n/a | n/a |\n', j, n_g1(j).h, accA, accB);
        end
    end
    fprintf(fid, '\n**Isolated-vs-swept summary @ %s**: matched-accepted pairs=%d, agree=%d, disagree=%d\n', tag, n_pairs2, n_pass2, n_fail2);

    % also report: nodes where accepted-status ITSELF differs between swept and isolated_cold
    % (the literal Sec. G item-1 question: "does the SAME node converge alone but not in continuation,
    % or vice versa" -- a 3A-relevant signature independent of invz_task2_agree)
    diff_accept = find([n_g1.accepted] ~= [n_cold.accepted]);
    fprintf(fid, '\nNodes where accepted(swept) ~= accepted(isolated_cold): %s\n', mat2str(diff_accept));
    diff_class = find(~strcmp({n_g1.class}, {n_cold.class}));
    fprintf(fid, 'Nodes where class(swept) ~= class(isolated_cold): %s\n', mat2str(diff_class));

    agree_summary(iF) = struct('field', F, 'cvm_n_pairs', n_pairs, 'cvm_n_pass', n_pass, ...
        'cvm_n_fail', n_fail, 'ivs_n_pairs', n_pairs2, 'ivs_n_pass', n_pass2, 'ivs_n_fail', n_fail2);
end
fclose(fid);
fprintf('wrote frag_agree.md\n');

%% =================================================================================================
%% Sec. E -- existence bar. Per field, per node j: stable+accepted in g1_swept AND isolated_cold
%% AND isolated_seed2; invz_task2_agree passes for (cold,seed2) and (swept,cold); AND >=1 of the
%% 3 G6 offset rungs matches (class-identical + Sec.-D numeric tolerance) the g1 baseline AT that
%% node (giving >=2 mesh offsets total, per D).
fid = fopen(fullfile(outdir, 'frag_existence.md'), 'w');
existence_met = false(1, numel(FIELDS));
for iF = 1:numel(FIELDS)
    F = FIELDS(iF); tag = FTAGS{iF};
    id_g1   = sprintf('g1_swept_%s', tag);
    id_cold = sprintf('g1_isolated_cold_%s', tag);
    id_sd2  = sprintf('g1_isolated_seed2_%s', tag);
    id_u40  = sprintf('g6_offunshifted_dp40_swept_%s', tag);
    id_h30  = sprintf('g6_offhalfstep_dp30_swept_%s', tag);
    id_h40  = sprintf('g6_offhalfstep_dp40_swept_%s', tag);
    k_g1 = idxc(results, id_g1); k_cold = idxc(results, id_cold); k_sd2 = idxc(results, id_sd2);
    k_u40 = idxc(results, id_u40); k_h30 = idxc(results, id_h30); k_h40 = idxc(results, id_h40);
    J0eff = results(k_g1).out.meta.J0eff;
    aopts = struct('J0eff', J0eff);

    n_g1 = results(k_g1).out.nodes; n_cold = results(k_cold).out.nodes; n_sd2 = results(k_sd2).out.nodes;
    n_u40 = results(k_u40).out.nodes; n_h30 = results(k_h30).out.nodes; n_h40 = results(k_h40).out.nodes;

    good_nodes = [];
    detail_rows = {};
    for j = 1:34
        stable3 = strcmp(n_g1(j).class, 'stable') && strcmp(n_cold(j).class, 'stable') && strcmp(n_sd2(j).class, 'stable');
        if ~stable3, continue; end
        [ok_cvm, ~] = invz_task2_agree(n_cold(j).state, n_sd2(j).state, aopts);
        [ok_ivs, ~] = invz_task2_agree(n_g1(j).state, n_cold(j).state, aopts);
        seeds_ok = ok_cvm && ok_ivs;
        offs = [n_u40(j), n_h30(j), n_h40(j)];
        offlabels = {'u40', 'h30', 'h40'};
        off_ok_list = false(1, 3);
        for o = 1:3
            r1 = mk_rung1(n_g1(j)); r2 = mk_rung1(offs(o));
            [rok, ~] = invz_task2_ladder_ok([r1, r2]);
            off_ok_list(o) = rok;
        end
        offsets_ok = any(off_ok_list);
        if seeds_ok && offsets_ok
            good_nodes(end+1) = j; %#ok<AGROW>
        end
        detail_rows{end+1} = sprintf('node %d (h=%.6g): stable3=%d seeds_ok=%d(cvm=%d,ivs=%d) offsets_ok=%d [u40=%d,h30=%d,h40=%d]', ...
            j, n_g1(j).h, stable3, seeds_ok, ok_cvm, ok_ivs, offsets_ok, off_ok_list(1), off_ok_list(2), off_ok_list(3)); %#ok<AGROW>
    end
    existence_met(iF) = ~isempty(good_nodes);
    fprintf(fid, '\n### field = %.6f T (%s)\n\n', F, tag);
    fprintf(fid, 'Existence-bar candidate nodes (stable in g1_swept+isolated_cold+isolated_seed2, seeds agree, >=1 offset agrees): %s\n', mat2str(good_nodes));
    fprintf(fid, 'Existence bar MET at this field: %d\n\n', existence_met(iF));
    fprintf(fid, 'All stable-in-all-3-seed-variants nodes examined:\n');
    for r = 1:numel(detail_rows)
        fprintf(fid, '  - %s\n', detail_rows{r});
    end
end
fclose(fid);
n_fields_met = nnz(existence_met);
fprintf('Existence bar met at %d/4 fields\n', n_fields_met);
fprintf('wrote frag_existence.md\n');

%% =================================================================================================
%% Step 4 -- density/distribution discriminator table (G1 vs G3-stride2/4/8 vs G4 vs G5), per field.
fid = fopen(fullfile(outdir, 'frag_density.md'), 'w');
fprintf(fid, '| field | source | n_accepted | stable | marginal | unstable | unconverged | note |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
for iF = 1:numel(FIELDS)
    F = FIELDS(iF); tag = FTAGS{iF};
    rows = {sprintf('g1_swept_%s', tag), 'G1 baseline (full 16384)'; ...
            sprintf('g3_ds2_swept_%s', tag), 'G3 stride2 (density 1/2, 8192)'; ...
            sprintf('g3_ds4_swept_%s', tag), 'G3 stride4 (density 1/4, 4096)'; ...
            sprintf('g3_ds8_swept_%s', tag), 'G3 stride8 (density 1/8, 2048)'; ...
            sprintf('g4_cardsynth_swept_%s', tag), 'G4 cardinality-synth (16384, non-real shape)'; ...
            sprintf('g5_histmatch_swept_%s', tag), 'G5 histmatch-synth (16384, real shape)'};
    for r = 1:size(rows, 1)
        k = find(strcmp(cfg_ids, rows{r,1}), 1);
        if isempty(k), continue; end
        su = results(k).out.summary;
        fprintf(fid, '| %.6f | %s | %d | %d | %d | %d | %d | %s |\n', ...
            F, rows{r,1}, su.n_accepted, su.n_stable, su.n_marginal, su.n_unstable, su.n_unconverged, rows{r,2});
    end
end
fclose(fid);
fprintf('wrote frag_density.md\n');

fprintf('\nPART2_DONE\n');

% =================================================================================================
function k = idxc(results, cfg_id)
k = find(strcmp({results.cfg_id}, cfg_id), 1);
if isempty(k), error('cfg_id not found: %s', cfg_id); end
end

% =================================================================================================
function rung = mk_rung(nodes, label)
rung = struct();
rung.class    = {nodes.class};
rung.s        = [nodes.s];
rung.D_uni    = [nodes.D_uni];
rung.Dq_min   = [nodes.Dq_min];
rung.accepted = [nodes.accepted];
rung.label    = label;
end

% =================================================================================================
function rung = mk_rung1(node)
rung = struct();
rung.class  = {node.class};
rung.s      = node.s;
rung.D_uni  = node.D_uni;
rung.Dq_min = node.Dq_min;
end

% =================================================================================================
function s = tf2str(tf)
if tf, s = 'YES'; else, s = 'NO'; end
end

% =================================================================================================
function f = safe_field(name)
if isempty(name), f = 'Sigma'; else, f = name; end
end
