% task2b_part1.m -- complete 40-cell matrix table, full node dump, contiguity analysis,
% q/branch resolution for key nodes. Read-only against task2_matrix_results.mat.
% NOTE: MATLAB's run() changes the current folder to this script's own folder before
% executing -- ALL paths below are therefore absolute, not relative, to avoid the
% double-nesting bug that produces (repo)/.superpowers/sdd/.superpowers/sdd/... .
REPO = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
outdir = fullfile(REPO, '.superpowers', 'sdd', 'task2b_out');
if ~exist(outdir, 'dir'), mkdir(outdir); end

tic;
S = load(fullfile(REPO, '.superpowers', 'sdd', 'task2_matrix_results.mat'), 'results');
results = S.results;
fprintf('loaded %d cells in %.1f s\n', numel(results), toc);

cfg_ids = {results.cfg_id};

%% ---------- full node dump CSV ----------
fid = fopen(fullfile(outdir, 'full_node_dump.csv'), 'w');
fprintf(fid, 'cell_idx,cfg_id,group,variant,field_T,node_id,h,phase,seed_kind,accepted,term_reason,D_uni,Dq_min,Dq_max,s,class\n');
for k = 1:numel(results)
    o = results(k).out;
    for j = 1:numel(o.nodes)
        n = o.nodes(j);
        fprintf(fid, '%d,%s,%s,%s,%.6f,%d,%.8g,%s,%s,%d,%s,%.8g,%.8g,%.8g,%.8g,%s\n', ...
            k, results(k).cfg_id, results(k).group, results(k).variant, results(k).field_T, ...
            n.id, n.h, n.phase, n.seed_kind, n.accepted, n.term_reason, n.D_uni, n.Dq_min, n.Dq_max, n.s, n.class);
    end
end
fclose(fid);
fprintf('wrote full_node_dump.csv\n');

%% ---------- complete matrix markdown table ----------
fid = fopen(fullfile(outdir, 'frag_complete_matrix.md'), 'w');
fprintf(fid, '| # | cfg_id | group | variant | field_T (T) | hmf_status | n_nodes | n_accepted | stable | marginal | unstable | unconverged | s range (accepted) | replay_mismatch |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(results)
    o = results(k).out;
    su = o.summary;
    acc_mask = false(1, numel(o.nodes));
    if ~isempty(o.nodes), acc_mask = [o.nodes.accepted]; end
    if any(acc_mask)
        s_acc = [o.nodes(acc_mask).s];
        s_range = sprintf('[%.4g, %.4g]', min(s_acc), max(s_acc));
    else
        s_range = 'n/a (0 accepted)';
    end
    fprintf(fid, '| %d | %s | %s | %s | %.6f | %s | %d | %d | %d | %d | %d | %d | %s | %d |\n', ...
        k, results(k).cfg_id, results(k).group, results(k).variant, results(k).field_T, ...
        su.hmf_status, su.n_nodes, su.n_accepted, su.n_stable, su.n_marginal, su.n_unstable, ...
        su.n_unconverged, s_range, su.replay_mismatch_any);
end
fclose(fid);
fprintf('wrote frag_complete_matrix.md\n');

%% ---------- global checks: replay_mismatch, class consistency, id==index, h monotonic in id order ----------
fid = fopen(fullfile(outdir, 'frag_global_checks.md'), 'w');
n_mismatch_replay = 0;
n_class_bad = 0;
n_id_bad = 0;
n_nonmonotonic_cells = 0;
n_nodes_total = 0;
for k = 1:numel(results)
    o = results(k).out;
    if o.summary.replay_mismatch_any, n_mismatch_replay = n_mismatch_replay + 1; end
    hprev = -Inf; mono = true;
    for j = 1:numel(o.nodes)
        n = o.nodes(j);
        n_nodes_total = n_nodes_total + 1;
        should_unconv = ~n.accepted || ~isfinite(n.D_uni) || ~isfinite(n.Dq_min);
        is_unconv = strcmp(n.class, 'unconverged');
        if should_unconv ~= is_unconv, n_class_bad = n_class_bad + 1; end
        if n.accepted && isfinite(n.D_uni) && isfinite(n.Dq_min)
            se = min(n.D_uni, n.Dq_min);
            if abs(se - n.s) > 1e-12, n_class_bad = n_class_bad + 1; end
        end
        if n.id ~= j, n_id_bad = n_id_bad + 1; end
        if n.h < hprev - 1e-15, mono = false; end
        hprev = n.h;
    end
    if ~mono, n_nonmonotonic_cells = n_nonmonotonic_cells + 1; end
end
fprintf(fid, '- Total cells: %d\n', numel(results));
fprintf(fid, '- Total nodes across all cells: %d\n', n_nodes_total);
fprintf(fid, '- Cells with replay_mismatch_any=true: %d (must be 0)\n', n_mismatch_replay);
fprintf(fid, '- Nodes violating the accepted-only class rule (Sec. A) or s=min(D_uni,Dq_min) identity: %d (must be 0)\n', n_class_bad);
fprintf(fid, '- Nodes where node.id ~= array position j: %d (must be 0 for id-based position matching to be valid)\n', n_id_bad);
fprintf(fid, '- Cells where h is NOT monotonically non-decreasing in id/visit order: %d (out of %d)\n', n_nonmonotonic_cells, numel(results));
fclose(fid);
type(fullfile(outdir, 'frag_global_checks.md'));

%% ---------- h-grid index-alignment check across ALL cells at each field (not just field 1) ----------
fields_all = unique([results.field_T]);
fid = fopen(fullfile(outdir, 'frag_hgrid_check.md'), 'w');
for fF = 1:numel(fields_all)
    F = fields_all(fF);
    idxF = find(abs([results.field_T] - F) < 1e-6);
    fprintf(fid, '\n### field = %.6f T\n\n', F);
    fprintf(fid, '| cfg_id | n_nodes | max|Δh| vs field reference | id==index? |\n|---|---|---|---|\n');
    href = [];
    for jj = 1:numel(idxF)
        k = idxF(jj);
        o = results(k).out;
        if isempty(o.nodes), fprintf(fid, '| %s | 0 | n/a (no nodes) | n/a |\n', results(k).cfg_id); continue; end
        hvec = [o.nodes.h];
        idok = all([o.nodes.id] == (1:numel(o.nodes)));
        if isempty(href)
            href = hvec;
            fprintf(fid, '| %s | %d | REFERENCE | %d |\n', results(k).cfg_id, numel(hvec), idok);
        else
            if numel(hvec) == numel(href)
                md = max(abs(hvec - href));
            else
                md = NaN;
            end
            fprintf(fid, '| %s | %d | %.3g | %d |\n', results(k).cfg_id, numel(hvec), md, idok);
        end
    end
end
fclose(fid);
fprintf('wrote frag_hgrid_check.md\n');
type(fullfile(outdir, 'frag_hgrid_check.md'));

%% ---------- contiguity analysis (h-sorted) for every SWEPT cell ----------
fid = fopen(fullfile(outdir, 'frag_contiguity.md'), 'w');
for k = 1:numel(results)
    o = results(k).out;
    if ~strcmp(o.summary.mode, 'swept') || isempty(o.nodes), continue; end
    hh = [o.nodes.h];
    cc = {o.nodes.class};
    [hs, ord] = sort(hh);
    cs = cc(ord);
    runs = local_runs(cs);
    fprintf(fid, '\n### %s (field=%.6f T)\n\n', results(k).cfg_id, results(k).field_T);
    fprintf(fid, '| run # | class | n_nodes | h_start | h_end |\n|---|---|---|---|---|\n');
    i0 = 1;
    for r = 1:numel(runs)
        i1 = i0 + runs(r).n - 1;
        fprintf(fid, '| %d | %s | %d | %.6g | %.6g |\n', r, runs(r).class, runs(r).n, hs(i0), hs(i1));
        i0 = i1 + 1;
    end
end
fclose(fid);
fprintf('wrote frag_contiguity.md\n');

%% ---------- physical q resolution helpers, key-node q/branch table ----------
% qc_for_field(F): the field's own g1_swept sibling's out.extra.trc.meta.qc (real 16^3/dp30 lattice)
fid = fopen(fullfile(outdir, 'frag_qbranch.md'), 'w');
fprintf(fid, '| cfg_id | most-marginal ACCEPTED node (min s) | h | s | D_uni | Dq_min | closest +Dq (q_idx,branch,Jq,[h k l]) | closest -Dq (q_idx,branch,Jq,[h k l]) |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
for k = 1:numel(results)
    o = results(k).out;
    grp = results(k).group;
    if isempty(o.nodes), continue; end
    if ~strcmp(o.summary.mode, 'swept'), continue; end  % focus on swept profiles for this table
    accmask = [o.nodes.accepted];
    if ~any(accmask), continue; end
    svals = [o.nodes.s];
    svals(~accmask) = Inf;  % only consider accepted
    [smin, jmin] = min(svals);
    node = o.nodes(jmin);

    % resolve qc source
    qc = [];
    switch grp
        case 'G1G2'
            % use this field's own g1_swept sibling's trc.meta.qc (works even if this IS that cell)
            sib_id = sprintf('g1_swept_F%07d', round(results(k).field_T*1e6));
            ksib = find(strcmp(cfg_ids, sib_id), 1);
            if ~isempty(ksib) && isfield(results(ksib).out.extra, 'trc') && ~isempty(results(ksib).out.extra.trc)
                qc = results(ksib).out.extra.trc.meta.qc;
            end
        case {'G3', 'G6'}
            lp = results(k).lattice_provenance;
            if isfield(lp, 'qc'), qc = lp.qc; end
        otherwise
            qc = [];  % G4/G5 synthetic: no physical q
    end

    pos_str = 'n/a'; neg_str = 'n/a';
    if isfinite(node.idx_pos_flat)
        qi = node.qbranch_pos.q_idx; bi = node.qbranch_pos.branch_idx; Jq = node.qbranch_pos.Jq;
        if isnan(qi) && strcmp(grp, {'G3'}) || (isnan(qi) && strcmp(grp,'G6'))
            % resolve via provenance sidecar (G3/G6 native qbranch is NaN by design)
            lp = results(k).lattice_provenance;
            try
                [qi, bi, Jq] = invz_ordered_trace_resolve(struct('is_synthetic', false, ...
                    'Jnu_unflat', lp.Jnu_unflat, 'nq', lp.nq), node.idx_pos_flat);
            catch ME_DBG
                fprintf('DEBUG pos-resolve FAILED for %s node %d: %s : %s\n', results(k).cfg_id, node.id, ME_DBG.identifier, ME_DBG.message);
            end
        end
        if ~isnan(qi) && ~isempty(qc) && qi <= size(qc,1)
            pos_str = sprintf('(%d,%d,%.4g,[%.3f %.3f %.3f])', qi, bi, Jq, qc(qi,1), qc(qi,2), qc(qi,3));
        elseif ~isnan(qi)
            pos_str = sprintf('(%d,%d,%.4g)', qi, bi, Jq);
        end
    end
    if isfinite(node.idx_neg_flat)
        qi = node.qbranch_neg.q_idx; bi = node.qbranch_neg.branch_idx; Jq = node.qbranch_neg.Jq;
        if isnan(qi) && (strcmp(grp,'G3') || strcmp(grp,'G6'))
            lp = results(k).lattice_provenance;
            try
                [qi, bi, Jq] = invz_ordered_trace_resolve(struct('is_synthetic', false, ...
                    'Jnu_unflat', lp.Jnu_unflat, 'nq', lp.nq), node.idx_neg_flat);
            catch
            end
        end
        if ~isnan(qi) && ~isempty(qc) && qi <= size(qc,1)
            neg_str = sprintf('(%d,%d,%.4g,[%.3f %.3f %.3f])', qi, bi, Jq, qc(qi,1), qc(qi,2), qc(qi,3));
        elseif ~isnan(qi)
            neg_str = sprintf('(%d,%d,%.4g)', qi, bi, Jq);
        end
    end
    fprintf(fid, '| %s | node id=%d | %.6g | %.4g | %.4g | %.4g | %s | %s |\n', ...
        results(k).cfg_id, node.id, node.h, smin, node.D_uni, node.Dq_min, pos_str, neg_str);
end
fclose(fid);
fprintf('wrote frag_qbranch.md\n');

fprintf('\nPART1_DONE\n');

% =================================================================================================
function runs = local_runs(cs)
%LOCAL_RUNS Maximal contiguous runs of identical class labels in cs (cell array of char, in
% h-sorted order). runs(r).class, runs(r).n.
runs = struct('class', {}, 'n', {});
if isempty(cs), return; end
cur = cs{1}; n = 1;
for i = 2:numel(cs)
    if strcmp(cs{i}, cur)
        n = n + 1;
    else
        runs(end+1) = struct('class', cur, 'n', n); %#ok<AGROW>
        cur = cs{i}; n = 1;
    end
end
runs(end+1) = struct('class', cur, 'n', n);
end
