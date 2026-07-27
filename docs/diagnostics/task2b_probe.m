% task2b_probe.m -- Stage A structural probe of task2_matrix_results.mat (read-only).
% Run: matlab -batch, foreground. Prints structure only; makes no claims used in the report
% without a follow-up direct computation.
tic;
S = load('.superpowers/sdd/task2_matrix_results.mat', 'results');
results = S.results;
fprintf('n cells = %d\n', numel(results));
fprintf('load time = %.1f s\n', toc);

fprintf('\n--- fieldnames(results) ---\n');
disp(fieldnames(results));

fprintf('\n--- cfg_id list (all 40) ---\n');
for k = 1:numel(results)
    fprintf('%2d  %-32s group=%-5s variant=%-20s field_T=%.6f\n', k, results(k).cfg_id, ...
        results(k).group, results(k).variant, results(k).field_T);
end

fprintf('\n--- out.summary per cell ---\n');
for k = 1:numel(results)
    o = results(k).out;
    su = o.summary;
    status = '';
    if isfield(o, 'status'), status = o.status; end
    fprintf('%2d %-32s hmf_status=%-10s n_nodes=%3d n_accepted=%3d stable=%2d marginal=%2d unstable=%2d unconverged=%2d hstar=%.6g mismatch=%d status=%s\n', ...
        k, results(k).cfg_id, su.hmf_status, su.n_nodes, su.n_accepted, su.n_stable, su.n_marginal, ...
        su.n_unstable, su.n_unconverged, su.hstar, su.replay_mismatch_any, status);
end

% Check h-grid identity across cells at same field: pick field 1173192 group
fprintf('\n--- h-grid identity check @ field=1.173192 ---\n');
idxF = find(abs([results.field_T] - 1.173192) < 1e-6);
h_ref = [];
for jj = 1:numel(idxF)
    k = idxF(jj);
    o = results(k).out;
    if isempty(o.nodes)
        fprintf('%-32s NO NODES (failed?)\n', results(k).cfg_id);
        continue;
    end
    hvec = [o.nodes.h];
    if isempty(h_ref)
        h_ref = hvec;
        fprintf('%-32s n=%d h(1:5)=%s  [REFERENCE]\n', results(k).cfg_id, numel(hvec), mat2str(hvec(1:min(5,end)),6));
    else
        same_n = numel(hvec) == numel(h_ref);
        if same_n
            maxdiff = max(abs(hvec - h_ref));
        else
            maxdiff = NaN;
        end
        fprintf('%-32s n=%d same_n=%d maxdiff_vs_ref=%.3g\n', results(k).cfg_id, numel(hvec), same_n, maxdiff);
    end
end

% Check node phase labels + term_reason distribution for one physical cell
fprintf('\n--- g1_swept_F1173192 node detail (all 34) ---\n');
k = find(strcmp({results.cfg_id}, 'g1_swept_F1173192'), 1);
o = results(k).out;
for j = 1:numel(o.nodes)
    n = o.nodes(j);
    fprintf('id=%2d h=%10.6f phase=%-10s seed=%-6s accepted=%d term=%-20s D_uni=%10.4g Dq_min=%10.4g s=%10.4g class=%-12s\n', ...
        n.id, n.h, n.phase, n.seed_kind, n.accepted, n.term_reason, n.D_uni, n.Dq_min, n.s, n.class);
end

fprintf('\n--- accepted-only class consistency spot check (ALL 40 cells) ---\n');
n_bad = 0; n_nodes_total = 0;
for k = 1:numel(results)
    o = results(k).out;
    for j = 1:numel(o.nodes)
        n = o.nodes(j);
        n_nodes_total = n_nodes_total + 1;
        is_unconv = strcmp(n.class, 'unconverged');
        should_be_unconv = ~n.accepted || ~isfinite(n.D_uni) || ~isfinite(n.Dq_min);
        if is_unconv ~= should_be_unconv
            n_bad = n_bad + 1;
            fprintf('MISMATCH cell=%d(%s) node=%d accepted=%d D_uni=%g Dq_min=%g class=%s\n', ...
                k, results(k).cfg_id, j, n.accepted, n.D_uni, n.Dq_min, n.class);
        end
        % also check accepted nodes' class matches s vs +-1e-3 threshold
        if n.accepted && isfinite(n.D_uni) && isfinite(n.Dq_min)
            s_expect = min(n.D_uni, n.Dq_min);
            if abs(s_expect - n.s) > 1e-12
                fprintf('S-MISMATCH cell=%d(%s) node=%d s=%g expect=%g\n', k, results(k).cfg_id, j, n.s, s_expect);
                n_bad = n_bad + 1;
            end
            if s_expect > 1e-3 && ~strcmp(n.class,'stable')
                fprintf('CLASS-MISMATCH(stable) cell=%d node=%d s=%g class=%s\n', k, j, s_expect, n.class); n_bad=n_bad+1;
            elseif abs(s_expect) <= 1e-3 && ~strcmp(n.class,'marginal')
                fprintf('CLASS-MISMATCH(marginal) cell=%d node=%d s=%g class=%s\n', k, j, s_expect, n.class); n_bad=n_bad+1;
            elseif s_expect < -1e-3 && ~strcmp(n.class,'unstable')
                fprintf('CLASS-MISMATCH(unstable) cell=%d node=%d s=%g class=%s\n', k, j, s_expect, n.class); n_bad=n_bad+1;
            end
        end
    end
end
fprintf('total nodes=%d, total mismatches=%d\n', n_nodes_total, n_bad);

fprintf('\n--- replay_mismatch_any across all 40 ---\n');
any_true = false;
for k = 1:numel(results)
    if results(k).out.summary.replay_mismatch_any
        any_true = true;
        fprintf('TRUE for cell %d (%s)\n', k, results(k).cfg_id);
    end
end
if ~any_true, fprintf('all false (good)\n'); end

fprintf('\n--- lattice_provenance presence ---\n');
for k = 1:numel(results)
    lp = results(k).lattice_provenance;
    has_qc = isfield(lp,'qc') && ~isempty(lp.qc);
    fprintf('%2d %-32s has_provenance=%d\n', k, results(k).cfg_id, has_qc);
end

fprintf('\nDONE\n');
