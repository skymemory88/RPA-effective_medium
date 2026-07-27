% Task 18 Step 3/5 support: merge every chunked invz_gate0_report() call into ONE final rep,
% call invz_gate0_aggregate ONCE on the union (the authoritative verdict, not any chunk's own
% partial one), and print every table the brief/verdict doc requires. Run via
% matlab -batch "run('/abs/path/task18_merge_report.m')" -- absolute addpath below because
% `run` executes with THIS script's own directory as context (bit Task 16's implementer).
ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
SDD  = fullfile(ROOT, '.superpowers', 'sdd');
addpath(fullfile(ROOT, 'invz_projected'));
addpath(fullfile(ROOT, 'invz_common'));
addpath(ROOT);

ordered_fields = [0.05 0.25 0.5 1 2 2.5 2.9 3.0];
nH_list        = [33 65 129];
pm_fields      = [3.1 3.5];

chunk_files = {'task18_chunk_B0p05.mat', 'task18_chunk_B0p25.mat', 'task18_chunk_B0p5.mat', ...
               'task18_calib_1p0_nH33.mat', 'task18_calib_1p0_nH65_129.mat', ...
               'task18_chunk_B2p0.mat', 'task18_chunk_B2p5.mat', 'task18_chunk_B2p9.mat', ...
               'task18_chunk_B3p0.mat'};

all_ordered = [];
for k = 1:numel(chunk_files)
    L = load(fullfile(SDD, chunk_files{k}));
    if isempty(all_ordered)
        all_ordered = L.rep.ordered;
    else
        all_ordered = [all_ordered, L.rep.ordered]; %#ok<AGROW>
    end
end

Lpm = load(fullfile(SDD, 'task18_chunk_PM_B0.mat'));
all_pm = Lpm.rep.pm;
b0     = Lpm.rep.b0;
digest_pm_chunk = Lpm.rep.digest;

fprintf('=== coverage check ===\n');
fprintf('numel(all_ordered) = %d (expect %d = %d fields x %d nH)\n', numel(all_ordered), ...
    numel(ordered_fields)*numel(nH_list), numel(ordered_fields), numel(nH_list));
fprintf('numel(all_pm) = %d (expect %d)\n', numel(all_pm), numel(pm_fields));
fprintf('digest (from PM/B0 chunk) = %s\n', digest_pm_chunk);

agg = invz_gate0_aggregate(all_ordered, all_pm, ordered_fields, nH_list);

fprintf('\n=== Table 1: ordered (B, nH) rows ===\n');
fprintf(['%6s %4s %-20s %12s %10s %10s %10s %6s %6s %5s %5s %5s %5s %5s %12s %10s %12s %12s %-12s %-12s\n'], ...
    'B', 'nH', 'status', 'hstar', 'crit*', 'D_uni*', 'Dq_min*', 'D_tol*', 'n_nod', 'n_ok', 'n_mod', ...
    'n_dd', 'n_sf', 'n_unrec', 'min_refmar', 'maxomitL', 'int_Sig0', 'int_rm1', 'g17_r', 'g17_crit');
for k = 1:numel(all_ordered)
    r = all_ordered(k);
    n_ok  = nnz(strcmp(r.bucket, 'ok'));
    n_mod = nnz(strcmp(r.bucket, 'medium_out_of_domain'));
    n_dd  = nnz(strcmp(r.bucket, 'degenerate_doublet'));
    n_sf  = nnz(strcmp(r.bucket, 'solver_failed'));
    n_un  = nnz(strcmp(r.bucket, 'unrecognized'));
    fprintf(['%6g %4d %-20s %12.6g %10.4g %10.4g %10.4g %10.3g %6d %6d %5d %5d %5d %5d %12.5g %10.5g %12.6g %12.6g %-12s %-12s\n'], ...
        r.B, r.nH, r.status, r.hstar, r.crit_star, r.D_uni_star, r.Dq_min_star, r.D_tol_star, ...
        r.n_nodes, n_ok, n_mod, n_dd, n_sf, n_un, r.min_ref_margin, r.max_omit_ledger, ...
        r.int_Sigma0, r.int_r_minus_1, r.g17_r.status, r.g17_crit.status);
end

fprintf('\n=== Table 2: PM controls ===\n');
fprintf('%6s %9s %12s %10s %-20s\n', 'B', 'converged', 'crit', 'crit_tol', 'medium_status');
for k = 1:numel(all_pm)
    r = all_pm(k);
    fprintf('%6g %9d %12.6g %10.3g %-20s\n', r.B, r.converged, r.crit, r.crit_tol, r.medium_status);
end

fprintf('\n=== Table 3: exact B=0 hard-domain control ===\n');
fprintf('status=%s hstar=%.6g n_nodes=%d n_accounted=%d expected_degenerate=%d\n', ...
    b0.status, b0.hstar, b0.n_nodes, b0.n_accounted, b0.expected_degenerate);

fprintf('\n=== Table 4: G5 path integrals -- raw values at each nH (see Table 1 too) ===\n');
fprintf('%6s %6s %16s %16s\n', 'B', 'nH', 'int_Sigma0', 'int_r_minus_1');
for iB = 1:numel(ordered_fields)
    B = ordered_fields(iB);
    rows_B = all_ordered([all_ordered.B] == B);
    [~, order] = sort([rows_B.nH]);
    rows_B = rows_B(order);
    for kk = 1:3
        fprintf('%6g %6d %16.8g %16.8g\n', B, rows_B(kk).nH, rows_B(kk).int_Sigma0, rows_B(kk).int_r_minus_1);
    end
end

fprintf(['\n=== Table 4b: G5 convergence diagnostic -- |I_fine-I_prev| vs I_atol+1e-3*max(|I_fine|,' ...
         '|I_prev|) ===\n']);
fprintf(['(33->65 is the APPROACH diagnostic only, prereg SS5 -- NOT gated. 65->129 is the ' ...
         'frozen G5 criterion.)\n']);
fprintf('%6s | %14s %14s %10s | %14s %14s %10s\n', 'B', 'd33->65(Sig0)', 'd33->65(rm1)', '(diag)', ...
    'd65->129(Sig0)', 'd65->129(rm1)', 'g5_pass');
I_atol = 1e-10;
for iB = 1:numel(ordered_fields)
    B = ordered_fields(iB);
    rows_B = all_ordered([all_ordered.B] == B);
    [~, order] = sort([rows_B.nH]);
    rows_B = rows_B(order);
    S = [rows_B.int_Sigma0];  Rm1 = [rows_B.int_r_minus_1];
    d1S = abs(S(2) - S(1));    d1R = abs(Rm1(2) - Rm1(1));
    d2S = abs(S(3) - S(2));    d2R = abs(Rm1(3) - Rm1(2));
    tolS = I_atol + 1e-3*max(abs(S(3)), abs(S(2)));
    tolR = I_atol + 1e-3*max(abs(Rm1(3)), abs(Rm1(2)));
    okS = isfinite(d2S) && d2S <= tolS;
    okR = isfinite(d2R) && d2R <= tolR;
    fprintf('%6g | %14.6g %14.6g %10s | %14.6g %14.6g %10s   (tolS=%.3g tolR=%.3g)\n', ...
        B, d1S, d1R, '--', d2S, d2R, mat2str(okS && okR), tolS, tolR);
end

fprintf('\n=== Table 5: G17 (d) per-field nH=65 vs nH=129 ===\n');
fprintf('%6s | %10s %12s %12s | %10s %12s %12s | %10s %12s %12s | %10s %12s %12s\n', ...
    'B', 'r@65', 'jump@65', '', 'r@129', 'jump@129', '', 'crit@65', 'jump@65', '', 'crit@129', 'jump@129', '');
for iB = 1:numel(ordered_fields)
    B = ordered_fields(iB);
    r65  = all_ordered([all_ordered.B] == B & [all_ordered.nH] == 65);
    r129 = all_ordered([all_ordered.B] == B & [all_ordered.nH] == 129);
    fprintf('%6g | %10s %12.4g %12s | %10s %12.4g %12s | %10s %12.4g %12s | %10s %12.4g %12s\n', ...
        B, r65.g17_r.status, r65.g17_r.max_jump, '', r129.g17_r.status, r129.g17_r.max_jump, '', ...
        r65.g17_crit.status, r65.g17_crit.max_jump, '', r129.g17_crit.status, r129.g17_crit.max_jump, '');
end

fprintf('\n=== Verdict ===\n');
fprintf('fail_a=%d fail_b=%d fail_c=%d fail_d=%d fail_e=%d\n', agg.fail_a, agg.fail_b, agg.fail_c, agg.fail_d, agg.fail_e);
fprintf('rep.pass = %d\n', agg.pass);
fprintf('detail.fail_a_rows   = %s\n', mat2str(agg.detail.fail_a_rows));
fprintf('detail.fail_b_rows   = %s\n', mat2str(agg.detail.fail_b_rows));
fprintf('detail.fail_c_rows   = %s\n', mat2str(agg.detail.fail_c_rows));
fprintf('detail.fail_d_fields = %s\n', mat2str(cell2mat(agg.detail.fail_d_fields)));
fprintf('detail.fail_e_rows   = %s\n', mat2str(agg.detail.fail_e_rows));
fprintf('detail.fail_e_pm     = %s\n', mat2str(agg.detail.fail_e_pm));
fprintf('detail.missing_rows: %d entries\n', numel(agg.detail.missing_rows));
for k = 1:numel(agg.detail.missing_rows), fprintf('  %s\n', agg.detail.missing_rows{k}); end

% rep.g5 / rep.g17: same reduction invz_gate0_report.m's local_build_g5/local_build_g17 apply
% (those are file-local and not externally callable, so replicated here verbatim over the
% UNIONED rows -- this merge script is a support artifact, not the committed deliverable, whose
% own copy of this logic was already validated byte-for-byte against this replica on B=2 below).
I_atol = 1e-10;
g5 = [];
for iB = 1:numel(ordered_fields)
    B = ordered_fields(iB);
    rB = all_ordered(abs([all_ordered.B] - B) < 1e-9);
    g = struct('B', B, 'int_Sigma0_33', NaN, 'int_Sigma0_65', NaN, 'int_Sigma0_129', NaN, ...
        'int_r_minus_1_33', NaN, 'int_r_minus_1_65', NaN, 'int_r_minus_1_129', NaN, ...
        'd_33_65_Sigma0', NaN, 'd_33_65_r_minus_1', NaN, ...
        'd_65_129_Sigma0', NaN, 'd_65_129_r_minus_1', NaN, ...
        'tol_65_129_Sigma0', NaN, 'tol_65_129_r_minus_1', NaN, ...
        'pass_Sigma0', false, 'pass_r_minus_1', false);
    for k = 1:numel(rB)
        switch rB(k).nH
            case 33,  g.int_Sigma0_33  = rB(k).int_Sigma0;  g.int_r_minus_1_33  = rB(k).int_r_minus_1;
            case 65,  g.int_Sigma0_65  = rB(k).int_Sigma0;  g.int_r_minus_1_65  = rB(k).int_r_minus_1;
            case 129, g.int_Sigma0_129 = rB(k).int_Sigma0;  g.int_r_minus_1_129 = rB(k).int_r_minus_1;
        end
    end
    g.d_33_65_Sigma0     = abs(g.int_Sigma0_65  - g.int_Sigma0_33);
    g.d_33_65_r_minus_1  = abs(g.int_r_minus_1_65  - g.int_r_minus_1_33);
    g.d_65_129_Sigma0    = abs(g.int_Sigma0_129 - g.int_Sigma0_65);
    g.d_65_129_r_minus_1 = abs(g.int_r_minus_1_129 - g.int_r_minus_1_65);
    g.tol_65_129_Sigma0    = I_atol + 1e-3*max(abs(g.int_Sigma0_129),    abs(g.int_Sigma0_65));
    g.tol_65_129_r_minus_1 = I_atol + 1e-3*max(abs(g.int_r_minus_1_129), abs(g.int_r_minus_1_65));
    g.pass_Sigma0    = isfinite(g.d_65_129_Sigma0)    && g.d_65_129_Sigma0    <= g.tol_65_129_Sigma0;
    g.pass_r_minus_1 = isfinite(g.d_65_129_r_minus_1) && g.d_65_129_r_minus_1 <= g.tol_65_129_r_minus_1;
    if isempty(g5), g5 = g; else, g5(end+1) = g; end %#ok<AGROW>
end

fail_d_fields = agg.detail.fail_d_fields;
if isempty(fail_d_fields), fail_d_set = []; else, fail_d_set = cell2mat(fail_d_fields); end
g17 = [];
for iB = 1:numel(ordered_fields)
    B = ordered_fields(iB);
    r65  = all_ordered(abs([all_ordered.B] - B) < 1e-9 & [all_ordered.nH] == 65);
    r129 = all_ordered(abs([all_ordered.B] - B) < 1e-9 & [all_ordered.nH] == 129);
    g = struct('B', B, 'g17_r_65', [], 'g17_r_129', [], 'g17_crit_65', [], 'g17_crit_129', [], ...
        'pass', ~any(abs(fail_d_set - B) < 1e-9));
    if ~isempty(r65),  g.g17_r_65 = r65.g17_r;    g.g17_crit_65 = r65.g17_crit;  end
    if ~isempty(r129), g.g17_r_129 = r129.g17_r;  g.g17_crit_129 = r129.g17_crit; end
    if isempty(g17), g17 = g; else, g17(end+1) = g; end %#ok<AGROW>
end

rep = struct();
rep.ordered = all_ordered;  rep.pm = all_pm;  rep.b0 = b0;
rep.g5 = g5;  rep.g17 = g17;
rep.fail_a = agg.fail_a;  rep.fail_b = agg.fail_b;  rep.fail_c = agg.fail_c;
rep.fail_d = agg.fail_d;  rep.fail_e = agg.fail_e;  rep.pass = agg.pass;
rep.detail = agg.detail;
rep.digest = digest_pm_chunk;
save(fullfile(SDD, 'task18_final_rep.mat'), 'rep');
fprintf('\nSAVED final merged rep -> %s\n', fullfile(SDD, 'task18_final_rep.mat'));

fprintf('\n=== cross-check: rep.g5/rep.g17 (B=2) vs the live driver (task18_verify_g5g17.m) ===\n');
g5_2  = g5([g5.B] == 2);
g17_2 = g17([g17.B] == 2);
fprintf('g5(B=2):  pass_Sigma0=%d pass_r_minus_1=%d d6529S=%.6g tol=%.6g d6529R=%.6g tol=%.6g\n', ...
    g5_2.pass_Sigma0, g5_2.pass_r_minus_1, g5_2.d_65_129_Sigma0, g5_2.tol_65_129_Sigma0, ...
    g5_2.d_65_129_r_minus_1, g5_2.tol_65_129_r_minus_1);
fprintf('g17(B=2): pass=%d r65=%s/%.6g r129=%s/%.6g crit65=%s/%.6g crit129=%s/%.6g\n', ...
    g17_2.pass, g17_2.g17_r_65.status, g17_2.g17_r_65.max_jump, g17_2.g17_r_129.status, ...
    g17_2.g17_r_129.max_jump, g17_2.g17_crit_65.status, g17_2.g17_crit_65.max_jump, ...
    g17_2.g17_crit_129.status, g17_2.g17_crit_129.max_jump);
