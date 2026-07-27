function tests = test_invz_gate0_aggregate
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% ============================================================================================
% Synthetic fixture builders. NO physics: these are hand-built rows in the exact schema
% invz_gate0_report reduces trc/prof into (see invz_gate0_aggregate.m's header), so every one
% of (a)-(e) can be forced independently and cheaply -- the whole point of Task 18 brief Step 2's
% mandatory aggregation fixture (no Boolean merely printed while never load-bearing).
% ============================================================================================
function rows = good_rows(B)
g17ok = struct('status', 'no_crossing', 'max_jump', NaN, 'crossings', struct([]));
nHs = [33 65 129];
rows = repmat(local_blank_row(), 1, 3);
for k = 1:3
    r = local_blank_row();
    r.B = B;  r.nH = nHs(k);
    r.status = 'ok';  r.hstar = 0.4;
    r.crit_star = 0.5;  r.D_uni_star = 0.5;  r.Dq_min_star = 0.5;
    r.crit_tol = 1e-6;  r.D_tol_star = 1e-6;  r.Dq_tol_star = 1e-6;
    r.n_nodes = 5;
    r.phase         = {'predictor', 'sweep', 'sweep', 'sweep', 'root'};
    r.bucket        = {'ok', 'ok', 'ok', 'ok', 'ok'};
    r.medium_status = {'ok', 'ok', 'ok', 'ok', 'ok'};
    r.omit_max      = [0.01 0.02 0.02 0.02 0.01];
    r.ok_final      = true(1, 5);
    r.h             = [0 0.1 0.2 0.3 0.4];
    r.g17_r = g17ok;  r.g17_crit = g17ok;  r.g17_anomaly = false;
    rows(k) = r;
end
end

function r = local_blank_row()
r = struct('B', NaN, 'nH', NaN, 'status', '', 'hstar', NaN, 'crit_star', NaN, ...
    'D_uni_star', NaN, 'Dq_min_star', NaN, 'crit_tol', NaN, 'D_tol_star', NaN, ...
    'Dq_tol_star', NaN, 'n_nodes', 0, 'phase', {{}}, 'bucket', {{}}, ...
    'medium_status', {{}}, 'omit_max', [], 'ok_final', [], 'h', [], ...
    'g17_r', [], 'g17_crit', [], 'g17_anomaly', false);
end

function pm = good_pm(Bs)
pm = repmat(struct('B', NaN, 'converged', true, 'crit', 0.5, 'crit_tol', 1e-6, ...
    'medium_status', 'ok'), 1, numel(Bs));
for k = 1:numel(Bs), pm(k).B = Bs(k); end
end

function verify_only(testCase, rep, which)
%VERIFY_ONLY Assert exactly the named clause(s) fired and no other, then rep.pass is false.
names = {'a', 'b', 'c', 'd', 'e'};
for k = 1:numel(names)
    fld = ['fail_' names{k}];
    expected = any(strcmp(which, names{k}));
    verifyEqual(testCase, rep.(fld), expected, ...
        sprintf('%s: expected %d, got %d', fld, expected, rep.(fld)));
end
verifyFalse(testCase, rep.pass);
end

% ============================================================================================
% Sanity: the baseline fixture itself passes clean -- otherwise every "isolated mutation"
% test below would be meaningless.
% ============================================================================================
function test_baseline_passes(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verifyTrue(testCase, rep.pass);
verifyFalse(testCase, rep.fail_a);  verifyFalse(testCase, rep.fail_b);
verifyFalse(testCase, rep.fail_c);  verifyFalse(testCase, rep.fail_d);
verifyFalse(testCase, rep.fail_e);
end

function test_exact_pass_formula(testCase)
% rep.pass must be EXACTLY ~(a|b|c|d|e) -- pin the formula itself, not just spot values.
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verifyEqual(testCase, rep.pass, ...
    ~(rep.fail_a || rep.fail_b || rep.fail_c || rep.fail_d || rep.fail_e));
end

% ============================================================================================
% (a): a non-'ok' reference denominator status on a required solved-path node.
% ============================================================================================
function test_fail_a_from_bad_medium_status(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).medium_status{2} = 'ref_denom_nonpositive';
rows(1).bucket{2} = 'medium_out_of_domain';   % realistic pairing; bucket stays RECOGNIZED
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'a'});
verifyEqual(testCase, rep.detail.fail_a_rows, 1);
end

function test_not_applicable_medium_status_does_not_fire_a(testCase)
% 'not_applicable' (e.g. a degenerate_doublet node that never reached reference evaluation)
% is excluded from (a) -- this is the same rule that keeps the B=0 hard-domain control from
% firing (a) when invz_gate0_report evaluates it (it is kept entirely out of the rows passed
% here, but the underlying medium_status semantics are the same and must be tested).
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).medium_status{2} = 'not_applicable';
rows(1).bucket{2} = 'degenerate_doublet';
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verifyFalse(testCase, rep.fail_a);
end

% ============================================================================================
% (b): coverage identity -- unrecognized bucket, missing predictor, missing root entry, and a
% dropped chunk (coverage-completeness across the expected field/nH grid).
% ============================================================================================
function test_fail_b_from_unrecognized_bucket(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).bucket{3} = 'unrecognized';   % an empty/unknown term_reason must not be silently folded
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'b'});
end

function test_fail_b_from_missing_predictor(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).phase{1} = 'sweep';   % predictor phase entry absent
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'b'});
end

function test_fail_b_from_missing_root_entry(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).phase{5} = 'sweep';   % hstar is finite but no 'root' phase entry backs it up
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'b'});
end

function test_fail_b_from_dropped_chunk(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(3) = [];   % simulate a lost nH=129 chunk during aggregation
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verifyTrue(testCase, rep.fail_b);
verifyFalse(testCase, rep.fail_a);  verifyFalse(testCase, rep.fail_c);
verifyFalse(testCase, rep.fail_d);  verifyFalse(testCase, rep.fail_e);
verifyNotEmpty(testCase, rep.detail.missing_rows);
end

% ============================================================================================
% (c): omitted-order ratio -- exceeding omit_promote, and a missing/non-finite ratio at an
% 'ok'-labelled node (never dropped by isfinite filtering).
% ============================================================================================
function test_fail_c_from_exceeding_omit_promote(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).omit_max(2) = 0.50;   % >> frozen omit_promote = 0.10, node stays medium_status='ok'
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'c'});
end

function test_fail_c_from_nonfinite_omit_at_ok_node(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).omit_max(2) = NaN;   % node still labelled 'ok' -- NOT dropped by isfinite filtering
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'c'});
end

function test_fail_c_absent_when_bad_node_already_excluded_by_a(testCase)
% A node whose medium is already non-'ok' contributes no omit_max opinion to (c) -- (a) alone
% covers it; (c) must not double-fire from the same node once it is no longer 'ok'-labelled.
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).medium_status{2} = 'ref_denom_nonpositive';
rows(1).bucket{2} = 'medium_out_of_domain';
rows(1).omit_max(2) = NaN;   % consistent with a domain-failed node (see invz_medium_moment_closure)
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'a'});
end

% ============================================================================================
% (d): G17 actual-path crossing continuity -- unresolved/jump-exceeded status, an increasing
% jump from nH=65 to 129, and the "finite integrands" anomaly flag.
% ============================================================================================
function test_fail_d_from_jump_exceeded_status(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(3).g17_r.status = 'jump_exceeded';  rows(3).g17_r.max_jump = 5.0;   % the nH=129 row
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'d'});
verifyEqual(testCase, rep.detail.fail_d_fields, {1.0});
end

function test_fail_d_from_unresolved_status(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(2).g17_crit.status = 'unresolved';   % the nH=65 row
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'d'});
end

function test_fail_d_from_increasing_jump(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(2).g17_r.status = 'ok';  rows(2).g17_r.max_jump = 1e-5;    % nH=65
rows(3).g17_r.status = 'ok';  rows(3).g17_r.max_jump = 1e-4;    % nH=129, INCREASED -> forbidden
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'d'});
end

function test_fail_d_from_finite_integrand_anomaly(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(3).g17_anomaly = true;   % an 'ok'-labelled node with non-finite r/crit/gstat_local_denom
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'d'});
end

function test_pass_when_jump_decreases_or_holds(testCase)
% the "must not increase" rule must not misfire on a DEcreasing or held jump.
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(2).g17_crit.status = 'ok';  rows(2).g17_crit.max_jump = 1e-4;
rows(3).g17_crit.status = 'ok';  rows(3).g17_crit.max_jump = 1e-4;   % held exactly
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verifyTrue(testCase, rep.pass);
end

% ============================================================================================
% (e): endpoint stability -- the ordered half and the PM-control half, independently.
% ============================================================================================
function test_fail_e_from_bad_endpoint(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(3).D_uni_star = 0;   % <= D_tol_star = 1e-6
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'e'});
end

function test_fail_e_from_nonok_status(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(2).status = 'medium_out_of_domain';   % hstar/crit_star stay finite, but status != 'ok'
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'e'});
end

function test_fail_e_from_pm_control(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
pm(2).converged = false;
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'e'});
verifyEqual(testCase, rep.detail.fail_e_pm, 2);
end

function test_fail_e_from_pm_medium_status(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
pm(1).medium_status = 'ref_denom_small';
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'e'});
end

% ============================================================================================
% Combined: two independent clauses firing at once must both be reported (no short-circuit).
% ============================================================================================
function test_two_clauses_fire_together(testCase)
rows = good_rows(1.0);  pm = good_pm([3.1 3.5]);
rows(1).omit_max(2) = 0.50;         % (c)
rows(3).D_uni_star = 0;             % (e)
rep = invz_gate0_aggregate(rows, pm, 1.0, [33 65 129]);
verify_only(testCase, rep, {'c', 'e'});
end
