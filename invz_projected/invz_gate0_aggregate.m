function rep = invz_gate0_aggregate(ordered_rows, pm_rows, expected_fields, expected_nH)
%INVZ_GATE0_AGGREGATE Pure Gate-0 pass/fail reduction (Task 18 Step 2). The frozen
% predicate -- the "Gate-0 negative-outcome rule" -- is implemented below and summarized in
% invzp_convg_diagnosis.md Section 8.2. G = -chi (meV^-1), ferromagnetic positive J.
%
% This is the ONLY place the five clause Booleans (a)-(e) are decided. invz_gate0_report.m
% reduces raw invz_hmf_ordered/invz_solve_point output (trc/prof) into the fixed-schema ROW
% structs consumed here; this function contains NO solving and NO I/O, so every clause is
% independently forceable from a synthetic fixture (Task 18 brief Step 2's mandatory
% aggregation unit test) without running any physics.
%
% ordered_rows (struct array, one entry per required ordered (B, nH) combination; may be 0x0):
%   B, nH                       field (T) and profile resolution
%   status                      prof.status (char)
%   hstar, crit_star, D_uni_star, Dq_min_star     the accepted root's own record (NaN if none)
%   crit_tol, D_tol_star, Dq_tol_star              the tolerances AT the root state (prereg SS1;
%                                                   D_tol/Dq_tol are state-dependent, so the
%                                                   CALLER computes them from Gstat_star/Jscale)
%   n_nodes                     numel(trc.nodes) for this (B,nH)
%   bucket                      1xn_nodes cellstr, ALREADY mapped by the caller's frozen
%                                term_reason table into {'ok','medium_out_of_domain',
%                                'degenerate_doublet','solver_failed','unrecognized'} -- an
%                                'unrecognized' entry is what an empty/unknown term_reason maps
%                                to; it is never silently folded into another bucket.
%   medium_status                1xn_nodes cellstr, trc.nodes.medium_status verbatim
%   omit_max                     1xn_nodes double, trc.nodes.omit_max verbatim
%   phase                        1xn_nodes cellstr, trc.nodes.phase verbatim
%   ok_final                     1xn_nodes logical, trc.nodes.ok_final verbatim
%   h                            1xn_nodes double, trc.nodes.h verbatim
%   g17_r, g17_crit               invz_pole_continuity(...) result structs (.status/.max_jump),
%                                computed by the caller on THIS (B,nH)'s final adapted
%                                prof.hgrid, after dropping any non-finite (h,d,y) triple
%   g17_anomaly                  logical: true iff any node with prof.medium_status=='ok' at
%                                this (B,nH) had a non-finite r, crit, or gstat_local_denom (the
%                                "require finite integrands everywhere" sub-clause of (d) -- an
%                                'ok'-labelled node must never be silently dropped from G17
%                                without explanation)
% pm_rows (struct array, one entry per PM control field; may be 0x0):
%   B, converged, crit, crit_tol, medium_status
% expected_fields (numeric vector): every required ordered field (T), for the coverage-
%   completeness check (a required field missing ANY of expected_nH's rows fires (b) -- a lost
%   chunk during aggregation must be loud, never silently short a field).
% expected_nH (numeric vector, e.g. [33 65 129]).
%
% rep.fail_a .. rep.fail_e (logical scalars) and rep.pass = ~(fail_a|fail_b|fail_c|fail_d|fail_e)
% -- EXACTLY that formula, no extra term. rep.detail.* records which row(s)/field(s) tripped
% each clause, so a failure is always attributable, never merely asserted.
%
% Frozen constants (user-approved 2026-07-25; see invzp_convg_diagnosis.md Section 8.2):
%   omit_promote = 0.10.  crit_tol is READ per-row (it is threaded through, not re-derived, so
%   this function never needs to guess the caller's Ecut/T-independent constant); D_tol/Dq_tol
%   are likewise read per-row (state-dependent). This function hardcodes ONLY omit_promote,
%   the one clause-(c) constant it evaluates directly.
if nargin < 4
    error('invz:gate0AggregateArgs', 'invz_gate0_aggregate requires (ordered_rows, pm_rows, expected_fields, expected_nH).');
end
omit_promote = 0.10;

detail = struct('fail_a_rows', [], 'fail_b_rows', [], 'fail_c_rows', [], ...
                 'fail_d_fields', [], 'fail_e_rows', [], 'fail_e_pm', [], ...
                 'missing_rows', {{}});

n = numel(ordered_rows);
row_fail_a = false(1, n);  row_fail_b = false(1, n);  row_fail_c = false(1, n);
row_fail_e = false(1, n);

for k = 1:n
    row = ordered_rows(k);
    row_fail_a(k) = local_fail_a(row);
    row_fail_b(k) = local_fail_b(row);
    row_fail_c(k) = local_fail_c(row, omit_promote);
    row_fail_e(k) = local_fail_e(row);
end
detail.fail_a_rows = find(row_fail_a);
detail.fail_b_rows = find(row_fail_b);
detail.fail_c_rows = find(row_fail_c);
detail.fail_e_rows = find(row_fail_e);

% ---- coverage-completeness: every expected (field, nH) combination must have exactly one row.
missing = {};
if n == 0
    have = zeros(0, 2);
else
    have = [[ordered_rows.B]', [ordered_rows.nH]'];
end
for fB = expected_fields(:).'
    for fnH = expected_nH(:).'
        cnt = nnz(abs(have(:,1) - fB) < 1e-9 & have(:,2) == fnH);
        if cnt ~= 1
            missing{end+1} = sprintf('B=%g nH=%d: %d rows (expected 1)', fB, fnH, cnt); %#ok<AGROW>
        end
    end
end
detail.missing_rows = missing;
fail_b_coverage = ~isempty(missing);

% ---- PM controls: contribute to (e) only (prereg SS3(e) names PM controls there explicitly).
np = numel(pm_rows);
row_fail_e_pm = false(1, np);
for k = 1:np
    row_fail_e_pm(k) = local_fail_e_pm(pm_rows(k));
end
detail.fail_e_pm = find(row_fail_e_pm);

% ---- (d): G17 crossing continuity, decided PER REQUIRED FIELD from its nH=65/129 rows.
fail_d_fields = {};
for fB = expected_fields(:).'
    r65  = local_find_row(ordered_rows, fB, 65);
    r129 = local_find_row(ordered_rows, fB, 129);
    if isempty(r65) || isempty(r129)
        continue;   % already reported via fail_b_coverage; (d) needs both rows to compare
    end
    bad = r65.g17_anomaly || r129.g17_anomaly;
    for yfield = {'g17_r', 'g17_crit'}
        g65 = r65.(yfield{1});  g129 = r129.(yfield{1});
        if any(strcmp(g65.status, {'unresolved', 'jump_exceeded'})), bad = true; end
        if any(strcmp(g129.status, {'unresolved', 'jump_exceeded'})), bad = true; end
        if strcmp(g65.status, 'ok') && strcmp(g129.status, 'ok') && ...
                g129.max_jump > g65.max_jump + 1e-15
            bad = true;
        end
    end
    if bad, fail_d_fields{end+1} = fB; end %#ok<AGROW>
end
detail.fail_d_fields = fail_d_fields;

rep = struct();
rep.fail_a = any(row_fail_a);
rep.fail_b = any(row_fail_b) || fail_b_coverage;
rep.fail_c = any(row_fail_c);
rep.fail_d = ~isempty(fail_d_fields);
rep.fail_e = any(row_fail_e) || any(row_fail_e_pm);
rep.pass   = ~(rep.fail_a || rep.fail_b || rep.fail_c || rep.fail_d || rep.fail_e);
rep.detail = detail;
end

% ---------------------------------------------------------------------------------------------
function row = local_find_row(rows, B, nH)
row = [];
for k = 1:numel(rows)
    if abs(rows(k).B - B) < 1e-9 && rows(k).nH == nH
        row = rows(k);  return;
    end
end
end

% ---------------------------------------------------------------------------------------------
function f = local_fail_a(row)
%LOCAL_FAIL_A prereg SS3(a): any required solved-path node has a non-'ok' REFERENCE
% denominator status. 'not_applicable' (a node that never reached reference evaluation, e.g.
% degenerate_doublet) is NOT a reference-status opinion and is excluded here, matching
% invz_hmf_status.m's own precedent for what counts as a medium domain failure.
ms = row.medium_status;
f = any(~cellfun(@(s) any(strcmp(s, {'ok', 'not_applicable'})), ms));
end

% ---------------------------------------------------------------------------------------------
function f = local_fail_b(row)
%LOCAL_FAIL_B prereg SS3(b): coverage identity. Every trc.nodes entry must classify into
% exactly one of the four buckets (n_accounted == n_nodes), the predictor phase must be
% present, and -- when a finite root was actually returned -- a matching 'root' phase entry
% with ok_final=true must exist in the ledger (a missing final root ledger entry).
n_accounted = nnz(~strcmp(row.bucket, 'unrecognized'));
coverage_ok = (n_accounted == row.n_nodes) && (row.n_nodes == numel(row.bucket));
predictor_ok = row.n_nodes >= 1 && strcmp(row.phase{1}, 'predictor');
if isfinite(row.hstar)
    ridx = find(strcmp(row.phase, 'root'));
    root_ok = ~isempty(ridx) && row.ok_final(ridx(end)) && ...
              abs(row.h(ridx(end)) - row.hstar) <= 1e-9 * max(1, abs(row.hstar));
else
    root_ok = true;   % no root claimed: nothing to cross-check
end
f = ~(coverage_ok && predictor_ok && root_ok);
end

% ---------------------------------------------------------------------------------------------
function f = local_fail_c(row, omit_promote)
%LOCAL_FAIL_C prereg SS3(c): max(omit_max) over the actual solved ledger's 'ok'-labelled nodes
% must not exceed omit_promote, AND a missing/non-finite ratio at an 'ok' node fails outright
% (never dropped by isfinite filtering / MATLAB's NaN-omitting max()).
ok_mask = strcmp(row.medium_status, 'ok');
ok_omit = row.omit_max(ok_mask);
any_nonfinite_ok = any(~isfinite(ok_omit));
if isempty(ok_omit)
    exceeds = false;
else
    finite_omit = ok_omit(isfinite(ok_omit));
    exceeds = ~isempty(finite_omit) && max(finite_omit) > omit_promote;
end
f = any_nonfinite_ok || exceeds;
end

% ---------------------------------------------------------------------------------------------
function f = local_fail_e(row)
%LOCAL_FAIL_E prereg SS3(e), ordered-field half: status='ok', a finite nonzero root, and a
% stable endpoint under the frozen crit/D_uni/Dq margins (state-dependent D_tol/Dq_tol, read
% from the row -- see the file header).
f = ~(strcmp(row.status, 'ok') && isfinite(row.hstar) && row.hstar > 0 && ...
      isfinite(row.crit_star) && row.crit_star > row.crit_tol && ...
      isfinite(row.D_uni_star) && row.D_uni_star > row.D_tol_star && ...
      isfinite(row.Dq_min_star) && row.Dq_min_star > row.Dq_tol_star);
end

% ---------------------------------------------------------------------------------------------
function f = local_fail_e_pm(row)
%LOCAL_FAIL_E_PM prereg SS3(e), PM-control half: "either PM control fails to return a
% converged finite positive-mass PM state."
f = ~(row.converged && isfinite(row.crit) && row.crit > row.crit_tol && ...
      strcmp(row.medium_status, 'ok'));
end
