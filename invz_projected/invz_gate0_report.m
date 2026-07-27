function rep = invz_gate0_report(ion, T, ordered_fields, pm_fields, Jnu_flat, opts)
%INVZ_GATE0_REPORT Gate-0 domain/omitted-order diagnostic driver (Task 18). The frozen
% predicate it measures is quoted verbatim in
% docs/invzp_strict_medium_gate0_report.md SS1. A DIAGNOSTIC DRIVER, NOT A TEST: it measures
% and reports the frozen promotion predicate on the real production coupling multiset. It
% never widens a tolerance, switches scheme, or filters a failed node out of the gate. G = -chi
% (meV^-1), ferromagnetic positive J.
%
% rep = invz_gate0_report(ion, T, ordered_fields, pm_fields, Jnu_flat, opts)
%
% ordered_fields, pm_fields (row vectors, Tesla, ALONG [1 0 0] per the production convention):
%   May each be a SUBSET of the full preregistered lists, or []. The 60-70 minute full sweep
%   does not fit one foreground MATLAB call (600 s cap), so the intended usage is to invoke
%   this function ONCE PER CHUNK (one field, or one (field,nH) pair) and merge the returned
%   `rep.ordered` / `rep.pm` rows across calls before a FINAL invz_gate0_aggregate call on the
%   union -- the chunk plan actually used is recorded in
%   docs/invzp_strict_medium_gate0_report.md. The digest
%   check below runs on EVERY invocation regardless: reloading Jnu_flat from a .mat does not
%   exempt a caller from it.
% Jnu_flat: the exact production coupling column (invz_bz_couplings, frozen tuple, prereg SS8).
%   Verified here against the frozen exact-byte digest BEFORE any solve; a mismatch aborts.
% opts.J0eff (REQUIRED, no default -- caller-owned, matching invz_hmf_ordered's own contract).
% opts.Jxx0 (default ion.Jxx0), opts.Ecut (default 40), opts.nH_list (default [33 65 129] --
%   pass a SINGLE value, e.g. 33, to chunk one nH at a time), opts.include_b0 (default false --
%   the exact-B=0 hard-domain CONTROL; run it in exactly one small dedicated call, never
%   repeated per chunk), opts.verbose (default true -- prints one line per solved field/nH/PM
%   point, satisfying the brief's "the run command must print" requirement directly).
% The static-medium scheme (strict_1z_dyson_ref) and every judgement-bearing constant
% (crit_tol=1e-6, omit_promote=0.10, pole_cont_tol=1e-3, ref_margin=1e-6) are FROZEN and
% hardcoded below, never exposed as a caller override -- widening any of them is exactly what
% prereg SS3 forbids on a Gate-0 failure.
%
% rep.ordered (one row per (B,nH) processed THIS call), rep.pm (one row per PM field processed
% THIS call), rep.b0 ([] unless opts.include_b0) -- see invz_gate0_aggregate.m's header for the
% exact row schema (this function is what BUILDS those rows from trc/prof). rep.fail_a .. e and
% rep.pass are invz_gate0_aggregate's verdict OVER ONLY THE ROWS THIS CALL PRODUCED -- printed
% as PARTIAL, never the final verdict unless this one call covered every required field/nH/PM
% point. rep.digest, rep.Jscale, rep.Jmom, rep.covered_fields/covered_nH/covered_pm are threaded
% through for provenance and for the caller's own merge/aggregate step.
%
% No broad catch: an invz:* identifier or any other exception from invz_static_domain_scan,
% invz_hmf_ordered, invz_solve_point, or invz_pole_continuity propagates to the caller
% unchanged and aborts this report, exactly like an unclassified/wiring error should.
if nargin < 6, opts = struct(); end

% ---- FROZEN constants (user approval 2026-07-25; quoted in
% docs/invzp_strict_medium_gate0_report.md SS1).
% Never derived from output, never widened on a failure. ------------------------------------
FROZEN_DIGEST  = 'ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17';
CRIT_TOL       = 1e-6;
OMIT_PROMOTE   = 0.10;    %#ok<NASGU> -- consumed inside invz_gate0_aggregate, recorded here too
POLE_CONT_TOL  = 1e-3;
REF_MARGIN     = 1e-6;
STATIC_MEDIUM  = 'strict_1z_dyson_ref';

% ---- digest check FIRST, before any solve, on EVERY invocation (prereg SS8; brief Step 2). -
if ~(isnumeric(Jnu_flat) && isreal(Jnu_flat) && ~isempty(Jnu_flat))
    error('invz:gate0Config', 'Jnu_flat must be a nonempty real numeric array.');
end
digest = invz_exact_numeric_digest(Jnu_flat(:));
if ~strcmp(digest, FROZEN_DIGEST)
    error('invz:gate0DigestMismatch', ...
        'coupling digest mismatch: got %s, frozen preregistered digest is %s. Aborting before any solve.', ...
        digest, FROZEN_DIGEST);
end
% ---- strict scalar-vector configuration (brief Step 2): the ordered leg never accepts a
% [nJ,nw] retardation matrix under a strict scheme (invz_hmf_ordered.m / invz:staticMedium).
if ~isvector(Jnu_flat)
    error('invz:gate0Config', ...
        'strict scalar-vector configuration required: Jnu_flat must be a vector; got size %s.', ...
        mat2str(size(Jnu_flat)));
end
Jnu_flat = Jnu_flat(:);

if ~isfield(opts, 'J0eff')
    error('invz:gate0Config', 'opts.J0eff is required (no default -- caller-owned).');
end
J0eff   = opts.J0eff;
Jxx0    = getf(opts, 'Jxx0', ion.Jxx0);
Ecut    = getf(opts, 'Ecut', 40);
nH_list = getf(opts, 'nH_list', [33 65 129]);
include_b0 = getf(opts, 'include_b0', false);
verbose = getf(opts, 'verbose', true);

Jmom   = invz_coupling_moments(Jnu_flat);
Jscale = max(abs(Jnu_flat));

ordered_fields = ordered_fields(:).';
pm_fields      = pm_fields(:).';
nH_list        = nH_list(:).';

obase = struct('J0eff', J0eff, 'Jxx0', Jxx0, 'hyp', true, ...
               'static_medium', STATIC_MEDIUM, 'ref_margin', REF_MARGIN, ...
               'Jmom', Jmom, 'Ecut', Ecut);

ordered_rows = repmat(local_blank_row(), 1, 0);
for iB = 1:numel(ordered_fields)
    B = ordered_fields(iB);
    for inH = 1:numel(nH_list)
        nH = nH_list(inH);
        row = local_run_one_ordered(ion, T, B, Jnu_flat, obase, nH, CRIT_TOL, POLE_CONT_TOL, ...
                                     Jscale, verbose);
        ordered_rows(end+1) = row; %#ok<AGROW>
    end
end

pm_rows = repmat(local_blank_pm_row(), 1, 0);
for iB = 1:numel(pm_fields)
    B = pm_fields(iB);
    row = local_run_one_pm(ion, T, B, Jnu_flat, obase, CRIT_TOL, verbose);
    pm_rows(end+1) = row; %#ok<AGROW>
end

b0 = [];
if include_b0
    b0 = local_run_b0(ion, T, Jnu_flat, obase, verbose);
end

agg = invz_gate0_aggregate(ordered_rows, pm_rows, ordered_fields, nH_list);

rep = struct();
rep.digest = digest;
rep.Jscale = Jscale;
rep.Jmom   = Jmom;
rep.ordered = ordered_rows;
rep.pm      = pm_rows;
rep.b0      = b0;
rep.fail_a = agg.fail_a;  rep.fail_b = agg.fail_b;  rep.fail_c = agg.fail_c;
rep.fail_d = agg.fail_d;  rep.fail_e = agg.fail_e;  rep.pass   = agg.pass;
rep.detail = agg.detail;
rep.covered_fields = ordered_fields;
rep.covered_nH     = nH_list;
rep.covered_pm     = pm_fields;

% rep.g5 / rep.g17: explicit per-required-field SUMMARY tables (interface contract, brief
% "Produces" list) -- built from the SAME ordered_rows ingredients above, over whatever
% (field,nH) subset THIS call covers. rep.g17's per-field `pass` is read directly from
% invz_gate0_aggregate's OWN fail_d decision (never re-derived), so the two can never disagree.
rep.g5  = local_build_g5(ordered_rows, ordered_fields);
rep.g17 = local_build_g17(ordered_rows, ordered_fields, agg.detail.fail_d_fields);

if verbose
    fprintf(['[gate0] THIS CALL ONLY -- fields=%s nH=%s pm=%s b0=%d => PARTIAL pass=%d ' ...
             '(a=%d b=%d c=%d d=%d e=%d)\n'], mat2str(ordered_fields), mat2str(nH_list), ...
        mat2str(pm_fields), include_b0, rep.pass, rep.fail_a, rep.fail_b, rep.fail_c, ...
        rep.fail_d, rep.fail_e);
end
end

% ---------------------------------------------------------------------------------------------
function row = local_run_one_ordered(ion, T, B, Jnu_flat, obase, nH, crit_tol, pole_cont_tol, ...
                                      Jscale, verbose)
%LOCAL_RUN_ONE_ORDERED One (B,nH) ordered-field pass: prospective scan (metadata only), the
% actual solved-path trace, ledger reduction into invz_gate0_aggregate's row schema, the root
% cross-check, and the G17 pole-continuity evaluation on the final adapted prof.hgrid.
t0 = tic;
o = obase;  o.nH = nH;
scan = invz_static_domain_scan(ion, T, [B 0 0], Jnu_flat, o);   % prospective only (brief 2.1)

o.trace = true;
[hstar, prof, trc] = invz_hmf_ordered(ion, T, [B 0 0], Jnu_flat, o);

n_nodes = numel(trc.nodes);
if n_nodes == 0
    phase = {};  term_reason = {};  medium_status = {};  omit_max = [];
    ref_margin = [];  ok_final = false(1, 0);  h = [];
else
    phase         = {trc.nodes.phase};
    term_reason   = {trc.nodes.term_reason};
    medium_status = {trc.nodes.medium_status};
    omit_max      = [trc.nodes.omit_max];
    ref_margin    = [trc.nodes.ref_margin];
    ok_final      = [trc.nodes.ok_final];
    h             = [trc.nodes.h];
end
bucket = cellfun(@local_classify_term_reason, term_reason, 'UniformOutput', false);

% Root cross-check (brief Step 2.7): trc.nodes' 'root' entry and prof.*_star are built from the
% SAME eval_node call (invz_hmf_ordered.m eval_node -> append_trace_node, same `rec`), so they
% are bit-identical by construction; any disagreement is a wiring defect, never a numeric one,
% and must abort rather than be silently absorbed into a row flag.
if isfinite(hstar)
    ridx = find(strcmp(phase, 'root'));
    if isempty(ridx)
        error('invz:gate0Wiring', ...
            'B=%g nH=%d: hstar is finite (%.6g) but no ''root'' phase entry exists in trc.nodes.', ...
            B, nH, hstar);
    end
    ridx = ridx(end);
    if trc.nodes(ridx).crit ~= prof.crit_star || trc.nodes(ridx).D_uni ~= prof.D_uni_star || ...
            trc.nodes(ridx).Dq_min ~= prof.Dq_min_star
        error('invz:gate0Wiring', ...
            ['B=%g nH=%d: root trace record (crit=%.17g D_uni=%.17g Dq_min=%.17g) disagrees ' ...
             'with prof.crit_star/D_uni_star/Dq_min_star (%.17g/%.17g/%.17g).'], B, nH, ...
            trc.nodes(ridx).crit, trc.nodes(ridx).D_uni, trc.nodes(ridx).Dq_min, ...
            prof.crit_star, prof.D_uni_star, prof.Dq_min_star);
    end
end

D_tol_star  = crit_tol * max(1, abs(prof.Gstat_star) * Jscale);   % prereg SS1, at the root state
Dq_tol_star = D_tol_star;                                          % same construction (prereg SS1)

% G17 (d): the FINAL ADAPTED prof.hgrid, filtered to invz_pole_continuity's aligned-finite
% precondition. g17_anomaly: an 'ok'-labelled PROFILE node (prof.medium_status, aligned with
% prof.hgrid/r/crit/gstat_local_denom -- NOT the broader trc.nodes ledger) whose r, crit, or
% gstat_local_denom is non-finite anyway -- "require finite integrands everywhere" (brief 2.6).
pms          = prof.medium_status;
ok_mask_prof = strcmp(pms, 'ok');
finite_r     = isfinite(prof.hgrid) & isfinite(prof.gstat_local_denom) & isfinite(prof.r);
finite_crit  = isfinite(prof.hgrid) & isfinite(prof.gstat_local_denom) & isfinite(prof.crit);
g17_anomaly  = any(ok_mask_prof & ~finite_r) || any(ok_mask_prof & ~finite_crit);
g17_r    = invz_pole_continuity(prof.hgrid(finite_r), prof.gstat_local_denom(finite_r), ...
                                 prof.r(finite_r), pole_cont_tol);
g17_crit = invz_pole_continuity(prof.hgrid(finite_crit), prof.gstat_local_denom(finite_crit), ...
                                 prof.crit(finite_crit), pole_cont_tol);

row = local_blank_row();
row.B = B;  row.nH = nH;
row.status = prof.status;  row.hstar = hstar;
row.crit_star = prof.crit_star;  row.D_uni_star = prof.D_uni_star;  row.Dq_min_star = prof.Dq_min_star;
row.crit_tol = crit_tol;  row.D_tol_star = D_tol_star;  row.Dq_tol_star = Dq_tol_star;
row.n_nodes = n_nodes;
row.phase = phase;  row.bucket = bucket;  row.medium_status = medium_status;
row.omit_max = omit_max;  row.ok_final = logical(ok_final);  row.h = h;
row.g17_r = g17_r;  row.g17_crit = g17_crit;  row.g17_anomaly = g17_anomaly;
row.int_Sigma0 = prof.int_Sigma0;  row.int_r_minus_1 = prof.int_r_minus_1;
row.min_ref_margin = local_safe_min(ref_margin);
ok_ledger = strcmp(medium_status, 'ok');
row.max_omit_ledger = local_safe_max(omit_max(ok_ledger));
row.scan_n_valid = scan.n_valid;  row.scan_n_degenerate = scan.n_degenerate;
row.scan_n_out_of_domain = scan.n_out_of_domain;  row.scan_n_skipped = scan.n_skipped;
row.scan_max_omit = local_safe_max(scan.omit_max);
row.elapsed_s = toc(t0);
if verbose
    fprintf(['  B=%-6g nH=%-4d status=%-22s hstar=%-13.6g crit*=%-11.4g D_uni*=%-11.4g ' ...
             'Dq_min*=%-11.4g n_nodes=%-3d max_omit_ledger=%-11.5g int_rm1=%-13.6g ' ...
             'g17r=%-13s g17c=%-13s (%.1fs)\n'], ...
        B, nH, row.status, hstar, row.crit_star, row.D_uni_star, row.Dq_min_star, n_nodes, ...
        row.max_omit_ledger, row.int_r_minus_1, g17_r.status, g17_crit.status, row.elapsed_s);
end
end

% ---------------------------------------------------------------------------------------------
function row = local_run_one_pm(ion, T, B, Jnu_flat, obase, crit_tol, verbose)
%LOCAL_RUN_ONE_PM PM control field: invz_solve_point ONLY -- never the Jensen HMF solver.
t0 = tic;
o = obase;
pt = invz_solve_point(ion, T, [B 0 0], Jnu_flat, o);
row = local_blank_pm_row();
row.B = B;  row.converged = pt.converged;  row.crit = pt.crit;  row.crit_tol = crit_tol;
row.medium_status = pt.medium_status;  row.elapsed_s = toc(t0);
if verbose
    fprintf('  PM B=%-6g converged=%-5d crit=%-13.6g medium_status=%-20s (%.1fs)\n', ...
        B, pt.converged, pt.crit, pt.medium_status, row.elapsed_s);
end
end

% ---------------------------------------------------------------------------------------------
function b0 = local_run_b0(ion, T, Jnu_flat, obase, verbose)
%LOCAL_RUN_B0 Exact B=0 hard-domain CONTROL (prereg SS3 note; brief trap #2): the return-mode
% two-level/domain path is expected to report 'degenerate_doublet' at the h=0 predictor (no
% field at all to lift the zero-field doublet degeneracy). Reported separately; NEVER passed
% into invz_gate0_aggregate, so it structurally cannot fire (a) or (e).
t0 = tic;
o = obase;  o.nH = 33;  o.trace = true;
[hstar, prof, trc] = invz_hmf_ordered(ion, T, [0 0 0], Jnu_flat, o);
n_nodes = numel(trc.nodes);
if n_nodes == 0
    term_reason = {};
else
    term_reason = {trc.nodes.term_reason};
end
bucket = cellfun(@local_classify_term_reason, term_reason, 'UniformOutput', false);
n_accounted = nnz(~strcmp(bucket, 'unrecognized'));
b0 = struct('status', prof.status, 'hstar', hstar, 'n_nodes', n_nodes, ...
    'n_accounted', n_accounted, 'expected_degenerate', strcmp(prof.status, 'degenerate_doublet'), ...
    'elapsed_s', toc(t0));
if verbose
    fprintf('  B=0 (hard-domain CONTROL): status=%-22s hstar=%-11.6g n_nodes=%d n_accounted=%d expected_degenerate=%d (%.1fs)\n', ...
        b0.status, hstar, n_nodes, n_accounted, b0.expected_degenerate, b0.elapsed_s);
end
end

% ---------------------------------------------------------------------------------------------
function g5 = local_build_g5(ordered_rows, fields)
%LOCAL_BUILD_G5 Explicit per-required-field G5 summary (prereg SS5): int_Sigma0 and
% int_r_minus_1 at each nH THIS call covered (NaN for an nH not covered), the 33->65 approach
% diagnostic (NOT gated), and the frozen 65->129 criterion |I_fine-I_prev| <= I_atol +
% 1e-3*max(|I_fine|,|I_prev|) with I_atol = 1e-10 -- reported per integral, never folded into
% rep.pass (G5 is its own gate, prereg SS5, distinct from the (a)-(e) Gate-0 clauses).
I_atol = 1e-10;
g5 = [];
for iB = 1:numel(fields)
    B = fields(iB);
    rB = ordered_rows(abs([ordered_rows.B] - B) < 1e-9);
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
end

% ---------------------------------------------------------------------------------------------
function g17 = local_build_g17(ordered_rows, fields, fail_d_fields)
%LOCAL_BUILD_G17 Explicit per-required-field G17 summary: the nH=65 and nH=129
% invz_pole_continuity results for both y=r and y=crit. `pass` is read directly from
% invz_gate0_aggregate's OWN fail_d decision for this field (never re-derived here), so the
% summary and the gate can never silently disagree.
if isempty(fail_d_fields), fail_d_set = []; else, fail_d_set = cell2mat(fail_d_fields); end
g17 = [];
for iB = 1:numel(fields)
    B = fields(iB);
    r65  = ordered_rows(abs([ordered_rows.B] - B) < 1e-9 & [ordered_rows.nH] == 65);
    r129 = ordered_rows(abs([ordered_rows.B] - B) < 1e-9 & [ordered_rows.nH] == 129);
    g = struct('B', B, 'g17_r_65', [], 'g17_r_129', [], 'g17_crit_65', [], 'g17_crit_129', [], ...
        'pass', ~any(abs(fail_d_set - B) < 1e-9));
    if ~isempty(r65),  g.g17_r_65 = r65.g17_r;    g.g17_crit_65 = r65.g17_crit;  end
    if ~isempty(r129), g.g17_r_129 = r129.g17_r;  g.g17_crit_129 = r129.g17_crit; end
    if isempty(g17), g17 = g; else, g17(end+1) = g; end %#ok<AGROW>
end
end

% ---------------------------------------------------------------------------------------------
function b = local_classify_term_reason(tr)
%LOCAL_CLASSIFY_TERM_REASON THE frozen term_reason -> {ok, medium_out_of_domain,
% degenerate_doublet, solver_failed} mapping (brief Step 2.3), enumerated from source
% (invz_hmf_ordered.m's local_term_reason and its eval_node; invz_ordered_node_solve.m's
% run_attempt): 'converged' (<- accepted) and 'bare_shortcut' (the force_bare caller shortcut,
% rec.accepted=true unconditionally; unreached in this driver since force_bare is never set,
% but mapped for completeness) both terminate ok. 'refresh_failed' (<- loop_converged_not_
% accepted) and 'max_iter' both terminate a genuine SOLVER failure -- the map ran but did not
% produce an accepted fixed point. 'medium_out_of_domain' and 'degenerate_doublet' pass through
% verbatim (both are labelled domain events, never a generic failure -- invz_hmf_status.m's own
% binding precedence). Anything else (including the blank-record placeholder 'not_evaluated',
% which should never survive to a completed trc.nodes entry) is 'unrecognized' -- NOT folded
% into any other bucket, so an empty/unknown status is always visible to the coverage identity
% (brief Step 2.3's n_accounted == numel(trc.nodes)) rather than silently passing through an
% `otherwise` branch.
switch tr
    case 'converged',                    b = 'ok';
    case 'bare_shortcut',                 b = 'ok';
    case 'medium_out_of_domain',          b = 'medium_out_of_domain';
    case 'degenerate_doublet',            b = 'degenerate_doublet';
    case {'refresh_failed', 'max_iter'},  b = 'solver_failed';
    otherwise,                            b = 'unrecognized';
end
end

% ---------------------------------------------------------------------------------------------
function v = local_safe_min(x)
%LOCAL_SAFE_MIN min over the finite entries of x; NaN when none exist (never errors on empty).
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = min(x); end
end

function v = local_safe_max(x)
%LOCAL_SAFE_MAX max over the finite entries of x; NaN when none exist (never errors on empty).
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end

% ---------------------------------------------------------------------------------------------
function r = local_blank_row()
%LOCAL_BLANK_ROW Fixed schema for one ordered (B,nH) row -- see invz_gate0_aggregate.m's header
% for the fields invz_gate0_aggregate itself reads; the remainder are report-only.
r = struct('B', NaN, 'nH', NaN, 'status', '', 'hstar', NaN, 'crit_star', NaN, ...
    'D_uni_star', NaN, 'Dq_min_star', NaN, 'crit_tol', NaN, 'D_tol_star', NaN, ...
    'Dq_tol_star', NaN, 'n_nodes', 0, 'phase', {{}}, 'bucket', {{}}, ...
    'medium_status', {{}}, 'omit_max', [], 'ok_final', false(1, 0), 'h', [], ...
    'g17_r', [], 'g17_crit', [], 'g17_anomaly', false, ...
    'int_Sigma0', NaN, 'int_r_minus_1', NaN, 'min_ref_margin', NaN, ...
    'max_omit_ledger', NaN, 'scan_n_valid', NaN, 'scan_n_degenerate', NaN, ...
    'scan_n_out_of_domain', NaN, 'scan_n_skipped', NaN, 'scan_max_omit', NaN, ...
    'elapsed_s', NaN);
end

function r = local_blank_pm_row()
r = struct('B', NaN, 'converged', false, 'crit', NaN, 'crit_tol', NaN, ...
    'medium_status', '', 'elapsed_s', NaN);
end
