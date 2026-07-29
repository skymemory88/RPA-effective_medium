function H = invzf_cluster_coeff_harness(fixture, opts)
%INVZF_CLUSTER_COEFF_HARNESS Exact-cluster coefficient gate for a candidate functional.
%
%   H = INVZF_CLUSTER_COEFF_HARNESS(FIXTURE, OPTS) implements the four steps of
%   invzp_convg_fix.md Sec. 6 for the two-site X-X bond:
%     1. compute each coupling-order coefficient by exact finite-cluster differentiation
%        (invzf_cluster_exact + invzf_taylor_extract);
%     2. compute it independently, from local data only (invzf_pair_local_manifest);
%     3. compare BEFORE any resummation or numerical fitting;
%     4. return the exact coefficients in a form that can be frozen as a regression
%        fixture.
%
%   A candidate functional is graded by passing OPTS.candidate, a handle J -> F_pair(J)
%   built from that candidate's own equations. The SAME extractor is applied to it, so
%   the comparison is between coefficients, never between converged end results --
%   "agreement of a final critical field cannot compensate for a wrong low-order
%   coefficient" (fix doc Sec. 6).
%
%   FIXTURE fields: Delta, M, h, beta  (a two-level local model; the exact cluster
%   oracle invzf_cluster_exact is two-level by construction).
%   OPTS: .K (4) highest coupling order graded; .extract (struct forwarded to
%   invzf_taylor_extract); .manifest (struct forwarded to invzf_pair_local_manifest);
%   .candidate ([]) handle J -> F; .candidate_name (''); .rel_tol (1e-8) and
%   .abs_floor (1e-12) for the verdicts.
%
%   VERDICT RULE. An order is graded only when BOTH sides have it: the manifest must
%   have derived it, and the extraction error bar must be smaller than the tolerance
%   band. Orders that fail either condition are returned as 'ungraded' with the
%   reason, never as a pass. A mismatch on a graded order is a REJECTION of the
%   candidate at that order, not an invitation to fit.
if nargin < 2, opts = struct(); end
K        = getf_local(opts, 'K', 4);
exopts   = getf_local(opts, 'extract', struct());
manopts  = getf_local(opts, 'manifest', struct());
cand     = getf_local(opts, 'candidate', []);
candname = getf_local(opts, 'candidate_name', 'candidate');
rel_tol  = getf_local(opts, 'rel_tol', 1e-8);
abs_floor= getf_local(opts, 'abs_floor', 1e-12);

D = fixture.Delta;  M = fixture.M;  h = fixture.h;  bt = fixture.beta;

% --- 1. exact cluster coefficients ---------------------------------------------------------
Jmat = @(j) [0 j; j 0];
fexact = @(j) getfield(invzf_cluster_exact(D, M, h, Jmat(j), bt, 0, 1), 'F'); %#ok<GFLD>
ex = invzf_taylor_extract(fexact, K, exopts);

% --- 2. independent local manifest ---------------------------------------------------------
man = invzf_pair_local_manifest(D, M, h, bt, manopts);

% --- 3. compare, order by order ------------------------------------------------------------
% A coefficient that the manifest predicts to be ZERO cannot be graded against its own
% magnitude, so tolerances are set from the PROBLEM scale -- the largest coefficient in
% the expansion -- and only floored, never scaled, by the individual coefficient. Without
% this, a correctly-vanishing coefficient is reported as ungradeable purely because its
% relative band collapses to zero.
a_scale = max(abs(ex.a(isfinite(ex.a))));
if isempty(a_scale) || a_scale == 0, a_scale = 1; end
grade_frac = getf_local(opts, 'grade_frac', 1e-6);   % extraction error above this fraction
rel_floor  = getf_local(opts, 'rel_floor', 1e-12);   % of a_scale means the order is ungradeable

rows = cell(K+1, 1);
for k = 0:K
    i = k + 1;
    ae = ex.a(i);  ee = ex.err(i);
    scale = max([abs(ae), abs_floor, rel_floor*a_scale]);
    band  = max(rel_tol*scale, ee);
    if i <= numel(man.derived) && man.derived(i)
        am = man.a(i);
        d  = abs(am - ae);
        if ee > grade_frac*a_scale
            verdict = 'ungraded';
            reason  = sprintf(['extraction error %.3g exceeds %.3g = %.0e of the ' ...
                'expansion scale %.3g'], ee, grade_frac*a_scale, grade_frac, a_scale);
        elseif d <= band
            verdict = 'match';   reason = '';
        else
            verdict = 'MISMATCH'; reason = sprintf('|manifest - exact| = %.6g > band %.3g', d, band);
        end
    else
        am = NaN;  d = NaN;  verdict = 'ungraded';
        reason = 'manifest does not derive this order (see invzf_pair_local_manifest)';
    end
    rows{i} = struct('order', k, 'exact', ae, 'exact_err', ee, ...
        'manifest', am, 'abs_diff', d, 'band', band, ...
        'verdict', verdict, 'reason', reason);
end
cmp = [rows{:}];

% --- candidate grading (same extractor) ----------------------------------------------------
candrows = [];
if ~isempty(cand)
    cx = invzf_taylor_extract(cand, K, exopts);
    crow = cell(K+1,1);
    for k = 0:K
        i = k+1;
        ae = ex.a(i);  ee = ex.err(i);  ac = cx.a(i);  ec = cx.err(i);
        scale = max(abs(ae), abs_floor);
        band = max([rel_tol*scale, ee, ec]);
        d = abs(ac - ae);
        if d <= band, v = 'match'; else, v = 'MISMATCH'; end
        crow{i} = struct('order', k, 'exact', ae, 'candidate', ac, ...
            'abs_diff', d, 'band', band, 'verdict', v, ...
            'candidate_err', ec, 'exact_err', ee);
    end
    candrows = [crow{:}];
end

anyMismatch = any(strcmp({cmp.verdict}, 'MISMATCH'));
nGraded = nnz(strcmp({cmp.verdict}, 'match'));

H = struct('fixture', fixture, 'K', K, 'exact', ex, 'manifest', man, ...
    'comparison', cmp, 'candidate_name', candname, 'candidate', candrows, ...
    'manifest_consistent', ~anyMismatch, 'n_orders_graded', nGraded, ...
    'rel_tol', rel_tol, 'abs_floor', abs_floor, ...
    'regression_fixture', struct('schema', 'invzf_pair_coeff_fixture/v1', ...
        'Delta', D, 'M', M, 'h', h, 'beta', bt, ...
        'a_exact', ex.a, 'a_exact_err', ex.err, ...
        'a_manifest', man.a, 'manifest_derived', man.derived));
end

function v = getf_local(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
