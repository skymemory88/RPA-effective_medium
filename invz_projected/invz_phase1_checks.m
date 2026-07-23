function res = invz_phase1_checks(ion, g, c, J_ref)
%INVZ_PHASE1_CHECKS Items 1, 3, 4, and item 5's per-configuration summary (docs/invzp_phase1_
% quadrature_prereg.md "Frozen pass/fail rules"), evaluated on an ALREADY-BUILT grid (g,
% invz_phase1_qgrid.m) and its ALREADY-EVALUATED couplings (c, invz_phase1_couplings.m, SAME g)
% for ONE (convention,offset,N,dpRng,gammaPolicy) configuration.
%
% Item 2 (periodicity) is convention/dpRng-level, not offset-specific -- see
% invz_phase1_check_periodicity.m. Item 5's REFINEMENT gate and item 6's OFFSET-sensitivity gate
% compare MULTIPLE configs' res.item5 summaries (via invz_phase1_gate.m) -- that cross-config
% assembly is owned by the driver, invz_phase1_run.m, not by this per-config function.
%
% INPUTS
%   ion    invz_ion()-shaped struct (needs .tau, for item 3's P_drop post-filter Gamma re-check).
%   g      struct from invz_phase1_qgrid (this config's grid).
%   c      struct from invz_phase1_couplings (this config's couplings, built from the SAME g).
%   J_ref  frozen reference scale (0.006424435656 meV), normalizes item 5's shape summaries.
%
% OUTPUT res (struct): res.item1 / .item3 / .item4 / .item5, each a sub-struct carrying both the
% frozen numbers and a logical .pass (item5 carries no single .pass -- its gate is cross-config).
if ~(isnumeric(J_ref) && isscalar(J_ref) && isfinite(J_ref) && J_ref > 0)
    error('invz:phase1Config', 'invz_phase1_checks: J_ref must be a finite positive scalar; got %s.', mat2str(J_ref));
end

% --- item 1: point uniqueness (tol_uniq = 1e-12, wrapped coordinates) ----------------------
TOL_UNIQ = 1e-12;
qw = mod(g.qvec + 0.5, 1) - 0.5;    % g.qvec is already wrapped by invz_phase1_qgrid; re-wrapping
                                     % here is an idempotent no-op that keeps this check
                                     % self-contained if ever handed a foreign qvec.
qdistinct = uniquetol(qw, TOL_UNIQ, 'ByRows', true, 'DataScale', 1);   % DataScale=1: ABSOLUTE tol
res.item1.tol_uniq = TOL_UNIQ;
res.item1.nominal  = size(qw, 1);
res.item1.distinct = size(qdistinct, 1);
res.item1.n_dup    = res.item1.nominal - res.item1.distinct;
res.item1.pass     = (res.item1.n_dup == 0);

% --- item 3: cardinality + Gamma count -------------------------------------------------------
res.item3.nominal     = g.nominal;
res.item3.n_gamma     = g.n_gamma;
res.item3.rows        = size(g.qvec, 1);
res.item3.gammaPolicy = g.gammaPolicy;
res.item3.total_weight = sum(g.w);
switch g.gammaPolicy
    case 'P_complete'
        res.item3.expected_rows = g.nominal;
        res.item3.weight_rule   = sprintf('uniform 1/%d over all nominal rows (Gamma kept)', g.nominal);
        res.item3.pass = (res.item3.rows == res.item3.expected_rows);
    case 'P_drop'
        res.item3.expected_rows = g.nominal - g.n_gamma;
        res.item3.weight_rule   = sprintf('uniform 1/%d over the %d Gamma-filtered rows', ...
            res.item3.expected_rows, res.item3.expected_rows);
        n_gamma_remaining = nnz_gamma(ion, g.qvec);
        res.item3.n_gamma_after_drop = n_gamma_remaining;
        res.item3.pass = (res.item3.rows == res.item3.expected_rows) && (n_gamma_remaining == 0);
    otherwise
        error('invz:phase1Config', 'invz_phase1_checks: unknown gammaPolicy ''%s''.', g.gammaPolicy);
end

% --- item 4: weight normalization (|sum(w)-1| <= 1e-12) --------------------------------------
res.item4.total_weight = sum(g.w);
res.item4.abs_err      = abs(res.item4.total_weight - 1);
res.item4.tol          = 1e-12;
res.item4.pass         = res.item4.abs_err <= res.item4.tol;

% --- item 5: multiset summary (single-config; the refinement/offset GATE is cross-config, see
% invz_phase1_run.m) ---------------------------------------------------------------------------
qs = [0.05 0.25 0.5 0.75 0.95];
x  = c.Jnu_flat;
raw.mean = mean(x);  raw.var = var(x);  raw.min = min(x);  raw.max = max(x);  raw.q = local_quantile(x, qs);
res.item5.raw   = raw;
res.item5.J_ref = J_ref;
res.item5.norm.mean = raw.mean / J_ref;
res.item5.norm.var  = raw.var  / J_ref^2;
res.item5.norm.min  = raw.min  / J_ref;
res.item5.norm.max  = raw.max  / J_ref;
res.item5.norm.q    = raw.q    / J_ref;
res.item5.energy.J0eff  = c.J0eff;
res.item5.energy.Jcc0   = c.Jcc0;
res.item5.energy.maxJnu = c.maxJnu;
res.item5.quantile_levels = qs;
end

function n = nnz_gamma(ion, qvec)
n = 0;
for i = 1:size(qvec, 1)
    if invz_is_gamma_equiv(qvec(i,:), ion.tau), n = n + 1; end
end
end

function v = local_quantile(x, p)
%LOCAL_QUANTILE Linear-interpolation quantile (the "R-7"/numpy-default convention), toolbox-free
% (Statistics and Machine Learning Toolbox not assumed, matching invz_robust_pct.m's own
% toolbox-free stance -- this is a DIFFERENT, more standard interpolated estimator than that
% function's nearest-rank clip helper; the two are not interchangeable and serve different
% purposes). p may be a vector; v is returned the same size as p.
x = sort(x(:));
n = numel(x);
v = zeros(size(p));
for k = 1:numel(p)
    if n == 1
        v(k) = x(1);
        continue;
    end
    idx  = 1 + p(k) * (n - 1);      % 1-based fractional rank
    lo   = floor(idx);  hi = ceil(idx);
    frac = idx - lo;
    v(k) = x(lo) + frac * (x(hi) - x(lo));
end
end
