function M = invz_ewald_metrics()
%INVZ_EWALD_METRICS Frozen Gate-A comparison metrics (Step-4 Task 3).
%   M = INVZ_EWALD_METRICS() returns function-handle fields
%   M.mt / M.mj / M.mid / M.mhfd / M.mfd implementing the tolerance
%   definitions frozen in docs/invzp_ewald_prereg.md sec 3. Plain
%   fixture/helper function (NOT a test; runtests on tests/ does not collect
%   it -- the invz_odd_anchors.m precedent).
%
% Every metric returns a struct with AT LEAST:
%   pass            logical scalar; true iff worst_margin <= 0 (M_HFD also
%                   requires zero sign mismatches on gated components).
%   worst_margin    worst (largest) (normalized-error - threshold); the
%                   comparison PASSES everywhere iff worst_margin <= 0.
%   worst_ratio     worst normalized-error reading (dimensionless "how many
%                   multiples of the natural scale" diagnostic used for
%                   logging; for M_HFD this is exactly max|A-B|/H_scale, the
%                   quantity the prereg calibration reports and compares to
%                   the 1e-7 threshold).
%   max_abs_error   max_component |A-B|, RAW (unnormalized) units.
% plus metric-specific extras (M_HFD adds sign_ok/H_scale/n_gated; every
% metric tags its own name in .metric so a result can never be silently
% mislabeled, e.g. an M_FD result reported as M_T).
%
% Every metric REJECTS NaN/Inf in its numeric inputs with a hard error (not
% a soft fail), and requires the two compared operands to match in size.
% Comparisons are symmetric in the two operands, per design sec 5.2 / prereg
% sec 3. Metrics operate elementwise and are shape-agnostic (3x3 blocks,
% full [3,3,ntau,ntau(,nq)] tensors, or coupling arrays/scalars all work).
M = struct();
M.mt   = @metric_mt;
M.mj   = @metric_mj;
M.mid  = @metric_mid;
M.mhfd = @metric_mhfd;
M.mfd  = @metric_mfd;
end

% =====================================================================
% shared guards
% =====================================================================
function check_finite(X, label)
if ~all(isfinite(X(:)))
    error('invz:ewaldMetricNonFinite', 'invz_ewald_metrics: %s contains NaN/Inf.', label);
end
end

function check_same_size(A, B)
if ~isequal(size(A), size(B))
    error('invz:ewaldMetricShape', ...
        'invz_ewald_metrics: operands must match in size (got [%s] vs [%s]).', ...
        num2str(size(A)), num2str(size(B)));
end
end

function r = safe_ratio(numv, denv)
% Elementwise numv./denv with a realmin floor on denv to avoid 0/0 -> NaN.
% The frozen definitions guarantee denv==0 implies numv==0 too (a zero
% "allowed" only arises from an all-zero comparison), so flooring denv is
% exact bookkeeping, not an approximation.
r = numv ./ max(denv, realmin);
end

function out = empty_pass_result(tag)
% Vacuous PASS for a zero-element comparison (nothing to disagree on).
out.pass = true; out.worst_margin = -Inf; out.worst_ratio = 0; out.max_abs_error = 0;
out.metric = tag;
end

% =====================================================================
% M_T -- tensor agreement (design sec 5.2 / prereg sec 3)
%   |A-B| <= AbsTol_T + RelTol_T*max(|A|,|B|),  AbsTol_T = 1e-8*T_scale,
%   RelTol_T = 1e-8,  T_scale = max over the COMPLETE tensor of max(|A|,|B|).
% =====================================================================
function out = metric_mt(A, B)
check_finite(A, 'A'); check_finite(B, 'B'); check_same_size(A, B);
if isempty(A), out = empty_pass_result('M_T'); return; end
Tscale = max(max(abs(A(:))), max(abs(B(:))));
AbsTol_T = 1e-8*Tscale;
RelTol_T = 1e-8;
diffv = abs(A - B);
allowed = AbsTol_T + RelTol_T*max(abs(A), abs(B));
margin = diffv - allowed;
out.pass          = max(margin(:)) <= 0;
out.worst_margin  = max(margin(:));
out.worst_ratio   = max(safe_ratio(diffv(:), allowed(:)));
out.max_abs_error = max(diffv(:));
out.T_scale       = Tscale;
out.metric        = 'M_T';
end

% =====================================================================
% M_J -- coupling agreement (Jnu/Jcc0/Jaa0/Juni; fixed J_ref scale)
%   |Ja-Jb| <= AbsTol_J + RelTol_J*max(|Ja|,|Jb|), AbsTol_J = 1e-8*J_ref,
%   RelTol_J = 1e-8, J_ref = 0.006424435656 (frozen prereg sec 3 literal).
% =====================================================================
function out = metric_mj(a, b)
check_finite(a, 'a'); check_finite(b, 'b'); check_same_size(a, b);
if isempty(a), out = empty_pass_result('M_J'); return; end
J_ref = 0.006424435656;             % frozen literal, docs/invzp_ewald_prereg.md sec 3
AbsTol_J = 1e-8*J_ref;
RelTol_J = 1e-8;
diffv = abs(a - b);
allowed = AbsTol_J + RelTol_J*max(abs(a), abs(b));
margin = diffv - allowed;
out.pass          = max(margin(:)) <= 0;
out.worst_margin  = max(margin(:));
out.worst_ratio   = max(safe_ratio(diffv(:), allowed(:)));
out.max_abs_error = max(diffv(:));
out.J_ref         = J_ref;
out.metric        = 'M_J';
end

% =====================================================================
% M_id -- exact machine-level structural identity (prereg sec 3)
%   max_components|A-B| <= 1e-12*T_scale (ONE global scale, no per-component
%   relative term).
% =====================================================================
function out = metric_mid(A, B)
check_finite(A, 'A'); check_finite(B, 'B'); check_same_size(A, B);
if isempty(A), out = empty_pass_result('M_id'); return; end
Tscale = max(max(abs(A(:))), max(abs(B(:))));
allowed = 1e-12*Tscale;
diffv = abs(A - B);
me = max(diffv(:));
out.pass          = (me - allowed) <= 0;
out.worst_margin  = me - allowed;
out.worst_ratio   = me / max(allowed, realmin);
out.max_abs_error = me;
out.T_scale       = Tscale;
out.metric        = 'M_id';
end

% =====================================================================
% M_HFD -- screened-Hessian finite-difference metric (prereg sec 3)
%   H_scale(x;A,B) = max(|x|^-3, max_ab|A_ab|, max_ab|B_ab|).
%   PASS iff max_ab|A_ab-B_ab|/H_scale <= 1e-7 AND every component with
%   max(|A_ab|,|B_ab|) >= 1e-8*H_scale has sign(A_ab) == sign(B_ab).
% =====================================================================
function out = metric_mhfd(A, B, x)
check_finite(A, 'A'); check_finite(B, 'B'); check_finite(x, 'x'); check_same_size(A, B);
if numel(x) ~= 3
    error('invz:ewaldMetricShape', 'invz_ewald_metrics: x must have exactly 3 components (got %d).', numel(x));
end
if isempty(A), out = empty_pass_result('M_HFD'); return; end
rx = norm(x(:));
Hscale = max([rx^-3, max(abs(A(:))), max(abs(B(:)))]);
diffv = abs(A - B);
me = max(diffv(:));
normerr = me / Hscale;
gate = max(abs(A), abs(B)) >= 1e-8*Hscale;
signA = sign(A); signB = sign(B);
sign_ok = all(signA(gate) == signB(gate));
out.pass          = (normerr <= 1e-7) && sign_ok;
out.worst_margin  = normerr - 1e-7;
out.worst_ratio   = normerr;
out.max_abs_error = me;
out.H_scale       = Hscale;
out.sign_ok       = sign_ok;
out.n_gated       = nnz(gate);
out.metric        = 'M_HFD';
end

% =====================================================================
% M_FD -- scalar-oracle finite-difference metric (fixed dimensioned floor;
%   deliberately looser than M_T to cover calibrated FD error). NEVER label
%   an M_FD pass as M_T.
%   |A-B| <= AbsTol_FD + RelTol_FD*max(|A|,|B|), AbsTol_FD = 2e-8 Angstrom^-3,
%   RelTol_FD = 1e-7.
% =====================================================================
function out = metric_mfd(A, B)
check_finite(A, 'A'); check_finite(B, 'B'); check_same_size(A, B);
if isempty(A), out = empty_pass_result('M_FD'); return; end
AbsTol_FD = 2e-8;
RelTol_FD = 1e-7;
diffv = abs(A - B);
allowed = AbsTol_FD + RelTol_FD*max(abs(A), abs(B));
margin = diffv - allowed;
out.pass          = max(margin(:)) <= 0;
out.worst_margin  = max(margin(:));
out.worst_ratio   = max(safe_ratio(diffv(:), allowed(:)));
out.max_abs_error = max(diffv(:));
out.metric        = 'M_FD';
end
