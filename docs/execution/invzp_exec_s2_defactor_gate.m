function G = invzp_exec_s2_defactor_gate(save_path)
%INVZP_EXEC_S2_DEFACTOR_GATE Gate for the opt-in defactored static closure coordinate.
%
% Execution packet S2 of docs/execution/invzp_plan_execution_diary.md.
% The fixture is a REAL node (T = 0.1 K, Bx = 3.825 T, h = 0.00241 -- the node that
% binds the column under mix_outer = 0.40) built through the same node constructor
% the solver uses, so every weight, lambda and coupling is the production one.
%
%   G1 EQUIVALENCE. Over the K0 range the failing node's iteration actually visits,
%      the factored and reciprocal q-averages must agree numerically. The bar is
%      rel_err <= 1e-11, NOT a few ulps: mean_q Gq and mean_q(J*Gq) are signed sums
%      over 16384 modes whose summands change sign (D_q crosses zero inside the
%      grid), and d0 = 1 + Sigma0 + K0*G0inel0 is itself formed by cancellation, so
%      both arrangements inherit an amplified rounding error that is common-mode and
%      not attributable to the reassociation. Measured behaviour on this fixture:
%      typical rel_err ~ 1e-15, worst ~8e-13, and -- decisively for the claim being
%      made -- the agreement is BEST, not worst, as d0 -> 0 (rel_err ~ 3e-15 at
%      d0 = 0.028). A pole-tracking error growth would falsify the reassociation;
%      its absence is what this gate checks.
%
%   G2 FINITENESS AT THE POLE. At the K0 that zeroes d0 = 1 + Sigma0 + K0*G0inel0 the
%      factored form must be non-finite while the reciprocal form stays finite and
%      reproduces the analytic limit Gq -> 1/(J(q) - K0).
%
%   G3 DEFAULT-PATH BIT-IDENTITY plus SAME FIXED POINT. Calling
%      invz_emt_static_ordered with no closure_coordinate option must return exactly
%      what the explicit 'factored' option returns (the historical arithmetic is
%      untouched), and the 'defactored' option must converge to the same K0.
if nargin < 1, save_path = ''; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));
addpath(fullfile(root, 'docs', 'diagnostics', 'invzp_solver_stability_2026-07-27'));

T = 0.1;  Bx = 3.825;  h = 0.0024100;
ion = invz_ion();
bz = struct('grid', [16 16 16], 'dpRng', 30, 'cache', true, 'dipole', 'bruteforce');
[Jf, ci, Jaa0] = invz_bz_couplings(ion, bz);
Jf = Jf(:);
ctx = invz_ordered_node_context(ion, T, [Bx 0 0], Jf, ...
    struct('J0eff', ci.Jcc0, 'Jxx0', Jaa0, 'hyp', true, 'transverse_mf', 'legacy_x'));
[node, nmeta] = invz_ordered_make_node(ctx, h);
if isempty(node)
    error('node construction failed at h = %g: %s', h, nmeta.status);
end
tl = node.tl;  beta = node.beta;  G0inel0 = node.G0inel0;  G0el0 = node.G0el0;
J0eff = node.J0eff;
% lam / Sigma0 frozen at the values the converged neighbour node carries; the static
% closure treats them as fixed inputs (invz_emt_static_ordered header).
Sigma0 = 0.0;
K = zeros(numel(node.wn), 1);
lam3 = invz_lambdas(K, node.g, node.wts, beta, [1 2 3]);
lam = lam3(1:2);

d0 = @(K0) 1 + Sigma0 + K0*G0inel0;
K0_pole = -(1 + Sigma0)/G0inel0;

% --- G1: equivalence over the visited K0 range ---------------------------------------------
K0_scan = linspace(-0.009, 0.019, 57).';    % the tail range of the failing node's iteration
n = numel(K0_scan);
relG = nan(n,1); relJ = nan(n,1); kap = nan(n,1); d0v = nan(n,1); bothfin = false(n,1);
for k = 1:n
    [gb_f, jl_f, ~, dk, Gq] = local_factored(tl, lam, K0_scan(k), Sigma0, beta, G0inel0, G0el0, Jf);
    rc = local_defactored(tl, lam, K0_scan(k), Sigma0, beta, G0inel0, G0el0, Jf);
    d0v(k) = dk;
    kap(k) = mean(abs(Gq))/max(abs(gb_f), realmin);
    bothfin(k) = isfinite(gb_f) && isfinite(jl_f) && isfinite(rc.Gbar) && isfinite(rc.Jloc);
    if bothfin(k)
        relG(k) = abs(rc.Gbar - gb_f)/max(abs(gb_f), realmin);
        relJ(k) = abs(rc.Jloc - jl_f)/max(abs(jl_f), realmin);
    end
end
sel = bothfin;
G1_TOL = 1e-11;
G1_worstG = max(relG(sel));
G1_worstJ = max(relJ(sel));
% Near-pole subset: the reassociation's whole claim is that accuracy does NOT degrade
% as d0 -> 0. Compare the decade closest to the pole against the rest.
near = sel & (abs(d0v) < 0.15);
far  = sel & (abs(d0v) >= 0.15);
G1_near_worst = max([relG(near); relJ(near)]);
G1_far_worst  = max([relG(far);  relJ(far)]);
G1_pass = isfinite(G1_worstG) && isfinite(G1_worstJ) && ...
          G1_worstG <= G1_TOL && G1_worstJ <= G1_TOL && ...
          G1_near_worst <= G1_far_worst;

% --- G2: at and adjacent to the pole -------------------------------------------------------
eps_list = [0, 1e-16, 1e-13, 1e-10, 1e-7];
fac_finite = false(size(eps_list));  def_finite = false(size(eps_list));
def_Jloc = nan(size(eps_list));      def_Gbar = nan(size(eps_list));
for k = 1:numel(eps_list)
    K0k = K0_pole + eps_list(k);
    [gb_f, jl_f] = local_factored(tl, lam, K0k, Sigma0, beta, G0inel0, G0el0, Jf);
    rc = local_defactored(tl, lam, K0k, Sigma0, beta, G0inel0, G0el0, Jf);
    fac_finite(k) = isfinite(gb_f) && isfinite(jl_f);
    def_finite(k) = isfinite(rc.Gbar) && isfinite(rc.Jloc);
    def_Jloc(k) = rc.Jloc;  def_Gbar(k) = rc.Gbar;
end
w_lim = 1./(Jf - K0_pole);
Jloc_lim = mean(Jf.*w_lim)/mean(w_lim);
Gbar_lim = mean(w_lim);
rel_lim_J = abs(def_Jloc(1) - Jloc_lim)/abs(Jloc_lim);
rel_lim_G = abs(def_Gbar(1) - Gbar_lim)/abs(Gbar_lim);
G2_pass = ~fac_finite(1) && all(def_finite) && rel_lim_J < 1e-10 && rel_lim_G < 1e-10;

% --- G3: default bit-identity + same fixed point -------------------------------------------
K0_seed = 0.0015;
eso_def = struct('static_medium', 'resummed', 'warn', false);
eso_fac = eso_def;  eso_fac.closure_coordinate = 'factored';
eso_dfc = eso_def;  eso_dfc.closure_coordinate = 'defactored';
[Ka, Ga, oa] = invz_emt_static_ordered(tl, lam, Sigma0, Jf, K0_seed, beta, ...
    J0eff, G0inel0, G0el0, eso_def);
[Kb, Gb, ob] = invz_emt_static_ordered(tl, lam, Sigma0, Jf, K0_seed, beta, ...
    J0eff, G0inel0, G0el0, eso_fac);
[Kc, Gc, oc] = invz_emt_static_ordered(tl, lam, Sigma0, Jf, K0_seed, beta, ...
    J0eff, G0inel0, G0el0, eso_dfc);
bit_same = isequal(Ka, Kb) && isequal(Ga, Gb) && isequal(oa.resid, ob.resid) && ...
           isequal(oa.iters, ob.iters);
K_rel = abs(Kc - Ka)/max(abs(Ka), realmin);
both_conv = oa.converged && oc.converged;
G3_pass = bit_same && both_conv && K_rel < 1e-8;

G = struct('fixture', struct('T', T, 'Bx', Bx, 'h', h, 'G0inel0', G0inel0, ...
        'G0el0', G0el0, 'Sigma0', Sigma0, 'lam', lam, 'K0_pole', K0_pole, 'nq', numel(Jf)), ...
    'G1_worst_rel_Gbar', G1_worstG, 'G1_worst_rel_Jloc', G1_worstJ, 'G1_tol', G1_TOL, ...
    'G1_near_pole_worst', G1_near_worst, 'G1_far_from_pole_worst', G1_far_worst, ...
    'G1_kappa_range', [min(kap(sel)) max(kap(sel))], 'G1_d0_range', [min(d0v) max(d0v)], ...
    'G1_n_compared', nnz(sel), 'G1_pass', G1_pass, ...
    'G2_factored_finite_at_pole', fac_finite(1), 'G2_defactored_finite', def_finite, ...
    'G2_rel_err_vs_limit_Jloc', rel_lim_J, 'G2_rel_err_vs_limit_Gbar', rel_lim_G, ...
    'G2_pass', G2_pass, ...
    'G3_bit_identical_default', bit_same, 'G3_both_converged', both_conv, ...
    'G3_K0_rel_diff', K_rel, 'G3_pass', G3_pass, ...
    'K0_factored', Ka, 'K0_defactored', Kc, 'Gstat_factored', Ga, 'Gstat_defactored', Gc, ...
    'resid_factored', oa.resid, 'resid_defactored', oc.resid, ...
    'iters_factored', oa.iters, 'iters_defactored', oc.iters);

fprintf('=== S2 defactored-coordinate gate (real node T=%.2f K, Bx=%.3f T, h=%.6g) ===\n', T, Bx, h);
fprintf('fixture: G0inel0 = %.6g, G0el0 = %.6g, nq = %d, K0 pole at %.6g\n', ...
    G0inel0, G0el0, numel(Jf), K0_pole);
fprintf('G1 equivalence: %d points compared, d0 in [%.3g %.3g], kappa in [%.3g %.3g]\n', ...
    nnz(sel), min(d0v), max(d0v), min(kap(sel)), max(kap(sel)));
fprintf('   worst rel_err: Gbar %.3g, Jloc %.3g (tol %.g)\n', G1_worstG, G1_worstJ, G1_TOL);
fprintf('   worst rel_err near the pole (|d0|<0.15) %.3g  vs  away %.3g  -> %s\n', ...
    G1_near_worst, G1_far_worst, pf(G1_pass));
fprintf('G2 at the pole: factored finite = %d; defactored finite = %s\n', ...
    fac_finite(1), mat2str(def_finite));
fprintf('   defactored vs analytic limit: rel err Jloc = %.3g, Gbar = %.3g  -> %s\n', ...
    rel_lim_J, rel_lim_G, pf(G2_pass));
fprintf('G3 default == explicit ''factored'' bitwise: %d; both converged: %d; K0 rel diff %.3g -> %s\n', ...
    bit_same, both_conv, K_rel, pf(G3_pass));
fprintf('   factored   K0 = %.16g  Gstat = %.10g  resid = %.3g  iters = %d  conv = %d\n', ...
    Ka, Ga, oa.resid, oa.iters, oa.converged);
fprintf('   defactored K0 = %.16g  Gstat = %.10g  resid = %.3g  iters = %d  conv = %d\n', ...
    Kc, Gc, oc.resid, oc.iters, oc.converged);
fprintf('OVERALL: %s\n', pf(G1_pass && G2_pass && G3_pass));

if ~isempty(save_path), save(save_path, 'G'); end
end

function s = pf(t), if t, s = 'PASS'; else, s = 'FAIL'; end, end

function [Gbar, Jloc, Gs, d0, Gq] = local_factored(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, Jf)
[Gs, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0);
d0 = go.gstat_local_denom;
Gq = Gs ./ (1 + (Jf - K0).*Gs);
Gbar = mean(Gq);
Jloc = mean(Jf .* Gq) / Gbar;
end

function rc = local_defactored(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, Jf)
[~, go] = invz_gstat_ordered(tl, lam, K0, Sigma0, beta, G0inel0, G0el0, ...
    struct('stable_form', true));
rc = invz_reciprocal_static_closure(go, Jf, K0);
end
