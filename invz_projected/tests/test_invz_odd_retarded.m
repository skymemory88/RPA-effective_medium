function tests = test_invz_odd_retarded
%TEST_INVZ_ODD_RETARDED T2.1 + T2.2 retarded ODD-mediated coupling.
% T2.1: invz_emt_scalar generalized to per-frequency mode spectra Jnu_flat
% [nJ,nw] (vector path byte-identical, constant columns bitwise-equal);
% invz_solve_point opts.odd_retarded builds deltaJ(q, i*wn) via the scalar
% surrogate chi_perp(i*wn) ~ r_n*chi_perp(0) + first-order eigenvalue
% perturbation, cross-checked against the full per-(q,n) eig behind
% opts.odd_retarded_exact. Physics numbers are REPORTED, never tuned.
% T2.2 (INVZ_SLOW=1): the static-vs-retarded decision measurement -- Tc(0.1 T)
% and Bc(1.2 K) both flags; the 10 mK rule picks the default (see ODD-LOG T2).
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
end

function test_emt_matrix_column_constant_bitwise(testCase)
% [nJ,1] path untouched; [nJ,nw] with identical columns reproduces it exactly.
G0 = -linspace(30, 5, 12).';  Sigma = 0.05*ones(12,1);
Jf = linspace(-2e-3, 6e-3, 24).';
m1 = invz_emt_scalar(G0, Sigma, Jf, struct('debug', true));
m2 = invz_emt_scalar(G0, Sigma, repmat(Jf, 1, 12), struct('debug', true));
verifyEqual(testCase, m2.K, m1.K);           % bitwise
verifyEqual(testCase, m2.G, m1.G);
end

function test_retarded_rn_unity_bitwise(testCase)
% r_n forced to 1 (test hook opts.odd_rn_override) == static ODD solve, bitwise.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
p1 = invz_solve_point(ion, 1.8, 0.5, [], o);
o2 = o;  o2.odd_retarded = true;  o2.odd_rn_override = 1;    % force r_n = 1
p2 = invz_solve_point(ion, 1.8, 0.5, [], o2);
verifyEqual(testCase, p2.Sigma, p1.Sigma);   % bitwise (same EMT inputs)
verifyEqual(testCase, p2.crit,  p1.crit);
end

function test_retarded_physics_and_surrogate_error(testCase)
% Retarded weakens the n~=0 mediated coupling: Sigma0 decreases vs static.
% Also logs the scalar-surrogate + first-order-perturbation error vs exact eig.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
ps = invz_solve_point(ion, 1.8, 0.5, [], o);
o.odd_retarded = true;
pr = invz_solve_point(ion, 1.8, 0.5, [], o);
verifyTrue(testCase, ps.converged && pr.converged);
verifyLessThanOrEqual(testCase, pr.Sigma0, ps.Sigma0 + 1e-12);
verifyTrue(testCase, all(pr.odd.r_n > 0 & pr.odd.r_n <= 1 + 1e-12));
verifyLessThan(testCase, pr.sumrule_rel, ps.sumrule_rel + 1e-6);   % T2.1 acceptance
o.odd_retarded_exact = true;
pe = invz_solve_point(ion, 1.8, 0.5, [], o);
fprintf('retarded: Sigma0 static %.5f -> ret %.5f (exact %.5f); surrogate err %.2e\n', ...
    ps.Sigma0, pr.Sigma0, pe.Sigma0, abs(pr.Sigma0 - pe.Sigma0));
verifyEqual(testCase, pr.Sigma0, pe.Sigma0, 'AbsTol', 1e-3);       % surrogate sanity
% Scalar-surrogate smallness diagnostics (chi_ab tiny; chi_bb ~ chi_aa) + r_n
% profile head -- REPORTED for ODD-LOG T2, not gated:
fprintf('r_n(1:4) = %.6f %.6f %.6f %.6f; chi_ab/chi_aa max %.2e; |chi_bb-chi_aa|/chi_aa max %.2e\n', ...
    pr.odd.r_n(1:4), pr.odd.ab_ratio_max, pr.odd.bb_aa_dev_max);
end

function test_t2_decision_measurement_slow(testCase)
% T2.2: static-vs-retarded decision measurement + the 10 mK rule.
% Tc leg AMENDED from the brief's Tc(0.1 T)-via-invz_critical_T route
% (measurement + rationale in ODD-LOG T2): with ODD on the ordered side of
% the boundary NEVER converges (probed at 0.1-3 T; 16^3/dpRng 30), so
% invz_critical_T's converged-sign-change bracket is structurally unreachable
% at any field -- invz_critical survives only through its para_edge fallback,
% which the T-cut finder lacks (pre-existing gap, reported, out of this task's
% file scope). Instead the Tc leg extrapolates the converged PM-side crit(T)
% to zero at the 2 T proxy (well-split doublet; transversal cut) with the SAME
% deterministic estimator for both flags, so their difference isolates the
% retardation shift (the estimator's own PM-side bias cancels). Note Tc(0)
% itself is retardation-INVARIANT by construction: the zero-field closed form
% lives entirely in the elastic n = 0 sector where r_1 = 1 (plan section 4).
% The Bc(1.2 K) leg is per the brief.
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
tS = tc_pm_extrap(ion, 2.0, o);
o.odd_retarded = true;
tR = tc_pm_extrap(ion, 2.0, o);
bS = invz_critical(ion, 1.2, [], rmfield(o, 'odd_retarded'));
bR = invz_critical(ion, 1.2, [], o);
fprintf('T2.2: Tc*(2T) static %.6f vs retarded %.6f (d = %+.4f mK); Bc(1.2K) %.4f vs %.4f T\n', ...
    tS, tR, 1e3*(tR - tS), bS, bR);
verifyTrue(testCase, all(isfinite([tS tR bS bR])));
end

function tc = tc_pm_extrap(ion, B, o)
% PM-side Tc estimator at fixed B (T2.2 decision leg): sample crit on the house
% grid step (1/30 K), keep converged points only (the classifier convention of
% invz_critical_T), and linearly extrapolate the two LOWEST converged PM points
% to crit = 0. The FIXED grid makes the estimator deterministic: the static and
% retarded runs vote on identical T points, so tc(static) - tc(retarded) is a
% pure crit-shift readout, free of grid/bisection noise.
Tg = 1.20:1/30:1.50;                             % spans the 2 T slowing band
c  = nan(size(Tg));  okv = false(size(Tg));
for i = 1:numel(Tg)
    [c(i), okv(i)] = invz_crit_at(ion, Tg(i), B, [], o);
end
Tv = Tg(okv);  cv = c(okv);
assert(numel(cv) >= 2 && all(cv > 0), ...
    'tc_pm_extrap: need >= 2 converged PM points in the scan window (got %d).', numel(cv));
tc = Tv(1) - cv(1)*(Tv(2) - Tv(1))/(cv(2) - cv(1));   % two lowest converged points
end
