function tests = test_invz_odd_tier2
%TEST_INVZ_ODD_TIER2 T3.3 + T3.4 Tier-2 outer self-consistency in invz_solve_point.
% T3.3: opts.odd_tier2 closes the Dollberg variable-moments loop -- converged
% Tier-1 solve -> internal-field covariance C (invz_odd_fieldvar, E3) ->
% Gauss-Hermite-dressed doublet + disorder-averaged propagator
% (invz_twolevel_avg) -> re-run the inner EMT<->Sigma loop -> re-converge C
% (mix 0.5, tol 1e-3 rel, cap 8). Physics gates are DIRECTIONS only (plan
% sign contract): Tier 2 ADDS suppression (crit increases at a fixed PM
% point), the dressed splitting shows level repulsion (Delta_eff >= Delta).
% IR safety (slow): C stays bounded as crit -> 0 (ODD blocks vanish at q = 0).
% Combined Tc measurement (slow): PM-side crit(T)-extrapolation estimator at
% the 0.5 T proxy for off / Tier-1 / Tier-1+2 (AMENDED route; shared fixture
% tests/invz_odd_tc_pm_extrap.m).
% T3.4 (characterization, ODD-LOG SS T3.4): the exact-B = 0 behaviour is
% REPORTED, not gated -- the Tier-1 scaffold itself (bare invz_twolevel)
% throws invz:degenerateDoublet at B = 0 before any Tier-2 machinery runs.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));   % repo root: MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_tier2_flag_off_byte_identical(testCase)
% Flag off (absent OR explicit false): output isequaln to the plain ODD solve
% -- every Tier-2 field sits behind the flag. Plus the two structural guards:
% odd_tier2 requires odd, and the tier2 x retarded combination is refused
% (invz:oddArgs, "not yet validated" -- invz_odd_fieldvar assembles the
% equal-time Scc from the STATIC mode spectrum; a retarded solve's
% per-frequency modes have no validated E3 counterpart).
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
p1 = invz_solve_point(ion, 1.6, 0.5, [], o);
o2 = o;  o2.odd_tier2 = false;
p2 = invz_solve_point(ion, 1.6, 0.5, [], o2);
verifyTrue(testCase, isequaln(p1, p2));
% structural guards (both error BEFORE any lattice work)
oBad = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd_tier2', true);
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [], oBad), 'invz:oddArgs');
oRet = o;  oRet.odd_tier2 = true;  oRet.odd_retarded = true;
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [], oRet), 'invz:oddArgs');
oRex = o;  oRex.odd_tier2 = true;  oRex.odd_retarded_exact = true;
verifyError(testCase, @() invz_solve_point(ion, 1.6, 0.5, [], oRex), 'invz:oddArgs');
end

function test_tier2_converges_and_suppresses(testCase)
% T3.3 convergence gate at a GUARANTEED-PM point (1.80 K, 0.05 T): the plan's
% (1.55 K, 0.05 T proxy) sits below the no-ODD Tc and is REPORTED separately
% below (whether it converges tells which side of the ODD-shifted boundary it
% lands on). Variable moments suppress ordering: crit increases vs Tier 1.
ion = invz_ion();
n = 6;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 10, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
p1 = invz_solve_point(ion, 1.80, 0.05, [], o);
o.odd_tier2 = true;
p2 = invz_solve_point(ion, 1.80, 0.05, [], o);
verifyTrue(testCase, p1.converged && p2.converged);
verifyLessThanOrEqual(testCase, p2.tier2_iters, 8);
verifyGreaterThan(testCase, p2.crit, p1.crit - 1e-12);     % suppression direction
verifyGreaterThanOrEqual(testCase, p2.tla.Delta, p1.tl.Delta - 1e-12);  % level repulsion
fprintf('T3.3: tier2 iters = %d, crit %.5f -> %.5f, Delta %.5g -> %.5g, C_aa = %.3g meV^2\n', ...
    p2.tier2_iters, p1.crit, p2.crit, p1.tl.Delta, p2.tla.Delta, p2.C(1,1));
% REPORT the plan's T3.3 point (1.55 K, 0.05 T), no gate:
r = invz_solve_point(ion, 1.55, 0.05, [], o);
fprintf('plan point (1.55 K, 0.05 T) Tier1+2: converged=%d, crit=%.4g\n', r.converged, r.crit);
end

function test_tier2_C_bounded_near_boundary_slow(testCase)
% T3.3 IR safety: C stays finite as crit -> 0 (ODD blocks vanish at q = 0).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();  T = 1.2;
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
o = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'odd', true, 'odd_blocks', S);
Bc = testCase.verifyWarning(@() invz_critical(ion, T, [], o), ...
    'invz:orderedSideNoConverge');                           % para-edge estimate expected:
Cn = zeros(3,1);  cr = zeros(3,1);  db = [0.5 0.2 0.05];    % ordered side never converges
for i = 1:3                                                 % with ODD on (T2.2 finding) --
                                                              % now ASSERTED, not loose noise
    pt = invz_solve_point(ion, T, Bc + db(i), [], o);
    assumeTrue(testCase, pt.converged);
    C = invz_odd_fieldvar(ion, pt, S, T, struct());
    Cn(i) = norm(C);  cr(i) = pt.crit;
end
verifyLessThan(testCase, Cn(3)/Cn(1), 20);                  % grows but saturates
fprintf('T3.3 IR: Bc(1.2 K) = %.4f T; |C| at Bc+0.5/0.2/0.05 T = %.3g / %.3g / %.3g meV^2 (crit %.4g / %.4g / %.4g)\n', ...
    Bc, Cn, cr);
end

function test_tier2_combined_measurement_slow(testCase)
% Combined dTc and the Tier1 : Tier2 split (REPORT), at the Bx = 0.5 T proxy
% (NOT lower: small-B speckle, plan SS8). AMENDED (Task-7 finding, 2026-07-17):
% invz_critical_T cannot bracket with ODD on -- no metastable converged-PM
% window below the boundary -- so Tc here uses the deterministic PM-side
% crit(T)-extrapolation estimator introduced in test_invz_odd_retarded.m
% (T2.2 leg): same estimator for all three configurations, so the
% off/T1/T1+T2 DIFFERENCES are meaningful (common bias cancels).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW only');
ion = invz_ion();
n = 16;  [h, k, l] = ndgrid((0:n-1)/n);  qvec = [h(:) k(:) l(:)];  qvec(1,:) = [];
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, struct('dpRng', 30, 'cache', true));
S = struct('Vca', Vca, 'Vcb', Vcb, 'Vcc', Vcc, 'Jcc0', infoB.Jcc0);
Jnu0 = invz_odd_modes(Vcc, []);
Tc0f = invz_odd_zero_field(ion, struct('mode', 'off'));
o0 = struct('J0eff', infoB.Jcc0, 'Jxx0', infoB.Jaa0, 'Tc0', Tc0f);
% Shared PM-side estimator (tests/invz_odd_tc_pm_extrap.m fixture; Jf lets the
% ODD-off config pass its explicit baseline modes, ODD configs pass []).
% ADAPTATION (Task 10, documented per the amended brief): T2.2's cross-run
% "same grid points" guard compared two nearly identical configurations
% (static vs retarded, dTc ~ 0.02 mK), where a converged-set mismatch could
% masquerade as the measured shift. Here the three configurations (off / T1 /
% T1+T2) have Tc differing by ~0.2 K, so each config MUST extrapolate near
% ITS own boundary and a cross-config same-points guard is meaningless --
% only the two ODD configs are compared point-for-point (verify below); the
% fixture's internal guards, used-points report and distance FLAG cover the
% rest. Common window bracketing all three Tc. Top at 1.90 K, NOT 1.80: the
% off-config PM-side slowing band at 0.5 T is ~0.06 K wide (probed 2026-07-17:
% 1.7333/1.7667 non-converged, first converged PM point 1.8000 with crit
% +0.0253), so a 1.80 top leaves the off config a single converged point.
Tg = 1.30:1/30:1.90;
t_off = invz_odd_tc_pm_extrap(ion, 0.5, Jnu0(:), o0, Tg);
o1 = o0;  o1.odd = true;  o1.odd_blocks = S;
[t_t1, Tused_t1] = invz_odd_tc_pm_extrap(ion, 0.5, [], o1, Tg);
o2 = o1;  o2.odd_tier2 = true;
[t_t12, Tused_t12] = invz_odd_tc_pm_extrap(ion, 0.5, [], o2, Tg);
fprintf('T3.3 combined (0.5 T): Tc off %.4f, +T1 %.4f, +T1+T2 %.4f K; split %.1f%% : %.1f%%\n', ...
    t_off, t_t1, t_t12, 100*(t_off - t_t1)/max(t_off - t_t12, 1e-9), ...
    100*(t_t1 - t_t12)/max(t_off - t_t12, 1e-9));
verifyLessThanOrEqual(testCase, t_t1, t_off + 5e-3);
verifyLessThanOrEqual(testCase, t_t12, t_t1 + 5e-3);
% The Tier-2 split (6 mK) is a pure dressed-crit readout only when both ODD
% configs extrapolate from IDENTICAL grid points; the off config's grid
% points legitimately differ (0.2 K away).
verifyEqual(testCase, Tused_t12, Tused_t1);
end
