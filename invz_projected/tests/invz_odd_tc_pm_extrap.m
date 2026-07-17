function [tc, Tused] = invz_odd_tc_pm_extrap(ion, B, Jf, o, Tg)
%INVZ_ODD_TC_PM_EXTRAP Shared PM-side Tc estimator at fixed B (test fixture).
%   [tc, Tused] = INVZ_ODD_TC_PM_EXTRAP(ion, B, Jf, o, Tg) samples crit on the
%   caller-supplied FIXED grid Tg (house step 1/30 K), keeps converged points
%   only (the classifier convention of invz_critical_T, voted via invz_crit_at),
%   and linearly extrapolates the two LOWEST converged PM points to crit = 0.
%   Plain fixture function, NOT a test (the invz_odd_anchors precedent --
%   runtests on tests/ does not collect it). It consolidates the previously
%   duplicated tc_pm_extrap (test_invz_odd_retarded, T2.2 decision leg) and
%   odd_tc_extrap (test_invz_odd_tier2, T3.3 combined measurement): the
%   estimator MATH is character-identical to what both shared; the union of
%   their surfaces is the Jf passthrough, the used-points report and the
%   distance FLAG.
%
%   Jf: baseline mode spectrum forwarded to invz_crit_at ([] on ODD
%   configurations; an ODD-off configuration passes its explicit modes).
%   The FIXED grid makes the estimator deterministic: two configurations voting
%   on identical T points give a pure crit-shift readout in their difference --
%   Tused (the two grid points actually used) is returned so callers can gate
%   or report same-points agreement where that is meaningful (see each caller's
%   guard/adaptation notes). Each extrapolation is internally guarded (>= 2
%   converged PM points, all crit > 0, PM-side monotonicity on the two lowest
%   points); the points used are REPORTED, and an extrapolation reaching > 2
%   grid steps below the lowest converged point is FLAGGED in the log
%   (distance sanity, report not gate).
c  = nan(size(Tg));  okv = false(size(Tg));
for i = 1:numel(Tg)
    [c(i), okv(i)] = invz_crit_at(ion, Tg(i), B, Jf, o);
end
Tv = Tg(okv);  cv = c(okv);
assert(numel(cv) >= 2, ...
    'invz_odd_tc_pm_extrap: need >= 2 converged PM points in the scan window (got %d).', numel(cv));
assert(all(cv > 0), ...
    'invz_odd_tc_pm_extrap: converged PM crit values must all be positive (got %s).', mat2str(cv));
assert(cv(2) > cv(1), ...
    ['invz_odd_tc_pm_extrap: PM-side crit must increase with T on the two lowest ' ...
     'converged points (got cv(1) = %.6g at T = %.6g, cv(2) = %.6g at T = %.6g).'], ...
    cv(1), Tv(1), cv(2), Tv(2));
Tused = Tv(1:2);                                 % two lowest converged points
tc = Tv(1) - cv(1)*(Tv(2) - Tv(1))/(cv(2) - cv(1));
fprintf('invz_odd_tc_pm_extrap(B = %.2g T): used T = [%.4f %.4f] K, crit = [%.4g %.4g], tc = %.4f K\n', ...
    B, Tused(1), Tused(2), cv(1), cv(2), tc);
if Tused(1) - tc > 2*(1/30)
    fprintf(['invz_odd_tc_pm_extrap: FLAG -- extrapolation reaches %.1f grid steps below the ' ...
        'lowest converged point (%.4f K -> %.4f K).\n'], (Tused(1) - tc)/(1/30), Tused(1), tc);
end
end
