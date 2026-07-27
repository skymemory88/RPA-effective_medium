function [tc, used] = invzt_tc_pm_extrap(critfun, Tg)
%INVZT_TC_PM_EXTRAP Handle-based PM-side Tc estimator by two-point crit extrapolation.
%   [tc, used] = INVZT_TC_PM_EXTRAP(critfun, Tg) samples the criticality handle
%   critfun on the caller-supplied FIXED grid Tg (house step 1/30 K), keeps the
%   CONVERGED PARAMAGNETIC points (ok AND crit > 0), and linearly extrapolates the
%   two LOWEST such points (in T) to crit = 0 -- the small-B-proxy ordering
%   temperature. This is the projected side's small-B-proxy extrapolation MATH,
%   RE-OWNED here with a pure handle interface (no ion/B/Jf/opts, no invz_crit_at
%   dependency), so any tensor-side criticality source can drive it.
%
%   critfun : function handle T -> [crit, ok]. crit is the criticality scalar
%             (crit > 0 in the paramagnet, < 0 ordered); ok is the sample's
%             VALIDITY (converged AND finite). The A1 usage is
%             critfun = @(T) crit_ok(invzt_solve_point(ion, T, [0.05 0 0], lat, o))
%             with the small-Bx proxy (mode 'a1' forbids B = 0, invzt:a1ZeroField).
%   Tg      : ascending temperature grid (K); the FIXED grid makes the estimator
%             deterministic -- two configurations voting on identical T points give
%             a pure crit-shift readout in their difference.
%
%   FILTER, NOT ASSERT (the key tensor re-owning): the projected fixture keeps ALL
%   ok points and ASSERTS every one has crit > 0. That is wrong for the tensor A1
%   solver, whose Anderson map CONVERGES metastable PM fixed points in the ORDERED
%   phase (ok = true, crit < 0) at low T (Task-6 finding). Here the crit > 0 filter
%   selects the paramagnetic points and DISCARDS those converged-ordered samples,
%   so the two lowest PM points bracket the boundary from above.
%
%   tc   = Tv(1) - cv(1)*(Tv(2) - Tv(1))/(cv(2) - cv(1))   (crit -> 0 line through
%          the two lowest PM points (Tv(1),cv(1)), (Tv(2),cv(2))).
%   used = Tv(1:2), the two grid temperatures actually extrapolated (reported so a
%          caller can gate/report same-points agreement between configurations).
%
%   ERRORS invzt:tcNoWindow when a valid PM extrapolation window cannot be formed:
%     - fewer than 2 converged PM points on Tg, or
%     - the two lowest PM points are not crit-monotonic (cv(2) <= cv(1)), which
%       would extrapolate the wrong way (PM-side crit increases with T).
%   A window whose crit = 0 target sits more than 2 grid steps below the lowest
%   converged PM point is FLAGGED (distance sanity, reported not gated).
%
%   See also INVZT_SOLVE_POINT, INVZT_CRITICAL, INVZ_ODD_TC_PM_EXTRAP (projected
%   fixture whose math this re-owns), INVZ_ODD_ZERO_FIELD (true-B=0 closed form).
Tg = Tg(:).';
n  = numel(Tg);
c   = nan(1, n);
okv = false(1, n);
for i = 1:n
    [c(i), okv(i)] = critfun(Tg(i));
end

% Converged PARAMAGNETIC points: valid sample AND crit > 0. The crit > 0 filter is
% what re-owns the projected math for the tensor solver (see header): it drops the
% ok-but-crit<0 metastable-PM-in-ordered-phase samples instead of erroring on them.
pm = okv & (c > 0);
Tv = Tg(pm);
cv = c(pm);

if numel(cv) < 2
    error('invzt:tcNoWindow', ...
        ['invzt_tc_pm_extrap: need >= 2 converged PM points (ok AND crit > 0) on ' ...
         'the grid to extrapolate Tc (got %d; %d converged, %d PM). Widen/lower Tg ' ...
         'into the paramagnet or check convergence.'], ...
        numel(cv), sum(okv), numel(cv));
end

% Two LOWEST PM points (Tg ascending => pm-filtered Tv/cv stay ascending in T).
Tused = Tv(1:2);
cused = cv(1:2);
if ~(cused(2) > cused(1))
    error('invzt:tcNoWindow', ...
        ['invzt_tc_pm_extrap: the two lowest converged PM points are not ' ...
         'crit-monotonic (cv(1) = %.6g at T = %.6g, cv(2) = %.6g at T = %.6g); ' ...
         'PM-side crit must increase with T for the extrapolation to be valid.'], ...
        cused(1), Tused(1), cused(2), Tused(2));
end

tc   = Tused(1) - cused(1)*(Tused(2) - Tused(1))/(cused(2) - cused(1));
used = Tused;

% --- report (stdout only; never writes a tracked file) ------------------------
fprintf('invzt_tc_pm_extrap: used T = [%.4f %.4f] K, crit = [%.4g %.4g], tc = %.4f K\n', ...
    Tused(1), Tused(2), cused(1), cused(2), tc);
if n >= 2
    step = Tg(2) - Tg(1);                       % house step (1/30 K); measured from Tg
    if isfinite(step) && step > 0 && (Tused(1) - tc) > 2*step
        fprintf(['invzt_tc_pm_extrap: FLAG -- extrapolation reaches %.1f grid steps below ' ...
            'the lowest converged PM point (%.4f K -> %.4f K).\n'], ...
            (Tused(1) - tc)/step, Tused(1), tc);
    end
end
end
