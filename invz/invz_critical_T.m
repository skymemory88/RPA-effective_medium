function tc = invz_critical_T(ion, Bx, Jnu_flat, opts)
%INVZ_CRITICAL_T Critical temperature at fixed transverse field Bx (paramagnetic side).
% Mirror of invz_critical with the roles of T and Bx swapped: bisection on
% pt.crit = 1 + Sigma(0) - J(0)*chi0_cc(0) over temperature at fixed field.
% crit is positive in the paramagnet (high T); inside the ordered phase (low T)
% the paramagnetic EMT fixed point does not exist and invz_solve_point returns
% non-finite crit; non-finite or <=0 is classified as the ordered side.
% Use where the boundary is nearly parallel to the field axis (T near Tc0 =
% 1.74 K): a fixed-B (horizontal) cut crosses it transversally, exactly where
% the fixed-T cut of invz_critical is ill-conditioned.
% opts.window = [Tlo Thi] (K, default [1.0 2.0]): Tlo must be on the ordered
% side (i.e. Bx < Bc(Tlo)) and Thi paramagnetic. opts.tol (K, default 0.01).
% Remaining opts pass through to invz_solve_point.
% Small-B reliability floor: keep Bx >= 0.5 T. At 0.2-0.3 T (16^3 grid,
% default solver opts) the paramagnetic solve has non-convergence patches
% near the boundary; the classifier reads them as ordered and biases tc
% upward by ~0.04-0.05 K (at Bx -> 0 invz_twolevel additionally raises
% invz:degenerateDoublet). The 0 < B < 0.5 T boundary segment spans only
% ~4 mK below the zero-field Tc and is represented by the closed-form
% invz_critical_T0field endpoint instead.
if nargin < 4, opts = struct(); end
win = [1.0 2.0]; if isfield(opts,'window'), win = opts.window; end
tol = 0.01;      if isfield(opts,'tol'),    tol = opts.tol;    end
is_ordered = @(c) ~isfinite(c) || c <= 0;
f = @(T) crit_of(ion, T, Bx, Jnu_flat, opts);
flo = f(win(1));  fhi = f(win(2));
assert(is_ordered(flo) && isfinite(fhi) && fhi > 0, 'invz:bracket', ...
    ['Boundary not bracketed: crit(%.3f K) = %.3g, crit(%.3f K) = %.3g. ' ...
     'Likely Bx = %.2f T exceeds Bc(%.3f K): lower the field or the window low edge.'], ...
    win(1), flo, win(2), fhi, Bx, win(1));
lo = win(1);  hi = win(2);
while hi - lo > tol
    mid = 0.5*(lo + hi);
    if is_ordered(f(mid)), lo = mid; else, hi = mid; end
end
tc = 0.5*(lo + hi);
end

function c = crit_of(ion, T, Bx, Jf, opts)
pt = invz_solve_point(ion, T, Bx, Jf, opts);
c = pt.crit;
end
