function bx = invz_critical(ion, T, Jnu_flat, opts)
%INVZ_CRITICAL Critical transverse field at temperature T (paramagnetic side).
% Bisection on pt.crit = 1 + Sigma(0) - J(0)*chi0_cc(0), which is positive in the
% paramagnet and crosses zero at the boundary. Inside the ordered phase the
% paramagnetic EMT fixed point does not exist and invz_solve_point returns
% non-finite crit; non-finite or <=0 is classified as the ordered side.
if nargin < 4, opts = struct(); end
win = [2 7];  if isfield(opts,'window'), win = opts.window; end
tol = 0.02;   if isfield(opts,'tol'),    tol = opts.tol;    end
is_ordered = @(c) ~isfinite(c) || c <= 0;
f = @(B) crit_of(ion, T, B, Jnu_flat, opts);
flo = f(win(1));  fhi = f(win(2));
assert(is_ordered(flo) && isfinite(fhi) && fhi > 0, 'invz:bracket', ...
    'Boundary not bracketed: crit(%.2fT)=%.3g, crit(%.2fT)=%.3g', win(1), flo, win(2), fhi);
lo = win(1);  hi = win(2);
while hi - lo > tol
    mid = 0.5*(lo + hi);
    if is_ordered(f(mid)), lo = mid; else, hi = mid; end
end
bx = 0.5*(lo + hi);
end

function c = crit_of(ion, T, B, Jf, opts)
pt = invz_solve_point(ion, T, B, Jf, opts);
c = pt.crit;
end
