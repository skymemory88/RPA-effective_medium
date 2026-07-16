function [c, ok] = invz_crit_at(ion, T, B, Jf, opts)
%INVZ_CRIT_AT crit = pt.crit at (T,B) and whether it is a trustworthy (converged, finite) verdict.
% Shared by invz_critical (fixed-T field scan) and invz_critical_T (fixed-B temperature scan).
try
    pt = invz_solve_point(ion, T, B, Jf, opts);
    c  = pt.crit;
    ok = pt.converged && isfinite(c);
catch
    c = NaN;  ok = false;                     % e.g. invz:degenerateDoublet
end
end
