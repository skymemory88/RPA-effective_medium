function [c, ok] = invz_crit_at(ion, T, B, Jf, opts)
%INVZ_CRIT_AT crit = pt.crit at (T,B) and whether it is a trustworthy (converged, finite) verdict.
% Shared by invz_critical (fixed-T field scan) and invz_critical_T (fixed-B temperature scan).
% opts (and Jf) forward VERBATIM to invz_solve_point, so the ODD flags (T1.4:
% opts.odd + opts.odd_blocks, with Jf = []) reach the solver unfiltered.
% ODD configuration errors are structural, not phase verdicts, and rethrow;
% genuine non-convergence signals may be classified as the ordered side.
try
    pt = invz_solve_point(ion, T, B, Jf, opts);
    c  = pt.crit;
    ok = pt.converged && isfinite(c);
catch err
    if strncmp(err.identifier, 'invz:odd', 8)
        rethrow(err);
    end
    c = NaN;  ok = false;                     % e.g. invz:degenerateDoublet
end
end
