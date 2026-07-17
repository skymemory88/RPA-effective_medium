function [c, ok] = invz_crit_at(ion, T, B, Jf, opts)
%INVZ_CRIT_AT crit = pt.crit at (T,B) and whether it is a trustworthy (converged, finite) verdict.
% Shared by invz_critical (fixed-T field scan) and invz_critical_T (fixed-B temperature scan).
% opts (and Jf) forward VERBATIM to invz_solve_point, so the ODD flags (T1.4:
% opts.odd + opts.odd_blocks, with Jf = []) reach the solver unfiltered — no
% threading needed here. ODD misuse errors (invz:odd*) and Tier-2 structural
% errors (invz:tlavg* from invz_twolevel_avg, invz:fieldvar* from
% invz_odd_fieldvar) are STRUCTURAL config/machinery errors, not an
% ordered-phase signal, and rethrow instead of being absorbed: structural
% config/machinery errors must not read as a phase verdict — absorbing them
% as c=NaN, ok=false would silently bias Bc/Tc in a production Tier-2 sweep
% by reading a misconfiguration as "ordered side". Only genuine
% non-convergence/physics signals (e.g. invz:degenerateDoublet) may read as
% the ordered side. They can only fire with opts.odd (+ opts.odd_tier2) on
% (flag off: behavior identical, the catch is unchanged).
try
    pt = invz_solve_point(ion, T, B, Jf, opts);
    c  = pt.crit;
    ok = pt.converged && isfinite(c);
catch err
    if strncmp(err.identifier, 'invz:odd', 8) || ...
       strncmp(err.identifier, 'invz:tlavg', 10) || ...
       strncmp(err.identifier, 'invz:fieldvar', 13)
        rethrow(err);
    end
    c = NaN;  ok = false;                     % e.g. invz:degenerateDoublet
end
end
