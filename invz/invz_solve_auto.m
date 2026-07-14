function [pt, phase, di] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Ordered-first 1/z solve at one (T, Bx) point, paramagnetic fallback.
% Tries invz_solve_point_ordered; a converged spontaneous-moment solution returns with
% phase = 1 (ferromagnet). Otherwise invz_solve_point is attempted: phase = 2 on a
% converged paramagnetic solve. phase = 0 means no usable 1/z solution: pt is then the
% non-converged paramagnetic pt when one was produced (its Sigma0 may still be of
% diagnostic value), or [] when a solver raised an expected numerical condition.
% opts passes through to both solvers (hyp, J0eff, Jxx0, ...).
%
% Error policy: ONLY invz:* identifiers (expected numerical conditions, e.g.
% invz:degenerateDoublet near Bx -> 0) are absorbed; their identifiers are returned in
% di.ordered_err / di.para_err. Any other exception is a programming/API defect and
% is RETHROWN rather than converted into "no solution".
%
% Option-A caveat (see invz_solve_point_ordered): the FM/PM handoff happens at the bare
% MEAN-FIELD boundary, which sits slightly above the 1/z critical field.
if nargin < 5, opts = struct(); end
pt = [];  phase = 0;  di = struct('ordered_err', '', 'para_err', '');
try
    pto = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts);
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
        pt = pto;  phase = 1;  return;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    di.ordered_err = err.identifier;
end
try
    pt = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
    if pt.converged && isfinite(pt.Sigma0)
        phase = 2;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    di.para_err = err.identifier;
    pt = [];
end
end
