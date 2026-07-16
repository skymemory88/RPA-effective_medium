function [pt, phase, di] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Ordered-first 1/z solve at one (T, Bx) point, paramagnetic fallback.
% Tries invz_solve_point_ordered; a converged spontaneous-moment solution returns with
% phase = 1 (ferromagnet). Otherwise invz_solve_point is attempted: phase = 2 on a
% converged paramagnetic solve. phase = 0 means no usable 1/z solution: pt is the
% non-converged paramagnetic pt if one was produced, or [] otherwise.
% opts passes through to both solvers (hyp, J0eff, Jxx0, ...).
%
% Error policy: only invz:* identifiers (expected numerical conditions, e.g.
% invz:degenerateDoublet near Bx -> 0) are absorbed and recorded in
% di.ordered_err / di.para_err; any other exception is rethrown.
%
% Option-A caveat (see invz_solve_point_ordered): the FM/PM handoff happens at the
% bare mean-field boundary, which sits slightly above the 1/z critical field.
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
