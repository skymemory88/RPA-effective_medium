function [pt, phase, di] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Ordered-first 1/z solve at one (T, B) point, paramagnetic fallback.
% Bx: scalar (transverse, historical) or [Bx By Bz] vector (T).
%
% Transverse route (|Bz| <= opts.bz_tol, default 1e-9 T; the component is ZEROED in
% the dead band): tries invz_solve_point_ordered (spontaneous moment, phase = 1),
% then invz_solve_point (paramagnet, phase = 2) -- identical to the historical logic.
%
% Longitudinal route (|Bz| > bz_tol): the moment is field-induced at every (T, B), so
% the solve goes EXCLUSIVELY through invz_solve_point_ordered with forced_moment = true
% (the strict-paramagnet solver would reject m ~= 0; invz_sigma_ordered reduces to the
% paramagnet form as m -> 0). phase = 1 on acceptance; on failure phase = 0 and pt is
% the failed pto WHENEVER it carries a nonempty si (RPA-overlay fallback for the map),
% pt = [] only when no usable single-ion state exists.
%
% Error policy: only invz:* identifiers are absorbed into di; anything else rethrows.
if nargin < 5, opts = struct(); end
bz_tol = getf(opts, 'bz_tol', 1e-9);
B = invz_field_vec(Bx);
if abs(B(3)) <= bz_tol, B(3) = 0; end          % dead band: exactly transverse below tolerance
pt = [];  phase = 0;  di = struct('ordered_err', '', 'para_err', '');

if B(3) ~= 0
    oo = opts;  oo.forced_moment = true;
    try
        pto = invz_solve_point_ordered(ion, T, B, Jnu_flat, oo);
        if pto.is_ordered && pto.converged && isfinite(pto.Sigma0) && pto.si.mf_converged
            pt = pto;  phase = 1;
        elseif ~isempty(pto.si)
            pt = pto;                          % failed, but si (and maybe tl) support an overlay
        end
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        di.ordered_err = err.identifier;
    end
    return;
end

try
    pto = invz_solve_point_ordered(ion, T, B, Jnu_flat, opts);
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
        pt = pto;  phase = 1;  return;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    di.ordered_err = err.identifier;
end
try
    pt = invz_solve_point(ion, T, B, Jnu_flat, opts);
    if pt.converged && isfinite(pt.Sigma0)
        phase = 2;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    di.para_err = err.identifier;
    pt = [];
end
end
