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
% Error policy (spec SS5.1, narrowed from the historical 'invz:' PREFIX match): the outer
% dispatcher boundary invz_try_solver_call absorbs ONLY the reviewed recoverable identifiers
% invz_is_recoverable_solver_error whitelists, and records the exact identifier in di. Every
% other identifier -- including an unseen invz:* one -- is a wiring or programming error and
% rethrows, because a prefix match silently downgrades those to a masked "no solution".
%
% opts.transverse_mf ('legacy_x' | 'none' | 'vector_ab', default 'legacy_x'): no code change
% here -- opts is forwarded wholesale to invz_solve_point_ordered / invz_solve_point on both
% routes, so the MF model reaches every invz_single_ion call transitively.
if nargin < 5, opts = struct(); end
bz_tol = getf(opts, 'bz_tol', 1e-9);
B = invz_field_vec(Bx);
if abs(B(3)) <= bz_tol, B(3) = 0; end          % dead band: exactly transverse below tolerance
pt = [];  phase = 0;  di = struct('ordered_err', '', 'para_err', '');

if B(3) ~= 0
    oo = opts;  oo.forced_moment = true;
    [pto, completed, error_id] = invz_try_solver_call( ...
        @() invz_solve_point_ordered(ion, T, B, Jnu_flat, oo));
    if ~completed
        di.ordered_err = error_id;
    elseif pto.is_ordered && pto.converged && isfinite(pto.Sigma0) && pto.si.mf_converged
        pt = pto;  phase = 1;
    elseif ~isempty(pto.si)
        pt = pto;                              % failed, but si (and maybe tl) support an overlay
    end
    return;
end

% The acceptance tests now sit OUTSIDE the boundary: only the solver call itself is eligible
% for absorption, so a defect in reading the returned point can no longer masquerade as one.
[pto, completed, error_id] = invz_try_solver_call( ...
    @() invz_solve_point_ordered(ion, T, B, Jnu_flat, opts));
if ~completed
    di.ordered_err = error_id;
elseif pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
    pt = pto;  phase = 1;  return;
end
[ptp, completed, error_id] = invz_try_solver_call( ...
    @() invz_solve_point(ion, T, B, Jnu_flat, opts));
if ~completed
    di.para_err = error_id;
    pt = [];
else
    pt = ptp;
    if pt.converged && isfinite(pt.Sigma0)
        phase = 2;
    end
end
end
