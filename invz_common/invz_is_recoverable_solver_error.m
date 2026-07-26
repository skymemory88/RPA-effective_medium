function tf = invz_is_recoverable_solver_error(id)
%INVZ_IS_RECOVERABLE_SOLVER_ERROR Single classifier for every catch on the strict-medium path
% (spec SS5.1). Returns true ONLY for identifiers that are genuine physics/branch/domain
% signals; everything else -- including unseen invz:* identifiers -- is treated as a wiring or
% programming error and must be rethrown.
%
% WHY A WHITELIST, NOT A PREFIX MATCH. invz_ordered_node_solve.m:213-216 currently absorbs
% every 'invz:*' id. Its docstring justifies that by the premise that the only error() site in
% its whole chain (invz_emt_scalar / invz_emt_static_ordered / invz_gstat_ordered /
% invz_lambdas / invz_sigma_ordered / invz_sigma) is invz:emtJnu. The strict scheme adds throw
% sites in that same chain, so a prefix match would silently downgrade a wiring error to "node
% not accepted" -- a masked column. That is exactly the failure mode that let the original
% masking defect hide for a whole stage. And there are at least three such absorbers on the path
% (invz_spectra_map x2, invz_solve_auto) -- invz_ordered_residual's own per-block safe_eval
% absorber was REMOVED by task 9 (exceptions now escape that checker unconditionally) -- so
% narrowing one only relocates the swallow: every catch must use THIS predicate.
%
% ADDING AN IDENTIFIER HERE IS A REVIEWED CONTRACT CHANGE, never a convenience response to a
% failing run. Strict-medium domain outcomes deliberately return STATUSES rather than throwing
% (spec SS5.2), so the inner node map and residual checker should have no expected recoverable
% throw at all.
%
% Whitelist rationale:
%   invz:orderedPhase      the strict-PM solver's m = 0 branch gate under a longitudinal tilt
%                          (invz_spectra_map.m:291) -- a branch signal, not a defect.
%   invz:degenerateDoublet the two-level domain floor Delta < 1e-4 meV
%                          (invz_twolevel_ordered.m:19) -- a domain signal, not a defect.
recoverable = {'invz:orderedPhase', 'invz:degenerateDoublet'};
tf = (ischar(id) || isstring(id)) && isscalar(string(id)) && ...
     any(strcmp(char(id), recoverable));
end
