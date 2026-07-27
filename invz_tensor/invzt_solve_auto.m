function [pt, phase, di] = invzt_solve_auto(ion, T, B, lat, opts)
%INVZT_SOLVE_AUTO Stability-based tensor a1 phase dispatcher at one (T, B) point.
%   [pt, phase, di] = INVZT_SOLVE_AUTO(ion, T, B, lat, opts) assigns the phase
%   by STABILITY, not by moment onset (review P0-1, 2026-07-19): the bare-MF
%   ordered leg keeps m0 > 0 well onto the PM side (measured 0.1 K: m0 = 1.17
%   at 4.8 T vs the corrected tensor QPT between 4.64 and 4.65 T), so the projected
%   ordered-first pattern would misassign [Bc, ~5.0 T]. Here the PM leg decides
%   FIRST via the three-part validity rule -- its crit > 0 IS the tensor QPT
%   criterion -- and the ordered leg is consulted only when the PM sample is
%   invalid:
%     phase 2 : PM   accepted -- converged, finite crit > 0, Sigma0 >= sigma_floor.
%     phase 1 : PM invalid AND ordered accepted -- is_ordered, converged,
%               si.mf_converged, finite Sigma0, crit > crit_tol_ordered
%               (default -1e-3): INVZT_CRIT_STATIC at the ordered renormalized
%               ctil0 measures local stability of the broken-symmetry state; a
%               converged ordered point failing it is NOT labeled a phase.
%     phase 0 : neither accepted -- near-Bc Option-A window or solver failure;
%               pt carries the last usable attempt ([] if none), di says why.
%
%   TRANSVERSE ONLY: invzt:orderedLongitudinal is RETHROWN, not absorbed (no
%   tensor forced-moment route -- caller error, not a phase outcome).
%
%   di (review P2-1): STRUCTURED per-leg diagnostics -- di.para / di.ordered,
%   each with attempted, accepted, converged, m0 (NaN for PM), crit, Sigma0,
%   err (exception identifier, '' if returned normally). A leg that returned
%   but was rejected is fully described -- rejection reasons are read off the
%   numbers (e.g. converged = true, crit < 0 = the spurious below-Bc PM point).
%
%   INPUT POLICY (validated AT ENTRY, before either leg -- round-2 P1-5):
%     - B is normalized (invz_field_vec); |Bz| > bz_tol raises
%       invzt:orderedLongitudinal HERE, so a PM-valid point can never silently
%       accept a longitudinal field.
%     - opts.mode must be 'a1' (invzt:autoMode otherwise): the phase is NEVER
%       classified by a3d -- evaluate a3d separately on accepted ordered points.
%
%   Error policy: ENUMERATED ALLOWLIST -- only the recoverable physics outcomes
%   invz:degenerateDoublet, invzt:a1ZeroField, invz:orderedPhase are absorbed
%   into .err; every other identifier (config, shape, unsupported mode,
%   invariant violations) RETHROWS.
%
%   See also INVZT_SOLVE_POINT, INVZT_SOLVE_POINT_ORDERED, INVZ_SOLVE_AUTO
%   (projected reference -- NOTE its ordered-first rule is NOT ported, P0-1).
if nargin < 5, opts = struct(); end
sfloor = getf(opts, 'sigma_floor', -0.5);            % single-sourced validity floor
ctolo  = getf(opts, 'crit_tol_ordered', -1e-3);      % ordered stability tolerance
bztol  = getf(opts, 'bz_tol', 1e-9);

% --- entry validation (round-2 P1-5): field and mode, BEFORE any solve ----------
B = invz_field_vec(B);
if abs(B(3)) > bztol
    error('invzt:orderedLongitudinal', ['invzt_solve_auto is the TRANSVERSE ' ...
        '(spontaneous-moment) route only: got Bz = %.3g T > bz_tol = %.3g T. ' ...
        'No tensor forced-moment route exists (2026-07-19 scope decision).'], B(3), bztol);
end
B(3) = 0;                                            % dead band: exactly transverse
if ~strcmp(char(getf(opts, 'mode', 'a1')), 'a1')
    error('invzt:autoMode', ['invzt_solve_auto classifies the phase with mode ' ...
        '''a1'' ONLY (got ''%s''): a3d is an after-classification evaluation on ' ...
        'accepted ordered points, never the phase criterion.'], ...
        char(getf(opts, 'mode', 'a1')));
end
RECOVERABLE = {'invz:degenerateDoublet', 'invzt:a1ZeroField', 'invz:orderedPhase'};
leg0 = struct('attempted', false, 'accepted', false, 'converged', false, ...
              'm0', NaN, 'crit', NaN, 'Sigma0', NaN, 'err', '');
di = struct('para', leg0, 'ordered', leg0);
pt = [];  ptp = [];  phase = 0;

% --- PM leg first: its crit > 0 is the QPT criterion ----------------------------
di.para.attempted = true;
try
    ptp = invzt_solve_point(ion, T, B, lat, opts);
    di.para.converged = ptp.converged;
    di.para.crit = ptp.crit;  di.para.Sigma0 = ptp.Sigma0;
    if ptp.converged && isfinite(ptp.crit) && ptp.crit > 0 && ptp.Sigma0 >= sfloor
        di.para.accepted = true;
        pt = ptp;  phase = 2;  return;
    end
    pt = ptp;                                        % keep for diagnostics
catch err
    if ~ismember(err.identifier, RECOVERABLE), rethrow(err); end
    di.para.err = err.identifier;
end

% --- ordered leg: only when the PM sample is invalid ----------------------------
% The ordered single-ion state must be expanded about Jensen's modified molecular
% field, not the bare MF field.  The full H_MF <-> H integral is non-local in the
% longitudinal field; at a continuous transition its linearized form is exact:
% J0z_mf = J0z*chi_tilde_PM(0)/chi_full_PM(0). This uses the PM solver's ACTUAL
% local convention, chi_tilde = cdom/(1+Sigma0)+crest, rather than incorrectly
% assuming that the full electronuclear response is divided by 1+Sigma0.
% Passing the rejected PM point supplies that boundary-matched factor (and an
% excellent Sigma seed) to the ordered leg.
% This removes the old mismatch in which the PM instability occurred near 4.65 T
% while the ordered eigenstates still carried the large bare-MF moment associated
% with the ~5 T MF boundary.
di.ordered.attempted = true;
try
    oo = opts;
    if ~isempty(ptp) && ptp.converged && isfinite(ptp.Sigma0) && (1 + ptp.Sigma0) > 0
        oo.hmf_sigma0 = ptp.Sigma0;
        split0 = struct('elastic', true);
        if isfield(ptp, 'mspec') && strcmp(getf(ptp.mspec, 'selection', ''), 'fixed_rank')
            split0.dominant_count = ptp.mspec.ndom;
        elseif isfield(ptp, 'mspec') && strcmp(getf(ptp.mspec, 'selection', ''), 'energy')
            split0.Esplit = ptp.mspec.Esplit;
        end
        [cdom0, crest0] = invzt_chi0_split(ptp.si, T, 0, split0);
        if ~ptp.chi_rest, crest0 = zeros(size(crest0)); end
        cfull0 = invz_chi0z(ptp.si, T, 0, struct('elastic', true));
        ctilde0 = cdom0/(1 + ptp.Sigma0) + crest0;
        oo.hmf_J0z = lat.info.Jcc0 * real(ctilde0(3,3,1)) / real(cfull0(3,3,1));
        if isfield(ptp, 'Sigma') && all(isfinite(ptp.Sigma))
            oo.Sigma_seed = ptp.Sigma;
        end
    end
    pto = invzt_solve_point_ordered(ion, T, B, lat, oo);
    di.ordered.converged = pto.converged;
    di.ordered.m0 = pto.m0;  di.ordered.crit = pto.crit;  di.ordered.Sigma0 = pto.Sigma0;
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0) ...
            && pto.si.mf_converged && pto.crit > ctolo
        di.ordered.accepted = true;
        pt = pto;  phase = 1;  return;
    end
    if ~isempty(pto.si), pt = pto; end               % last usable attempt wins the pt slot
catch err
    if ~ismember(err.identifier, RECOVERABLE), rethrow(err); end
    di.ordered.err = err.identifier;
end
end
