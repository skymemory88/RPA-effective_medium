function [ok, detail] = invz_task2_agree(stateA, stateB, opts)
%INVZ_TASK2_AGREE Full-state reproducibility/agreement test (stage-2c task 2a; prereg
% docs/invzp_task2_prereg.md Sec. C, FROZEN 2026-07-23 -- values below are cited VERBATIM).
%
% Compares the COMPLETE exported state -- not merely Sigma0/D_uni -- across seeds
% (cold/warm/multi) and isolated-vs-swept: Sigma(:), K(:), lambda(:), K0, Gstat, D_uni, and
% the endpoint/root quantities hstar, m_star, r_star. Per-quantity test is a COMBINED
% absolute + relative tolerance -- relative tolerance ALONE is invalid near D_uni = 0
% (prereg Sec. C prose, verbatim):
%     |a - b| <= AbsTol_q + RelTol*max(|a|,|b|),   RelTol = 1e-6 (FROZEN)
%     AbsTol_q = 1e-8 * scale_q                    (FROZEN floor magnitude)
% RelTol-alone is invalid near zero because a pure-relative bound COLLAPSES toward zero
% exactly where D_uni/Sigma/etc. legitimately live near a stability boundary, so it would
% WRONGLY FLAG two re-solves of the SAME fixed point (differing only by ordinary
% tol_outer=1e-8-level solver noise) as disagreeing -- the standard reason a bare relative
% tolerance is unsafe near zero (cf. numpy.isclose's own documented rtol-near-zero caveat).
% The AbsTol_q floor is what allows such a pair to correctly agree; test 3 below proves
% this is load-bearing, not decorative, by showing the SAME pair disagrees under a
% RelTol-alone (AbsTol_q forced to 0) bound and agrees only once the floor is restored.
%
% scale_q (prereg: "the implementer fixes each scale_q in code from the quantity's units and
% documents it; no post-hoc adjustment" -- fixed HERE, from the physics, not from the data
% being compared):
%   Sigma, D_uni   : scale = 1           (dimensionless self-energy / static pole observable)
%   K, lambda, K0  : scale = |J0eff|     (energy-dimensioned quantities, same family as J0eff)
%   Gstat          : scale = 1/|J0eff|   (inverse-energy; D_uni = 1+(J0eff-K0)*Gstat is O(1)
%                    exactly when Gstat ~ 1/J0eff -- the natural, problem-native unit, NOT a
%                    per-comparison/data-dependent quantity)
%   hstar          : scale = |J0eff|     (hz_fixed / the H_MF root is energy-dimensioned:
%                    SS9.3's spontaneous condition h0(hmf) = J0eff*<Jz>(hmf))
%   m_star         : scale = 1           (dimensionless <Jz> moment, O(1) for this ion)
%   r_star         : scale = 1           (dimensionless Green's-function ratio r = G0/Gtil0)
%
% stateA, stateB (struct, REQUIRED fields, all of): Sigma, K, lambda, K0, Gstat, D_uni,
%   hstar, m_star, r_star. Sigma/K/lambda may be vectors of matching length (elementwise
%   combined-tolerance test: EVERY element must satisfy the bound -- the natural vector
%   generalization, the same worst-case convention invz_ordered_residual.m's own
%   robust_max_abs uses); K0/Gstat/D_uni/hstar/m_star/r_star are scalars.
% opts.J0eff (REQUIRED): sets the |J0eff|-derived scale_q's above.
%
% ok (logical): true iff every quantity passes.
% detail (struct): one field per quantity, struct('diff','AbsTol','RelTol','pass'); plus
%   detail.ok (== ok) and detail.worst (name of the quantity with the largest fail margin,
%   '' when all pass).
%
% NaN/Inf convention: a quantity that is BIT-IDENTICAL between A and B (via isequaln, which
% treats matching NaN patterns as equal) passes trivially -- this is the common, legitimate
% case for the endpoint fields (hstar/m_star/r_star) when neither state is at a root, both
% consistently NaN. Any OTHER non-finite pattern (one side NaN and the other not, an Inf,
% etc.) is an automatic disagreement -- never silently compared as if finite.
if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'J0eff')
    error('invz:task2Agree', 'opts.J0eff is required (sets the energy-dimensioned scale_q terms).');
end
J0eff = opts.J0eff;

RelTol   = 1e-6;   % prereg Sec. C, FROZEN
AbsFloor = 1e-8;   % prereg Sec. C, FROZEN (AbsTol_q = 1e-8*scale_q)

scale = struct('Sigma', 1, 'D_uni', 1, 'K', abs(J0eff), 'lambda', abs(J0eff), ...
               'K0', abs(J0eff), 'Gstat', 1/abs(J0eff), 'hstar', abs(J0eff), ...
               'm_star', 1, 'r_star', 1);
qnames = {'Sigma', 'K', 'lambda', 'K0', 'Gstat', 'D_uni', 'hstar', 'm_star', 'r_star'};

detail = struct();
ok = true;
worst_name = '';  worst_margin = -Inf;
for i = 1:numel(qnames)
    q = qnames{i};
    if ~isfield(stateA, q) || ~isfield(stateB, q)
        error('invz:task2Agree', 'both states must carry field ''%s''.', q);
    end
    a = stateA.(q)(:);  b = stateB.(q)(:);
    if numel(a) ~= numel(b)
        error('invz:task2Agree', 'field ''%s'' size mismatch (%d vs %d).', q, numel(a), numel(b));
    end
    AbsTol_q = AbsFloor * scale.(q);

    if isequaln(a, b)
        pass_q = true;  dmax = 0;  margin = -Inf;
    elseif all(isfinite(a)) && all(isfinite(b))
        d   = abs(a - b);
        tol = AbsTol_q + RelTol*max(abs(a), abs(b));
        pass_q = all(d <= tol);
        dmax   = max(d);
        margin = max(d - tol);
    else
        pass_q = false;  dmax = NaN;  margin = Inf;   % any other non-finite pattern: disagree
    end

    detail.(q) = struct('diff', dmax, 'AbsTol', AbsTol_q, 'RelTol', RelTol, 'pass', pass_q);
    if ~pass_q, ok = false; end
    if margin > worst_margin, worst_margin = margin; worst_name = q; end
end
detail.ok = ok;
detail.worst = worst_name;
end
