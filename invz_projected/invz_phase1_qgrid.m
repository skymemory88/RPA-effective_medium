function g = invz_phase1_qgrid(ion, N, offsetFlags, convention, gammaPolicy)
%INVZ_PHASE1_QGRID Phase-1 BZ quadrature/offset/Gamma-policy grid builder (stage-2c coupling-only
% audit, ADDITIVE; docs/invzp_phase1_quadrature_prereg.md "Conventions"/"Offsets"/"Gamma policy",
% FROZEN). Calls qVec_generator (the SAME grid function invz_bz_couplings.m /
% invz_task2_couplings_shifted_grid.m call) and invz_is_gamma_equiv (the SAME Gamma-equivalence
% test invz_jq_modes.m / invz_task2_couplings_shifted_grid.m use); does NOT modify either file.
% NO field argument anywhere (Phase 1 is field-independent, prereg "Scope").
%
% CONSTRUCTION ("built by constructing the refined 2N grid once and partitioning it into the
% eight sub-offsets", prereg "Offsets"). Per axis (h,k,l) this builds exactly TWO 1-D coordinate
% arrays:
%   ax0 = this convention's own DIRECT N-point axis (a single qVec_generator call at grid=N),
%         i.e. bit-identical to what invz_bz_couplings.m/invz_task2_couplings_shifted_grid.m
%         themselves call for the "no shift" case, for BOTH conventions.
%   ax1 = ax0 + halfstep, where halfstep = (ax0(2)-ax0(1))/2 -- HALF of this convention's own
%         N-axis spacing -- the literal "{0,1/2}" phase-offset language.
% Each of the 3 axes independently selects ax0 (flag false) or ax1 (flag true); the eight
% offsetFlags combinations are exactly the eight {0,1/2}^3 offsets. The three per-axis arrays are
% combined via meshgrid in the SAME (h,k,l) -> (QX,QY,QZ) -> [QX(:) QY(:) QZ(:)] order
% qVec_generator's own generate_grid uses (grid_size(1)/(2)/(3) -> qx/qy/qz -> meshgrid(qx,qy,qz)
% -> [QX(:),QY(:),QZ(:)]), so offset [0 0 0] (ax0 on every axis) reproduces a DIRECT
% qVec_generator(ion.a,'grid',[N N N],...) call bit-for-bit -- SAME VALUES AND SAME ROW ORDER,
% since ax0 itself (from a [N 1 1] call) is computed by the identical internal formula
% (depending only on grid_size(1)=N, not on grid_size(2)/(3)) that a direct NxNxN call would use
% for its own qx. Verified in test_invz_phase1_quadrature.m (test_offset_000_matches_direct_call)
% -- no special-casing is needed in this file for the baseline offset.
%
% WHY THIS IS "the refined 2N grid, partitioned" for the convention that matters, and an honest,
% documented choice for the legacy-parity diagnostic:
%   HALF-OPEN bisects EXACTLY: ax0 = -0.5+(0:N-1)/N, and a literal refined grid at 2N points,
%   -0.5+(0:2N-1)/(2N), splits by index parity into EXACTLY ax0 (even) and ax0+1/(2N) (odd) --
%   proven algebraically and confirmed numerically against qVec_generator at N=16 (max abs diff
%   0.0 across all N points). So for halfopen, ax1 as defined above IS the literal 2N-grid's
%   odd-parity partition slice; this builder's output is bit-identical (as a value set) to
%   "build qVec_generator(2N,'endpoint',false) once and partition by index parity" for every one
%   of the eight offsets.
%   LEGACY_INCLUSIVE does NOT bisect exactly: linspace(-0.5,0.5,N) has step 1/(N-1), and
%   linspace(-0.5,0.5,2N)'s own index-parity partition is NOT a constant shift of the direct
%   N-point grid (confirmed numerically: max abs diff ~0.032 r.l.u. at N=16 between the direct
%   16-point legacy grid and the even-parity half of a literal 32-point legacy grid). Using that
%   literal 2N-grid parity partition for legacy would make the [0 0 0] slot fail to reproduce the
%   frozen pre-registration's own worked numeric example (item 1: "N^3 nominal vs (N-1)^3
%   distinct -- for N=16, 4096 vs 3375") and would make every OTHER legacy offset spuriously PASS
%   item 1 (0 duplicates) -- the historical +-0.5 duplicate faces land in different parity classes
%   of a literal 2N grid and never co-occur within one partition slot, silently HIDING, not
%   demonstrating, the defect prereg item 1 documents as "expected to FAIL". This builder instead
%   applies the identical ax0/ax1-by-constant-shift rule to both conventions. For legacy this has
%   the opposite, and on reflection more informative, property: because the direct N-point legacy
%   axis spans exactly one full BZ period (ax0(N)-ax0(1) == 1 identically, endpoint-inclusive),
%   ANY constant translation of it is STILL exactly one full period wide, so ax1 ALSO wrap-
%   collides its own first and last point. Consequence (confirmed numerically): under this
%   builder, EVERY ONE of the eight legacy offsets exhibits the endpoint-duplicate defect, not
%   just [0 0 0] -- correctly demonstrating that the defect is intrinsic to the endpoint-inclusive
%   CONVENTION itself, not an artifact of one particular offset choice. Reported per-offset,
%   honestly, in docs/invzp_phase1_report.md (item 1 table); legacy_inclusive remains "retained
%   ONLY as a labeled legacy-parity diagnostic" (prereg) throughout -- this file makes no
%   convention selection.
%
% INPUTS
%   ion          invz_ion()-shaped struct (needs .a, .tau).
%   N            positive integer grid size per axis (12, 16, 20 in the frozen ladder).
%   offsetFlags  [1x3] logical/0-1, [oh ok ol] -- true means that axis sits at the ax1 (half-step
%                shifted) phase. Use invz_phase1_offsets() to enumerate the canonical eight.
%   convention   'halfopen' | 'legacy_inclusive' (qVec_generator's 'endpoint' false/true resp.).
%   gammaPolicy  'P_complete' (keep Gamma, uniform weight 1/N^3) | 'P_drop' (remove exact-Gamma
%                rows via invz_is_gamma_equiv, renormalize to 1/(N^3-n_gamma)).
%
% OUTPUT g (struct):
%   g.qvec         [Npts x 3] reduced-q points, wrapped into the BZ via mod(q+0.5,1)-0.5.
%                  Duplicate periodic points (legacy_inclusive only) are NOT deduplicated -- kept
%                  as literal duplicate rows, exactly reproducing the over-weighted quadrature
%                  under test (this IS what item 1 measures downstream).
%   g.w            [Npts x 1] uniform weights, sum(g.w) == 1 to machine precision (P_complete:
%                  1/N^3 each, including any duplicate rows; P_drop: 1/(N^3-n_gamma) over the
%                  Gamma-filtered rows).
%   g.n_gamma      exact-Gamma row count found (via invz_is_gamma_equiv, scanning ALL nominal
%                  rows) BEFORE any P_drop filtering -- report-only under P_complete, equals the
%                  number of rows P_drop removed under P_drop.
%   g.nominal      N^3 (pre-Gamma-filter cardinality; the P_complete weight denominator always).
%   g.N, g.offsetFlags, g.convention, g.gammaPolicy   echoed provenance.
%   g.note         human-readable provenance string (mirrors invzt_qgrid.m's g.note contract).
%
% No RNG/Date anywhere: calling this twice with identical inputs yields an identical g.qvec/g.w.
if ~(isnumeric(N) && isscalar(N) && isfinite(N) && N == round(N) && N >= 2)
    error('invz:phase1Config', 'invz_phase1_qgrid: N must be a scalar integer >= 2; got %s.', mat2str(N));
end
offsetFlags = logical(offsetFlags(:).');
if ~isequal(size(offsetFlags), [1 3])
    error('invz:phase1Config', 'invz_phase1_qgrid: offsetFlags must be a 1x3 logical/0-1 vector; got size %s.', mat2str(size(offsetFlags)));
end
if ~(ischar(convention) || isstring(convention)) || ~ismember(convention, {'halfopen','legacy_inclusive'})
    error('invz:phase1Config', 'invz_phase1_qgrid: convention must be ''halfopen'' or ''legacy_inclusive''.');
end
if ~(ischar(gammaPolicy) || isstring(gammaPolicy)) || ~ismember(gammaPolicy, {'P_complete','P_drop'})
    error('invz:phase1Config', 'invz_phase1_qgrid: gammaPolicy must be ''P_complete'' or ''P_drop''.');
end

endpointFlag = strcmp(convention, 'legacy_inclusive');   % qVec_generator's own 'endpoint' name
qax = qVec_generator(ion.a, 'mode', 'grid', 'grid', [N 1 1], 'range', [-0.5 0.5], ...
    'endpoint', endpointFlag, 'verbose', false);
ax0 = qax(:,1).';                          % [1xN], this convention's own direct N-axis
halfstep = (ax0(2) - ax0(1)) / 2;          % half of this convention's own N-axis spacing
ax1 = ax0 + halfstep;

axByFlag = {ax0, ax1};
axh = axByFlag{double(offsetFlags(1)) + 1};
axk = axByFlag{double(offsetFlags(2)) + 1};
axl = axByFlag{double(offsetFlags(3)) + 1};

[QH, QK, QL] = meshgrid(axh, axk, axl);    % SAME order as qVec_generator's own generate_grid
qraw = [QH(:), QK(:), QL(:)];
qvec = mod(qraw + 0.5, 1) - 0.5;           % wrap into one BZ (prereg "Offsets")

nominal = size(qvec, 1);                   % == N^3 always
is_g = false(nominal, 1);
for i = 1:nominal
    is_g(i) = invz_is_gamma_equiv(qvec(i,:), ion.tau);
end
n_gamma = nnz(is_g);

switch gammaPolicy
    case 'P_complete'
        w = ones(nominal, 1) / nominal;
        note_g = sprintf('P_complete: Gamma kept (found %d exact-Gamma row(s) of %d nominal); uniform weight 1/%d.', ...
            n_gamma, nominal, nominal);
    case 'P_drop'
        qvec = qvec(~is_g, :);
        ndrop = nominal - n_gamma;
        if ndrop <= 0
            error('invz:phase1Config', ['invz_phase1_qgrid: P_drop left 0 points (N=%d, all %d ' ...
                'nominal rows were Gamma-equivalent).'], N, nominal);
        end
        w = ones(ndrop, 1) / ndrop;
        note_g = sprintf('P_drop: removed %d exact-Gamma row(s) of %d nominal via invz_is_gamma_equiv; renormalized to uniform weight 1/%d.', ...
            n_gamma, nominal, ndrop);
end

g.qvec        = qvec;
g.w           = w;
g.n_gamma     = n_gamma;
g.nominal     = nominal;
g.N           = N;
g.offsetFlags = offsetFlags;
g.convention  = convention;
g.gammaPolicy = gammaPolicy;
g.note = sprintf('%s, N=%d, offset=[%d %d %d]: %d nominal points (ax0/ax1 per-axis half-step construction, wrapped into one BZ). %s', ...
    convention, N, offsetFlags(1), offsetFlags(2), offsetFlags(3), nominal, note_g);
end
