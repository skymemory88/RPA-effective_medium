function [K, chi_bar, info] = invzt_emt_matrix(ctil, lat, opts)
%INVZT_EMT_MATRIX A2 DIRECT matrix effective-medium closure (K = ctil^-1 - chibar^-1).
%   [K, chi_bar, info] = INVZT_EMT_MATRIX(ctil, lat, opts) is the tensor
%   generalization of the scalar effective medium (invz_emt_scalar), computed as
%   a DIRECT per-frequency closure -- NOT an iteration. It takes the site-local
%   renormalized susceptibility ctil [3,3,nz] (one Cartesian 3x3 block per
%   frequency) and the LOCKED lattice struct lat (from INVZT_JQ_TENSOR; the
%   [12,12,nq] coupling tensor and BZ weights lat.w) and returns the medium
%   coupling K [3,3,nz] and the effective-medium response chi_bar [3,3,nz].
%
%   WHY DIRECT (the K-cancellation derivation). The effective-medium closure pairs
%   a lattice equation with an impurity (single-site embedded) equation:
%       chi_lat(q) = (ctil^-1 - J(q))^-1                          (lattice RPA of ctil)
%       chi_imp    = (ctil^-1 - K)^-1                             (impurity in the medium K)
%       chi_bar    = <chi_lat(q)>_BZ,sublattice   (site-diagonal 3x3 average)
%       chi_imp    == chi_bar                     (self-consistency: impurity = medium)
%   Substituting the impurity relation into the self-consistency condition,
%       (ctil^-1 - K)^-1 = chi_bar   =>   K = ctil^-1 - chi_bar^-1,
%   and the medium coupling K does NOT appear inside chi_lat/chi_bar: it CANCELS
%   exactly. So chi_bar is just the BARE tensor RPA of ctil (K absent), and K is
%   recovered afterward from the impurity relation. No fixed point to iterate --
%   one BZ average and one 3x3 (active-subspace) solve per frequency. This mirrors
%   invz_emt_scalar's closed form (the scalar A -> K -> G algebra likewise never
%   iterates), of which this is the exact 1-channel reduction (see below).
%
%   ALGORITHM.
%     (1) chi_lat(q) = the bare tensor RPA of ctil, via INVZT_CHI_RPA (page solves
%         X = ctil*(I - J*ctil)^-1 in the numerically regular left-division form
%         (I - X0*J)\X0, X0 = kron(eye(4), ctil) -- NEVER a literal inv(ctil), which
%         would be singular for a rank-deficient local response).
%     (2) chi_bar(n) = weighted BZ + sublattice average of the site-diagonal 3x3
%         blocks X(3(s-1)+(1:3), 3(s-1)+(1:3), q) of chi_lat, RAW (NOT Hermitized):
%             chi_bar(n) = sum_q lat.w(q) * (1/4) sum_s X_block(s,q).
%         Per-sublattice blocks are also kept and their S4 spread reported
%         (info.persub_spread); the sublattice-resolved cc medium is info.diag4_cc.
%     (3) K(n) = ctil(:,:,n)^-1 - chi_bar(:,:,n)^-1, computed by SOLVES on ONE
%         COMMON ACTIVE SUBSPACE. ctil may be rank-deficient (spectator/near-null
%         Cartesian channels), and a bare inverse would be ill-posed and could
%         rotate discontinuously between frequencies. So a FREQUENCY-CONSISTENT
%         projector P [3 x r] (orthonormal columns) is built ONCE as the union of
%         the active channels over the whole frequency grid -- the range of
%         sum_n Cn^H*Cn (Cn = Hermitized ctil(:,:,n)), whose null space is exactly
%         the intersection of the per-frequency null spaces -- rank-revealed at
%         opts.rank_tol (relative to the largest active singular value). Then per
%         frequency, restricting to the active subspace,
%             ctil_P = P'*ctil*P;  chibar_P = P'*chi_bar*P;
%             K = P*( (ctil_P \ I_r) - (chibar_P \ I_r) )*P';   (K = 0 on the null complement)
%         Using ONE common P (not a per-frequency rank reveal) keeps K continuous
%         across frequencies when a channel skims the rank threshold.
%
%   SCALAR REDUCTION -- SIGN LOCKED POSITIVE (v2). With one active channel, ctil =
%   -G0/D (D = 1+Sigma) and chi_bar = -G (the chi = -G convention), so
%       K = ctil^-1 - chi_bar^-1 = -D/G0 - (-1/G) = 1/G - D/G0,
%   which is EXACTLY invz_emt_scalar's med.K (G = G0/(D + K*G0) rearranged gives
%   1/G = D/G0 + K, i.e. K = 1/G - D/G0). The SAME POSITIVE relation -- NO sign
%   flip. Consequently invzt_solve_point mode 'a2' takes K_scalar(n) = +K_cc(n)
%   (positive) into invz_lambdas/invz_sigma. And -chi_bar_cc = G0*<1/(D+Jnu*G0)> =
%   med.G identically (the site-diagonal average equals the eigenmode average since
%   trace is basis invariant), so -chi_bar(3,3) matches med.G in the decoupled cc
%   limit.
%
%   HERMITICITY / TRANSPOSE RELATION (LOCKED constraint 9). Neither chi_bar nor K
%   is Hermitized. The physical single-ion tensor obeys chi(-iwn) = chi(iwn).' (the
%   TRANSPOSE, not componentwise evenness -- plan hard-math constraint 9), so at
%   wn != 0 it carries a PHYSICAL gyrotropic (~B) anti-Hermitian cross-Cartesian
%   part (the symmetry-allowed yz term; measured rel anti-Hermitian magnitude
%   ~6-17% at Bx = 0.5 T). chi_bar and K INHERIT the transpose relation exactly:
%       chi_bar(-iwn) = chi_bar(iwn).'  and  K(-iwn) = K(iwn).'
%   (because ctil(-iwn) = ctil(iwn).', the lattice satisfies J(q).' = J(-q), and a
%   grid symmetric under q -> -q lets the BZ sum relabel; K then transposes through
%   the inverse since K = ctil^-1 - chi_bar^-1). Both are Hermitian ONLY at the
%   static wn = 0 slot (where the gyrotropic part vanishes). The transpose relation
%   is EXACT only under ONE shared, frequency-consistent active-subspace projector
%   over the passed +/-wn set -- hence the union-of-ranges projector below. The cc
%   diagonal K_cc(n) stays REAL for all n (mode 'a2' uses +K_cc). info records the
%   per-frequency single-wn Hermiticity residuals as diagnostics; single-wn
%   Hermiticity is deliberately NOT asserted (it would contradict constraint 9).
%
%   OPTIONS (getf defaults):
%     rank_tol   1e-12 : active-subspace threshold, RELATIVE to the largest active
%                        singular value of the frequency-summed Gram (see above).
%     static_idx []    : index of a wn = 0 slot (where K IS Hermitian) to gate; []
%                        (default) asserts nothing on K -- correct for an arbitrary
%                        frequency grid, since K's invariant is the transpose relation
%                        (constraint 9), NOT single-frequency Hermiticity. Pass a real
%                        static index only when the caller knows that slot is wn = 0.
%
%   OUTPUTS.
%     K       [3,3,nz] : medium coupling; K(-iwn)=K(iwn).' (constraint 9), K_cc real.
%     chi_bar [3,3,nz] : RAW effective-medium response; chi_bar(-iwn)=chi_bar(iwn).'.
%     info    : struct with fields
%         projector      [3 x r] : orthonormal active-subspace basis (common to all n).
%         active_rank    r       : number of active channels.
%         rank_tol               : the threshold used.
%         persub_spread          : max_n max_s ||chi_s(n) - <chi_s(n)>_s||_F, the
%                                  sublattice (S4) spread of the per-sublattice
%                                  BZ-averaged site-diagonal blocks (a report; ~0 on
%                                  a symmetry-complete BZ grid, noisy for explicit-q lat).
%         diag4_cc       [4,nz]  : per-sublattice BZ-averaged site-diagonal cc medium
%                                  (parallels invzt_gcc_lattice's diag4; real).
%         herm_resid_chibar [1,nz]: per-frequency single-wn chi_bar Hermiticity
%                                  residual (~0 at static, O(gyrotropic) off it).
%         herm_resid_K   [1,nz]  : per-frequency single-wn K Hermiticity residual
%                                  (nonzero at wn != 0 by constraint 9). DIAGNOSTIC.
%         antisym_K      [1,nz]  : per-frequency anti-Hermitian magnitude of K.
%         Kcc            [1,nz]  : real cc diagonal medium coupling K(3,3,n).
%
%   See also INVZT_CHI_RPA, INVZT_GCC_LATTICE, INVZT_JQ_TENSOR, INVZT_SOLVE_POINT,
%   INVZ_EMT_SCALAR (the scalar-limit parity target).
if nargin < 3, opts = struct(); end
rank_tol   = getf(opts, 'rank_tol', 1e-12);
static_idx = getf(opts, 'static_idx', []);

nz = size(ctil, 3);
nq = size(lat.Jt, 3);
w  = lat.w(:);

% --- (1)+(2): chi_bar = symmetrized weighted BZ + sublattice average of the site-
%     diagonal 3x3 blocks of the bare tensor RPA of ctil. The accumulation order
%     MIRRORS the defining identity (sublattice mean, then w-weighted BZ sum).
%     Per-sublattice blocks are kept for the S4 spread report and the sublattice-
%     resolved cc medium.
chi_bar = zeros(3, 3, nz);
persub  = zeros(3, 3, 4, nz);
for k = 1:nz
    X = invzt_chi_rpa(ctil(:,:,k), lat.Jt);             % [12,12,nq], bare RPA (K cancels)
    acc = zeros(3, 3, nq);
    for s = 1:4
        blk = X(3*(s-1)+(1:3), 3*(s-1)+(1:3), :);       % [3,3,nq] site-diagonal block
        acc = acc + blk / 4;                             % sublattice mean
        persub(:,:,s,k) = sum(blk .* reshape(w, 1, 1, nq), 3);   % per-sublattice BZ average
    end
    bz = sum(acc .* reshape(w, 1, 1, nq), 3);           % w-weighted BZ sum of the sublattice mean
    chi_bar(:,:,k) = bz;                                  % RAW average (v3): do NOT Hermitize --
    % chi0(iwn) is non-Hermitian off the static slot (physical gyrotropic ~B term),
    % so chi_bar obeys the LOCKED transpose relation chi_bar(-iwn) = chi_bar(iwn).'
    % (constraint 9), Hermitian ONLY at wn=0. Hermitizing would corrupt that relation.
end

% --- frequency-consistent active-subspace projector (union over frequencies) ----
%     range(sum_n Cn^H*Cn) (Cn = Hermitized ctil) is exactly the UNION of the per-
%     frequency active subspaces, rank-revealed relative to the largest active channel
%     (shared policy with INVZT_SIGMA_TENSOR); see INVZT_ACTIVE_PROJECTOR.
P = invzt_active_projector(ctil, rank_tol);
r = size(P, 2);

% --- (3): K = ctil^-1 - chi_bar^-1 via solves on the common active subspace ------
%     K is NOT Hermitized: it carries the physical gyrotropic antisymmetric part
%     (constraint 9). K = 0 on the null complement.
K  = zeros(3, 3, nz);
Ir = eye(r);
for k = 1:nz
    if r > 0
        ctil_P   = P' * ctil(:,:,k) * P;                % [r,r] restriction (raw ctil)
        chibar_P = P' * chi_bar(:,:,k) * P;             % [r,r]
        Kp = (ctil_P \ Ir) - (chibar_P \ Ir);           % subspace solves, never inv(full ctil)
        K(:,:,k) = P * Kp * P';                          % lift to 3x3
    end
end

% --- info + diagnostics ---------------------------------------------------------
info.projector   = P;
info.active_rank = r;
info.rank_tol    = rank_tol;
% per-sublattice S4 spread of the site-diagonal blocks (report; ~0 on symmetric BZ)
sp = 0;
for k = 1:nz
    mb = mean(persub(:,:,:,k), 3);
    for s = 1:4
        sp = max(sp, norm(persub(:,:,s,k) - mb, 'fro'));
    end
end
info.persub_spread = sp;
info.diag4_cc = real(reshape(persub(3,3,:,:), 4, nz));  % sublattice-resolved cc medium
% Per-frequency single-frequency Hermiticity residuals (DIAGNOSTIC, NOT asserted):
% ~0 at the static wn=0 slot, O(|gyrotropic|) off it (constraint 9). K and chi_bar
% obey the transpose relation X(-iwn) = X(iwn).' across +/-wn slots (checked by the
% caller, which knows the frequency pairing).
info.herm_resid_chibar = reshape(max(max(abs(chi_bar - pagectranspose(chi_bar)), [], 1), [], 2), 1, nz);
info.herm_resid_K      = reshape(max(max(abs(K - pagectranspose(K)), [], 1), [], 2), 1, nz);
info.antisym_K         = reshape(max(max(abs((K - pagetranspose(K))/2), [], 1), [], 2), 1, nz);
info.Kcc               = real(reshape(K(3,3,:), 1, nz));
% NO single-wn Hermiticity assertion: chi_bar and K are legitimately non-Hermitian
% off the static slot (LOCKED constraint 9). The static_idx option is accepted for
% callers that know a slot is wn=0 and want a static Hermitian sanity gate.
if ~isempty(static_idx) && static_idx >= 1 && static_idx <= nz
    hs = max(info.herm_resid_K(static_idx), info.herm_resid_chibar(static_idx));
    assert(hs < 1e-8, 'invzt:emtHermStatic', ...
        'static-slot Hermiticity residual %.3e exceeds 1e-8.', hs);
end
end
