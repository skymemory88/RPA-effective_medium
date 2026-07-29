function [dJ, d, dinfo] = invz_odd_deltaJ(Vca, Vcb, Xp)
%INVZ_ODD_DELTAJ ODD-mediated cc coupling deltaJ^{cc}(q): E1 + E4/E5 (meV).
%   [dJ, d, dinfo] = INVZ_ODD_DELTAJ(Vca, Vcb, Xp) contracts the geometric
%   off-diagonal dipolar (ODD) blocks with the transverse single-ion
%   susceptibility and enforces the self-site subtraction:
%
%     Vca, Vcb [4,4,nq] : J^{ca}(q), J^{cb}(q) from invz_odd_blocks (meV;
%                         individually NOT Hermitian -- never Hermitize them).
%     Xp       [2,2]    : real symmetric (a,b)x(a,b) transverse susceptibility
%                         block from invz_chiperp (meV^-1).
%
%   E1 (per q; mediated coupling):
%     dJpre(q) = Vca*Xp(1,1)*Vca' + Vca*Xp(1,2)*Vcb'
%              + Vcb*Xp(2,1)*Vca' + Vcb*Xp(2,2)*Vcb'
%     then Hermitize (M+M')/2. Since Xp is symmetric PSD, dJpre =
%     [Vca Vcb]*(Xp kron I4)*[Vca Vcb]' is PSD BY CONSTRUCTION -- the
%     Hermitize only cleans rounding (dinfo.presub_min_eig verifies this).
%   E4 (self-site subtraction, DIAGONAL ONLY):
%     dJ(s,s,q) = dJpre(s,s,q) - mean_q dJpre(s,s,q);  off-diagonal untouched.
%   E5 (uniform reduction, meV, >= 0):
%     d = mean_s mean_q dJpre(s,s,q).  The four per-sublattice means agree to
%     1e-10 relative (S4 symmetry) -- asserted here ('invz:oddS4'); the
%     diagonal of the Hermitized dJpre is asserted real to machine precision
%     before real() is taken ('invz:oddImag').
%
%   CALLER CONTRACT (grid): the qvec behind Vca/Vcb must be a FULL UNIFORM BZ
%   MESH (e.g. n = 16; [h,k,l] = ndgrid((0:n-1)/n); qvec = [h(:) k(:) l(:)]) --
%   the E4/E5 "BZ averages" are plain means over the supplied grid. Driver
%   grids EXCLUDING Gamma are fine: deltaJ(q -> 0) -> 0 on-axis and Gamma is a
%   single point of the mesh, so including/excluding it changes the averages
%   only at O(1/Nq) -- the same discretization order as the mesh itself. Do
%   NOT feed q-paths or partial meshes: the subtraction would then remove a
%   path average, not a BZ average, and E4/E5 lose their meaning.
%
%   BOOKKEEPING -- NO DOUBLE COUNTING (ODD main-body T1.3; the theory is retained in
%   jensen_1z_framework.html Section 11 and the operational contract below):
%     The on-site constant subtracted by E4 multiplies (sigma^z_i)^2 = 1 in
%     the strict two-level limit -- it is a PURE ENERGY SHIFT with no Tc
%     content, which is why it is removed from the grid matrices here and
%     re-applied EXACTLY ONCE as the explicit -d on info.Jcc0 (E5, the
%     invz_jq_modes opts.odd path). Its PHYSICAL residue -- the internal
%     transverse field renormalizing the two-level parameters (Delta, M^2,
%     n01) -- is exactly what Tier 2 owns (plan section 1.2). Tier 1 must
%     never add an on-site deltaJ(ii) term, and Tier 2 must never re-subtract
%     d: cite (E4)/(E5) and this plan when touching either side.
%
%   OUTPUTS
%     dJ    [4,4,nq] Hermitian per q: POST-subtraction deltaJ (its diagonal
%           already carries -d plus the q-structure; BZ-average of each
%           sublattice diagonal is 0 to machine precision).
%     d     scalar (meV): E5 uniform reduction, applied by the caller to
%           Jcc0/J0eff exactly once. d = 0 exactly when Xp = 0.
%     dinfo struct (computed only when requested, nargout >= 3 -- the per-q eig
%           sweep dominates this function's cost; the invz:oddImag / invz:oddS4
%           guards above stay unconditional):
%       .d_per_sublattice [4,1] pre-subtraction per-sublattice BZ means (E5)
%       .presub_min_eig   min over q of min eig of dJpre (PSD check, >~ -eps)
%       .postsub_diag_bzavg max_s |mean_q dJ(s,s,:)| (E4 exactness, ~ 0)
%       .dJ_max           max |element| of the PRE-subtraction dJpre (meV)
%
%   See also INVZ_ODD_BLOCKS, INVZ_CHIPERP, INVZ_JQ_MODES.

% --- input guards (cheap; prevent silent garbage entering E1) ---
if ~isnumeric(Vca) || ~isnumeric(Vcb) ...
        || size(Vca,1) ~= 4 || size(Vca,2) ~= 4 || ~isequal(size(Vca), size(Vcb))
    error('invz:oddBlocks', ...
        'Vca/Vcb must both be [4,4,nq] ODD blocks from invz_odd_blocks (got %s vs %s).', ...
        mat2str(size(Vca)), mat2str(size(Vcb)));
end
if ~isnumeric(Xp) || ~isequal(size(Xp), [2 2]) || ~isreal(Xp) ...
        || ~all(isfinite(Xp(:))) || ~issymmetric(Xp)
    error('invz:oddXp', ...
        'Xp must be a [2 2] real symmetric finite matrix (meV^-1, invz_chiperp output).');
end
nq = size(Vca, 3);

% --- E1: mediated coupling per q, Hermitized (PSD by construction, see header)
dJpre = Xp(1,1)*pagemtimes(Vca, 'none', Vca, 'ctranspose') ...
      + Xp(1,2)*pagemtimes(Vca, 'none', Vcb, 'ctranspose') ...
      + Xp(2,1)*pagemtimes(Vcb, 'none', Vca, 'ctranspose') ...
      + Xp(2,2)*pagemtimes(Vcb, 'none', Vcb, 'ctranspose');
dJpre = (dJpre + pagectranspose(dJpre))/2;

% --- E5: per-sublattice BZ averages of the pre-subtraction diagonal ---
dg = zeros(4, nq);              % complex-safe staging of the diagonals
for s = 1:4
    dgs = squeeze(dJpre(s, s, :)).';
    dgi = max(abs(imag(dgs)));  % Hermitized diagonal is exactly real
    assert(dgi <= 1e-12*max(max(abs(real(dgs))), 1e-30), 'invz:oddImag', ...
        'E5: imaginary content %.3g on the dJpre diagonal (sublattice %d) is not machine-small.', ...
        dgi, s);
    dg(s, :) = real(dgs);
end
m_s = mean(dg, 2);              % [4,1] pre-subtraction per-sublattice means
spread = max(m_s) - min(m_s);
scale  = max(abs(m_s));
assert(spread <= 1e-10*max(scale, 1e-30), 'invz:oddS4', ...
    ['E5: per-sublattice BZ-averaged diagonals disagree beyond 1e-10 relative ' ...
     '(spread %.3g, scale %.3g meV) -- S4 symmetry violated; qvec must be a ' ...
     'full uniform BZ mesh (see caller contract in the header).'], spread, scale);
d = mean(m_s);                  % E5 uniform reduction (>= 0 since dJpre is PSD)

% --- E4: self-site subtraction, diagonal only ---
dJ = dJpre;
for s = 1:4
    dJ(s, s, :) = dJ(s, s, :) - m_s(s);
end

% --- diagnostics (dinfo only; gated so the hot [dJ, d] callers skip the per-q
% eig sweep -- the guards above are unconditional and use none of these) ---
if nargout >= 3
    dinfo.d_per_sublattice = m_s;
    dinfo.dJ_max = max(abs(dJpre(:)));
    mn = Inf;
    for iq = 1:nq
        mn = min(mn, min(real(eig(dJpre(:,:,iq)))));   % Hermitian -> real eigs
    end
    dinfo.presub_min_eig = mn;
    pv = 0;
    for s = 1:4
        pv = max(pv, abs(mean(real(squeeze(dJ(s, s, :))))));
    end
    dinfo.postsub_diag_bzavg = pv;
end
end
