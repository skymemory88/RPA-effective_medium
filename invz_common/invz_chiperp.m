function [Xp, info] = invz_chiperp(ion, T, B, opts)
%INVZ_CHIPERP Static (and Matsubara-ready) transverse single-ion susceptibility.
%   [Xp, info] = INVZ_CHIPERP(ion, T, B, opts) returns Xp, the real symmetric
%   2x2 (a,b)x(a,b) block (meV^-1) of the full electronuclear single-ion
%   susceptibility invz_chi0z(si, T, z, ...) at complex frequency z (default 0).
%   The (a,b) axes are the Cartesian x (field/a) and y (b) directions.
%
%   opts fields (all optional):
%     hyp           (true)         include the nuclear I = 7/2 manifold in si
%     transverse_mf ('legacy_x')   transverse mean-field mode passed to si
%     elastic       (true)         include the elastic (Curie/beta) term in chi0z
%     Jxx0          (ion.Jxx0)     transverse MF coupling forwarded to the si MF
%     z             (0)            scalar OR vector of complex frequencies. A
%                                  vector returns Xp [2,2,nz] (symmetrized and
%                                  real-parted per slice) -- the retardation form
%                                  used by the optional retarded ODD path. At z = i*wn the
%                                  elastic term is dropped automatically by chi0z's
%                                  |z| < ztol gating; wn = 0 keeps it.
%
%   info fields:
%     asym          max abs antisymmetric (gyrotropic) part of the (a,b) block,
%                   over all z slices. This chi_ab part is gyrotropic and drops
%                   out of the static Gaussian elimination -- it is DISCARDED
%                   here (Xp is symmetrized); callers assert it is small.
%     elastic_share (Xp_elastic - Xp_noelastic)/Xp_elastic on the z = 0, (1,1)
%                   element -- the fraction of chi_aa carried by the elastic
%                   (Curie/beta) term (small => Van Vleck dominated).
%     si            the converged single-ion state, so downstream callers (the
%                   EMT<->Sigma solve, Tasks 4/5) can reuse it without a second
%                   diagonalization.
%
%   DESIGN DECISION (ODD main-body plan T1.2 -- recorded verbatim): chi_perp is
%   Van Vleck-dominated and insensitive to the cc order parameter and to K -- it
%   is computed ONCE per (T, Bx) BEFORE the EMT<->Sigma loop, never
%   self-consistently inside it. This function does NOT route through
%   invz_twolevel: the transverse response crosses the ~10 K CF gap and is
%   regular where the two-level object is not; the invz:degenerateDoublet guard
%   must never fire here -- invz_chiperp works at Bx = 0.
%
%   Units: energies meV, chi in meV^-1 (invz_const conventions). No extra
%   g-factors: chi0 is meV^-1 and the ODD blocks it contracts with are meV.

if nargin < 4, opts = struct(); end
hyp   = getf(opts, 'hyp', true);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
elast = getf(opts, 'elastic', true);
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
z     = getf(opts, 'z', 0);

% Converged single-ion state ONCE, with the SAME options the rest of the point
% solve uses (T1.2: chi_perp is evaluated at that shared single-ion state).
si = invz_single_ion(ion, T, invz_field_vec(B), ...
    struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf));

% Full electronuclear chi0 at the requested frequencies; keep the (a,b) block.
c0  = invz_chi0z(si, T, z, struct('elastic', elast));    % [3,3,nz]
blk = c0(1:2, 1:2, :);                                   % [2,2,nz] complex

% Symmetrize per slice; the antisymmetric (gyrotropic) part is discarded, but
% measured and reported so callers can assert it is small.
blkT  = permute(blk, [2 1 3]);
Xsym  = (blk + blkT)/2;
Xasym = (blk - blkT)/2;
info.asym = max(abs(Xasym(:)));

% Realness: the symmetrized transverse block is real on the static and Matsubara
% axes; assert the residual imaginary content is negligible, then drop it.
mxr = max(abs(real(Xsym(:))));
mxi = max(abs(imag(Xsym(:))));
assert(mxi < 1e-10*mxr, 'invz:chiperpImag', ...
    'chi_perp imaginary part %.3g exceeds 1e-10*max|real| = %.3g', mxi, 1e-10*mxr);
Xp = real(Xsym);
if size(Xp, 3) == 1, Xp = Xp(:, :, 1); end               % clean 2x2 for scalar z

% Elastic share on the STATIC (z = 0) (1,1) element -- a fixed Van-Vleck-dominance
% diagnostic, independent of opts.z / opts.elastic.
if isscalar(z) && z == 0 && elast
    % Default case: the block computed above IS static_ab_block(si, T, true)
    % (same chi0z call, same symmetrize/real ops) -- reuse it bit-identically.
    Xel = Xp;
else
    Xel = static_ab_block(si, T, true);
end
Xnel = static_ab_block(si, T, false);
info.elastic_share = (Xel(1,1) - Xnel(1,1)) / Xel(1,1);

info.si = si;
end

function Xs = static_ab_block(si, T, elast)
%STATIC_AB_BLOCK Symmetrized real (a,b) block of chi0z at z = 0 (elastic on/off).
c  = invz_chi0z(si, T, 0, struct('elastic', elast));
b  = c(1:2, 1:2, 1);
Xs = real((b + b.')/2);
end
