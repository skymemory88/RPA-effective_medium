function P = invz_jq_path(ion, qpath, opts)
%INVZ_JQ_PATH Coupling branches along a q-path with a direction-aware Gamma-limit guard.
%   Evaluates the 4 cc coupling branches (invz_jq_modes convention, ascending sort) at
%   each row of qpath (r.l.u.), then REPLACES values inside a trust radius of any
%   Gamma-equivalent point G by the exact directional long-wavelength limit
%       eig( J_reg(Gamma) + gfac*(4*pi/Vc)*(1/3 - khat_z^2) )   [scalar sublattice broadcast]
%   with khat the Cartesian direction of q - G (Cohen-Keffer nonanalytic term). This
%   avoids the raw branch collapsing on approach to G and jumping at the exact endpoint,
%   since MF_dipole's real-space cutoff cannot resolve |q - G| below ~1/(dpRng*a). For
%   in-plane approach (khat_z = 0, e.g. the (h,0,0) R2007 Fig-3 path) the limit equals the
%   uniform-mode Lorentz value (info.Jcc0 convention), giving a continuous endpoint.
%   Points at exactly G use the LOCAL PATH DIRECTION; a single-point path at G defaults
%   to the in-plane (uniform-mode) convention.
%
%   Scope: the guard covers Gamma-EQUIVALENT points only. Near non-Gamma-equivalent
%   integer points (e.g. (1,0,0), structure factor 0) the branches retain truncated-sum
%   quality. Branch index = sorted-eigenvalue position per q; branch identity is NOT
%   tracked through crossings.
%
%   P.Juni is the PHYSICAL single-mode dispersion: the uniform ferromagnetic-mode
%   projection v'*Jcc(q)*v (see invz_jq_modes), guard applied. Prefer it over any P.Jnu
%   column for a q-path dispersion: max(eig) selects the wrong sublattice branch for
%   h < 1.5 and mirrors the curve about h = 1.5. The four P.Jnu branches remain
%   available for exploratory branch-resolved views.
%
%   Returns:
%     P.Juni    [nq x 1]  uniform FM-mode coupling (meV), guard applied -- physical mode
%     P.Jnu     [nq x 4]  sorted branch couplings (meV), guard applied (exploratory)
%     P.snapped [nq x 1]  logical: true where the directional limit replaced the raw sum
%     P.s       [1 x nq]  cumulative path distance in INDEX (r.l.u.) coordinates
%     P.s_cart  [1 x nq]  cumulative path distance in Cartesian reciprocal Ang^-1
%     P.ksnap   scalar    trust radius used (Ang^-1)
%
%   opts: .dpRng (30), .cache (true)  -- forwarded to invz_jq_modes for the raw sums;
%         .snapfac (2.5)              -- trust radius = snapfac*2*pi/(dpRng*min ||a_i||).
if nargin < 3, opts = struct(); end
dpRng    = 30;  if isfield(opts,'dpRng'),   dpRng    = opts.dpRng;   end
useCache = ~isfield(opts,'cache') || opts.cache;
snapfac  = 2.5; if isfield(opts,'snapfac'), snapfac  = opts.snapfac; end
C  = invz_const();
nq = size(qpath, 1);

[Jnu, info, Juni] = invz_jq_modes(ion, qpath, struct('dpRng', dpRng, 'cache', useCache));

v     = ones(4, 1) / 2;                     % uniform ferromagnetic mode (see invz_jq_modes)
Brec  = 2*pi*inv(ion.a).';                 % reciprocal basis rows: k_cart = q_rlu * Brec
ksnap = snapfac * 2*pi / (dpRng * min(vecnorm(ion.a, 2, 2)));
snapped = false(nq, 1);
Greg = [];                                  % regular Gamma-point matrix, computed lazily once
for iq = 1:nq
    G = round(qpath(iq, :));
    if ~invz_is_gamma_equiv(G, ion.tau), continue; end
    k = (qpath(iq, :) - G) * Brec;
    if norm(k) >= ksnap, continue; end
    if norm(k) < 1e-12                      % exactly at G: use the local path direction
        if iq < nq, dq = qpath(iq+1, :) - qpath(iq, :);
        elseif iq > 1, dq = qpath(iq, :) - qpath(iq-1, :);
        else, dq = [0 0 0];
        end
        k = dq * Brec;
        if norm(k) < 1e-12, k = [1 0 0] * Brec; end   % single point at G: in-plane default
    end
    kz2 = (k(3) / norm(k))^2;
    if isempty(Greg)
        % Reuse the q-independent geometry invz_jq_modes already built (info.geomD/geomX)
        % instead of re-deriving it from scratch: bit-identical (MF_dipole/exchange 5-arg
        % reuse form), pinned by test_invz_dipole_geometry_reuse.
        dip0 = MF_dipole([0 0 0], dpRng, ion.a, ion.tau, info.geomD);
        ex0  = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau, info.geomX);
        Greg = -squeeze(C.gfac*dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));
    end
    Jm = Greg + C.gfac*(4*pi/ion.Vc)*(1/3 - kz2);     % directional nonanalytic broadcast
    Jm = (Jm + Jm')/2;
    Jnu(iq, :) = sort(real(eig(Jm))).';
    Juni(iq)   = real(v.'*Jm*v);                      % same directional limit, uniform mode
    snapped(iq) = true;
end

P.Jnu = Jnu;  P.Juni = Juni(:);  P.snapped = snapped;  P.ksnap = ksnap;
P.s      = [0 cumsum(vecnorm(diff(qpath, 1, 1),        2, 2)).'];
P.s_cart = [0 cumsum(vecnorm(diff(qpath, 1, 1) * Brec, 2, 2)).'];
end
