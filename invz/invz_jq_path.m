function P = invz_jq_path(ion, qpath, opts)
%INVZ_JQ_PATH Coupling branches along a q-path with a direction-aware Gamma-limit guard.
%   P = invz_jq_path(ion, qpath) evaluates the 4 cc coupling branches (invz_jq_modes
%   convention, ascending sort) at each row of qpath (r.l.u.), then REPLACES values inside
%   a trust radius of any Gamma-equivalent point G by the exact directional long-wavelength
%   limit. Rationale: MF_dipole's sharp real-space cutoff cannot resolve |q - G| below
%   ~1/(dpRng*a) -- the raw branch collapses on approach to G and then jumps at the exact
%   endpoint where the Lorentz term is added (verified: max branch 0.0016 meV at h = 1.999
%   vs the correct 0.0064 at dpRng 30). The directional limit is
%       eig( J_reg(Gamma) + gfac*(4*pi/Vc)*(1/3 - khat_z^2) )   [scalar sublattice broadcast]
%   with khat the Cartesian direction of q - G (Cohen-Keffer nonanalytic term). For any
%   in-plane approach (khat_z = 0, e.g. the (h,0,0) R2007 Fig-3 path) this equals the
%   uniform-mode Lorentz value, i.e. the info.Jcc0 convention -- continuous endpoint by
%   construction. Points at exactly G use the LOCAL PATH DIRECTION; a single-point path
%   at G defaults to the in-plane (uniform-mode) convention.
%
%   Scope: the guard covers Gamma-EQUIVALENT points only. Near non-Gamma-equivalent
%   integer points (e.g. (1,0,0), structure factor 0) the staggered branches retain
%   truncated-sum quality -- exploratory. Branch index = sorted-eigenvalue position per q;
%   branch identity is NOT tracked through crossings.
%
%   Returns:
%     P.Jnu     [nq x 4]  branch couplings (meV), guard applied
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

Jnu = invz_jq_modes(ion, qpath, struct('dpRng', dpRng, 'cache', useCache));

Brec  = 2*pi*inv(ion.a).';                 % reciprocal basis rows: k_cart = q_rlu * Brec
ksnap = snapfac * 2*pi / (dpRng * min(vecnorm(ion.a, 2, 2)));
snapped = false(nq, 1);
Greg = [];                                  % regular Gamma-point matrix, computed lazily once
for iq = 1:nq
    G = round(qpath(iq, :));
    if ~is_gamma_equiv(G, ion.tau), continue; end
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
        dip0 = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
        ex0  = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
        Greg = -squeeze(C.gfac*dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));
    end
    Jm = Greg + C.gfac*(4*pi/ion.Vc)*(1/3 - kz2);     % directional nonanalytic broadcast
    Jm = (Jm + Jm')/2;
    Jnu(iq, :) = sort(real(eig(Jm))).';
    snapped(iq) = true;
end

P.Jnu = Jnu;  P.snapped = snapped;  P.ksnap = ksnap;
P.s      = [0 cumsum(vecnorm(diff(qpath, 1, 1),        2, 2)).'];
P.s_cart = [0 cumsum(vecnorm(diff(qpath, 1, 1) * Brec, 2, 2)).'];
end

function tf = is_gamma_equiv(q, tau)
tf = abs(real(sum(exp(2i*pi*(tau*q.'))))/size(tau,1) - 1) < 1e-9;
end
