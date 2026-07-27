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
%   Dipolar backend (Step-5 Task 6, opt-in; parameters frozen 2026-07-24):
%   opts.dipole = absent | 'bruteforce'
%   (both resolve to the unchanged brute-force path) | 'ewald' (opts.ewald must then
%   be a complete {alpha,r_cut,g_cut,boundary} struct; not synthesized here). The
%   production default remains bruteforce. Either way the Gamma-limit base Greg is
%   read from invz_jq_modes' info.Jpath_base_cc (no local MF_dipole/exchange
%   reconstruction lives in this file anymore). Under bruteforce the trust-radius
%   snap band is UNCHANGED (ksnap from dpRng/snapfac). Under ewald there is no trust
%   radius -- invz_jq_modes already evaluates every genuinely nonzero q directly, so
%   only the exact-Gamma[-equivalent] point itself (norm(k) < 1e-12) is replaced, by
%   the local path direction; P.ksnap is NaN and opts.snapfac/opts.dpRng become
%   documented no-ops for the snap logic (opts.dpRng is still accepted, but is not
%   forwarded into an Ewald invz_jq_modes call and does not affect it).
%
%   Returns:
%     P.Juni    [nq x 1]  uniform FM-mode coupling (meV), guard applied -- physical mode
%     P.Jnu     [nq x 4]  sorted branch couplings (meV), guard applied (exploratory)
%     P.snapped [nq x 1]  logical: true where the directional limit replaced the raw sum
%     P.s       [1 x nq]  cumulative path distance in INDEX (r.l.u.) coordinates
%     P.s_cart  [1 x nq]  cumulative path distance in Cartesian reciprocal Ang^-1
%     P.ksnap   scalar    trust radius used (Ang^-1); NaN under the ewald backend
%     P.dipole  struct    invz_jq_modes' info.dipole (backend/ewald/q_reduction/
%                         primitive_schema provenance), so callers can verify which
%                         backend actually produced this P
%
%   opts: .dpRng (30), .cache (true)  -- forwarded to invz_jq_modes for the raw sums
%         (bruteforce only; not forwarded under ewald, see above);
%         .snapfac (2.5)              -- trust radius = snapfac*2*pi/(dpRng*min ||a_i||),
%                                        bruteforce only, a documented no-op under ewald;
%         .dipole, .ewald             -- dipolar backend selection (see above).
if nargin < 3, opts = struct(); end
dpRng    = 30;  if isfield(opts,'dpRng'),   dpRng    = opts.dpRng;   end
useCache = ~isfield(opts,'cache') || opts.cache;
snapfac  = 2.5; if isfield(opts,'snapfac'), snapfac  = opts.snapfac; end
C  = invz_const();
nq = size(qpath, 1);

% Forward dipole/ewald by PRESENCE (design Sec.6.4). Absent backend keeps the
% EXACT legacy call shape below (byte-identical); explicit 'bruteforce' adds
% only .dipole; 'ewald' forwards cache/dipole/ewald but NOT dpRng (dpRng does
% not affect the Ewald calculation or its cache identity). invz_jq_path does
% no backend validation of its own -- invz_jq_modes is the single source of
% truth and raises invz:jqModes* on any unrecognized/incomplete request.
if ~isfield(opts, 'dipole') || isempty(opts.dipole)
    jqOpts = struct('dpRng', dpRng, 'cache', useCache);
else
    dipRaw = opts.dipole;
    if isstring(dipRaw) && isscalar(dipRaw), dipRaw = char(dipRaw); end
    if ischar(dipRaw) && isrow(dipRaw) && strcmp(dipRaw, 'ewald')
        jqOpts = struct('cache', useCache, 'dipole', opts.dipole);
        if isfield(opts, 'ewald'), jqOpts.ewald = opts.ewald; end
    else
        jqOpts = struct('dpRng', dpRng, 'cache', useCache, 'dipole', opts.dipole);
    end
end
[Jnu, info, Juni] = invz_jq_modes(ion, qpath, jqOpts);
backend = info.dipole.backend;              % 'bruteforce' | 'ewald' (design Sec.6.5)

v     = ones(4, 1) / 2;                     % uniform ferromagnetic mode (see invz_jq_modes)
Brec  = 2*pi*inv(ion.a).';                 % reciprocal basis rows: k_cart = q_rlu * Brec
if strcmp(backend, 'ewald')
    ksnap = NaN;                             % no trust radius under ewald (design Sec.6.4)
else
    ksnap = snapfac * 2*pi / (dpRng * min(vecnorm(ion.a, 2, 2)));
end
snapped = false(nq, 1);
Greg = info.Jpath_base_cc;                  % backend-agnostic Gamma-limit base (invz_jq_modes metadata)
for iq = 1:nq
    G = round(qpath(iq, :));
    if ~invz_is_gamma_equiv(G, ion.tau), continue; end
    k = (qpath(iq, :) - G) * Brec;
    if strcmp(backend, 'ewald')
        trigger = norm(k) < 1e-12;          % ewald: only exact Gamma[-equivalent] snaps
    else
        trigger = norm(k) < ksnap;           % bruteforce: unchanged trust-radius band
    end
    if ~trigger, continue; end
    if norm(k) < 1e-12                      % exactly at G: use the local path direction
        if iq < nq, dq = qpath(iq+1, :) - qpath(iq, :);
        elseif iq > 1, dq = qpath(iq, :) - qpath(iq-1, :);
        else, dq = [0 0 0];
        end
        k = dq * Brec;
        if norm(k) < 1e-12, k = [1 0 0] * Brec; end   % single point at G: in-plane default
    end
    kz2 = (k(3) / norm(k))^2;
    Jm = Greg + C.gfac*(4*pi/ion.Vc)*(1/3 - kz2);     % directional nonanalytic broadcast
    Jm = (Jm + Jm')/2;
    Jnu(iq, :) = sort(real(eig(Jm))).';
    Juni(iq)   = real(v.'*Jm*v);                      % same directional limit, uniform mode
    snapped(iq) = true;
end

P.Jnu = Jnu;  P.Juni = Juni(:);  P.snapped = snapped;  P.ksnap = ksnap;
P.s      = [0 cumsum(vecnorm(diff(qpath, 1, 1),        2, 2)).'];
P.s_cart = [0 cumsum(vecnorm(diff(qpath, 1, 1) * Brec, 2, 2)).'];
P.dipole = info.dipole;   % additive: backend/ewald/q_reduction/primitive_schema provenance
end
