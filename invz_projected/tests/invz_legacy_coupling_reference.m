function ref = invz_legacy_coupling_reference()
%INVZ_LEGACY_COUPLING_REFERENCE Frozen, cache-free, test-only oracle for the
% pre-Step-5 brute-force invz_jq_modes (non-ODD branch) and invz_jq_path
% arithmetic (Ewald Step-5 Task 1 -- docs/superpowers/plans/2026-07-24-
% ewald-step5-integration.md).
%
%   ref = INVZ_LEGACY_COUPLING_REFERENCE() returns function-handle fields
%     ref.jq_modes(ion, qvec, opts) -> [Jnu, info, Juni]
%     ref.jq_path(ion, qpath, opts) -> P
%   replicating, statement-for-statement, invz_projected/invz_jq_modes.m's
%   non-ODD code path and invz_projected/invz_jq_path.m as they exist at the
%   Step-5 base commit (a58bdb4): same call order to MF_dipole/exchange, same
%   Hermitization order (Lorentz/broadcast added BEFORE Hermitization, exactly
%   as in the sources), same sort(real(eig(.)))/uniform-mode-projection
%   convention, same directional Gamma-limit guard and trust-radius logic.
%
%   This file is the FROZEN ORACLE that later Step-5 tasks regress the real
%   (soon-to-be-modified) invz_jq_modes/invz_jq_path against. Consequently:
%     * it MUST NOT call the production invz_jq_modes or invz_jq_path (that
%       would make the oracle track whatever those files are changed to,
%       defeating its purpose as an independent regression reference);
%     * it MAY call, and does call, the UNCHANGED production primitives
%       MF_dipole, exchange, invz_const, invz_is_gamma_equiv -- Step 5 does
%       not modify any of those;
%     * it is FULLY CACHE-FREE: it never reads or writes a jq4_*/jq5_* cache
%       file under any circumstances. An opts.cache field, if supplied, is
%       accepted only for call-signature parity with the production
%       functions and is otherwise ignored (ref.jq_modes/ref.jq_path always
%       recompute from the primitives).
%     * it covers ONLY the non-ODD branch of invz_jq_modes (opts.odd
%       absent/false): Step 5 does not touch the ODD diversion, and this
%       reference is not a general-purpose invz_jq_modes reimplementation.
%       ref.jq_modes raises invz:legacyRefOddUnsupported if opts.odd is
%       present and not false/empty, rather than silently ignoring it.
%
% Plain fixture/helper function (NOT a test; runtests on tests/ does not
% collect it -- the invz_ewald_metrics.m/invz_ewald_fixtures.m precedent).
ref = struct();
ref.jq_modes = @local_jq_modes_ref;
ref.jq_path  = @local_jq_path_ref;
end

% =====================================================================
% frozen reference: invz_jq_modes.m non-ODD branch (cache-free)
% =====================================================================
function [Jnu, info, Juni] = local_jq_modes_ref(ion, qvec, opts)
% Cache-free copy of invz_projected/invz_jq_modes.m's non-ODD arithmetic
% (its lines below the opts.odd diversion). Call order, the +lorz-before-
% Hermitization order, and the eig/uniform-projection convention are copied
% exactly; only cache read/write is removed.
if nargin < 3, opts = struct(); end
if isfield(opts, 'odd') && ~isempty(opts.odd) && ~isequal(opts.odd, false)
    error('invz:legacyRefOddUnsupported', ...
        ['invz_legacy_coupling_reference covers only the non-ODD invz_jq_modes ' ...
         'branch (Step 5 does not touch the ODD diversion); opts.odd must be ' ...
         'absent or false.']);
end
dpRng = 30;  if isfield(opts,'dpRng'), dpRng = opts.dpRng; end
% useCache is accepted for call-signature parity only; this reference never
% reads or writes a cache file regardless of its value.

C = invz_const();
demag = 0;   if isfield(ion,'demag')  && ~isempty(ion.demag),  demag = ion.demag;  end
if isfield(opts,'demag') && ~isempty(opts.demag), demag = opts.demag; end
alpha = 1;   if isfield(ion,'alpha')  && ~isempty(ion.alpha),  alpha = ion.alpha;  end
if isfield(opts,'alpha') && ~isempty(opts.alpha), alpha = opts.alpha; end
if demag ~= 0
    Nd    = ellipsoid_demagn(alpha);                  % trace-1 demag tensor (sphere -> 1/3 each)
    dm_cc = C.gfac*(4*pi/ion.Vc)*demag*Nd(3,3);       % c-axis demag share (exported as info.Jshape_cc)
    dm_aa = C.gfac*(4*pi/ion.Vc)*demag*Nd(1,1);       % a-axis demag share (folded into info.Jaa0)
else
    dm_cc = 0;  dm_aa = 0;                             % off: byte-identical to the pre-demag code
end

% --- CACHE-FREE: no cacheDir/pkey/key/cacheFile/load/save at all ----------
v = ones(4,1)/2;                 % uniform (all-sublattices-in-phase) ferromagnetic mode
nq = size(qvec,1);
Jnu  = zeros(nq, 4);
Juni = zeros(nq, 1);
lorz = 4*pi/(3*ion.Vc)*C.gfac;   % scalar; broadcasts to ones(4,4)-type Lorentz block
% Build the q-independent lattice geometry ONCE and reuse it for every q below.
% This priming call is itself at q=[0 0 0], so capture its dip0 for the Gamma-
% point info block instead of recomputing it. MF_dipole/exchange are otherwise
% bit-identical whether the geometry is rebuilt or passed in.
[dip0, ~, geomD] = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
[~,       geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
for iq = 1:nq
    q = qvec(iq,:);
    dip = MF_dipole(q, dpRng, ion.a, ion.tau, geomD);       % [3,3,4,4], Å^-3
    ex  = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX); % [3,3,4,4], carries |J12|
    Jcc = -squeeze(C.gfac*dip(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:));
    if invz_is_gamma_equiv(q, ion.tau)
        Jcc = Jcc + lorz;                            % uniform-mode Lorentz cavity (demag-invariant)
    end
    Jcc = (Jcc + Jcc')/2;
    Jnu(iq,:) = sort(real(eig(Jcc))).';
    Juni(iq)  = real(v.'*Jcc*v);                     % uniform FM-mode coupling (physical dispersion)
end
% Gamma-point info block (dip0 from the priming call above), uniform-mode projection:
Jcc0d = -squeeze(C.gfac*dip0(3,3,:,:)) + lorz;
Jaa0d = -squeeze(C.gfac*dip0(1,1,:,:)) + lorz - dm_aa;
Jcc0d = (Jcc0d + Jcc0d')/2;
Jaa0d = (Jaa0d + Jaa0d')/2;
info.Jcc0_dipole = real(v.'*Jcc0d*v);
info.Jaa0_dipole = real(v.'*Jaa0d*v);
info.Jcc0 = info.Jcc0_dipole + 4*ion.J12;
info.Jaa0      = info.Jaa0_dipole + 4*ion.J12;   % transverse J(0), demag-aware (meV)
info.Jshape_cc = 4*dm_cc;                        % strict-uniform observable correction (meV); 0 when demag = 0
info.dpRng = dpRng;
info.geomD = geomD;   % q-independent lattice geometry (MF_dipole/exchange 5-arg reuse form)
info.geomX = geomX;
% --- CACHE-FREE: no save at all, regardless of any opts.cache value -------
end

% =====================================================================
% frozen reference: invz_jq_path.m (cache-free; calls local_jq_modes_ref,
% never the production invz_jq_modes)
% =====================================================================
function P = local_jq_path_ref(ion, qpath, opts)
% Cache-free copy of invz_projected/invz_jq_path.m. Identical directional
% Gamma-limit guard, trust radius, and Hermitization order (Jm Hermitized
% AFTER the broadcast term is added, never Greg alone -- matching the
% production source exactly).
if nargin < 3, opts = struct(); end
dpRng    = 30;  if isfield(opts,'dpRng'),   dpRng    = opts.dpRng;   end
% useCache accepted for call-signature parity only; forwarded into the
% cache-free local_jq_modes_ref call below, which ignores it.
useCache = ~isfield(opts,'cache') || opts.cache;
snapfac  = 2.5; if isfield(opts,'snapfac'), snapfac  = opts.snapfac; end
C  = invz_const();
nq = size(qpath, 1);

[Jnu, info, Juni] = local_jq_modes_ref(ion, qpath, struct('dpRng', dpRng, 'cache', useCache));

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
        % Reuse the q-independent geometry local_jq_modes_ref already built
        % (info.geomD/geomX) instead of re-deriving it from scratch:
        % bit-identical (MF_dipole/exchange 5-arg reuse form), mirroring the
        % production invz_jq_path's own reuse of invz_jq_modes' info.
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
