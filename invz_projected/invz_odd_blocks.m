function [Vca, Vcb, Vcc, info] = invz_odd_blocks(ion, qvec, opts)
%INVZ_ODD_BLOCKS Geometric off-diagonal dipolar (ODD) sublattice blocks (meV).
%   [Vca, Vcb, Vcc, info] = INVZ_ODD_BLOCKS(ion, qvec, opts) returns, per q in
%   qvec (rows, r.l.u.), the three 4x4 sublattice coupling matrices that feed the
%   ODD-mediated coupling deltaJ^{cc} (E1/E4/E5, Tasks 3-4). These are the pure
%   geometric lattice sums: no chi_perp, no deltaJ, no field.
%
%   OUTPUTS
%     Vca [4,4,nq] complex : J^{ca}(q)_{ss'} = -C.gfac*dip(3,1,s,s')  (DIPOLE-ONLY)
%     Vcb [4,4,nq] complex : J^{cb}(q)_{ss'} = -C.gfac*dip(3,2,s,s')  (DIPOLE-ONLY)
%     Vcc [4,4,nq] Hermitian: the SAME assembly invz_jq_modes eigendecomposes,
%                  -C.gfac*dip(3,3,:,:) + sign(ion.J12)*ex(3,3,:,:), plus the
%                  Lorentz cavity `lorz` at Gamma-equivalent q, Hermitized per q.
%     info : struct with fields dpRng, Jcc0, Jaa0, Jcc0_dipole, Jaa0_dipole, lorz
%            (the same Gamma info block invz_jq_modes reports; verified equal by
%            the parity test test_vcc_parity_with_jq_modes).
%
%   CONVENTIONS & UNITS (must match the codebase; ODD plan §1.3)
%     - Energies meV. `dip` from MF_dipole is [3,3,4,4] in Å^-3, pre-gfac;
%       Cartesian index mu,nu = 1..3 = (a,b,c) = (x,y,z); sublattice s,s' = 1..4.
%       `C.gfac` carries mu0/4pi*(gL*muB)^2, so J^{calpha} = -gfac*dip is in meV.
%     - J^{ca}(q) is NOT Hermitian and the individual ODD blocks are NEVER
%       Hermitized. The pair identities are J^{ca}(-q) = conj(J^{ca}(q)) and
%       J^{ac}(q) = J^{ca}(q)' (real-space couplings are real). Vcc, being a
%       cc-diagonal Cartesian block, IS Hermitized per q (as in invz_jq_modes).
%     - Exchange is isotropic (Cartesian-diagonal) => ex(3,alpha)=0 for alpha~=3,
%       so it contributes NOTHING to Vca/Vcb; Lorentz/demag are Cartesian-diagonal
%       => no ODD shape terms. Vca/Vcb are therefore pure dipole lattice sums and
%       receive NO Lorentz. Vcc is the full ordering-channel assembly.
%     - info.Jcc0 = info.Jcc0_dipole + 4*ion.J12; info.Jcc0_dipole = v'*Jcc0d*v
%       with v = [1 1 1 1]/2, Jcc0d = -gfac*dip0(3,3) + lorz (from the priming
%       call at q=[0 0 0]). Likewise info.Jaa0 from dip0(1,1). Matches
%       invz_jq_modes at demag = 0 (this layer is intrinsic-only, see below).
%
%   INPUTS / OPTS
%     opts.dpRng (default 30) : MF_dipole/exchange lattice-sum range (unit cells).
%     opts.cache (default true): read/write the geometric blocks to the
%       invz_projected/cache directory under the `odd1_` namespace.
%
%   CACHE CONTRACT (Global Constraints; own namespace, never touches jq4_)
%     Key   : odd1_<dpRng>_<hash(qvec)>_<hash(pkey)>.mat  (schema tag odd1)
%     pkey  : [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]  (trailing 1)
%     Stores: Vca, Vcb, Vcc, info, pkey, qvec in ONE file; the loader
%             isequal-verifies BOTH pkey AND qvec before trusting a hit (stale or
%             legacy entries fall through, recompute, and overwrite).
%
%   DEMAG GUARD
%     The ODD layer is intrinsic-only (demag handling stays in invz_jq_modes,
%     where the shape term cancels from the critical condition per R2007). If
%     ion.demag ~= 0 this errors with id 'invz:oddDemag'.
%
%   Mirrors invz_jq_modes' geometry priming, per-q loop and Gamma handling so the
%   cc channel is bit-for-bit the matrix invz_jq_modes eigendecomposes.
%   See also INVZ_JQ_MODES, MF_DIPOLE, EXCHANGE, INVZ_IS_GAMMA_EQUIV.
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng'), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;
C = invz_const();

% Intrinsic-only layer: demag must be off. The shape term is Cartesian-diagonal
% and cancels from the critical condition (R2007); demag stays in invz_jq_modes.
demag = 0;  if isfield(ion,'demag') && ~isempty(ion.demag), demag = ion.demag; end
if demag ~= 0
    error('invz:oddDemag', ...
        ['invz_odd_blocks is intrinsic-only (ion.demag must be 0; got %g). ' ...
         'Demagnetization handling stays in invz_jq_modes.'], demag);
end

cacheDir = fullfile(fileparts(mfilename('fullpath')), 'cache');
pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1];   % trailing 1 = odd1 schema
key  = sprintf('odd1_%d_%s_%s.mat', dpRng, hash_vec(qvec(:)), hash_vec(pkey));
cacheFile = fullfile(cacheDir, key);
if useCache && exist(cacheFile, 'file')
    S = load(cacheFile);
    if isfield(S,'pkey') && isfield(S,'qvec') && isequal(S.pkey, pkey) && isequal(S.qvec, qvec)
        Vca = S.Vca;  Vcb = S.Vcb;  Vcc = S.Vcc;  info = S.info;  return;
    end
    % stale or legacy cache entry: fall through and recompute (file overwritten).
end

nq   = size(qvec,1);
v    = ones(4,1)/2;                 % uniform (all-sublattices-in-phase) FM mode
lorz = 4*pi/(3*ion.Vc)*C.gfac;      % scalar Lorentz cavity, broadcast at Gamma
Vca  = zeros(4,4,nq);
Vcb  = zeros(4,4,nq);
Vcc  = zeros(4,4,nq);
% Build the q-independent lattice geometry ONCE and reuse it for every q. The
% priming call is at q=[0 0 0]; capture its dip0 for the Gamma info block instead
% of recomputing (MF_dipole/exchange are bit-identical whether geom is rebuilt).
[dip0, ~, geomD] = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
[~,       geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
for iq = 1:nq
    q   = qvec(iq,:);
    dip = MF_dipole(q, dpRng, ion.a, ion.tau, geomD);       % [3,3,4,4], Å^-3
    ex  = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX); % [3,3,4,4], carries |J12|
    Vca(:,:,iq) = -C.gfac*squeeze(dip(3,1,:,:));            % J^{ca} (NOT Hermitized)
    Vcb(:,:,iq) = -C.gfac*squeeze(dip(3,2,:,:));            % J^{cb} (NOT Hermitized)
    Jcc = -squeeze(C.gfac*dip(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:));
    if invz_is_gamma_equiv(q, ion.tau)
        Jcc = Jcc + lorz;                                   % uniform-mode Lorentz cavity
    end
    Vcc(:,:,iq) = (Jcc + Jcc')/2;                           % cc Hermitian per q
end
% Gamma-point info block (dip0 from the priming call), uniform-mode projection.
% demag = 0 here (guarded above) => bit-identical to invz_jq_modes' Jcc0/Jaa0.
Jcc0d = -squeeze(C.gfac*dip0(3,3,:,:)) + lorz;
Jaa0d = -squeeze(C.gfac*dip0(1,1,:,:)) + lorz;
Jcc0d = (Jcc0d + Jcc0d')/2;
Jaa0d = (Jaa0d + Jaa0d')/2;
info.dpRng       = dpRng;
info.Jcc0_dipole = real(v.'*Jcc0d*v);
info.Jaa0_dipole = real(v.'*Jaa0d*v);
info.Jcc0        = info.Jcc0_dipole + 4*ion.J12;
info.Jaa0        = info.Jaa0_dipole + 4*ion.J12;
info.lorz        = lorz;
if useCache
    if ~exist(cacheDir,'dir'), mkdir(cacheDir); end
    save(cacheFile, 'Vca', 'Vcb', 'Vcc', 'info', 'pkey', 'qvec');
end
end

function h = hash_vec(v)
h = sprintf('%dv_%08x', numel(v), ...
    typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
