function [Jnu, info, Juni] = invz_jq_modes(ion, qvec, opts)
%INVZ_JQ_MODES Eigenvalue branches of the 4x4 cc sublattice coupling matrix (meV).
% J_cc(q)_{rs} = -gfac*dip_cc_{rs}(q) [+ Lorentz at q≡0] + sign(J12)*|J12|*ex_cc_{rs}(q).
% Convention: ferromagnetic-positive; criticality when J(0)*chi0 = 1+Sigma(0).
%
% Third output Juni [nq x 1] is the UNIFORM ferromagnetic-mode projection
% v'*Jcc(q)*v (v=[1 1 1 1]/2), the physical single-mode dispersion under the
% mean-field/RPA approximation (Juni([0 0 0]) == info.Jcc0; matches
% MF_RPA_Yikai.m's J_avg exactly). Use Juni, NOT max(eig), to trace a
% dispersion along a q-path: away from Gamma the uniform mode stops being an
% eigenvector, so sorted branches cross and max(eig) picks the wrong branch
% (mirrors the (1,0,0)->(2,0,0) dispersion about h=1.5). See invz_jq_path.
%
% The demagnetizing/shape term is excluded from Jnu/info.Jcc0 (cancels in the
% critical condition per Ronnow, PRB 75, 054426 (2007)); exported separately
% as info.Jshape_cc and inside info.Jaa0. dpRng=30 is the production default
% (grid-convergence checked against R2007 targets to <0.3%).
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng'), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;
C = invz_const();
% Lorentz cavity (+4pi/(3Vc)) is always added at the uniform mode (mandatory
% split term, not a toggle). Demagnetization (ion.demag/ion.alpha, default
% off) cancels from the critical condition per R2007, so Jnu/info.Jcc0/Tc(B=0)
% are demag-invariant; it is exported instead as info.Jshape_cc (applied
% downstream via chi_meas = chi/(1+Jshape_cc*chi) in invz_chi_realaxis) and
% folded into demag-aware info.Jaa0, through which Bc(T) vs applied field can
% still shift.
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
cacheDir = fullfile(fileparts(mfilename('fullpath')), 'cache');
pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; demag; alpha; 4];   % trailing 4 = cache schema v4 (adds info.geomD/geomX)
key = sprintf('jq4_%d_%s_%s.mat', dpRng, hash_vec(qvec(:)), hash_vec(pkey));
cacheFile = fullfile(cacheDir, key);
if useCache && exist(cacheFile, 'file')
    S = load(cacheFile);
    if isfield(S,'pkey') && isfield(S,'qvec') && isfield(S,'Juni') && isequal(S.pkey, pkey) && isequal(S.qvec, qvec)
        Jnu = S.Jnu;  info = S.info;  Juni = S.Juni;  return;
    end
    % stale or legacy cache entry: fall through and recompute (file will be overwritten)
end
v = ones(4,1)/2;                 % uniform (all-sublattices-in-phase) ferromagnetic mode
nq = size(qvec,1);
Jnu  = zeros(nq, 4);
Juni = zeros(nq, 1);
lorz = 4*pi/(3*ion.Vc)*C.gfac;   % scalar; broadcasts to ones(4,4)-type Lorentz block (see header)
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
% Γ-point info block (dip0 from the priming call above), uniform-mode projection:
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
info.geomD = geomD;   % q-independent lattice geometry (MF_dipole/exchange 5-arg reuse form),
info.geomX = geomX;   % exposed so callers (e.g. invz_jq_path's Gamma-limit Greg) can rebuild a
                       % q=0 dip/ex matrix without re-deriving the geometry (bit-identical either way).
if useCache
    if ~exist(cacheDir,'dir'), mkdir(cacheDir); end
    save(cacheFile, 'Jnu', 'info', 'Juni', 'pkey', 'qvec');
end
end

function h = hash_vec(v)
h = sprintf('%dv_%08x', numel(v), ...
    typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
