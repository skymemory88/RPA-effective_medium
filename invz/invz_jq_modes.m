function [Jnu, info, Juni] = invz_jq_modes(ion, qvec, opts)
%INVZ_JQ_MODES Eigenvalue branches of the 4x4 cc sublattice coupling matrix (meV).
% J_cc(q)_{rs} = -gfac*dip_cc_{rs}(q) [+ Lorentz at q≡0] + sign(J12)*|J12|*ex_cc_{rs}(q)
% Convention: ferromagnetic-positive; criticality when J(0)*chi0 = 1+Sigma(0).
%
% Third output Juni [nq x 1] is the UNIFORM ferromagnetic-mode projection of the cc
% coupling at each q,  Juni(q) = v'*Jcc(q)*v  with v = [1 1 1 1]/2 (all sublattices in
% phase). This is the q-resolved generalization of info.Jcc0 (Juni([0 0 0]) == Jcc0),
% and equals MF_RPA_Yikai.m's J_avg = sum(sum(J_int,4),3)/unitN exactly. It is the
% PHYSICAL single-mode dispersion under the "all four moments equivalent" mean-field/RPA
% approximation both codes make (mode pole set by v'*Jcc*v). Use Juni -- NOT max(eig) of
% the 4x4 block -- to trace a dispersion along a q-path: away from Gamma the uniform mode
% stops being an eigenvector, so the sorted branches cross and max(eig) selects the wrong
% branch (it mirrors the (1,0,0)->(2,0,0) dispersion about h=1.5). See invz_jq_path.
%
% ---------------------------------------------------------------------------
% Step-4 physics checkpoint (Rønnow et al., PRB 75, 054426 (2007), eq. 4):
% target J_D*D_aa(0) = 3.912 ueV, J_D*D_cc(0) = 6.821 ueV. Candidates computed
% at dpRng=30, q=[0 0 0] (gfac=0.08388 meV*Ang^3, Vc=287.9 Ang^3):
%
%   Quantity                                            cc (ueV)   aa (ueV)
%   max-eig(+gfac*dip), no Lorentz [brief Step-3 draft]    5.206      0.971
%   max-eig(+gfac*dip) + lorz (draft lorz, no x4 fix)      6.426*     2.192   (*coincidental near-miss)
%   uniform-mode(+gfac*dip), no Lorentz (row sum)         -1.943      0.971
%   uniform-mode(-gfac*dip), no Lorentz                    1.943     -0.971
%   uniform-mode(-gfac*dip) + lorz  <-- RESOLVED           6.824      3.910
%   max-eig(-gfac*dip) + lorz                              6.824      6.113   (aa WRONG: max != uniform)
%
% Resolution (both targets matched simultaneously, <0.1% error, within the
% 3% test tolerance):
%   1. SIGN: MF_dipole.m returns dip(n,m,i,j) = -sum_r exp(-iq.r)(3 r_n r_m/r^5
%      - delta_nm/r^3) (note its own built-in leading minus). Our convention
%      needs the OPPOSITE overall sign from the brief's literal Step-3 draft,
%      i.e. use -gfac*dip (not +gfac*dip), so the c-axis dipole sum comes out
%      ferromagnetic-positive at Gamma.
%   2. LORENTZ MAGNITUDE: lorz = 4*pi/(3*Vc)*gfac is added as a MATLAB scalar
%      broadcast to the full 4x4 sublattice matrix (Jcc = Jcc + lorz), NOT as
%      lorz*eye(4). Because the uniform (ferromagnetic) mode is the constant
%      vector v=[1 1 1 1], a scalar broadcast is equivalent to adding
%      lorz*ones(4,4), whose only nonzero eigenvalue is 4*lorz on the uniform
%      mode (ones(4,4)*v = 4*v) -- reproducing the expected Lorentz share
%      J_D*(4*pi/3) = 4.882 ueV with J_D = 4*gfac/Vc, WITHOUT needing an
%      explicit extra factor of 4 in the `lorz` formula itself. (This matches
%      how MF_RPA_Yikai.m broadcasts eins=repmat(eye(3),1,1,4,4), whose (3,3)
%      slice is 1 for every sublattice pair (i,j) -- i.e. also a ones(4,4)-type
%      correction, not eye(4).)
%   3. EXTRACTION: info.Jcc0_dipole / info.Jaa0_dipole are the eigenvalue of
%      the UNIFORM (all-sublattices-equal, k=0 ferromagnetic) eigenvector
%      v=[1 1 1 1]/2, computed as the quadratic form v.'*J*v (equivalently the
%      row sum, since lattice translational symmetry guarantees v is an exact
%      eigenvector -- verified numerically to residual ~1e-18). This is NOT
%      always max(eig(...)): for the cc axis the uniform mode happens to BE
%      the largest eigenvalue (consistent with LiHoF4 being a true q=0
%      ferromagnet in the Ising/c channel), but for the aa axis the uniform
%      mode is the SECOND-SMALLEST eigenvalue -- taking max(eig(...)) for aa
%      would give 6.113 ueV (56% off target). The brief's Step-3 draft code
%      (`max(real(eig(...)))`) is therefore corrected to an explicit uniform-
%      mode projection for both info.Jcc0_dipole and info.Jaa0_dipole.
%   The demagnetizing/sample-shape term is deliberately EXCLUDED from Jnu and info.Jcc0
%   (it cancels in the critical condition per R2007); it is exported separately as
%   info.Jshape_cc (strict-uniform observable correction) and inside info.Jaa0
%   (transverse channel). See the knob block below.
%
% The per-q branch matrix Jcc(q) (used for Jnu, all 4 modes) uses the SAME
% sign-flip and Lorentz-broadcast convention as the Gamma-point diagnostic,
% so the full dispersion is internally consistent with info.Jcc0.
%
% dpRng = 30 chosen as production default: grid-convergence spot check (Task 5
% Step 5) gave Jcc0_dipole = 6.8394 (R=10), 6.8230 (R=20), 6.8244 (R=30),
% 6.8218 ueV (R=40) -- all within ~0.3% of each other (target 6.821 ueV);
% Jaa0_dipole = 3.9030, 3.9111, 3.9104, 3.9117 ueV over the same R's (target
% 3.912 ueV). The sharp real-space sphere cutoff in MF_dipole gives small
% (<0.3%) non-monotonic fluctuations rather than smooth decay, but R=30 is
% already well inside that noise band, so it is used as the default.
% ---------------------------------------------------------------------------
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng'), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;
C = invz_const();
% Sample-shape terms: the Lorentz cavity +4pi/(3Vc) is ALWAYS added at the uniform mode (it is
% the mandatory cavity term of the dipole-sum split, not a physical toggle). The demagnetization
% correction (ion.demag/ion.alpha, default off) NEVER touches the ordering channel: per R2007 it
% cancels from the critical condition (ordering occurs at q -> 0+, not strict q = 0), so Jnu,
% info.Jcc0, and Tc(B=0) are demag-invariant. The shape is exported instead as (a) info.Jshape_cc,
% the strict-uniform OBSERVABLE correction applied downstream in invz_chi_realaxis
% (chi_meas = chi/(1 + Jshape_cc*chi)), and (b) demag-aware info.Jaa0, the transverse
% (non-critical) mean-field coupling -- through which Bc(T) vs APPLIED field can still shift
% (internal-vs-applied transverse field relation).
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
pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; demag; alpha; 3];   % trailing 3 = cache schema v3 (adds Juni)
key = sprintf('jq3_%d_%s_%s.mat', dpRng, hash_vec(qvec(:)), hash_vec(pkey));
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
for iq = 1:nq
    q = qvec(iq,:);
    dip = MF_dipole(q, dpRng, ion.a, ion.tau);       % [3,3,4,4], Å^-3
    ex  = exchange(q, abs(ion.J12), ion.a, ion.tau); % [3,3,4,4], carries |J12|
    Jcc = -squeeze(C.gfac*dip(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:));
    if is_gamma_equiv(q, ion.tau)
        Jcc = Jcc + lorz;                            % uniform-mode Lorentz cavity (demag-invariant)
    end
    Jcc = (Jcc + Jcc')/2;
    Jnu(iq,:) = sort(real(eig(Jcc))).';
    Juni(iq)  = real(v.'*Jcc*v);                     % uniform FM-mode coupling (physical dispersion)
end
% Γ-point info block (always computed at q=[0 0 0]), uniform-mode projection:
dip0 = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
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
if useCache
    if ~exist(cacheDir,'dir'), mkdir(cacheDir); end
    save(cacheFile, 'Jnu', 'info', 'Juni', 'pkey', 'qvec');
end
end

function tf = is_gamma_equiv(q, tau)
tf = abs(real(sum(exp(2i*pi*(tau*q.'))))/size(tau,1) - 1) < 1e-9;
end

function h = hash_vec(v)
h = sprintf('%dv_%08x', numel(v), ...
    typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end
