function [dip, counts, geom] = invz_dipole_ewald(q, a, tau, eopts, geom)
%INVZ_DIPOLE_EWALD  Opt-in Ewald dipolar-tensor primitive (Step-4, Tasks 1-2).
%
%   [dip, counts, geom] = invz_dipole_ewald(q, a, tau, eopts, geom)
%
% Returns the q-dependent dipolar tensor in MF_dipole's quantity/sign/units
% (Angstrom^-3), i.e. dip = -sum'_R T(R+d) e^{-i q.(R+d)}, split into the three
% convergent Ewald parts with the exact k=0 reciprocal term omitted
% (boundary = 'conducting_k0_omitted'; there is NO primitive surface term):
%
%   T_ab(x)  = (3 x_a x_b - |x|^2 delta_ab)/|x|^5,   x = R + d_nm,  d_nm = tau_m - tau_n
%   dip^(r)  = -sum'_{|x|<=r_cut} g_ab(x) e^{-i q_cart.x}
%              g_ab = P(r)(3 x_a x_b - r^2 delta)/r^5 + Q(r) x_a x_b
%              P(r) = erfc(a r) + (2 a r/sqrt(pi)) e^{-a^2 r^2}
%              Q(r) = (4 a^3/sqrt(pi)) e^{-a^2 r^2}/r^2
%   dip^(G)  = +(4 pi/Vc) sum_{k~=0,|k|<=g_cut} (k_a k_b/|k|^2) e^{-|k|^2/4a^2} e^{+i G.d_nm}
%   dip^(self) = -delta_nm delta_ab 4 a^3/(3 sqrt(pi))
%
% with a = eopts.alpha, k = q_cart + G_cart, G_cart = G_hkl * B, B = 2*pi*inv(a)'.
%
% INPUTS
%   q      [nq,3] reduced reciprocal (Miller) coordinates, finite real.
%   a      [3,3]  direct lattice, rows = lattice vectors, finite real nonsingular.
%   tau    [ntau,3] basis positions in fractional coordinates, finite real.
%   eopts  struct with EXACTLY {alpha, r_cut, g_cut, boundary}. alpha/r_cut/g_cut
%          are finite positive scalars; boundary must be 'conducting_k0_omitted'.
%   geom   (optional) reusable q-independent geometry payload from a prior call.
%          Reuse is exact-fingerprint checked (a, tau, alpha, r_cut, g_cut,
%          boundary, q-reduction convention, schema); a mismatch raises
%          invz:ewaldGeomReuse.
%
% OUTPUTS
%   dip    [3,3,ntau,ntau] for nq==1, else [3,3,ntau,ntau,nq]. dip(:,:,n,m) is the
%          block for pair (n,m); dip(:,:,m,n) = conj(dip(:,:,n,m)) is an EMERGENT
%          property (every ordered pair is built independently).
%   counts diagnostic struct (frozen schema, see invzp_ewald_prereg.md sec 1):
%            .real_pair                       [ntau,ntau] retained real vectors
%            .recip_candidates                scalar cached candidate count
%            .recip_used                      [nq,1] per-q retained reciprocal count
%            .preflight.real_cube_bound       [ntau,ntau] conservative integer-mesh bound
%            .preflight.recip_cube_bound      scalar conservative integer-mesh bound
%            .preflight.estimated_peak_bytes  scalar, includes the 1.25 margin
%            .preflight.array_manifest        struct array (name/class/is_complex/
%                                             size/numel/bytes per planned array)
%   geom   q-independent geometry payload (pass back as the 5th argument to skip
%          the rebuild across a q-sweep).
%
% The preflight is a hard pre-ALLOCATION gate: conservative real/reciprocal cube
% bounds and a byte manifest are enforced against the frozen caps BEFORE any
% ndgrid / large allocation (invz:ewaldRealCap / invz:ewaldRecipCap /
% invz:ewaldMemoryCap). See docs/invzp_ewald_prereg.md (FROZEN 2026-07-24).

% ---- frozen constants -------------------------------------------------------
SCHEMA    = 'invz_dipole_ewald/v1';
QCONV     = 'K=floor(q+0.5); qbar=q-K in [-0.5,0.5); gauge=exp(-i*2*pi*K.(tau_m-tau_n))';
CAP_REAL  = 3.0e6;      % max real-space pair vectors per sublattice pair
CAP_RECIP = 3.0e6;      % max reciprocal candidates
CAP_BYTES = 4*2^30;     % estimated peak working bytes (4 GiB)
MARGIN    = 1.25;       % allocator margin on the byte estimate
GUARD     = 4.5;        % accuracy guards: alpha*r_cut>=4.5 AND g_cut/(2 alpha)>=4.5

% ---- strict input validation ------------------------------------------------
if nargin < 4
    error('invz:ewaldArgs', 'invz_dipole_ewald requires q, a, tau, eopts.');
end
local_validate_q(q);
local_validate_a(a);
local_validate_tau(tau);
[alpha, r_cut, g_cut, boundary] = local_validate_eopts(eopts);

if ~strcmp(boundary, 'conducting_k0_omitted')
    error('invz:ewaldBoundary', ...
        'unsupported boundary ''%s''; the primitive only supports ''conducting_k0_omitted''.', ...
        boundary);
end
if ~(alpha*r_cut >= GUARD)
    error('invz:ewaldGuard', ...
        'accuracy guard alpha*r_cut = %.6g < %.3g.', alpha*r_cut, GUARD);
end
if ~(g_cut/(2*alpha) >= GUARD)
    error('invz:ewaldGuard', ...
        'accuracy guard g_cut/(2*alpha) = %.6g < %.3g.', g_cut/(2*alpha), GUARD);
end

nq   = size(q, 1);
ntau = size(tau, 1);
fp_now = local_fingerprint(a, tau, alpha, r_cut, g_cut, boundary, QCONV, SCHEMA);

% ---- preflight step 2: dimension-only lower bound (scalar; no ntau^2 array) --
% Unavoidable pair-indexed metadata + complex output. This early-out runs before
% any ntau^2 displacement/count array is formed.
lower_bytes = (9*ntau^2*nq)*16 ...   % complex output
            + (ntau^2)*8 ...         % counts.real_pair
            + (nq)*8;                % counts.recip_used
if MARGIN*lower_bytes > CAP_BYTES
    error('invz:ewaldMemoryCap', ...
        ['dimension-only lower bound %.4g bytes (x%.2f margin) exceeds cap %.4g ' ...
         'bytes before allocation.'], lower_bytes, MARGIN, CAP_BYTES);
end

% ---- obtain conservative bounds + geometry ---------------------------------
reuse = (nargin >= 5) && ~isempty(geom);
if reuse
    local_check_reuse(geom, fp_now);           % invz:ewaldGeomReuse on mismatch
    real_cube_bound  = geom.real_cube_bound;
    recip_cube_bound = geom.recip_cube_bound;
else
    [real_cube_bound, recip_cube_bound, B, Vc, taucart, sa, sb] = ...
        local_scalar_bounds(a, tau, r_cut, g_cut);
end

% ---- preflight step 4: count caps + byte manifest + memory cap --------------
if any(real_cube_bound(:) > CAP_REAL)
    error('invz:ewaldRealCap', ...
        'real-space candidate cube bound %.4g exceeds cap %.4g before allocation.', ...
        max(real_cube_bound(:)), CAP_REAL);
end
if recip_cube_bound > CAP_RECIP
    error('invz:ewaldRecipCap', ...
        'reciprocal candidate cube bound %.4g exceeds cap %.4g before allocation.', ...
        recip_cube_bound, CAP_RECIP);
end
[manifest, est_bytes] = local_manifest(real_cube_bound, recip_cube_bound, ntau, nq, MARGIN);
if est_bytes > CAP_BYTES
    error('invz:ewaldMemoryCap', ...
        'estimated peak %.4g bytes (incl. %.2f margin) exceeds cap %.4g bytes before allocation.', ...
        est_bytes, MARGIN, CAP_BYTES);
end

% ---- preflight step 5: allocate + build geometry (fresh only) ---------------
if ~reuse
    geom = local_build_geom(a, tau, alpha, r_cut, g_cut, boundary, ...
        B, Vc, taucart, sa, sb, real_cube_bound, recip_cube_bound, fp_now);
end

% ---- q-work + assembly ------------------------------------------------------
[dip, recip_used] = local_assemble(geom, q, alpha, g_cut);

% ---- counts (frozen schema) -------------------------------------------------
counts.real_pair                     = geom.real_pair_count;
counts.recip_candidates              = size(geom.Gcart, 1);
counts.recip_used                    = recip_used;
counts.preflight.real_cube_bound     = real_cube_bound;
counts.preflight.recip_cube_bound    = recip_cube_bound;
counts.preflight.estimated_peak_bytes = est_bytes;
counts.preflight.array_manifest      = manifest;
end

% =====================================================================
% validation
% =====================================================================
function local_validate_q(q)
if ~(isnumeric(q) && isreal(q) && ismatrix(q) && size(q,2)==3 && size(q,1)>=1 ...
        && all(isfinite(q(:))))
    error('invz:ewaldArgs', 'q must be a finite real [nq,3] matrix.');
end
end

function local_validate_a(a)
if ~(isnumeric(a) && isreal(a) && isequal(size(a),[3 3]) && all(isfinite(a(:))))
    error('invz:ewaldArgs', 'a must be a finite real [3,3] matrix.');
end
rc = rcond(a);
if ~isfinite(rc) || rc < 1e-12
    error('invz:ewaldArgs', 'a is singular or numerically singular (rcond = %.3g).', rc);
end
end

function local_validate_tau(tau)
if ~(isnumeric(tau) && isreal(tau) && ismatrix(tau) && size(tau,2)==3 ...
        && size(tau,1)>=1 && all(isfinite(tau(:))))
    error('invz:ewaldArgs', 'tau must be a finite real [ntau,3] matrix.');
end
end

function [alpha, r_cut, g_cut, boundary] = local_validate_eopts(eopts)
if ~(isstruct(eopts) && isscalar(eopts))
    error('invz:ewaldArgs', 'eopts must be a scalar struct.');
end
want = {'alpha','g_cut','r_cut','boundary'};
have = sort(fieldnames(eopts));
if ~isequal(have, sort(want(:)))
    error('invz:ewaldArgs', ...
        'eopts must have EXACTLY the fields {alpha, r_cut, g_cut, boundary}; got {%s}.', ...
        strjoin(fieldnames(eopts).', ', '));
end
alpha = eopts.alpha; r_cut = eopts.r_cut; g_cut = eopts.g_cut; boundary = eopts.boundary;
if ~local_isposscalar(alpha)
    error('invz:ewaldArgs', 'eopts.alpha must be a finite positive scalar.');
end
if ~local_isposscalar(r_cut)
    error('invz:ewaldArgs', 'eopts.r_cut must be a finite positive scalar.');
end
if ~local_isposscalar(g_cut)
    error('invz:ewaldArgs', 'eopts.g_cut must be a finite positive scalar.');
end
if isstring(boundary) && isscalar(boundary)
    boundary = char(boundary);
elseif ~(ischar(boundary) && isrow(boundary))
    error('invz:ewaldArgs', 'eopts.boundary must be a character-row/string scalar.');
end
end

function tf = local_isposscalar(x)
tf = isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0;
end

% =====================================================================
% geometry fingerprint / reuse
% =====================================================================
function fp = local_fingerprint(a, tau, alpha, r_cut, g_cut, boundary, qconv, schema)
fp = struct('a', a, 'tau', tau, 'alpha', alpha, 'r_cut', r_cut, 'g_cut', g_cut, ...
    'boundary', char(boundary), 'qconv', qconv, 'schema', schema);
end

function local_check_reuse(geom, fp_now)
if ~(isstruct(geom) && isscalar(geom) && isfield(geom, 'fingerprint'))
    error('invz:ewaldGeomReuse', 'reused geom is malformed (no fingerprint).');
end
if ~isequaln(geom.fingerprint, fp_now)
    error('invz:ewaldGeomReuse', ...
        ['reused geom fingerprint does not match the current (a, tau, alpha, ' ...
         'r_cut, g_cut, boundary, q-reduction, schema).']);
end
end

% =====================================================================
% conservative scalar bounds (no large allocation)
% =====================================================================
function [real_cube_bound, recip_cube_bound, B, Vc, taucart, sa, sb] = ...
        local_scalar_bounds(a, tau, r_cut, g_cut)
B  = 2*pi*inv(a).';                 %#ok<MINV> reciprocal rows, B = 2*pi*inv(a)'
Vc = abs(det(a));
taucart = tau*a;
ntau = size(tau,1);
sa = min(svd(a));
sb = min(svd(B));
% per-pair conservative real integer-mesh bound
real_cube_bound = zeros(ntau,ntau);
for n = 1:ntau
    for m = 1:ntau
        d = taucart(m,:) - taucart(n,:);
        nmax_r = ceil((r_cut + norm(d))/sa);
        real_cube_bound(n,m) = (2*nmax_r + 1)^3;
    end
end
% conservative reciprocal integer-mesh bound
[c1,c2,c3] = ndgrid([-0.5 0.5],[-0.5 0.5],[-0.5 0.5]);
corners = [c1(:) c2(:) c3(:)];
qmax = max(vecnorm(corners*B, 2, 2));
nmax_G = ceil((g_cut + qmax)/sb);
recip_cube_bound = (2*nmax_G + 1)^3;
end

% =====================================================================
% conservative byte manifest (prereg sec 2, literal: "every planned retained
% and temporary numeric/logical array" -- no size-dependence qualifier).
% Enumerates: real+reciprocal geometry INCLUDING the raw ndgrid triples
% [H,Kk,L] that are live simultaneously with the concatenated hkl/Ghkl_all
% meshes built from them; local_gab's internal work arrays (gab_*); the
% local_boxmin_dist temporaries (boxmin_*), INCLUDING the fixed [3,3] metric
% scratch M and the conservative [3,3] MFF -- previously excluded as
% "O(1)/non-size-dependent", now INCLUDED per the frozen manifest-
% completeness decision: an array counts regardless of how small/fixed its
% shape is, so only genuinely scalar bookkeeping (tol, loop indices, and the
% <=3-element active-set index vectors s/F/Cc used purely for control flow,
% never stored as data) stays off this list; the reciprocal candidate-retain
% mask (recip_keepG); the per-q loop-body temporaries (qwork_*), including
% the complex w.*kK broadcast product; and the output. Rows are summed as if
% simultaneously live, which conservatively over-bounds the true peak
% working set. See test_manifest_names_are_complete's manifest_source_table
% for the authoritative row-by-row source mapping this function must stay in
% sync with.
% =====================================================================
function [manifest, est_bytes] = local_manifest(real_cube_bound, recip_cube_bound, ntau, nq, MARGIN)
RBmax = max(real_cube_bound(:));          % worst single-pair temporary mesh
RBsum = sum(real_cube_bound(:));          % conservative total retained real vectors
GC    = recip_cube_bound;                 % conservative candidate count
rows = struct('name',{},'class',{},'is_complex',{},'size',{},'numel',{},'bytes',{});
% -- real-space geometry: raw ndgrid triple + per-pair build mesh + retained union --
rows(end+1) = local_mkrow('real_ndgrid_raw','double',  false, [RBmax 3]);
rows(end+1) = local_mkrow('real_int_mesh',  'double',  false, [RBmax 3]);
rows(end+1) = local_mkrow('real_cart_mesh', 'double',  false, [RBmax 3]);
rows(end+1) = local_mkrow('real_radius',    'double',  false, [RBmax 1]);
rows(end+1) = local_mkrow('real_mask',      'logical', false, [RBmax 1]);
rows(end+1) = local_mkrow('real_x',         'double',  false, [RBsum 3]);
rows(end+1) = local_mkrow('real_gab',       'double',  false, [RBsum 9]);
% -- local_gab internal work arrays (one call per ordered pair, sized at
%    that pair's retained count <= RBmax) --
rows(end+1) = local_mkrow('gab_P',          'double',  false, [RBmax 1]);
rows(end+1) = local_mkrow('gab_Q',          'double',  false, [RBmax 1]);
rows(end+1) = local_mkrow('gab_r2',         'double',  false, [RBmax 1]);
rows(end+1) = local_mkrow('gab_r5',         'double',  false, [RBmax 1]);
rows(end+1) = local_mkrow('gab_a2r2',       'double',  false, [RBmax 1]);
rows(end+1) = local_mkrow('gab_Tab',        'double',  false, [RBmax 1]);
% -- reciprocal candidate build mesh: raw ndgrid triple + concatenated meshes --
rows(end+1) = local_mkrow('recip_ndgrid_raw','double', false, [GC 3]);
rows(end+1) = local_mkrow('recip_int_mesh', 'double',  false, [GC 3]);
rows(end+1) = local_mkrow('recip_Gcart',    'double',  false, [GC 3]);
rows(end+1) = local_mkrow('recip_dmin',     'double',  false, [GC 1]);
rows(end+1) = local_mkrow('recip_keepG',    'logical', false, [GC 1]);
% -- local_boxmin_dist temporaries (nG = GC candidate rows), including the
%    fixed [3,3] SPD metric M and the conservative [3,3] active-set MFF --
rows(end+1) = local_mkrow('boxmin_M',       'double',  false, [3 3]);
rows(end+1) = local_mkrow('boxmin_MFF',     'double',  false, [3 3]);
rows(end+1) = local_mkrow('boxmin_lo',      'double',  false, [GC 3]);
rows(end+1) = local_mkrow('boxmin_hi',      'double',  false, [GC 3]);
rows(end+1) = local_mkrow('boxmin_v',       'double',  false, [GC 3]);
rows(end+1) = local_mkrow('boxmin_RHS',     'double',  false, [GC 3]);
rows(end+1) = local_mkrow('boxmin_VF',      'double',  false, [GC 3]);
rows(end+1) = local_mkrow('boxmin_feas',    'logical', false, [GC 1]);
rows(end+1) = local_mkrow('boxmin_f',       'double',  false, [GC 1]);
rows(end+1) = local_mkrow('boxmin_fbest',   'double',  false, [GC 1]);
% -- retained reciprocal per-pair phase --
rows(end+1) = local_mkrow('recip_phase',    'double',  true,  [GC*ntau^2 1]);
% -- per-q assembly work: qbar/K vectorized over nq; loop-body temporaries at
%    single-iteration size, matching the existing qwork_k/qwork_kernel modelling --
rows(end+1) = local_mkrow('qwork_qbar',     'double',  false, [nq 3]);
rows(end+1) = local_mkrow('qwork_K',        'double',  false, [nq 3]);
rows(end+1) = local_mkrow('qwork_k',        'double',  false, [GC 3]);
rows(end+1) = local_mkrow('qwork_kernel',   'double',  true,  [GC 1]);
rows(end+1) = local_mkrow('qwork_kk',       'double',  false, [GC 1]);
rows(end+1) = local_mkrow('qwork_keep',     'logical', false, [GC 1]);
rows(end+1) = local_mkrow('qwork_kK',       'double',  false, [GC 3]);
rows(end+1) = local_mkrow('qwork_kk2',      'double',  false, [GC 1]);
rows(end+1) = local_mkrow('qwork_ph',       'double',  true,  [RBmax 1]);
rows(end+1) = local_mkrow('qwork_w',        'double',  true,  [GC 1]);
rows(end+1) = local_mkrow('qwork_w_kK',     'double',  true,  [GC 3]);
% -- output + per-q retained count --
rows(end+1) = local_mkrow('dip_output',     'double',  true,  [3 3 ntau ntau nq]);
rows(end+1) = local_mkrow('recip_used',     'double',  false, [nq 1]);
manifest  = rows(:);
est_bytes = MARGIN * sum([manifest.bytes]);
end

function row = local_mkrow(name, cls, iscx, sz)
n = prod(sz);
switch cls
    case 'double'
        if iscx, w = 16; else, w = 8; end
    case 'logical'
        w = 1;
    otherwise
        w = 8;
end
row = struct('name', name, 'class', cls, 'is_complex', logical(iscx), ...
    'size', sz, 'numel', n, 'bytes', n*w);
end

% =====================================================================
% geometry build (allocations happen only here, after all caps pass)
% =====================================================================
function geom = local_build_geom(a, tau, alpha, r_cut, g_cut, boundary, ...
        B, Vc, taucart, sa, ~, real_cube_bound, recip_cube_bound, fp)
ntau = size(tau,1);

% --- real-space geometry per ordered pair (independent for every (n,m)) ---
realc = cell(ntau,ntau);
real_pair_count = zeros(ntau,ntau);
for n = 1:ntau
    for m = 1:ntau
        d = taucart(m,:) - taucart(n,:);
        nmax_r = ceil((r_cut + norm(d))/sa);
        rng = -nmax_r:nmax_r;
        % [H,Kk,L] (raw ndgrid outputs) are live simultaneously with the
        % concatenated hkl/x built from them below -- accounted in
        % local_manifest as real_ndgrid_raw/real_int_mesh/real_cart_mesh;
        % keep the three in sync.
        [H,Kk,L] = ndgrid(rng, rng, rng);
        hkl = [H(:) Kk(:) L(:)];
        x = hkl*a + d;
        r = vecnorm(x, 2, 2);
        keep = (r <= r_cut);
        if n == m
            keep = keep & ~all(hkl == 0, 2);   % drop the self cell R=0 by cell index
        end
        x = x(keep,:);
        r = r(keep);
        realc{n,m}.x   = x;
        realc{n,m}.gab = local_gab(x, r, alpha);
        real_pair_count(n,m) = size(x,1);
    end
end

% --- reciprocal candidate union (contains G=[0,0,0]) ---------------------
[c1,c2,c3] = ndgrid([-0.5 0.5],[-0.5 0.5],[-0.5 0.5]);
corners = [c1(:) c2(:) c3(:)];
qmax = max(vecnorm(corners*B, 2, 2));
nmax_G = ceil((g_cut + qmax)/min(svd(B)));
rng = -nmax_G:nmax_G;
% [H,Kk,L] (raw ndgrid outputs) are live simultaneously with Ghkl_all/
% Gcart_all built from them below -- accounted in local_manifest as
% recip_ndgrid_raw/recip_int_mesh/recip_Gcart; keep the three in sync.
[H,Kk,L] = ndgrid(rng, rng, rng);
Ghkl_all  = [H(:) Kk(:) L(:)];
Gcart_all = Ghkl_all*B;
% retain every integer G whose min over the half-open cube of |(qbar+G)B| <= g_cut,
% treated conservatively over the CLOSED box (the upper face belongs to the union
% in the half-open limit), with a tiny numerical slack for boundary candidates.
dmin = local_boxmin_dist(Ghkl_all, B);
keepG = dmin <= g_cut*(1 + 1e-12);          % accounted in local_manifest as recip_keepG
Ghkl  = Ghkl_all(keepG,:);
Gcart = Gcart_all(keepG,:);

% --- per-pair reciprocal phase exp(+i G_cart . d_nm) ---------------------
recipc = cell(ntau,ntau);
for n = 1:ntau
    for m = 1:ntau
        d = taucart(m,:) - taucart(n,:);
        recipc{n,m}.phase = exp(1i*(Gcart*d.'));
    end
end

% --- pack ----------------------------------------------------------------
geom = struct();
geom.a = a; geom.tau = tau; geom.taucart = taucart;
geom.alpha = alpha; geom.r_cut = r_cut; geom.g_cut = g_cut; geom.boundary = boundary;
geom.B = B; geom.Vc = Vc; geom.ntau = ntau;
geom.real = realc;
geom.real_pair_count = real_pair_count;
geom.Ghkl = Ghkl; geom.Gcart = Gcart; geom.recip = recipc;
geom.real_cube_bound = real_cube_bound;
geom.recip_cube_bound = recip_cube_bound;
geom.fingerprint = fp;
end

% =====================================================================
% screened Hessian kernel g_ab
% =====================================================================
function gab = local_gab(x, r, alpha)
% x [K,3], r [K,1]. Returns gab [K,3,3]. The internal work arrays
% P,Q,r2,r5,a2r2,Tab (each ~[K,1], K<=RBmax, one call per ordered pair) are
% accounted in local_manifest as gab_P/gab_Q/gab_r2/gab_r5/gab_a2r2/gab_Tab;
% keep the two in sync.
K = size(x,1);
gab = zeros(K,3,3);
if K == 0
    return;
end
r2  = r.^2;
r5  = r2.^2 .* r;
a2r2 = (alpha*r).^2;
P = erfc(alpha*r) + (2*alpha*r/sqrt(pi)).*exp(-a2r2);
Q = (4*alpha^3/sqrt(pi))*exp(-a2r2)./r2;
for aa = 1:3
    for bb = 1:3
        xa = x(:,aa); xb = x(:,bb);
        Tab = (3*xa.*xb - r2*(aa==bb))./r5;
        gab(:,aa,bb) = P.*Tab + Q.*(xa.*xb);
    end
end
end

% =====================================================================
% deterministic box-constrained minimization of |(qbar+G)B| over qbar-box
% =====================================================================
function dmin = local_boxmin_dist(G, B)
% For each integer row G, returns min over qbar in [-0.5,0.5]^3 of |(qbar+G)*B|.
% Metric M = B*B' is SPD; minimize the convex quadratic f(v)=v*M*v' over the box
% v in [G-0.5, G+0.5] by enumerating the 27 free/lower/upper active sets (no
% Optimization Toolbox). The global box-min is the smallest feasible KKT value.
% The size-dependent temporaries below (lo/hi/v/RHS/VF/feas/f/fbest, ~[nG,*]),
% PLUS the fixed [3,3] metric M and the conservative [3,3] active-set MFF, are
% accounted in local_manifest as boxmin_* -- keep the two in sync.
M  = B*B.';
nG = size(G,1);
lo = G - 0.5;
hi = G + 0.5;
tol = 1e-12;
fbest = inf(nG,1);
for s1 = -1:1
    for s2 = -1:1
        for s3 = -1:1
            s = [s1 s2 s3];
            F  = find(s == 0);       % free coordinates
            Cc = find(s ~= 0);       % constrained coordinates (fixed at a bound)
            v = zeros(nG,3);
            for j = Cc
                v(:,j) = G(:,j) + 0.5*s(j);
            end
            feas = true(nG,1);
            if ~isempty(F)
                MFF = M(F,F);          % accounted in local_manifest as boxmin_MFF (conservative [3,3])
                if isempty(Cc)
                    RHS = zeros(numel(F), nG);
                else
                    RHS = -M(F,Cc) * v(:,Cc).';
                end
                VF = (MFF \ RHS).';           % [nG, |F|]
                v(:,F) = VF;
                for jj = 1:numel(F)
                    j = F(jj);
                    feas = feas & (v(:,j) >= lo(:,j) - tol) & (v(:,j) <= hi(:,j) + tol);
                end
            end
            f = sum((v*M).*v, 2);
            f(~feas) = inf;
            fbest = min(fbest, f);
        end
    end
end
dmin = sqrt(max(fbest, 0));
end

% =====================================================================
% assembly (real + reciprocal + self, per q, with extended-zone gauge)
% =====================================================================
function [dip, recip_used] = local_assemble(geom, q, alpha, g_cut)
ntau = geom.ntau;
nq   = size(q,1);
B    = geom.B;
Vc   = geom.Vc;
taufrac = geom.tau;
Gcart   = geom.Gcart;
fourpiVc = 4*pi/Vc;
selfval  = 4*alpha^3/(3*sqrt(pi));
inv4a2   = 1/(4*alpha^2);
gc2      = g_cut^2;

dip = complex(zeros(3,3,ntau,ntau,nq));
recip_used = zeros(nq,1);

% Per-q loop-body temporaries (k/kk/keep/kK/kk2/kernel and, per pair,
% ph/w/w.*kK) are accounted in local_manifest as qwork_* -- keep the two in
% sync.
for i = 1:nq
    qraw = q(i,:);
    K    = floor(qraw + 0.5);        % extended-zone translation
    qbar = qraw - K;                 % reduced into [-0.5,0.5)
    qcart = qbar*B;

    % per-q reciprocal filter (pair-independent): |q_cart+G_cart|<=g_cut, k~=0
    k  = Gcart + qcart;
    kk = sum(k.^2, 2);
    keep = (kk <= gc2) & (kk > 0);
    recip_used(i) = sum(keep);
    kK  = k(keep,:);
    kk2 = kk(keep);
    kernel = fourpiVc * exp(-kk2*inv4a2) ./ kk2;   % [Kk,1]

    for n = 1:ntau
        for m = 1:ntau
            % --- real part: -sum_x g_ab(x) e^{-i q_cart.x},  x = R+d_nm ---
            X   = geom.real{n,m}.x;
            gab = geom.real{n,m}.gab;
            Kr  = size(X,1);
            ph  = exp(-1i*(X*qcart.'));            % [Kr,1]
            gmat = reshape(gab, Kr, 9);
            block = reshape(-(ph.'*gmat), 3, 3);   % [3,3]

            % --- reciprocal part: +(4pi/Vc) sum k_a k_b/|k|^2 ... e^{+iG.d} ---
            w = kernel .* geom.recip{n,m}.phase(keep);   % [Kk,1]
            % w.*kK [Kk,3] complex is accounted in local_manifest as qwork_w_kK.
            block = block + (kK.' * (w .* kK));          % [3,3]

            % --- self part (same-site diagonal only) ---
            if n == m
                block = block - selfval*eye(3);
            end

            % --- extended-zone gauge restore ---
            gph = exp(-1i*2*pi*(K*(taufrac(m,:) - taufrac(n,:)).'));
            dip(:,:,n,m,i) = gph*block;
        end
    end
end
end
