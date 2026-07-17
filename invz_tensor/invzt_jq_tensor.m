function lat = invzt_jq_tensor(ion, g, opts)
%INVZT_JQ_TENSOR Cached [12,12,nq] full-tensor sublattice x Cartesian coupling (meV).
%   lat = INVZT_JQ_TENSOR(ion, g, opts) returns the LOCKED lattice struct:
%     lat.Jt      [12,12,nq] : Hermitian per page (mu,nu = 1(a),2(b),3(c);
%                 s,s' = 1..4; composite index i = 3*(s-1)+mu). Slices:
%                 cc = Jt(3:3:12,3:3:12,:), ca = Jt(3:3:12,1:3:12,:), etc.
%     lat.qvec, lat.w, lat.conv : echoed from g (see below).
%     lat.JtGamma [12,12]     : the SAME assembly evaluated at q = [0 0 0]
%                 internally (independent of whether qvec contains Gamma).
%     lat.info    : struct(Jcc0, Jaa0, Jcc0_dipole, Jaa0_dipole, lorz, dpRng)
%                 -- the Gamma uniform-mode info block, bit-identical (at
%                 demag = 0) to invz_odd_blocks/invz_jq_modes.
%
%   g is either:
%     - a grid struct from INVZT_QGRID (fields qvec, w, conv), or
%     - a raw [nq,3] qvec array, for single-point/test use: then
%       lat.conv = 'explicit' and lat.w is uniform (1/nq each).
%
%   opts.dpRng (default 30)  : MF_dipole/exchange lattice-sum range (unit cells).
%   opts.cache  (default true): read/write invz_tensor/cache/jqt1_*.mat. Self-
%       verifying (isequal-checks stored pkey + qvec before trusting a hit).
%   opts.parts  (default 'full') : 'full' | 'dipole'. 'dipole' omits the
%       exchange term AND the Lorentz cavity term entirely (full - dipole =
%       exchange + Lorentz exactly); 'dipole' results are NEVER cached.
%
%   ASSEMBLY (Global Constraints, LOCKED; mirrors invz_jq_modes/invz_odd_blocks
%   block-for-block so cc/ca/cb parity is bit-exact, not approximate):
%     J^{mu nu}_{s s'}(q) = -C.gfac*dip(mu,nu,s,s') + sign(ion.J12)*ex(mu,nu,s,s')
%     (parts == 'dipole' drops the second term entirely). At Gamma-equivalent q
%     (invz_is_gamma_equiv), the Lorentz cavity lorz = 4*pi/(3*ion.Vc)*C.gfac is
%     ADDED to the Cartesian-DIAGONAL (mu == nu) entries of ALL 16 sublattice
%     pairs (aa, bb, and cc blocks alike) -- a straight generalisation of the
%     projected cc-only `Jcc = Jcc + lorz` step. Exchange is Cartesian-diagonal
%     by construction (exchange.m), so it never touches mu ~= nu entries.
%     Each Cartesian-diagonal principal sub-block (aa/bb/cc, [4,4] in
%     sublattice space) is Hermitized as (J+J')/2 per q, exactly mirroring
%     invz_odd_blocks' treatment of Vcc; the off-Cartesian-diagonal blocks
%     (ab/ac/bc and mirrors, e.g. ca/cb) are left RAW/un-Hermitized, exactly
%     mirroring Vca/Vcb -- MF_dipole/exchange already enforce the exact
%     (bit-for-bit) sublattice-conjugate relation T(mu,nu,s,s') =
%     conj(T(nu,mu,s',s)) internally (their own `d(:,:,mt,nt) = conj(...)`
%     assignments), so the FULL 12x12 page comes out Hermitian to machine
%     precision without needing a global symmetrization step that would risk
%     perturbing the raw ca/cb values away from bit-exact parity.
%
%   DEMAG GUARD: this layer is intrinsic-only. Errors 'invzt:demag' if
%   ion.demag ~= 0 (checked FIRST, before any lattice work).
%
%   CACHE CONTRACT (Global Constraints; own namespace, never touches
%   invz_projected/cache or odd1_/jq4_ keys):
%     key   : invz_cache_key('jqt1', dpRng, [qvec(:); convcode], pkey)
%     pkey  : [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1]
%     Stores lat, pkey, qvec in ONE file; the loader isequal-verifies BOTH
%     pkey AND qvec before trusting a hit (stale/legacy entries fall through,
%     recompute, and overwrite). 'dipole' parts never read or write this cache.
%
%   See also INVZT_QGRID, INVZ_ODD_BLOCKS, INVZ_JQ_MODES, MF_DIPOLE, EXCHANGE,
%   INVZ_IS_GAMMA_EQUIV, INVZ_CACHE_KEY.
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng') && ~isempty(opts.dpRng), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;
parts = 'full';  if isfield(opts,'parts') && ~isempty(opts.parts), parts = opts.parts; end

% --- demag guard FIRST (tensor layer is intrinsic-only; checked before the
%     parts-string validation below and before any lattice work) -----------
demag = 0;  if isfield(ion,'demag') && ~isempty(ion.demag), demag = ion.demag; end
if demag ~= 0
    error('invzt:demag', ['invzt_jq_tensor is intrinsic-only (ion.demag must be ' ...
        '0; got %g). Demagnetization handling is not part of the tensor layer.'], demag);
end

if ~(ischar(parts) || isstring(parts)) || ~ismember(char(parts), {'full','dipole'})
    error('invzt:parts', 'opts.parts must be ''full'' or ''dipole''; got %s.', ...
        local_str(parts));
end
parts = char(parts);
isFull = strcmp(parts, 'full');

% --- resolve grid input: struct (invzt_qgrid) or raw [nq,3] qvec ------------
if isstruct(g)
    qvec = g.qvec;
    w    = g.w(:);
    convstr = g.conv;
else
    qvec = g;
    nq0  = size(qvec, 1);
    w    = ones(nq0, 1) / nq0;
    convstr = 'explicit';
end
nq = size(qvec, 1);

C = invz_const();
pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 1];
useCacheEff = useCache && isFull;      % dipole-only never cached
cacheDir  = fullfile(fileparts(mfilename('fullpath')), 'cache');
cacheFile = '';
if useCacheEff
    key = invz_cache_key('jqt1', dpRng, [qvec(:); conv_code(convstr)], pkey);
    cacheFile = fullfile(cacheDir, key);
    if exist(cacheFile, 'file')
        S = load(cacheFile);
        if isfield(S,'pkey') && isfield(S,'qvec') && isequal(S.pkey, pkey) && isequal(S.qvec, qvec)
            lat = S.lat;
            return;
        end
        % stale or legacy cache entry: fall through and recompute (file overwritten).
    end
end

% --- geometry priming (once), mirrors invz_jq_modes/invz_odd_blocks ---------
[dip0, ~, geomD] = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
if isFull
    [ex0, geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
else
    ex0 = [];  geomX = [];
end
lorz = 4*pi/(3*ion.Vc)*C.gfac;

Jt = zeros(12, 12, nq);
for iq = 1:nq
    q = qvec(iq,:);
    dip = MF_dipole(q, dpRng, ion.a, ion.tau, geomD);       % [3,3,4,4], Å^-3
    if isFull
        ex = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX); % [3,3,4,4], carries |J12|
    else
        ex = [];
    end
    isGamma = invz_is_gamma_equiv(q, ion.tau);
    Jt(:,:,iq) = assemble_page(dip, ex, isGamma, isFull, C, ion, lorz);
end
JtGamma = assemble_page(dip0, ex0, true, isFull, C, ion, lorz);

% --- Gamma-point info block (bit-identical to invz_odd_blocks/invz_jq_modes
%     at demag = 0; parts-INDEPENDENT -- always the canonical reference
%     values, using the closed-form 4*ion.J12 uniform-exchange shortcut, not
%     an actual q=0 exchange sum) -----------------------------------------
v = ones(4,1)/2;
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

lat.Jt      = Jt;
lat.qvec    = qvec;
lat.w       = w;
lat.conv    = convstr;
lat.JtGamma = JtGamma;
lat.info    = info;

if useCacheEff
    if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
    save(cacheFile, 'lat', 'pkey', 'qvec');
end
end

% ------------------------------- local helpers ------------------------------

function Jb = assemble_page(dip, ex, isGamma, isFull, C, ion, lorz)
% Build one [12,12] page from the [3,3,4,4] dip (and, if isFull, ex) arrays
% for a single q. See the file header for the exact per-block rule.
Jb = zeros(12, 12);
for mu = 1:3
    for nu = 1:3
        Blk = -C.gfac*squeeze(dip(mu,nu,:,:));            % [4,4], always present
        if isFull
            Blk = Blk + sign(ion.J12)*squeeze(ex(mu,nu,:,:));
        end
        if mu == nu && isFull && isGamma
            Blk = Blk + lorz;                              % uniform-mode Lorentz cavity
        end
        if mu == nu
            Blk = (Blk + Blk')/2;                          % Hermitize Cartesian-diagonal block
        end
        rows = mu:3:12;   % s = 1..4 -> row 3*(s-1)+mu
        cols = nu:3:12;
        Jb(rows, cols) = Blk;
    end
end
end

function c = conv_code(convstr)
% Numeric tag folded into the cache key so distinct grid conventions at the
% same q-mesh content never alias (independent copy of invzt_qgrid's own
% conv_code -- matches the codebase's existing per-file hash-helper duplication
% convention; see invz_cache_key.m header).
switch convstr
    case 'explicit',         c = 0;
    case 'halfopen',         c = 1;
    case 'legacy_inclusive', c = 2;
    otherwise
        error('invzt:conv', 'Unknown grid convention %s.', local_str(convstr));
end
end

function s = local_str(x)
if ischar(x) || (isstring(x) && isscalar(x))
    s = char(x);
else
    s = mat2str(x);
end
end
