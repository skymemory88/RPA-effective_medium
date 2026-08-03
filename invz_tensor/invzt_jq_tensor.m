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
%   opts.dpRng (default 30)  : brute-force MF_dipole lattice-sum range (unit cells).
%   opts.dipole (default 'bruteforce'): 'bruteforce' | 'ewald'.
%   opts.ewald                : with dipole='ewald', an exact scalar struct
%       {alpha,r_cut,g_cut,boundary}; boundary='conducting_k0_omitted'.
%   opts.cache  (default true): read/write invz_tensor/cache/jqt2_*.mat. Hits
%       require exact structured-metadata equality, not only a weak filename hash.
%   opts.parts  (default 'full') : 'full' | 'dipole'. 'dipole' omits the
%       exchange term. For brute force it also omits the wrapper Lorentz cavity
%       (full-dipole = exchange+Lorentz); Ewald has no wrapper cavity term
%       (full-dipole = exchange). 'dipole' results are NEVER cached.
%
%   ASSEMBLY (Global Constraints, LOCKED; mirrors invz_jq_modes/invz_odd_blocks
%   block-for-block so cc/ca/cb parity is bit-exact, not approximate):
%     J^{mu nu}_{s s'}(q) = -C.gfac*dip(mu,nu,s,s') + sign(ion.J12)*ex(mu,nu,s,s')
%     (parts == 'dipole' drops the second term entirely). On the brute backend,
%     at Gamma-equivalent q (invz_is_gamma_equiv), the Lorentz cavity
%     lorz = 4*pi/(3*ion.Vc)*C.gfac is
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
%     perturbing the raw ca/cb values away from bit-exact parity. On the Ewald
%     backend invz_dipole_ewald supplies the regularized Gamma tensor under its
%     declared conducting/k=0-omitted convention, so NO additional +lorz is made.
%
%   DEMAG GUARD: this layer is intrinsic-only. Errors 'invzt:demag' if
%   ion.demag ~= 0 (checked FIRST, before any lattice work).
%
%   CACHE CONTRACT: schema invzt_jq_tensor/v2 in the jqt2_<backend> namespace.
%   The stored cacheMeta exactly records qvec/weights/convention, lattice/basis,
%   couplings, backend, exact Ewald controls (or a canonical empty value), and
%   brute dpRng/Ewald NaN sentinel. Legacy jqt1 files are never consumed.
%
%   See also INVZT_QGRID, INVZ_ODD_BLOCKS, INVZ_JQ_MODES, INVZ_DIPOLE_EWALD,
%   MF_DIPOLE, EXCHANGE, INVZ_IS_GAMMA_EQUIV, INVZ_CACHE_KEY.
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng') && ~isempty(opts.dpRng), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;
parts = 'full';  if isfield(opts,'parts') && ~isempty(opts.parts), parts = opts.parts; end
[backend, eopts] = resolve_dipole_backend(opts);

% --- demag guard FIRST (tensor layer is intrinsic-only; checked before the
%     parts-string validation below and before any lattice work) -----------
demag = 0;  if isfield(ion,'demag') && ~isempty(ion.demag), demag = ion.demag; end
if demag ~= 0
    error('invzt:demag', ['invzt_jq_tensor is intrinsic-only (ion.demag must be ' ...
        '0; got %g). Demagnetization handling is not part of the tensor layer.'], demag);
end

if ~(ischar(parts) || isstring(parts)) || ~ismember(char(parts), {'full','dipole'})
    error('invzt:parts', 'opts.parts must be ''full'' or ''dipole''; got %s.', ...
        invzt_str(parts));
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
cacheMeta = build_cache_meta(backend, eopts, dpRng, qvec, w, convstr, ion, C, parts);
if strcmp(backend, 'ewald')
    cacheDp = 0;
    pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 2; ...
            eopts.alpha; eopts.r_cut; eopts.g_cut];
else
    cacheDp = dpRng;
    pkey = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; 2; dpRng];
end
useCacheEff = useCache && isFull;      % dipole-only never cached
cacheDir  = fullfile(fileparts(mfilename('fullpath')), 'cache');
cacheFile = '';
if useCacheEff
    key = invz_cache_key(['jqt2_' backend], cacheDp, ...
        [qvec(:); w(:); conv_code(convstr)], pkey);
    cacheFile = fullfile(cacheDir, key);
    if exist(cacheFile, 'file')
        S = load(cacheFile);
        if cache_hit_valid(S, cacheMeta)
            lat = S.lat;
            return;
        end
        % stale or legacy cache entry: fall through and recompute (file overwritten).
    end
end

% --- geometry priming (once), mirrors invz_jq_modes/invz_odd_blocks ---------
geomD = [];
geomE = [];
if strcmp(backend, 'ewald')
    [dip0, ~, geomE] = invz_dipole_ewald([0 0 0], ion.a, ion.tau, eopts);
else
    [dip0, ~, geomD] = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
end
if isFull
    [ex0, geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
else
    ex0 = [];  geomX = [];
end
lorz = 4*pi/(3*ion.Vc)*C.gfac;

Jt = zeros(12, 12, nq);
for iq = 1:nq
    q = qvec(iq,:);
    if strcmp(backend, 'ewald')
        dip = invz_dipole_ewald(q, ion.a, ion.tau, eopts, geomE);
    else
        dip = MF_dipole(q, dpRng, ion.a, ion.tau, geomD);   % [3,3,4,4], Ang^-3
    end
    if isFull
        ex = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX); % [3,3,4,4], carries |J12|
    else
        ex = [];
    end
    isGamma = invz_is_gamma_equiv(q, ion.tau);
    Jt(:,:,iq) = assemble_page(dip, ex, isGamma, isFull, C, ion, lorz, backend);
end
JtGamma = assemble_page(dip0, ex0, true, isFull, C, ion, lorz, backend);

% --- Gamma-point info block (bit-identical to invz_odd_blocks/invz_jq_modes
%     at demag = 0; parts-INDEPENDENT -- always the canonical reference
%     values, using the closed-form 4*ion.J12 uniform-exchange shortcut, not
%     an actual q=0 exchange sum) -----------------------------------------
v = ones(4,1)/2;
gammaLorz = 0;
if strcmp(backend, 'bruteforce'), gammaLorz = lorz; end
Jcc0d = -squeeze(C.gfac*dip0(3,3,:,:)) + gammaLorz;
Jaa0d = -squeeze(C.gfac*dip0(1,1,:,:)) + gammaLorz;
Jcc0d = (Jcc0d + Jcc0d')/2;
Jaa0d = (Jaa0d + Jaa0d')/2;
if strcmp(backend, 'ewald'), info.dpRng = NaN; else, info.dpRng = dpRng; end
info.Jcc0_dipole = real(v.'*Jcc0d*v);
info.Jaa0_dipole = real(v.'*Jaa0d*v);
info.Jcc0        = info.Jcc0_dipole + 4*ion.J12;
info.Jaa0        = info.Jaa0_dipole + 4*ion.J12;
info.lorz        = lorz;
if strcmp(backend, 'ewald')
    qReduction = geomE.fingerprint.qconv;
    primitiveSchema = geomE.fingerprint.schema;
else
    qReduction = ['bruteforce: q used directly as MF_dipole/exchange Miller ' ...
        'indices; no canonical q-domain reduction'];
    primitiveSchema = 'MF_dipole+exchange (legacy, unversioned)';
end
info.dipole = struct('backend', backend, 'ewald', eopts, ...
    'q_reduction', qReduction, 'primitive_schema', primitiveSchema);

lat.Jt      = Jt;
lat.qvec    = qvec;
lat.w       = w;
lat.conv    = convstr;
lat.JtGamma = JtGamma;
lat.info    = info;

if useCacheEff
    if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
    save(cacheFile, 'lat', 'cacheMeta');
end
end

% ------------------------------- local helpers ------------------------------

function [backend, eopts] = resolve_dipole_backend(opts)
% Resolve the opt-in backend without inventing Ewald controls.  Keeping this
% contract identical to invz_jq_modes prevents the projected and tensor paths
% from silently using different regularizations for the same user options.
if ~isfield(opts, 'dipole') || isempty(opts.dipole)
    backend = 'bruteforce';
else
    raw = opts.dipole;
    if isstring(raw) && isscalar(raw), raw = char(raw); end
    if ~(ischar(raw) && isrow(raw))
        error('invzt:jqTensorBackend', ...
            'opts.dipole must be a scalar string/char (''bruteforce''|''ewald'').');
    end
    if ~ismember(raw, {'bruteforce','ewald'})
        error('invzt:jqTensorBackend', ...
            'Unknown opts.dipole backend ''%s''.', raw);
    end
    backend = raw;
end

hasEwald = isfield(opts, 'ewald') && ~isempty(opts.ewald);
if hasEwald && ~strcmp(backend, 'ewald')
    error('invzt:jqTensorEwaldOptsUnexpected', ...
        'opts.ewald is accepted only with opts.dipole = ''ewald''.');
end

if strcmp(backend, 'ewald')
    if ~hasEwald || ~isstruct(opts.ewald) || ~isscalar(opts.ewald)
        error('invzt:jqTensorEwaldOptsFields', ...
            ['opts.dipole=''ewald'' requires a scalar opts.ewald with EXACTLY ' ...
             '{alpha, r_cut, g_cut, boundary}.']);
    end
    want = sort({'alpha','r_cut','g_cut','boundary'});
    have = sort(reshape(fieldnames(opts.ewald), 1, []));
    if ~isequal(have, want)
        error('invzt:jqTensorEwaldOptsFields', ...
            'opts.ewald must have EXACTLY {alpha, r_cut, g_cut, boundary}.');
    end
    b = opts.ewald.boundary;
    if isstring(b) && isscalar(b), b = char(b); end
    if ~(ischar(b) && isrow(b) && strcmp(b, 'conducting_k0_omitted'))
        error('invzt:jqTensorEwaldBoundary', ...
            'opts.ewald.boundary must be ''conducting_k0_omitted''.');
    end
    eopts = opts.ewald;
    eopts.boundary = b;       % canonicalize string scalar to character row
else
    eopts = empty_ewald();
end
end

function eopts = empty_ewald()
eopts = struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', '');
end

function meta = build_cache_meta(backend, eopts, dpRng, qvec, w, convstr, ion, C, parts)
% Full structured identity is the safety check; the filename digest is only a
% lookup convenience and may collide without compromising correctness.
meta = struct();
meta.schema = 'invzt_jq_tensor/v2';
meta.backend = backend;
meta.ewald = eopts;
if strcmp(backend, 'ewald'), meta.dpRng = NaN; else, meta.dpRng = dpRng; end
meta.qvec = qvec;
meta.w = w;
meta.conv = convstr;
meta.lattice = ion.a;
meta.basis = ion.tau;
meta.Vc = ion.Vc;
meta.J12 = ion.J12;
meta.gfac = C.gfac;
meta.parts = parts;
meta.output_shape = [12 12 size(qvec,1)];
end

function ok = cache_hit_valid(S, meta)
% Treat legacy, incomplete, or corrupted payloads as misses.  Exact metadata
% equality protects identity; these checks protect the promised output shape
% and the backend provenance consumed by real-axis reconstruction.
ok = false;
if ~isfield(S, 'cacheMeta') || ~isfield(S, 'lat') || ~isequaln(S.cacheMeta, meta)
    return;
end
lat = S.lat;
want = {'Jt','qvec','w','conv','JtGamma','info'};
if ~(isstruct(lat) && isscalar(lat) && all(isfield(lat, want)))
    return;
end
if size(lat.Jt,1) ~= 12 || size(lat.Jt,2) ~= 12 ...
        || size(lat.Jt,3) ~= meta.output_shape(3) ...
        || ~isequal(size(lat.JtGamma), [12 12]) ...
        || ~isequaln(lat.qvec, meta.qvec) || ~isequaln(lat.w, meta.w) ...
        || ~isequaln(lat.conv, meta.conv)
    return;
end
if ~(isstruct(lat.info) && isscalar(lat.info) && isfield(lat.info, 'dipole'))
    return;
end
d = lat.info.dipole;
if ~(isstruct(d) && isscalar(d) && isfield(d, 'backend') && isfield(d, 'ewald') ...
        && isequaln(d.backend, meta.backend) && isequaln(d.ewald, meta.ewald))
    return;
end
ok = true;
end

function Jb = assemble_page(dip, ex, isGamma, isFull, C, ion, lorz, backend)
% Build one [12,12] page from the [3,3,4,4] dip (and, if isFull, ex) arrays
% for a single q. See the file header for the exact per-block rule.
Jb = zeros(12, 12);
for mu = 1:3
    for nu = 1:3
        Blk = -C.gfac*squeeze(dip(mu,nu,:,:));            % [4,4], always present
        if isFull
            Blk = Blk + sign(ion.J12)*squeeze(ex(mu,nu,:,:));
        end
        if mu == nu && isFull && isGamma && strcmp(backend, 'bruteforce')
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
        error('invzt:conv', 'Unknown grid convention %s.', invzt_str(convstr));
end
end
