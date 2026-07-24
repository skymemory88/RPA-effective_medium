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
%
% ODD extension (T1.3, opt-in): opts.odd = false (default) | struct('Xp', [2x2])
% with Xp the real symmetric transverse susceptibility block from invz_chiperp
% (meV^-1). With opts.odd a struct, the cc matrices gain the ODD-mediated
% coupling deltaJ(q) (invz_odd_deltaJ, E1 + self-site subtraction E4) built
% from the cached geometric blocks (invz_odd_blocks, odd1_ namespace), and
% info.Jcc0 carries the explicit uniform shift -d (E5) -- the single line that
% feeds the DS2023 mean-field mechanism into the MF, the RPA denominator and
% the 1/z criticality. Diagnostics land in info.odd (see jq_modes_odd below).
% With opts.odd false/absent the default path below runs untouched (bitwise
% regression-gated); the ODD path never reads or writes the jq4_/jq5_ cache.
% opts.odd remains brute-force-only: an active ODD request combined with the
% Ewald dipolar backend (below) is rejected before either diversion runs.
%
% Dipolar backend (Step-5 Task 2, opt-in; docs/invzp_ewald_prereg.md FROZEN,
% docs/invzp_ewald_design.md Sec.2.2/2.3/4.2, docs/invzp_ewald_integration_map.md
% Sec.6.3): opts.dipole = absent | 'bruteforce' (both resolve to the unchanged
% brute-force MF_dipole path; identical cache identity) | 'ewald' (opt-in
% invz_dipole_ewald primitive; opts.ewald must then be a scalar struct with
% EXACTLY {alpha,r_cut,g_cut,boundary} -- this function does not synthesize
% frozen defaults, see the higher-level drivers for that). The production
% default remains bruteforce; no default flip happens here.
%
% Both backends additively export (appended AFTER every legacy field):
%   info.dipole         = struct('backend',...,'ewald',...,'q_reduction',...,
%                          'primitive_schema',...) -- full backend provenance;
%                          bruteforce reports a canonical empty 'ewald' value
%                          and a documented legacy q-convention string.
%   info.Jpath_base_cc  = [4x4], NOT pre-Hermitized, the backend-agnostic
%                          q-path reconstruction base invz_jq_path consumes:
%                            bruteforce = -gfac*dip_sphere_cc(0) + Jex_cc(0)
%                            ewald      = -gfac*dip_reg_cc(0) + Jex_cc(0) - lorz*ones(4)
%   info.Jgamma_cc      = info.Jpath_base_cc + lorz*ones(4), the exact-Gamma
%                          backend-agnostic production matrix (this single
%                          formula reduces to both frozen Gate-C decompositions
%                          -- see docs/invzp_ewald_prereg.md Sec.5 Gate-C).
% Under Ewald, the regularized dipolar tensor already contains the isotropic
% Lorentz term (design Sec.4.2: "Ewald adds 0 at Gamma"), so the per-q loop and
% the Gamma-point Jcc0_dipole/Jaa0_dipole priming block add NO extra +lorz;
% info.dpRng is NaN (dpRng does not affect the Ewald calculation or its cache
% identity) and info.geomD is never populated (info.geomX, the UNCHANGED
% exchange geometry, remains present under both backends).
%
% Cache: schema 'invz_jq_modes/v5', filenames 'jq5_<backend>_...'. A hit is
% accepted only after an exact isequaln validation of a structured cacheMeta
% payload (qvec, lattice, basis, Vc, J12, gfac, demag, top-level aspect ratio,
% backend, exact Ewald controls or a canonical empty brute-force value, brute
% dpRng or the Ewald NaN sentinel, the BZ cacheContext if the caller supplies
% one via opts.cacheContext [Task 4] or a canonical direct-call context,
% schema version, and the required info-field/output-shape contract); any
% missing/legacy/malformed/mismatched payload is a miss and is recomputed.
% Absent backend and explicit 'bruteforce' resolve to the identical backend
% string and therefore share one canonical cache identity. This replaces the
% v4 'jq4_' scheme (not extended in place); the ODD path's separate 'odd1_'
% cache is untouched.
if nargin < 3, opts = struct(); end
dpRng = 30;  if isfield(opts,'dpRng'), dpRng = opts.dpRng; end
useCache = ~isfield(opts,'cache') || opts.cache;

% --- Backend dispatch (Step-5 Task 2): resolved/validated BEFORE the ODD
% diversion below, per docs/invzp_ewald_design.md Sec.1.1/2.2. ---
[backend, eopts] = local_resolve_dipole_backend(opts);

activeOdd = isfield(opts, 'odd') && ~isempty(opts.odd) && ~isequal(opts.odd, false);
if strcmp(backend, 'ewald') && activeOdd
    error('invz:jqModesOddEwald', ['the ODD extension (opts.odd) is not supported together with ' ...
        'the Ewald dipolar backend (opts.dipole=''ewald''); opts.odd must be false/absent when ' ...
        'requesting Ewald. The ODD path remains brute-force-only (docs/invzp_ewald_design.md Sec.1.1).']);
end

% --- ODD diversion (T1.3): strictly additive and opt-in. Everything below this
% block is the pre-ODD code path, byte-untouched (regression test
% test_jq_modes_odd_off_bitwise gates isequaln on all three outputs).
if activeOdd
    [Jnu, info, Juni] = jq_modes_odd(ion, qvec, opts, dpRng, useCache);
    return
end
C = invz_const();
% Lorentz cavity (+4pi/(3Vc)) is always added at the uniform mode (mandatory
% split term, not a toggle) for the BRUTE-FORCE backend. Demagnetization
% (ion.demag/ion.alpha, default off) cancels from the critical condition per
% R2007, so Jnu/info.Jcc0/Tc(B=0) are demag-invariant; it is exported instead
% as info.Jshape_cc (applied downstream via chi_meas = chi/(1+Jshape_cc*chi) in
% invz_chi_realaxis) and folded into demag-aware info.Jaa0, through which
% Bc(T) vs applied field can still shift. This block is backend-independent.
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
cacheContext = local_resolve_cache_context(opts);
reqInfoFields = local_required_info_fields(backend);
cacheMeta = local_build_cache_meta(backend, eopts, dpRng, qvec, ion, C, demag, alpha, cacheContext, reqInfoFields);
pkeyNum = local_pkey_numeric(backend, eopts, dpRng, ion, C, demag, alpha);
dpTag = 'NaN';  if strcmp(backend,'bruteforce'), dpTag = sprintf('%d', dpRng); end
key = sprintf('jq5_%s_%s_%s_%s.mat', backend, dpTag, hash_vec(qvec(:)), hash_vec(pkeyNum));
cacheFile = fullfile(cacheDir, key);
if useCache && exist(cacheFile, 'file')
    S = load(cacheFile);
    if local_cache_hit_valid(S, cacheMeta, reqInfoFields)
        Jnu = S.Jnu;  info = S.info;  Juni = S.Juni;  return;
    end
    % missing/legacy/malformed/mismatched cache entry: fall through and recompute (file will be overwritten)
end
v = ones(4,1)/2;                 % uniform (all-sublattices-in-phase) ferromagnetic mode
nq = size(qvec,1);
Jnu  = zeros(nq, 4);
Juni = zeros(nq, 1);
lorz = 4*pi/(3*ion.Vc)*C.gfac;   % scalar; broadcasts to ones(4,4)-type Lorentz block (see header)

if strcmp(backend, 'bruteforce')
    % ================= BRUTE-FORCE PATH (operation order byte-preserved) =================
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
    % ---- Additive Gamma metadata (Step-5 Task 2), appended AFTER all legacy fields. ----
    % Recovers the Gamma exchange tensor the priming call above discards (bit-identical
    % reuse of the already-built geomX -- docs/invzp_ewald_integration_map.md Sec.6.1).
    ex0 = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau, geomX);
    info.dipole = struct('backend', 'bruteforce', 'ewald', local_empty_ewald(), ...
        'q_reduction', ['bruteforce: q used directly as MF_dipole/exchange Miller indices ' ...
                         '(q*geom.b); no canonical q-domain reduction applied'], ...
        'primitive_schema', 'MF_dipole+exchange (legacy, unversioned)');
    info.Jpath_base_cc = -C.gfac*squeeze(dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:));
    info.Jgamma_cc     = info.Jpath_base_cc + lorz*ones(4);
else
    % ================= EWALD PATH (additive; opt-in) =================
    % Build/reuse ONE invz_dipole_ewald geometry across every q below (same
    % priming-call pattern as the brute-force branch). exchange is UNCHANGED
    % and used identically to the brute-force branch (design Sec.9: "exchange
    % is out of scope"). No +lorz is added anywhere here: the regularized
    % Ewald tensor already contains the isotropic Lorentz term at Gamma
    % (design Sec.4.2; prereg Sec.5 Gate-C).
    [dip0, ~, geomE] = invz_dipole_ewald([0 0 0], ion.a, ion.tau, eopts);
    [~,       geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
    for iq = 1:nq
        q = qvec(iq,:);
        dip = invz_dipole_ewald(q, ion.a, ion.tau, eopts, geomE);
        ex  = exchange(q, abs(ion.J12), ion.a, ion.tau, geomX);
        Jcc = -squeeze(C.gfac*dip(3,3,:,:)) + sign(ion.J12)*squeeze(ex(3,3,:,:));
        Jcc = (Jcc + Jcc')/2;
        Jnu(iq,:) = sort(real(eig(Jcc))).';
        Juni(iq)  = real(v.'*Jcc*v);
    end
    Jcc0d = -squeeze(C.gfac*dip0(3,3,:,:));            % NO +lorz: already inside dip_reg
    Jaa0d = -squeeze(C.gfac*dip0(1,1,:,:)) - dm_aa;    % NO +lorz: already inside dip_reg
    Jcc0d = (Jcc0d + Jcc0d')/2;
    Jaa0d = (Jaa0d + Jaa0d')/2;
    info.Jcc0_dipole = real(v.'*Jcc0d*v);
    info.Jaa0_dipole = real(v.'*Jaa0d*v);
    info.Jcc0 = info.Jcc0_dipole + 4*ion.J12;
    info.Jaa0      = info.Jaa0_dipole + 4*ion.J12;     % transverse J(0), demag-aware (meV); same semantics as bruteforce
    info.Jshape_cc = 4*dm_cc;                          % same caller-level demag semantics as bruteforce
    info.dpRng = NaN;                  % Ewald: dpRng does not affect the calculation or cache identity
    info.geomX = geomX;                % UNCHANGED exchange geometry; info.geomD is intentionally NEVER set
    ex0 = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau, geomX);
    info.dipole = struct('backend', 'ewald', 'ewald', eopts, ...
        'q_reduction', geomE.fingerprint.qconv, 'primitive_schema', geomE.fingerprint.schema);
    info.Jpath_base_cc = -C.gfac*squeeze(dip0(3,3,:,:)) + sign(ion.J12)*squeeze(ex0(3,3,:,:)) - lorz*ones(4);
    info.Jgamma_cc     = info.Jpath_base_cc + lorz*ones(4);   % = -gfac*dip_reg_cc(0)+Jex_cc(0): adds 0 extra
end

if useCache
    if ~exist(cacheDir,'dir'), mkdir(cacheDir); end
    save(cacheFile, 'Jnu', 'info', 'Juni', 'cacheMeta');
end
end

function h = hash_vec(v)
h = sprintf('%dv_%08x', numel(v), ...
    typecast(single(sum(v.*(1:numel(v))')), 'uint32'));
end

% =====================================================================
% Step-5 Task 2: backend dispatch validation, cache metadata, cache-hit
% validation. Every raised identifier is stable and invz:jqModes*-namespaced.
% =====================================================================
function [backend, eopts] = local_resolve_dipole_backend(opts)
% Resolves opts.dipole (absent or 'bruteforce' -> 'bruteforce' | 'ewald' -> 'ewald')
% and strictly validates opts.ewald against the resolved backend. jq_modes does
% NOT synthesize frozen Ewald defaults -- opts.ewald must already be a
% complete {alpha,r_cut,g_cut,boundary} struct for the Ewald backend; default
% derivation is a higher-layer concern (docs/invzp_ewald_design.md Sec.2.2).
if ~isfield(opts, 'dipole') || isempty(opts.dipole)
    backend = 'bruteforce';
else
    raw = opts.dipole;
    if isstring(raw) && isscalar(raw)
        raw = char(raw);
    end
    if ~(ischar(raw) && isrow(raw))
        error('invz:jqModesBackend', ...
            ['opts.dipole must be a scalar string/char naming a backend ' ...
             '(''bruteforce''|''ewald''); got class %s.'], class(opts.dipole));
    end
    if ~(strcmp(raw,'bruteforce') || strcmp(raw,'ewald'))
        error('invz:jqModesBackend', ...
            'unknown opts.dipole backend ''%s''; supported backends are ''bruteforce'' and ''ewald''.', raw);
    end
    backend = raw;
end

hasEwaldOpts = isfield(opts, 'ewald') && ~isempty(opts.ewald);
if hasEwaldOpts && ~strcmp(backend, 'ewald')
    error('invz:jqModesEwaldOptsUnexpected', ...
        ['opts.ewald was supplied but the resolved opts.dipole backend is ''%s'', not ''ewald''; ' ...
         'Ewald controls are only accepted with an explicit opts.dipole=''ewald'' request.'], backend);
end

eopts = [];
if strcmp(backend, 'ewald')
    if ~hasEwaldOpts || ~isstruct(opts.ewald) || ~isscalar(opts.ewald)
        error('invz:jqModesEwaldOptsFields', ...
            ['opts.dipole=''ewald'' requires a complete scalar struct opts.ewald with EXACTLY the ' ...
             'fields {alpha, r_cut, g_cut, boundary}; jq_modes does not synthesize frozen defaults.']);
    end
    want = sort({'alpha','r_cut','g_cut','boundary'});
    have = sort(reshape(fieldnames(opts.ewald), 1, []));
    if ~isequal(have, want)
        error('invz:jqModesEwaldOptsFields', ...
            'opts.ewald must have EXACTLY the fields {alpha, r_cut, g_cut, boundary}; got {%s}.', ...
            strjoin(reshape(fieldnames(opts.ewald), 1, []), ', '));
    end
    b = opts.ewald.boundary;
    if isstring(b) && isscalar(b), b = char(b); end
    if ~(ischar(b) && isrow(b) && strcmp(b, 'conducting_k0_omitted'))
        error('invz:jqModesEwaldBoundary', ...
            'opts.ewald.boundary must be ''conducting_k0_omitted''; got %s.', local_describe_value(opts.ewald.boundary));
    end
    eopts = opts.ewald;
end
end

function s = local_describe_value(x)
if (ischar(x) && isrow(x)) || (isstring(x) && isscalar(x))
    s = ['''' char(x) ''''];
else
    s = class(x);
end
end

function ewaldEmpty = local_empty_ewald()
% Canonical empty Ewald value: same field names as a real eopts struct (for a
% predictable, backend-independent info.dipole.ewald shape), sentinel-empty
% values. Used for BOTH info.dipole.ewald and the bruteforce cacheMeta.ewald.
ewaldEmpty = struct('alpha', [], 'r_cut', [], 'g_cut', [], 'boundary', '');
end

function ctx = local_resolve_cache_context(opts)
% The BZ layer (Task 4) will supply/own a private, validated opts.cacheContext
% (route, grid dimensions, convention, offset, Gamma policy, q digest); stored
% and isequaln-validated as-is here. Direct/anchor calls that supply none get
% a canonical direct-call sentinel so their cache identity stays well-defined
% and distinct from any future BZ-layer context 'kind'.
if isfield(opts, 'cacheContext') && ~isempty(opts.cacheContext)
    ctx = opts.cacheContext;
else
    ctx = struct('kind', 'direct_call');
end
end

function fn = local_required_info_fields(backend)
% Backend-specific fixed info field-name contract (part of the v5 cacheMeta
% "required output field names/shapes"); info.geomD exists ONLY for bruteforce.
core = {'Jcc0_dipole','Jaa0_dipole','Jcc0','Jaa0','Jshape_cc','dpRng','geomX', ...
        'dipole','Jpath_base_cc','Jgamma_cc'};
if strcmp(backend, 'bruteforce')
    fn = sort(reshape([core, {'geomD'}], [], 1));
else
    fn = sort(reshape(core, [], 1));
end
end

function cacheMeta = local_build_cache_meta(backend, eopts, dpRng, qvec, ion, C, demag, alpha, cacheContext, reqInfoFields)
% Structured v5 cache-identity payload (global cache contract: exact qvec,
% lattice, basis, Vc, J12, gfac, demag, top-level aspect ratio, backend, exact
% Ewald controls or a canonical empty brute-force value, brute dpRng or the
% Ewald NaN sentinel, BZ cacheContext or a canonical direct-call context,
% schema version, required output field names/shapes). Accepted on a filename
% hit ONLY after an exact isequaln match (see local_cache_hit_valid).
cacheMeta = struct();
cacheMeta.schema  = 'invz_jq_modes/v5';
cacheMeta.qvec    = qvec;
cacheMeta.lattice = ion.a;
cacheMeta.basis   = ion.tau;
cacheMeta.Vc      = ion.Vc;
cacheMeta.J12     = ion.J12;
cacheMeta.gfac    = C.gfac;
cacheMeta.demag   = demag;
cacheMeta.aspect  = alpha;             % top-level ellipsoid aspect ratio (NOT eopts.alpha)
cacheMeta.backend = backend;
if strcmp(backend, 'ewald')
    cacheMeta.ewald = eopts;
    cacheMeta.dpRng = NaN;
else
    cacheMeta.ewald = local_empty_ewald();
    cacheMeta.dpRng = dpRng;
end
cacheMeta.cacheContext  = cacheContext;
cacheMeta.reqInfoFields = reqInfoFields;
cacheMeta.JnuCols  = 4;
cacheMeta.JuniCols = 1;
end

function pkeyNum = local_pkey_numeric(backend, eopts, dpRng, ion, C, demag, alpha)
% Compact numeric fingerprint for the cache FILENAME digest only (convenience/
% collision-reduction, NOT the safety mechanism -- see local_cache_hit_valid's
% isequaln check on the full cacheMeta, which is what actually gates a hit).
pkeyNum = [ion.a(:); ion.tau(:); ion.Vc; ion.J12; C.gfac; demag; alpha; 5];   % trailing 5 = schema v5
if strcmp(backend, 'ewald')
    pkeyNum = [pkeyNum; eopts.alpha; eopts.r_cut; eopts.g_cut];
else
    pkeyNum = [pkeyNum; dpRng];
end
end

function ok = local_cache_hit_valid(S, cacheMeta, reqInfoFields)
% A filename hit is accepted only after this exact structural + isequaln
% validation. Missing, legacy (pre-v5), malformed, or mismatched payloads are
% ALL treated as misses (recomputed), never trusted.
ok = false;
if ~isfield(S,'cacheMeta') || ~isfield(S,'Jnu') || ~isfield(S,'info') || ~isfield(S,'Juni')
    return   % missing/legacy (e.g. jq4_-style) payload
end
if ~isequaln(S.cacheMeta, cacheMeta)
    return   % any mismatch anywhere in the structured payload
end
nq = size(cacheMeta.qvec, 1);
if ~isequal(size(S.Jnu), [nq 4]) || ~isequal(size(S.Juni), [nq 1])
    return   % malformed shapes
end
if ~isequal(sort(reshape(fieldnames(S.info), [], 1)), reqInfoFields)
    return   % malformed/legacy info field set
end
ok = true;
end

function [Jnu, info, Juni] = jq_modes_odd(ion, qvec, opts, dpRng, useCache)
%JQ_MODES_ODD opts.odd branch of invz_jq_modes (T1.3; E1/E4/E5).
% Rebuilds the cc channel as M(q) = Vcc(q) + deltaJ(q) from the geometric ODD
% blocks (invz_odd_blocks; odd1_ cache when useCache) and the supplied
% transverse susceptibility Xp = opts.odd.Xp, then eigendecomposes per q
% exactly as the default path does (same sort(real(eig(.))) and uniform-mode
% Juni = v'*M*v). NEVER reads or writes the jq4_/jq5_ cache: the ODD Jnu depend
% on Xp, which is not part of that key; the heavy geometric blocks are cached
% under odd1_ by invz_odd_blocks itself. Always brute-force (Ewald+ODD is
% rejected before this function is ever reached -- see the dispatch above).
%
% info mirrors the default path field-for-field:
%   Jcc0            = infoB.Jcc0 - d  <- THE single line that carries the
%                     MF-level DS2023 mechanism into the mean field, the RPA
%                     denominator and the 1/z criticality. Bookkeeping (plan
%                     T1.3): the grid matrices carry the POST-subtraction
%                     deltaJ -- whose diagonal already contains -d (E4) --
%                     and Jcc0 carries the explicit -d (E5); there is NO
%                     other q = 0 handling anywhere.
%   Jcc0_dipole, Jaa0_dipole, Jaa0, dpRng: from invz_odd_blocks' Gamma info
%                     block, bit-identical to the default path at demag = 0
%                     (T1.1 parity test). info.lorz rides along (extra,
%                     harmless).
%   Jshape_cc       = 0: the ODD layer is intrinsic-only (demag == 0 enforced
%                     below), and the default path's 4*dm_cc is 0 at demag = 0.
%   geomD/geomX     : invz_odd_blocks does not export the primed geometry, so
%                     it is rebuilt here once (~0.16 s at dpRng 30; MF_dipole's
%                     documented reuse form is bit-identical either way) to
%                     keep the info contract whole for callers (invz_jq_path
%                     reads info.geomD for its Gamma-limit rebuild).
%   odd             = struct(d, dJ_mean_diag [4,1] pre-subtraction diagonal
%                     means, dJ_max, uniform_residual, Xp). uniform_residual =
%                     max |PRE-subtraction deltaJ| element on the smallest-|q|
%                     (r.l.u.) non-Gamma shell of qvec (E1 recomputed exactly
%                     on that shell) -- report-only; on-axis c* shells sit at
%                     machine zero, near-a* shells carry the linear-in-q ODD
%                     residual (ODD-LOG P0.3/T1.3), and tilted rays keep a
%                     direction-dependent macroscopic term, so never gate
%                     small-q decay off-axis.
odd = opts.odd;
if ~isstruct(odd) || ~isfield(odd, 'Xp')
    error('invz:oddXp', ['opts.odd must be false (default) or a struct with a ' ...
        'field Xp = [2x2] real symmetric chi_perp block (invz_chiperp); got %s.'], ...
        class(odd));
end
Xp = odd.Xp;
if ~isnumeric(Xp) || ~isequal(size(Xp), [2 2]) || ~isreal(Xp) ...
        || ~all(isfinite(Xp(:))) || ~issymmetric(Xp)
    error('invz:oddXp', ['opts.odd.Xp must be a [2 2] real symmetric finite ' ...
        'matrix (meV^-1, invz_chiperp output).']);
end
demag = 0;
if isfield(ion,'demag')  && ~isempty(ion.demag),  demag = ion.demag;  end
if isfield(opts,'demag') && ~isempty(opts.demag), demag = opts.demag; end
if demag ~= 0
    error('invz:oddDemag', ['the ODD path of invz_jq_modes is intrinsic-only ' ...
        '(demag must be 0; got %g). Demag handling stays on the default path/' ...
        'invz_chi_realaxis (R2007: shape term cancels from the critical ' ...
        'condition).'], demag);
end
[Vca, Vcb, Vcc, infoB] = invz_odd_blocks(ion, qvec, ...
    struct('dpRng', dpRng, 'cache', useCache));
[dJ, d, dinfo] = invz_odd_deltaJ(Vca, Vcb, Xp);
nq   = size(qvec, 1);
[Jnu, Juni] = invz_odd_modes(Vcc, dJ);   % shared values-only kernel + uniform projection
info = infoB;
info.Jcc0 = infoB.Jcc0 - d;      % E5 explicit uniform shift (see header block)
info.Jshape_cc = 0;              % demag == 0 enforced above (intrinsic-only)
% Primed q-independent geometry, mirroring the default path's info contract:
[~, ~, info.geomD] = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);
[~,    info.geomX] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
% uniform_residual: E1 recomputed exactly on the smallest-|q| non-Gamma shell.
qn = sqrt(sum(qvec.^2, 2));
for iq = 1:nq
    if invz_is_gamma_equiv(qvec(iq,:), ion.tau), qn(iq) = Inf; end
end
qmin = min(qn);
if isfinite(qmin)
    ur = 0;
    for iq = reshape(find(qn <= qmin*(1 + 1e-9)), 1, [])
        A = Vca(:,:,iq);  B = Vcb(:,:,iq);
        P = Xp(1,1)*(A*A') + Xp(1,2)*(A*B') + Xp(2,1)*(B*A') + Xp(2,2)*(B*B');
        ur = max(ur, max(abs(P(:))));
    end
else
    ur = NaN;                    % all-Gamma grid: no finite-q shell to measure
end
info.odd = struct('d', d, 'dJ_mean_diag', dinfo.d_per_sublattice, ...
    'dJ_max', dinfo.dJ_max, 'uniform_residual', ur, 'Xp', Xp);
end
