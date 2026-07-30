function [Jnu, info, Juni, Jcc_pages] = invz_jq_modes(ion, qvec, opts)
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
% Optional fourth output Jcc_pages [4 x 4 x nq] returns the exact Hermitian
% longitudinal matrices immediately before the existing eigenvalue call.
% It is nargout-gated: ordinary one- through three-output calls retain the
% existing cache/control path and do not allocate the pages.  A fourth-output
% request intentionally recomputes rather than accepting an eigenvalue-only
% cache hit.
%
% The demagnetizing/shape term is excluded from Jnu/info.Jcc0 (cancels in the
% critical condition per Ronnow, PRB 75, 054426 (2007)); exported separately
% as info.Jshape_cc and inside info.Jaa0. dpRng=30 is the production default
% (grid-convergence checked against R2007 targets to <0.3%).
%
% Dipolar backend (Step-5 Task 2, opt-in; parameters frozen 2026-07-24; operative
% contract retained here):
% opts.dipole = absent | 'bruteforce' (both resolve to the unchanged
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
%                          formula reduces to both frozen Gate-C decompositions).
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
useCache = isfield(opts,'cache') && opts.cache;
wantJccPages = nargout >= 4;
Jcc_pages = [];

% --- Backend dispatch. ---
[backend, eopts] = local_resolve_dipole_backend(opts);
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
if ~wantJccPages && useCache && exist(cacheFile, 'file')
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
if wantJccPages
    Jcc_pages = complex(zeros(4,4,nq));
end
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
        if wantJccPages, Jcc_pages(:,:,iq) = Jcc; end %#ok<AGROW>
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
    % reuse of the already-built geomX).
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
        if wantJccPages, Jcc_pages(:,:,iq) = Jcc; end %#ok<AGROW>
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
% derivation is a higher-layer concern.
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
