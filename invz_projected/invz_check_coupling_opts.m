function H = invz_check_coupling_opts()
%INVZ_CHECK_COUPLING_OPTS Shared precomputed-coupling backend/grid provenance
% and conflict validation, used by BOTH invz_spectra_map.m and
% invz_spectra_qpath.m (Ewald Step-5 Task 7). Extracted (task-7 review, dedup fix)
% from a ~116-line block that was previously duplicated byte-for-byte in both
% drivers -- same shared-helper precedent as invz_check_solve_opts.m /
% invz_common/invz_check_transverse_mf.m, structured as a struct-of-handles
% because several of the five functions are called directly from the drivers'
% own control flow, not only from one another.
%
%   H = INVZ_CHECK_COUPLING_OPTS() returns a struct of function handles:
%     H.validate_dipole_opts(opts)                           -> [backend, eopts]
%     H.has_complete_dipole_provenance(info)                  -> tf
%     H.has_complete_grid_provenance(info)                    -> tf
%     H.check_backend_provenance(info, backendReq, eoptsReq)
%     H.check_grid_provenance(info, opts, grid)
%
%   H.validate_dipole_opts intentionally replicates invz_jq_modes.m's own
%   local_resolve_dipole_backend checks AND error identifiers, so a malformed
%   opts.dipole/opts.ewald request is rejected identically whether couplings
%   are computed for real or a precomputed opts.Jnu/opts.info bypasses
%   invz_jq_modes entirely.
%
%   Pure refactor: every identifier, message string, and comparison rule below
%   is unchanged from the pre-dedup driver copies (verified byte-identical
%   before extraction).
H = struct();
H.validate_dipole_opts           = @local_validate_dipole_opts;
H.has_complete_dipole_provenance = @local_has_complete_dipole_provenance;
H.has_complete_grid_provenance   = @local_has_complete_grid_provenance;
H.check_backend_provenance       = @local_check_backend_provenance;
H.check_grid_provenance          = @local_check_grid_provenance;
end

function [backend, eopts] = local_validate_dipole_opts(opts)
if ~isfield(opts, 'dipole') || isempty(opts.dipole)
    backend = 'bruteforce';
else
    raw = opts.dipole;
    if isstring(raw) && isscalar(raw), raw = char(raw); end
    if ~(ischar(raw) && isrow(raw))
        error('invz:jqModesBackend', ...
            ['opts.dipole must be a scalar string/char naming a backend ' ...
             '(''bruteforce''|''ewald''); got class %s.'], class(opts.dipole));
    end
    if ~(strcmp(raw, 'bruteforce') || strcmp(raw, 'ewald'))
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
    want = sort({'alpha', 'r_cut', 'g_cut', 'boundary'});
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
            'opts.ewald.boundary must be ''conducting_k0_omitted''; got %s.', class(opts.ewald.boundary));
    end
    eopts = opts.ewald;
end
end

function tf = local_has_complete_dipole_provenance(info)
%LOCAL_HAS_COMPLETE_DIPOLE_PROVENANCE True iff info.dipole is a complete backend-provenance
% struct (invz_jq_modes' info.dipole contract: backend/ewald/q_reduction/primitive_schema).
tf = isfield(info, 'dipole') && isstruct(info.dipole) && isscalar(info.dipole) && ...
     all(isfield(info.dipole, {'backend', 'ewald', 'q_reduction', 'primitive_schema'}));
end

function tf = local_has_complete_grid_provenance(info)
%LOCAL_HAS_COMPLETE_GRID_PROVENANCE True iff info.grid is a complete BZ grid-policy provenance
% struct (invz_bz_couplings' info.grid contract).
tf = isfield(info, 'grid') && isstruct(info.grid) && isscalar(info.grid) && ...
     all(isfield(info.grid, {'schema', 'convention', 'offset', 'gammaPolicy', 'requested', ...
                             'nominal', 'retained', 'n_gamma', 'qhash'}));
end

function local_check_backend_provenance(info, backendReq, eoptsReq)
%LOCAL_CHECK_BACKEND_PROVENANCE An explicit opts.dipole/opts.ewald request was made alongside
% precomputed opts.Jnu/opts.info: require complete info.dipole and compare EXACTLY.
if ~local_has_complete_dipole_provenance(info)
    error('invz:spectraBackendProvenanceMissing', ['an explicit opts.dipole/opts.ewald request ' ...
        'was given together with precomputed opts.Jnu/opts.info, but opts.info lacks complete ' ...
        'info.dipole provenance (a struct with backend/ewald/q_reduction/primitive_schema); cannot ' ...
        'verify the precomputed couplings actually used the requested backend.']);
end
conflict = ~strcmp(info.dipole.backend, backendReq);
if ~conflict && strcmp(backendReq, 'ewald')
    conflict = ~isequaln(info.dipole.ewald, eoptsReq);
end
if conflict
    error('invz:spectraBackendConflict', ['the requested opts.dipole/opts.ewald backend/controls ' ...
        'conflict with the precomputed opts.info.dipole provenance (requested ''%s'', precomputed ' ...
        'info carries ''%s'').'], backendReq, info.dipole.backend);
end
end

function local_check_grid_provenance(info, opts, grid)
%LOCAL_CHECK_GRID_PROVENANCE An explicit gridConvention/gridOffset/gammaPolicy request was made
% alongside precomputed opts.Jnu/opts.info: require complete info.grid, resolve any omitted
% request member to invz_bz_couplings' own defaults, and compare EXACTLY (convention, offset,
% gammaPolicy, and the requested grid dimensions).
if ~local_has_complete_grid_provenance(info)
    error('invz:spectraGridProvenanceMissing', ['an explicit gridConvention/gridOffset/gammaPolicy ' ...
        'request was given together with precomputed opts.Jnu/opts.info, but opts.info lacks ' ...
        'complete info.grid provenance; cannot verify the precomputed couplings actually used the ' ...
        'requested grid policy.']);
end
convention = getf(opts, 'gridConvention', 'legacy_inclusive');
if isstring(convention) && isscalar(convention), convention = char(convention); end
gridOffset  = logical(reshape(getf(opts, 'gridOffset', [0 0 0]), 1, 3));
gammaPolicy = getf(opts, 'gammaPolicy', 'P_drop');
if isstring(gammaPolicy) && isscalar(gammaPolicy), gammaPolicy = char(gammaPolicy); end
gr = info.grid;
conflict = ~isequal(gr.convention, convention) || ~isequal(gr.offset, gridOffset) || ...
    ~isequal(gr.gammaPolicy, gammaPolicy) || ~isequal(gr.requested, reshape(grid, 1, 3));
if conflict
    error('invz:spectraGridConflict', ['the requested gridConvention/gridOffset/gammaPolicy/grid ' ...
        'dimensions conflict with the precomputed opts.info.grid provenance.']);
end
end
