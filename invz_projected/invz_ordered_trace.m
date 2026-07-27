function trace = invz_ordered_trace(ion, T, Bx, Jnu_flat, opts)
%INVZ_ORDERED_TRACE Diagnostic wrapper: runs invz_hmf_ordered (the jensen ordered leg's
% production node loop, invz_solve_point_ordered's ordered_mode='jensen') with its
% behaviour-neutral trace hook (opts.trace, stage-2c task 0) and packages a versioned,
% self-contained struct -- run-level metadata stored ONCE + per-node/per-iteration records
% -- intended for .mat storage (schema_version-tagged; NOT printed text parsed back).
%
% This helper calls invz_hmf_ordered DIRECTLY (not invz_solve_point_ordered): the trace hook
% lives entirely inside invz_hmf_ordered's own node loop (see its header for the full field
% contract), and calling it directly keeps this helper decoupled from invz_solve_point_ordered's
% richer (untouched) acceptance logic -- pt.is_ordered/pt.converged/pt.Sigma0 are NOT
% reproduced or re-derived here, only trace.result.hstar/hmf_status (invz_hmf_ordered's own
% root + profile status). For the pt-level fields, call invz_solve_point_ordered separately
% on the SAME (ion, T, Bx, Jnu_flat, solve opts): both functions are deterministic, so the two
% calls agree by construction -- this is exactly the behaviour-neutrality property required
% here: a traced invz_hmf_ordered call must reproduce the untraced production result
% bit-for-bit.
%
% ion, T, Bx, Jnu_flat: EXACTLY the production positional arguments (flat Jnu_flat -- see
%   invz_bz_couplings.m:17, Jnu = Jnu(:)). This interface is UNCHANGED (stage-2c task 0 does
%   not touch it); q/branch provenance is reconstructed here from the CALLER-supplied
%   qc/Jnu_unflat below, never re-derived from Jnu_flat alone.
%
% opts (all fields optional):
%   .qc          [nq x 3]  q-grid used to build Jnu_flat (invz_jq_modes/invz_bz_couplings),
%                          recorded for PROVENANCE only -- never consumed by the solve.
%   .Jnu_unflat  [nq x 4]  UNFLATTENED branch matrix, i.e. invz_jq_modes' own Jnu output
%                          BEFORE invz_bz_couplings.m:17's Jnu = Jnu(:). Supplied together
%                          with qc, trace.meta can resolve a flat index to (q-index,
%                          branch-index, J(q)) -- see invz_ordered_trace_resolve.m. Omit both
%                          for a synthetic (non-lattice) Jnu_flat; trace.meta.is_synthetic is
%                          then true and qc/Jnu_unflat are recorded empty.
%   .grid, .dpRng          recorded verbatim into trace.meta (this helper does not itself
%                          call invz_jq_modes/invz_bz_couplings -- the caller already did, to
%                          get Jnu_flat/qc/Jnu_unflat in the first place).
%   .lattice_info  struct  the invz_jq_modes/invz_bz_couplings info struct, recorded verbatim
%                          (e.g. carries Jcc0/Jaa0/dpRng) when available.
%   .solve         struct  forwarded to invz_hmf_ordered AS ITS OWN opts (J0eff -- REQUIRED,
%                          no default, exactly as invz_hmf_ordered itself demands -- Jxx0,
%                          hyp, transverse_mf, nH, ...). opts.solve.trace is driver-owned
%                          here (forced on) and must not be set by the caller.
%   .save_path     char    if given, `save(save_path, 'trace')` -- the .mat IS the schema;
%                          nothing here is printed text meant to be parsed back.
%
% Returns trace (schema_version = 2):
%   .schema_version  (double) versioning for downstream (Task 2) consumers.
%   .meta            run-level, stored ONCE: T, Bx, J0eff (the value invz_hmf_ordered
%                    actually used), solve_opts (verbatim), qc, Jnu_unflat, nq,
%                    is_synthetic, grid, dpRng, lattice_info, Jnu_hash, lattice_hash.
%   .nodes, .iters   invz_hmf_ordered's OWN per-node / per-outer-iteration trace records,
%                    copied verbatim -- see invz_hmf_ordered.m's header for the full field
%                    list (phase/seed provenance/term_reason per node; raw map + static
%                    residuals, dynamic/static pole coordinates, K0, Gstat, D_uni, Dq
%                    min/max/absolute-min/neg-count + closest-mode flat indices per iteration).
%   .result          struct('hstar', ..., 'hmf_status', ...) -- invz_hmf_ordered's own
%                    verdict for this call.
if nargin < 5, opts = struct(); end
qc         = getf(opts, 'qc', []);
Jnu_unflat = getf(opts, 'Jnu_unflat', []);
grid_      = getf(opts, 'grid', []);
dpRng_     = getf(opts, 'dpRng', NaN);
latinfo    = getf(opts, 'lattice_info', struct());
solve_opts = getf(opts, 'solve', struct());
save_path  = getf(opts, 'save_path', '');
if isfield(solve_opts, 'trace')
    error('invz:orderedTrace', ...
        'opts.solve.trace is driver-owned (invz_ordered_trace forces it on) -- do not set it.');
end
Jnu_flat = Jnu_flat(:);

% --- the ONLY solve call: invz_hmf_ordered directly, tracing forced on ---------------------
hopts = solve_opts;  hopts.trace = true;
[hstar, hprof, trcRaw] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, hopts);

% --- run-level metadata, stored ONCE -------------------------------------------------------
is_synthetic = isempty(qc) || isempty(Jnu_unflat);
nq = NaN;
if ~is_synthetic
    nq = size(Jnu_unflat, 1);
    if size(Jnu_unflat, 2) ~= 4 || size(qc, 1) ~= nq
        error('invz:orderedTrace', ['qc/Jnu_unflat size mismatch: qc is [%d x %d], ' ...
            'Jnu_unflat is [%d x %d] (need [nq x 3] / [nq x 4]).'], ...
            size(qc,1), size(qc,2), size(Jnu_unflat,1), size(Jnu_unflat,2));
    end
    if numel(Jnu_flat) ~= numel(Jnu_unflat)
        error('invz:orderedTrace', ['Jnu_flat (numel %d) does not match numel(Jnu_unflat) ' ...
            '(%d): not the same coupling set.'], numel(Jnu_flat), numel(Jnu_unflat));
    end
end
meta = struct('T', T, 'Bx', Bx, 'J0eff', hopts.J0eff, 'solve_opts', solve_opts, ...
    'qc', qc, 'Jnu_unflat', Jnu_unflat, 'nq', nq, 'is_synthetic', is_synthetic, ...
    'grid', grid_, 'dpRng', dpRng_, 'lattice_info', latinfo, ...
    'Jnu_hash', weak_hash(Jnu_flat), 'lattice_hash', '');
if ~is_synthetic
    meta.lattice_hash = weak_hash([qc(:); dpRng_(:)]);
end

trace = struct('schema_version', 2, 'meta', meta, ...
    'nodes', trcRaw.nodes, 'iters', trcRaw.iters, ...
    'result', struct('hstar', hstar, 'hmf_status', hprof.status));

if ~isempty(save_path)
    save(save_path, 'trace');
end
end

function h = weak_hash(v)
% Weak, deterministic digest -- mirrors invz_common/invz_cache_key.m's hash_vec /
% invz_jq_modes.m's own local hash_vec (32-bit truncation of a single-precision,
% index-weighted sum). Collisions are possible: this is a cache-style run fingerprint for
% provenance/identification, not a security hash.
v = v(:);
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v)).')), 'uint32'));
end
