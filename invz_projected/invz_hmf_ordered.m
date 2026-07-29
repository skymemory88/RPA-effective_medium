function [hmf_star, prof, trc] = invz_hmf_ordered(ion, T, Bx, Jnu_flat, opts)
%INVZ_HMF_ORDERED Jensen applied-field/H_MF self-consistency, spontaneous root (SS9.3, J 2.31-2.33).
%   h0(hmf) = int_0^hmf r(h') dh',   r = G0(0;h')/Gtil0(0;h')
% with Gtil0 built on the STATIC-CLOSURE K0 (invz_emt_static_ordered, P0-2), evaluated on
% fixed-field single-ion states (hz_fixed WITHOUT order -- P0-1, invz:hzFixed asserted).
% Spontaneous condition (zero applied longitudinal field): h0(hmf) = J0eff*<Jz>(hmf); the
% nonzero root is bracketed on a GEOMETRIC profile clustered at 0 (P1-4) and refined by
% bisection with DIRECT node evaluations to opts.tol_root. F(h)/h -> crit as h -> 0+
% (SS5), returned as prof.slope0. Returns NaN when no nonzero root exists, or when the
% separate bare (order=true) bracketing solve does not order.
%
% Status contract (round-5 P2, binding): 'unresolved' means ordering was PREDICTED
% (slope0 < 0) but no bracket was found above hmin_abs, OR the bracket did not refine
% to tol_root. 'node_failed' means ANY node evaluated along the way -- the h=0
% predictor, a profile node, a bisection iterate, or the final root evaluation --
% failed to converge/close. Both map to a NaN hmf_star and MUST be read by callers as
% converged = false, never as a PM label.
%
% Task 13 (binding precedence, the pure reducer invz_hmf_status.m): 'degenerate_doublet'
% means some evaluated node's two-level splitting Delta fell below the 1e-4 meV domain
% floor (invz_twolevel_ordered opts.domain_policy = 'return', spec SS5.3) -- a labelled
% domain outcome, not a solver failure. PUBLIC CONTRACT CHANGE: this function previously
% let that construction's invz:degenerateDoublet throw ESCAPE unchanged (so Bx = 0, in
% particular, threw); it now returns normally with hmf_star = NaN and this status instead.
% 'medium_out_of_domain' means some evaluated node's strict-scheme reference/closure event
% left medium_status outside {'ok','not_applicable'}. Both are prepended ABOVE
% 'node_failed'/'unresolved' in the binding precedence -- degenerate_doublet >
% medium_out_of_domain > node_failed > unresolved > ok -- so a domain/degenerate reason on
% ANY evaluated node (predictor, profile, extension, redensification, bisection iterate,
% or the final root) is never masked as a generic 'node_failed'/'unresolved' verdict. Both
% new statuses still map to a NaN hmf_star and MUST be read by callers as converged = false.
% prof.status_detail is [] on 'ok'; otherwise it carries the compact failed-node census and
% deterministic binding node returned by invz_hmf_status.
%
% Third output trc (stage-2c task 0, diagnostic-only; RELOCATED task 1b-ii-A): opts.trace --
% absent/false/empty by default -- gates a BEHAVIOUR-NEUTRAL trace of this function's node
% loop (eval_node/run_sweep, including the h=0 predictor call). Every trace statement is
% guarded by `if tracing`; hmf_star/prof are bit-identical whether or not opts.trace is
% set. eval_node's own per-node solve (Sigma<->EMT loop + post-loop static refresh) is now
% ONE call to the shared invz_ordered_node_solve (sopts.trace=tracing forwarded through);
% trc.nodes/trc.iters are rebuilt from that call's info/info.iters, field-for-field, in
% eval_node itself. The one UNCONDITIONAL change (now inside the solver's own per-iteration
% step) is capturing Gstat from the per-iteration static call instead of discarding it via
% `~`: that costs nothing regardless of tracing (invz_emt_static_ordered always computes
% it, independent of nargout). When
% tracing is off, trc is a small fixed placeholder (schema_version, enabled=false, empty
% meta/nodes/iters); callers that never request the 3rd output (all existing callers) are
% unaffected either way. When on: trc.nodes has one record per eval_node CALL (id, h,
% phase in {predictor,sweep,extend,redensify,bisect,root}, seed_kind cold/warm +
% seed_from node id, outer_iters, term_reason in {converged,max_iter,refresh_failed,
% bare_shortcut,medium_out_of_domain,degenerate_doublet} (the second-to-last is a
% strict-scheme reference/closure domain event, passed through verbatim by
% local_term_reason -- see its own docstring; the last is task 13's two-level Delta
% domain screen, set directly by eval_node BEFORE local_term_reason ever runs -- see the
% Status contract paragraph above), K0, D_uni,
% resid_static, PLUS every field of that call's per-node record -- see the diagnostics
% paragraph below; ONE nested finalizer, append_trace_node, builds every
% entry, so no exit path can append a shorter schema); trc.iters has one record per OUTER Sigma<->K
% iteration across all nodes (node_id back-link, the RAW map residual `dS` and the RAW
% static-closure residual `sout.resid`, dynamic x/denominator proximity, static y/coupling
% interval, K0, Gstat/local denominator/xi, D_uni, min/max/absolute-min/negative-count of
% Dq = 1+(J(q)-K0)*Gstat, closest-to-zero positive/negative Dq flat indices + values, and
% the static-closed/outer-open discriminator). trc.meta carries this call's own
% (T, Bx, J0eff, Jnu_flat, opts) -- q/branch
% lattice provenance is NOT available here (flat-Jnu interface, unchanged) and is
% reconstructed by the caller; see invz_projected/invz_ordered_trace.m, the packaging
% helper that adds qc/unflattened-Jnu/hashes and documents the versioned schema in full.
%
% Per-node diagnostics (task 12). G = -chi (meV^-1), ferromagnetic positive J. eval_node
% returns ONE fixed-schema record (blank_node_record) instead of a positional output list, and
% the solved path exports it wholesale, so prof's numeric arrays and its cell arrays are
% derived from the SAME record array in the SAME order at ONE point and cannot drift apart:
%   prof.crit              r + J0eff*G0bare -- the dimensionless mass, the per-node
%                          generalization of prof.slope0 (= crit at h = 0, SS5)
%   prof.r_minus_1         r - 1. NOT Sigma0: at finite ordered moment the hybrid elastic
%                          factor xi makes r depend on K0 and lambda(1:2) as well, so the two
%                          coincide only at m = 0 (spec SSA Consequence 2)
%   prof.Delta             the node's own two-level splitting (invz_twolevel_ordered)
%   prof.Dq_min            res.stability.Dq_min from the independent residual checker
%   prof.ref_denom         strict reference denominator (1 + Sigma0 for the Dyson reference)
%   prof.ref_margin        the ACTUAL distance to the configured floor (denom - ref_margin),
%                          never the denominator and never the floor itself
%   prof.gstat_local_denom 1 + Sigma0 + K0*G0inel0, Gstat's own signed local denominator
%   prof.omit_mu3/_cubic/_max   the strict closure's omitted-term ratios
%   prof.medium_status, prof.node_term_reason   cell arrays, one entry per node
%   prof.r_pm0/.G0bare_pm0/.Sigma0_pm0/.K0_pm0  the h = 0 predictor seeds
%   prof.G0bare_star/.crit_star/.Dq_min_star    the accepted root's own record
%   prof.int_Sigma0, prof.int_r_minus_1   path integrals over h, seeded on the SAME first
%                          panel as h0 so all three quadratures are consistent. BOTH are
%                          needed (spec SSA binding caution, gate G5): int Sigma0 dh is a
%                          component diagnostic, not the whole correction.
% All of these are diagnostics: none of them gates the root search or the status contract.
%
% opts.static_medium / opts.ref_margin are resolved ONCE per call by invz_check_static_medium,
% the sole scheme authority, which stamps both leg option structs (opts.emt -> the dynamic
% slot, opts.emt_static -> the ordered static sector). The scheme is never re-derived here.
% opts.Jmom (optional): invz_coupling_moments of the coupling multiset, resolved once per call
% and threaded into every node solve; a direct caller may omit it and it is derived here.
% opts.hmf_seed (optional): struct('Sigma',full Matsubara vector,'K0s',scalar), used only as
% the numerical initial guess for the h = 0 predictor. It changes neither the equations nor
% any acceptance gate. A failed warm attempt retains invz_ordered_node_solve's existing
% one-shot cold retry. prof.hmf_seed_out is populated only when the predictor passes the
% complete residual checker, allowing a caller to continue between physical fields.
% opts.cold_acceleration (optional, default 'none'): 'none' | 'signed_aitken1'. The latter
% is the safeguarded, default-off resummed cold-predictor experiment implemented by
% invz_ordered_node_solve. A warm attempt is never accelerated; its existing fresh-cold
% retry is eligible. The full summary for the h=0 predictor is returned in
% prof.predictor_acceleration.
if nargin < 5, opts = struct(); end
J0eff = opts.J0eff;                                  % required, no default (caller-owned)
Jxx0  = getf(opts, 'Jxx0', ion.Jxx0);
hyp   = getf(opts, 'hyp', true);
tmf   = getf(opts, 'transverse_mf', 'legacy_x');
nH    = getf(opts, 'nH', 33);
hfac  = getf(opts, 'hmax_fac', 1.25);
hfrac = getf(opts, 'hmin_frac', 1e-3);
trt   = getf(opts, 'tol_root', 1e-3);
fbare = getf(opts, 'force_bare', false);
mixo  = getf(opts, 'mix_outer', 0.7);
tolo  = getf(opts, 'tol_outer', 1e-8);
maxo  = getf(opts, 'max_outer', 200);                % ALIGNED with both solvers' default
cold_accel = getf(opts, 'cold_acceleration', 'none');
% Static-medium scheme, resolved ONCE by the SOLE authority (invz_check_static_medium, spec
% SS4.2) and stamped into BOTH leg option structs -- eopts (invz_emt_scalar's PM static slot)
% and eso (invz_emt_static_ordered's ordered static sector) -- so the two sectors can never run
% different truncation orders. Resolved HERE, ahead of the bare-order early return below, so
% prof.static_medium is stamped on every exit path. The scheme is NEVER re-derived locally, and
% eopts is NOT fetched a second time near the Matsubara block.
eopts = getf(opts, 'emt', struct());                 % dynamic-sector opts, threaded (P1-6)
eso   = getf(opts, 'emt_static', struct());          % static-closure opts, threaded (P1-F)
[sm, eopts, eso] = invz_check_static_medium(opts, eopts, eso);
eso.warn = false;   % node loop gates on so.converged; suppress the per-node console flood

% Coupling moments, resolved ONCE per call (spec SS4.3): recomputing them inside the outer
% iteration of every node would repeat an O(nJ) pass up to max_outer x nNodes times per field.
% A caller that already resolved the point's coupling spectrum (invz_solve_point_ordered) passes
% opts.Jmom straight through; a direct call derives it here as a compatibility fallback.
Jmom = getf(opts, 'Jmom', []);
if isempty(Jmom), Jmom = invz_coupling_moments(Jnu_flat); end
if sm.is_strict && ~isvector(Jnu_flat)
    error('invz:staticMedium', 'strict ordered HMF does not support [nJ,nw] Jnu_flat.');
end

% single-ion opts for FIXED-FIELD nodes: NO 'order' (P0-1); forward mf knobs (P1-6)
sibase = struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), sibase.(f{1}) = opts.(f{1}); end
end
hmin_abs = getf(opts, 'hmin_abs', NaN);              % resolved after hmax below (P1-C)

prof = struct('hgrid', [], 'r', [], 'h0', [], 'm', [], 'Sigma0', [], 'K0', [], ...
              'D_uni', [], 'G0bare', [], 'Gstat', [], 'node_conv', [], 'F', [], ...
              'crit', [], 'r_minus_1', [], 'Delta', [], 'Dq_min', [], ...
              'ref_denom', [], 'ref_margin', [], 'gstat_local_denom', [], ...
              'omit_mu3', [], 'omit_cubic', [], ...
              'omit_max', [], 'medium_status', {{}}, 'node_term_reason', {{}}, ...
              'slope0', NaN, 'r_pm0', NaN, 'G0bare_pm0', NaN, ...
              'Sigma0_pm0', NaN, 'K0_pm0', NaN, 'J0eff', J0eff, ...
              'hmf_seed_out', [], ...
              'cold_acceleration', cold_accel, 'predictor_acceleration', struct(), ...
              'n_extend', 0, 'hmin_initial', NaN, 'status', 'no_bare_order', ...
              'status_detail', [], ...
              'redensified', false, 'int_Sigma0', NaN, 'int_r_minus_1', NaN, ...
              'm_star', NaN, 'D_uni_star', NaN, 'r_star', NaN, 'Gstat_star', NaN, ...
              'G0bare_star', NaN, 'crit_star', NaN, 'Dq_min_star', NaN, ...
              'static_medium', sm.scheme);
hmf_star = NaN;

% --- Trace hook init (stage-2c task 0; diagnostic-only, DEFAULT OFF -- see the header
% block above). tracing/trc are assigned unconditionally, BEFORE the first possible early
% return below, so trc is always well-formed regardless of which exit path fires.
% node_seq/cur_phase are bookkeeping-only counters shared with eval_node/run_sweep below
% via ordinary nested-function workspace sharing -- the SAME mechanism this file already
% uses (read-only) for mixo/tolo/maxo/Jxx0/tmf etc.; no new production dependency.
tracing = isfield(opts, 'trace') && ~isempty(opts.trace) && ~isequal(opts.trace, false);
% schema_version 4: trc.nodes/iters gained the cold-acceleration summary/proposal fields.
% Version 3 added the compact failure-ledger and Phase-0 pole-proximity fields. Bumped because this field
% exists precisely to signal a trace-record shape change. No consumer in the repo reads/checks
% THIS struct's schema_version equality
% (grepped repo-wide): invz_ordered_trace.m never forwards trcRaw.schema_version into its
% own (separately-versioned) wrapper, so bumping it is inert everywhere today.
trc = struct('schema_version', 4, 'enabled', tracing, 'meta', struct(), ...
             'nodes', struct([]), 'iters', struct([]));
if tracing
    optsRec = opts;
    if isfield(optsRec, 'trace'), optsRec = rmfield(optsRec, 'trace'); end
    trc.meta = struct('T', T, 'Bx', Bx, 'J0eff', J0eff, 'Jnu_flat', Jnu_flat(:), 'opts', optsRec);
end
node_seq  = 0;    % running eval_node CALL counter (assigns nrec.id; bookkeeping only)
cur_phase = '';   % phase tag for the NEXT eval_node call, set at each call site below
% id/h/phase/seed provenance of the eval_node call currently in flight, captured at its ENTRY
% (BEFORE Sigma is touched) and read back by the single append_trace_node finalizer. Declared
% here so it is genuinely SHARED across both nested functions -- the SAME bookkeeping-only
% nested-workspace idiom node_seq/cur_phase already use; never read on the production path.
cur_node_meta = struct('id', NaN, 'h', NaN, 'phase', '', 'seed_kind', '', 'seed_from', NaN);

% Bracket ceiling from the BARE ordered fixed point: SAME MF option base plus
% order=true and J0z (P1-F -- the bracket runs under the caller's MF knobs too).
sibo = sibase;  sibo.order = true;  sibo.J0z = J0eff;
sib = invz_single_ion(ion, T, Bx, sibo);
if ~(sib.mf_converged && abs(sib.Jexp(3)) > 1e-6), return; end     % bare does not order
% Ceiling selection (kept SEPARATE from the bare-order check above -- third §7 review
% P1): opts.hmax_abs, when supplied, is an EXACT OVERRIDE (hmax = hmax_abs), never a
% `min` cap -- a cap would not give a COMMON cutoff across states with different bare
% hz (Task 6b Step 1). A nonpositive/nonfinite value is a caller error, not silently
% clamped.
if isfield(opts, 'hmax_abs')
    hmax_abs = opts.hmax_abs;
    if ~(isscalar(hmax_abs) && isfinite(hmax_abs) && hmax_abs > 0)
        error('invz:hmfOpts', 'hmax_abs must be a positive finite scalar: got %s', mat2str(hmax_abs));
    end
    hmax = hmax_abs;
else
    hmax = hfac * abs(sib.hz);
end
if isnan(hmin_abs), hmin_abs = 1e-10*hmax; end

% --- Matsubara grid, weights, beta: MIRROR invz_solve_point_ordered's setup block
% verbatim (wn, wts, beta -- honors Ecut and EMT options, P1-6; eopts was resolved and
% scheme-stamped once at the top of this function, so it is NOT re-fetched here).
Ecut  = getf(opts, 'Ecut', 40);
[wn, wts, beta] = invz_matsubara(T, Ecut);

% Independent h = 0 PM predictor node (round-3 P0-3; doubles as Gate 6b's comparator):
% ONE node solve at hz_fixed = 0 gives THIS machinery's PM fixed point. Its mass
%   slope_pred = r(0) + J0eff*G0bare(0) = 1 + Sigma0(0) - J0eff*chi_path(0)   (= crit, SS5)
% predicts root existence INDEPENDENTLY of any sampled profile value.
hmf_seed = getf(opts, 'hmf_seed', []);
if isempty(hmf_seed)
    Sigma = [];  K0s = 0;
elseif isstruct(hmf_seed) && isscalar(hmf_seed) && ...
        isfield(hmf_seed, 'Sigma') && isfield(hmf_seed, 'K0s') && ...
        isnumeric(hmf_seed.Sigma) && isreal(hmf_seed.Sigma) && ...
        numel(hmf_seed.Sigma) == numel(wn) && all(isfinite(hmf_seed.Sigma(:))) && ...
        isnumeric(hmf_seed.K0s) && isreal(hmf_seed.K0s) && isscalar(hmf_seed.K0s) && ...
        isfinite(hmf_seed.K0s)
    Sigma = hmf_seed.Sigma(:);  K0s = hmf_seed.K0s;
else
    error('invz:hmfOpts', ...
        ['hmf_seed must be empty or a finite real struct with Sigma matching the ' ...
         'Matsubara grid and scalar K0s.']);
end
if tracing, cur_phase = 'predictor'; end             % stage-2c task 0: node phase tag (bookkeeping only)
[pred, Sigma, K0s] = eval_node(0, Sigma, K0s);
pred.is_predictor = true;
prof.predictor_acceleration = pred.acceleration;
if pred.accepted
    prof.hmf_seed_out = struct('Sigma', Sigma, 'K0s', K0s);
end
% pred is the predictor's own fixed-schema record (task 13): passed to invz_hmf_status
% below/at every later early return, alongside whichever other records that site has in
% hand, so a degenerate/domain reason on THIS node is never masked as a generic failure.
r0n = pred.r;  S0pm = pred.Sigma0;  K0pm = pred.K0;  Gb0 = pred.G0bare; %#ok<NASGU> (S0pm/K0pm: parity locals)
ok0 = pred.accepted;
% A predictor-node convergence failure is NOT one of the three enumerated
% 'node_failed' triggers (round-5 P2: profile/bisection/final-evaluation), and it is
% NOT 'unresolved' either (that label presupposes a computed slope_pred). It is its
% own case: h0(h) = int_0^h r seeds on r(0), so without a converged predictor the
% cumulative integral (hence F, hence the root search) is categorically undefined --
% but the per-node grid quantities below (si/c0-derived G0bare in particular) are
% direct diagonalizations, INDEPENDENT of Sigma<->K convergence, so the grid is still
% built and exported for diagnostic use (measured: invz_solve_point's OWN default
% mix_outer = 0.7 oscillates rather than converging for transverse_mf = 'none' at a
% field where the m -> 0 reduction is genuinely marginal -- a pre-existing outer-loop
% conditioning fact shared by both solvers, not a Task-3 algebra issue). The overall
% verdict is still forced to 'node_failed' / NaN below, never silently 'ok'.
if ok0
    slope_pred = pred.crit;              % pred.crit IS r(0) + J0eff*G0bare(0) (task 12
                                          % record schema); read it rather than re-deriving
                                          % the same formula from r0n/Gb0 a second time.
    prof.Sigma0_pm0 = pred.Sigma0;  prof.K0_pm0 = pred.K0;  prof.slope0 = slope_pred;
    % Exported predictor seeds: r_pm0/G0bare_pm0 are the SAME two values slope0 is built from
    % (so prof.slope0 == prof.r_pm0 + J0eff*prof.G0bare_pm0 is a literal identity), and
    % r_pm0/Sigma0_pm0 are the h = 0 seeds of int_r_minus_1/int_Sigma0 below.
    prof.r_pm0 = r0n;  prof.G0bare_pm0 = Gb0;
else
    slope_pred = NaN;
end

[hgrid, ratio] = invz_hmf_grid(hmax, nH, hfrac);     % geometric, clustered at 0 (P1-4)
prof.hmin_initial = hgrid(1);

if tracing, cur_phase = 'sweep'; end                 % stage-2c task 0: node phase tag (bookkeeping only)
[nodes, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s);
[h0, F] = path_from_nodes(nodes, hgrid);

% ADAPTIVE lower extension (round-3 P0-3): predictor-driven, NOT self-referential.
% slope_pred < 0 predicts an ordered root; extend geometrically downward until a
% negative F sample appears or the absolute floor is reached.
n_extend = 0;
while ok0 && slope_pred < 0 && all(F >= 0) && hgrid(1) > hmin_abs
    n_extend = n_extend + 1;
    hext = hgrid(1) * ratio.^(3:-1:1);                % three more decades-fraction nodes
    erec = repmat(blank_node_record(), 1, 3);         % SAME schema as every other node
    if tracing, cur_phase = 'extend'; end             % stage-2c task 0: node phase tag (bookkeeping only)
    for k = 1:3
        [erec(k), Sigma, K0s] = eval_node(hext(k), Sigma, K0s);
    end
    hgrid = [hext hgrid];  nodes = [erec nodes];      % ONE array, prepended in grid order
    [h0, F] = path_from_nodes(nodes, hgrid);          % grid changed -> recompute h0/F
end

if n_extend > 0 && any(F < 0)
    % RE-DENSIFY (execution amendment 3, 2026-07-22): the extension's sparse geometric
    % panels feed O(coarse-grid) quadrature error into h0 exactly where F is a small
    % difference of large terms (measured: 11% root error at Bc_1z - 0.01 on a
    % deliberately coarse grid vs the fine default). Rebuild the profile at FULL nH
    % resolution anchored to the discovered bracket scale, so adaptive-path roots match
    % default-path quality. Cost: one extra nH-sweep, only when extension fired.
    idx0 = find(F < 0, 1, 'first');
    hfrac_eff = max(hmin_abs/hmax, 0.25*hgrid(idx0)/hmax);
    hgrid = invz_hmf_grid(hmax, nH, hfrac_eff);
    if tracing, cur_phase = 'redensify'; end          % stage-2c task 0: node phase tag (bookkeeping only)
    [nodes, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s);
    [h0, F] = path_from_nodes(nodes, hgrid);          % grid changed -> recompute h0/F
    prof.redensified = true;
end
prof.n_extend = n_extend;

% --- SINGLE derivation point for every solved-path array (task 12 step 8) -------------------
% `nodes` is now FINAL (extension and redensification are both behind us), so the numeric
% vectors AND the two cell arrays below are all read off THIS one struct array in THIS one
% order. There is deliberately no second, separately-maintained set of per-node vectors:
% that is exactly how prof.medium_status / prof.node_term_reason would drift out of alignment
% with prof.r / prof.m. h0/F come from path_from_nodes on the SAME array and grid.
rv = [nodes.r];      mv = [nodes.m];       S0v = [nodes.Sigma0];  K0v = [nodes.K0];
Dv = [nodes.D_uni];  Gbv = [nodes.G0bare]; Gsv = [nodes.Gstat];   cnv = [nodes.accepted];
prof.hgrid = hgrid;  prof.r = rv;  prof.h0 = h0;  prof.m = mv;
prof.Sigma0 = S0v;   prof.K0 = K0v;  prof.D_uni = Dv;  prof.node_conv = cnv;  prof.F = F;
prof.G0bare = Gbv;   prof.Gstat = Gsv;
prof.crit = [nodes.crit];  prof.r_minus_1 = rv - 1;  prof.Delta = [nodes.Delta];
prof.Dq_min = [nodes.Dq_min];  prof.ref_denom = [nodes.ref_denom];
prof.ref_margin = [nodes.ref_margin];
prof.gstat_local_denom = [nodes.gstat_local_denom];
prof.omit_mu3 = [nodes.omit_mu3];
prof.omit_cubic = [nodes.omit_cubic];  prof.omit_max = [nodes.omit_max];
prof.medium_status = {nodes.medium_status};  prof.node_term_reason = {nodes.term_reason};
% Path corrections (spec SSA binding caution, gate G5). BOTH are needed: at finite ordered
% moment r - 1 is NOT Sigma0 -- the hybrid elastic factor xi makes r depend on K0 and
% lambda(1:2) -- so integral Sigma0 dh is a component diagnostic, not the whole correction.
% The ~0.3% PM boundary shift bounds NEITHER of them deep in the ordered phase.
% Same first-panel seeding as h0 above, so all three integrals are quadrature-consistent, and
% computed here from the FINAL adapted grid so no superseded-grid integral is ever published.
if ok0
    prof.int_Sigma0    = trapz([0 hgrid], [prof.Sigma0_pm0, S0v]);
    prof.int_r_minus_1 = trapz([0 hgrid], [r0n - 1, rv - 1]);
end

% Task 13: ONE pure precedence reducer replaces the old ~ok0 / slope_pred<0 / any(~cnv)
% chain -- degenerate_doublet > medium_out_of_domain > node_failed > unresolved > ok. A
% predictor failure (~ok0) is still correctly labelled 'node_failed' here (pred.accepted
% is false, and neither of the two domain reasons applies), UNLESS pred itself carries a
% degenerate/domain reason, in which case that label now correctly takes precedence
% instead of masking as a generic 'node_failed' verdict.
[prof.status, prof.status_detail] = invz_hmf_status(pred, nodes, slope_pred, F);
if strcmp(prof.status, 'unresolved')
    warning('invz:hmfUnresolved', ...
        'ordering predicted (slope_pred = %.3g) but no negative F above hmin_abs = %.3g', ...
        slope_pred, hmin_abs);
end
if ~strcmp(prof.status, 'ok'), return; end
% ROOT-SELECTION CONVENTION (documented 2026-07-29, execution packet S5; behaviour
% UNCHANGED). Two separate choices are made on this line, and only the first is derived:
%   (1) sign pattern (- -> >= 0): the INCREASING crossing of F, i.e. a Landau minimum
%       rather than a maximum. This is the derived requirement; see the crit_star
%       comment near the end of this function.
%   (2) 'last': among ALL increasing crossings present on the h grid, the LARGEST-h one
%       is taken. This is a CONVENTION, not a derived equilibrium rule. The projected
%       ordered equations are known to admit several admissible roots at one operating
%       point -- seven distinct h = 0 roots with folds at 1.5 T, and two states that both
%       pass the full A-D contract at 4.05 T (invzp_convg_diagnosis.md SS6.1, SS6.2) --
%       so when more than one increasing crossing exists this line silently selects one
%       of several stationary states, with no thermodynamic ranking behind the choice.
% It must NOT be replaced by argmin of any integrated-potential score (Phi = int F dm or
% similar): that swaps one uncertified convention for another, and the integrated route
% is outside its validated domain as an absolute selector (blind_convg_plan.md SS2.1-2.3,
% SS5). The correct replacement is the stationary functional's own selection rule
% (invzp_convg_fix.md WP7), and nothing before it.
s = sign(F);  idx = find(s(1:end-1) < 0 & s(2:end) >= 0, 1, 'last');
if isempty(idx), return; end                          % no nonzero root: PM side

% --- Root refinement by DIRECT evaluation (P1-4): bisection on F between the
% bracketing nodes, fresh node solve per iterate, cumulative h0 via local trapezoid
% panel from the bracket's left node.
a = hgrid(idx);  b = hgrid(idx+1);  Fa = F(idx);  h0a = h0(idx);  ra = rv(idx);
if tracing, cur_phase = 'bisect'; end                % stage-2c task 0: node phase tag (bookkeeping only)
for it = 1:12
    c = 0.5*(a + b);
    [recc, Sigma, K0s] = eval_node(c, Sigma, K0s);
    if ~recc.accepted                                 % round-5 P1-A: a failed bisection node
        % Task 13: this iterate may itself be the first node to carry a domain/degenerate
        % reason (Delta is re-screened on every eval_node call, not only the profile) --
        % the reducer promotes that above a generic 'node_failed' verdict when it applies.
        recc.in_bracket = true;
        [prof.status, prof.status_detail] = invz_hmf_status(pred, [nodes recc], slope_pred, F);
        hmf_star = NaN;                               % TERMINATES the solve -- never a root
        return;                                       % from a partial bracket
    end
    h0c = h0a + 0.5*(ra + recc.r)*(c - a);
    Fc  = h0c - J0eff*recc.m;
    if sign(Fc) == sign(Fa), a = c; Fa = Fc; h0a = h0c; ra = recc.r; else, b = c; end
    if (b - a) < trt*b, break; end
end
if (b - a) >= trt*b                                   % round-5 P1-A: tol_root not reached --
    % Task 13: every bisection iterate reaching this line was accepted (else the branch
    % above already returned), so the reducer can only confirm 'ok' here -- this site's
    % own default, 'unresolved' (a distinct refinement failure, NOT a node failure), then
    % applies. Routed through the same reducer regardless, so a future domain/degenerate
    % reason on an accepted node would still correctly override it.
    recc.in_bracket = true;
    [prof.status, prof.status_detail] = invz_hmf_status( ...
        pred, [nodes recc], slope_pred, F, 'unresolved');
    hmf_star = NaN;
    warning('invz:hmfUnresolved', 'root bracket not refined to tol_root: (b-a)/b = %.3g', (b-a)/b);
    return;
end
hmf_star = 0.5*(a + b);
if tracing, cur_phase = 'root'; end                  % stage-2c task 0: node phase tag (bookkeeping only)
[root, ~, ~] = eval_node(hmf_star, Sigma, K0s);
root.in_bracket = true;
if ~root.accepted
    % Task 13: the final root evaluation is also screened -- a domain/degenerate reason
    % here must not be masked as a generic 'node_failed' verdict either.
    [prof.status, prof.status_detail] = invz_hmf_status(pred, [nodes root], slope_pred, F);
    hmf_star = NaN;
    return;
end
% the root's COMPLETE record, G0bare included (it was previously discarded here, which is
% why crit_star was not computable): crit_star = r_star + J0eff*G0bare_star must be > 0 for
% an accepted root, i.e. the INCREASING crossing of F -- a Landau minimum, not a maximum.
prof.m_star = root.m;  prof.D_uni_star = root.D_uni;  prof.Dq_min_star = root.Dq_min;
prof.r_star = root.r;  prof.Gstat_star = root.Gstat;  prof.G0bare_star = root.G0bare;
prof.crit_star = root.crit;

    function rec = blank_node_record()
    % The FIXED per-node record schema (task 12 step 7) -- the ONE initializer every eval_node
    % exit path starts from, and the one run_sweep/the extension preallocate with. Because the
    % field set, order and defaults are identical on every exit (normal, fbare shortcut, and any
    % future domain / degenerate-doublet exit), the records concatenate into one struct array and
    % prof's arrays can be read off it field-for-field. NEVER build a partial record by hand.
    % G = -chi (meV^-1), ferromagnetic positive J.
    rec = struct('h',NaN,'r',NaN,'m',NaN,'Sigma0',NaN,'K0',NaN, ...
        'D_uni',NaN,'Dq_min',NaN,'Dq_abs_min',NaN, ...
        'G0bare',NaN,'Gstat',NaN,'accepted',false,'crit',NaN,'Delta',NaN, ...
        'medium_status','not_applicable','term_reason','not_evaluated', ...
        'ref_denom',NaN,'ref_margin',NaN,'gstat_local_denom',NaN, ...
        'omit_mu3',NaN,'omit_cubic',NaN,'omit_max',NaN, ...
        'outer_iters',0,'resid_static',NaN, ...
        'resid_A',NaN,'resid_B',NaN,'resid_C',NaN,'resid_D',NaN, ...
        'resid_norm',NaN,'is_predictor',false,'in_bracket',false, ...
        'acceleration',blank_acceleration_record());
    end

    function [h0, F] = path_from_nodes(nodes, hgrid)
    % h0(h) = int_0^h r dh' and F = h0 - J0eff*m on the CURRENT grid, always recomputed from
    % the SAME record array the export block reads. Called after every grid change (initial
    % sweep, each extension prepend, redensification) so a published integral can never
    % describe a superseded grid. First panel seeded with the h = 0 predictor's r (P0-3).
    if ok0
        h0 = cumtrapz([0 hgrid], [r0n [nodes.r]]);  h0 = h0(2:end);
    else
        h0 = nan(1, numel(hgrid));                 % undefined without a real r(0)
    end
    F = h0 - J0eff*[nodes.m];
    end

    function append_trace_node(rec, info)
    % THE single trc.nodes finalizer (task 12 step 7). EVERY eval_node exit routes through here
    % when tracing -- normal, fbare shortcut, and any future domain / degenerate exit -- so no
    % early return can hand-append a shorter schema. Each entry is: this call's id/h/phase/seed
    % provenance (cur_node_meta, captured at eval_node ENTRY before Sigma is touched), the outer
    % loop provenance read from `info`, then EVERY field of the per-node record verbatim. Gate 0
    % reads this ledger to account for the predictor, extension/redensification nodes, bisection
    % iterates and the final root -- a missing or short entry silently corrupts that gate.
    % info = [] for the fbare shortcut: no solver call was made, so there is no loop provenance.
    nrec = cur_node_meta;                          % id, h, phase, seed_kind, seed_from
    nrec.outer_iters = 0;  nrec.outer_hit_max = false;  nrec.dS_break = false;
    nrec.resid_static = NaN;
    if ~isempty(info)
        nrec.outer_iters   = info.outer_iters;
        nrec.dS_break      = info.loop_converged;
        nrec.outer_hit_max = ~info.loop_converged && info.outer_iters == maxo;
        nrec.resid_static  = info.so.resid;
    end
    nrec.ok_final    = rec.accepted;               % == info.accepted on the solver path
    for fn = fieldnames(rec).'                     % the full record, field for field
        nrec.(fn{1}) = rec.(fn{1});                % (NOT `f`: that name is a parent loop index,
    end                                            % and nested functions SHARE parent variables)
    if isempty(trc.nodes), trc.nodes = nrec; else, trc.nodes(end+1) = nrec; end
    end

    function [recs, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s)
    % Behavior-neutral factoring of the per-node nH sweep (execution amendment 3):
    % evaluate eval_node at every point of hgrid in ascending order, warm-starting
    % Sigma/K0s across calls. Shared by the initial profile sweep and the re-densify
    % pass; the extension's 3-node prepends stay inline (unchanged). Returns the struct
    % ARRAY of per-node records plus the two continuation carriers.
    recs = repmat(blank_node_record(), 1, numel(hgrid));
    for is = 1:numel(hgrid)
        [recs(is), Sigma, K0s] = eval_node(hgrid(is), Sigma, K0s);
    end
    end

    function [rec, Sigma, K0s] = eval_node(hp, Sigma, K0s)
    % One fixed-field node: si (hz_fixed, NO order), tl, c0/G0, then the ordered
    % Sigma<->EMT loop WITH the static-sector closure each pass (Interfaces bullet).
    % Returns ONE fixed-schema record (blank_node_record) plus the continuation carriers --
    % NOT a positional output list; task 12 collapsed the previous ten outputs precisely
    % because every new diagnostic grew it and every call site destructured a different subset.
    rec = blank_node_record();
    rec.h = hp;
    if tracing
        rec.is_predictor = strcmp(cur_phase, 'predictor');
        rec.in_bracket = any(strcmp(cur_phase, {'bisect', 'root'}));
    end
    if tracing                                     % stage-2c task 0: node identity + seed
        node_seq = node_seq + 1;                   % provenance, captured BEFORE Sigma is
        cur_node_meta = struct('id', node_seq, 'h', hp, 'phase', cur_phase, ...     % touched
            'seed_kind', 'warm', 'seed_from', node_seq - 1);   % single sequential thread
        if isempty(Sigma)                          % this file's OWN cold-start criterion
            cur_node_meta.seed_kind = 'cold';  cur_node_meta.seed_from = NaN;
        end
    end
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bx, sio);
    if abs(si.hz - hp) > 1e-12
        error('invz:hzFixed', 'hz_fixed not held: si.hz = %.6g vs %.6g', si.hz, hp);
    end
    % Delta domain screen (spec SS5.3), SINGLE evaluation: return mode reuses the one
    % diagonalization the constructor already performs, instead of pre-screening with a
    % duplicate one and then calling it again. Delta is measured at THIS node's molecular field
    % hp, and the geometric grid clusters at 0, so the predictor and lowest nodes are the ones
    % at risk whenever Bx is small -- not only at exactly Bx = 0. Previously the constructor's
    % throw escaped this function entirely and the column masked as a solver failure.
    tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', Jxx0, 'transverse_mf', tmf, ...
                                                      'domain_policy', 'return'));
    if ~tl.valid
        rec.m = si.Jexp(3);  rec.Delta = tl.Delta;
        rec.accepted = false;  rec.term_reason = 'degenerate_doublet';
        if tracing, append_trace_node(rec, []); end  % Task-12 single schema finalizer
        return;
    end
    rec.m = si.Jexp(3);
    rec.Delta = tl.Delta;                          % this node's own doublet splitting
    if fbare
        rec.r = 1;  rec.Sigma0 = 0;  rec.K0 = 0;  rec.accepted = true;
        rec.term_reason = 'bare_shortcut';
        rec.crit = rec.r + J0eff*rec.G0bare;       % G0bare is not computed here -> NaN crit
        if tracing, append_trace_node(rec, []); end % degenerate node, same full schema
        return;
    end
    c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
    G0 = -real(squeeze(c0(3,3,:)));
    c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));   % static inelastic only
    G0inel0 = -real(c0i(3,3,1));                                   % fixed-Hamiltonian slot
    % PATH-CONSISTENT bare static (SS4a, round-3 P0-2; round-4 P1-A): -dm/dh ALONG the
    % node path. Transverse-MF feedback by mode: 'none' forces hx = hy = 0 in
    % invz_single_ion, so the path derivative is X_zz alone; 'legacy_x' feeds back the
    % x channel; 'vector_ab' the 2x2 {x,y} block.
    X = real(c0(:, :, 1));                                         % static chi tensor (chi = -G)
    switch tmf
        case 'none'
            fb = 0;
        case 'legacy_x'
            fb = X(3, 1) * (Jxx0 / (1 - Jxx0*X(1, 1))) * X(1, 3);
        case 'vector_ab'
            t = [1 2];
            fb = X(3, t) * (Jxx0 * ((eye(2) - Jxx0*X(t, t)) \ X(t, 3)));
        otherwise
            error('invz:transverseMF', 'unknown transverse_mf ''%s''', tmf);
    end
    G0bare0 = -(X(3, 3) + fb);
    G0el0   = G0bare0 - G0inel0;                                   % elastic + feedback (SS4a)
    g  = real(invz_g(tl, 1i*wn));
    % --- stage-2c task 1b-ii-A: ONE call to the shared, checker-gated node solver ----------
    % (replaces the old inline cold-init + 7-step loop + post-loop refresh + ad hoc ctol
    % gate: invz_ordered_node_solve runs that SAME map verbatim and gates acceptance on the
    % complete residual checker, invz_ordered_residual -- see both files' headers.)
    Sigma_in = Sigma;  K0s_in = K0s;      % P0-3: capture the INPUT seed BEFORE the call --
                                          % the last-good continuation state to roll back to
                                          % if this node is not accepted (a failed node must
                                          % NOT export its own iterate as the next warm start).
    node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
        'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
        'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu_flat, 'Jmom', Jmom);
    % Jmom (resolved ONCE at the top of this function) is REQUIRED by both the node solver and
    % the residual checker under a strict scheme and IGNORED by both under 'resummed' -- so it
    % is threaded unconditionally rather than branching on the scheme at every node.
    if isempty(Sigma_in)                 % this file's OWN cold-start criterion (unchanged)
        seed = [];
    else
        seed = struct('Sigma', Sigma_in, 'K0s', K0s_in);
    end
    sopts = struct('mix_outer', mixo, 'max_outer', maxo, 'tol_outer', tolo, ...
        'cold_retry', true, 'trace', tracing, 'cold_acceleration', cold_accel);
    [state, info] = invz_ordered_node_solve(node, seed, sopts);

    % checker-gated acceptance: info.accepted is the COMPLETE invz_ordered_residual verdict,
    % replacing the old in-loop `dS<tolo && sout.converged` + post-loop `so.converged &&
    % so.resid<ctol` gate -- that logic now lives entirely in the solver/checker.
    rec.r = info.so.r;  rec.Sigma0 = state.Sigma(1);  rec.K0 = state.K0s;
    rec.D_uni = info.so.D_uni;  rec.G0bare = G0bare0;  rec.Gstat = info.so.Gstat;
    rec.accepted = info.accepted;
    rec.outer_iters = info.outer_iters;  rec.resid_static = info.so.resid;
    rec.resid_A = info.res.blockA.resid;  rec.resid_B = info.res.blockB.resid;
    rec.resid_C = info.res.blockC.resid;  rec.resid_D = info.res.blockD.resid;
    block_resid = [rec.resid_A rec.resid_B rec.resid_C rec.resid_D];
    block_scale = [info.res.blockA.scale_abs info.res.blockB.scale_abs ...
                   info.res.blockC.scale_abs info.res.blockD.scale_abs];
    if all(isfinite(block_resid)) && all(isfinite(block_scale)) && all(block_scale > 0)
        rec.resid_norm = max(abs(block_resid)./block_scale);
    end
    % crit = r + J0eff*G0bare (SS5): the SAME expression prof.slope0 uses at h = 0, evaluated
    % per node. Diagnostic only -- intermediate path nodes legitimately sit in the unstable
    % Landau interval (crit < 0); only the accepted root is held to positivity.
    rec.crit = rec.r + J0eff*rec.G0bare;
    rec.Dq_min = info.res.stability.Dq_min;        % the CHECKER's independent recomputation
    if ~rec.accepted
        Dq_final = 1 + (Jnu_flat(:) - state.K0s).*info.so.Gstat;
        if all(isfinite(Dq_final)), rec.Dq_abs_min = min(abs(Dq_final)); end
    end
    rec.gstat_local_denom = info.so.gstat_local_denom;
    rec.omit_mu3 = info.so.omit_mu3;  rec.omit_cubic = info.so.omit_cubic;
    rec.omit_max = info.so.omit_max;
    rec.medium_status = info.medium_status;
    rec.term_reason   = local_term_reason(info.term_reason);
    rec.acceleration = info.acceleration;
    % Reference denominator + ACTUAL distance-to-floor (denom - ref_margin), not the floor and
    % not the denominator repeated. Both stay NaN under 'resummed', where medium.ref is [].
    if isstruct(info.medium) && isfield(info.medium, 'ref') && isstruct(info.medium.ref)
        rec.ref_denom  = info.medium.ref.denom;
        rec.ref_margin = info.medium.ref.margin;
    end

    % P0-3 seed-safety: thread the ACCEPTED state forward; a non-accepted node instead
    % returns the INPUT Sigma/K0s UNCHANGED (rollback -- a cold input, i.e. Sigma_in = [],
    % rolls back to cold, matching the pre-existing cold criterion for the next node).
    if rec.accepted
        Sigma = state.Sigma;  K0s = state.K0s;
    else
        Sigma = Sigma_in;  K0s = K0s_in;
    end

    if tracing                                     % stage-2c task 1b-ii-A: node record
        append_trace_node(rec, info);              % task 12: the ONE finalizer, full schema

        for ii = 1:numel(info.iters)                % per-outer-iteration records, relocated
            src = info.iters(ii);                   % from info.iters (the solver's own trace)
            irec = struct('node_id', cur_node_meta.id, 'outer', src.outer, 'resid_map', src.dS, ...
                'resid_static', src.resid_static, 'K0', src.K0, 'D_uni', src.D_uni, ...
                'Dq_min', src.Dq_min, 'Dq_max', src.Dq_max, 'Dq_neg_count', src.Dq_neg_count, ...
                'Dq_abs_min', src.Dq_abs_min, ...
                'idx_pos_flat', src.idx_pos_flat, 'Dq_pos_val', src.Dq_pos_val, ...
                'idx_neg_flat', src.idx_neg_flat, 'Dq_neg_val', src.Dq_neg_val, ...
                'converged_flag', src.converged_flag, ...
                'x', src.x, 'x_minus_Jmax', src.x_minus_Jmax, ...
                'Ddyn_abs_min', src.Ddyn_abs_min, ...
                'Ddyn_nonpositive_count', src.Ddyn_nonpositive_count, ...
                'y', src.y, 'y_rank', src.y_rank, ...
                'y_interval_lo', src.y_interval_lo, 'y_interval_hi', src.y_interval_hi, ...
                'Gstat', src.Gstat, 'gstat_local_denom', src.gstat_local_denom, ...
                'xi', src.xi, 'static_closed_outer_open', src.static_closed_outer_open, ...
                'accel_attempted', src.accel_attempted, ...
                'accel_accepted', src.accel_accepted, ...
                'accel_lambda', src.accel_lambda, ...
                'accel_lambda_spread', src.accel_lambda_spread, ...
                'accel_mode_fit_rel', src.accel_mode_fit_rel, ...
                'accel_resid_ordinary', src.accel_resid_ordinary, ...
                'accel_resid_proposal', src.accel_resid_proposal, ...
                'accel_interval_rank', src.accel_interval_rank, ...
                'accel_reject_reason', src.accel_reject_reason);
            if isempty(trc.iters), trc.iters = irec; else, trc.iters(end+1) = irec; end
        end
    end
    end
end

% =============================================================================================
function out = blank_acceleration_record()
out = struct('mode', 'none', 'enabled', false, ...
    'attempted', 0, 'accepted', 0, 'accepted_outers', [], ...
    'lambda', NaN, 'lambda_spread', NaN, 'mode_fit_rel', NaN, ...
    'resid_ordinary', NaN, 'resid_proposal', NaN, ...
    'interval_rank', NaN, 'last_reject_reason', '');
end

% =============================================================================================
function tr = local_term_reason(info_reason)
%LOCAL_TERM_REASON invz_ordered_node_solve's vocabulary -> THIS file's own per-node vocabulary
% {converged, max_iter, refresh_failed, bare_shortcut} (an unmapped solver reason, e.g.
% 'medium_out_of_domain', passes through verbatim). ONE mapping serves both consumers:
% prof.node_term_reason and trc.nodes.term_reason, so the profile and the trace ledger can
% never disagree about a node's outcome.
switch info_reason
    case 'accepted',                    tr = 'converged';
    case 'loop_converged_not_accepted', tr = 'refresh_failed';
    case 'max_iter',                    tr = 'max_iter';
    otherwise,                          tr = info_reason;
end
end
