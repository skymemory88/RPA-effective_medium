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
% bare_shortcut}, K0, D_uni, resid_static); trc.iters has one record per OUTER Sigma<->K
% iteration across all nodes (node_id back-link, the RAW map residual `dS` and the RAW
% static-closure residual `sout.resid`, K0, D_uni, min/max/neg-count of the static
% Dq = 1+(J(q)-K0)*Gstat, and the closest-to-zero positive/negative Dq flat indices +
% values). trc.meta carries this call's own (T, Bx, J0eff, Jnu_flat, opts) -- q/branch
% lattice provenance is NOT available here (flat-Jnu interface, unchanged) and is
% reconstructed by the caller; see invz_projected/invz_ordered_trace.m, the packaging
% helper that adds qc/unflattened-Jnu/hashes and documents the versioned schema in full.
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
eso   = getf(opts, 'emt_static', struct());          % static-closure opts, threaded (P1-F)
eso.warn = false;   % node loop gates on so.converged; suppress the per-node console flood

% single-ion opts for FIXED-FIELD nodes: NO 'order' (P0-1); forward mf knobs (P1-6)
sibase = struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mf_maxit', 'mf_mix'}
    if isfield(opts, f{1}), sibase.(f{1}) = opts.(f{1}); end
end
hmin_abs = getf(opts, 'hmin_abs', NaN);              % resolved after hmax below (P1-C)

prof = struct('hgrid', [], 'r', [], 'h0', [], 'm', [], 'Sigma0', [], 'K0', [], ...
              'D_uni', [], 'G0bare', [], 'Gstat', [], 'node_conv', [], 'F', [], ...
              'slope0', NaN, 'Sigma0_pm0', NaN, 'K0_pm0', NaN, 'J0eff', J0eff, ...
              'n_extend', 0, 'hmin_initial', NaN, 'status', 'no_bare_order', ...
              'redensified', false, ...
              'm_star', NaN, 'D_uni_star', NaN, 'r_star', NaN, 'Gstat_star', NaN);
hmf_star = NaN;

% --- Trace hook init (stage-2c task 0; diagnostic-only, DEFAULT OFF -- see the header
% block above). tracing/trc are assigned unconditionally, BEFORE the first possible early
% return below, so trc is always well-formed regardless of which exit path fires.
% node_seq/cur_phase are bookkeeping-only counters shared with eval_node/run_sweep below
% via ordinary nested-function workspace sharing -- the SAME mechanism this file already
% uses (read-only) for mixo/tolo/maxo/Jxx0/tmf etc.; no new production dependency.
tracing = isfield(opts, 'trace') && ~isempty(opts.trace) && ~isequal(opts.trace, false);
trc = struct('schema_version', 1, 'enabled', tracing, 'meta', struct(), ...
             'nodes', struct([]), 'iters', struct([]));
if tracing
    optsRec = opts;
    if isfield(optsRec, 'trace'), optsRec = rmfield(optsRec, 'trace'); end
    trc.meta = struct('T', T, 'Bx', Bx, 'J0eff', J0eff, 'Jnu_flat', Jnu_flat(:), 'opts', optsRec);
end
node_seq  = 0;    % running eval_node CALL counter (assigns nrec.id; bookkeeping only)
cur_phase = '';   % phase tag for the NEXT eval_node call, set at each call site below

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
% verbatim (wn, wts, beta, eopts from opts -- honors Ecut and EMT options, P1-6).
Ecut  = getf(opts, 'Ecut', 40);
eopts = getf(opts, 'emt', struct());
[wn, wts, beta] = invz_matsubara(T, Ecut);

% Independent h = 0 PM predictor node (round-3 P0-3; doubles as Gate 6b's comparator):
% ONE node solve at hz_fixed = 0 gives THIS machinery's PM fixed point. Its mass
%   slope_pred = r(0) + J0eff*G0bare(0) = 1 + Sigma0(0) - J0eff*chi_path(0)   (= crit, SS5)
% predicts root existence INDEPENDENTLY of any sampled profile value.
Sigma = [];  K0s = 0;                                % warm-start carriers across nodes
if tracing, cur_phase = 'predictor'; end             % stage-2c task 0: node phase tag (bookkeeping only)
[r0n, ~, S0pm, K0pm, ~, Gb0, ~, ok0, Sigma, K0s] = eval_node(0, Sigma, K0s);
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
    slope_pred = r0n + J0eff*Gb0;
    prof.Sigma0_pm0 = S0pm;  prof.K0_pm0 = K0pm;  prof.slope0 = slope_pred;
else
    slope_pred = NaN;
end

ratio = hfrac^(1/(nH-1));
hgrid = hmax * ratio.^((nH-1):-1:0);                 % geometric, clustered at 0 (P1-4)
prof.hmin_initial = hgrid(1);

if tracing, cur_phase = 'sweep'; end                 % stage-2c task 0: node phase tag (bookkeeping only)
[rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s);

if ok0
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end); % first panel seeded with r(0)
else
    h0 = nan(1, nH);                                     % undefined without a real r(0)
end
F  = h0 - J0eff*mv;

% ADAPTIVE lower extension (round-3 P0-3): predictor-driven, NOT self-referential.
% slope_pred < 0 predicts an ordered root; extend geometrically downward until a
% negative F sample appears or the absolute floor is reached.
n_extend = 0;
while ok0 && slope_pred < 0 && all(F >= 0) && hgrid(1) > hmin_abs
    n_extend = n_extend + 1;
    hext = hgrid(1) * ratio.^(3:-1:1);                % three more decades-fraction nodes
    [re, me, S0e, K0e, De, Gbe, Gse] = deal(nan(1, 3));  ce = false(1, 3);
    if tracing, cur_phase = 'extend'; end             % stage-2c task 0: node phase tag (bookkeeping only)
    for k = 1:3
        [re(k), me(k), S0e(k), K0e(k), De(k), Gbe(k), Gse(k), ce(k), Sigma, K0s] = ...
            eval_node(hext(k), Sigma, K0s);
    end
    hgrid = [hext hgrid];  rv = [re rv];  mv = [me mv];  cnv = [ce cnv];
    S0v = [S0e S0v];  K0v = [K0e K0v];  Dv = [De Dv];  Gbv = [Gbe Gbv];  Gsv = [Gse Gsv];
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
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
    ratio2 = hfrac_eff^(1/(nH-1));
    hgrid = hmax * ratio2.^((nH-1):-1:0);
    if tracing, cur_phase = 'redensify'; end          % stage-2c task 0: node phase tag (bookkeeping only)
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s);
    h0 = cumtrapz([0 hgrid], [r0n rv]);  h0 = h0(2:end);
    F  = h0 - J0eff*mv;
    prof.redensified = true;
end
prof.n_extend = n_extend;

prof.hgrid = hgrid;  prof.r = rv;  prof.h0 = h0;  prof.m = mv;
prof.Sigma0 = S0v;   prof.K0 = K0v;  prof.D_uni = Dv;  prof.node_conv = cnv;  prof.F = F;
prof.G0bare = Gbv;   prof.Gstat = Gsv;

if ~ok0                          % predictor never converged: h0/F are undefined (NaN)
    prof.status = 'node_failed'; % above -- report the honest verdict now that the grid's
    return;                      % own (convergence-independent) diagnostics are exported;
end                              % NEVER fall through to the F-based search on NaN data
if slope_pred < 0 && all(F >= 0)                      % floor hit without a bracket:
    prof.status = 'unresolved';                       % NEVER silently PM (round-3 P0-3)
    warning('invz:hmfUnresolved', ...
        'ordering predicted (slope_pred = %.3g) but no negative F above hmin_abs = %.3g', ...
        slope_pred, hmin_abs);
    return;                                           % hmf_star stays NaN; the jensen solver
end                                                   % must return converged = false here
if any(~cnv)                                          % round-4 P1-C: status must be truthful
    prof.status = 'node_failed';                      % on node failure -- never 'ok'
    return;
end
prof.status = 'ok';
s = sign(F);  idx = find(s(1:end-1) < 0 & s(2:end) >= 0, 1, 'last');
if isempty(idx), return; end                          % no nonzero root: PM side

% --- Root refinement by DIRECT evaluation (P1-4): bisection on F between the
% bracketing nodes, fresh node solve per iterate, cumulative h0 via local trapezoid
% panel from the bracket's left node.
a = hgrid(idx);  b = hgrid(idx+1);  Fa = F(idx);  h0a = h0(idx);  ra = rv(idx);
if tracing, cur_phase = 'bisect'; end                % stage-2c task 0: node phase tag (bookkeeping only)
for it = 1:12
    c = 0.5*(a + b);
    [rc, mc, ~, ~, ~, ~, ~, okc, Sigma, K0s] = eval_node(c, Sigma, K0s);
    if ~okc                                           % round-5 P1-A: a failed bisection node
        prof.status = 'node_failed';  hmf_star = NaN; % TERMINATES the solve -- never a root
        return;                                       % from a partial bracket
    end
    h0c = h0a + 0.5*(ra + rc)*(c - a);
    Fc  = h0c - J0eff*mc;
    if sign(Fc) == sign(Fa), a = c; Fa = Fc; h0a = h0c; ra = rc; else, b = c; end
    if (b - a) < trt*b, break; end
end
if (b - a) >= trt*b                                   % round-5 P1-A: tol_root not reached --
    prof.status = 'unresolved';  hmf_star = NaN;      % a distinct refinement failure
    warning('invz:hmfUnresolved', 'root bracket not refined to tol_root: (b-a)/b = %.3g', (b-a)/b);
    return;
end
hmf_star = 0.5*(a + b);
if tracing, cur_phase = 'root'; end                  % stage-2c task 0: node phase tag (bookkeeping only)
[r_s, m_s, ~, ~, D_s, ~, Gs_s, ok_s] = eval_node(hmf_star, Sigma, K0s);
if ~ok_s
    prof.status = 'node_failed';  hmf_star = NaN;  return;
end
prof.m_star = m_s;  prof.D_uni_star = D_s;  prof.r_star = r_s;  prof.Gstat_star = Gs_s;

    function [rv, mv, S0v, K0v, Dv, Gbv, Gsv, cnv, Sigma, K0s] = run_sweep(hgrid, Sigma, K0s)
    % Behavior-neutral factoring of the per-node nH sweep (execution amendment 3):
    % evaluate eval_node at every point of hgrid in ascending order, warm-starting
    % Sigma/K0s across calls. Shared by the initial profile sweep and the re-densify
    % pass; the extension's 3-node prepends stay inline (unchanged).
    n = numel(hgrid);
    [rv, mv, S0v, K0v, Dv, Gbv, Gsv] = deal(nan(1, n));  cnv = false(1, n);
    for is = 1:n
        [rv(is), mv(is), S0v(is), K0v(is), Dv(is), Gbv(is), Gsv(is), cnv(is), Sigma, K0s] = ...
            eval_node(hgrid(is), Sigma, K0s);
    end
    end

    function [rk, mk, S0k, K0k, Dk, Gbk, Gsk, ok, Sigma, K0s] = eval_node(hp, Sigma, K0s)
    % One fixed-field node: si (hz_fixed, NO order), tl, c0/G0, then the ordered
    % Sigma<->EMT loop WITH the static-sector closure each pass (Interfaces bullet).
    if tracing                                     % stage-2c task 0: node identity + seed
        node_seq = node_seq + 1;                   % provenance, captured BEFORE Sigma is
        nrec = struct('id', node_seq, 'h', hp, 'phase', cur_phase, ...             % touched
            'seed_kind', 'warm', 'seed_from', node_seq - 1, ...    % single sequential thread
            'outer_iters', NaN, 'outer_hit_max', false, 'dS_break', false, ...
            'ok_final', false, 'term_reason', '', 'K0', NaN, 'D_uni', NaN, ...
            'resid_static', NaN);
        if isempty(Sigma)                          % this file's OWN cold-start criterion
            nrec.seed_kind = 'cold';  nrec.seed_from = NaN;
        end
    end
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bx, sio);
    if abs(si.hz - hp) > 1e-12
        error('invz:hzFixed', 'hz_fixed not held: si.hz = %.6g vs %.6g', si.hz, hp);
    end
    tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
    mk = si.Jexp(3);
    if fbare
        rk = 1;  S0k = 0;  K0k = 0;  Dk = NaN;  Gbk = NaN;  Gsk = NaN;  ok = true;
        if tracing                                 % stage-2c task 0: degenerate node record
            nrec.outer_iters = 0;  nrec.outer_hit_max = false;  nrec.dS_break = false;
            nrec.ok_final = true;  nrec.term_reason = 'bare_shortcut';
            if isempty(trc.nodes), trc.nodes = nrec; else, trc.nodes(end+1) = nrec; end
        end
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
        'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu_flat);
    if isempty(Sigma_in)                 % this file's OWN cold-start criterion (unchanged)
        seed = [];
    else
        seed = struct('Sigma', Sigma_in, 'K0s', K0s_in);
    end
    sopts = struct('mix_outer', mixo, 'max_outer', maxo, 'tol_outer', tolo, ...
        'cold_retry', true, 'trace', tracing);
    [state, info] = invz_ordered_node_solve(node, seed, sopts);

    % checker-gated acceptance: info.accepted is the COMPLETE invz_ordered_residual verdict,
    % replacing the old in-loop `dS<tolo && sout.converged` + post-loop `so.converged &&
    % so.resid<ctol` gate -- that logic now lives entirely in the solver/checker.
    rk = info.so.r;  S0k = state.Sigma(1);  K0k = state.K0s;  Dk = info.so.D_uni;
    Gbk = G0bare0;   Gsk = info.so.Gstat;   ok = info.accepted;

    % P0-3 seed-safety: thread the ACCEPTED state forward; a non-accepted node instead
    % returns the INPUT Sigma/K0s UNCHANGED (rollback -- a cold input, i.e. Sigma_in = [],
    % rolls back to cold, matching the pre-existing cold criterion for the next node).
    if ok
        Sigma = state.Sigma;  K0s = state.K0s;
    else
        Sigma = Sigma_in;  K0s = K0s_in;
    end

    if tracing                                     % stage-2c task 1b-ii-A: node record
        nrec.outer_iters   = info.outer_iters;      % finalize -- relocated to source info/state
        nrec.dS_break      = info.loop_converged;
        nrec.outer_hit_max = ~info.loop_converged && info.outer_iters == maxo;
        nrec.ok_final      = info.accepted;
        switch info.term_reason                    % solver vocabulary -> this file's own
            case 'accepted',                    nrec.term_reason = 'converged';
            case 'loop_converged_not_accepted', nrec.term_reason = 'refresh_failed';
            case 'max_iter',                    nrec.term_reason = 'max_iter';
            otherwise,                          nrec.term_reason = info.term_reason;
        end
        nrec.K0 = state.K0s;  nrec.D_uni = info.so.D_uni;  nrec.resid_static = info.so.resid;
        if isempty(trc.nodes), trc.nodes = nrec; else, trc.nodes(end+1) = nrec; end

        for ii = 1:numel(info.iters)                % per-outer-iteration records, relocated
            src = info.iters(ii);                   % from info.iters (the solver's own trace)
            irec = struct('node_id', nrec.id, 'outer', src.outer, 'resid_map', src.dS, ...
                'resid_static', src.resid_static, 'K0', src.K0, 'D_uni', src.D_uni, ...
                'Dq_min', src.Dq_min, 'Dq_max', src.Dq_max, 'Dq_neg_count', src.Dq_neg_count, ...
                'idx_pos_flat', src.idx_pos_flat, 'Dq_pos_val', src.Dq_pos_val, ...
                'idx_neg_flat', src.idx_neg_flat, 'Dq_neg_val', src.Dq_neg_val, ...
                'converged_flag', src.converged_flag);
            if isempty(trc.iters), trc.iters = irec; else, trc.iters(end+1) = irec; end
        end
    end
    end
end
