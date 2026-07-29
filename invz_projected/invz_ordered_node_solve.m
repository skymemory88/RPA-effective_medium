function [state, info] = invz_ordered_node_solve(node, seed, sopts)
%INVZ_ORDERED_NODE_SOLVE Shared, pure, seed-safe ordered "jensen" node solver (stage-2c,
% task 1b-i). Consolidates the two duplicated per-node damped-Picard loops --
% invz_solve_point_ordered.m:206-235 (jensen closure loop + post-loop static refresh) and
% invz_hmf_ordered.m's eval_node (:313-368, excluding its fbare shortcut -- see the note
% below) -- into ONE function, run VERBATIM (same helper calls, same argument order, no
% physics change), gating node ACCEPTANCE on the Task-1a complete residual checker
% (invz_ordered_residual.m; contract summarized in invzp_convg_diagnosis.md Section 2.4)
% instead of the
% current incomplete in-loop `dS < tol_outer && sout.converged` alone (that verdict is still
% computed and exposed as info.loop_converged, a diagnostic -- see that contract summary for
% why it under-certifies a fixed point: the dynamic K(2:end) is never recomputed from the
% final mixed Sigma).
%
% THIS IS AN ADDITIVE, STANDALONE, PURE FUNCTION (stage-2c task 1b-i scope fence): neither
% invz_hmf_ordered.m nor invz_solve_point_ordered.m is modified or rewired to call this --
% that is task 1b-ii. This file and its test are the entire diff.
%
% node  (struct, REQUIRED fields, the SAME 13-field bundle invz_ordered_residual consumes,
%        constant across one node's outer loop):
%   tl, G0, g, wts, wn, beta, J0eff, G0inel0, G0el0, G0bare0, eso, eopts, Jnu_flat
%   (read-only: never mutated here; G0bare0 is accepted for provenance/parity with the
%   checker's own node schema and is not read by the map below, exactly as in
%   invz_ordered_residual.m.)
%   Jmom (14th field, invz_coupling_moments of the static Jnu_flat column) is REQUIRED once
%   node.eso.static_medium selects a strict_1z_* scheme (spec SS4.3); harmless-absent under
%   'resummed'. Threaded verbatim into both EMT leaves (eopts.Jmom for the dynamic slot,
%   eso_local.Jmom for the static slot) so neither leaf silently re-derives it.
% seed  (struct('Sigma', size(node.wn), 'K0s', scalar), OPTIONAL): the warm-start
%        continuation state -- the SAME two quantities BOTH current loops actually thread
%        across node calls (confirmed directly at invz_hmf_ordered.m's eval_node call sites,
%        run_sweep/the extension loop/the bisection loop, all of which pass Sigma/K0s through
%        as nested-function in/out arguments). EMPTY ([] or struct() with no Sigma/K0s
%        fields) => COLD start: Sigma = zeros(size(node.wn)), K0s = 0. `lam` is ALWAYS reset
%        to [0;0;0] at the start of the loop, WARM OR COLD -- this is not a new convention
%        invented here, it is what BOTH existing loops actually do: invz_hmf_ordered.m's
%        eval_node resets `lam = [0;0;0]` unconditionally at the top of every call
%        (:312), even when Sigma/K0s are warm-threaded from the previous node/profile point;
%        invz_solve_point_ordered.m only ever cold-starts its single per-point solve
%        (:205), so it is trivially consistent with the same rule. seed/node are never
%        mutated (read-only throughout; only local copies are used).
% sopts (struct, OPTIONAL, solver knobs, all defaulting to current production values):
%   .mix_outer  (default 0.7,   = the current mixo)
%   .max_outer  (default 200,   = the current maxo)
%   .tol_outer  (default 1e-8,  = the current tolo; also threaded into the checker's own
%               opts.tol_outer, so the acceptance scales and the loop's in-loop gate never
%               disagree about which tolerance is "the" outer tolerance)
%   .cold_retry (default true): see the retry semantics below.
%   .trace      (default false): gates info.iters (see below); default path stays cheap.
%   .cold_acceleration (default 'none'): 'none' | 'signed_aitken1'. The opt-in pilot is
%               applied ONLY to a genuinely cold attempt under the resummed static medium.
%               It fits one signed scalar factor to four successive full Sigma-vector
%               increments (three ratios), requires a stable negative oscillatory factor,
%               a common static coupling interval, closed inner solves, and a small
%               one-mode fit error. A proposal is retained only when fresh unmixed
%               full-Sigma residual evaluations at the ordinary and accelerated candidates
%               stay in that interval and the accelerated residual is strictly smaller.
%               Each accepted/rejected proposal restarts the history, so another proposal
%               must independently re-earn every gate; ordinary Picard iteration runs
%               between proposals. Final acceptance remains the unchanged A--D checker.
%
% state (struct('Sigma','K','lam','K0s'), FIELD-COMPATIBLE with invz_ordered_residual's
%        `state` input): the exported tuple of the WINNING attempt (see retry semantics
%        below) -- Sigma/K on the node.wn grid, lam = invz_lambdas(...,[1 2 3]), K0s the
%        scalar continuation seed (also written into state.K(1), exactly as both current
%        loops do).
%
% info (struct, ALL fields ALWAYS present; NaN/empty when not applicable):
%   .accepted       (logical) = the winning attempt's res.accepted.
%   .res            the full invz_ordered_residual(...) result for the RETURNED state
%                   (blockA-D, D_uni, Dq_min/Dq_max, finite, stall).
%   .loop_converged (logical) the WINNING attempt's in-loop `dS<tol_outer && sout.converged`
%                   verdict -- diagnostic only, NEVER the acceptance criterion.
%   .so             the WINNING attempt's final (post-loop-refresh) static closure struct
%                   (.r, .D_uni, .resid, .converged, plus .Gstat merged in -- Gstat is
%                   invz_emt_static_ordered's sibling 2nd output, not a field of its own 3rd
%                   output struct, so it is folded in here for a single self-contained bundle).
%   .med            the WINNING attempt's final dynamic medium (the LAST outer iteration's
%                   invz_emt_scalar return, UNTOUCHED: med.G(1) is left as the ordinary-Dyson
%                   value here -- the CALLER overwrites med.G(1) = Gstat at pack time in
%                   1b-ii, matching invz_solve_point_ordered.m:235; this solver does not
%                   perform that substitution itself).
%   .outer_iters    (integer) the WINNING attempt's outer-loop iteration count.
%   .seed_kind      'cold' | 'warm' | 'cold_after_warm_fail'.
%   .term_reason    'accepted' | 'medium_out_of_domain' | 'loop_converged_not_accepted' |
%                   'max_iter' -- the WINNING attempt's own classification
%                   ('medium_out_of_domain': a strict-scheme reference/closure domain event
%                   stopped the attempt before invz_lambdas/invz_sigma_ordered ran, spec SS4.4;
%                   'bare_shortcut' is reserved for 1b-ii's caller-level fbare branch, see the
%                   note below; this solver never emits it).
%   .iters          struct array, one record per OUTER iteration of the WINNING attempt,
%                   populated ONLY when sopts.trace: fields `outer`, `dS` (raw pre-mix Sigma-
%                   map residual), `resid_static` (sout.resid), `K0` (=K0s), `D_uni`
%                   (sout.D_uni), `converged_flag` (sout.converged), `Dq_min`, `Dq_max`,
%                   `Dq_abs_min`, `Dq_neg_count`, `idx_pos_flat`, `Dq_pos_val`, `idx_neg_flat`,
%                   `Dq_neg_val`, dynamic `x`/denominator proximity, static `y`/coupling
%                   interval, Gstat/local denominator/xi, and `static_closed_outer_open`
%                   (closest-to-zero positive/negative Dq = 1+(Jnu_flat-K0s).*Gstat, flat
%                   indices + values) -- mirroring invz_hmf_ordered.m eval_node's own per-
%                   iteration trace record (:329-349) field-for-field, except `dS` (named
%                   `resid_map` there) which uses the name THIS brief specifies. Empty (0x0,
%                   no fields) when sopts.trace is false, so the default path pays nothing
%                   for it. NOT wired to invz_hmf_ordered's trc schema here -- that
%                   relocation is 1b-ii's job.
%   .medium         the WINNING attempt's final medium struct (.scheme/.status/.ref/.closure),
%                   from the LAST invz_emt_static_ordered call that actually ran (the post-loop
%                   refresh when it ran; the last in-loop call when a domain event pre-empted
%                   it). status='not_applicable'/scheme='resummed' placeholder under the legacy
%                   scheme (task 8's convention).
%   .medium_status  the same call's out.medium_status ('not_applicable' | 'ok' |
%                   'ref_denom_nonpositive' | 'ref_denom_small' | 'nonfinite'), exposed
%                   top-level so a caller does not need to reach into .medium.status.
%   .acceleration   fixed-schema summary of the cold accelerator: mode/enabled, proposal
%                   attempt/accept counts, accepted outer iterations and latest signed
%                   factor, fresh ordinary/proposal residuals, interval rank, mode-fit
%                   error, and last rejection reason. This is diagnostic only.
%
% Retry semantics (deterministic; point 3 of the brief): a WARM-seeded attempt that is not
% accepted, with sopts.cold_retry true, is retried EXACTLY ONCE from a cold start (fresh
% Sigma = zeros/K0s = 0/lam = [0;0;0] -- the failed warm attempt's state never leaks in).
% Whichever cold-retry result comes back (accepted or not) is what is returned, with
% info.seed_kind = 'cold_after_warm_fail'. If the ORIGINAL seed was already cold, or
% sopts.cold_retry is false, no retry is attempted: retrying an already-cold start against
% itself is a deterministic no-op (identical inputs, identical map => identical result), so
% it can never change the outcome and is not attempted.
%
% No exception absorber (task 10 removed the previous per-attempt try/catch): ordinary
% numerical non-convergence never throws anywhere in the call chain this map exercises
% (invz_emt_scalar / invz_emt_static_ordered / invz_gstat_ordered / invz_lambdas /
% invz_sigma_ordered / invz_sigma) -- it always returns non-finite/non-converged outputs,
% which the checker below then correctly scores as not accepted. Under a strict scheme, a
% reference/closure domain event is likewise never an exception: it is a returned
% medium_status, and this map halts the attempt at that point (before invz_lambdas /
% invz_sigma_ordered can consume the invalid reference, spec SS4.4), scored
% term_reason='medium_out_of_domain' above. Any invz:* identifier that IS thrown and reaches
% this function -- e.g. invz_emt_scalar's 'invz:emtJnu' (a malformed [nJ,nw] Jnu_flat matrix
% shape) or invz_static_medium_reference's 'invz:staticMedium' (an unrecognized scheme) -- is
% therefore always a wiring/programming defect, never a numerical one, and escapes to the
% caller unchanged, exactly like any non-invz:* exception.
%
% Note on `fbare` (invz_hmf_ordered.m eval_node's force-bare shortcut: rk=1, ok=true, NO
% loop at all): that is a CALLER-level branch on the fixed-field h, not part of the per-node
% outer map itself, so it is intentionally NOT implemented here. Task 1b-ii keeps `fbare`
% handling in the caller and simply does not invoke this solver on a force-bare node.
%
% Purity: node/seed are read-only throughout (only local copies are used, and nothing is
% ever assigned back into a field of either); no persistent/global state; deterministic
% given identical inputs.
if nargin < 2, seed = []; end
if nargin < 3 || isempty(sopts), sopts = struct(); end

req_node = {'tl', 'G0', 'g', 'wts', 'wn', 'beta', 'J0eff', 'G0inel0', 'G0el0', 'G0bare0', ...
            'eso', 'eopts', 'Jnu_flat'};
for k = 1:numel(req_node)
    if ~isfield(node, req_node{k})
        error('invz:nodeSolveNode', 'node is missing required field ''%s''.', req_node{k});
    end
end

% Jmom is required only once a STRICT scheme is selected (spec SS4.3).
smid_node = getf(node.eso, 'static_medium', 'resummed');
if ~strcmp(smid_node, 'resummed') && (~isfield(node, 'Jmom') || isempty(node.Jmom))
    error('invz:nodeSolveNode', ['node.Jmom is required under static_medium ''%s'' ' ...
        '(invz_coupling_moments of the static coupling column).'], smid_node);
end

cold_acceleration = getf(sopts, 'cold_acceleration', 'none');
if ~(ischar(cold_acceleration) && isrow(cold_acceleration)) || ...
        ~any(strcmp(cold_acceleration, {'none', 'signed_aitken1'}))
    error('invz:coldAcceleration', ...
        'cold_acceleration must be ''none'' or ''signed_aitken1''.');
end
sopts = struct('mix_outer',       getf(sopts, 'mix_outer', 0.7), ...
               'max_outer',       getf(sopts, 'max_outer', 200), ...
               'tol_outer',       getf(sopts, 'tol_outer', 1e-8), ...
               'cold_retry',      getf(sopts, 'cold_retry', true), ...
               'trace',           getf(sopts, 'trace', false), ...
               'cold_acceleration', cold_acceleration);
accel_requested = strcmp(sopts.cold_acceleration, 'signed_aitken1');
if accel_requested && ~strcmp(smid_node, 'resummed')
    error('invz:coldAcceleration', ...
        'signed_aitken1 is currently restricted to static_medium ''resummed''.');
end
if sopts.trace || accel_requested
    Jflat = node.Jnu_flat(:);
    jdiag = struct('flat', Jflat, 'max', max(Jflat), 'sorted', sort(Jflat));
else
    jdiag = struct([]);
end

is_warm = isstruct(seed) && isscalar(seed) && isfield(seed, 'Sigma') && isfield(seed, 'K0s') ...
          && ~isempty(seed.Sigma);
if is_warm
    Sigma0 = seed.Sigma(:);  K0s0 = seed.K0s;  seed_kind0 = 'warm';
else
    Sigma0 = zeros(size(node.wn));  K0s0 = 0;  seed_kind0 = 'cold';
end

[state1, info1] = run_attempt(node, Sigma0, K0s0, sopts, jdiag, ...
    accel_requested && strcmp(seed_kind0, 'cold'));

if info1.res.accepted || ~(sopts.cold_retry && strcmp(seed_kind0, 'warm'))
    state = state1;  info = info1;  seed_kind = seed_kind0;
else
    % deterministic cold retry: fresh cold start, the failed warm attempt is discarded
    % entirely (not merged/blended) -- point 3 of the brief.
    [state2, info2] = run_attempt(node, zeros(size(node.wn)), 0, sopts, jdiag, ...
        accel_requested);
    state = state2;  info = info2;  seed_kind = 'cold_after_warm_fail';
end
info.seed_kind = seed_kind;
info.accepted  = info.res.accepted;
end

% =============================================================================================
function [state, info] = run_attempt(node, Sigma0, K0s0, sopts, jdiag, accel_enabled)
%RUN_ATTEMPT One full damped-Picard attempt: the outer Sigma<->EMT loop
% (invz_solve_point_ordered.m:206-221 / invz_hmf_ordered.m eval_node:313-328, VERBATIM) plus
% the post-loop static refresh (invz_solve_point_ordered.m:225-227, VERBATIM), from a given
% (Sigma0, K0s0) start (lam ALWAYS [0;0;0] -- see the file header). No exception absorber: an
% invz:* identifier reaching here is a wiring defect and escapes; ordinary non-convergence and
% strict-scheme domain events are never exceptions (see the file header). A domain event halts
% the attempt early and `state` reflects the last fully-computed values at that point. Builds
% `state`, then scores it with the complete checker (invz_ordered_residual).
wn = node.wn(:);  G0 = node.G0(:);  g = node.g(:);  wts = node.wts(:);
Jnu_flat = node.Jnu_flat;  if isvector(Jnu_flat), Jnu_flat = Jnu_flat(:); end

Sigma = Sigma0(:);
K     = zeros(size(wn));           % matches both loops' per-call reset; inert in practice
                                    % (med.K overwrites it on iteration 1 regardless --
                                    % invz_emt_scalar's own docstring: K0 is accepted for
                                    % backward compatibility but unused, the solve is direct).
K0s   = K0s0;
lam   = [0; 0; 0];
eopts = node.eopts;                 % LOCAL copy: .K0 is mutated below; node's own untouched.
if isfield(node, 'Jmom') && ~isempty(node.Jmom)
    eopts.Jmom = node.Jmom;                 % strict slot reads it; legacy path ignores it
    eso_local  = node.eso;  eso_local.Jmom = node.Jmom;
else
    eso_local  = node.eso;
end
medium_status = 'not_applicable';  medium = struct('scheme', getf(node.eso, ...
    'static_medium', 'resummed'), 'status', 'not_applicable', 'ref', [], 'closure', []);

dS = NaN;  loop_converged = false;  outer_used = 0;  iters = struct([]);
med = struct('G', nan(size(wn)), 'K', nan(size(wn)), 'converged', false, 'closure', NaN, ...
             'iters', 0, 'dynamic_converged', false);
so    = struct('D_uni', NaN, 'resid', NaN, 'converged', false, 'iters', 0);
Gstat = NaN;
acceleration = blank_acceleration_summary(sopts.cold_acceleration, accel_enabled);
accel_history = struct('deltas', {{}}, 'interval_rank', [], ...
                       'inner_closed', [], 'pole_clear', []);

for outer = 1:sopts.max_outer
    % (1) dynamic sector -- MIRRORS both loops' emt call verbatim
    eopts.K0 = K;
    med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
    K   = med.K;
    % (2) static sector (threaded node.eso opts, verbatim argument order):
    [K0s, Gstat_it, sout] = invz_emt_static_ordered(node.tl, lam(1:2), Sigma(1), ...
        Jnu_flat, K0s, node.beta, node.J0eff, node.G0inel0, node.G0el0, eso_local);
    % Domain event: stop BEFORE invz_lambdas / invz_sigma_ordered consume an invalid
    % reference (spec SS4.4). Propagating a non-finite K0 through the outer map would turn a
    % reportable domain outcome into an unexplained non-convergence.
    medium_status = getf(sout, 'medium_status', 'not_applicable');
    medium        = getf(sout, 'medium', medium);
    if ~any(strcmp(medium_status, {'not_applicable', 'ok'}))
        if sopts.trace
            iters = append_iter(iters, outer, NaN, sout, K0s, Gstat_it, ...
                Sigma(1), G0(1), jdiag, sopts.tol_outer, blank_acceleration_iter());
        end
        outer_used = outer;
        break;
    end
    K(1) = K0s;
    % (3)-(5) lambdas, ordered Sigma map, pre-mix step -- MIRRORS both loops verbatim
    lam = invz_lambdas(K, g, wts, node.beta, [1 2 3]);
    sg  = invz_sigma_ordered(node.tl, lam, K, g, node.beta);
    dS  = invz_finite_max_abs(sg.Sigma, Sigma);
    % (6) damped mix
    if accel_enabled
        Sigma_before = Sigma;
        Sigma_picard = Sigma_before + sopts.mix_outer*(sg.Sigma - Sigma_before);
        [Sigma, accel_history, acceleration, accel_iter] = try_cold_acceleration( ...
            node, Sigma_before, Sigma_picard, dS, sout, K0s, Gstat_it, lam, K, ...
            eopts, eso_local, jdiag, outer, accel_history, acceleration);
        if sopts.trace
            iters = append_iter(iters, outer, dS, sout, K0s, Gstat_it, ...
                Sigma_before(1), G0(1), jdiag, sopts.tol_outer, accel_iter);
        end
    else
        if sopts.trace
            iters = append_iter(iters, outer, dS, sout, K0s, Gstat_it, ...
                Sigma(1), G0(1), jdiag, sopts.tol_outer, blank_acceleration_iter());
        end
        % Keep the default-off arithmetic statement byte-for-byte identical to the
        % pre-experiment path (I1).
        Sigma = Sigma + sopts.mix_outer*(sg.Sigma - Sigma);
    end
    outer_used = outer;
    % (7) in-loop verdict -- DIAGNOSTIC ONLY (info.loop_converged); NEVER the acceptance
    % test (that is the incomplete gate this whole task replaces).
    if dS < sopts.tol_outer && sout.converged
        loop_converged = true;
        break;
    end
end
% post-loop static refresh (invz_solve_point_ordered.m:225-227, verbatim): recomputed on
% the just-mixed Sigma(1), continuing from the loop's own last K0s. Skipped after a domain
% event -- that would feed the invalid reference into the static function a second time and
% contradict "stop before consumption" (spec SS4.4).
if any(strcmp(medium_status, {'not_applicable', 'ok'}))
    [K0s, Gstat, so] = invz_emt_static_ordered(node.tl, lam(1:2), Sigma(1), Jnu_flat, K0s, ...
        node.beta, node.J0eff, node.G0inel0, node.G0el0, eso_local);
    medium_status = getf(so, 'medium_status', medium_status);
    medium        = getf(so, 'medium', medium);
    if any(strcmp(medium_status, {'not_applicable', 'ok'})), K(1) = K0s; end
else
    so = sout;  Gstat = Gstat_it;       % preserve the classified domain record
end

state = struct('Sigma', Sigma, 'K', K, 'lam', lam, 'K0s', K0s);

% ---- acceptance = the complete checker, NOT the bare in-loop step (point 2 of the brief) --
res = invz_ordered_residual(node, state, struct('tol_outer', sopts.tol_outer, 'dS', dS));

if res.accepted
    term_reason = 'accepted';
elseif ~any(strcmp(medium_status, {'not_applicable', 'ok'}))
    term_reason = 'medium_out_of_domain';   % distinct from a convergence failure
elseif loop_converged
    term_reason = 'loop_converged_not_accepted';
else
    term_reason = 'max_iter';
end

so_out = so;  so_out.Gstat = Gstat;

info = struct('res', res, 'loop_converged', loop_converged, 'so', so_out, 'med', med, ...
              'outer_iters', outer_used, 'term_reason', term_reason, 'iters', iters, ...
              'medium', medium, 'medium_status', medium_status, ...
              'acceleration', acceleration);
end

% =============================================================================================
function [Sigma_out, hist, summary, diag] = try_cold_acceleration( ...
        node, Sigma_old, Sigma_picard, dS, sout, K0s, Gstat_it, lam, K, ...
        eopts, eso_local, jdiag, outer, hist, summary)
%TRY_COLD_ACCELERATION One default-off, residual-decreasing signed Aitken-1 proposal.
%
% The thresholds below are fixed experiment safeguards, not production tuning:
%   * four full Sigma-vector increments / three signed ratios;
%   * at least eight ordinary outer iterations;
%   * lambda in [-0.99,-0.50], ratio spread <= 0.02;
%   * pooled scalar-mode fit error <= 0.10;
%   * one common finite pole interval and closed inner solve over the window;
%   * fresh proposal residual strictly below both the ordinary-candidate residual and
%     the current residual;
%   * a fresh four-increment history after every accepted or rejected proposal.
policy = struct('history', 4, 'min_outer', 8, ...
                'lambda_min', -0.99, 'lambda_max', -0.50, ...
                'lambda_spread', 0.02, 'mode_fit_rel', 0.10);

Sigma_out = Sigma_picard;
diag = blank_acceleration_iter();
delta = Sigma_picard - Sigma_old;
[interval_rank, ~, ~] = coupling_interval(jdiag.sorted, K0s - 1/Gstat_it);
Dq = 1 + (jdiag.flat - K0s).*Gstat_it;
pole_clear = all(isfinite(Dq)) && ~isempty(Dq) && min(abs(Dq)) > 0;
inner_closed = sout.converged && isfinite(sout.resid);

hist.deltas{end+1} = delta;
hist.interval_rank(end+1) = interval_rank;
hist.inner_closed(end+1) = inner_closed;
hist.pole_clear(end+1) = pole_clear;
if numel(hist.deltas) > policy.history
    hist.deltas = hist.deltas(end-policy.history+1:end);
    hist.interval_rank = hist.interval_rank(end-policy.history+1:end);
    hist.inner_closed = hist.inner_closed(end-policy.history+1:end);
    hist.pole_clear = hist.pole_clear(end-policy.history+1:end);
end

if numel(hist.deltas) < policy.history || outer < policy.min_outer
    return;
end

[lambda, ratio_spread, fit_rel, fit_ok] = signed_increment_factor(hist.deltas);
diag.lambda = lambda;
diag.lambda_spread = ratio_spread;
diag.mode_fit_rel = fit_rel;
if ~fit_ok || ~(lambda >= policy.lambda_min && lambda <= policy.lambda_max) || ...
        ratio_spread > policy.lambda_spread || fit_rel > policy.mode_fit_rel
    return;
end
target_rank = hist.interval_rank(end);
if ~isfinite(target_rank) || any(hist.interval_rank ~= target_rank) || ...
        ~all(hist.inner_closed) || ~all(hist.pole_clear)
    return;
end

coeff = lambda/(1-lambda);
Sigma_proposal = Sigma_picard + coeff*delta;
diag.attempted = true;
diag.interval_rank = target_rank;
summary.attempted = summary.attempted + 1;

ordinary = fresh_outer_residual(node, Sigma_picard, K0s, lam, K, ...
    eopts, eso_local, jdiag);
proposal = fresh_outer_residual(node, Sigma_proposal, K0s, lam, K, ...
    eopts, eso_local, jdiag);
diag.resid_ordinary = ordinary.resid;
diag.resid_proposal = proposal.resid;

same_interval = ordinary.interval_rank == target_rank && ...
                proposal.interval_rank == target_rank;
residual_decrease = isfinite(proposal.resid) && isfinite(ordinary.resid) && ...
                    isfinite(dS) && proposal.resid < ordinary.resid && ...
                    proposal.resid < dS;
if ordinary.valid && proposal.valid && ordinary.inner_closed && ...
        proposal.inner_closed && ordinary.pole_clear && proposal.pole_clear && ...
        same_interval && residual_decrease
    Sigma_out = Sigma_proposal;
    diag.accepted = true;
    diag.reject_reason = '';
    summary.accepted = summary.accepted + 1;
    summary.accepted_outers(end+1) = outer;
    summary.lambda = lambda;
    summary.lambda_spread = ratio_spread;
    summary.mode_fit_rel = fit_rel;
    summary.resid_ordinary = ordinary.resid;
    summary.resid_proposal = proposal.resid;
    summary.interval_rank = target_rank;
    summary.last_reject_reason = '';
else
    if ~(ordinary.valid && proposal.valid)
        reason = 'nonfinite_candidate';
    elseif ~(ordinary.inner_closed && proposal.inner_closed)
        reason = 'inner_not_closed';
    elseif ~(ordinary.pole_clear && proposal.pole_clear)
        reason = 'pole_not_clear';
    elseif ~same_interval
        reason = 'interval_changed';
    else
        reason = 'residual_not_decreased';
    end
    diag.reject_reason = reason;
    summary.last_reject_reason = reason;
end

% A rejected proposal must earn a fresh four-increment window before another attempt.
% An accepted proposal is also a restart, so ordinary Picard must establish a new local
% history before any later proposal. Final A--D acceptance remains independent.
hist = struct('deltas', {{}}, 'interval_rank', [], ...
              'inner_closed', [], 'pole_clear', []);
end

function [lambda, spread, fit_rel, ok] = signed_increment_factor(deltas)
%SIGNED_INCREMENT_FACTOR Pooled signed scalar fit d_{k+1} ~= lambda*d_k.
lambda = NaN;  spread = NaN;  fit_rel = NaN;  ok = false;
n = numel(deltas);
if n < 2, return; end
num = 0;  den = 0;  ratios = nan(1, n-1);
for k = 2:n
    a = deltas{k-1}(:);  b = deltas{k}(:);
    da = real(a'*a);
    if ~(all(isfinite(a)) && all(isfinite(b)) && isfinite(da) && da > 0)
        return;
    end
    cross = real(a'*b);
    ratios(k-1) = cross/da;
    num = num + cross;
    den = den + da;
end
if ~(isfinite(den) && den > 0 && all(isfinite(ratios))), return; end
lambda = num/den;
err2 = 0;  out2 = 0;
for k = 2:n
    a = deltas{k-1}(:);  b = deltas{k}(:);
    err2 = err2 + real((b-lambda*a)'*(b-lambda*a));
    out2 = out2 + real(b'*b);
end
if ~(isfinite(lambda) && isfinite(err2) && isfinite(out2) && out2 > 0), return; end
spread = max(abs(ratios-lambda));
fit_rel = sqrt(max(err2,0)/out2);
ok = isfinite(spread) && isfinite(fit_rel);
end

function out = fresh_outer_residual(node, Sigma, K0_seed, lam_seed, K_seed, ...
                                    eopts, eso_local, jdiag)
%FRESH_OUTER_RESIDUAL Unmixed full-Sigma residual at one candidate state.
out = struct('valid', false, 'inner_closed', false, 'pole_clear', false, ...
             'resid', NaN, 'interval_rank', NaN);
eopts_local = eopts;
eopts_local.K0 = K_seed;
med = invz_emt_scalar(node.G0(:), Sigma(:), node.Jnu_flat, eopts_local);
Kc = med.K;
[K0c, Gstatc, soc] = invz_emt_static_ordered(node.tl, lam_seed(1:2), Sigma(1), ...
    node.Jnu_flat, K0_seed, node.beta, node.J0eff, node.G0inel0, node.G0el0, eso_local);
medium_status = getf(soc, 'medium_status', 'not_applicable');
if ~any(strcmp(medium_status, {'not_applicable', 'ok'})), return; end
Kc(1) = K0c;
lamc = invz_lambdas(Kc, node.g(:), node.wts(:), node.beta, [1 2 3]);
sgc = invz_sigma_ordered(node.tl, lamc, Kc, node.g(:), node.beta);
out.resid = invz_finite_max_abs(sgc.Sigma, Sigma(:));
[out.interval_rank, ~, ~] = coupling_interval(jdiag.sorted, K0c - 1/Gstatc);
Dq = 1 + (jdiag.flat-K0c).*Gstatc;
out.pole_clear = all(isfinite(Dq)) && ~isempty(Dq) && min(abs(Dq)) > 0;
out.inner_closed = soc.converged && isfinite(soc.resid);
out.valid = all(isfinite(Sigma(:))) && all(isfinite(Kc(:))) && ...
    all(isfinite(lamc(:))) && isfinite(K0c) && isfinite(Gstatc) && ...
    isfinite(out.resid) && isfinite(out.interval_rank);
end

function out = blank_acceleration_summary(mode, enabled)
out = struct('mode', mode, 'enabled', logical(enabled), ...
    'attempted', 0, 'accepted', 0, 'accepted_outers', [], ...
    'lambda', NaN, 'lambda_spread', NaN, 'mode_fit_rel', NaN, ...
    'resid_ordinary', NaN, 'resid_proposal', NaN, ...
    'interval_rank', NaN, 'last_reject_reason', '');
end

function out = blank_acceleration_iter()
out = struct('attempted', false, 'accepted', false, ...
    'lambda', NaN, 'lambda_spread', NaN, 'mode_fit_rel', NaN, ...
    'resid_ordinary', NaN, 'resid_proposal', NaN, ...
    'interval_rank', NaN, 'reject_reason', '');
end

% =============================================================================================
function iters = append_iter(iters, outer, dS, sout, K0s, Gstat_it, ...
                             Sigma0, G00, jdiag, tol_outer, accel)
%APPEND_ITER One per-outer-iteration trace record (sopts.trace only), mirroring
% invz_hmf_ordered.m eval_node's own per-iteration record (:329-349) field-for-field except
% `dS` (named `resid_map` there; `dS` is the name this brief specifies).
Jflat = jdiag.flat;
Dq = 1 + (Jflat - K0s) .* Gstat_it;
Ddyn = 1 + Sigma0 + Jflat.*G00;
x = -(1 + Sigma0)/G00;
y = K0s - 1/Gstat_it;
Dq_abs_min = NaN;
Ddyn_abs_min = NaN;
if all(isfinite(Dq)), Dq_abs_min = min(abs(Dq)); end
if all(isfinite(Ddyn)), Ddyn_abs_min = min(abs(Ddyn)); end
[y_rank, y_lo, y_hi] = coupling_interval(jdiag.sorted, y);
irec = struct('outer', outer, 'dS', dS, 'resid_static', sout.resid, 'K0', K0s, ...
              'D_uni', sout.D_uni, 'converged_flag', sout.converged, ...
              'Dq_min', min(Dq), 'Dq_max', max(Dq), 'Dq_neg_count', nnz(Dq <= 0), ...
              'Dq_abs_min', Dq_abs_min, ...
              'idx_pos_flat', NaN, 'Dq_pos_val', NaN, ...
              'idx_neg_flat', NaN, 'Dq_neg_val', NaN, ...
              'x', x, 'x_minus_Jmax', x-jdiag.max, ...
              'Ddyn_abs_min', Ddyn_abs_min, ...
              'Ddyn_nonpositive_count', nnz(Ddyn <= 0), ...
              'y', y, 'y_rank', y_rank, ...
              'y_interval_lo', y_lo, 'y_interval_hi', y_hi, ...
              'Gstat', Gstat_it, ...
              'gstat_local_denom', getf(sout, 'gstat_local_denom', NaN), ...
              'xi', getf(sout, 'xi', NaN), ...
              'static_closed_outer_open', sout.converged && ...
                  ~(isfinite(dS) && dS < tol_outer), ...
              'accel_attempted', accel.attempted, ...
              'accel_accepted', accel.accepted, ...
              'accel_lambda', accel.lambda, ...
              'accel_lambda_spread', accel.lambda_spread, ...
              'accel_mode_fit_rel', accel.mode_fit_rel, ...
              'accel_resid_ordinary', accel.resid_ordinary, ...
              'accel_resid_proposal', accel.resid_proposal, ...
              'accel_interval_rank', accel.interval_rank, ...
              'accel_reject_reason', accel.reject_reason);
posmask = Dq > 0;
if any(posmask)
    ix = find(posmask);  [vmin, jm] = min(Dq(ix));
    irec.idx_pos_flat = ix(jm);  irec.Dq_pos_val = vmin;
end
negmask = ~posmask;                    % Dq <= 0: the non-positive/unstable side
if any(negmask)
    ix = find(negmask);  [vmax, jm] = max(Dq(ix));   % closest to zero from below
    irec.idx_neg_flat = ix(jm);  irec.Dq_neg_val = vmax;
end
if isempty(iters), iters = irec; else, iters(end+1) = irec; end
end

function [rank_lower, lo, hi] = coupling_interval(sorted_J, y)
rank_lower = NaN;  lo = NaN;  hi = NaN;
if ~isfinite(y) || isempty(sorted_J), return; end
rank_lower = nnz(sorted_J < y);
if rank_lower == 0
    lo = -Inf;  hi = sorted_J(1);
elseif rank_lower == numel(sorted_J)
    lo = sorted_J(end);  hi = Inf;
else
    lo = sorted_J(rank_lower);  hi = sorted_J(rank_lower + 1);
end
end
