function [state, info] = invz_ordered_node_solve(node, seed, sopts)
%INVZ_ORDERED_NODE_SOLVE Shared, pure, seed-safe ordered "jensen" node solver (stage-2c,
% task 1b-i). Consolidates the two duplicated per-node damped-Picard loops --
% invz_solve_point_ordered.m:206-235 (jensen closure loop + post-loop static refresh) and
% invz_hmf_ordered.m's eval_node (:313-368, excluding its fbare shortcut -- see the note
% below) -- into ONE function, run VERBATIM (same helper calls, same argument order, no
% physics change), gating node ACCEPTANCE on the Task-1a complete residual checker
% (invz_ordered_residual.m; contract: docs/invz_ordered_residual_contract.md) instead of the
% current incomplete in-loop `dS < tol_outer && sout.converged` alone (that verdict is still
% computed and exposed as info.loop_converged, a diagnostic -- see the contract Sec 0/2 for
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
%                   `Dq_neg_count`, `idx_pos_flat`, `Dq_pos_val`, `idx_neg_flat`, `Dq_neg_val`
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

sopts = struct('mix_outer',  getf(sopts, 'mix_outer', 0.7), ...
               'max_outer',  getf(sopts, 'max_outer', 200), ...
               'tol_outer',  getf(sopts, 'tol_outer', 1e-8), ...
               'cold_retry', getf(sopts, 'cold_retry', true), ...
               'trace',      getf(sopts, 'trace', false));

is_warm = isstruct(seed) && isscalar(seed) && isfield(seed, 'Sigma') && isfield(seed, 'K0s') ...
          && ~isempty(seed.Sigma);
if is_warm
    Sigma0 = seed.Sigma(:);  K0s0 = seed.K0s;  seed_kind0 = 'warm';
else
    Sigma0 = zeros(size(node.wn));  K0s0 = 0;  seed_kind0 = 'cold';
end

[state1, info1] = run_attempt(node, Sigma0, K0s0, sopts);

if info1.res.accepted || ~(sopts.cold_retry && strcmp(seed_kind0, 'warm'))
    state = state1;  info = info1;  seed_kind = seed_kind0;
else
    % deterministic cold retry: fresh cold start, the failed warm attempt is discarded
    % entirely (not merged/blended) -- point 3 of the brief.
    [state2, info2] = run_attempt(node, zeros(size(node.wn)), 0, sopts);
    state = state2;  info = info2;  seed_kind = 'cold_after_warm_fail';
end
info.seed_kind = seed_kind;
info.accepted  = info.res.accepted;
end

% =============================================================================================
function [state, info] = run_attempt(node, Sigma0, K0s0, sopts)
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
        outer_used = outer;
        break;
    end
    K(1) = K0s;
    % (3)-(5) lambdas, ordered Sigma map, pre-mix step -- MIRRORS both loops verbatim
    lam = invz_lambdas(K, g, wts, node.beta, [1 2 3]);
    sg  = invz_sigma_ordered(node.tl, lam, K, g, node.beta);
    dS  = max(abs(sg.Sigma - Sigma));
    if sopts.trace
        iters = append_iter(iters, outer, dS, sout, K0s, Gstat_it, Jnu_flat);
    end
    % (6) damped mix
    Sigma = Sigma + sopts.mix_outer*(sg.Sigma - Sigma);
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
              'medium', medium, 'medium_status', medium_status);
end

% =============================================================================================
function iters = append_iter(iters, outer, dS, sout, K0s, Gstat_it, Jnu_flat)
%APPEND_ITER One per-outer-iteration trace record (sopts.trace only), mirroring
% invz_hmf_ordered.m eval_node's own per-iteration record (:329-349) field-for-field except
% `dS` (named `resid_map` there; `dS` is the name this brief specifies).
Dq = 1 + (Jnu_flat(:) - K0s) .* Gstat_it;
irec = struct('outer', outer, 'dS', dS, 'resid_static', sout.resid, 'K0', K0s, ...
              'D_uni', sout.D_uni, 'converged_flag', sout.converged, ...
              'Dq_min', min(Dq), 'Dq_max', max(Dq), 'Dq_neg_count', nnz(Dq <= 0), ...
              'idx_pos_flat', NaN, 'Dq_pos_val', NaN, ...
              'idx_neg_flat', NaN, 'Dq_neg_val', NaN);
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
