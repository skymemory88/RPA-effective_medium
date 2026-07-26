function res = invz_ordered_residual(node, state, opts)
%INVZ_ORDERED_RESIDUAL Pure, non-mutating complete-residual checker for a jensen ordered
% node's exported (Sigma, K, lam, K0s) tuple (stage-2c task 1a/9; contract:
% docs/invz_ordered_residual_contract.md -- read that file for the full derivation, units,
% and scale justification of every number used below; NOTHING here is tuned to a run).
% G = -chi (meV^-1), ferromagnetic positive J.
%
% Implements four independently-recomputed residual blocks:
%   A -- outer Sigma map: one full independent pass of the map, seeded from the exported
%        lam/K0s, starting from the exported Sigma; rA = max|F(Sigma)-Sigma|.
%   B -- static medium (REVISED IN PLACE for the strict scheme, task 9): under 'resummed'
%        (unchanged), a fresh invz_emt_static_ordered call at the exported Sigma(1)/lam(1:2),
%        rB = so.resid, the q-average closure residual. Under a strict_1z_* scheme, the
%        load-bearing residual instead becomes the ALGEBRAIC check |K0s - Kstrict(Gref)|,
%        independently recomputed from the exported state -- see the strict-only domain
%        preflight below and the contract's dated Block-B subsection.
%   C -- Sigma self-consistency of the DERIVED lambda/Sigma chain: lam and sg recomputed
%        FRESH from the exported K (never from state.lam -- lam is derived); this is
%        production's existing final_resid, named and independently scaled here.
%   D -- dynamic EMT identity (the block final_resid OMITS): a fresh invz_emt_scalar call
%        at the exported Sigma; compares the exported K(2:end) against the fresh med.K(2:end)
%        (K(1) excluded BY DESIGN -- it is the elastic-hybrid static value, a different
%        physical quantity from invz_emt_scalar's own ordinary-Dyson K(1)). Gated on
%        med.dynamic_converged (slots 2:end), never med.converged: ordered callers replace
%        the discarded PM static slot before lambdas, so the whole-PM flag must not vote here.
%
% A strict scheme (node.eso.static_medium in {'strict_1z_dyson_ref','strict_1z_bare_ref'})
% additionally runs an independent reference/closure domain preflight BEFORE Blocks A/C/D or
% either EMT leaf are touched (invz_static_medium_reference + invz_medium_moment_closure on
% node.G0bare0/node.Jmom). A non-'ok' status returns a complete, non-accepted res immediately
% (blockA/C/D carry resid=NaN/pass=false; res.stability is 'undefined') -- an invalid static
% medium must never be fed into Jensen's local denominators.
%
% node  (struct, REQUIRED fields, constant across one node's outer loop -- see
%        invz_solve_point_ordered.m:177-205 / invz_hmf_ordered.m eval_node preamble):
%   tl, G0, g, wts, wn, beta, J0eff, G0inel0, G0el0, G0bare0, eso, eopts, Jnu_flat
%   (G0bare0 is accepted for provenance/documentation and is not read by any block under a
%   'resummed' scheme -- see the contract Sec. 2 -- but IS read by the strict-only domain
%   preflight and by res.stability.crit, both below, once node.eso.static_medium selects a
%   strict_1z_* scheme.)
%   Jmom (invz_coupling_moments of the static Jnu_flat column) is REQUIRED once
%   node.eso.static_medium selects a strict_1z_* scheme (spec SS4.3); harmless-absent under
%   'resummed'.
% state (struct, REQUIRED fields, the EXPORTED node state under test):
%   Sigma (size(node.wn)), K (size(node.wn)), lam ([3x1], invz_lambdas p=[1 2 3]), K0s (scalar)
% opts (OPTIONAL):
%   .tol_outer (default 1e-8): the outer Sigma-loop tolerance (matches the production
%       loops' own tolo default); gates blocks A/C directly and, non-dimensionalized by
%       max|Jnu_flat|, sets block D's absolute floor.
%   .dS (default NaN): the calling loop's own last PRE-MIX step (max|sg.Sigma-Sigma| from
%       the live iteration) -- used ONLY for the .stall diagnostic. Omitted/NaN disables
%       .stall (reported as NaN, "undecidable without dS").
%   .K_atol (default 1e-14), .K_rtol (default 1e-12): block B's strict-scheme gate,
%       gate = K_atol + K_rtol*max(|K0s|,|Kstrict|,Jscale) (prereg SS2) -- a mis-wiring
%       floor, not a physics tolerance: a correctly-wired one-shot call is exact.
%   .crit_tol (default 1e-6): the dimensionless res.stability classifier's crit floor
%       (prereg SS1).
%   .debug_resummed (default false): under a strict scheme only, additionally runs the
%       discarded resummed q-average closure as a nullable diagnostic
%       (res.blockB.resid_resummed); never affects .finite/.accepted/.pass.
%
% res (struct):
%   .blockA/.blockB/.blockC/.blockD -- each .resid, .scale_abs, .scale_rel, .pass (block B
%       additionally carries .converged, .status, .scheme, .ref_denom, .omit_mu3,
%       .omit_cubic, and -- only when opts.debug_resummed -- .resid_resummed; see the
%       contract Sec. 4 and its dated Block-B subsection).
%   .D_uni, .Dq_min, .Dq_max -- block B's independent recomputation, exposed per the
%       contract (never folded away: a sign-changing Dq/D_uni may be PHYSICAL instability,
%       not noise -- stage2c-context.md).
%   .stability -- struct('crit','D_uni','Dq_min','class','pass'), class in {'stable',
%       'unstable','boundary_band','undefined'} (prereg SS1). Computed for EVERY node but
%       folded into .accepted for NONE: intermediate path nodes are the unstable Landau
%       interval by construction, so per-node stability is a caller-owned endpoint check,
%       never a per-node gate.
%   .finite -- state (Sigma/K/lam/K0s) AND all four raw residuals finite.
%   .stall -- diagnostic only (see opts.dS above); NEVER used to set .accepted.
%   .accepted -- res.finite && blockA.pass && blockB.pass && blockC.pass && blockD.pass
%       (res.stability.pass is NOT a term -- see .stability above).
%
% Non-mutating (contract invariant): node/state are read-only throughout; only local copies
% are used, and nothing is ever assigned back into a field of either input. Every block
% recomputes independently from those local copies -- never from another block's
% intermediate result, never from any cached solve output.
%
% Error policy: this checker contains NO exception absorber (task 9 removed the previous
% per-block safe_eval). invz:residualNode/invz:residualState are raised directly by this
% function's own input validation; every other exception -- including any invz:* identifier
% raised by a called helper -- escapes to the caller unchanged. Physics non-convergence and
% strict-medium domain outcomes already surface as return values/statuses (out.converged,
% out.medium_status, the preflight below, ...), so an exception reaching this layer is a
% wiring/programming defect, not a numerical one, and must not be silently downgraded into a
% failed block. Each block's .err field is retained as '' for schema compatibility only.
if nargin < 3 || isempty(opts), opts = struct(); end
tol_outer = getf(opts, 'tol_outer', 1e-8);
dS_in     = getf(opts, 'dS', NaN);

req_node  = {'tl','G0','g','wts','wn','beta','J0eff','G0inel0','G0el0','G0bare0','eso','eopts','Jnu_flat'};
for k = 1:numel(req_node)
    if ~isfield(node, req_node{k})
        error('invz:residualNode', 'node is missing required field ''%s''.', req_node{k});
    end
end

% Jmom is required only once a STRICT scheme is selected (spec SS4.3): legacy/resummed
% direct-node fixtures may omit it without changing their numerical path.
smid_node = getf(node.eso, 'static_medium', 'resummed');
if ~strcmp(smid_node, 'resummed') && (~isfield(node, 'Jmom') || isempty(node.Jmom))
    error('invz:residualNode', ['node.Jmom is required under static_medium ''%s'' ' ...
        '(invz_coupling_moments of the static coupling column).'], smid_node);
end

req_state = {'Sigma','K','lam','K0s'};
for k = 1:numel(req_state)
    if ~isfield(state, req_state{k})
        error('invz:residualState', 'state is missing required field ''%s''.', req_state{k});
    end
end

% ---- local, read-only copies -- node/state are NEVER assigned into below. --------------
tl = node.tl;  G0 = node.G0(:);  g = node.g(:);  wts = node.wts(:);
beta = node.beta;  J0eff = node.J0eff;  G0inel0 = node.G0inel0;  G0el0 = node.G0el0;
eso = node.eso;  eopts = node.eopts;  Jnu_flat = node.Jnu_flat;
if isvector(Jnu_flat), Jnu_flat = Jnu_flat(:); end

Sigma = state.Sigma(:);  K = state.K(:);  lam = state.lam(:);  K0s = state.K0s;

if numel(Sigma) ~= numel(K) || numel(Sigma) ~= numel(node.wn)
    error('invz:residualState', ...
        'size mismatch: numel(Sigma)=%d, numel(K)=%d, numel(node.wn)=%d must all agree.', ...
        numel(Sigma), numel(K), numel(node.wn));
end
if numel(lam) ~= 3
    error('invz:residualState', 'state.lam must have 3 elements (p = [1 2 3]); got %d.', numel(lam));
end

Jscale = max(abs(Jnu_flat(:)));                  % problem-native coupling scale (meV)

% =========================================================================================
% Strict-only domain preflight (independent reference/closure recomputation BEFORE Blocks
% A/C/D or either EMT leaf are ever touched). Mirrors, precedence for precedence, the
% ref/closure status combination invz_emt_static_ordered's own strict branch uses internally,
% so this independently-recomputed verdict can never disagree with the node's own. A non-'ok'
% status means the static medium itself is out of domain: feeding it into Jensen's local
% denominators would manufacture a residual, not measure one, so this returns a complete,
% non-accepted res immediately -- local_F/local_blockB/local_blockC/invz_emt_scalar/
% invz_lambdas/invz_sigma_ordered are never called on this path.
% =========================================================================================
if ~strcmp(smid_node, 'resummed')
    [pf_Gref, pf_ref] = invz_static_medium_reference(node.G0bare0, Sigma(1), smid_node, eso);
    [~, pf_clo] = invz_medium_moment_closure(pf_Gref, node.Jmom, smid_node);
    if strcmp(pf_ref.status, 'ok') && strcmp(pf_clo.status, 'ok')
        pf_status = 'ok';
    elseif ~strcmp(pf_ref.status, 'ok')
        pf_status = pf_ref.status;
    else
        pf_status = pf_clo.status;
    end
    if ~strcmp(pf_status, 'ok')
        pfB_abs = getf(opts, 'K_atol', 1e-14);  pfB_rel = getf(opts, 'K_rtol', 1e-12);
        res.blockA = local_dead_block(tol_outer, tol_outer);
        res.blockB = struct('resid', NaN, 'scale_abs', pfB_abs, 'scale_rel', pfB_rel, ...
                             'pass', false, 'converged', false, 'err', '', 'status', pf_status, ...
                             'scheme', smid_node, 'ref_denom', pf_ref.denom, ...
                             'omit_mu3', pf_clo.omit_mu3, 'omit_cubic', pf_clo.omit_cubic);
        res.blockC = local_dead_block(tol_outer, tol_outer);
        res.blockD = local_dead_block(tol_outer*Jscale, tol_outer);
        res.D_uni  = NaN;  res.Dq_min = NaN;  res.Dq_max = NaN;
        res.stability = local_undefined_stability();
        res.finite = false;
        if isnan(dS_in)
            res.stall = NaN;
        else
            res.stall = isfinite(dS_in) && (dS_in < tol_outer);   % res.finite is always false here
        end
        res.accepted = false;
        return;
    end
end

% =========================================================================================
% Block A -- outer Sigma map (contract Sec. 4): one independent full pass of the map,
% seeded from state.lam/state.K0s, starting from state.Sigma.
% =========================================================================================
sgA_Sigma = local_F(tl, G0, Jnu_flat, eopts, g, wts, beta, J0eff, ...
                     G0inel0, G0el0, eso, Sigma, lam, K0s);
rA = robust_max_abs(sgA_Sigma, Sigma);
scaleA_abs = tol_outer;  scaleA_rel = tol_outer;   % Sigma dimensionless O(1): abs === rel
passA = isfinite(rA) && (rA < scaleA_abs);
res.blockA = struct('resid', rA, 'scale_abs', scaleA_abs, 'scale_rel', scaleA_rel, ...
                     'pass', passA, 'err', '');

% =========================================================================================
% Block B -- static medium (contract Sec. 4, REVISED IN PLACE for the strict scheme).
%   'resummed'          : unchanged -- the q-average closure residual from a fresh
%                         invz_emt_static_ordered at the exported Sigma(1)/lam(1:2).
%   strict_1z_*         : the load-bearing residual becomes the ALGEBRAIC check
%                         |K0s - Kstrict(Gref)|, independently recomputed from the exported
%                         state. It is identically zero for a correctly wired one-shot call,
%                         so it exists to catch MIS-WIRING, which is exactly the class of
%                         defect a prefix-matching error catch used to hide.
% The discarded resummed closure is NOT run under a strict scheme unless
% opts.debug_resummed: doing so would restore the inner iteration and pole exposure this
% design removes, and the analytic-continuation path deliberately crosses that pole.
% =========================================================================================
outB = local_blockB(tl, lam, Sigma, Jnu_flat, K0s, beta, J0eff, ...
                    G0inel0, G0el0, eso);
strictB = ~strcmp(smid_node, 'resummed');
if strictB
    scaleB_abs = getf(opts, 'K_atol', 1e-14);
    scaleB_rel = getf(opts, 'K_rtol', 1e-12);
else
    rtolB = getf(eso, 'resid_tol', 1e-10);
    scaleB_abs = rtolB;  scaleB_rel = rtolB;
end
statusB = 'nonfinite';  omit3 = NaN;  omit4 = NaN;  refdenB = NaN;
Gstat_b = outB.Gstat;
D_uni   = outB.so.D_uni;
Dq      = 1 + (Jnu_flat(:) - outB.K0) .* Gstat_b;
Dq_min  = min(Dq);  Dq_max = max(Dq);
convB   = outB.so.converged;
if strictB
    statusB = outB.so.medium_status;
    omit3 = outB.so.omit_mu3;  omit4 = outB.so.omit_cubic;
    refdenB = outB.so.medium.ref.denom;
    if strcmp(statusB, 'ok')
        Kstrict = outB.so.medium.closure.Kstrict;
        rB = abs(K0s - Kstrict);
        gate = scaleB_abs + scaleB_rel*max([abs(K0s), abs(Kstrict), Jscale]);
        passB = isfinite(rB) && (rB < gate);
    else
        rB = NaN;  passB = false;  convB = false;    % domain status, not an exception
    end
else
    statusB = 'not_applicable';
    rB = outB.so.resid;
    passB = isfinite(rB) && (rB < scaleB_abs + scaleB_rel*abs(Gstat_b));
end
res.blockB = struct('resid', rB, 'scale_abs', scaleB_abs, 'scale_rel', scaleB_rel, ...
                     'pass', passB, 'converged', convB, 'err', '', 'status', statusB, ...
                     'scheme', smid_node, 'ref_denom', refdenB, ...
                     'omit_mu3', omit3, 'omit_cubic', omit4);
if strictB && getf(opts, 'debug_resummed', false)
    esoR = eso;  esoR.static_medium = 'resummed';  esoR.warn = false;
    outR = local_blockB(tl, lam, Sigma, Jnu_flat, K0s, beta, J0eff, ...
                        G0inel0, G0el0, esoR);
    res.blockB.resid_resummed = outR.so.resid;
end

% =========================================================================================
% Block C -- Sigma self-consistency of the DERIVED lambda/Sigma chain (contract Sec. 4):
% lam and sg recomputed FRESH from the exported K (never from state.lam).
% =========================================================================================
rC = local_blockC(tl, K, g, wts, beta, Sigma);
scaleC_abs = tol_outer;  scaleC_rel = tol_outer;
passC = isfinite(rC) && (rC < scaleC_abs);
res.blockC = struct('resid', rC, 'scale_abs', scaleC_abs, 'scale_rel', scaleC_rel, ...
                     'pass', passC, 'err', '');

% =========================================================================================
% Block D -- dynamic EMT identity (contract Sec. 4; the block final_resid OMITS): fresh
% invz_emt_scalar at the exported Sigma; compares the exported K(2:end) to the fresh
% med.K(2:end). K(1) EXCLUDED by design (elastic-hybrid vs ordinary-Dyson -- see contract).
% Gated on med.dynamic_converged (slots 2:end), never med.converged (whole-PM, slot 1
% included): ordered callers replace the discarded PM static slot before lambdas, so letting
% that slot vote here would reintroduce a static-sector veto through a dynamic residual.
% =========================================================================================
medD = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);
finite_D = all(isfinite(K)) && all(isfinite(lam)) && all(isfinite(Sigma));
scaleD_abs = tol_outer * Jscale;  scaleD_rel = tol_outer;
if numel(K) > 1
    rD    = robust_max_abs(K(2:end), medD.K(2:end));
    refD  = max(abs(medD.K(2:end)));
    passD = isfinite(rD) && (rD < scaleD_abs + scaleD_rel*refD) && finite_D && medD.dynamic_converged;
else                                          % degenerate-size guard (contract Sec. 4)
    rD    = 0;
    passD = finite_D && medD.dynamic_converged;
end
res.blockD = struct('resid', rD, 'scale_abs', scaleD_abs, 'scale_rel', scaleD_rel, ...
                     'pass', passD, 'err', '');

% ---- top-level diagnostics / finite / stall / aggregate (contract Sec. 6-7) -------------
res.D_uni  = D_uni;
res.Dq_min = Dq_min;
res.Dq_max = Dq_max;
% ---- STABILITY TIER (spec SS1): computed for every node, folded into .accepted for NONE.
% Intermediate path nodes may legitimately be unstable -- they are the analytic continuation
% through the unstable Landau interval, and requiring per-node positivity would re-mask the
% ordered phase. Only an ENDPOINT root is held to this tier, by the caller.
crit_tol = getf(opts, 'crit_tol', 1e-6);
if all(isfinite([D_uni, Dq_min]))
    critv = outB.so.r + J0eff*node.G0bare0;
    D_tol = 1e-6*max(1, abs(Gstat_b)*Jscale);        % prereg SS1: noise scales with |Gstat|
    Dq_tol = D_tol;
    if ~isfinite(critv)
        cls = 'undefined';
    elseif critv > crit_tol && D_uni > D_tol && Dq_min > D_tol
        cls = 'stable';
    elseif critv < -crit_tol || D_uni < -D_tol || Dq_min < -Dq_tol
        cls = 'unstable';
    else
        cls = 'boundary_band';
    end
    res.stability = struct('crit', critv, 'D_uni', D_uni, 'Dq_min', Dq_min, ...
                           'class', cls, 'pass', strcmp(cls, 'stable'));
else
    res.stability = local_undefined_stability();
end
res.finite = all(isfinite(Sigma)) && all(isfinite(K)) && all(isfinite(lam)) && isfinite(K0s) ...
             && isfinite(rA) && isfinite(rB) && isfinite(rC) && isfinite(rD);
if isnan(dS_in)
    res.stall = NaN;
else
    res.stall = isfinite(dS_in) && (dS_in < tol_outer) && ...
                ~(passA && passB && passC && passD && res.finite);
end
res.accepted = res.finite && passA && passB && passC && passD;
end

% =============================================================================================
function b = local_dead_block(scale_abs, scale_rel)
%LOCAL_DEAD_BLOCK A/C/D placeholder for the strict-only domain preflight short circuit above:
% an invalid static medium is never fed into Jensen's local denominators, so these blocks are
% never evaluated on that path -- resid/pass report that honestly (NaN/false) rather than a
% fabricated number, while still keeping the exact schema a normal return has.
b = struct('resid', NaN, 'scale_abs', scale_abs, 'scale_rel', scale_rel, 'pass', false, 'err', '');
end

% =============================================================================================
function s = local_undefined_stability()
%LOCAL_UNDEFINED_STABILITY The 'undefined' res.stability struct, shared by the strict-only
% domain preflight short circuit and the ordinary D_uni/Dq_min-non-finite case below -- ONE
% definition, so the two call sites can never drift apart.
s = struct('crit', NaN, 'D_uni', NaN, 'Dq_min', NaN, 'class', 'undefined', 'pass', false);
end

% =============================================================================================
function sgSigma = local_F(tl, G0, Jnu_flat, eopts, g, wts, beta, J0eff, G0inel0, G0el0, ...
                            eso, Sigma, lam_seed, K0s_seed)
%LOCAL_F One full independent pass of the outer map's body (contract Sec. 3), mirroring
% invz_solve_point_ordered.m:206-221 / invz_hmf_ordered.m eval_node:313-328 verbatim.
med = invz_emt_scalar(G0, Sigma, Jnu_flat, eopts);              % (1) dynamic sector
Kf  = med.K;
[K0s_new, ~, ~] = invz_emt_static_ordered(tl, lam_seed(1:2), Sigma(1), Jnu_flat, K0s_seed, ...
                                          beta, J0eff, G0inel0, G0el0, eso);   % (2) static
Kf(1) = K0s_new;
lam_next = invz_lambdas(Kf, g, wts, beta, [1 2 3]);             % (3) derived lambdas
sg = invz_sigma_ordered(tl, lam_next, Kf, g, beta);             % (4) ordered Sigma map
sgSigma = sg.Sigma;
end

% =============================================================================================
function out = local_blockB(tl, lam, Sigma, Jnu_flat, K0s, beta, J0eff, G0inel0, G0el0, eso)
%LOCAL_BLOCKB Fresh static-closure recomputation at the exported Sigma(1)/lam(1:2), seeded
% from the exported K0s (contract Sec. 4, Block B).
[K0_b, Gstat_b, so] = invz_emt_static_ordered(tl, lam(1:2), Sigma(1), Jnu_flat, K0s, beta, ...
                                              J0eff, G0inel0, G0el0, eso);
out = struct('K0', K0_b, 'Gstat', Gstat_b, 'so', so);
end

% =============================================================================================
function rC = local_blockC(tl, K, g, wts, beta, Sigma)
%LOCAL_BLOCKC Fresh lam/sg recomputation from the exported K (contract Sec. 4, Block C;
% mirrors invz_solve_point_ordered.m:232-234's final_resid exactly).
lam_check = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg = invz_sigma_ordered(tl, lam_check, K, g, beta);
rC = robust_max_abs(sg.Sigma, Sigma);
end

% =============================================================================================
function r = robust_max_abs(a, b)
%ROBUST_MAX_ABS max(abs(a-b)), but NaN-PROPAGATING. MATLAB's max() ignores individual NaN
% entries by default (unlike sum/mean), so a bare max(abs(a-b)) would silently ignore a single
% non-finite component in an otherwise-finite vector rather than reporting the block as failed
% -- exactly the kind of silent masking this diagnostic stage must not do. Any non-finite entry
% in either operand makes the whole residual NaN (isfinite(NaN) is false, so the block's own
% .pass then correctly reads false, without relying on the separate top-level .finite flag to
% catch it).
if any(~isfinite(a(:))) || any(~isfinite(b(:)))
    r = NaN;
else
    r = max(abs(a(:) - b(:)));
end
end
