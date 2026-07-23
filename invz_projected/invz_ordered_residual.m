function res = invz_ordered_residual(node, state, opts)
%INVZ_ORDERED_RESIDUAL Pure, non-mutating complete-residual checker for a jensen ordered
% node's exported (Sigma, K, lam, K0s) tuple (stage-2c task 1a; contract:
% docs/invz_ordered_residual_contract.md -- read that file for the full derivation, units,
% and scale justification of every number used below; NOTHING here is tuned to a run).
%
% Implements four independently-recomputed residual blocks:
%   A -- outer Sigma map: one full independent pass of the map, seeded from the exported
%        lam/K0s, starting from the exported Sigma; rA = max|F(Sigma)-Sigma|.
%   B -- static EMT closure: a fresh invz_emt_static_ordered call at the exported
%        Sigma(1)/lam(1:2), seeded from the exported K0s; rB = so.resid.
%   C -- Sigma self-consistency of the DERIVED lambda/Sigma chain: lam and sg recomputed
%        FRESH from the exported K (never from state.lam -- lam is derived); this is
%        production's existing final_resid, named and independently scaled here.
%   D -- dynamic EMT identity (the block final_resid OMITS): a fresh invz_emt_scalar call
%        at the exported Sigma; compares the exported K(2:end) against the fresh med.K(2:end)
%        (K(1) excluded BY DESIGN -- it is the elastic-hybrid static value, a different
%        physical quantity from invz_emt_scalar's own ordinary-Dyson K(1)).
%
% node  (struct, REQUIRED fields, constant across one node's outer loop -- see
%        invz_solve_point_ordered.m:177-205 / invz_hmf_ordered.m eval_node preamble):
%   tl, G0, g, wts, wn, beta, J0eff, G0inel0, G0el0, G0bare0, eso, eopts, Jnu_flat
%   (G0bare0 is accepted for provenance/documentation; it is not read by any block below --
%   see the contract Sec. 2.)
% state (struct, REQUIRED fields, the EXPORTED node state under test):
%   Sigma (size(node.wn)), K (size(node.wn)), lam ([3x1], invz_lambdas p=[1 2 3]), K0s (scalar)
% opts (OPTIONAL):
%   .tol_outer (default 1e-8): the outer Sigma-loop tolerance (matches the production
%       loops' own tolo default); gates blocks A/C directly and, non-dimensionalized by
%       max|Jnu_flat|, sets block D's absolute floor.
%   .dS (default NaN): the calling loop's own last PRE-MIX step (max|sg.Sigma-Sigma| from
%       the live iteration) -- used ONLY for the .stall diagnostic. Omitted/NaN disables
%       .stall (reported as NaN, "undecidable without dS").
%
% res (struct):
%   .blockA/.blockB/.blockC/.blockD -- each .resid, .scale_abs, .scale_rel, .pass (block B
%       additionally carries .converged = so.converged, the closure's OWN pure-absolute
%       flag, which can legitimately differ from .pass -- see the contract Sec. 4).
%   .D_uni, .Dq_min, .Dq_max -- block B's independent recomputation, exposed per the
%       contract (never folded away: a sign-changing Dq/D_uni may be PHYSICAL instability,
%       not noise -- stage2c-context.md).
%   .finite -- state (Sigma/K/lam/K0s) AND all four raw residuals finite.
%   .stall -- diagnostic only (see opts.dS above); NEVER used to set .accepted.
%   .accepted -- res.finite && blockA.pass && blockB.pass && blockC.pass && blockD.pass.
%
% Non-mutating (contract invariant): node/state are read-only throughout; only local copies
% are used, and nothing is ever assigned back into a field of either input. Every block
% recomputes independently from those local copies -- never from another block's
% intermediate result, never from any cached solve output.
%
% Error policy: only 'invz:*' identifiers are absorbed per block (that block is then
% reported as non-finite/failing, never silently masked, and its error id is recorded in
% res.blockX.err); anything else rethrows immediately.
if nargin < 3 || isempty(opts), opts = struct(); end
tol_outer = getf(opts, 'tol_outer', 1e-8);
dS_in     = getf(opts, 'dS', NaN);

req_node  = {'tl','G0','g','wts','wn','beta','J0eff','G0inel0','G0el0','G0bare0','eso','eopts','Jnu_flat'};
for k = 1:numel(req_node)
    if ~isfield(node, req_node{k})
        error('invz:residualNode', 'node is missing required field ''%s''.', req_node{k});
    end
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
% Block A -- outer Sigma map (contract Sec. 4): one independent full pass of the map,
% seeded from state.lam/state.K0s, starting from state.Sigma.
% =========================================================================================
[sgA_Sigma, errA] = safe_eval(@() local_F(tl, G0, Jnu_flat, eopts, g, wts, beta, J0eff, ...
                                           G0inel0, G0el0, eso, Sigma, lam, K0s));
if isempty(errA)
    rA = robust_max_abs(sgA_Sigma, Sigma);
else
    rA = NaN;
end
scaleA_abs = tol_outer;  scaleA_rel = tol_outer;   % Sigma dimensionless O(1): abs === rel
passA = isfinite(rA) && (rA < scaleA_abs);
res.blockA = struct('resid', rA, 'scale_abs', scaleA_abs, 'scale_rel', scaleA_rel, ...
                     'pass', passA, 'err', errA);

% =========================================================================================
% Block B -- static EMT closure (contract Sec. 4): fresh invz_emt_static_ordered at the
% exported Sigma(1)/lam(1:2), seeded from the exported K0s.
% =========================================================================================
rtolB = getf(eso, 'resid_tol', 1e-10);
[outB, errB] = safe_eval(@() local_blockB(tl, lam, Sigma, Jnu_flat, K0s, beta, J0eff, ...
                                           G0inel0, G0el0, eso));
scaleB_abs = rtolB;  scaleB_rel = rtolB;
if isempty(errB)
    rB      = outB.so.resid;
    Gstat_b = outB.Gstat;
    passB   = isfinite(rB) && (rB < scaleB_abs + scaleB_rel*abs(Gstat_b));
    D_uni   = outB.so.D_uni;
    Dq      = 1 + (Jnu_flat(:) - outB.K0) .* Gstat_b;
    Dq_min  = min(Dq);  Dq_max = max(Dq);
    convB   = outB.so.converged;
else
    rB = NaN;  passB = false;  D_uni = NaN;  Dq_min = NaN;  Dq_max = NaN;  convB = false;
end
res.blockB = struct('resid', rB, 'scale_abs', scaleB_abs, 'scale_rel', scaleB_rel, ...
                     'pass', passB, 'converged', convB, 'err', errB);

% =========================================================================================
% Block C -- Sigma self-consistency of the DERIVED lambda/Sigma chain (contract Sec. 4):
% lam and sg recomputed FRESH from the exported K (never from state.lam).
% =========================================================================================
[rC, errC] = safe_eval(@() local_blockC(tl, K, g, wts, beta, Sigma));
if ~isempty(errC), rC = NaN; end
scaleC_abs = tol_outer;  scaleC_rel = tol_outer;
passC = isfinite(rC) && (rC < scaleC_abs);
res.blockC = struct('resid', rC, 'scale_abs', scaleC_abs, 'scale_rel', scaleC_rel, ...
                     'pass', passC, 'err', errC);

% =========================================================================================
% Block D -- dynamic EMT identity (contract Sec. 4; the block final_resid OMITS): fresh
% invz_emt_scalar at the exported Sigma; compares the exported K(2:end) to the fresh
% med.K(2:end). K(1) EXCLUDED by design (elastic-hybrid vs ordinary-Dyson -- see contract).
% =========================================================================================
[medD, errD] = safe_eval(@() invz_emt_scalar(G0, Sigma, Jnu_flat, eopts));
finite_D = all(isfinite(K)) && all(isfinite(lam)) && all(isfinite(Sigma));
scaleD_abs = tol_outer * Jscale;  scaleD_rel = tol_outer;
if isempty(errD)
    if numel(K) > 1
        rD    = robust_max_abs(K(2:end), medD.K(2:end));
        refD  = max(abs(medD.K(2:end)));
        passD = isfinite(rD) && (rD < scaleD_abs + scaleD_rel*refD) && finite_D && medD.converged;
    else                                          % degenerate-size guard (contract Sec. 4)
        rD    = 0;
        passD = finite_D && medD.converged;
    end
else
    rD = NaN;  passD = false;
end
res.blockD = struct('resid', rD, 'scale_abs', scaleD_abs, 'scale_rel', scaleD_rel, ...
                     'pass', passD, 'err', errD);

% ---- top-level diagnostics / finite / stall / aggregate (contract Sec. 6-7) -------------
res.D_uni  = D_uni;
res.Dq_min = Dq_min;
res.Dq_max = Dq_max;
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
function [out, errid] = safe_eval(fn)
%SAFE_EVAL Evaluate a zero-arg function handle; absorb only 'invz:*' errors (per-block, into
% errid, with out = []); anything else rethrows immediately (repo convention, e.g.
% invz_solve_auto.m: "catch err; if ~strncmp(err.identifier,'invz:',5), rethrow(err); end").
errid = '';
try
    out = fn();
catch err
    if ~strncmp(err.identifier, 'invz:', 5)
        rethrow(err);
    end
    out = [];
    errid = err.identifier;
end
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
