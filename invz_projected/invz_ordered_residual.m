function res = invz_ordered_residual(node, state, opts)
%INVZ_ORDERED_RESIDUAL Pure, non-mutating complete-residual checker for a jensen ordered
% node's exported (Sigma, K, lam, K0s) tuple (stage-2c task 1a/9; contract summarized
% in invzp_convg_diagnosis.md Section 2.4; NOTHING here is tuned to a run).
% G = -chi (meV^-1), ferromagnetic positive J.
%
% Implements four independently-recomputed residual blocks:
%   A -- Sigma map. opts.formulation='nested' (default) replays one full legacy outer
%        iteration, including the inner static Picard solve. 'coupled' instead holds the
%        exported K0 fixed while independently rebuilding dynamic K, lambda, and Sigma;
%        the static equation is then graded separately by defactored Block B. The latter
%        is the simultaneous [Sigma,K0] formulation and does not ask an ambiguous inner
%        solver to select a second static root while auditing the supplied one. No nested
%        static solve is run on the coupled path unless opts.debug_legacy_nested is true.
%   B -- static medium (REVISED IN PLACE for the strict scheme, task 9): under the default
%        resummed audit coordinate, a fresh invz_emt_static_ordered call at the exported
%        Sigma(1)/lam(1:2), rB = so.resid, the q-average closure residual. The diagnostic
%        node option eso.audit_coordinate='defactored' instead independently measures the
%        equivalent, well-conditioned |K0-Jloc|/Jscale equation while retaining the raw
%        closure as blockB.resid_raw. Under a strict_1z_* scheme, the
%        load-bearing residual instead becomes the ALGEBRAIC check |K0s - Kstrict(Gref)|,
%        independently recomputed from the exported state -- see the strict-only domain
%        preflight below and the contract's dated Block-B subsection.
%   C -- Sigma self-consistency of the DERIVED lambda/Sigma chain: lam and sg recomputed
%        FRESH from the exported K (never from state.lam -- lam is derived); this is
%        production's existing final_resid, named and independently scaled here. In the
%        coupled formulation it also verifies the exact derived-state wiring
%        state.lam == invz_lambdas(state.K,...).
%   D -- dynamic EMT identity (the block final_resid OMITS): a fresh invz_emt_scalar call
%        at the exported Sigma; compares the exported K(2:end) against the fresh med.K(2:end)
%        (K(1) excluded BY DESIGN -- it is the elastic-hybrid static value, a different
%        physical quantity from invz_emt_scalar's own ordinary-Dyson K(1)). Gated on
%        med.dynamic_converged (slots 2:end), never med.converged: ordered callers replace
%        the discarded PM static slot before lambdas, so the whole-PM flag must not vote here.
%        Coupled mode separately verifies the exact wiring identity state.K(1)==state.K0s;
%        it still never compares state.K(1) with the ordinary-Dyson med.K(1).
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
%   .formulation ('nested' default | 'coupled'): Block-A definition. 'coupled' is
%       diagnostic-only, supported only for resummed static medium with
%       node.eso.audit_coordinate='defactored'.
%   .debug_legacy_nested (default false): with formulation='coupled', additionally run and
%       report the legacy nested replay. This can be branch-ambiguous and expensive, and
%       never votes on coupled acceptance.
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
%       not noise; see invzp_convg_diagnosis.md).
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
formulation = getf(opts,'formulation','nested');
debugLegacyNested = getf(opts,'debug_legacy_nested',false);
if ~(islogical(debugLegacyNested) && isscalar(debugLegacyNested))
    error('invz:residualNode', ...
        'opts.debug_legacy_nested must be a logical scalar.');
end
if isstring(formulation)
    if ~isscalar(formulation) || ismissing(formulation)
        error('invz:residualNode', ...
            'opts.formulation must be ''nested'' or ''coupled''.');
    end
    formulation = char(formulation);
end
if ~(ischar(formulation) && isrow(formulation) && ...
        any(strcmp(formulation,{'nested','coupled'})))
    error('invz:residualNode', ...
        'opts.formulation must be ''nested'' or ''coupled''.');
end

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
auditCoordinate = getf(eso,'audit_coordinate','raw_closure');
if ~strcmp(smid_node,'resummed') && strcmp(formulation,'coupled')
    error('invz:residualNode', ...
        'opts.formulation=''coupled'' currently supports static_medium=''resummed'' only.');
end
if strcmp(smid_node,'resummed') && ...
        ~any(strcmp(auditCoordinate,{'raw_closure','defactored'}))
    error('invz:residualNode', ...
        'node.eso.audit_coordinate must be ''raw_closure'' or ''defactored''.');
end
if strcmp(formulation,'coupled') && ~strcmp(auditCoordinate,'defactored')
    error('invz:residualNode', ...
        ['opts.formulation=''coupled'' requires ', ...
         'node.eso.audit_coordinate=''defactored''.']);
end
res = struct('formulation',formulation);

% =========================================================================================
% Strict-only domain preflight (independent reference/closure recomputation BEFORE Blocks
% A/C/D or either EMT leaf are ever touched). This calls invz_static_medium_reference /
% invz_medium_moment_closure on node.G0bare0 and node.Jmom, in the SAME precedence order
% invz_emt_static_ordered's own strict branch uses internally (:86-91) -- but that branch reads
% G0inel0+G0el0 and eso.Jmom (falling back to invz_coupling_moments(Jf) when eso.Jmom is
% absent), a DIFFERENT provenance path for nominally the same quantities. The two verdicts
% therefore agree -- and this preflight's status matches invz_emt_static_ordered's own -- only
% when node.G0bare0 is bitwise-equal to node.G0inel0+node.G0el0 and node.Jmom matches what
% eso.Jmom (or its fallback) would produce. Measured true at every fixture to date (bitwise,
% task-9 review), but NOT enforced here and not guaranteed for an arbitrary caller. Block B's
% defensive `else rB = NaN; passB = false; convB = false;` branch below exists precisely for a
% caller where that agreement fails; it is unreachable today only because the equality above
% happens to hold everywhere it has been measured. A non-'ok' status here means the static
% medium itself is out of domain: feeding it into Jensen's local denominators would manufacture
% a residual, not measure one, so this returns a complete, non-accepted res immediately --
% local_F/local_blockB/local_blockC/invz_emt_scalar/invz_lambdas/invz_sigma_ordered are never
% called on this path.
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
        if getf(opts, 'debug_resummed', false)
            % Schema parity with the normal path (below, near the debug_resummed block after
            % local_blockB): a flag-keyed Task 10/12 consumer must see the same
            % fieldnames(res.blockB) regardless of which path produced res. NaN, not computed:
            % the whole point of this preflight is that no EMT leaf runs on an invalid medium,
            % so there is no resummed closure to evaluate here even as an opt-in diagnostic.
            res.blockB.resid_resummed = NaN;
        end
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
% Block A -- Sigma map. The nested formulation replays the legacy outer map, including
% its inner static solver. The coupled formulation independently rebuilds dynamic K,
% replaces only its static slot by the supplied K0, and then rebuilds lambda/Sigma. In
% the latter formulation Block B, not an inner iteration, grades the static equation.
% =========================================================================================
if strcmp(formulation,'coupled')
    sgA_Sigma = local_F_coupled( ...
        tl,G0,Jnu_flat,eopts,g,wts,beta,Sigma,K0s);
else
    sgA_Sigma = local_F(tl, G0, Jnu_flat, eopts, g, wts, beta, J0eff, ...
                         G0inel0, G0el0, eso, Sigma, lam, K0s);
end
rA = invz_finite_max_abs(sgA_Sigma, Sigma);
scaleA_abs = tol_outer;  scaleA_rel = tol_outer;   % Sigma dimensionless O(1): abs === rel
passA = isfinite(rA) && (rA < scaleA_abs);
res.blockA = struct('resid', rA, 'scale_abs', scaleA_abs, 'scale_rel', scaleA_rel, ...
                     'pass', passA, 'err', '');
if strcmp(formulation,'coupled')
    res.blockA.legacy_nested_computed = debugLegacyNested;
    res.blockA.legacy_nested_resid = NaN;
    res.blockA.legacy_nested_pass = false;
    if debugLegacyNested
        nestedSigma = local_F(tl,G0,Jnu_flat,eopts,g,wts,beta,J0eff, ...
            G0inel0,G0el0,eso,Sigma,lam,K0s);
        nestedResidual = invz_finite_max_abs(nestedSigma,Sigma);
        res.blockA.legacy_nested_resid = nestedResidual;
        res.blockA.legacy_nested_pass = ...
            isfinite(nestedResidual) && nestedResidual < scaleA_abs;
    end
end

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
strictB = ~strcmp(smid_node, 'resummed');
defactoredB = ~strictB && strcmp(auditCoordinate,'defactored');
if strcmp(formulation,'coupled')
    % The coupled formulation requires this coordinate by preflight above. Build a
    % schema-compatible local record without calling the branch-selecting nested solver.
    directB = local_blockB_coupled(tl,G0,Sigma,K0s,Jnu_flat,eopts, ...
        g,wts,beta,G0inel0,G0el0);
    directDUni = 1+(J0eff-K0s)*directB.Gstat;
    outB = struct('K0',NaN,'Gstat',directB.Gstat, ...
        'so',struct('r',directB.r,'D_uni',directDUni, ...
        'converged',false,'resid',directB.resid_raw));
else
    outB = local_blockB(tl,lam,Sigma,Jnu_flat,K0s,beta,J0eff, ...
        G0inel0,G0el0,eso);
end
if strictB
    scaleB_abs = getf(opts, 'K_atol', 1e-14);
    scaleB_rel = getf(opts, 'K_rtol', 1e-12);
elseif defactoredB
    scaleB_abs = tol_outer;
    scaleB_rel = tol_outer;
else
    rtolB = getf(eso, 'resid_tol', 1e-10);
    scaleB_abs = rtolB;  scaleB_rel = rtolB;
end
omit3 = NaN;  omit4 = NaN;  refdenB = NaN;
Gstat_b = outB.Gstat;
rStatic = outB.so.r;
D_uni   = outB.so.D_uni;
Dq      = 1 + (Jnu_flat(:) - outB.K0) .* Gstat_b;
Dq_min  = min(Dq);  Dq_max = max(Dq);
% so.converged means different things by scheme (contract Sec. 4's base Block-B text describes
% only the 'resummed' meaning): under 'resummed' it is so.resid < rtol_B, a PURE-ABSOLUTE
% closure flag independent of blockB.pass's own combined gate; under strict it is a
% domain-validity flag (strcmp(medium.status,'ok')), forced false in the domain-failure branch
% below. Tasks 10/12 reading res.blockB.converged must branch on res.blockB.scheme, not assume
% either meaning holds unconditionally.
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
elseif defactoredB
    if ~strcmp(formulation,'coupled')
        directB = local_blockB_defactored(tl,lam,Sigma,K0s,beta, ...
            G0inel0,G0el0,Jnu_flat);
    end
    rB = directB.resid;
    passB = isfinite(rB) && (rB < scaleB_abs);
    statusB = 'not_applicable';
    Gstat_b = directB.Gstat;
    rStatic = directB.r;
    D_uni = 1+(J0eff-K0s)*Gstat_b;
    Dq = 1+(Jnu_flat(:)-K0s).*Gstat_b;
    Dq_min = min(Dq);  Dq_max = max(Dq);
    convB = passB;
else
    statusB = 'not_applicable';
    rB = outB.so.resid;
    passB = isfinite(rB) && (rB < scaleB_abs + scaleB_rel*abs(Gstat_b));
end
% blockB.pass is likewise scheme-dependent: under 'resummed' it is the q-average closure gate
% (isfinite(rB) && rB < scale_abs+scale_rel*|Gstat_b|, set above); under strict it is the
% algebraic |K0s-Kstrict| mis-wiring check AND status=='ok' (both already folded into passB by
% the if/else above) -- see the contract's dated Block-B subsection for the full scheme split.
res.blockB = struct('resid', rB, 'scale_abs', scaleB_abs, 'scale_rel', scaleB_rel, ...
                     'pass', passB, 'converged', convB, 'err', '', 'status', statusB, ...
                     'scheme', smid_node, 'ref_denom', refdenB, ...
                     'omit_mu3', omit3, 'omit_cubic', omit4);
if defactoredB
    res.blockB.coordinate = 'defactored_K0';
    res.blockB.resid_raw = directB.resid_raw;
    res.blockB.raw_gate = getf(eso,'resid_tol',1e-10)* ...
        (1+abs(directB.Gstat));
    res.blockB.raw_pass = isfinite(res.blockB.resid_raw) && ...
        res.blockB.resid_raw < res.blockB.raw_gate;
    if strcmp(formulation,'coupled')
        res.blockB.K0_seed_drift = NaN;
    else
        res.blockB.K0_seed_drift = abs(outB.K0-K0s);
    end
end
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
if strcmp(formulation,'coupled')
    lamCheck = invz_lambdas(K,g,wts,beta,[1 2 3]);
    res.blockC.lambda_consistent = isequaln(lam,lamCheck);
    res.blockC.lambda_abs_resid = abs(lam-lamCheck);
    res.blockC.pass = res.blockC.pass && res.blockC.lambda_consistent;
    passC = res.blockC.pass;
end

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
    rD    = invz_finite_max_abs(K(2:end), medD.K(2:end));
    refD  = max(abs(medD.K(2:end)));
    passD = isfinite(rD) && (rD < scaleD_abs + scaleD_rel*refD) && finite_D && medD.dynamic_converged;
else                                          % degenerate-size guard (contract Sec. 4)
    rD    = 0;
    passD = finite_D && medD.dynamic_converged;
end
res.blockD = struct('resid', rD, 'scale_abs', scaleD_abs, 'scale_rel', scaleD_rel, ...
                     'pass', passD, 'err', '');
if strcmp(formulation,'coupled')
    res.blockD.static_slot_resid = abs(K(1)-K0s);
    res.blockD.static_slot_consistent = isequaln(K(1),K0s);
    res.blockD.pass = res.blockD.pass && res.blockD.static_slot_consistent;
    passD = res.blockD.pass;
end

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
    critv = rStatic + J0eff*node.G0bare0;
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
function sgSigma = local_F_coupled( ...
        tl,G0,Jnu_flat,eopts,g,wts,beta,Sigma,K0)
%LOCAL_F_COUPLED Independent simultaneous-map recomputation. K0 is an
% independent coordinate whose closure is graded by Block B; no nested
% static solver is called here.
med = invz_emt_scalar(G0,Sigma,Jnu_flat,eopts);
K = med.K;
K(1) = K0;
lam = invz_lambdas(K,g,wts,beta,[1 2 3]);
sg = invz_sigma_ordered(tl,lam,K,g,beta);
sgSigma = sg.Sigma;
end

% =============================================================================================
function out = local_blockB_coupled( ...
        tl,G0,Sigma,K0,Jnu_flat,eopts,g,wts,beta,G0inel0,G0el0)
%LOCAL_BLOCKB_COUPLED Independently rebuild the canonical derived K/lambda tuple before
% evaluating the simultaneous static equation. This makes Block B a function of the
% independent [Sigma,K0] coordinates, exactly like INVZ_ORDERED_NODE_EQUATIONS, rather than
% trusting the exported state's derived lambda.
med = invz_emt_scalar(G0,Sigma,Jnu_flat,eopts);
K = med.K;
K(1) = K0;
lam = invz_lambdas(K,g,wts,beta,[1 2 3]);
out = local_blockB_defactored(tl,lam,Sigma,K0,beta,G0inel0,G0el0,Jnu_flat);
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
function out = local_blockB_defactored(tl,lam,Sigma,K0,beta,G0inel0,G0el0,Jnu_flat)
%LOCAL_BLOCKB_DEFACTORED Independent reciprocal-coordinate static equation.
[Gstat,go] = invz_gstat_ordered(tl,lam(1:2),K0,Sigma(1),beta, ...
    G0inel0,G0el0,struct('stable_form',true));
d0 = go.gstat_local_denom;
Hz = G0inel0+go.xi*G0el0*d0;
z = d0/Hz;
Jf = Jnu_flat(:);
Jscale = max(abs(Jf));
if isinf(z)
    Jloc = mean(Jf);
    Gbar = 0;
elseif isfinite(z)
    scale = max(abs(z),Jscale);
    weights = scale./(z+Jf-K0);
    meanWeights = mean(weights);
    Gbar = meanWeights/scale;
    Jloc = mean(Jf.*weights)/meanWeights;
else
    Gbar = NaN;
    Jloc = NaN;
end
out = struct('Gstat',Gstat,'Gbar',Gbar,'Jloc',Jloc,'r',go.r, ...
    'resid',abs(K0-Jloc)/Jscale, ...
    'resid_raw',abs(Gbar-Gstat));
end

% =============================================================================================
function rC = local_blockC(tl, K, g, wts, beta, Sigma)
%LOCAL_BLOCKC Fresh lam/sg recomputation from the exported K (contract Sec. 4, Block C;
% mirrors invz_solve_point_ordered.m:232-234's final_resid exactly).
lam_check = invz_lambdas(K, g, wts, beta, [1 2 3]);
sg = invz_sigma_ordered(tl, lam_check, K, g, beta);
rC = invz_finite_max_abs(sg.Sigma, Sigma);
end
