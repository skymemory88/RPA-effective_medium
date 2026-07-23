function out = invz_task2_run_config(cfg)
%INVZ_TASK2_RUN_CONFIG Per-config discriminator harness (stage-2c task 2a; spec:
% .superpowers/sdd/task2a-brief.md deliverable 2, docs/invzp_task2_prereg.md).
%
% Runs ONE discriminator configuration through the SAME production ordered path
% (invz_hmf_ordered / invz_solve_point_ordered's ordered_mode='jensen') and returns
% per-NODE records -- checker-accept?, the stability scalar/class (prereg Sec. A), the
% FULL exported state (Sigma/K/lambda/K0/Gstat/D_uni + endpoint hstar/m_star/r_star, for
% Sec. C comparisons), seed provenance, and q/branch of the closest +/-Dq modes -- plus a
% per-config summary. ADDITIVE ONLY: this file calls invz_ordered_trace / invz_hmf_ordered
% / invz_ordered_node_solve / invz_ordered_residual (Task 0/1) and invz_single_ion /
% invz_twolevel_ordered / invz_chi0z / invz_g / invz_matsubara (the same public per-node
% construction helpers invz_hmf_ordered.m's own eval_node calls, in the SAME order, with
% the SAME arguments -- see task2_build_node below); it modifies NONE of them.
%
% ---------------------------------------------------------------------------------------
% WHY a harness-level node reconstruction is needed (read before changing this file):
% invz_hmf_ordered's own trace hook (opts.trace, Task 0) records per-node scalars (h,
% phase, seed_kind/seed_from, term_reason, ok_final = the checker's accepted verdict, K0,
% D_uni, resid_static) and per-OUTER-ITERATION diagnostics (Dq_min/max, closest-mode flat
% indices) -- but it does NOT export the full per-node (Sigma(:), K(:), lambda(:), Gstat)
% state Sec. C needs (by design: those are Matsubara-length arrays, not part of the Task-0
% trace schema, and this task's scope fence forbids adding them there). The ONLY way to
% recover that full state without touching invz_hmf_ordered.m is to re-solve each node
% ourselves via the SAME shared, pure, additive Task-1 solver (invz_ordered_node_solve),
% fed a `node` struct built by the SAME recipe eval_node uses (task2_build_node below,
% copied statement-for-statement from invz_hmf_ordered.m's eval_node, :276-325, and
% cross-checked against test_invz_ordered_node_solve.m / test_invz_ordered_residual.m's own
% independent fixture reconstruction of the same recipe). Because invz_ordered_node_solve
% is a PURE, deterministic function of (node, seed, sopts), replaying it on the identical
% (node, seed) pair invz_hmf_ordered's eval_node used internally reproduces that node's
% state bit-for-bit -- this is verified, not merely assumed: every swept-mode node's
% replayed (accepted, D_uni) is cross-checked against the production trace's own recorded
% (ok_final, D_uni) at the SAME node (see run_swept_mode's `replay_ok`); a mismatch is
% surfaced via warning('invz:task2ReplayMismatch', ...) and out.summary.replay_mismatch_any,
% never silently absorbed. test_invz_task2_harness.m asserts replay_mismatch_any is false
% on both validation configs.
%
% cfg fields:
%   .ion         (default invz_ion())
%   .T, .Bx      REQUIRED. Bx: scalar or [Bx By Bz] (invz_field_vec convention).
%   .couplings   struct, one of:
%                  synthetic: struct('Jnu_flat', vec, 'J0eff', scalar[, 'Jxx0', scalar])
%                  real lattice (default when omitted): struct('grid',[16 16 16],
%                    'dpRng', 30) -- resolved via invz_bz_couplings/invz_jq_modes, exactly
%                    as invz_run_phase_diagram.m / test_invz_ordered_trace.m do.
%   .J0eff, .Jxx0  optional overrides (else taken from .couplings / ion defaults).
%   .solve_opts  struct forwarded to the production numerics (Jxx0, hyp, transverse_mf,
%                mix_outer, tol_outer, max_outer, Ecut, emt, emt_static, nH, hmax_fac,
%                hmin_frac, tol_root, force_bare, mf_maxit, mf_mix) -- production defaults
%                apply when omitted (mirrored here, never re-tuned).
%   .mode        'swept' (default) | 'isolated'.
%     'swept':   drives invz_ordered_trace (-> invz_hmf_ordered, tracing on): the FULL
%                H_MF profile (predictor/sweep/extend/redensify/bisect/root), exactly as
%                the production ordered leg does. One record per eval_node call.
%   .h_list      (isolated only, REQUIRED) field value(s) hz_fixed to solve, independently.
%   .seed_list   (isolated only, default {[]}, cell array of seeds: [] = cold, or
%                struct('Sigma',...,'K0s',...) = warm) -- one record per (h, seed) pair.
%   .label       free-form tag, stored in out.meta.label.
%
% out.meta     : run-level (T, Bx, J0eff, Jxx0, mode, label, grid, dpRng, is_synthetic,
%                solve_opts, Jnu_hash).
% out.nodes    : struct array, one per node (see task2_empty_record for the full field
%                list): id, h, phase, seed_kind, seed_from, accepted, term_reason,
%                outer_iters, D_uni, Dq_min, Dq_max, s, class, state (Sigma/K/lambda/K0/
%                Gstat/D_uni/hstar/m_star/r_star), idx_pos_flat/Dq_pos_val/idx_neg_flat/
%                Dq_neg_val, qbranch_pos/qbranch_neg (q_idx/branch_idx/Jq; NaN on synthetic
%                couplings), r, m, replay_ok.
% out.summary  : per-class counts, n_accepted, hstar, hmf_status, replay_mismatch_any.
% out.extra    : raw invz_ordered_trace output (swept mode) for provenance/debugging.
%
% Error policy (review-fix pass, stage-2c task 2a): the per-config node-production calls
% (the switch(mode) block below -- run_swept_mode / run_isolated_mode, i.e.
% invz_ordered_trace -> invz_hmf_ordered, and the task2_build_node -> invz_ordered_node_solve
% replay chain) are wrapped in ONE try/catch. An identifier matching 'invz:*' but NOT
% 'invz:task2*' -- a genuine production physics/numerics exception, e.g.
% invz:degenerateDoublet, invz:transverseMF, invz:hzFixed, invz:orderedTrace,
% invz:nodeSolveNode (none currently thrown by the two Task-2a validation configs) -- is
% ABSORBED into a structured failed-config record: out.status = 'failed', out.err_id/
% out.err_msg carry the caught identifier/message, out.nodes = struct([]) (empty, this
% file's existing "no nodes" convention -- see run_swept_mode/run_isolated_mode), and
% out.summary reflects zero counts with hmf_status = 'error' -- returned NORMALLY, so a
% caller (Task 2b's matrix loop) continues to the next cell instead of crashing. Anything
% else RETHROWS unchanged: this file's OWN 'invz:task2Config'/'invz:task2Node'
% argument-validation errors, this deliverable's sibling checker files' own 'invz:task2*'
% errors (invz_task2_classify/invz_task2_agree/invz_task2_ladder_ok -- a harness programming
% defect, not a physics signal, see task2_is_absorbable), and any non-'invz' MATLAB error.
% The two validation configs never throw, so their output is byte-for-behaviour unchanged by
% this policy.
if nargin < 1 || isempty(cfg), cfg = struct(); end
ion = getf(cfg, 'ion', invz_ion());
if ~isfield(cfg, 'T'), error('invz:task2Config', 'cfg.T is required.'); end
if ~isfield(cfg, 'Bx'), error('invz:task2Config', 'cfg.Bx is required.'); end
T  = cfg.T;
Bx = invz_field_vec(cfg.Bx);
mode  = getf(cfg, 'mode', 'swept');
label = getf(cfg, 'label', '');

[Jnu_flat, J0eff, Jxx0, qc, Jnu_unflat, latinfo, grid_, dpRng_, is_synth] = ...
    resolve_couplings(ion, cfg);

cfg_opts = getf(cfg, 'solve_opts', struct());
cfg_opts.J0eff = J0eff;                                  % required downstream, no default
if ~isfield(cfg_opts, 'Jxx0'), cfg_opts.Jxx0 = Jxx0; end
ro = resolve_numerics(ion, T, cfg_opts);

failed = false;  err_id = '';  err_msg = '';
try
    switch mode
        case 'swept'
            [nodes, extra] = run_swept_mode(ion, T, Bx, Jnu_flat, cfg_opts, ro, qc, ...
                Jnu_unflat, latinfo, grid_, dpRng_);
        case 'isolated'
            [nodes, extra] = run_isolated_mode(ion, T, Bx, Jnu_flat, cfg, ro, Jnu_unflat, is_synth);
        otherwise
            error('invz:task2Config', 'cfg.mode must be ''swept'' or ''isolated''; got ''%s''.', mode);
    end
catch ME
    if task2_is_absorbable(ME.identifier)
        failed  = true;
        err_id  = ME.identifier;
        err_msg = ME.message;
        nodes = struct([]);
        extra = struct('hstar', NaN, 'hmf_status', 'error', 'replay_mismatch_any', false, 'trc', []);
    else
        rethrow(ME);
    end
end

out.meta = struct('T', T, 'Bx', Bx, 'J0eff', J0eff, 'Jxx0', Jxx0, 'mode', mode, ...
    'label', label, 'grid', grid_, 'dpRng', dpRng_, 'is_synthetic', is_synth, ...
    'solve_opts', cfg_opts, 'Jnu_hash', weak_hash(Jnu_flat));
out.nodes   = nodes;
out.summary = build_summary(nodes, extra, mode);
out.extra   = extra;
if failed
    out.status  = 'failed';
    out.err_id  = err_id;
    out.err_msg = err_msg;
end
end

% =================================================================================================
function tf = task2_is_absorbable(id)
%TASK2_IS_ABSORBABLE True for an 'invz:*' identifier that is NOT this deliverable's own
% 'invz:task2*' namespace. 'invz:task2*' covers BOTH this file's own argument-validation
% errors (invz:task2Config/invz:task2Node, above and in task2_build_node) AND the sibling
% Task-2 checker files' own argument-validation errors (invz_task2_classify.m's
% invz:task2Classify, invz_task2_agree.m's invz:task2Agree, invz_task2_ladder_ok.m's
% invz:task2LadderOk) -- all of these signal a harness/programming defect (a malformed node/
% state bundle this file itself assembled), never a physics exception, so they always
% rethrow. Everything else under 'invz:*' is a genuine production physics/numerics exception
% (e.g. invz_twolevel_ordered.m/invz_twolevel_avg.m's invz:degenerateDoublet,
% invz_hmf_ordered.m/invz_single_ion.m/invz_solve_point_ordered.m's invz:transverseMF/
% invz:hzFixed, invz_ordered_trace.m's invz:orderedTrace, invz_ordered_node_solve.m's
% invz:nodeSolveNode) and is absorbed into a failed-config record by the caller's try/catch.
% Any non-'invz' MATLAB error (id = '' or some other component) is never absorbable either.
tf = startsWith(id, 'invz:') && ~startsWith(id, 'invz:task2');
end

% =================================================================================================
function [Jnu_flat, J0eff, Jxx0, qc, Jnu_unflat, latinfo, grid_, dpRng_, is_synth] = ...
    resolve_couplings(ion, cfg)
%RESOLVE_COUPLINGS Synthetic (explicit Jnu_flat) or real-lattice (invz_bz_couplings, the
% SAME shared BZ-grid helper invz_run_phase_diagram.m / test_invz_ordered_trace.m use).
cpl = getf(cfg, 'couplings', struct());
if isfield(cpl, 'Jnu_flat') && ~isempty(cpl.Jnu_flat)
    Jnu_flat = cpl.Jnu_flat(:);
    J0eff = getf(cpl, 'J0eff', getf(cfg, 'J0eff', ion.J0eff));
    Jxx0  = getf(cpl, 'Jxx0',  getf(cfg, 'Jxx0',  ion.Jxx0));
    qc = [];  Jnu_unflat = [];  latinfo = struct();  grid_ = [];  dpRng_ = NaN;
    is_synth = true;
else
    grid_  = getf(cpl, 'grid',  [16 16 16]);
    dpRng_ = getf(cpl, 'dpRng', 30);
    [Jnu_flat, info, Jaa0] = invz_bz_couplings(ion, struct('grid', grid_, 'dpRng', dpRng_));
    J0eff = getf(cfg, 'J0eff', info.Jcc0);
    Jxx0  = getf(cfg, 'Jxx0',  Jaa0);
    qc = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid_, 'range', [-0.5 0.5], 'verbose', false);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    [Jnu_unflat, info4] = invz_jq_modes(ion, qc, struct('dpRng', dpRng_, 'cache', true));
    latinfo = info4;
    is_synth = false;
end
end

% =================================================================================================
function ro = resolve_numerics(ion, T, cfg_opts)
%RESOLVE_NUMERICS Mirrors invz_hmf_ordered.m's own top-of-function option resolution
% (:44-64, :112-116) verbatim, so this harness's replay solves under IDENTICAL numerical
% knobs to whatever invz_ordered_trace/invz_hmf_ordered used for the same cfg_opts.
J0eff = cfg_opts.J0eff;
Jxx0  = getf(cfg_opts, 'Jxx0', ion.Jxx0);
hyp   = getf(cfg_opts, 'hyp', true);
tmf   = getf(cfg_opts, 'transverse_mf', 'legacy_x');
mixo  = getf(cfg_opts, 'mix_outer', 0.7);
tolo  = getf(cfg_opts, 'tol_outer', 1e-8);
maxo  = getf(cfg_opts, 'max_outer', 200);
eso   = getf(cfg_opts, 'emt_static', struct());
eso.warn = false;
sibase = struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf);
for f = {'mf_maxit', 'mf_mix'}
    if isfield(cfg_opts, f{1}), sibase.(f{1}) = cfg_opts.(f{1}); end
end
Ecut  = getf(cfg_opts, 'Ecut', 40);
eopts = getf(cfg_opts, 'emt', struct());
[wn, wts, beta] = invz_matsubara(T, Ecut);
ro = struct('J0eff', J0eff, 'Jxx0', Jxx0, 'hyp', hyp, 'tmf', tmf, 'mixo', mixo, ...
    'tolo', tolo, 'maxo', maxo, 'eso', eso, 'sibase', sibase, 'Ecut', Ecut, ...
    'eopts', eopts, 'wn', wn, 'wts', wts, 'beta', beta);
end

% =================================================================================================
function [node, mk, si] = task2_build_node(ion, T, Bx, hp, Jnu_flat, ro)
%TASK2_BUILD_NODE One fixed-field node's input bundle, mirroring invz_hmf_ordered.m's
% eval_node (:276-325) STATEMENT FOR STATEMENT (si/tl/c0/G0/G0inel0/X/fb/G0bare0/G0el0/g,
% the SAME production helper calls in the SAME order/arguments) -- the identical recipe
% test_invz_ordered_node_solve.m / test_invz_ordered_residual.m's own build_fixture()
% reconstructs from a converged pt (there, from pt.si/pt.tl; here, freshly at an arbitrary
% hz_fixed = hp, exactly as eval_node itself does for every profile/bisection/root node).
sio = ro.sibase;  sio.hz_fixed = hp;
si = invz_single_ion(ion, T, Bx, sio);
if abs(si.hz - hp) > 1e-12
    error('invz:task2Node', 'hz_fixed not held at h=%.6g (si.hz=%.6g).', hp, si.hz);
end
tl = invz_twolevel_ordered(ion, T, Bx, hp, struct('Jxx0', ro.Jxx0, 'transverse_mf', ro.tmf));
mk = si.Jexp(3);
c0  = invz_chi0z(si, T, 1i*ro.wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*ro.wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X = real(c0(:, :, 1));
switch ro.tmf
    case 'none'
        fb = 0;
    case 'legacy_x'
        fb = X(3,1) * (ro.Jxx0 / (1 - ro.Jxx0*X(1,1))) * X(1,3);
    case 'vector_ab'
        t = [1 2];
        fb = X(3,t) * (ro.Jxx0 * ((eye(2) - ro.Jxx0*X(t,t)) \ X(t,3)));
    otherwise
        error('invz:task2Node', 'unknown transverse_mf ''%s''.', ro.tmf);
end
G0bare0 = -(X(3,3) + fb);
G0el0   = G0bare0 - G0inel0;
g = real(invz_g(tl, 1i*ro.wn));
node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', ro.wts, 'wn', ro.wn, 'beta', ro.beta, ...
    'J0eff', ro.J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
    'eso', ro.eso, 'eopts', ro.eopts, 'Jnu_flat', Jnu_flat);
end

% =================================================================================================
function [nodes, extra] = run_swept_mode(ion, T, Bx, Jnu_flat, cfg_opts, ro, qc, ...
    Jnu_unflat, latinfo, grid_, dpRng_)
%RUN_SWEPT_MODE Drives invz_ordered_trace (-> invz_hmf_ordered, tracing on) for the
% GROUND-TRUTH node sequence (h, phase, seed_kind/seed_from, ok_final = accepted, D_uni,
% hstar/hmf_status), then REPLAYS every node via task2_build_node + invz_ordered_node_solve
% (warm/cold threaded EXACTLY per eval_node's own P0-3 rollback rule: warm from the last
% ACCEPTED node's state, rollback to the prior carried state on a non-accepted node -- never
% the failed node's own iterate) to recover the full (Sigma, K, lambda, Gstat) state and the
% complete Dq array. Cross-checks every replayed node's (accepted, D_uni) against the
% production trace's own recorded (ok_final, D_uni) -- see the file header.
opts_trace = struct('qc', qc, 'Jnu_unflat', Jnu_unflat, 'grid', grid_, 'dpRng', dpRng_, ...
    'lattice_info', latinfo, 'solve', cfg_opts);
trc = invz_ordered_trace(ion, T, Bx, Jnu_flat, opts_trace);

n = numel(trc.nodes);
recs = cell(1, n);
carried = [];             % carried (Sigma, K0s); [] = cold (matches eval_node's own state)
mismatch_any = false;
sopts = struct('mix_outer', ro.mixo, 'max_outer', ro.maxo, 'tol_outer', ro.tolo, ...
    'cold_retry', true, 'trace', false);

for k = 1:n
    ntr = trc.nodes(k);
    if strcmp(ntr.term_reason, 'bare_shortcut')
        % force_bare caller-level shortcut (rk=1, ok=true, NO Sigma/K solve at all --
        % invz_hmf_ordered.m header note); not exercised by either Task-2a validation
        % config, handled defensively. D_uni is undefined (NaN) -> prereg Sec. A's
        % "or non-finite" clause correctly forces 'unconverged' regardless of ok=true.
        rec = task2_empty_record();
        rec.id = ntr.id;  rec.h = ntr.h;  rec.phase = ntr.phase;
        rec.seed_kind = ntr.seed_kind;  rec.seed_from = ntr.seed_from;
        rec.accepted = true;  rec.term_reason = 'bare_shortcut';  rec.outer_iters = 0;
        [rec.class, rec.s, rec.D_uni, rec.Dq_min] = invz_task2_classify( ...
            struct('accepted', true, 'D_uni', NaN, 'Dq_min', NaN));
        rec.r = 1;  rec.replay_ok = true;
        recs{k} = rec;
        continue;   % carried UNCHANGED: eval_node's fbare branch returns before touching Sigma/K0s
    end

    [nd, mk, ~] = task2_build_node(ion, T, Bx, ntr.h, Jnu_flat, ro);
    seed = carried;
    [state, info] = invz_ordered_node_solve(nd, seed, sopts);

    if info.accepted
        carried = struct('Sigma', state.Sigma, 'K0s', state.K0s);
    end   % else: rollback -- carried stays whatever it already was (eval_node's P0-3 rule)

    Dcheck = (isnan(info.so.D_uni) && isnan(ntr.D_uni)) || ...
             (isfinite(info.so.D_uni) && isfinite(ntr.D_uni) && ...
              abs(info.so.D_uni - ntr.D_uni) <= 1e-9);
    replay_ok = isequaln(info.accepted, ntr.ok_final) && Dcheck;
    if ~replay_ok
        mismatch_any = true;
        warning('invz:task2ReplayMismatch', ['node id %d (h=%.6g, phase=%s): harness replay ' ...
            'disagrees with the production trace (accepted %d vs %d, D_uni %.6g vs %.6g) -- ' ...
            'see invz_task2_run_config.m header.'], ntr.id, ntr.h, ntr.phase, ...
            info.accepted, ntr.ok_final, info.so.D_uni, ntr.D_uni);
    end

    Dq = 1 + (Jnu_flat(:) - state.K0s) .* info.so.Gstat;
    [class_, s_, D_uni_, Dq_min_] = invz_task2_classify(struct('accepted', info.accepted, ...
        'D_uni', info.so.D_uni, 'Dq', Dq));
    [idx_pos, Dq_pos, idx_neg, Dq_neg] = closest_modes(Dq);
    [qbp, qbn] = resolve_qbranch(trc.meta, idx_pos, idx_neg);

    rec = task2_empty_record();
    rec.id = ntr.id;  rec.h = ntr.h;  rec.phase = ntr.phase;
    rec.seed_kind = info.seed_kind;  rec.seed_from = ntr.seed_from;
    rec.accepted = info.accepted;  rec.term_reason = info.term_reason;
    rec.outer_iters = info.outer_iters;
    rec.D_uni = D_uni_;  rec.Dq_min = Dq_min_;
    if any(~isfinite(Dq)), rec.Dq_max = NaN; else, rec.Dq_max = max(Dq); end   % NaN-propagating,
                                                       % matching invz_task2_classify's robust_min
    rec.s = s_;  rec.class = class_;
    rec.state = struct('Sigma', state.Sigma, 'K', state.K, 'lambda', state.lam, ...
        'K0', state.K0s, 'Gstat', info.so.Gstat, 'D_uni', info.so.D_uni, ...
        'hstar', NaN, 'm_star', NaN, 'r_star', NaN);
    rec.idx_pos_flat = idx_pos;  rec.Dq_pos_val = Dq_pos;
    rec.idx_neg_flat = idx_neg;  rec.Dq_neg_val = Dq_neg;
    rec.qbranch_pos = qbp;  rec.qbranch_neg = qbn;
    rec.r = info.so.r;  rec.m = mk;
    rec.replay_ok = replay_ok;
    recs{k} = rec;
end
if isempty(recs), nodes = struct([]); else, nodes = [recs{:}]; end

% Root-node endpoint quantities (hstar/m_star/r_star, "where applicable" -- brief
% deliverable 2): populated ONLY on the phase=='root' node, matching production's own
% prof.m_star/D_uni_star/r_star/Gstat_star (invz_hmf_ordered.m:247) -- here independently
% recomputed from this harness's own replay rather than reused from prof (which
% invz_ordered_trace does not forward).
if ~isempty(nodes)
    root_ix = find(strcmp({nodes.phase}, 'root'), 1);
    if ~isempty(root_ix) && nodes(root_ix).accepted
        nodes(root_ix).state.hstar  = trc.result.hstar;
        nodes(root_ix).state.m_star = nodes(root_ix).m;
        nodes(root_ix).state.r_star = nodes(root_ix).r;
    end
end
extra = struct('hstar', trc.result.hstar, 'hmf_status', trc.result.hmf_status, ...
    'replay_mismatch_any', mismatch_any, 'trc', trc);
end

% =================================================================================================
function [nodes, extra] = run_isolated_mode(ion, T, Bx, Jnu_flat, cfg, ro, Jnu_unflat, is_synth)
%RUN_ISOLATED_MODE Solves cfg.h_list independently (no continuation): one record per
% (h, seed) pair in cfg.seed_list (default {[]}, i.e. cold only). Task-2 experiment-matrix
% cell 1 ("isolated vs swept ... node converges alone but not in continuation") re-solves
% the SAME h a swept run visited, isolated/cold, via exactly this path.
h_list = getf(cfg, 'h_list', []);
if isempty(h_list)
    error('invz:task2Config', 'cfg.mode=''isolated'' requires cfg.h_list (>=1 field value).');
end
seed_list = getf(cfg, 'seed_list', {[]});
if ~iscell(seed_list), seed_list = {seed_list}; end

nq = NaN;
if ~is_synth, nq = size(Jnu_unflat, 1); end
meta = struct('is_synthetic', is_synth, 'Jnu_unflat', Jnu_unflat, 'nq', nq);

sopts = struct('mix_outer', ro.mixo, 'max_outer', ro.maxo, 'tol_outer', ro.tolo, ...
    'cold_retry', true, 'trace', false);
recs = {};
node_id = 0;
for ih = 1:numel(h_list)
    hp = h_list(ih);
    [nd, mk, ~] = task2_build_node(ion, T, Bx, hp, Jnu_flat, ro);
    for is = 1:numel(seed_list)
        node_id = node_id + 1;
        seed = seed_list{is};
        [state, info] = invz_ordered_node_solve(nd, seed, sopts);
        Dq = 1 + (Jnu_flat(:) - state.K0s) .* info.so.Gstat;
        [class_, s_, D_uni_, Dq_min_] = invz_task2_classify(struct('accepted', info.accepted, ...
            'D_uni', info.so.D_uni, 'Dq', Dq));
        [idx_pos, Dq_pos, idx_neg, Dq_neg] = closest_modes(Dq);
        [qbp, qbn] = resolve_qbranch(meta, idx_pos, idx_neg);

        rec = task2_empty_record();
        rec.id = node_id;  rec.h = hp;  rec.phase = 'isolated';
        rec.seed_kind = info.seed_kind;  rec.seed_from = NaN;
        rec.accepted = info.accepted;  rec.term_reason = info.term_reason;
        rec.outer_iters = info.outer_iters;
        rec.D_uni = D_uni_;  rec.Dq_min = Dq_min_;
        if any(~isfinite(Dq)), rec.Dq_max = NaN; else, rec.Dq_max = max(Dq); end   % NaN-propagating,
                                                       % matching invz_task2_classify's robust_min
        rec.s = s_;  rec.class = class_;
        rec.state = struct('Sigma', state.Sigma, 'K', state.K, 'lambda', state.lam, ...
            'K0', state.K0s, 'Gstat', info.so.Gstat, 'D_uni', info.so.D_uni, ...
            'hstar', NaN, 'm_star', NaN, 'r_star', NaN);
        rec.idx_pos_flat = idx_pos;  rec.Dq_pos_val = Dq_pos;
        rec.idx_neg_flat = idx_neg;  rec.Dq_neg_val = Dq_neg;
        rec.qbranch_pos = qbp;  rec.qbranch_neg = qbn;
        rec.r = info.so.r;  rec.m = mk;
        rec.replay_ok = true;   % no production ground truth to cross-check in isolated mode
        recs{end+1} = rec; %#ok<AGROW>
    end
end
if isempty(recs), nodes = struct([]); else, nodes = [recs{:}]; end
extra = struct('hstar', NaN, 'hmf_status', 'isolated', 'replay_mismatch_any', false, 'trc', []);
end

% =================================================================================================
function rec = task2_empty_record()
%TASK2_EMPTY_RECORD Full field template (ensures every node record -- across bare_shortcut/
% normal/isolated branches -- shares an identical field set for struct-array concatenation).
rec = struct('id', NaN, 'h', NaN, 'phase', '', 'seed_kind', '', 'seed_from', NaN, ...
    'accepted', false, 'term_reason', '', 'outer_iters', NaN, ...
    'D_uni', NaN, 'Dq_min', NaN, 'Dq_max', NaN, 's', NaN, 'class', 'unconverged', ...
    'state', struct('Sigma', [], 'K', [], 'lambda', [], 'K0', NaN, 'Gstat', NaN, ...
        'D_uni', NaN, 'hstar', NaN, 'm_star', NaN, 'r_star', NaN), ...
    'idx_pos_flat', NaN, 'Dq_pos_val', NaN, 'idx_neg_flat', NaN, 'Dq_neg_val', NaN, ...
    'qbranch_pos', struct('q_idx', NaN, 'branch_idx', NaN, 'Jq', NaN), ...
    'qbranch_neg', struct('q_idx', NaN, 'branch_idx', NaN, 'Jq', NaN), ...
    'r', NaN, 'm', NaN, 'replay_ok', true);
end

% =================================================================================================
function [idx_pos, Dq_pos, idx_neg, Dq_neg] = closest_modes(Dq)
%CLOSEST_MODES Closest-to-zero positive/negative Dq mode (flat index + value); mirrors
% invz_ordered_node_solve.m's append_iter convention exactly.
Dq = Dq(:);
posmask = Dq > 0;
if any(posmask)
    ix = find(posmask);  [vmin, jm] = min(Dq(ix));  idx_pos = ix(jm);  Dq_pos = vmin;
else
    idx_pos = NaN;  Dq_pos = NaN;
end
negmask = ~posmask;
if any(negmask)
    ix = find(negmask);  [vmax, jm] = max(Dq(ix));  idx_neg = ix(jm);  Dq_neg = vmax;
else
    idx_neg = NaN;  Dq_neg = NaN;
end
end

% =================================================================================================
function [qbp, qbn] = resolve_qbranch(meta, idx_pos, idx_neg)
%RESOLVE_QBRANCH q/branch provenance of the closest +/-Dq modes via invz_ordered_trace_
% resolve.m (real-lattice grids only; NaN triple on synthetic couplings/missing mode).
qbp = struct('q_idx', NaN, 'branch_idx', NaN, 'Jq', NaN);
qbn = qbp;
if meta.is_synthetic, return; end
if isfinite(idx_pos)
    [q, b, J] = invz_ordered_trace_resolve(meta, idx_pos);
    qbp = struct('q_idx', q, 'branch_idx', b, 'Jq', J);
end
if isfinite(idx_neg)
    [q, b, J] = invz_ordered_trace_resolve(meta, idx_neg);
    qbn = struct('q_idx', q, 'branch_idx', b, 'Jq', J);
end
end

% =================================================================================================
function summary = build_summary(nodes, extra, mode)
%BUILD_SUMMARY Per-class node counts + config-level verdict.
summary = struct('n_nodes', numel(nodes), 'hstar', extra.hstar, 'hmf_status', extra.hmf_status, ...
    'replay_mismatch_any', extra.replay_mismatch_any, 'mode', mode, ...
    'n_stable', 0, 'n_marginal', 0, 'n_unstable', 0, 'n_unconverged', 0, 'n_accepted', 0);
if isempty(nodes), return; end
cls = {nodes.class};
summary.n_stable      = nnz(strcmp(cls, 'stable'));
summary.n_marginal    = nnz(strcmp(cls, 'marginal'));
summary.n_unstable    = nnz(strcmp(cls, 'unstable'));
summary.n_unconverged = nnz(strcmp(cls, 'unconverged'));
summary.n_accepted    = nnz([nodes.accepted]);
end

% =================================================================================================
function h = weak_hash(v)
%WEAK_HASH Deterministic run-fingerprint digest (mirrors invz_ordered_trace.m's own local
% weak_hash / invz_common/invz_cache_key.m's hash_vec convention) -- provenance only.
v = v(:);
h = sprintf('%dv_%08x', numel(v), typecast(single(sum(v.*(1:numel(v)).')), 'uint32'));
end
