function tests = test_invz_ordered_trace
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_production_trace_physical_defect_and_synthetic_control(testCase)
% Stage-2c task 0 pinning test. A versioned trace of the PRODUCTION jensen ordered-leg
% node loop (invz_hmf_ordered.m's opts.trace hook, via invz_ordered_trace.m) freezes:
%   (1) the physical-coupling masking defect on the real 16^3-grid BZ couplings, with the
%       node-failure / negative-Dq signal Task 2 will classify;
%   (2) a synthetic control where the SAME jensen machinery is accepted; and
%   (3) the behaviour-neutrality the trace hook depends on (traced-off reproduces the
%       untraced production result bit-for-bit, both through invz_solve_point_ordered and
%       directly at the invz_hmf_ordered nargout boundary).
% INVZ_SLOW-gated: the physical leg alone measures ~40 s (16^3 grid, 34 fixed-field nodes,
% T=0.1 K/B=2.85 T -- see .superpowers/sdd/task0-report.md), well above every other test in
% this file's suite; this repo already gates comparable real-grid/jensen costs identically
% (test_invz_ordered_phase.m, test_invz_qcp_closure.m).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (16^3-grid jensen solve).');
ion = invz_ion();

% =========================================================================================
% (1) PHYSICAL FAILING CASE: real BZ couplings mask the whole jensen ordered leg (memory:
% invzp-jensen-realcoupling-nonconvergence). T=0.1, B=[2.85 0 0], invz_bz_couplings 16^3/30.
% =========================================================================================
Tp = 0.1;  Bp = [2.85 0 0];
gridOpts = struct('grid', [16 16 16], 'dpRng', 30);
[Jnu_p, info_p, Jaa0_p] = invz_bz_couplings(ion, gridOpts);

% Unflattened q/branch provenance, obtained by calling invz_jq_modes directly on the SAME
% qc/dpRng invz_bz_couplings used internally (stage2c-context.md key-anchor: invz_bz_couplings
% .m:17 does Jnu=Jnu(:) column-major over invz_jq_modes' [nq x 4] branches).
qc = qVec_generator(ion.a, 'mode', 'grid', 'grid', gridOpts.grid, 'range', [-0.5 0.5], 'verbose', false);
qc = qc(any(abs(qc) > 1e-12, 2), :);
[Jnu4_p, info4_p] = invz_jq_modes(ion, qc, struct('dpRng', gridOpts.dpRng, 'cache', true));
verifyEqual(testCase, Jnu_p, Jnu4_p(:));   % same cached coupling set (exact: plain reshape)

opts_phys = struct('qc', qc, 'Jnu_unflat', Jnu4_p, 'grid', gridOpts.grid, 'dpRng', gridOpts.dpRng, ...
    'lattice_info', info4_p, 'solve', struct('J0eff', info_p.Jcc0, 'Jxx0', Jaa0_p, 'hyp', true));
trc_phys = invz_ordered_trace(ion, Tp, Bp, Jnu_p, opts_phys);

% masks: no accepted root, honest non-'ok' status (never silently PM/converged)
verifyTrue(testCase, ~isfinite(trc_phys.result.hstar));
verifyTrue(testCase, ~strcmp(trc_phys.result.hmf_status, 'ok'));
verifyTrue(testCase, ~trc_phys.meta.is_synthetic);
verifyEqual(testCase, trc_phys.meta.nq, size(Jnu4_p, 1));

% the trace captured the D_uni<0 / node-failure pattern (stage2c-context.md's binding
% physics: Dq = D_uni + (J0eff-J(q))*chi_stat, so Dq >= D_uni EVERYWHERE when chi_stat =
% -Gstat > 0 and max(Jnu) <= J0eff hold -- the GENERIC/near-convergence regime. This is a
% precondition, not a universal invariant: a minority of early, not-yet-converged outer
% iterations (checked separately, not asserted away) can carry a transiently POSITIVE
% Gstat before a node's own Sigma<->K loop settles, so this is pinned as an "at least one"
% co-occurrence below, not a blanket per-iteration identity.
verifyGreaterThanOrEqual(testCase, numel(trc_phys.nodes), 1);
verifyGreaterThanOrEqual(testCase, numel(trc_phys.iters), 1);
Duni_all = [trc_phys.iters.D_uni];  Dqmin_all = [trc_phys.iters.Dq_min];
neg_iters = [trc_phys.iters.Dq_neg_count] > 0;
verifyGreaterThanOrEqual(testCase, nnz(neg_iters), 1, 'expected >= 1 iteration with a negative Dq mode');
verifyLessThan(testCase, min(Dqmin_all(neg_iters)), 0);
verifyTrue(testCase, any(Duni_all(neg_iters) < 0), ...
    'expected >= 1 negative-Dq iteration to also show D_uni < 0 (the unstable-branch signature)');
% at least one node genuinely failed to close (never silently exported as 'ok')
verifyGreaterThanOrEqual(testCase, nnz(~[trc_phys.nodes.ok_final]), 1);
% Closed set of reachable per-node term_reason values (task 12 fix round 1: made
% 'medium_out_of_domain' reachable -- a strict-scheme reference/closure domain event,
% local_term_reason's own pass-through case -- so it belongs in this set alongside the
% four this file's own vocabulary defines. Still a closed set, not a vacuous check.
verifyTrue(testCase, all(ismember({trc_phys.nodes.term_reason}, ...
    {'converged', 'max_iter', 'refresh_failed', 'bare_shortcut', 'medium_out_of_domain'})));

% resolve the closest-to-zero negative Dq mode back to (q, branch) via invz_ordered_trace_
% resolve.m's flat-index formula (stage2c-context.md: q=mod(k-1,nq)+1, branch=floor((k-1)/nq)+1)
neg_with_idx = trc_phys.iters(neg_iters);
neg_with_idx = neg_with_idx(~isnan([neg_with_idx.idx_neg_flat]));
verifyGreaterThanOrEqual(testCase, numel(neg_with_idx), 1);
sample = neg_with_idx(1);
[q_idx, branch_idx, Jq] = invz_ordered_trace_resolve(trc_phys.meta, sample.idx_neg_flat);
verifyGreaterThanOrEqual(testCase, q_idx, 1);
verifyLessThanOrEqual(testCase, q_idx, trc_phys.meta.nq);
verifyTrue(testCase, ismember(branch_idx, 1:4));
verifyEqual(testCase, Jq, trc_phys.meta.Jnu_unflat(q_idx, branch_idx));

% =========================================================================================
% (2) SYNTHETIC PASSING CONTROL: same jensen machinery, ACCEPTED (repo's own 24-pt Jnu).
% =========================================================================================
Ts = 0.31;  Bs = [2.85 0 0];
Jnu_s = linspace(-2e-3, 6e-3, 24).';
J0eff_s = 6.42e-3;                    % = Jcc0 (brief-pinned value)
Jxx0_s  = ion.Jxx0;
solve_s = struct('J0eff', J0eff_s, 'Jxx0', Jxx0_s, 'hyp', true);

pt_s = invz_solve_point_ordered(ion, Ts, Bs, Jnu_s, ...
    struct('J0eff', J0eff_s, 'Jxx0', Jxx0_s, 'hyp', true, 'ordered_mode', 'jensen'));
verifyTrue(testCase, pt_s.is_ordered);
verifyTrue(testCase, pt_s.converged);
verifyTrue(testCase, isfinite(pt_s.Sigma0));

trc_syn = invz_ordered_trace(ion, Ts, Bs, Jnu_s, struct('solve', solve_s));
verifyTrue(testCase, isequaln(trc_syn.result.hstar, pt_s.hmf));   % same deterministic root
verifyEqual(testCase, trc_syn.result.hmf_status, 'ok');
verifyGreaterThanOrEqual(testCase, numel(trc_syn.nodes), 1);
verifyTrue(testCase, all([trc_syn.nodes.ok_final]), ...
    'expected every node to converge on the synthetic control');
verifyTrue(testCase, all(strcmp({trc_syn.nodes.term_reason}, 'converged')));
verifyTrue(testCase, trc_syn.meta.is_synthetic);

% =========================================================================================
% (3) BEHAVIOUR-NEUTRALITY: traced-off invz_solve_point_ordered reproduces the untraced
% production result bit-for-bit on the synthetic control (opts.trace forwards through
% invz_solve_point_ordered's own hopts = opts copy for free -- see invz_hmf_ordered.m header);
% and directly at invz_hmf_ordered's nargout boundary (requesting the 3rd/trace output must
% not perturb hmf_star/prof relative to a 2-output call).
% =========================================================================================
base_opts   = struct('J0eff', J0eff_s, 'Jxx0', Jxx0_s, 'hyp', true, 'ordered_mode', 'jensen');
traced_opts = base_opts;  traced_opts.trace = true;
pt_off = invz_solve_point_ordered(ion, Ts, Bs, Jnu_s, base_opts);
pt_on  = invz_solve_point_ordered(ion, Ts, Bs, Jnu_s, traced_opts);
verifyTrue(testCase, isequaln(pt_on.is_ordered, pt_off.is_ordered));
verifyTrue(testCase, isequaln(pt_on.converged,  pt_off.converged));
verifyTrue(testCase, isequaln(pt_on.Sigma0,     pt_off.Sigma0));
verifyTrue(testCase, isequaln(pt_on.hmf,        pt_off.hmf));
verifyTrue(testCase, isequaln(pt_on.m0,         pt_off.m0));

hmf_base = rmfield(base_opts, 'ordered_mode');  hmf_base.J0eff = J0eff_s;
hmf_traced = hmf_base;  hmf_traced.trace = true;
[h2a, p2a]       = invz_hmf_ordered(ion, Ts, Bs, Jnu_s, hmf_base);
[h2b, p2b, trc3] = invz_hmf_ordered(ion, Ts, Bs, Jnu_s, hmf_traced);
verifyTrue(testCase, isequaln(h2a, h2b));
verifyTrue(testCase, isequaln(p2a, p2b));
verifyTrue(testCase, trc3.enabled);
[h2c, p2c, trc3_off] = invz_hmf_ordered(ion, Ts, Bs, Jnu_s, hmf_base); % 3rd output requested, trace OFF
verifyFalse(testCase, trc3_off.enabled);
verifyEqual(testCase, numel(trc3_off.nodes), 0);
verifyEqual(testCase, numel(trc3_off.iters), 0);
verifyTrue(testCase, isequaln(h2c, h2a));
verifyTrue(testCase, isequaln(p2c, p2a));

% --- save/load round trip (the schema IS the .mat struct, never printed text) --
save_path = [tempname() '.mat'];
invz_ordered_trace(ion, Ts, Bs, Jnu_s, struct('solve', solve_s, 'save_path', save_path));
verifyEqual(testCase, exist(save_path, 'file'), 2);
S = load(save_path);
verifyTrue(testCase, isfield(S, 'trace') && isequal(S.trace.schema_version, 1));
if exist(save_path, 'file'), delete(save_path); end
end
