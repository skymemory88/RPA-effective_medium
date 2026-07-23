function tests = test_invz_task2_harness
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% ===========================================================================================
% (1) PHYSICAL DEFECT ANCHOR (stage-2c task 2a, deliverable 4, bullet 1) -- THE single most
% important assertion in this file: the harness runs on the real 16^3-grid physical config
% (T=0.1, B=[2.85 0 0], one of the prereg Sec. E frozen four physical fields -- the "existing
% 2.85 T defect anchor"), the overall solve masks (hstar NaN, matching the pinned defect --
% memory: invzp-jensen-realcoupling-nonconvergence; test_invz_ordered_trace.m's own physical
% leg), and its nodes classify as a MIX including 'unconverged' -- with the ACCEPTED-ONLY rule
% (prereg Sec. A) holding for EVERY node: a non-accepted node is 'unconverged', NEVER
% 'unstable', no matter how negative its transient D_uni/Dq went. This is the user's Task-2
% amendment; it is pinned directly, not merely implied by other assertions.
% INVZ_SLOW-gated: ~40 s for the production profile (invz_ordered_trace) plus a comparable
% replay cost (this harness's own full-state recovery -- see invz_task2_run_config.m's
% header) -- ~1-2 min total, above every other test in this file's suite, matching this
% repo's existing convention for the same physical config (test_invz_ordered_trace.m,
% test_invz_ordered_phase.m, test_invz_qcp_closure.m all gate comparable real-grid/jensen
% costs identically).
% ===========================================================================================
function test_physical_anchor_unconverged_not_unstable(testCase)
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), ...
    'INVZ_SLOW=1 to run (16^3-grid jensen profile + full-state replay, ~1-2 min).');
ion = invz_ion();
cfg = struct('ion', ion, 'T', 0.1, 'Bx', 2.85, 'mode', 'swept', ...
    'couplings', struct('grid', [16 16 16], 'dpRng', 30), 'label', 'physical_defect_anchor');
out = invz_task2_run_config(cfg);

% masks: honest non-'ok' verdict, never silently PM/converged (invz_hmf_ordered's own contract)
verifyTrue(testCase, ~isfinite(out.summary.hstar));
verifyFalse(testCase, strcmp(out.summary.hmf_status, 'ok'));
verifyGreaterThanOrEqual(testCase, out.summary.n_nodes, 1);

accepted = [out.nodes.accepted];
classes  = {out.nodes.class};
verifyGreaterThanOrEqual(testCase, nnz(~accepted), 1, ...
    'expected at least one non-accepted (max_iter) node, matching the pinned physical defect');

% THE central assertion (prereg Sec. A; stage2c-context.md): every non-accepted node is
% 'unconverged', full stop -- regardless of how negative its transient D_uni/Dq_min is.
verifyTrue(testCase, all(strcmp(classes(~accepted), 'unconverged')), ...
    'every non-accepted node must classify unconverged, never unstable/stable/marginal');
verifyTrue(testCase, ~any(strcmp(classes, 'unstable') & ~accepted), ...
    'no node may be labelled ''unstable'' unless it is checker-accepted (task-2 amendment)');
verifyTrue(testCase, any(strcmp(classes, 'unconverged')));

% Harness self-consistency (invz_task2_run_config.m header): every replayed node's own
% (accepted, D_uni) must match the production trace's recorded values -- proves this
% harness's node reconstruction is not silently diverging from what invz_hmf_ordered
% actually computed internally.
verifyFalse(testCase, out.summary.replay_mismatch_any, ...
    'harness replay disagreed with the production trace on >=1 node (see warnings above)');
verifyTrue(testCase, all([out.nodes.replay_ok]));
end

% ===========================================================================================
% (2) SYNTHETIC CONTROL (deliverable 4, bullet 2): same jensen machinery, ACCEPTED (repo's
% own 24-pt Jnu fixture, T=0.31/B=[2.85 0 0]/J0eff=6.42e-3 -- test_invz_ordered_trace.m:89-93,
% test_invz_ordered_node_solve.m/test_invz_ordered_residual.m's shared pinned fixture). The
% root node is checker-accepted and classifies stable or marginal with a well-defined s.
% ===========================================================================================
function test_synthetic_control_accepted_and_stable(testCase)
ion = invz_ion();
T = 0.31;  Bx = 2.85;  Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
cfg = struct('ion', ion, 'T', T, 'Bx', Bx, 'mode', 'swept', ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'label', 'synthetic_control');
out = invz_task2_run_config(cfg);

verifyTrue(testCase, isfinite(out.summary.hstar));
verifyEqual(testCase, out.summary.hmf_status, 'ok');
verifyGreaterThanOrEqual(testCase, out.summary.n_accepted, 1);

root_ix = find(strcmp({out.nodes.phase}, 'root'), 1);
verifyTrue(testCase, ~isempty(root_ix), 'expected a root-phase node in the swept profile');
rootrec = out.nodes(root_ix);
verifyTrue(testCase, rootrec.accepted);
verifyTrue(testCase, ismember(rootrec.class, {'stable', 'marginal'}), ...
    sprintf('root node classified ''%s'' (s=%.3g)', rootrec.class, rootrec.s));
verifyTrue(testCase, isfinite(rootrec.s));
verifyTrue(testCase, isfinite(rootrec.state.D_uni));
verifyEqual(testCase, rootrec.state.hstar, out.summary.hstar);
verifyTrue(testCase, isfinite(rootrec.state.m_star) && isfinite(rootrec.state.r_star));

verifyFalse(testCase, out.summary.replay_mismatch_any);
verifyTrue(testCase, all([out.nodes.replay_ok]));

% Cross-check against invz_solve_point_ordered's OWN (COLD-started) state at the SAME hstar
% (invz_solve_point_ordered.m:214 comment: "COLD: this leg always starts from Sigma=0/K0s=0"):
% the harness's WARM-threaded swept root node and the production point-solver's COLD root
% state are two different seeds converging on the SAME node -- exactly the isolated(cold)-
% vs-swept(warm) comparison prereg Sec. C / experiment-matrix cell 2 wants, dog-fooded here
% on real (not hand-built) data via invz_task2_agree itself.
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, ...
    struct('J0eff', J0eff, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen'));
verifyTrue(testCase, pt.converged);
stateA = rootrec.state;
stateB = struct('Sigma', pt.Sigma, 'K', pt.K, 'lambda', pt.lambda, 'K0', pt.K(1), ...
    'Gstat', pt.G(1), 'D_uni', pt.D_uni, 'hstar', pt.hmf, 'm_star', pt.m0, ...
    'r_star', rootrec.state.r_star);   % pt has no standalone r_star field; not compared
[ok_cw, detail_cw] = invz_task2_agree(stateA, stateB, struct('J0eff', J0eff));
verifyTrue(testCase, ok_cw, sprintf('cold-vs-warm root disagreement, worst=%s', detail_cw.worst));
end

% ===========================================================================================
% (2b) ISOLATED MODE, cross-checked against invz_solve_point_ordered's own cold solve at the
% SAME h=hstar -- both are COLD at the identical node, so this is the strongest available
% correctness check of task2_build_node's node-reconstruction fidelity (not required by the
% brief's minimum bar, but cheap and directly de-risks the "full exported state" gate item).
% ===========================================================================================
function test_isolated_mode_cold_matches_solve_point_ordered(testCase)
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, ...
    struct('J0eff', J0eff, 'Jxx0', ion.Jxx0, 'hyp', true, 'ordered_mode', 'jensen'));
verifyTrue(testCase, pt.converged);

cfg = struct('ion', ion, 'T', T, 'Bx', Bx, 'mode', 'isolated', 'h_list', pt.hmf, ...
    'couplings', struct('Jnu_flat', Jnu, 'J0eff', J0eff), 'seed_list', {{[]}});
out = invz_task2_run_config(cfg);
verifyEqual(testCase, numel(out.nodes), 1);
rec = out.nodes(1);
verifyTrue(testCase, rec.accepted);
verifyEqual(testCase, rec.seed_kind, 'cold');
verifyEqual(testCase, rec.state.Sigma, pt.Sigma,  'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, rec.state.K,     pt.K,      'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, rec.state.lambda, pt.lambda,'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, rec.state.D_uni, pt.D_uni,  'AbsTol', 1e-8, 'RelTol', 1e-6);
end

% ===========================================================================================
% (3) AGREEMENT HELPER unit checks (deliverable 4, bullet 3): identical states agree; a
% large perturbation disagrees; the AbsTol_q floor is proven load-bearing near D_uni = 0.
% ===========================================================================================
function test_agree_identical_states_agree(testCase)
J0eff = 6.42e-3;
s = struct('Sigma', [0.1; 0.2], 'K', [1e-3; 2e-3], 'lambda', [0.01; 0.02; 0.03], ...
    'K0', 1e-3, 'Gstat', -300, 'D_uni', 0.05, 'hstar', 0.02, 'm_star', 3.5, 'r_star', 0.8);
[ok, detail] = invz_task2_agree(s, s, struct('J0eff', J0eff));
verifyTrue(testCase, ok);
verifyTrue(testCase, detail.ok);
verifyEqual(testCase, detail.worst, '');
end

function test_agree_large_perturbation_disagrees(testCase)
J0eff = 6.42e-3;
a = struct('Sigma', [0.1; 0.2], 'K', [1e-3; 2e-3], 'lambda', [0.01; 0.02; 0.03], 'K0', 1e-3, ...
    'Gstat', -300, 'D_uni', 0.05, 'hstar', 0.02, 'm_star', 3.5, 'r_star', 0.8);
b = a;  b.K0 = a.K0 + 10*abs(J0eff);   % >> AbsTol_q + RelTol*scale_q -- genuine disagreement
[ok, detail] = invz_task2_agree(a, b, struct('J0eff', J0eff));
verifyFalse(testCase, ok);
verifyFalse(testCase, detail.K0.pass);
verifyEqual(testCase, detail.worst, 'K0');
end

% (3c) THE AbsTol_q FLOOR IS LOAD-BEARING (prereg Sec. C: "relative tolerance ALONE is
% invalid near D_uni = 0"). The standard reason a pure-relative bound is unsafe near zero is
% that it COLLAPSES toward zero exactly where D_uni legitimately lives near a stability
% boundary, so it would WRONGLY FLAG two re-solves of the SAME fixed point -- differing only
% by ordinary tol_outer=1e-8-level solver noise -- as disagreeing (the same caveat numpy's
% own isclose documents for rtol near zero). Constructed here: two D_uni values of O(1e-9),
% i.e. well inside the outer solver's own precision, differing by 6e-9 (< AbsTol_q=1e-8*1,
% so the COMBINED test correctly calls them the SAME state) but by MORE than the (collapsed)
% pure-relative bound RelTol*max(|a|,|b|) ~ 3e-15 (so RelTol ALONE would flag a spurious
% disagreement). This proves the AbsTol_q term changes the verdict -- it is not decorative.
function test_agree_abstol_floor_is_load_bearing_near_Duni_zero(testCase)
J0eff = 6.42e-3;
base  = struct('Sigma', 0.1, 'K', 1e-3, 'lambda', [0.01; 0.02; 0.03], 'K0', 1e-3, ...
    'Gstat', -300, 'D_uni', 3e-9, 'hstar', 0.02, 'm_star', 3.5, 'r_star', 0.8);
other = base;  other.D_uni = -3e-9;    % |a-b| = 6e-9: same fixed point up to solver noise

[ok, detail] = invz_task2_agree(base, other, struct('J0eff', J0eff));
verifyTrue(testCase, ok, 'the combined test must recognize this as the SAME state');
verifyTrue(testCase, detail.D_uni.pass);

% Explicit, load-bearing proof: the SAME pair, under RelTol alone (AbsTol_q forced to 0),
% must give the OPPOSITE verdict for D_uni -- i.e. deleting the floor changes the outcome.
a = base.D_uni;  b = other.D_uni;
relOnlyTol = 1e-6 * max(abs(a), abs(b));           % RelTol*max(|a|,|b|), NO AbsTol_q floor
verifyGreaterThan(testCase, abs(a - b), relOnlyTol, ...
    'RelTol alone must be strictly tighter than |a-b| here, or the floor would not be load-bearing');
scale_Duni  = 1;                                    % invz_task2_agree.m: scale_q(D_uni) = 1
combinedTol = 1e-8*scale_Duni + relOnlyTol;
verifyLessThanOrEqual(testCase, abs(a - b), combinedTol, ...
    'the documented AbsTol_q + RelTol*max(.) bound must accept this pair');
end

% ===========================================================================================
% (4) BONUS: light, cheap, pure-numeric self-checks of the classifier (accepted-only rule
% and the +-1e-3 thresholds in isolation) and the mesh-ladder helper -- not part of the
% brief's minimum bar, but derisking Sec. A/D fidelity at near-zero cost (no solver calls).
% ===========================================================================================
function test_classify_accepted_only_and_thresholds(testCase)
% not accepted, arbitrarily negative Dq/D_uni -> unconverged, NEVER unstable (the crux rule).
[c, s, ~, ~] = invz_task2_classify(struct('accepted', false, 'D_uni', -50, 'Dq_min', -50));
verifyEqual(testCase, c, 'unconverged');
verifyTrue(testCase, isnan(s) || isfinite(s));   % s is diagnostic-only here; never asserted stable/unstable

% accepted, s > 1e-3 -> stable.
[c, s] = invz_task2_classify(struct('accepted', true, 'D_uni', 0.5, 'Dq_min', 0.2));
verifyEqual(testCase, c, 'stable');  verifyEqual(testCase, s, 0.2);

% accepted, |s| <= 1e-3 (boundary, both signs) -> marginal.
[c] = invz_task2_classify(struct('accepted', true, 'D_uni', 1e-3, 'Dq_min', 0.5));
verifyEqual(testCase, c, 'marginal');
[c] = invz_task2_classify(struct('accepted', true, 'D_uni', -1e-3, 'Dq_min', 0.5));
verifyEqual(testCase, c, 'marginal');

% accepted, s < -1e-3 -> unstable.
[c, s] = invz_task2_classify(struct('accepted', true, 'D_uni', -0.5, 'Dq_min', 0.2));
verifyEqual(testCase, c, 'unstable');  verifyEqual(testCase, s, -0.5);

% accepted but non-finite Dq -> unconverged (prereg Sec. A: "checker NOT accepted, OR
% non-finite" -- the OR is independent of the accepted flag).
[c] = invz_task2_classify(struct('accepted', true, 'D_uni', 0.5, 'Dq', [0.1, NaN, 0.3]));
verifyEqual(testCase, c, 'unconverged');

% Dq vector reduction is NaN-propagating (mirrors invz_ordered_residual's robust_max_abs).
[c, s] = invz_task2_classify(struct('accepted', true, 'D_uni', 0.5, 'Dq', [0.6, 0.4, 0.9]));
verifyEqual(testCase, c, 'stable');  verifyEqual(testCase, s, 0.4);
end

function test_ladder_ok_class_and_numeric_convergence(testCase)
r1.class = {'stable', 'marginal'};  r1.s = [2e-3, 5e-4];
r1.D_uni = [3e-3, 6e-4];  r1.Dq_min = [2e-3, 5e-4];

r2 = r1;   % identical rung -> resolved
[resolved, ~] = invz_task2_ladder_ok([r1, r2]);
verifyTrue(testCase, resolved);

r3 = r1;  r3.class{1} = 'unstable';   % field 1 class disagrees across rungs
[resolved2, detail2] = invz_task2_ladder_ok([r1, r3]);
verifyFalse(testCase, resolved2);
verifyTrue(testCase, detail2.class_disagree(1));
verifyFalse(testCase, detail2.class_disagree(2));

r4 = r1;  r4.s(1) = r1.s(1) + 0.1;    % field 1 s fails |Delta| <= 1e-6+1e-4*max(.) badly
[resolved3, detail3] = invz_task2_ladder_ok([r1, r4]);
verifyFalse(testCase, resolved3);
verifyTrue(testCase, detail3.numeric_disagree(1));
verifyFalse(testCase, detail3.numeric_disagree(2));
end
