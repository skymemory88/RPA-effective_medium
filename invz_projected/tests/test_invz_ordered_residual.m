function tests = test_invz_ordered_residual
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% ===========================================================================================
% Shared fixture: the Task-0-pinned synthetic control (test_invz_ordered_trace.m:89-93),
% T=0.31, B=[2.85 0 0], Jnu=linspace(-2e-3,6e-3,24).', J0eff=6.42e-3, Jxx0=ion.Jxx0,
% ordered_mode='jensen'. invz_solve_point_ordered (the CURRENT, unmodified production
% solver) does not itself export the node-input bundle invz_ordered_residual needs (tl and
% si are exported as pt.tl/pt.si; G0/g/wts/wn/G0inel0/G0el0/G0bare0/eso/eopts are not) --
% so this fixture reconstructs them from the SAME public calls the solver's own preamble
% uses (invz_solve_point_ordered.m:177-205), documented step by step below, rather than
% reaching into the solver's internals.
% ===========================================================================================
function [node, state, pt] = build_fixture()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;                      % = Jcc0, the Task-0-pinned synthetic-control value
Jxx0 = ion.Jxx0;
tmf  = 'legacy_x';                    % production default (opts.transverse_mf), not overridden
Ecut = 40;                            % production default (opts.Ecut), not overridden

o = struct('J0eff', J0eff, 'Jxx0', Jxx0, 'hyp', true, 'ordered_mode', 'jensen');
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, o);

% --- reconstruct node inputs (invz_solve_point_ordered.m:110,177-204 verbatim) -----------
[wn, wts, beta] = invz_matsubara(T, Ecut);
c0  = invz_chi0z(pt.si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(pt.si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X = real(c0(:, :, 1));
switch tmf                                    % invz_solve_point_ordered.m:192-202, verbatim
    case 'legacy_x'
        fb = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
    case 'none'
        fb = 0;
    otherwise
        error('invz:fixtureTmf', 'fixture only reconstructs none/legacy_x; got ''%s''.', tmf);
end
G0bare0 = -(X(3,3) + fb);
G0el0   = G0bare0 - G0inel0;
g = real(invz_g(pt.tl, 1i*wn));

node = struct('tl', pt.tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
              'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
              'eso', struct('warn', false), 'eopts', struct(), 'Jnu_flat', Jnu);
state = struct('Sigma', pt.Sigma, 'K', pt.K, 'lam', pt.lambda, 'K0s', pt.K(1));
end

% ===========================================================================================
% (i) Accepted state passes all blocks.
% ===========================================================================================
function test_accepted_state_passes_all_blocks(testCase)
[node, state, pt] = build_fixture();
verifyTrue(testCase, pt.is_ordered && pt.converged, 'fixture solve did not converge');

res = invz_ordered_residual(node, state);

verifyTrue(testCase, res.blockA.pass, sprintf('block A: resid=%.3e gate=%.3e', ...
    res.blockA.resid, res.blockA.scale_abs));
verifyTrue(testCase, res.blockB.pass, sprintf('block B: resid=%.3e gate=%.3e+%.3e*|Gstat|', ...
    res.blockB.resid, res.blockB.scale_abs, res.blockB.scale_rel));
verifyTrue(testCase, res.blockC.pass, sprintf('block C: resid=%.3e gate=%.3e', ...
    res.blockC.resid, res.blockC.scale_abs));
verifyTrue(testCase, res.blockD.pass, sprintf('block D: resid=%.3e gate=%.3e+%.3e*ref', ...
    res.blockD.resid, res.blockD.scale_abs, res.blockD.scale_rel));
verifyTrue(testCase, res.finite);
verifyTrue(testCase, res.accepted);

% Block C reproduces production's existing final_resid quantity exactly (contract Sec. 4:
% "this is the current final_resid, keep it but as one named block").
verifyEqual(testCase, res.blockC.resid, pt.final_resid, 'AbsTol', 1e-14);

% Block D's excluded-K(1) design (contract Sec. 4): the elastic-hybrid exported K(1) must
% differ from the ordinary-Dyson K(1) invz_emt_scalar alone produces at the SAME Sigma --
% confirms block D is comparing the physically right thing (matches the existing
% test_invz_ordered_jensen.m pin that pt.G(1) differs from the ordinary-Dyson value).
medD = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts);
verifyGreaterThan(testCase, abs(medD.K(1) - state.K(1))/abs(state.K(1)), 1e-3);

% diagnostics exposed, never folded away (contract Sec. 4/7).
verifyTrue(testCase, isfinite(res.D_uni));
verifyTrue(testCase, isfinite(res.Dq_min) && isfinite(res.Dq_max));
verifyTrue(testCase, res.blockB.converged);
end

% ===========================================================================================
% (ii-a) Perturbed state fails the RIGHT block: a Sigma-space corruption hits A (the outer
% map is no longer at its own fixed point) and, as a natural consequence, C (the derived
% chain, re-evaluated against the now-wrong Sigma, disagrees too) and D (K is Sigma-dependent
% through invz_emt_scalar's own closed-form solve, so a Sigma corruption cascades into a
% K-space disagreement as well). Block B is EXPECTED to keep passing: its own internal
% Picard loop re-converges around whatever Sigma(1)/lam it is handed, so it measures internal
% self-consistency, not global correctness (contract Sec. 4/9) -- this is asserted
% explicitly, not merely tolerated, so a future change that made B spuriously fail here would
% also be caught.
% ===========================================================================================
function test_sigma_perturbation_fails_A_and_C(testCase)
[node, state] = build_fixture();
res0 = invz_ordered_residual(node, state);
verifyTrue(testCase, res0.accepted, 'fixture must be a genuine accepted state to begin with');

pstate = state;
pstate.Sigma = state.Sigma + 0.01;         % uniform shift, several orders above every gate
res = invz_ordered_residual(node, pstate);

verifyFalse(testCase, res.blockA.pass);
verifyFalse(testCase, res.blockC.pass);
verifyTrue(testCase, res.blockB.pass);
verifyFalse(testCase, res.accepted);
end

% ===========================================================================================
% (ii-b) Perturbed state fails the RIGHT block: a K(2:end)-only corruption (Sigma/lam/K0s
% untouched) hits D (the direct target -- "the dynamic-medium perturbation") and, as a
% natural consequence, C (its own lam/sg recomputation reads the now-perturbed K). Blocks A
% and B are proven -- not merely assumed -- INDEPENDENT of state.K: their residuals come back
% BIT-IDENTICAL to the unperturbed run, because neither block's recomputation ever reads
% state.K (contract Sec. 3/4). This is the cleanest single-block-discriminating signature
% the four-block design offers, and is asserted at full strength (equality, not just "still
% passes").
% ===========================================================================================
function test_dynamic_medium_perturbation_fails_D_not_A_or_B(testCase)
[node, state] = build_fixture();
res0 = invz_ordered_residual(node, state);
verifyTrue(testCase, res0.accepted, 'fixture must be a genuine accepted state to begin with');

pstate = state;
pstate.K(2:end) = state.K(2:end) + 1e-4;   % dynamic entries only, several orders above gate
res = invz_ordered_residual(node, pstate);

verifyFalse(testCase, res.blockD.pass);
verifyFalse(testCase, res.blockC.pass);
verifyFalse(testCase, res.accepted);
verifyTrue(testCase, res.blockA.pass);
verifyTrue(testCase, res.blockB.pass);
verifyEqual(testCase, res.blockA.resid, res0.blockA.resid);   % bit-identical: proves A never
verifyEqual(testCase, res.blockB.resid, res0.blockB.resid);   % reads state.K; same for B.
end

% ===========================================================================================
% (iii) Non-mutation: node/state must be untouched by the checker call, on BOTH an accepted
% and a failing evaluation.
% ===========================================================================================
function test_non_mutation(testCase)
[node, state] = build_fixture();
node_before = node;  state_before = state;
res = invz_ordered_residual(node, state); %#ok<NASGU>
verifyTrue(testCase, isequaln(node, node_before));
verifyTrue(testCase, isequaln(state, state_before));

pstate = state;  pstate.Sigma = state.Sigma + 0.01;
node_before2 = node;  pstate_before = pstate;
res2 = invz_ordered_residual(node, pstate); %#ok<NASGU>
verifyTrue(testCase, isequaln(node, node_before2));
verifyTrue(testCase, isequaln(pstate, pstate_before));
end

% ===========================================================================================
% (iv) stall / finite semantics.
% ===========================================================================================
function test_stall_is_diagnostic_only_never_gates_accepted(testCase)
[node, state] = build_fixture();

% Accepted state + a small caller-supplied dS: NOT stalled (the complete residual already
% agrees with the old step-size signal).
res_ok = invz_ordered_residual(node, state, struct('dS', 1e-10));
verifyFalse(testCase, res_ok.stall);
verifyTrue(testCase, res_ok.accepted);

% opts.dS omitted: undecidable, reported as NaN, never coerced to true/false.
res_nodS = invz_ordered_residual(node, state);
verifyTrue(testCase, isnan(res_nodS.stall));

% Perturbed (non-accepted) state with a small caller-supplied dS: THIS is the failure mode
% the contract exists to surface -- the OLD gate (small dS) would have accepted the node,
% but the COMPLETE residual says it should not have. stall is a diagnostic; accepted must
% still be false regardless of what stall reports.
pstate = state;  pstate.Sigma = state.Sigma + 0.01;
res_bad = invz_ordered_residual(node, pstate, struct('dS', 1e-10));
verifyTrue(testCase, res_bad.stall);
verifyFalse(testCase, res_bad.accepted);
end

function test_finite_flag_catches_partial_nan_corruption(testCase)
% A single non-finite entry, deep inside an otherwise-finite Sigma, must not be silently
% ignored: MATLAB's max() ignores individual NaNs by default (unlike sum/mean), so this
% additionally exercises the checker's explicit NaN-propagation safeguard, not just the
% top-level isfinite(state.Sigma) check.
[node, state] = build_fixture();
pstate = state;  pstate.Sigma(3) = NaN;
res = invz_ordered_residual(node, pstate);
verifyFalse(testCase, res.finite);
verifyFalse(testCase, res.accepted);
verifyTrue(testCase, isnan(res.blockA.resid));
verifyFalse(testCase, res.blockA.pass);
verifyTrue(testCase, isnan(res.blockC.resid));
verifyFalse(testCase, res.blockC.pass);
end

% ===========================================================================================
% Hard requirement: only 'invz:*' identifiers are ever the checker's own doing; malformed
% inputs raise a clearly-identified invz: error rather than an opaque MATLAB one.
% ===========================================================================================
function test_malformed_inputs_raise_invz_errors(testCase)
[node, state] = build_fixture();

verifyError(testCase, @() invz_ordered_residual(rmfield(node, 'tl'), state), ...
    'invz:residualNode');
verifyError(testCase, @() invz_ordered_residual(node, rmfield(state, 'K0s')), ...
    'invz:residualState');

bad_state = state;  bad_state.lam = [0; 0];
verifyError(testCase, @() invz_ordered_residual(node, bad_state), 'invz:residualState');

bad_state2 = state;  bad_state2.K = state.K(1:end-1);
verifyError(testCase, @() invz_ordered_residual(node, bad_state2), 'invz:residualState');
end
