function tests = test_invz_ordered_residual_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Build a fixed-h node directly. Do NOT call invz_solve_point_ordered here: that public solver
% does not become scheme-aware until Task 14, and committing a test that knowingly fails until
% then violates the per-task green-suite rule.
function [node, state] = build_strict_fixture()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;  Jxx0 = ion.Jxx0;  tmf = 'legacy_x';  Ecut = 40;
[wn, wts, beta] = invz_matsubara(T, Ecut);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', Jxx0, ...
                                        'transverse_mf', tmf));
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
G0el0   = G0bare0 - G0inel0;
g = real(invz_g(tl, 1i*wn));
mom = invz_coupling_moments(Jnu);
eso = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom);
eopts = struct('static_medium', 'strict_1z_dyson_ref', 'Jmom', mom);
node = struct('tl', tl, 'G0', G0, 'g', g, 'wts', wts, 'wn', wn, 'beta', beta, ...
              'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0el0, 'G0bare0', G0bare0, ...
              'eso', eso, 'eopts', eopts, 'Jnu_flat', Jnu, 'Jmom', mom);
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
assert(info.accepted, 'invz:testFixture', ...
       'strict residual fixture did not reach an accepted state; choose a deterministic algebra fixture');
end

% Block B is revised IN PLACE: same field name, new load-bearing residual. No fifth block.
function test_blockB_is_the_algebraic_k0_check(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyTrue(testCase, isfield(res, 'blockB'));
verifyFalse(testCase, isfield(res, 'strict'));        % NOT a parallel block
verifyEqual(testCase, res.blockB.scheme, 'strict_1z_dyson_ref');
verifyEqual(testCase, res.blockB.status, 'ok');
verifyLessThan(testCase, res.blockB.resid, res.blockB.scale_abs + ...
    res.blockB.scale_rel*max([abs(state.K0s), 6.0e-3]));
verifyTrue(testCase, res.blockB.pass);
end

% The gate must not be scaled by a vanishing correction or by eps(Jbar) (prereg §2).
function test_gate_has_a_problem_native_floor(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyGreaterThan(testCase, res.blockB.scale_abs + res.blockB.scale_rel*max(abs(node.Jnu_flat)), 1e-15);
end

% Mis-wiring is exactly what this residual exists to catch.
function test_mis_wired_k0_fails_blockB(testCase)
[node, state] = build_strict_fixture();
state.K0s = state.K0s * 1.01;                          % 1% off the strict formula
state.K(1) = state.K0s;
res = invz_ordered_residual(node, state);
verifyFalse(testCase, res.blockB.pass);
verifyFalse(testCase, res.accepted);
end

% The discarded resummed closure is opt-in ONLY, and never feeds .finite/.accepted.
function test_resummed_diagnostic_is_opt_in(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyFalse(testCase, isfield(res.blockB, 'resid_resummed'));
dbg = invz_ordered_residual(node, state, struct('debug_resummed', true));
verifyTrue(testCase, isfield(dbg.blockB, 'resid_resummed'));
verifyEqual(testCase, dbg.accepted, res.accepted);     % diagnostic changes no verdict
verifyEqual(testCase, dbg.finite, res.finite);
end

% THE TWO-TIER SEPARATION (spec §1, G4). Stability is computed but never gates acceptance --
% intermediate path nodes ARE the unstable Landau interval by construction.
function test_stability_is_separate_and_never_gates_acceptance(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
verifyTrue(testCase, isfield(res, 'stability'));
for f = {'crit', 'D_uni', 'Dq_min', 'class', 'pass'}
    verifyTrue(testCase, isfield(res.stability, f{1}), f{1});
end
verifyTrue(testCase, any(strcmp(res.stability.class, ...
    {'stable', 'unstable', 'boundary_band', 'undefined'})));
% forcing the stability verdict negative must NOT change res.accepted
node2 = node;  node2.J0eff = 10*node.J0eff;            % drives crit/D_uni negative
res2 = invz_ordered_residual(node2, state);
verifyEqual(testCase, res2.stability.class, 'unstable');
verifyEqual(testCase, res2.accepted, res.accepted, ...
    'stability must not feed acceptance (spec §1 two-tier verdict)');
end

% crit is the dimensionless mass r + J0eff*G0bare, not an inverse susceptibility.
function test_crit_definition(testCase)
[node, state] = build_strict_fixture();
res = invz_ordered_residual(node, state);
[~, ~, so] = invz_emt_static_ordered(node.tl, state.lam(1:2), state.Sigma(1), node.Jnu_flat, ...
    state.K0s, node.beta, node.J0eff, node.G0inel0, node.G0el0, node.eso);
verifyEqual(testCase, res.stability.crit, so.r + node.J0eff*node.G0bare0, 'RelTol', 1e-12);
end

% Wiring errors must escape the residual checker, not become a failed block (spec §5.1).
function test_wiring_error_is_not_absorbed(testCase)
[node, state] = build_strict_fixture();
node.eso.static_medium = 'not_a_scheme';
verifyError(testCase, @() invz_ordered_residual(node, state), 'invz:staticMedium');
end

% Missing Jmom is a wiring error under strict, and harmless under resummed.
function test_jmom_required_only_under_strict(testCase)
[node, state] = build_strict_fixture();
bad = rmfield(node, 'Jmom');
verifyError(testCase, @() invz_ordered_residual(bad, state), 'invz:residualNode');
leg = bad;  leg.eso.static_medium = 'resummed';  leg.eopts.static_medium = 'resummed';
r = invz_ordered_residual(leg, state);
verifyTrue(testCase, isstruct(r));                     % no throw on the legacy path
end

% A strict reference-domain event returns a complete nonaccepted residual WITHOUT evaluating
% A/C/D on an invalid medium.
function test_domain_preflight_short_circuits_checker(testCase)
[node, state] = build_strict_fixture();
state.Sigma(1) = -1;                                  % 1+Sigma0 = 0
res = invz_ordered_residual(node, state);
verifyEqual(testCase, res.blockB.status, 'ref_denom_nonpositive');
verifyFalse(testCase, res.accepted);
verifyFalse(testCase, res.finite);
verifyEqual(testCase, res.stability.class, 'undefined');
verifyTrue(testCase, all(isnan([res.blockA.resid,res.blockC.resid,res.blockD.resid])));
end

% Review-fix (Important finding #3): Block D's medD.dynamic_converged switch
% (invz_ordered_residual.m:273,276) had no discriminating test -- reverting it to whole-PM
% medD.converged left the whole suite green. Corrupt ONLY node.G0(1) (eso/eopts/state all
% untouched -- no scheme split, per the controller's explicit constraint) and Block D's own
% fresh invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts) call gets K(1)=G(1)=NaN
% straight out of the closed-form solve (invz_emt_scalar.m: A(1)/K(1)/G(1) alone go NaN; every
% other frequency is an independent elementwise computation of its own G0 element, so 2:end is
% untouched). medD.converged (invz_emt_scalar.m:128, ALL slots) is therefore false;
% medD.dynamic_converged (invz_emt_scalar.m:99, slots 2:end only) is true -- exactly the
% configuration the switch exists to handle. No other block reads node.G0(1) after Block A's
% own static overwrite (contract Sec. 3: local_F's step (2) clobbers Kf(1) before it is ever
% used), so A/B/C, and blockD.resid itself (which compares state.K(2:end) to medD.K(2:end)),
% come back BIT-IDENTICAL to the unperturbed fixture -- measured, not assumed (task9_fix_probe.m
% probe): only .converged vs .dynamic_converged differ. Reverting :273/:276 to medD.converged
% flips blockD.pass (and res.accepted) to false on this exact fixture -- the discriminating
% evidence review finding #3 asked for.
function test_blockD_gates_on_dynamic_converged_not_whole_pm(testCase)
[node, state] = build_strict_fixture();
res0 = invz_ordered_residual(node, state);
verifyTrue(testCase, res0.accepted, 'fixture must be a genuine accepted state to begin with');

node2 = node;
node2.G0(1) = NaN;                          % corrupt ONLY the discarded PM static slot

% Confirm the fixture actually discriminates BEFORE trusting the checker's own verdict.
medD = invz_emt_scalar(node2.G0, state.Sigma, node2.Jnu_flat, node2.eopts);
verifyFalse(testCase, medD.converged, ...
    'fixture no longer discriminates: medD.converged unexpectedly true');
verifyTrue(testCase, medD.dynamic_converged, ...
    'fixture no longer discriminates: medD.dynamic_converged unexpectedly false');
verifyTrue(testCase, all(isfinite(medD.K(2:end))) && all(isfinite(medD.G(2:end))));

res = invz_ordered_residual(node2, state);
% Every OTHER block is untouched -- isolates the corruption to exactly Block D's convergence
% gate, so a failure below can only be explained by :273/:276.
verifyEqual(testCase, res.blockA.resid, res0.blockA.resid);
verifyEqual(testCase, res.blockB.resid, res0.blockB.resid);
verifyEqual(testCase, res.blockC.resid, res0.blockC.resid);
verifyEqual(testCase, res.blockD.resid, res0.blockD.resid);   % same rD; only .pass can differ
verifyTrue(testCase, res.blockD.pass, ...
    'blockD.pass must gate on medD.dynamic_converged, not whole-PM medD.converged');
verifyTrue(testCase, res.accepted);
end
