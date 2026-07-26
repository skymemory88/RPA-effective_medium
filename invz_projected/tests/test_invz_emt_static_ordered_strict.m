function tests = test_invz_emt_static_ordered_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [tl, args] = fixture()
beta = 1/(0.0862*0.31);
tl = struct('Delta', 0.02, 'M2', 0.8, 'm', 0.6, 'n01', tanh(0.02*beta/2), 'g0', 1);
tl.g0  = 2*tl.n01/tl.Delta;
args = struct('lam', [0.01; 0.02], 'Sigma0', 0.05, ...
              'Jnu', linspace(-2e-3, 6e-3, 24).', 'K0_seed', 0, 'beta', beta, ...
              'J0eff', 6.42444e-3, 'G0inel0', -300, 'G0el0', -20);
end

function [K0, Gstat, out] = call(tl, a, o)
[K0, Gstat, out] = invz_emt_static_ordered(tl, a.lam, a.Sigma0, a.Jnu, a.K0_seed, a.beta, ...
                                           a.J0eff, a.G0inel0, a.G0el0, o);
end

% Absent field => legacy iteration, bit-identical.
function test_absent_scheme_is_bit_identical(testCase)
[tl, a] = fixture();
[K1, G1] = call(tl, a, struct('warn', false));
[K2, G2] = call(tl, a, struct('warn', false, 'static_medium', 'resummed'));
verifyEqual(testCase, K1, K2, 'AbsTol', 0);
verifyEqual(testCase, G1, G2, 'AbsTol', 0);
end

% Strict mode does not iterate, and K0 is exactly the primitive composition.
function test_strict_is_one_shot(testCase)
[tl, a] = fixture();
[K0, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.iters, 0);
verifyEqual(testCase, out.medium_status, 'ok');
mom  = invz_coupling_moments(a.Jnu);
Gref = invz_static_medium_reference(a.G0inel0 + a.G0el0, a.Sigma0, 'strict_1z_dyson_ref');
Kexp = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, K0, Kexp, 'AbsTol', 0);
end

% The strict residual is the algebraic K0 check, identically zero for a correct one-shot call,
% and out.converged reflects DOMAIN validity, not iteration success.
function test_strict_residual_is_the_algebraic_check(testCase)
[tl, a] = fixture();
[~, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.resid, 0, 'AbsTol', 0);
verifyTrue(testCase, out.converged);
end

% The K0_seed is IGNORED under strict: a one-shot medium has no warm start to inherit, so a
% contaminated seed can no longer propagate between nodes.
function test_seed_is_ignored_under_strict(testCase)
[tl, a] = fixture();
o = struct('warn', false, 'static_medium', 'strict_1z_dyson_ref');
[Ka] = call(tl, a, o);
a.K0_seed = 0.05;
[Kb] = call(tl, a, o);
verifyEqual(testCase, Ka, Kb, 'AbsTol', 0);
end

% Domain event: a status, no throw, no warning flood.
function test_out_of_domain_reference_is_a_status(testCase)
[tl, a] = fixture();
a.Sigma0 = -1;                                   % 1 + Sigma0 = 0
[K0, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.medium_status, 'ref_denom_nonpositive');
verifyTrue(testCase, isnan(K0));
verifyFalse(testCase, out.converged);
end

% Strict mode emits no invz:emtStatic warning even with warn = true: nothing iterates.
function test_no_closure_warning_under_strict(testCase)
[tl, a] = fixture();
w = warning('off', 'all');  restore = onCleanup(@() warning(w));
lastwarn('');
call(tl, a, struct('warn', true, 'static_medium', 'strict_1z_dyson_ref'));
[~, id] = lastwarn();
verifyNotEqual(testCase, id, 'invz:emtStatic');
end

% Both omitted-order ratios are surfaced for the caller's promotion gate.
function test_omitted_ratios_exposed(testCase)
[tl, a] = fixture();
[~, ~, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyTrue(testCase, isfinite(out.omit_mu3) && isfinite(out.omit_cubic));
verifyEqual(testCase, out.omit_max, max(out.omit_mu3, out.omit_cubic), 'AbsTol', 0);
end

% Dq / D_uni are still built from the physical Gstat and still reported in full (spec §0).
function test_collective_observables_still_reported(testCase)
[tl, a] = fixture();
[K0, Gstat, out] = call(tl, a, struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, out.D_uni, 1 + (a.J0eff - K0)*Gstat, 'RelTol', 1e-12);
Dq = 1 + (a.Jnu - K0).*Gstat;
verifyEqual(testCase, out.Dq_min, min(Dq), 'AbsTol', 0);
verifyEqual(testCase, out.Dq_max, max(Dq), 'AbsTol', 0);
verifyEqual(testCase, out.Dq_neg_count, nnz(Dq <= 0));
verifyTrue(testCase, isfinite(out.r));                       % Task 6 reassociation in effect
end
