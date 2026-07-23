function tests = test_invz_ordered_node_solve
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% ===========================================================================================
% Shared fixture: the Task-0-pinned synthetic control (test_invz_ordered_trace.m:89-93; the
% SAME fixture test_invz_ordered_residual.m's build_fixture() uses), T=0.31, B=[2.85 0 0],
% Jnu=linspace(-2e-3,6e-3,24).', J0eff=6.42e-3, Jxx0=ion.Jxx0, tmf='legacy_x', Ecut=40.
% Task 1b-i's brief permits factoring a shared fixture helper but forbids modifying
% test_invz_ordered_residual.m, and the task's own gate caps the diff at exactly two new
% files -- so this is a duplicate LOCAL copy (not a shared import), matching that file's
% build_fixture() statement-for-statement (invz_solve_point_ordered.m:177-205's own
% preamble), not a new construction.
% ===========================================================================================
function [node, seed_prod, pt] = build_fixture()
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
seed_prod = struct('Sigma', pt.Sigma, 'K0s', pt.K(1));   % the production converged state,
                                                          % reusable as a warm seed by (b)/(f)
end

% ===========================================================================================
% (a) Cold solve accepts and matches production.
% ===========================================================================================
function test_cold_solve_accepts_and_matches_production(testCase)
[node, ~, pt] = build_fixture();
verifyTrue(testCase, pt.is_ordered && pt.converged, 'fixture solve did not converge');

[st, inf] = invz_ordered_node_solve(node, []);

verifyTrue(testCase, inf.accepted);
verifyTrue(testCase, inf.res.accepted);
verifyEqual(testCase, inf.seed_kind, 'cold');
verifyEqual(testCase, inf.term_reason, 'accepted');

% Named-quantity match, explicit tolerance (brief: "NOT bitwise" -- assert closeness, not
% exact equality). Empirically (probe run against this exact fixture) the cold path
% reproduces production BIT-FOR-BIT (max|st.X - pt.X| = 0 for Sigma/K/lam/K0s): the cold
% branch here executes the IDENTICAL statement sequence invz_solve_point_ordered.m's jensen
% loop does, from the identical (Sigma=0, K0s=0, lam=[0;0;0]) start, with the identical
% mix_outer/tol_outer/max_outer defaults, so this is expected, not merely hoped for. The
% AbsTol/RelTol below is a declared, generous tolerance (matching the outer tolerance,
% tol_outer=1e-8) per the brief's instruction -- not a sign that a real discrepancy exists.
verifyEqual(testCase, st.Sigma, pt.Sigma,  'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, st.K,     pt.K,      'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, st.lam,   pt.lambda, 'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, st.K0s,   pt.K(1),   'AbsTol', 1e-8, 'RelTol', 1e-6);
verifyEqual(testCase, inf.outer_iters, pt.outer_iters);

% diagnostics present and sane (contract-aligned: res carries blockA-D, D_uni, Dq_min/max).
verifyTrue(testCase, isfinite(inf.res.D_uni));
verifyTrue(testCase, isfinite(inf.res.Dq_min) && isfinite(inf.res.Dq_max));
verifyTrue(testCase, inf.loop_converged);     % the in-loop gate ALSO fired here (diagnostic)
end

% ===========================================================================================
% (b) Warm solve from the production state is a near-fixed-point.
% ===========================================================================================
function test_warm_solve_from_production_state_is_near_fixed_point(testCase)
[node, seed_prod, pt] = build_fixture();

[st, inf] = invz_ordered_node_solve(node, seed_prod);

verifyTrue(testCase, inf.accepted);
verifyEqual(testCase, inf.seed_kind, 'warm');
% "very few" iterations: strictly fewer than the cold solve's own outer_iters (a warm start
% from an ALREADY-converged production state should need less work than starting from
% Sigma=0), AND well below max_outer=200. Empirically outer_iters=10 here (vs. the cold
% path's 14) -- the residual gap from the one-iteration lag intrinsic to lam always
% resetting to [0;0;0] (see invz_ordered_node_solve.m's header) costs a handful of extra
% iterations to re-settle, not a slow reconvergence.
verifyLessThan(testCase, inf.outer_iters, pt.outer_iters);
verifyLessThan(testCase, inf.outer_iters, 50);

res_check = invz_ordered_residual(node, st);
verifyTrue(testCase, res_check.accepted);
end

% ===========================================================================================
% (c) Cold retry triggers on a bad warm seed.
%
% The brief's OWN suggested example (Sigma = ones(...)*1e3, a uniform large shift) was tried
% FIRST (empirically, via a scratch probe) and does NOT fail the warm attempt on this
% fixture: a uniform shift keeps D = 1+Sigma huge relative to Jnu_flat.*G0 at every
% frequency, so invz_emt_scalar's A = mean(Jnu_flat./(D+Jnu_flat.*G0)) stays small and
% bounded rather than hitting any pole, and the damped mix (contraction ratio 1-mix_outer =
% 0.3 per iteration) relaxes a Sigma~1e3 perturbation back below tol_outer=1e-8 in ~23
% outer iterations -- comfortably inside max_outer=200, i.e. it is RESCUED WITHOUT any
% retry, not a failing seed here. (K0s alone perturbed to +-1e10 is similarly harmless: the
% static closure's own internal Picard loop re-converges K0 from that seed within its first
% call, independent of how far off the seed was.)
%
% The corruption below is chosen because it demonstrably fails: Sigma(1) = -1 (K0s = 0) puts
% the STATIC sector's own elastic denominator (invz_gstat_ordered.m: "Gstat = G0inel0/(1 +
% Sigma0 + K0*G0inel0) + xi*G0el0") at 1+Sigma0+K0*G0inel0 = 0 exactly -- a genuine pole the
% outer damped-Picard map cannot recover from within max_outer=200 (confirmed below by (i)).
% ===========================================================================================
function test_cold_retry_rescues_corrupt_warm_seed(testCase)
[node, ~, pt] = build_fixture();
bad_seed = struct('Sigma', [-1; zeros(numel(node.wn) - 1, 1)], 'K0s', 0);

% (i) non-vacuousness: the warm attempt ALONE (retry disabled) really does fail to close --
% this is not a test of a no-op perturbation.
[st_warm, inf_warm] = invz_ordered_node_solve(node, bad_seed, struct('cold_retry', false));
verifyFalse(testCase, inf_warm.accepted, 'the corrupted warm seed must NOT be accepted alone');
verifyEqual(testCase, inf_warm.seed_kind, 'warm');
verifyEqual(testCase, inf_warm.outer_iters, 200);
verifyEqual(testCase, inf_warm.term_reason, 'max_iter');
verifyFalse(testCase, all(isfinite(st_warm.Sigma)), ...
    'expected a genuine divergence (non-finite Sigma), not merely a missed tolerance');

% (ii) the default (cold_retry = true) path rescues it: seed_kind records the rescue, the
% accepted state passes the full checker, and (since the retry re-initialises from the same
% cold start as test (a)) reproduces that same cold fixed point.
[st, inf] = invz_ordered_node_solve(node, bad_seed);
verifyTrue(testCase, inf.accepted);
verifyEqual(testCase, inf.seed_kind, 'cold_after_warm_fail');
res_check = invz_ordered_residual(node, st);
verifyTrue(testCase, res_check.accepted);
verifyEqual(testCase, st.Sigma, pt.Sigma, 'AbsTol', 1e-8, 'RelTol', 1e-6);
end

% ===========================================================================================
% (d) A genuinely non-converging node returns accepted=false WITHOUT throwing.
% ===========================================================================================
function test_non_converging_node_returns_not_accepted_without_throwing(testCase)
[node, ~, ~] = build_fixture();
sopts = struct('max_outer', 1, 'cold_retry', false);   % cannot possibly close in one pass

[st, inf] = invz_ordered_node_solve(node, [], sopts);  % must return normally, not throw

verifyFalse(testCase, inf.accepted);
verifyEqual(testCase, inf.outer_iters, 1);
verifyEqual(testCase, inf.term_reason, 'max_iter');
verifyTrue(testCase, isstruct(inf.res));
verifyTrue(testCase, inf.res.finite, ...
    'one outer iteration from a cold start is still numerically well-behaved, just unclosed');
verifyTrue(testCase, all(isfinite(st.Sigma)) && all(isfinite(st.K)) && ...
    all(isfinite(st.lam)) && isfinite(st.K0s));
end

% ===========================================================================================
% (e) trace records populate only when requested.
% ===========================================================================================
function test_trace_records_populate_only_when_requested(testCase)
[node, ~, ~] = build_fixture();

[~, inf_on] = invz_ordered_node_solve(node, [], struct('trace', true));
verifyGreaterThan(testCase, numel(inf_on.iters), 0);
verifyEqual(testCase, numel(inf_on.iters), inf_on.outer_iters);
rec = inf_on.iters(1);
req_fields = {'dS', 'resid_static', 'K0', 'D_uni', 'Dq_min', 'Dq_max', 'Dq_neg_count', ...
              'idx_pos_flat', 'Dq_pos_val', 'idx_neg_flat', 'Dq_neg_val'};
for k = 1:numel(req_fields)
    verifyTrue(testCase, isfield(rec, req_fields{k}), ...
        sprintf('inf.iters record is missing field ''%s''', req_fields{k}));
end

[~, inf_off] = invz_ordered_node_solve(node, []);      % sopts.trace defaults to false
verifyTrue(testCase, isempty(inf_off.iters));
end

% ===========================================================================================
% (f) Purity: node/seed are never mutated, across the cold, warm, and cold-retry paths.
% ===========================================================================================
function test_purity_no_input_mutation(testCase)
[node, seed_prod, ~] = build_fixture();

node_before1 = node;
[~, ~] = invz_ordered_node_solve(node, []);
verifyTrue(testCase, isequaln(node, node_before1));

node_before2 = node;  seed_before2 = seed_prod;
[~, ~] = invz_ordered_node_solve(node, seed_prod);
verifyTrue(testCase, isequaln(node, node_before2));
verifyTrue(testCase, isequaln(seed_prod, seed_before2));

bad_seed = struct('Sigma', [-1; zeros(numel(node.wn) - 1, 1)], 'K0s', 0);
node_before3 = node;  bad_seed_before3 = bad_seed;
[~, ~] = invz_ordered_node_solve(node, bad_seed);      % exercises BOTH run_attempt calls
verifyTrue(testCase, isequaln(node, node_before3));
verifyTrue(testCase, isequaln(bad_seed, bad_seed_before3));
end
