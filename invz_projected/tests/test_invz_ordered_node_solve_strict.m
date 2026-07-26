function tests = test_invz_ordered_node_solve_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Minimal real node built from public calls at a field where the bare set orders.
function node = build_node(scheme)
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;  Jxx0 = ion.Jxx0;  tmf = 'legacy_x';
[wn, wts, beta] = invz_matsubara(T, 40);
si  = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', Jxx0, ...
                                         'transverse_mf', tmf));
tl  = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', Jxx0, 'transverse_mf', tmf));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (Jxx0 / (1 - Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
node = struct('tl', tl, 'G0', G0, 'g', real(invz_g(tl, 1i*wn)), 'wts', wts, 'wn', wn, ...
    'beta', beta, 'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0bare0 - G0inel0, ...
    'G0bare0', G0bare0, 'eso', struct('warn', false, 'static_medium', scheme), ...
    'eopts', struct('static_medium', scheme), 'Jnu_flat', Jnu, ...
    'Jmom', invz_coupling_moments(Jnu));
end

function test_strict_node_solves_and_reports_medium(testCase)
node = build_node('strict_1z_dyson_ref');
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
verifyEqual(testCase, info.medium_status, 'ok');
verifyEqual(testCase, info.medium.scheme, 'strict_1z_dyson_ref');
verifyTrue(testCase, isfinite(state.K0s));
verifyTrue(testCase, isfield(info.res, 'stability'));
end

% Jmom must actually be threaded into BOTH leaves, not silently re-derived.
function test_jmom_is_threaded_to_both_leaves(testCase)
node = build_node('strict_1z_dyson_ref');
bad = node;
bad.Jmom.mu2 = node.Jmom.mu2 * 2;          % a deliberately wrong moment
[sa] = invz_ordered_node_solve(node, [], struct('trace', false));
[sb] = invz_ordered_node_solve(bad,  [], struct('trace', false));
verifyNotEqual(testCase, sa.K0s, sb.K0s, ...
    'node.Jmom must reach the static leaf; identical K0s means it was re-derived');
end

% Missing Jmom under strict is a wiring error; under resummed it is harmless.
function test_missing_jmom_under_strict_throws(testCase)
node = rmfield(build_node('strict_1z_dyson_ref'), 'Jmom');
verifyError(testCase, @() invz_ordered_node_solve(node, [], struct()), 'invz:nodeSolveNode');
leg = rmfield(build_node('resummed'), 'Jmom');
[~, info] = invz_ordered_node_solve(leg, [], struct());
verifyEqual(testCase, info.medium_status, 'not_applicable');
end

% A domain event stops the attempt BEFORE lambdas/Sigma consume an invalid reference, and is
% reported as its own term_reason -- never as generic max_iter.
function test_domain_event_stops_before_lambdas(testCase)
node = build_node('strict_1z_dyson_ref');
node.eso.ref_margin = 1e9;                 % force every reference denominator out of domain
[state, info] = invz_ordered_node_solve(node, [], struct('trace', false));
verifyEqual(testCase, info.term_reason, 'medium_out_of_domain');
verifyEqual(testCase, info.medium_status, 'ref_denom_small');
verifyFalse(testCase, info.accepted);
verifyEqual(testCase, state.lam, [0; 0; 0], 'AbsTol', 0);   % lambdas never ran
end

% Wiring errors escape; they are not scored as a failed node.
function test_wiring_error_escapes(testCase)
node = build_node('strict_1z_dyson_ref');
node.eso.static_medium = 'not_a_scheme';
verifyError(testCase, @() invz_ordered_node_solve(node, [], struct()), 'invz:staticMedium');
end

% Legacy path unchanged.
function test_resummed_path_untouched(testCase)
node = build_node('resummed');
[~, info] = invz_ordered_node_solve(node, [], struct('trace', false));
verifyTrue(testCase, any(strcmp(info.term_reason, ...
    {'accepted', 'loop_converged_not_accepted', 'max_iter'})));
verifyEqual(testCase, info.medium_status, 'not_applicable');
end
