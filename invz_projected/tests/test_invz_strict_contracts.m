function tests = test_invz_strict_contracts
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% G11: a REAL-coupling ordered anchor, with its full provenance asserted. The original masking
% defect survived a whole stage because no test fed the jensen leg real bz_couplings densities.
function test_G11_real_coupling_ordered_anchor(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1'); end
ion = invz_ion();
prov = struct('grid', [16 16 16], 'dpRng', 30, 'dipole', 'bruteforce', 'cache', false);
[Jnu, info] = invz_bz_couplings(ion, prov);
mom = invz_coupling_moments(Jnu(:));
% provenance first: the hard-coded moments are valid ONLY for this tuple (spec §B)
verifyEqual(testCase, info.dipole.backend, 'bruteforce');
verifyFalse(testCase, isfield(info, 'grid'), ...
    'no grid-policy field means the bit-identical legacy route, whose info.grid is absent');
verifyEqual(testCase, invz_exact_numeric_digest(Jnu(:)), 'ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17');
verifyEqual(testCase, mom.n, 16384);
verifyEqual(testCase, mom.Jbar, 1.20766e-4, 'RelTol', 1e-4);
verifyEqual(testCase, mom.mu2,  5.48264e-6, 'RelTol', 1e-4);
% Build-blocking real anchor: the diagnosed 1 T production point must now have a solved root.
o = struct('J0eff', info.Jcc0, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', 'Jmom', mom, 'trace', true);
[hstar, p, trc] = invz_hmf_ordered(ion, 0.1, [1 0 0], Jnu(:), o);
verifyEqual(testCase, p.status, 'ok');
verifyTrue(testCase, isfinite(hstar) && hstar > 0);
verifyTrue(testCase, all(strcmp(p.medium_status, 'ok')));
verifyTrue(testCase, trc.enabled && ~isempty(trc.nodes));
verifyTrue(testCase, all(strcmp({trc.nodes.medium_status}, 'ok')));
verifyTrue(testCase, isfinite(p.crit_star));
% Dq diagnostics still need the FULL multiset, not just the two moments
verifyEqual(testCase, numel(Jnu(:)), 16384);
end

% G13: behavioural sentinel -- the PM slot must not leak into ordered lambdas. Not a
% source-text-order assertion (those are brittle; a prior test-regex bug on this branch came
% from exactly that style).
function test_G13_pm_slot_does_not_leak_into_ordered_lambdas(testCase)
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];  hz = 0.02;
Jnu = linspace(-2e-3, 6.0e-3, 24).';  J0eff = 6.42e-3;
[wn, wts, beta] = invz_matsubara(T, 40);
si = invz_single_ion(ion, T, Bx, struct('hyp', true, 'hz_fixed', hz, 'Jxx0', ion.Jxx0, ...
                                        'transverse_mf', 'legacy_x'));
tl = invz_twolevel_ordered(ion, T, Bx, hz, struct('Jxx0', ion.Jxx0, ...
                                                  'transverse_mf', 'legacy_x'));
c0  = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
G0  = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si, T, 1i*wn(1), struct('elastic', false));
G0inel0 = -real(c0i(3,3,1));
X  = real(c0(:, :, 1));
fb = X(3,1) * (ion.Jxx0/(1 - ion.Jxx0*X(1,1))) * X(1,3);
G0bare0 = -(X(3,3) + fb);
node = struct('tl', tl, 'G0', G0, 'g', real(invz_g(tl, 1i*wn)), 'wts', wts, 'wn', wn, ...
    'beta', beta, 'J0eff', J0eff, 'G0inel0', G0inel0, 'G0el0', G0bare0 - G0inel0, ...
    'G0bare0', G0bare0, 'eso', struct('warn', false, 'static_medium', 'strict_1z_dyson_ref'), ...
    'eopts', struct('static_medium', 'strict_1z_dyson_ref'), 'Jnu_flat', Jnu, ...
    'Jmom', invz_coupling_moments(Jnu));
[state, ~] = invz_ordered_node_solve(node, [], struct('trace', false));
% the exported K(1) must be the ORDERED static value, and the lambdas must be the ones derived
% from it -- recompute both and compare against invz_emt_scalar's own slot 1
med = invz_emt_scalar(node.G0, state.Sigma, node.Jnu_flat, node.eopts);
verifyNotEqual(testCase, state.K(1), med.K(1), ...
    'ordered K(1) must be the elastic-hybrid static value, not the PM slot');
lam_from_exported = invz_lambdas(state.K, node.g, node.wts, node.beta, [1 2 3]);
verifyEqual(testCase, state.lam, lam_from_exported, 'RelTol', 1e-10, ...
    'lambdas must derive from the exported K WITH the ordered K(1) substituted');
Kleak = state.K;  Kleak(1) = med.K(1);
lam_leaked = invz_lambdas(Kleak, node.g, node.wts, node.beta, [1 2 3]);
verifyGreaterThan(testCase, max(abs(lam_leaked - state.lam)), 1e-12, ...
    'the sentinel must be able to detect a leak at all');
end

% G15: fatal ids escape every layer; recoverable/domain outcomes keep their exact category.
function test_G15_fatal_ids_escape_every_layer(testCase)
verifyFalse(testCase, invz_is_recoverable_solver_error('invz:staticMedium'));
verifyTrue(testCase,  invz_is_recoverable_solver_error('invz:degenerateDoublet'));
[~, completed, rid] = invz_try_solver_call(@() local_degenerate_thrower());
verifyFalse(testCase, completed);
verifyEqual(testCase, rid, 'invz:degenerateDoublet');
% node solver
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'not_a_scheme');
verifyError(testCase, @() invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o), ...
            'invz:staticMedium');
% point solver
verifyError(testCase, @() invz_solve_point(ion, 0.31, [2.85 0 0], Jnu, o), ...
            'invz:staticMedium');
% auto dispatcher: inject a fatal error INSIDE its ordered/PM try blocks, not at its own entry.
obad = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
              'static_medium', 'strict_1z_dyson_ref');
Jbad = complex(Jnu, 1e-12);
verifyError(testCase, @() invz_solve_auto(ion, 0.31, [2.85 0 0], Jbad, obad), ...
            'invz:couplingMoments');
% spectra outer boundary: the same injected fatal must cross one_field/parfor unchanged.
sp = struct('Jnu', Jbad, 'info', struct('Jcc0', 6.42e-3), ...
    'static_medium', 'strict_1z_dyson_ref', 'parallel', false, 'verbose', false);
verifyError(testCase, @() invz_spectra_map(ion, 0.31, 2.85, (0.02:0.04:0.42).', sp), ...
            'invz:couplingMoments');
end

% G15b: a domain outcome must never be reported as generic node_failed.
function test_G15b_domain_outcome_keeps_its_category(testCase)
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref', 'ref_margin', 1e9);
[~, p] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o);
verifyEqual(testCase, p.status, 'medium_out_of_domain');
verifyNotEqual(testCase, p.status, 'node_failed');
end

function v = local_degenerate_thrower() %#ok<STOUT> -- v is intentionally never assigned: the
% function always throws before returning, so there is no normal-return path that could leave a
% caller holding an unset value. The declared output exists only so invz_try_solver_call's
% nargout probe sees 1 rather than -1. A bare `@() error(...)` surfaces as MATLAB:maxlhs and masks
% its own identifier -- see that function's docstring. Measured: bare form throws MATLAB:maxlhs;
% this form returns rid = the thrown id.
error('invz:degenerateDoublet', 'synthetic');
end
