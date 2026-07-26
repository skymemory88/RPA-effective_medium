function tests = test_invz_hmf_ordered_diagnostics
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, T, Bx, Jnu, o] = fixture()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true, ...
           'static_medium', 'strict_1z_dyson_ref');
end

% Every new per-node array is present and the same length as hgrid, on ANY exit path.
function test_per_node_arrays_are_present_and_aligned(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
n = numel(prof.hgrid);
% Without this, a no_bare_order/n==0 regression would make every numel(...)==n check
% below pass vacuously (0==0) instead of catching the break.
verifyGreaterThan(testCase, n, 0, 'fixture must produce a non-empty grid');
for f = {'crit', 'r_minus_1', 'Delta', 'Dq_min', 'ref_denom', 'ref_margin', ...
         'gstat_local_denom', ...
         'omit_mu3', 'omit_cubic', 'omit_max'}
    verifyEqual(testCase, numel(prof.(f{1})), n, f{1});
end
verifyEqual(testCase, numel(prof.medium_status), n);
verifyEqual(testCase, numel(prof.node_term_reason), n);
verifyTrue(testCase, iscell(prof.medium_status));
end

% Step 7 schema coverage, DOMAIN kind: opts.ref_margin=1e9 forces the strict reference
% denominator out of domain on every node (measured: status='node_failed', hstar=NaN,
% medium_status={ref_denom_small}, node_term_reason={medium_out_of_domain} on all 33
% nodes, ref_denom(1)=1, ref_margin(1)=-999999999). Same structural guarantee as
% test_per_node_arrays_are_present_and_aligned, plus the domain-specific field provenance
% (ref_margin is the ACTUAL distance to the configured floor, never the floor itself).
function test_schema_domain_record(testCase)
[ion, T, Bx, Jnu, o] = fixture();
o.ref_margin = 1e9;
[hstar, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.status, 'node_failed');
verifyTrue(testCase, isnan(hstar));
n = numel(prof.hgrid);
verifyGreaterThan(testCase, n, 0, 'fixture must produce a non-empty grid');
for f = {'crit', 'r_minus_1', 'Delta', 'Dq_min', 'ref_denom', 'ref_margin', ...
         'gstat_local_denom', ...
         'omit_mu3', 'omit_cubic', 'omit_max'}
    verifyEqual(testCase, numel(prof.(f{1})), n, f{1});
end
verifyEqual(testCase, numel(prof.medium_status), n);
verifyEqual(testCase, numel(prof.node_term_reason), n);
verifyTrue(testCase, iscell(prof.medium_status));
verifyTrue(testCase, all(strcmp(prof.medium_status, 'ref_denom_small')));
verifyTrue(testCase, all(strcmp(prof.node_term_reason, 'medium_out_of_domain')));
verifyEqual(testCase, prof.ref_margin, prof.ref_denom - 1e9, 'RelTol', 1e-9, ...
    'ref_margin must be the actual signed distance to the configured floor');
end

% Step 7 schema coverage, BARE-SHORTCUT kind: opts.force_bare=true routes every node
% through eval_node's fbare branch (measured: status='ok' with a finite hstar; r==1,
% Sigma0==0, K0==0, term_reason='bare_shortcut', medium_status='not_applicable' on all
% 33 nodes; crit stays NaN because G0bare is never computed on this path). Same
% structural guarantee as test_per_node_arrays_are_present_and_aligned, plus the
% bare-shortcut-specific field values the implementer's own probe first measured.
function test_schema_bare_shortcut_record(testCase)
[ion, T, Bx, Jnu, o] = fixture();
o.force_bare = true;
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
n = numel(prof.hgrid);
verifyGreaterThan(testCase, n, 0, 'fixture must produce a non-empty grid');
for f = {'crit', 'r_minus_1', 'Delta', 'Dq_min', 'ref_denom', 'ref_margin', ...
         'gstat_local_denom', ...
         'omit_mu3', 'omit_cubic', 'omit_max'}
    verifyEqual(testCase, numel(prof.(f{1})), n, f{1});
end
verifyEqual(testCase, numel(prof.medium_status), n);
verifyEqual(testCase, numel(prof.node_term_reason), n);
verifyTrue(testCase, iscell(prof.medium_status));
verifyTrue(testCase, all(strcmp(prof.node_term_reason, 'bare_shortcut')));
verifyTrue(testCase, all(strcmp(prof.medium_status, 'not_applicable')));
verifyTrue(testCase, all(prof.r == 1));
verifyTrue(testCase, all(prof.Sigma0 == 0));
verifyTrue(testCase, all(prof.K0 == 0));
verifyTrue(testCase, all(isnan(prof.crit)));
verifyTrue(testCase, all(isfinite(prof.Delta)));
end

% crit is the dimensionless mass r + J0eff*G0bare. The predictor identity is tested at h=0
% directly; do not compare it to the first nonzero geometric node with an arbitrary 5% band.
function test_crit_definition_and_slope0_consistency(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.crit, prof.r + o.J0eff*prof.G0bare, 'RelTol', 1e-12);
verifyEqual(testCase, prof.slope0, prof.r_pm0 + o.J0eff*prof.G0bare_pm0, 'RelTol', 1e-12);
end

% r_minus_1 is r - 1 -- NOT Sigma0. They coincide only at m = 0 (spec §A Consequence 2), so a
% finite-moment node must show them differing, or the whole (r-1) measurement is vacuous.
function test_r_minus_1_is_not_sigma0_at_finite_moment(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.r_minus_1, prof.r - 1, 'AbsTol', 0);
k = find(abs(prof.m) > 1e-3 & isfinite(prof.Sigma0) & isfinite(prof.r), 1, 'last');
verifyNotEmpty(testCase, k, 'fixture produced no finite-moment converged node');
verifyGreaterThan(testCase, abs(prof.r_minus_1(k) - prof.Sigma0(k)), 1e-9, ...
    'r - 1 must differ from Sigma0 at finite m (spec §A Consequence 2)');
end

% Both path integrals use the SAME first-panel seeding as h0, so the three are
% quadrature-consistent by construction.
function test_path_integrals_match_the_h0_quadrature(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[~, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, prof.status, 'ok', 'fixture must be a deterministic accepted profile');
verifyEqual(testCase, prof.int_r_minus_1, ...
    trapz([0 prof.hgrid], [prof.r_pm0 - 1, prof.r - 1]), ...
    'RelTol', 1e-9, 'int_r_minus_1 must use the actual h=0 r seed');
verifyEqual(testCase, prof.int_Sigma0, ...
    trapz([0 prof.hgrid], [prof.Sigma0_pm0, prof.Sigma0]), 'RelTol', 1e-9);
end

% crit_star requires the root call to stop discarding G0bare.
function test_crit_star_at_the_root(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[hstar, prof] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfinite(hstar), 'fixture must produce the preregistered root');
verifyTrue(testCase, isfinite(prof.G0bare_star));
verifyEqual(testCase, prof.crit_star, prof.r_star + o.J0eff*prof.G0bare_star, 'RelTol', 1e-12);
% the accepted root must be the INCREASING crossing of F, i.e. a Landau minimum
verifyGreaterThan(testCase, prof.crit_star, 0);
end

% Supplying Jmom must not change any result (hot-path optimisation, not a semantic change).
function test_supplied_jmom_is_behaviour_neutral(testCase)
[ion, T, Bx, Jnu, o] = fixture();
[h1, p1] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
o.Jmom = invz_coupling_moments(Jnu);
[h2, p2] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
verifyEqual(testCase, h2, h1, 'AbsTol', 0);
verifyEqual(testCase, p2.r, p1.r, 'AbsTol', 0);
verifyEqual(testCase, p2.crit, p1.crit, 'AbsTol', 0);
end
