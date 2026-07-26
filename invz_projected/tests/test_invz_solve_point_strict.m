function tests = test_invz_solve_point_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [ion, T, Bx, Jnu, o] = fx()
ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
end

% Provenance is mandatory on BOTH legs: an unlabelled result cannot be compared across schemes.
function test_pm_leg_carries_provenance(testCase)
[ion, T, Bx, Jnu, o] = fx();
pt = invz_solve_point(ion, T, Bx, Jnu, o);
verifyEqual(testCase, pt.static_medium, 'resummed');
verifyTrue(testCase, isfield(pt, 'Jmom') && isfield(pt, 'medium_margin'));
o.static_medium = 'strict_1z_dyson_ref';
pts = invz_solve_point(ion, T, Bx, Jnu, o);
verifyEqual(testCase, pts.static_medium, 'strict_1z_dyson_ref');
verifyEqual(testCase, pts.Jmom.mu2, invz_coupling_moments(Jnu).mu2, 'AbsTol', 0);
end

% Default absent => legacy numbers unchanged (G9 at point level).
function test_absent_scheme_reproduces_legacy_numbers(testCase)
[ion, T, Bx, Jnu, o] = fx();
a = invz_solve_point(ion, T, Bx, Jnu, o);
b = invz_solve_point(ion, T, Bx, Jnu, setfield(o, 'static_medium', 'resummed'));  %#ok<SFLD>
verifyEqual(testCase, b.Sigma0, a.Sigma0, 'AbsTol', 0);
verifyEqual(testCase, b.crit,   a.crit,   'AbsTol', 0);
end

% Strict changes the PM mass -- it must, since K(0) changed (spec §0.2).
function test_strict_shifts_the_pm_mass(testCase)
[ion, T, Bx, Jnu, o] = fx();
a = invz_solve_point(ion, T, Bx, Jnu, o);
b = invz_solve_point(ion, T, Bx, Jnu, setfield(o, 'static_medium', 'strict_1z_dyson_ref'));  %#ok<SFLD>
verifyNotEqual(testCase, b.crit, a.crit);
verifyTrue(testCase, isfinite(b.crit));
end

% The ordered leg exposes the two tiers SEPARATELY, and converged keeps its old meaning.
function test_ordered_leg_two_tier_fields(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.ordered_mode = 'jensen';  o.static_medium = 'strict_1z_dyson_ref';
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfield(pt, 'stable_1z') && isfield(pt, 'crit_1z'));
verifyTrue(testCase, isfield(pt, 'omit_max'));
verifyEqual(testCase, pt.static_medium, 'strict_1z_dyson_ref');
if pt.is_ordered && isfinite(pt.crit_1z)
    Dtol = 1e-6*max(1, abs(pt.G(1))*max(abs(Jnu)));
    expected = pt.converged && pt.crit_1z > 1e-6 && ...
               pt.D_uni > Dtol && pt.Dq_min > Dtol;
    verifyEqual(testCase, pt.stable_1z, expected);
end
end

% A per-leg override is rejected at the public entry, so the sectors cannot split.
function test_per_leg_override_rejected(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.emt = struct('static_medium', 'strict_1z_dyson_ref');
verifyError(testCase, @() invz_solve_point(ion, T, Bx, Jnu, o), 'invz:staticMedium');
end

% Review finding 1: the BARE ordered loop -- the DEFAULT ordered_mode, reachable from
% invz_solve_auto -- must halt at the FIRST strict-scheme domain event, like the PM loop.
% Pre-fix this ran to max_outer = 200 and exported 'nonfinite', because the NaN Sigma fed the
% next iteration's reference: provenance naming the WRONG condition. The reference denominator
% on iteration 1 is 1 + Sigma(1) = 1 (this leg always cold-starts at Sigma = 0), so a floor of
% 2 makes the true first cause 'ref_denom_small' on iteration 1, deterministically.
function test_bare_strict_halt_reports_the_first_cause(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.static_medium = 'strict_1z_dyson_ref';  o.ref_margin = 2;
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, pt.is_ordered);
verifyFalse(testCase, pt.converged);
verifyEqual(testCase, pt.medium_status, 'ref_denom_small');
verifyEqual(testCase, pt.outer_iters, 1);
end

% The halt path must export the SAME field set as the in-domain path (callers never probe a
% missing member) -- the property already measured for the PM solver.
function test_bare_halt_exports_the_full_field_set(testCase)
[ion, T, Bx, Jnu, o] = fx();
o.static_medium = 'strict_1z_dyson_ref';
ok   = invz_solve_point_ordered(ion, T, Bx, Jnu, o);
halt = invz_solve_point_ordered(ion, T, Bx, Jnu, setfield(o, 'ref_margin', 2));  %#ok<SFLD>
verifyEqual(testCase, ok.medium_status, 'ok');
verifyEqual(testCase, halt.medium_status, 'ref_denom_small');
verifyEqual(testCase, sort(fieldnames(halt)), sort(fieldnames(ok)));
end

% The guard is INERT under the default 'resummed' scheme: invz_emt_scalar's strict block is
% gated off there, so medium_status is always 'not_applicable' -- one of the two strings the
% guard whitelists -- and the halt can never be reached.
function test_resummed_bare_never_reaches_the_guard(testCase)
[ion, T, ~, Jnu, o] = fx();
for B = [0.5 2.85]
    pt = invz_solve_point_ordered(ion, T, [B 0 0], Jnu, o);
    verifyEqual(testCase, pt.medium_status, 'not_applicable');
end
end

% Review finding 2, the docstring's field-presence contract that task 15's spectra map reads:
% the jensen-only members are ABSENT ENTIRELY in bare mode (not NaN-valued), while the five
% static-medium provenance members are present on both modes.
function test_bare_mode_omits_the_jensen_only_members(testCase)
[ion, T, Bx, Jnu, o] = fx();
jonly = {'stable_1z', 'crit_1z', 'Dq_min', 'D_uni', ...
         'omit_mu3', 'omit_cubic', 'omit_max', 'path_omit_max'};
both = {'static_medium', 'medium_status', 'medium_denom', 'medium_margin', 'Jmom'};
pt = invz_solve_point_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, pt.is_ordered);
for f = jonly
    verifyFalse(testCase, isfield(pt, f{1}), ...
        sprintf('bare mode must not export the jensen-only member %s', f{1}));
end
for f = both
    verifyTrue(testCase, isfield(pt, f{1}), ...
        sprintf('%s must be present on both ordered modes', f{1}));
end
end
