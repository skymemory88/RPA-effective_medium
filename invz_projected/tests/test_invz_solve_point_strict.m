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
