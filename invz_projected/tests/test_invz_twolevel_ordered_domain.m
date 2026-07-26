function tests = test_invz_twolevel_ordered_domain
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Default is unchanged: exactly B = 0, hz = 0 still throws, deliberately (README.html:208).
function test_default_still_throws_at_zero_field(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_twolevel_ordered(ion, 0.31, [0 0 0], 0, ...
    struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x')), 'invz:degenerateDoublet');
end

function test_valid_flag_present_on_the_normal_path(testCase)
ion = invz_ion();
tl = invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, ...
    struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x'));
verifyTrue(testCase, tl.valid);
verifyGreaterThan(testCase, tl.Delta, 1e-4);
verifyTrue(testCase, isfield(tl, 'M2') && isfield(tl, 'g0'));
end

% Return mode reports the domain outcome instead of throwing, and reports the Delta it measured.
function test_return_mode_reports_instead_of_throwing(testCase)
ion = invz_ion();
o = struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x', 'domain_policy', 'return');
tl = invz_twolevel_ordered(ion, 0.31, [0 0 0], 0, o);
verifyFalse(testCase, tl.valid);
verifyLessThan(testCase, tl.Delta, 1e-4);
% the two-level fields are deliberately ABSENT: a consumer ignoring .valid must fail loudly
verifyFalse(testCase, isfield(tl, 'g0'));
end

% Return mode is behaviour-neutral where the doublet is resolved.
function test_return_mode_identical_when_valid(testCase)
ion = invz_ion();
base = struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x');
a = invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, base);
b = invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, ...
    setfield(base, 'domain_policy', 'return'));  %#ok<SFLD>
verifyEqual(testCase, b, a);
end

function test_unknown_policy_is_a_wiring_error(testCase)
ion = invz_ion();
verifyError(testCase, @() invz_twolevel_ordered(ion, 0.31, [2.85 0 0], 0.02, ...
    struct('Jxx0', ion.Jxx0, 'domain_policy', 'maybe')), 'invz:twolevelDomainPolicy');
end
