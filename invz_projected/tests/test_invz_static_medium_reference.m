function tests = test_invz_static_medium_reference
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_dyson_ref_divides_by_one_plus_sigma(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.denom, 1.25, 'AbsTol', 0);
verifyEqual(testCase, Gref, -0.5/1.25, 'AbsTol', 0);
verifyEqual(testCase, ref.status, 'ok');
verifyEqual(testCase, ref.scheme, 'strict_1z_dyson_ref');
end

function test_bare_ref_denominator_is_one(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_bare_ref');
verifyEqual(testCase, ref.denom, 1, 'AbsTol', 0);
verifyEqual(testCase, Gref, -0.5, 'AbsTol', 0);
verifyEqual(testCase, ref.status, 'ok');
end

% The two conventions coincide only as Sigma0 -> 0; that is the O(1/z^2) difference the spec
% labels a Dyson-improved scheme choice, and it must be visible, not hidden.
function test_conventions_differ_at_finite_sigma_and_agree_at_zero(testCase)
gD = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_dyson_ref');
gB = invz_static_medium_reference(-0.5, 0.25, 'strict_1z_bare_ref');
verifyGreaterThan(testCase, abs(gD - gB)/abs(gB), 0.1);
g0D = invz_static_medium_reference(-0.5, 0, 'strict_1z_dyson_ref');
g0B = invz_static_medium_reference(-0.5, 0, 'strict_1z_bare_ref');
verifyEqual(testCase, g0D, g0B, 'AbsTol', 0);
end

% Domain events RETURN a status; they never throw (spec §5.2).
function test_nonpositive_denominator_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, -1.5, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'ref_denom_nonpositive');
verifyTrue(testCase, isnan(Gref));
verifyEqual(testCase, ref.denom, -0.5, 'AbsTol', 1e-15);
end

function test_small_denominator_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, -1 + 1e-9, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'ref_denom_small');
verifyTrue(testCase, isnan(Gref));
verifyEqual(testCase, ref.floor, 1e-6, 'AbsTol', 0);
verifyLessThanOrEqual(testCase, ref.margin, 0);
end

function test_margin_is_configurable_and_boundary_is_inclusive(testCase)
[~, ref] = invz_static_medium_reference(-0.5, -1 + 1e-3, 'strict_1z_dyson_ref', ...
                                       struct('ref_margin', 1e-2));
verifyEqual(testCase, ref.status, 'ref_denom_small');
[~, ok] = invz_static_medium_reference(-0.5, -1 + 1e-1, 'strict_1z_dyson_ref', ...
                                       struct('ref_margin', 1e-2));
verifyEqual(testCase, ok.status, 'ok');
verifyGreaterThan(testCase, ok.margin, 0);
end

function test_nonfinite_input_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(NaN, 0.1, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'nonfinite');
verifyTrue(testCase, isnan(Gref));
end

% 'resummed' must never reach this primitive: that is a wiring error, not a domain event.
function test_resummed_scheme_is_a_wiring_error(testCase)
verifyError(testCase, @() invz_static_medium_reference(-0.5, 0.1, 'resummed'), ...
            'invz:staticMedium');
verifyError(testCase, @() invz_static_medium_reference(-0.5, 0.1, 'nonsense'), ...
            'invz:staticMedium');
end
