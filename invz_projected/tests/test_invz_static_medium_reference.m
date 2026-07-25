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

% Neither sub-case above puts denom exactly at opts.ref_margin (1e-3 is strictly below the
% 1e-2 floor; 1e-1 is strictly above it), so that test cannot detect a regression from '<=' to
% '<' in the ref_denom_small classification. This test constructs denom bit-exactly equal to
% opts.ref_margin using exactly-representable doubles (denom = 1 + Sigma0 = 1 + (-0.75) = 0.25,
% opts.ref_margin = 0.25) to test the boundary itself, and asserts the exactness of denom
% itself so the test cannot silently drift into a near-miss and stop testing the boundary.
function test_boundary_exactly_at_margin_is_ref_denom_small(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, -0.75, 'strict_1z_dyson_ref', ...
                                          struct('ref_margin', 0.25));
verifyEqual(testCase, ref.denom, 0.25, 'AbsTol', 0);
verifyEqual(testCase, ref.status, 'ref_denom_small');
verifyTrue(testCase, isnan(Gref));
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

% G0bare must be validated as scalar: an un-indexed [nJ,nw] retardation set is the same
% "forgot to index down to the static column" wiring error the mom-field guard in
% invz_medium_moment_closure guards against, and this is the sole documented producer of Gref.
function test_nonscalar_g0bare_rejected(testCase)
verifyError(testCase, @() invz_static_medium_reference([-0.5 -0.6], 0.25, ...
            'strict_1z_dyson_ref'), 'invz:staticMedium');
end

% Symmetric guard on Sigma0.
function test_nonscalar_sigma0_rejected(testCase)
verifyError(testCase, @() invz_static_medium_reference(-0.5, [0.25 0.3], ...
            'strict_1z_dyson_ref'), 'invz:staticMedium');
end

% Non-finite Sigma0 must still return a status rather than throw -- the new scalar guard on
% Sigma0 checks scalarity and numeric type only, never finiteness, so NaN must flow through to
% the existing 'nonfinite' handling (here via a non-finite denom, distinct from the existing
% G0bare-NaN test above which is non-finite via G0bare itself).
function test_nonfinite_sigma0_returns_status(testCase)
[Gref, ref] = invz_static_medium_reference(-0.5, NaN, 'strict_1z_dyson_ref');
verifyEqual(testCase, ref.status, 'nonfinite');
verifyTrue(testCase, isnan(Gref));
end
