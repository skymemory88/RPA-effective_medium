function tests = test_invz_medium_moment_closure
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function mom = fixture_mom()
% production-multiset values (spec §B), used as plain numbers so this stays a unit test
mom = struct('Jbar', 1.20766e-4, 'mu2', 5.48264e-6, 'mu3', -3.42228e-11, ...
             'mu4', 2.3894*5.48264e-6^2, 'n', 16384);
end

function test_one_shot_formula(testCase)
mom = fixture_mom();  Gref = -200;
[K0, info] = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, K0, mom.Jbar - mom.mu2*Gref, 'AbsTol', 0);
verifyEqual(testCase, info.Kstrict, K0, 'AbsTol', 0);   % identical for one-shot, by definition
verifyEqual(testCase, info.retained, 'mu2');
verifyEqual(testCase, info.status, 'ok');
verifyEqual(testCase, info.scheme, 'strict_1z_dyson_ref');
end

% Gref < 0 (G = -chi, chi > 0) so the mu2 term RAISES K0 above Jbar. A sign flip here would
% invert the whole medium correction, so pin the direction explicitly.
function test_correction_sign(testCase)
mom = fixture_mom();
K0 = invz_medium_moment_closure(-200, mom, 'strict_1z_dyson_ref');
verifyGreaterThan(testCase, K0, mom.Jbar);
end

% mu3*Gref^2 is the FIRST omitted term, before the cubic. Both ratios are always reported --
% mu3's near-zero value is a measured property of one multiset, never a general licence.
function test_both_omitted_ratios_reported(testCase)
mom = fixture_mom();  Gref = -200;
[~, info] = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, info.omit_mu3, abs(mom.mu3*Gref^2)/abs(mom.mu2*Gref), 'RelTol', 1e-12);
verifyEqual(testCase, info.omit_cubic, ...
    abs((2*mom.mu2^2 - mom.mu4)*Gref^3)/abs(mom.mu2*Gref), 'RelTol', 1e-12);
verifyEqual(testCase, info.omit_max, max(info.omit_mu3, info.omit_cubic), 'AbsTol', 0);
end

% A skewed multiset must make omit_mu3 DOMINATE, proving the ratio is not decoration.
function test_skewed_multiset_makes_mu3_dominate(testCase)
mom = fixture_mom();
mom.mu3 = 100*mom.mu2^1.5;                 % strongly skewed, unlike the production multiset
[~, info] = invz_medium_moment_closure(-200, mom, 'strict_1z_dyson_ref');
verifyGreaterThan(testCase, info.omit_mu3, info.omit_cubic);
verifyEqual(testCase, info.omit_max, info.omit_mu3, 'AbsTol', 0);
end

% Frozen thresholds (prereg §4) are the CALLER's gate, not this leaf's: the leaf never rejects.
function test_leaf_never_rejects_on_large_ratio(testCase)
mom = fixture_mom();
[K0, info] = invz_medium_moment_closure(-5000, mom, 'strict_1z_dyson_ref');
verifyGreaterThan(testCase, info.omit_max, 0.25);       % beyond omit_promote
verifyEqual(testCase, info.status, 'ok');               % still 'ok': the polynomial is defined
verifyTrue(testCase, isfinite(K0));
end

% Explicit zero convention (spec §4.1): 0/0 -> 0 only when the numerator vanishes too.
function test_zero_denominator_convention(testCase)
mom = struct('Jbar', 1e-4, 'mu2', 5e-6, 'mu3', 0, 'mu4', 0, 'n', 10);
[~, a] = invz_medium_moment_closure(0, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, a.omit_mu3, 0, 'AbsTol', 0);      % numerator also 0
verifyEqual(testCase, a.omit_cubic, 0, 'AbsTol', 0);
mom.mu3 = 1e-11;
[~, b] = invz_medium_moment_closure(0, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, b.omit_mu3, 0, 'AbsTol', 0);      % Gref = 0 => numerator 0 as well
mom2 = struct('Jbar', 1e-4, 'mu2', 0, 'mu3', 1e-11, 'mu4', 1e-14, 'n', 10);
[~, c] = invz_medium_moment_closure(-200, mom2, 'strict_1z_dyson_ref');
verifyEqual(testCase, c.omit_mu3, Inf);                 % mu2 = 0, numerator nonzero
end

function test_nonfinite_gref_returns_status(testCase)
mom = fixture_mom();
[K0, info] = invz_medium_moment_closure(NaN, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, info.status, 'nonfinite');
verifyTrue(testCase, isnan(K0));
end

function test_resummed_and_unknown_schemes_are_wiring_errors(testCase)
mom = fixture_mom();
verifyError(testCase, @() invz_medium_moment_closure(-200, mom, 'resummed'), ...
            'invz:staticMedium');
verifyError(testCase, @() invz_medium_moment_closure(-200, mom, 'nope'), 'invz:staticMedium');
end

% Non-scalar moment fields are a wiring error: the static slot must pass column 1, not the
% whole [nJ,nw] moment set (spec §4.3).
function test_nonscalar_moments_rejected(testCase)
mom = struct('Jbar', [1e-4 2e-4], 'mu2', [5e-6 6e-6], 'mu3', [0 0], 'mu4', [0 0], 'n', [10 10]);
verifyError(testCase, @() invz_medium_moment_closure(-200, mom, 'strict_1z_dyson_ref'), ...
            'invz:staticMedium');
end
