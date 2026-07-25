function tests = test_invz_coupling_moments
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% Population (N) normalization, NOT MATLAB's sample (N-1) var(): the difference is 4% at
% N = 24, exactly the size of the legacy synthetic fixtures.
function test_population_normalization_not_sample(testCase)
J = linspace(-2e-3, 6e-3, 24).';
mom = invz_coupling_moments(J);
Jb = mean(J);
verifyEqual(testCase, mom.Jbar, Jb, 'AbsTol', 0);
verifyEqual(testCase, mom.mu2, mean((J - Jb).^2), 'AbsTol', 0);
% and it is demonstrably NOT var()
verifyGreaterThan(testCase, abs(var(J) - mom.mu2)/mom.mu2, 0.03);
end

function test_third_and_fourth_moments(testCase)
J = [1e-3; -2e-3; 5e-3; 0; 3e-3];
Jb = mean(J);
mom = invz_coupling_moments(J);
verifyEqual(testCase, mom.mu3, mean((J - Jb).^3), 'AbsTol', 0);
verifyEqual(testCase, mom.mu4, mean((J - Jb).^4), 'AbsTol', 0);
verifyEqual(testCase, mom.n, numel(J));
end

% [nJ,nw] retardation interface: one moment set PER COLUMN, never one flattened multiset.
function test_matrix_input_is_per_column(testCase)
A = [1e-3; 2e-3; 3e-3];
B = [-4e-3; 0; 4e-3];
mom = invz_coupling_moments([A B]);
momA = invz_coupling_moments(A);
momB = invz_coupling_moments(B);
verifyEqual(testCase, size(mom.Jbar), [1 2]);
verifyEqual(testCase, mom.Jbar, [momA.Jbar momB.Jbar], 'AbsTol', 0);
verifyEqual(testCase, mom.mu2,  [momA.mu2  momB.mu2],  'AbsTol', 0);
verifyEqual(testCase, mom.n,    [momA.n    momB.n]);
% a flattened interpretation would give a single scalar equal to neither
flat = invz_coupling_moments([A; B]);
verifyNotEqual(testCase, flat.mu2, momA.mu2);
end

function test_rejects_nonfinite_and_complex(testCase)
verifyError(testCase, @() invz_coupling_moments([1e-3; NaN]), 'invz:couplingMoments');
verifyError(testCase, @() invz_coupling_moments([1e-3; Inf]), 'invz:couplingMoments');
verifyError(testCase, @() invz_coupling_moments([1e-3; 1i*1e-3]), 'invz:couplingMoments');
verifyError(testCase, @() invz_coupling_moments([]), 'invz:couplingMoments');
end

% Production-multiset regression, with its provenance asserted (spec §B: the numbers are valid
% ONLY for this tuple). INVZ_SLOW-gated like the other real-coupling anchors.
function test_production_multiset_moments(testCase)
if isempty(getenv('INVZ_SLOW')), assumeFail(testCase, 'set INVZ_SLOW=1 to run'); end
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid', [16 16 16], 'dpRng', 30, ...
    'dipole', 'bruteforce', 'cache', false));
mom = invz_coupling_moments(Jnu(:));
verifyEqual(testCase, mom.n, 16384);
verifyEqual(testCase, info.dipole.backend, 'bruteforce');
verifyFalse(testCase, isfield(info, 'grid'), ...
    'the frozen fixture uses the legacy absent-grid-policy route');
verifyEqual(testCase, info.Jcc0, 6.42444e-3, 'RelTol', 1e-4);
verifyEqual(testCase, mom.Jbar,  1.20766e-4, 'RelTol', 1e-4);
verifyEqual(testCase, mom.mu2,   5.48264e-6, 'RelTol', 1e-4);
verifyEqual(testCase, sqrt(mom.mu2), 2.3415e-3, 'RelTol', 1e-4);
end
