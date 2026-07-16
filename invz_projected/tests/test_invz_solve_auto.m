function tests = test_invz_solve_auto
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_expected_invz_errors_become_phase0(testCase)
% Known numerical conditions (invz:* identifiers) are absorbed as "no solution" with
% the identifier preserved in the diagnostic. At T=1.9 K, Bx=0.04 T the ordered branch
% does not converge, so invz_solve_auto falls through to the paramagnetic solve, whose
% two-level helper raises invz:degenerateDoublet at this small field.
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[pt, phase, di] = invz_solve_auto(ion, 1.9, 0.04, Jnu, struct('J0eff', 6.4e-3));
verifyEqual(testCase, phase, 0);
verifyEmpty(testCase, pt);
verifyEqual(testCase, di.para_err, 'invz:degenerateDoublet');
end

function test_unexpected_errors_rethrow(testCase)
% Programming/API defects must NOT be converted into phase = 0: a malformed ion
% struct (missing field) raises a non-invz error that must propagate.
ion = rmfield(invz_ion(), 'J');
Jnu = linspace(-2e-3, 6.0e-3, 24).';
verifyError(testCase, ...
    @() invz_solve_auto(ion, 0.31, 5.5, Jnu, struct('J0eff', 6.4e-3)), ...
    'MATLAB:nonExistentField');
end
