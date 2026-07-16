function tests = test_invz_field_vec
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
end

function test_scalar_maps_to_transverse_row(testCase)
verifyEqual(testCase, invz_field_vec(4), [4 0 0]);
verifyEqual(testCase, invz_field_vec(0), [0 0 0]);
verifyEqual(testCase, invz_field_vec(-2.5), [-2.5 0 0]);   % signed amplitude passes through
end

function test_vectors_normalize_to_row(testCase)
verifyEqual(testCase, invz_field_vec([1 2 3]), [1 2 3]);
verifyEqual(testCase, invz_field_vec([1; 2; 3]), [1 2 3]); % column -> row
verifyEqual(testCase, invz_field_vec(invz_field_vec(5)), [5 0 0]);  % idempotent
end

function test_invalid_inputs_error(testCase)
bad = {NaN, [1 Inf 0], 1+2i, [], [1 2], [1 2 3 4], 'abc', true};
for k = 1:numel(bad)
    verifyError(testCase, @() invz_field_vec(bad{k}), 'invz:fieldVec');
end
end
