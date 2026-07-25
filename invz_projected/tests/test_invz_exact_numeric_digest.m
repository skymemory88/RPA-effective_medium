function tests = test_invz_exact_numeric_digest
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_digest_is_64char_lowercase_hex_row(testCase)
d = invz_exact_numeric_digest([1 2 3]);
verifyEqual(testCase, size(d), [1 64]);
verifyTrue(testCase, ischar(d));
verifyTrue(testCase, all(ismember(d, '0123456789abcdef')));
end

function test_digest_is_stable_across_repeated_calls(testCase)
x = [1 2 3; 4 5 6];
verifyEqual(testCase, invz_exact_numeric_digest(x), invz_exact_numeric_digest(x));
end

% NOT shape-invariant by design: a reshaped multiset is a different input.
function test_digest_is_shape_sensitive(testCase)
verifyNotEqual(testCase, invz_exact_numeric_digest([1 2 3]), ...
    invz_exact_numeric_digest([1 2 3].'));
end

% NOT order-invariant by design: a reordered multiset is a different input.
function test_digest_is_order_sensitive(testCase)
verifyNotEqual(testCase, invz_exact_numeric_digest([1 2 3]), ...
    invz_exact_numeric_digest([3 2 1]));
end

function test_digest_is_class_sensitive(testCase)
verifyNotEqual(testCase, invz_exact_numeric_digest(int32([1 2 3])), ...
    invz_exact_numeric_digest(double([1 2 3])));
end

% REGRESSION on the "exact-byte" claim itself. These two int64 values are distinct but collapse
% to the SAME double (2^53 is the last integer with an exact double neighbour), so any
% implementation that hashes double(x) instead of the original class bytes digests them
% identically -- same class, same shape, same converted data. Hashing x(:) directly separates
% them. Without this test "exact-byte" would be true only for doubles.
function test_digest_hashes_original_class_bytes_not_doubles(testCase)
a = int64(9007199254740992);        % 2^53
b = int64(9007199254740993);        % 2^53 + 1; double(a) == double(b)
verifyEqual(testCase, double(a), double(b), 'AbsTol', 0);   % the collision is real
verifyNotEqual(testCase, invz_exact_numeric_digest(a), invz_exact_numeric_digest(b));
end

function test_nonnumeric_and_sparse_inputs_raise_exactDigest(testCase)
verifyError(testCase, @() invz_exact_numeric_digest('abc'), 'invz:exactDigest');
verifyError(testCase, @() invz_exact_numeric_digest(sparse([1 0 2])), 'invz:exactDigest');
verifyError(testCase, @() invz_exact_numeric_digest([1 2i]), 'invz:exactDigest');
end
