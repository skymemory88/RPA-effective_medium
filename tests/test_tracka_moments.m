function tests = test_tracka_moments
% TEST_TRACKA_MOMENTS Unit tests for Track-A lambda/X helpers.

tests = functiontests(localfunctions);

end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

function test_compute_lambdas_exact(testCase)
g = [1; 2; 3];
k = [2; 2; 2];
beta = 2;
p = [1 2];

lambda = emt_compute_lambdas(k, g, beta, p, []);

verifyEqual(testCase, lambda(1), 6, 'AbsTol', 1e-12);
verifyEqual(testCase, lambda(2), 14, 'AbsTol', 1e-12);
end

function test_compute_x_moment_ratio(testCase)
g = [1; 2];
lambda = [1, 2];
track = struct('model', 'moment_ratio', 'clamp_scale', 100);

[X, aux] = emt_compute_x_from_lambdas(g, lambda, track);

verifyEqual(testCase, X(1), 2, 'AbsTol', 1e-12);
verifyEqual(testCase, X(2), 3, 'AbsTol', 1e-12);
verifyEqual(testCase, aux.a, 1, 'AbsTol', 1e-12);
verifyEqual(testCase, aux.gamma, 1, 'AbsTol', 1e-12);
end

function test_compute_x_custom_callback(testCase)
g = [1; 2; 3];
lambda = [2, 4];
track = struct();
track.clamp_scale = 100;
track.custom_x_function = @custom_model;

[X, aux] = emt_compute_x_from_lambdas(g, lambda, track, []);

verifyEqual(testCase, X(:), [2; 2; 2], 'AbsTol', 1e-12);
verifyEqual(testCase, aux.a, 2, 'AbsTol', 1e-12);
verifyEqual(testCase, aux.model, 'custom');
end

function test_compute_x_jensen_216_219(testCase)
g = [1; 2];
k = [0.1; 0.2];
lambda = [3, 5];

track = struct();
track.model = 'jensen_216_219';
track.clamp_scale = 100;
track.beta = 2;
track.M2 = 4;
track.n01 = 0.5;

[X, aux] = emt_compute_x_from_lambdas(g, lambda, track, k);

% alpha = (4/0.25) * (5 - 0.5*(1 + 2*(1-0.25))*3)
%       = 16 * (5 - 3.75) = 20
verifyEqual(testCase, aux.a, 20, 'AbsTol', 1e-12);
verifyEqual(testCase, aux.lambda1, 3, 'AbsTol', 1e-12);
verifyEqual(testCase, aux.lambda2, 5, 'AbsTol', 1e-12);

% gamma(i) = (4/0.25) * (3 - 0.75*k(i)) = 16 * (3 - 0.75*k(i))
gamma_expected = 16 * ([3; 3] - 0.75 * k);
verifyEqual(testCase, aux.gamma(:), gamma_expected(:), 'AbsTol', 1e-12);

X_expected = 20 + gamma_expected .* g;
verifyEqual(testCase, X(:), X_expected(:), 'AbsTol', 1e-12);
end

function [X, aux] = custom_model(g, lambda, ~, ~)
X = 0 * g + lambda(1);
aux = struct('a', lambda(1), 'gamma', 0);
end
