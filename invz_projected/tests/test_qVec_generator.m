function tests = test_qVec_generator
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));
addpath(fullfile(here,'..','..'));
end

function lat = tetragonal()
lat = [5.175 0 0; 0 5.175 0; 0 0 10.75];
end

function test_default_grid_is_inclusive_and_unchanged(testCase)
% Default endpoint=true must reproduce the historical inclusive linspace grid
% exactly (the LiHoF4 benchmarks depend on it).
[q, ~] = evalc_grid(4, true);
axis = unique(q(:,1)).';
verifyEqual(testCase, axis, linspace(-0.5, 0.5, 4), 'AbsTol', 1e-14);
verifyEqual(testCase, size(q,1), 4^3);
% both boundary faces present (the redundancy finding #5 describes)
verifyTrue(testCase, any(abs(q(:,1) - (-0.5)) < 1e-12) && any(abs(q(:,1) - 0.5) < 1e-12));
end

function test_halfopen_grid_has_no_duplicate_face(testCase)
% endpoint=false gives N distinct points per axis with no +0.5 duplicate of -0.5.
N = 4;
[q, ~] = evalc_grid(N, false);
axis = unique(q(:,1)).';
verifyEqual(testCase, numel(axis), N);                       % N distinct, not N with a dup face
verifyEqual(testCase, axis, -0.5 + (0:N-1)/N, 'AbsTol', 1e-14);
verifyFalse(testCase, any(abs(q(:,1) - 0.5) < 1e-12));       % upper face excluded
verifyEqual(testCase, size(q,1), N^3);
end

function [q, qc] = evalc_grid(N, endpoint)
% Run the grid generator quietly (verbose=false) and capture q-vectors.
q  = qVec_generator(tetragonal(), 'mode', 'grid', 'grid', [N N N], ...
                    'range', [-0.5 0.5], 'verbose', false, 'endpoint', endpoint);
qc = [];
end

function test_getf_shared_helper(testCase)
% The centralized invz/getf.m returns the field value or the default (finding #8).
s = struct('a', 7);
verifyEqual(testCase, getf(s, 'a', 99), 7);
verifyEqual(testCase, getf(s, 'b', 99), 99);
end
