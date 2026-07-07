function tests = test_build_jq
% TEST_BUILD_JQ Tests for real-space to reciprocal J(q) builder.

tests = functiontests(localfunctions);

end

function setupOnce(~)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

function test_build_jq_two_point_phase(testCase)
qvec = [0 0 0; 0.5 0 0];
Rvec = [0 0 0; 1 0 0];

A = diag([1, 2, 3]);
B = [0.2 0.05 0; 0.05 0.1 0; 0 0 0.3];
J_of_R = zeros(3, 3, 2);
J_of_R(:,:,1) = A;
J_of_R(:,:,2) = B;

Jq = build_jq(qvec, Rvec, J_of_R, struct('symmetrize', true));

verifyEqual(testCase, Jq(:,:,1), A + B, 'AbsTol', 1e-12);
verifyEqual(testCase, Jq(:,:,2), A - B, 'AbsTol', 1e-12);
end

function test_validate_jq_input(testCase)
Jq = zeros(3, 3, 4);
for iq = 1:4
    J = [1 0.1 0; 0.1 2 0; 0 0 3] + 0.01 * iq * eye(3);
    Jq(:,:,iq) = J;
end

report = validate_jq_input(Jq);
verifyTrue(testCase, report.valid);
verifyLessThan(testCase, report.max_asymmetry, 1e-12);
end
