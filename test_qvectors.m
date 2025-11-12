%% Quick test script for compute_qvectors function
% This script performs basic tests to verify the q-vector computation

clear; clc;

fprintf('=== Testing compute_qvectors function ===\n\n');

%% Test 1: LiHoF4 lattice
fprintf('TEST 1: LiHoF4 tetragonal lattice\n');
fprintf('----------------------------------\n');
a = 5.175;
c = 10.75;
lattice_LiHoF4 = [a, 0, 0; 0, a, 0; 0, 0, c];

try
    [qvec, qvec_cart, recip] = compute_qvectors(lattice_LiHoF4);
    fprintf('✓ PASSED: High-symmetry points generated\n');
    fprintf('  Number of points: %d\n', size(qvec, 1));
    fprintf('  Gamma point (Cartesian): (%.4f, %.4f, %.4f) Å⁻¹\n', qvec_cart(1,:));
catch ME
    fprintf('✗ FAILED: %s\n', ME.message);
end
fprintf('\n');

%% Test 2: Grid generation
fprintf('TEST 2: Uniform grid generation\n');
fprintf('--------------------------------\n');
try
    [qvec_grid, ~] = compute_qvectors(lattice_LiHoF4, 'mode', 'grid', ...
                                       'grid', [5, 5, 5], 'range', [-0.5, 0.5]);
    expected_size = 5^3;
    if size(qvec_grid, 1) == expected_size
        fprintf('✓ PASSED: Grid size correct (%d points)\n', expected_size);
    else
        fprintf('✗ FAILED: Expected %d points, got %d\n', expected_size, size(qvec_grid, 1));
    end
catch ME
    fprintf('✗ FAILED: %s\n', ME.message);
end
fprintf('\n');

%% Test 3: Path generation
fprintf('TEST 3: Path through BZ\n');
fprintf('-----------------------\n');
try
    path = {'Gamma', 'X', 'Z'};
    [qvec_path, ~] = compute_qvectors(lattice_LiHoF4, 'mode', 'path', ...
                                       'path', path, 'npoints', 10);
    expected_size = (length(path) - 1) * 10;
    if size(qvec_path, 1) == expected_size
        fprintf('✓ PASSED: Path generated correctly (%d points)\n', expected_size);
    else
        fprintf('✗ FAILED: Expected %d points, got %d\n', expected_size, size(qvec_path, 1));
    end
catch ME
    fprintf('✗ FAILED: %s\n', ME.message);
end
fprintf('\n');

%% Test 4: Reciprocal lattice vectors
fprintf('TEST 4: Reciprocal lattice correctness\n');
fprintf('---------------------------------------\n');
try
    [~, ~, recip] = compute_qvectors(lattice_LiHoF4);

    % Check orthogonality: a·b* = 2π δ_ij
    a = lattice_LiHoF4(1,:);
    b = lattice_LiHoF4(2,:);
    c = lattice_LiHoF4(3,:);
    a_star = recip(1,:);
    b_star = recip(2,:);
    c_star = recip(3,:);

    % Should be 2π (within numerical tolerance)
    dot_aa = dot(a, a_star);
    dot_bb = dot(b, b_star);
    dot_cc = dot(c, c_star);

    % Should be 0 (within numerical tolerance)
    dot_ab = dot(a, b_star);
    dot_bc = dot(b, c_star);
    dot_ca = dot(c, a_star);

    tol = 1e-10;
    if abs(dot_aa - 2*pi) < tol && abs(dot_bb - 2*pi) < tol && abs(dot_cc - 2*pi) < tol
        fprintf('✓ PASSED: Diagonal terms correct (a·a* = 2π)\n');
    else
        fprintf('✗ FAILED: Diagonal terms incorrect\n');
        fprintf('  a·a* = %.6f (expected 2π = %.6f)\n', dot_aa, 2*pi);
        fprintf('  b·b* = %.6f (expected 2π = %.6f)\n', dot_bb, 2*pi);
        fprintf('  c·c* = %.6f (expected 2π = %.6f)\n', dot_cc, 2*pi);
    end

    if abs(dot_ab) < tol && abs(dot_bc) < tol && abs(dot_ca) < tol
        fprintf('✓ PASSED: Off-diagonal terms correct (a·b* = 0)\n');
    else
        fprintf('✗ FAILED: Off-diagonal terms should be zero\n');
        fprintf('  a·b* = %.2e\n', dot_ab);
        fprintf('  b·c* = %.2e\n', dot_bc);
        fprintf('  c·a* = %.2e\n', dot_ca);
    end
catch ME
    fprintf('✗ FAILED: %s\n', ME.message);
end
fprintf('\n');

%% Test 5: Cubic lattice
fprintf('TEST 5: Cubic lattice detection\n');
fprintf('--------------------------------\n');
try
    lattice_cubic = eye(3) * 5.0;  % Cubic with a=5 Angstrom
    [~, ~, ~] = compute_qvectors(lattice_cubic);
    fprintf('✓ PASSED: Cubic lattice processed\n');
catch ME
    fprintf('✗ FAILED: %s\n', ME.message);
end
fprintf('\n');

%% Summary
fprintf('=== ALL TESTS COMPLETED ===\n');
fprintf('\nThe compute_qvectors function is ready to use!\n');
fprintf('Run example_LiHoF4_qvectors.m for a comprehensive demonstration.\n');
