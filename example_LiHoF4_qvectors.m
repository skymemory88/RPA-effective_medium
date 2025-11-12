%% Example script: Computing q-vectors for LiHoF4
% This script demonstrates how to compute q-vectors for LiHoF4 using the
% compute_qvectors() function and visualize them in reciprocal space.
%
% LiHoF4 crystallizes in the scheelite structure (space group I41/a)
% - Tetragonal lattice
% - Body-centered tetragonal structure
% - 4 magnetic Ho³⁺ ions per unit cell
%
% Author: Claude (Anthropic)
% Date: 2025-11-12

clear; close all;

%% Define LiHoF4 lattice parameters
% From literature: J. Phys. C: Solid State Phys. 8, 3483 (1975)
% and other sources
a = 5.175;   % Å (tetragonal a-axis)
b = 5.175;   % Å (tetragonal b-axis, equal to a)
c = 10.75;   % Å (tetragonal c-axis)

% Construct real-space lattice matrix
% Each ROW is a lattice vector in Cartesian coordinates
lattice = [a,  0,  0;   % a-vector
           0,  b,  0;   % b-vector
           0,  0,  c];  % c-vector

fprintf('=== LiHoF4 Q-Vector Computation ===\n\n');

%% Example 1: High-symmetry points only (default)
fprintf('EXAMPLE 1: High-symmetry points\n');
fprintf('===============================\n');
[qvec_sym, qvec_sym_cart, recip_lattice] = compute_qvectors(lattice);
fprintf('\n');

%% Example 2: Path through reciprocal space
% Common path for tetragonal: Gamma -> X -> M -> Gamma -> Z -> A -> R
fprintf('EXAMPLE 2: Path through high-symmetry points\n');
fprintf('==============================================\n');
path_points = {'Gamma', 'X', 'M', 'Gamma', 'Z', 'A', 'R', 'Gamma'};
[qvec_path, qvec_path_cart] = compute_qvectors(lattice, ...
                                                'mode', 'path', ...
                                                'path', path_points, ...
                                                'npoints', 50);
fprintf('\n');

%% Example 3: Uniform grid for full q-space integration
fprintf('EXAMPLE 3: Uniform q-space grid\n');
fprintf('================================\n');
[qvec_grid, qvec_grid_cart] = compute_qvectors(lattice, ...
                                                'mode', 'grid', ...
                                                'grid', [11, 11, 11], ...
                                                'range', [-0.5, 0.5]);
fprintf('\n');

%% Example 4: Fine grid for RPA/effective medium calculations
fprintf('EXAMPLE 4: Fine grid for RPA calculations\n');
fprintf('=========================================\n');
[qvec_fine, qvec_fine_cart] = compute_qvectors(lattice, ...
                                                'mode', 'grid', ...
                                                'grid', [21, 21, 21], ...
                                                'range', [-1, 1]);
fprintf('\n');

%% Example 5: Line scan along [001] direction (c-axis)
fprintf('EXAMPLE 5: Line scan along [001] (c-axis)\n');
fprintf('==========================================\n');
qz_values = linspace(0, 1, 101)';
qvec_001_line = [zeros(length(qz_values), 1), zeros(length(qz_values), 1), qz_values];
qvec_001_cart = qvec_001_line * recip_lattice;
fprintf('Generated %d q-points along [001] direction\n', size(qvec_001_line, 1));
fprintf('Range: (0,0,0) to (0,0,1) in reduced coordinates\n');
fprintf('       (0,0,0) to (0,0,%.4f) Å⁻¹ in Cartesian\n', norm(recip_lattice(3,:)));
fprintf('\n');

%% Example 6: Transverse scan (perpendicular to c-axis)
fprintf('EXAMPLE 6: Transverse scan in (h,h,0) plane\n');
fprintf('============================================\n');
qh_values = linspace(0, 1, 101)';
qvec_hh0_line = [qh_values, qh_values, zeros(length(qh_values), 1)];
qvec_hh0_cart = qvec_hh0_line * recip_lattice;
fprintf('Generated %d q-points along (h,h,0) direction\n', size(qvec_hh0_line, 1));
fprintf('Range: (0,0,0) to (1,1,0) in reduced coordinates\n\n');

%% Visualization
fprintf('Generating visualization...\n');

figure('Position', [100, 100, 1400, 900]);

% Plot 1: High-symmetry points
subplot(2, 3, 1);
plot3(qvec_sym_cart(:,1), qvec_sym_cart(:,2), qvec_sym_cart(:,3), ...
      'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
grid on; axis equal;
xlabel('q_x (Å^{-1})'); ylabel('q_y (Å^{-1})'); zlabel('q_z (Å^{-1})');
title('High-Symmetry Points');
view(45, 30);

% Plot 2: Path
subplot(2, 3, 2);
plot3(qvec_path_cart(:,1), qvec_path_cart(:,2), qvec_path_cart(:,3), ...
      'b-', 'LineWidth', 2);
hold on;
plot3(qvec_sym_cart(:,1), qvec_sym_cart(:,2), qvec_sym_cart(:,3), ...
      'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on; axis equal;
xlabel('q_x (Å^{-1})'); ylabel('q_y (Å^{-1})'); zlabel('q_z (Å^{-1})');
title('Path Through BZ');
view(45, 30);

% Plot 3: Coarse grid
subplot(2, 3, 3);
scatter3(qvec_grid_cart(:,1), qvec_grid_cart(:,2), qvec_grid_cart(:,3), ...
         20, sqrt(sum(qvec_grid_cart.^2, 2)), 'filled');
colorbar; colormap('jet');
grid on; axis equal;
xlabel('q_x (Å^{-1})'); ylabel('q_y (Å^{-1})'); zlabel('q_z (Å^{-1})');
title('Coarse Grid (11³)');
view(45, 30);

% Plot 4: [001] line
subplot(2, 3, 4);
q_mag_001 = sqrt(sum(qvec_001_cart.^2, 2));
plot(qz_values, q_mag_001, 'b-', 'LineWidth', 2);
grid on;
xlabel('q_z (r.l.u.)'); ylabel('|q| (Å^{-1})');
title('[001] Direction');

% Plot 5: (h,h,0) line
subplot(2, 3, 5);
q_mag_hh0 = sqrt(sum(qvec_hh0_cart.^2, 2));
plot(qh_values, q_mag_hh0, 'r-', 'LineWidth', 2);
grid on;
xlabel('(h,h,0) (r.l.u.)'); ylabel('|q| (Å^{-1})');
title('(h,h,0) Plane');

% Plot 6: First Brillouin Zone boundaries
subplot(2, 3, 6);
% Draw first BZ for tetragonal lattice
a_star = recip_lattice(1,:);
b_star = recip_lattice(2,:);
c_star = recip_lattice(3,:);

% BZ corners
corners = [
     0.5,  0.5,  0.5;
     0.5,  0.5, -0.5;
     0.5, -0.5,  0.5;
     0.5, -0.5, -0.5;
    -0.5,  0.5,  0.5;
    -0.5,  0.5, -0.5;
    -0.5, -0.5,  0.5;
    -0.5, -0.5, -0.5;
];
corners_cart = corners * recip_lattice;

% Draw edges of BZ
edges = [1 2; 1 3; 1 5; 2 4; 2 6; 3 4; 3 7; 4 8; 5 6; 5 7; 6 8; 7 8];
hold on;
for i = 1:size(edges, 1)
    p1 = corners_cart(edges(i,1), :);
    p2 = corners_cart(edges(i,2), :);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', 1.5);
end

% Add high-symmetry points
plot3(qvec_sym_cart(:,1), qvec_sym_cart(:,2), qvec_sym_cart(:,3), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Add some q-points from fine grid
sample_indices = 1:50:size(qvec_fine_cart, 1);
scatter3(qvec_fine_cart(sample_indices,1), ...
         qvec_fine_cart(sample_indices,2), ...
         qvec_fine_cart(sample_indices,3), ...
         15, 'b', 'filled', 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeAlpha', 0.3);

grid on; axis equal;
xlabel('q_x (Å^{-1})'); ylabel('q_y (Å^{-1})'); zlabel('q_z (Å^{-1})');
title('First Brillouin Zone');
view(45, 30);

sgtitle('LiHoF4 Q-Space Sampling', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('Visualization complete!\n\n');

%% Save q-vectors to file (optional)
fprintf('Saving q-vectors to .mat file...\n');
save('LiHoF4_qvectors.mat', 'qvec_sym', 'qvec_sym_cart', ...
     'qvec_path', 'qvec_path_cart', ...
     'qvec_grid', 'qvec_grid_cart', ...
     'qvec_fine', 'qvec_fine_cart', ...
     'qvec_001_line', 'qvec_001_cart', ...
     'qvec_hh0_line', 'qvec_hh0_cart', ...
     'lattice', 'recip_lattice', ...
     'a', 'c');
fprintf('Saved to: LiHoF4_qvectors.mat\n\n');

%% Print usage instructions
fprintf('=== USAGE IN MF_RPA_Yikai.m ===\n');
fprintf('To use these q-vectors in your RPA calculation:\n\n');
fprintf('1. Load the q-vectors:\n');
fprintf('   load(''LiHoF4_qvectors.mat'', ''qvec_fine'');\n\n');
fprintf('2. Replace the qvec definition in MF_RPA_Yikai.m (around line 24) with:\n');
fprintf('   qvec = qvec_fine;  %% Use fine grid for integration\n\n');
fprintf('3. Or use path for dispersion plots:\n');
fprintf('   qvec = qvec_path;  %% Use path for plotting\n\n');
fprintf('4. Or create custom q-vectors on the fly:\n');
fprintf('   [qvec, ~, ~] = compute_qvectors(lattice, ''mode'', ''grid'', ''grid'', [31,31,31]);\n\n');

fprintf('=== INTEGRATION OVER Q-SPACE ===\n');
fprintf('For proper q-space integration in effective_medium.m:\n');
fprintf('- Use a dense uniform grid (e.g., 31³ or 41³ points)\n');
fprintf('- Weight by Brillouin zone volume: dq = (1/N_q)\n');
fprintf('- Example integration:\n');
fprintf('  G_local = (1/N_q) * sum(G_q, ''all'')  %% Average over q-points\n\n');

fprintf('Script complete!\n');
