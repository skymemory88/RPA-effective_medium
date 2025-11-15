function [qvec, qvec_cart, recip_lattice] = qVec_generator(lattice, varargin)
% qVec_generator Computes q-vectors in reciprocal space based on lattice structure
%
% USAGE:
%   [qvec, qvec_cart, recip_lattice] = qVec_generator(lattice)
%   [qvec, qvec_cart, recip_lattice] = qVec_generator(lattice, 'grid', [n1, n2, n3])
%   [qvec, qvec_cart, recip_lattice] = qVec_generator(lattice, 'path', path_points, n_points)
%
% INPUTS:
%   lattice       - Real space lattice vectors [3x3] matrix where each ROW is a lattice vector
%                   lattice(1,:) = a-vector
%                   lattice(2,:) = b-vector
%                   lattice(3,:) = c-vector
%
% OPTIONAL INPUTS (Name-Value pairs):
%   'mode'        - Computation mode: 'grid', 'path', or 'highsym' (default: 'highsym')
%
%   FOR 'grid' MODE:
%   'grid'        - Grid dimensions [n1, n2, n3] for uniform q-space mesh (default: [11, 11, 11])
%   'range'       - Q-space range in reduced coordinates (default: [-1, 1] for each direction)
%
%   FOR 'path' MODE:
%   'path'        - High-symmetry path definition as cell array of point labels
%                   e.g., {'Gamma', 'X', 'M', 'Gamma', 'Z'} or as [N x 3] matrix of points
%   'npoints'     - Number of points along each segment (default: 50)
%
%   FOR 'highsym' MODE:
%   'sympoints'   - Cell array of high symmetry point labels to include
%                   (default: all standard points for the lattice type)
%
% OUTPUTS:
%   qvec          - Q-vectors in reduced coordinates [n_q x 3]
%                   (in units of reciprocal lattice vectors: q = h*a* + k*b* + l*c*)
%   qvec_cart     - Q-vectors in Cartesian coordinates [n_q x 3] (in units of Å^-1)
%   recip_lattice - Reciprocal lattice vectors [3x3] where each ROW is a reciprocal vector
%                   recip_lattice(1,:) = a* vector (2π/a direction)
%                   recip_lattice(2,:) = b* vector
%                   recip_lattice(3,:) = c* vector
%
% EXAMPLES:
%   % For LiHoF4 (tetragonal scheelite structure)
%   a = 5.175;  % Angstrom
%   c = 10.75;  % Angstrom
%   lattice = [a, 0, 0; 0, a, 0; 0, 0, c];
%
%   % Get high-symmetry points
%   [qvec, qvec_cart, recip] = qVec_generator(lattice);
%
%   % Get a uniform grid
%   [qvec, qvec_cart] = qVec_generator(lattice, 'grid', [15, 15, 15]);
%
%   % Get a path through reciprocal space
%   [qvec, qvec_cart] = qVec_generator(lattice, 'mode', 'path', ...
%                                        'path', {'Gamma', 'X', 'M', 'Gamma', 'Z', 'R'}, ...
%                                        'npoints', 100);
%
% LATTICE TYPE DETECTION:
%   The function automatically detects:
%   - Cubic: a = b = c, all angles 90°
%   - Tetragonal: a = b ≠ c, all angles 90°  (LiHoF4)
%   - Orthorhombic: a ≠ b ≠ c, all angles 90°
%   - General: other cases
%
% HIGH-SYMMETRY POINTS FOR TETRAGONAL (LiHoF4):
%   Gamma = (0, 0, 0)      - Zone center
%   X     = (0.5, 0, 0)    - Zone boundary along a*
%   M     = (0.5, 0.5, 0)  - Zone boundary corner in ab* plane
%   Z     = (0, 0, 0.5)    - Zone boundary along c*
%   R     = (0.5, 0.5, 0.5)- Zone boundary corner
%   A     = (0.5, 0, 0.5)  - Zone boundary
%
% NOTES:
%   - Q-vectors are given in "reduced" coordinates (h,k,l) where
%     q_cart = h*a* + k*b* + l*c* with a* = 2π * (b×c)/(a·(b×c)), etc.
%   - For RPA calculations, these q-vectors can be directly used in MF_RPA_Yikai.m
%   - The function handles both primitive and conventional cells

%% Parse inputs
p = inputParser;
addRequired(p, 'lattice', @(x) isnumeric(x) && isequal(size(x), [3, 3]));
addParameter(p, 'mode', 'highsym', @(x) ischar(x) && ismember(x, {'grid', 'path', 'highsym'}));
addParameter(p, 'grid', [11, 11, 11], @(x) isnumeric(x) && length(x) == 3);
addParameter(p, 'range', [-1, 1], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'path', {}, @(x) iscell(x) || (isnumeric(x) && size(x,2) == 3));
addParameter(p, 'npoints', 50, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sympoints', {}, @iscell);
parse(p, lattice, varargin{:});

mode = p.Results.mode;
grid_size = p.Results.grid;
q_range = p.Results.range;
path_spec = p.Results.path;
n_path_points = p.Results.npoints;
sympoints_spec = p.Results.sympoints;

%% Compute reciprocal lattice vectors
% Real space lattice vectors (rows of lattice matrix)
a = lattice(1,:);
b = lattice(2,:);
c = lattice(3,:);

% Unit cell volume
V = abs(dot(a, cross(b, c)));

% Reciprocal lattice vectors (without 2π factor for reduced coordinates)
% These are defined so that a_i · b*_j = 2π δ_ij
a_star = 2*pi * cross(b, c) / V;
b_star = 2*pi * cross(c, a) / V;
c_star = 2*pi * cross(a, b) / V;

recip_lattice = [a_star; b_star; c_star];

%% Detect lattice type
lattice_type = detect_lattice_type(lattice);
fprintf('Detected lattice type: %s\n', lattice_type);
fprintf('Unit cell volume: %.4f Å³\n', V);
fprintf('Lattice parameters:\n');
a_len = norm(a);
b_len = norm(b);
c_len = norm(c);
fprintf('  a = %.4f Å\n', a_len);
fprintf('  b = %.4f Å\n', b_len);
fprintf('  c = %.4f Å\n', c_len);
fprintf('Reciprocal lattice parameters:\n');
fprintf('  a* = %.4f Å⁻¹ (= 2π/%.4f)\n', norm(a_star), 2*pi/norm(a_star));
fprintf('  b* = %.4f Å⁻¹ (= 2π/%.4f)\n', norm(b_star), 2*pi/norm(b_star));
fprintf('  c* = %.4f Å⁻¹ (= 2π/%.4f)\n', norm(c_star), 2*pi/norm(c_star));
fprintf('\n');

%% Generate q-vectors based on mode
switch mode
    case 'grid'
        [qvec, qvec_cart] = generate_grid(grid_size, q_range, recip_lattice);

    case 'path'
        [qvec, qvec_cart] = generate_path(path_spec, n_path_points, lattice_type, recip_lattice);

    case 'highsym'
        [qvec, qvec_cart] = generate_highsym(lattice_type, sympoints_spec, recip_lattice);
end

fprintf('Generated %d q-vectors\n\n', size(qvec, 1));

%% Display q-vectors
% display_qvectors(qvec, qvec_cart);

end

%% Helper function: Detect lattice type
function lattice_type = detect_lattice_type(lattice)
    a = lattice(1,:);
    b = lattice(2,:);
    c = lattice(3,:);

    a_len = norm(a);
    b_len = norm(b);
    c_len = norm(c);

    tol = 1e-3;  % Tolerance for comparing lengths

    % Check if angles are 90 degrees
    angle_ab = acos(dot(a,b)/(a_len*b_len)) * 180/pi;
    angle_bc = acos(dot(b,c)/(b_len*c_len)) * 180/pi;
    angle_ca = acos(dot(c,a)/(c_len*a_len)) * 180/pi;

    orthogonal = abs(angle_ab - 90) < tol && abs(angle_bc - 90) < tol && abs(angle_ca - 90) < tol;

    if orthogonal
        if abs(a_len - b_len) < tol && abs(b_len - c_len) < tol
            lattice_type = 'cubic';
        elseif abs(a_len - b_len) < tol && abs(b_len - c_len) > tol
            lattice_type = 'tetragonal';
        else
            lattice_type = 'orthorhombic';
        end
    else
        lattice_type = 'general';
    end
end

%% Helper function: Generate uniform grid
function [qvec, qvec_cart] = generate_grid(grid_size, q_range, recip_lattice)
    fprintf('Generating uniform q-space grid...\n');
    fprintf('Grid size: [%d, %d, %d]\n', grid_size(1), grid_size(2), grid_size(3));
    fprintf('Q-range: [%.2f, %.2f]\n', q_range(1), q_range(2));

    % Create grid in reduced coordinates
    qx = linspace(q_range(1), q_range(2), grid_size(1));
    qy = linspace(q_range(1), q_range(2), grid_size(2));
    qz = linspace(q_range(1), q_range(2), grid_size(3));

    [QX, QY, QZ] = meshgrid(qx, qy, qz);

    qvec = [QX(:), QY(:), QZ(:)];

    % Convert to Cartesian
    qvec_cart = qvec * recip_lattice;
end

%% Helper function: Generate path through high-symmetry points
function [qvec, qvec_cart] = generate_path(path_spec, n_points, lattice_type, recip_lattice)
    fprintf('Generating q-space path through high-symmetry points...\n');

    % Get high-symmetry points for this lattice type
    sym_points = get_symmetry_points(lattice_type);

    % Parse path specification
    if iscell(path_spec)
        % Convert labels to coordinates
        path_coords = zeros(length(path_spec), 3);
        for i = 1:length(path_spec)
            label = path_spec{i};
            if isfield(sym_points, label)
                path_coords(i,:) = sym_points.(label);
            else
                error('Unknown symmetry point: %s', label);
            end
        end
    else
        % Already in coordinate form
        path_coords = path_spec;
    end

    % Generate points along path
    n_segments = size(path_coords, 1) - 1;
    qvec = zeros(n_segments * n_points, 3);

    for i = 1:n_segments
        start_pt = path_coords(i,:);
        end_pt = path_coords(i+1,:);

        t = linspace(0, 1, n_points)';
        segment = start_pt + t * (end_pt - start_pt);

        idx_start = (i-1) * n_points + 1;
        idx_end = i * n_points;
        qvec(idx_start:idx_end, :) = segment;
    end

    % Convert to Cartesian
    qvec_cart = qvec * recip_lattice;

    fprintf('Path segments: %d\n', n_segments);
    fprintf('Points per segment: %d\n', n_points);
end

%% Helper function: Generate high-symmetry points only
function [qvec, qvec_cart] = generate_highsym(lattice_type, sympoints_spec, recip_lattice)
    fprintf('Generating high-symmetry q-points for %s lattice...\n', lattice_type);

    % Get high-symmetry points
    sym_points = get_symmetry_points(lattice_type);
    point_names = fieldnames(sym_points);

    % Filter if specific points requested
    if ~isempty(sympoints_spec)
        point_names = intersect(point_names, sympoints_spec, 'stable');
    end

    % Extract coordinates
    n_points = length(point_names);
    qvec = zeros(n_points, 3);

    for i = 1:n_points
        qvec(i,:) = sym_points.(point_names{i});
    end

    % Convert to Cartesian
    qvec_cart = qvec * recip_lattice;

    fprintf('High-symmetry points included:\n');
    for i = 1:n_points
        fprintf('  %s = (%.3f, %.3f, %.3f)\n', point_names{i}, qvec(i,1), qvec(i,2), qvec(i,3));
    end
end

%% Helper function: Get high-symmetry points for lattice type
function sym_points = get_symmetry_points(lattice_type)
    switch lattice_type
        case 'tetragonal'
            % High-symmetry points for tetragonal lattice (e.g., LiHoF4)
            sym_points.Gamma = [0.0, 0.0, 0.0];  % Zone center
            sym_points.X     = [0.5, 0.0, 0.0];  % Zone boundary along a*
            sym_points.Y     = [0.0, 0.5, 0.0];  % Zone boundary along b* (equiv to X)
            sym_points.M     = [0.5, 0.5, 0.0];  % Zone boundary corner in ab* plane
            sym_points.Z     = [0.0, 0.0, 0.5];  % Zone boundary along c*
            sym_points.A     = [0.5, 0.0, 0.5];  % Zone boundary
            sym_points.R     = [0.5, 0.5, 0.5];  % Zone boundary corner

        case 'cubic'
            % High-symmetry points for cubic lattice
            sym_points.Gamma = [0.0, 0.0, 0.0];  % Zone center
            sym_points.X     = [0.5, 0.0, 0.0];  % Zone boundary
            sym_points.M     = [0.5, 0.5, 0.0];  % Zone boundary corner (face)
            sym_points.R     = [0.5, 0.5, 0.5];  % Zone boundary corner

        case 'orthorhombic'
            % High-symmetry points for orthorhombic lattice
            sym_points.Gamma = [0.0, 0.0, 0.0];
            sym_points.X     = [0.5, 0.0, 0.0];
            sym_points.Y     = [0.0, 0.5, 0.0];
            sym_points.Z     = [0.0, 0.0, 0.5];
            sym_points.U     = [0.5, 0.5, 0.0];
            sym_points.T     = [0.5, 0.0, 0.5];
            sym_points.S     = [0.0, 0.5, 0.5];
            sym_points.R     = [0.5, 0.5, 0.5];

        otherwise
            % General case - just include zone center and boundary points
            sym_points.Gamma = [0.0, 0.0, 0.0];
            sym_points.X     = [0.5, 0.0, 0.0];
            sym_points.Y     = [0.0, 0.5, 0.0];
            sym_points.Z     = [0.0, 0.0, 0.5];
    end
end

%% Helper function: Display q-vectors
function display_qvectors(qvec, qvec_cart)
    fprintf('Q-vectors (showing first 20 and last 5):\n');
    fprintf('%-4s  %-30s  %-30s  %-10s\n', 'No.', 'Reduced (h, k, l)', 'Cartesian (Å⁻¹)', '|q| (Å⁻¹)');
    fprintf('%s\n', repmat('-', 1, 80));

    n_show = min(20, size(qvec, 1));
    for i = 1:n_show
        q_mag = norm(qvec_cart(i,:));
        fprintf('%-4d  (%7.4f, %7.4f, %7.4f)  (%7.4f, %7.4f, %7.4f)  %7.4f\n', ...
                i, qvec(i,1), qvec(i,2), qvec(i,3), ...
                qvec_cart(i,1), qvec_cart(i,2), qvec_cart(i,3), q_mag);
    end

    if size(qvec, 1) > 25
        fprintf('  ...\n');
        for i = size(qvec,1)-4:size(qvec,1)
            q_mag = norm(qvec_cart(i,:));
            fprintf('%-4d  (%7.4f, %7.4f, %7.4f)  (%7.4f, %7.4f, %7.4f)  %7.4f\n', ...
                    i, qvec(i,1), qvec(i,2), qvec(i,3), ...
                    qvec_cart(i,1), qvec_cart(i,2), qvec_cart(i,3), q_mag);
        end
    end
    fprintf('\n');
end
