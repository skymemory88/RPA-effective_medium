function [field, poles, sigma_map, modes, qvec] = RPA_lineplot(mion, dscrt_var, omega_grid, theta, phi, gama, hyp, varargin)
% RPA_LINEPLOT Find RPA poles via SVD with optional q-vector scan.
%
% Usage
%   [field, poles, sigma_map, modes] = RPA_lineplot(mion, T, omega_grid, theta, phi, gamma, hyp)
%   [field, poles, sigma_map, modes, qvec] = RPA_lineplot(..., 'qvec', qpts, 'method', 'page')
%
% Examples:
%   % Single q-point (default q=0)
%   [B, p, s, m] = RPA_lineplot('Er', 0.300, linspace(0,35,3001), 0, 0, 1e-3, 1.0);
%
%   % Multiple q-points with page-wise optimization
%   qpts = [[0 0 0]; [0.1 0 0]; [0.2 0 0]; [0.3 0 0]];
%   [B, p, s, m, q] = RPA_lineplot('Er', 0.3, omega, 0, 0, 1e-3, 1.0, 'qvec', qpts, 'method', 'page');
%
% Required Inputs:
%   mion        Element symbol, e.g. 'Er', 'Ho'.
%   dscrt_var   Temperature (K) for field scan file name.
%   omega_grid  Frequency grid (GHz) to scan for poles.
%   theta       Angle from c-axis in ac-plane (deg).
%   phi         Rotation around c-axis (deg).
%   gama        Homogeneous linewidth γ (GHz).
%   hyp         Hyperfine isotope proportion; 1 -> 'Hz_I=1', else 'Hz_I=0'.
%
% Optional Name-Value Pairs:
%   'qvec'      [N x 3] array of q-vectors in r.l.u. (default: [0 0 0])
%   'method'    Computation method: 'vec' | 'page' | 'sequential' (default: 'vec')
%               'vec'        - Vectorized, compatible with all MATLAB versions
%               'page'       - Page-wise operations, requires MATLAB R2020b+
%               'sequential' - No parallelization, for debugging
%
% Outputs:
%   field       Field magnitudes (T) from loaded file [1 x num_fields]
%   poles       [n_modes x num_fields x num_q] RPA pole frequencies (NaN padded)
%   sigma_map   [num_omega x num_fields x num_q] σ_min(I − χ0·Jq) amplitudes
%   modes       Struct with axis-resolved modes [n_modes x num_fields x num_q]:
%               .xx, .yy, .zz - peak positions of Im χ^RPA_{aa}(ω)
%   qvec        [num_q x 3] - q-vectors used
%
% Notes:
% - Requires eigen-state MAT from LiReF4_MF_Yikai matching (mion,T,theta,phi,hyp)
% - Poles found as local minima of σ_min(I − χ0·Jq) via findpeaks on −σ_min
% - For multi-q scans, 'vec' or 'page' methods provide 10-100x speedup

% Parse optional inputs
p = inputParser;
addParameter(p, 'qvec', [0 0 0], @(x) isnumeric(x) && size(x,2)==3);
addParameter(p, 'method', 'vec', @(x) ismember(x, {'vec', 'page', 'sequential'}));
parse(p, varargin{:});

qvec = p.Results.qvec;
method = p.Results.method;

fprintf('\n========================================\n');
fprintf('RPA SUSCEPTIBILITY POLE FINDER\n');
fprintf('Method: %s | Q-points: %d\n', upper(method), size(qvec,1));
fprintf('Using σ_min(I - χ0·Jq) pole criterion\n');
fprintf('========================================\n\n');

% Options
options.min_prom = 1e-2;
options.min_sep = 5;
options.method = method;

% Locate eigen-state data
if hyp > 0
    nZee_path = 'Hz_I=1';
else
    nZee_path = 'Hz_I=0';
end

if ispc
    eigen_dir = ['C:\Users\skyme\OneDrive - Nexus365\Postdoc\Research projects\Li', mion, ...
        'F4 project\Data\Simulations\mean field\eigen-states\', nZee_path, '\'];
else
    eigen_dir = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Postdoc/Research projects/Li', ...
        mion, 'F4 project/Data/Simulations/mean field/eigen-states/', nZee_path, '/'];
end

filename = strcat(['Hscan_Li',mion,'F4_'],...
    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var, theta, phi, hyp),'.mat');
file_load = fullfile(eigen_dir, filename);

fprintf('Loading: %s\n', filename);
load(file_load, '-mat', 'eee', 'vvv', 'fff', 'ttt', 'ion');
fprintf('Eigen-state data loaded\n');

eigenE = eee;
eigenW = vvv;
fields = fff;
temperatures = ttt;
field = vecnorm(fields, 2, 1);

% Setup constants and compute Jq for all q-vectors
const = setup_constants(ion);
Jq_all = compute_Jq_array(qvec, ion, const);

fprintf('\nNumber of q-points: %d\n', size(qvec, 1));
if size(qvec, 1) == 1
    fprintf('Jq(q=[%.2f %.2f %.2f]) [meV]:\n', qvec(1), qvec(2), qvec(3));
    disp(Jq_all(:,:,1));
end
fprintf('Field: %.3f - %.3f T (%d points)\n', min(field), max(field), length(field));
fprintf('Frequency: %.3f - %.3f GHz (%d pts)\n', min(omega_grid), max(omega_grid), numel(omega_grid));
tempK_print = temperatures; if numel(tempK_print) > 1, tempK_print = tempK_print(1); end
fprintf('Temperature: %.3f mK\n\n', tempK_print*1000);

% Store in options
options.eigenE = eigenE;
options.eigenW = eigenW;
options.temperatures = temperatures;
options.ion = ion;
options.const = const;
options.gamma = gama;

% Find poles using selected method
fprintf('Searching for RPA poles...\n');
[poles, sigma_map, modes] = RPA_poles_unified(Jq_all, qvec, omega_grid, field, options);

% Print statistics
num_q = size(qvec, 1);
for nq = 1:min(num_q, 10)  % Limit output for many q-points
    counts = zeros(1, numel(field));
    for j = 1:numel(field)
        counts(j) = sum(~isnan(poles(:, j, nq)));
    end
    fprintf('Q=[%.2f %.2f %.2f]: pole counts min=%d, median=%d, max=%d\n', ...
        qvec(nq,1), qvec(nq,2), qvec(nq,3), min(counts), median(counts), max(counts));
end
if num_q > 10
    fprintf('... (%d more q-points)\n', num_q - 10);
end

% Plot (only for single q-point or first q-point)
if size(qvec, 1) == 1
    plot_opts.field_label = 'Magnetic Field';
    plot_opts.field_units = 'T';
    plot_opts.omega_label = 'Excitation Energy';
    plot_opts.omega_units = 'GHz';
    plot_opts.title = sprintf('RPA Modes: Li%sF4 at %.2f K, Q=[%.2f %.2f %.2f]', ...
        mion, dscrt_var, qvec(1), qvec(2), qvec(3));

    plot_rpa_modes(field, poles(:,:,1), plot_opts);

    % σ_min color map
    try
        figure('Name','SigmaMin Map','Position',[120,120,900,650]);
        imagesc(field, omega_grid, sigma_map(:,:,1));
        set(gca, 'YDir','normal'); colorbar;
        xlabel('Magnetic Field (T)'); ylabel('Frequency (GHz)');
        title(sprintf('σ_{min}(I - χ0·Jq), Q=[%.2f %.2f %.2f]', qvec(1), qvec(2), qvec(3)));
    catch
    end
end
end


function const = setup_constants(ion)
const.hbar = 1.05457E-34;
const.J2meV = 6.24151e+21;
const.Gh2mV = const.hbar * 2*pi * 1e9 * const.J2meV;
const.muB = 9.27401e-24;
const.muN = 5.05078e-27;
const.elem = find(ion.prop);
const.mu0 = 4*pi*1e-7;
const.dpRng = 100;
const.ELEf = ion.gLande(const.elem) * const.muB * const.J2meV;
const.NUCf = ion.nLande(const.elem) * const.muN * const.J2meV;
const.gfac = ion.gLande(const.elem)^2 * (const.muB)^2 * (const.mu0/4/pi) * const.J2meV * 10^30;
const.kB = 8.61733e-2; % [meV/K]
end


function Jq_all = compute_Jq_array(qvec, ion, const)
% Compute interaction matrix Jq for multiple q-vectors
unitN = 4;
lattice = ion.abc{const.elem};
Vc = sum(lattice(1,:) .* cross(lattice(2,:), lattice(3,:)));

eins = diag([1 1 1]);
eins = repmat(eins, 1, 1, 4, 4);
demagn_t = ellipsoid_demagn(ion.alpha);
demagn = repmat(demagn_t, 1, 1, 4, 4);

num_q = size(qvec, 1);
Jq_all = zeros(3, 3, num_q);

for nq = 1:num_q
    q = qvec(nq, :);

    % Lorentz term applies only for q including all sublattices coherently
    lorz_on = 1;
    if abs(real(sum(exp(1i*2*pi*q*ion.tau'))/size(ion.tau,1))-1) > 1e-10
        lorz_on = 0;
    end

    % Dipole + exchange interaction
    D = const.gfac * (MF_dipole(q, const.dpRng, lattice, ion.tau) - ...
        lorz_on*4*pi/Vc*(eins/3 - ion.demag*demagn)) + ...
        exchange(q, ion.ex(const.elem), lattice, ion.tau);

    % Sublattice average and renormalization
    Jav = squeeze(sum(sum(D(:,:,:,:), 4), 3) / unitN);
    Jq_all(:,:,nq) = -diag(ion.renorm(const.elem,:)) .* Jav;
end
end


function [poles, sigma_map, modes] = RPA_poles_unified(Jq_all, qvec, omega_scan, field, options)
% Unified RPA pole finder - dispatches to appropriate method

switch options.method
    case 'vec'
        [poles, sigma_map, modes] = RPA_poles_vec(Jq_all, qvec, omega_scan, field, options);
    case 'page'
        [poles, sigma_map, modes] = RPA_poles_page(Jq_all, qvec, omega_scan, field, options);
    case 'sequential'
        [poles, sigma_map, modes] = RPA_poles_sequential(Jq_all, qvec, omega_scan, field, options);
    otherwise
        error('Unknown method: %s. Use ''vec'', ''page'', or ''sequential''.', options.method);
end
end


function [poles, sigma_map, modes] = RPA_poles_vec(Jq_all, qvec, omega_scan, field, options)
% Vectorized RPA pole finder - parallelizes over q-vectors

num_q = size(qvec, 1);
num_fields = size(field, 2);
num_omega = numel(omega_scan);

poles = [];
sigma_map = zeros(num_omega, num_fields, num_q);
mx_all = cell(num_fields, num_q);
my_all = cell(num_fields, num_q);
mz_all = cell(num_fields, num_q);

tempK = options.temperatures;
if numel(tempK) > 1, tempK = tempK(1); end
beta = 1 / max(options.const.kB * max(tempK, 0), eps);

JhT = spin_ops(options.ion, options.const);
gamma_meV = options.gamma * options.const.Gh2mV;
I3 = eye(3);

% Parallelize over q-vectors
parfor nq = 1:num_q
    Jq = Jq_all(:,:,nq);
    p_list = cell(1, num_fields);
    im_diag_local = zeros(num_omega, num_fields, 3);
    sigma_local_q = zeros(num_omega, num_fields);

    for ii = 1:num_fields
        en = squeeze(options.eigenE(ii, :))';
        v = squeeze(options.eigenW(ii, :, :));
        if isempty(en) || isempty(v), continue; end

        % Thermal populations
        if tempK > 0
            zn = exp(-beta * (en - min(en)));
            zn = zn / max(sum(zn), eps);
        else
            zn = zeros(size(en)); zn(1) = 1;
        end
        [n, np] = meshgrid(zn, zn);
        NN = n - np;

        % Transition energies
        [ee, eep] = meshgrid(en, en);
        omega_nm = eep - ee;

        % Precompute J_nm matrices
        Jnm = cell(3,1);
        J_ops = {JhT.x, JhT.y, JhT.z};
        for a = 1:3
            Jnm{a} = v' * J_ops{a} * v;
        end

        % Build stacked matrix
        N = size(omega_nm, 1);
        Mstack = zeros(9, N*N);
        row = 1;
        for b = 1:3
            for a = 1:3
                M_ab = Jnm{a} .* conj(Jnm{b}) .* NN;
                Mstack(row, :) = reshape(M_ab, 1, []);
                row = row + 1;
            end
        end

        sigma_local = zeros(num_omega, 1);
        im_diag_field = zeros(num_omega, 3);

        % Process frequency in chunks
        W_all = omega_scan(:).' * options.const.Gh2mV;
        max_elems = 5e6;
        chunk = max(1, floor(max_elems / (N*N)));

        for s = 1:chunk:num_omega
            e = min(s + chunk - 1, num_omega);
            W = W_all(s:e);
            DenInv = 1 ./ ((omega_nm(:) - 1i*gamma_meV) - W);
            chi0_flat = Mstack * DenInv;

            for jj = 1:(e - s + 1)
                idx = s + jj - 1;
                chi0 = reshape(chi0_flat(:, jj), [3, 3]);
                A = I3 - chi0 * Jq;
                svals = svd(A);
                sigma_local(idx) = min(svals);
                chi_rpa = A \ chi0;
                d = imag(diag(chi_rpa));
                im_diag_field(idx, :) = d(:);
            end
        end

        sigma_local_q(:, ii) = sigma_local;
        im_diag_local(:, ii, :) = im_diag_field;

        % Find poles
        [~, locs] = findpeaks(-sigma_local, ...
            'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        if ~isempty(locs)
            locs = sort(locs, 'ascend');
            p_list{ii} = omega_scan(locs(:));
        else
            p_list{ii} = [];
        end
    end

    sigma_map(:, :, nq) = sigma_local_q;

    % Convert poles to matrix
    max_modes = 0;
    for j = 1:num_fields
        max_modes = max(max_modes, numel(p_list{j}));
    end
    poles_q = nan(max_modes, num_fields);
    for j = 1:num_fields
        pj = p_list{j};
        if ~isempty(pj)
            poles_q(1:numel(pj), j) = pj(:);
        end
    end

    if nq == 1
        poles = nan(max_modes, num_fields, num_q);
    end
    poles(:, :, nq) = poles_q;

    % Extract axis-resolved modes
    for ii = 1:num_fields
        sigx = im_diag_local(:, ii, 1);
        sigy = im_diag_local(:, ii, 2);
        sigz = im_diag_local(:, ii, 3);
        [~, lx] = findpeaks(sigx, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        [~, ly] = findpeaks(sigy, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        [~, lz] = findpeaks(sigz, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        mx_all{ii, nq} = omega_scan(lx(:));
        my_all{ii, nq} = omega_scan(ly(:));
        mz_all{ii, nq} = omega_scan(lz(:));
    end
end

modes.xx = pad_cell_3d(mx_all);
modes.yy = pad_cell_3d(my_all);
modes.zz = pad_cell_3d(mz_all);
end


function [poles, sigma_map, modes] = RPA_poles_page(Jq_all, qvec, omega_scan, field, options)
% Page-wise RPA pole finder - requires MATLAB R2020b+

num_q = size(qvec, 1);
num_fields = size(field, 2);
num_omega = numel(omega_scan);

poles = [];
sigma_map = zeros(num_omega, num_fields, num_q);
mx_all = cell(num_fields, num_q);
my_all = cell(num_fields, num_q);
mz_all = cell(num_fields, num_q);

tempK = options.temperatures;
if numel(tempK) > 1, tempK = tempK(1); end
beta = 1 / max(options.const.kB * max(tempK, 0), eps);

JhT = spin_ops(options.ion, options.const);
gamma_meV = options.gamma * options.const.Gh2mV;

% Parallelize over q-vectors
parfor nq = 1:num_q
    Jq = Jq_all(:,:,nq);
    p_list = cell(1, num_fields);
    im_diag_local = zeros(num_omega, num_fields, 3);
    sigma_local_q = zeros(num_omega, num_fields);

    for ii = 1:num_fields
        en = squeeze(options.eigenE(ii, :))';
        v = squeeze(options.eigenW(ii, :, :));
        if isempty(en) || isempty(v), continue; end

        % Thermal populations
        if tempK > 0
            zn = exp(-beta * (en - min(en)));
            zn = zn / max(sum(zn), eps);
        else
            zn = zeros(size(en)); zn(1) = 1;
        end
        [n, np] = meshgrid(zn, zn);
        NN = n - np;

        % Transition energies
        [ee, eep] = meshgrid(en, en);
        omega_nm = eep - ee;

        % Precompute J_nm matrices
        Jnm = cell(3,1);
        J_ops = {JhT.x, JhT.y, JhT.z};
        for a = 1:3
            Jnm{a} = v' * J_ops{a} * v;
        end

        % Build stacked matrix
        N = size(omega_nm, 1);
        Mstack = zeros(9, N*N);
        row = 1;
        for b = 1:3
            for a = 1:3
                M_ab = Jnm{a} .* conj(Jnm{b}) .* NN;
                Mstack(row, :) = reshape(M_ab, 1, []);
                row = row + 1;
            end
        end

        sigma_local = zeros(num_omega, 1);
        im_diag_field = zeros(num_omega, 3);

        % Process frequency in chunks with page-wise operations
        W_all = omega_scan(:).' * options.const.Gh2mV;
        max_elems = 5e6;
        chunk = max(1, floor(max_elems / (N*N)));

        for s = 1:chunk:num_omega
            e = min(s + chunk - 1, num_omega);
            W = W_all(s:e);
            n_chunk = e - s + 1;

            DenInv = 1 ./ ((omega_nm(:) - 1i*gamma_meV) - W);
            chi0_flat = Mstack * DenInv;

            % Reshape to pages: [3, 3, n_chunk]
            chi0_pages = reshape(chi0_flat, 3, 3, n_chunk);

            % Page-wise operations
            Jq_expanded = repmat(Jq, 1, 1, n_chunk);
            MM_pages = pagemtimes(chi0_pages, Jq_expanded);
            I_pages = repmat(eye(3), 1, 1, n_chunk);
            A_pages = I_pages - MM_pages;
            chi_rpa_pages = pagemldivide(A_pages, chi0_pages);

            % Extract diagonal imaginary parts (vectorized)
            for a = 1:3
                im_diag_field(s:e, a) = squeeze(imag(chi_rpa_pages(a, a, :)));
            end

            % Compute singular values
            for jj = 1:n_chunk
                idx = s + jj - 1;
                A = A_pages(:, :, jj);
                svals = svd(A);
                sigma_local(idx) = min(svals);
            end
        end

        sigma_local_q(:, ii) = sigma_local;
        im_diag_local(:, ii, :) = im_diag_field;

        % Find poles
        [~, locs] = findpeaks(-sigma_local, ...
            'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        if ~isempty(locs)
            locs = sort(locs, 'ascend');
            p_list{ii} = omega_scan(locs(:));
        else
            p_list{ii} = [];
        end
    end

    sigma_map(:, :, nq) = sigma_local_q;

    % Convert poles to matrix
    max_modes = 0;
    for j = 1:num_fields
        max_modes = max(max_modes, numel(p_list{j}));
    end
    poles_q = nan(max_modes, num_fields);
    for j = 1:num_fields
        pj = p_list{j};
        if ~isempty(pj)
            poles_q(1:numel(pj), j) = pj(:);
        end
    end

    if nq == 1
        poles = nan(max_modes, num_fields, num_q);
    end
    poles(:, :, nq) = poles_q;

    % Extract axis-resolved modes
    for ii = 1:num_fields
        sigx = im_diag_local(:, ii, 1);
        sigy = im_diag_local(:, ii, 2);
        sigz = im_diag_local(:, ii, 3);
        [~, lx] = findpeaks(sigx, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        [~, ly] = findpeaks(sigy, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        [~, lz] = findpeaks(sigz, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        mx_all{ii, nq} = omega_scan(lx(:));
        my_all{ii, nq} = omega_scan(ly(:));
        mz_all{ii, nq} = omega_scan(lz(:));
    end
end

modes.xx = pad_cell_3d(mx_all);
modes.yy = pad_cell_3d(my_all);
modes.zz = pad_cell_3d(mz_all);
end


function [poles, sigma_map, modes] = RPA_poles_sequential(Jq_all, qvec, omega_scan, field, options)
% Sequential RPA pole finder - no parallelization (for debugging/comparison)

num_q = size(qvec, 1);
num_fields = size(field, 2);
num_omega = numel(omega_scan);

poles = [];
sigma_map = zeros(num_omega, num_fields, num_q);
mx_all = cell(num_fields, num_q);
my_all = cell(num_fields, num_q);
mz_all = cell(num_fields, num_q);

tempK = options.temperatures;
if numel(tempK) > 1, tempK = tempK(1); end
beta = 1 / max(options.const.kB * max(tempK, 0), eps);

JhT = spin_ops(options.ion, options.const);
gamma_meV = options.gamma * options.const.Gh2mV;
I3 = eye(3);

% Sequential loops
for nq = 1:num_q
    Jq = Jq_all(:,:,nq);
    p_list = cell(1, num_fields);
    im_diag = zeros(num_omega, num_fields, 3);

    for ii = 1:num_fields
        en = squeeze(options.eigenE(ii, :))';
        v = squeeze(options.eigenW(ii, :, :));
        if isempty(en) || isempty(v), continue; end

        % Thermal populations
        if tempK > 0
            zn = exp(-beta * (en - min(en)));
            zn = zn / max(sum(zn), eps);
        else
            zn = zeros(size(en)); zn(1) = 1;
        end
        [n, np] = meshgrid(zn, zn);
        NN = n - np;

        % Transition energies
        [ee, eep] = meshgrid(en, en);
        omega_nm = eep - ee;

        % Precompute J_nm matrices
        Jnm = cell(3,1);
        J_ops = {JhT.x, JhT.y, JhT.z};
        for a = 1:3
            Jnm{a} = v' * J_ops{a} * v;
        end

        % Build stacked matrix
        N = size(omega_nm, 1);
        Mstack = zeros(9, N*N);
        row = 1;
        for b = 1:3
            for a = 1:3
                M_ab = Jnm{a} .* conj(Jnm{b}) .* NN;
                Mstack(row, :) = reshape(M_ab, 1, []);
                row = row + 1;
            end
        end

        sigma_local = zeros(num_omega, 1);
        im_diag_local = zeros(num_omega, 3);

        % Process frequency in chunks
        W_all = omega_scan(:).' * options.const.Gh2mV;
        max_elems = 5e6;
        chunk = max(1, floor(max_elems / (N*N)));

        for s = 1:chunk:num_omega
            e = min(s + chunk - 1, num_omega);
            W = W_all(s:e);
            DenInv = 1 ./ ((omega_nm(:) - 1i*gamma_meV) - W);
            chi0_flat = Mstack * DenInv;

            for jj = 1:(e - s + 1)
                idx = s + jj - 1;
                chi0 = reshape(chi0_flat(:, jj), [3, 3]);
                A = I3 - chi0 * Jq;
                svals = svd(A);
                sigma_local(idx) = min(svals);
                chi_rpa = A \ chi0;
                d = imag(diag(chi_rpa));
                im_diag_local(idx, :) = d(:);
            end
        end

        sigma_map(:, ii, nq) = sigma_local;
        im_diag(:, ii, :) = im_diag_local;

        % Find poles
        [~, locs] = findpeaks(-sigma_local, ...
            'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        if ~isempty(locs)
            locs = sort(locs, 'ascend');
            p_list{ii} = omega_scan(locs(:));
        else
            p_list{ii} = [];
        end
    end

    % Convert to matrix for this q
    max_modes = 0;
    for j = 1:num_fields
        max_modes = max(max_modes, numel(p_list{j}));
    end
    poles_q = nan(max_modes, num_fields);
    for j = 1:num_fields
        pj = p_list{j};
        if ~isempty(pj)
            poles_q(1:numel(pj), j) = pj(:);
        end
    end

    if nq == 1
        poles = nan(max_modes, num_fields, num_q);
    end
    poles(:, :, nq) = poles_q;

    % Extract axis-resolved modes
    for ii = 1:num_fields
        sigx = im_diag(:, ii, 1);
        sigy = im_diag(:, ii, 2);
        sigz = im_diag(:, ii, 3);
        [~, lx] = findpeaks(sigx, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        [~, ly] = findpeaks(sigy, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        [~, lz] = findpeaks(sigz, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
        mx_all{ii, nq} = omega_scan(lx(:));
        my_all{ii, nq} = omega_scan(ly(:));
        mz_all{ii, nq} = omega_scan(lz(:));
    end
end

modes.xx = pad_cell_3d(mx_all);
modes.yy = pad_cell_3d(my_all);
modes.zz = pad_cell_3d(mz_all);
end


function JhT = spin_ops(ion, const)
ionJ = ion.J(const.elem);
Jz = diag(ionJ:-1:-ionJ);
Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ)).*(ionJ+1+((ionJ-1):-1:-ionJ))), 1);
Jm = Jp';

if ion.hyp(const.elem) > 0
    ionI = ion.I(const.elem);
    Iz = diag(ionI:-1:-ionI);
    IhT.z = kron(eye(2*ionJ+1), Iz);
    Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))), 1);
    Im = Ip';
    Iph = kron(eye(2*ionJ+1), Ip);
    Imh = kron(eye(2*ionJ+1), Im);
    IhT.x = (Iph + Imh)/2;
    IhT.y = (Iph - Imh)/2i;

    Jph = kron(Jp, eye(2*ionI+1));
    Jmh = kron(Jm, eye(2*ionI+1));
    JhT.x = (Jph + Jmh)/2;
    JhT.y = (Jph - Jmh)/2i;
    JhT.z = kron(Jz, eye(2*ionI+1));
else
    JhT.x = (Jp + Jm)/2;
    JhT.y = (Jp - Jm)/2i;
    JhT.z = Jz;
    IhT.x = 0; IhT.y = 0; IhT.z = 0;
end

JhT.x = JhT.x + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.x;
JhT.y = JhT.y + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.y;
JhT.z = JhT.z + ion.hyp(const.elem) * const.NUCf/const.ELEf * IhT.z;
end


function M = pad_cell_3d(C)
[nf, nq] = size(C);
k = 0;
for jj = 1:nf
    for qq = 1:nq
        k = max(k, numel(C{jj, qq}));
    end
end
M = nan(k, nf, nq);
for jj = 1:nf
    for qq = 1:nq
        v = C{jj, qq};
        if ~isempty(v)
            M(1:numel(v), jj, qq) = v(:);
        end
    end
end
end


function fig = plot_rpa_modes(field, pole_frequencies, plot_options)
if nargin < 3, plot_options = struct(); end
if ~isfield(plot_options, 'field_label'), plot_options.field_label = 'Magnetic Field'; end
if ~isfield(plot_options, 'field_units'), plot_options.field_units = 'T'; end
if ~isfield(plot_options, 'omega_label'), plot_options.omega_label = 'Mode Frequency'; end
if ~isfield(plot_options, 'omega_units'), plot_options.omega_units = 'GHz'; end
if ~isfield(plot_options, 'title'), plot_options.title = 'RPA Modes'; end

field = field(:);

% Accept either a numeric matrix [n_modes x n_fields] or a cell array per field
if iscell(pole_frequencies)
    n_fields = numel(pole_frequencies);
    max_modes = 0;
    for j = 1:n_fields
        max_modes = max(max_modes, numel(pole_frequencies{j}));
    end
    n_modes = max_modes;
    M = nan(n_modes, n_fields);
    for j = 1:n_fields
        f = sort(pole_frequencies{j}(:), 'ascend');
        M(1:numel(f), j) = f;
    end
else
    M = pole_frequencies;
    n_modes = size(M, 1);
end

if ~isfield(plot_options, 'mode_labels')
    plot_options.mode_labels = cell(n_modes, 1);
    for i = 1:n_modes
        plot_options.mode_labels{i} = sprintf('Mode %d', i);
    end
else
    if ~iscell(plot_options.mode_labels)
        plot_options.mode_labels = {plot_options.mode_labels};
    end
    n_labels = length(plot_options.mode_labels);
    if n_labels < n_modes
        for i = (n_labels+1):n_modes
            plot_options.mode_labels{i} = sprintf('Mode %d', i);
        end
    end
end

fig = figure('Position', [100, 100, 1000, 700]);
hold on;
colors = viridis(n_modes);

for i = 1:n_modes
    freq_row = M(i, :).';
    valid_mask = ~isnan(freq_row);
    field_valid = field(valid_mask);
    freq_valid = freq_row(valid_mask);

    if isempty(freq_valid), continue; end
    sizes = 50 * ones(size(freq_valid));

    scatter(field_valid, freq_valid, sizes, ...
        'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.5, 'MarkerFaceAlpha', 0.7, ...
        'DisplayName', plot_options.mode_labels{i});

    plot(field_valid, freq_valid, '-', 'Color', [colors(i, :), 0.3], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlabel(sprintf('%s (%s)', plot_options.field_label, plot_options.field_units), ...
    'FontSize', 13, 'FontWeight', 'bold');
ylabel(sprintf('%s (%s)', plot_options.omega_label, plot_options.omega_units), ...
    'FontSize', 13, 'FontWeight', 'bold');
title(plot_options.title, 'FontSize', 15, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'GridAlpha', 0.3, 'FontSize', 11);
box on;
hold off;
end


function cmap = viridis(n)
base_colors = [
    0.267004, 0.004874, 0.329415; 0.282623, 0.140926, 0.457517;
    0.253935, 0.265254, 0.529983; 0.206756, 0.371758, 0.553117;
    0.163625, 0.471133, 0.558148; 0.127568, 0.566949, 0.550556;
    0.134692, 0.658636, 0.517649; 0.266941, 0.748751, 0.440573;
    0.477504, 0.821444, 0.318195; 0.741388, 0.873449, 0.149561;
    0.993248, 0.906157, 0.143936
    ];
x_base = linspace(0, 1, size(base_colors, 1));
x_new = linspace(0, 1, n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(x_base, base_colors(:, i), x_new);
end
end
