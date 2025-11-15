function [field, poles, sigma_map, modes] = RPA_line(mion, dscrt_var, omega_grid, theta, phi, gama, hyp)
% RPA_LINEPLOT Find RPA poles via SVD and plot modes vs field.
%
% Usage
%   [field, poles, sigma_map] = RPA_line(mion, T, omega_grid, theta, phi, gamma, hyp)
%   Example:
%     [B, modes, sig] = RPA_line('Er', 0.300, linspace(0,35,3001), 0, 0, 1e-3, 1.0);
%
% Inputs
%   mion        Element symbol, e.g. 'Er', 'Ho'.
%   dscrt_var   Temperature (K) used in LiReF4_MF_Yikai field scan file name.
%   omega_grid  Frequency grid (GHz) to scan for poles.
%   theta       Angle from c-axis in ac-plane (deg).
%   phi         Rotation around c-axis (deg).
%   gama        Homogeneous linewidth γ (GHz).
%   hyp         Hyperfine isotope proportion; 1 -> use 'Hz_I=1' data, else 'Hz_I=0'.
%
% Outputs
%   field       Field magnitudes (T) from the loaded eigen-state file.
%   poles       [n_modes x num_fields] numeric; each column lists poles at a field (NaN padded).
%   sigma_map   [num_omega x num_fields] σ_min(I − χ0·Jq) amplitudes across ω and field.
%   modes       struct with axis-resolved mode positions (GHz):
%               modes.xx, modes.yy, modes.zz are [n_modes x num_fields] numeric (NaN padded),
%               giving peak positions of Im χ^RPA_{aa}(ω) per field.
%
% Notes
% - Requires eigen-state MAT produced by LiReF4_MF_Yikai with matching (mion,T,theta,phi,hyp).
% - J(q=0) uses sublattice-averaged dipole+Lorentz minus exchange, scaled by diag(ion.renorm).
% - Poles are local minima of σ_min(I − χ0·Jq) along ω found via findpeaks on −σ_min.
%
% Physics (Jensen & Mackintosh, ch. 3)
% - χ0_{αβ}(ω) = Σ_nm (ρ_n−ρ_m)⟨n|J_α|m⟩⟨m|J_β|n⟩/(E_m−E_n−ℏω−iγ)
% - χ^RPA = (I − χ0·Jq)\χ0; poles where det(I − χ0·Jq) = 0
% - Numerically track σ_min(I − χ0·Jq) instead of determinant for stability.

fprintf('\n========================================\n');
fprintf('RPA SUSCEPTIBILITY POLE FINDER\n');
fprintf('Using σ_min(I - χ0·Jq) pole criterion\n');
fprintf('========================================\n\n');

% Options
options.min_prom = 1e-2; % peak prominence on σ_min for robustness
options.min_sep = 5;     % min peak distance [grid points]

% Locate eigen-state data produced by LiReF4_MF_Yikai
if hyp > 0
    nZee_path = 'Hz_I=1';
else
    nZee_path = 'Hz_I=0';
end

if ispc
    eigen_dir = ['C:\Users\skyme\OneDrive - Nexus365\Postdoc\Research projects\Li', mion, 'F4 project\Data\Simulations\mean field\eigen-states\', nZee_path, '\'];
else
    eigen_dir = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Postdoc/Research projects/Li', mion, 'F4 project/Data/Simulations/mean field/eigen-states/', nZee_path, '/'];
end

filename = strcat(['Hscan_Li',mion,'F4_'],...
    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', dscrt_var, theta, phi, hyp),'.mat');
file_load = fullfile(eigen_dir, filename);

% Try to load eigenstate data for on-demand chi0 computation
fprintf('Loading data from: %s\n', filename);

% Load eigenstate data
load(file_load, '-mat', 'eee', 'vvv', 'fff', 'ttt', 'ion');
fprintf('Eigen-state data loaded\n');
eigenE = eee;
eigenW = vvv;
fields = fff;
temperatures = ttt;
field = vecnorm(fields, 2, 1);

% Setup constants and exchange
const = setup_constants(ion);
Jq = interac_q(ion, const);

fprintf('\nExchange matrix Jq [meV]:\n');
disp(Jq);
fprintf('\nField: %.3f - %.3f T (%d points)\n', min(field), max(field), length(field));
fprintf('Frequency grid: %.3f - %.3f GHz (%d pts)\n', min(omega_grid), max(omega_grid), numel(omega_grid));
tempK_print = temperatures; if numel(tempK_print) > 1, tempK_print = tempK_print(1); end
fprintf('Temperature: %.3f mK\n\n', tempK_print*1000);

% Store in options
options.eigenE = eigenE;
options.eigenW = eigenW;
options.temperatures = temperatures;
options.ion = ion;
options.const = const;
options.gamma = gama; % [GHz]

% Find poles via σ_min
fprintf('Searching for RPA poles\n');

[poles, sigma_map, modes] = RPA_poles(Jq, omega_grid, field, options);

counts = zeros(1, numel(field));
for j = 1:numel(field)
    counts(j) = sum(~isnan(poles(:, j))); 
end
fprintf('Per-field pole counts: min=%d, median=%d, max=%d\n', min(counts), median(counts), max(counts));

% Plot
plot_opts.field_label = 'Magnetic Field';
plot_opts.field_units = 'T';
plot_opts.omega_label = 'Excitation Energy';
plot_opts.omega_units = 'GHz';
plot_opts.title = sprintf('RPA Modes: Li%sF4 at %.2f K', mion, dscrt_var);

plot_rpa_modes(field, poles, plot_opts);

% Optional diagnostic: σ_min color map
try
    figure('Name','SigmaMin Map','Position',[120,120,900,650]);
    imagesc(field, omega_grid, sigma_map);
    set(gca, 'YDir','normal'); colorbar;
    xlabel('Magnetic Field (T)'); ylabel('Frequency (GHz)');
    title('σ_{min}(I - χ0·Jq)');
catch
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


function Jq = interac_q(ion, const)
unitN = 4;
lattice = ion.abc{const.elem};
Vc = sum(lattice(1,:) .* cross(lattice(2,:), lattice(3,:)));

eins = diag([1 1 1]);
eins = repmat(eins,1,1,4,4);
demagn_t = ellipsoid_demagn(ion.alpha);
demagn = repmat(demagn_t,1,1,4,4);

% J(q=0) from dipole - Lorentz (demag) + exchange; sublattice-averaged and renormalized
lorz_on = 1;
D = const.gfac * (MF_dipole([0 0 0], const.dpRng, lattice, ion.tau) - ...
    lorz_on*4*pi/Vc*(eins/3 - ion.demag*demagn)) + ...
    exchange([0 0 0], ion.ex(const.elem), lattice, ion.tau);
Jav = squeeze(sum(sum(D(:,:,:,:),4),3)/unitN);
Jq = -diag(ion.renorm(const.elem,:)) .* Jav;
end

function [poles, sigma_map, modes] = RPA_poles(Jq, omega_scan, field, options)
% Compute σ_min(I - χ0·Jq) across ω for each field; pick minima as poles
tempK = options.temperatures; % [K]
if numel(tempK) > 1, tempK = tempK(1); end
beta = 1 / max(options.const.kB * max(tempK, 0), eps); % [meV^-1]; T=0 -> β→∞ handled below

JhT = spin_ops(options.ion, options.const); % electron+nuclear spin operators (renormalized)

num_fields = size(field, 2);
num_omega = numel(omega_scan);
sigma_map = zeros(num_omega, num_fields);
p_list = cell(1, num_fields);
im_diag = zeros(num_omega, num_fields, 3); % store Im diag(χ^RPA)

gamma_meV = options.gamma * options.const.Gh2mV; % linewidth in meV

for ii = 1:num_fields
    % Eigen-system for this field
    en = squeeze(options.eigenE(ii, :))'; % [meV]
    v = squeeze(options.eigenW(ii, :, :)); % [N x N]
    if isempty(en) || isempty(v)
        continue
    end
    % Thermal populations ρ_n
    if tempK > 0
        zn = exp(-beta * (en - min(en)));
        zn = zn / max(sum(zn), eps);
    else
        zn = zeros(size(en)); zn(1) = 1;
    end
    [n, np] = meshgrid(zn, zn);
    NN = n - np; % ρ_n - ρ_m

    % Transition energies and matrix elements
    [ee, eep] = meshgrid(en, en);
    omega_nm = eep - ee; % [meV]

    % Precompute J_nm matrices once (3 costly ops instead of 18)
    Jnm = cell(3,1);
    J_ops = {JhT.x, JhT.y, JhT.z};
    for a = 1:3
        Jnm{a} = v' * J_ops{a} * v; % [N x N]
    end

    % Build stacked M rows for all αβ pairs: row k = vec(Jnm{α} .* conj(Jnm{β}) .* NN)
    % Row order chosen so reshape(vec,[3,3]) reconstructs χ0 by column-major (β fastest)
    % k = (β-1)*3 + α
    N = size(omega_nm,1);
    Mstack = zeros(9, N*N);
    row = 1;
    for b = 1:3
        for a = 1:3
            M_ab = Jnm{a} .* conj(Jnm{b}) .* NN;
            Mstack(row, :) = reshape(M_ab, 1, []);
            row = row + 1;
        end
    end

    % Process ω in chunks to limit memory: [N^2 x n_chunk]
    W_all = omega_scan(:).' * options.const.Gh2mV; % [1 x num_omega] meV
    max_elems = 5e6; % ~80 MB for complex doubles
    chunk = max(1, floor(max_elems / (N*N)));
    for s = 1:chunk:num_omega
        e = min(s + chunk - 1, num_omega);
        W = W_all(s:e);
        DenInv = 1 ./ ((omega_nm(:) - 1i*gamma_meV) - W); % [N^2 x n_chunk]
        chi0_flat = Mstack * DenInv; % [9 x n_chunk]

        % Small 3x3 ops per ω in chunk
        for jj = 1:(e - s + 1)
            idx = s + jj - 1;
            chi0 = reshape(chi0_flat(:, jj), [3, 3]);
            A = eye(3) - chi0 * Jq; % RPA denominator
            svals = svd(A);
            sigma_map(idx, ii) = min(svals);
            chi_rpa = A \ chi0;
            d = imag(diag(chi_rpa));
            im_diag(idx, ii, :) = d(:);
        end
    end

    % Identify local minima of σ_min along ω
    [~, locs] = findpeaks(-sigma_map(:, ii), ...
        'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
    if ~isempty(locs)
        locs = sort(locs, 'ascend');
        p_list{ii} = omega_scan(locs(:));
    else
        p_list{ii} = [];
    end
end

% Convert cell list to numeric matrix with NaN padding
max_modes = 0;
for j = 1:num_fields
    max_modes = max(max_modes, numel(p_list{j}));
end
poles = nan(max_modes, num_fields);
for j = 1:num_fields
    pj = p_list{j};
    if ~isempty(pj)
        poles(1:numel(pj), j) = pj(:);
    end
end

% Axis-resolved modes from Im χ^RPA_{aa} peaks
mx = cell(1, num_fields); my = cell(1, num_fields); mz = cell(1, num_fields);
for ii = 1:num_fields
    sigx = im_diag(:, ii, 1);
    sigy = im_diag(:, ii, 2);
    sigz = im_diag(:, ii, 3);
    [~, lx] = findpeaks(sigx, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
    [~, ly] = findpeaks(sigy, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
    [~, lz] = findpeaks(sigz, 'MinPeakProminence', options.min_prom, 'MinPeakDistance', options.min_sep);
    mx{ii} = omega_scan(lx(:));
    my{ii} = omega_scan(ly(:));
    mz{ii} = omega_scan(lz(:));
end

% helper: convert cell vector of varying-length vectors to NaN-padded matrix
function M = pad_cell_vectors(C)
    n = numel(C);
    k = 0; for j = 1:n, k = max(k, numel(C{j})); end
    M = nan(k, n);
    for j = 1:n
        v = C{j};
        if ~isempty(v)
            M(1:numel(v), j) = v(:);
        end
    end
end

modes.xx = pad_cell_vectors(mx);
modes.yy = pad_cell_vectors(my);
modes.zz = pad_cell_vectors(mz);
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


% (removed) compute_transition_matrices: superseded by vectorized build in RPA_poles


function fig = plot_rpa_modes(field, pole_frequencies, plot_options)
if nargin < 4, plot_options = struct(); end
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
