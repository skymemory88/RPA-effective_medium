function R = invzp_qcp_coupling_scan(opts)
%INVZP_QCP_COUPLING_SCAN Coupling-only predictor for the ordered-QCP grid gate.
%
% This diagnostic intentionally uses the production legacy grid route:
% opts.grid is present, but gridConvention/gridOffset/gammaPolicy are absent.
% Adding any of those policy fields changes invz_bz_couplings' grid builder and
% would confound the current 16^3 spectra baseline.
%
% Returned quantities are diagnostics, not acceptance gates:
%   W              = Jmax-Jmin
%   control_2_over_W = 2/W
%   S_edge         = mean(1./(J-J0)) at the excluded-Gamma continuum edge
%   S_fixed        = S_N(J0+delta) at fixed positive delta
%   S_scaled       = S_N(J0+c*dlevel), with dlevel the top distinct-level gap
%
% R.fit reports simple linear/quadratic 1/N intercepts. Four grids do not
% provide a rigorous continuum error bound; model sensitivity must be retained.
if nargin < 1, opts = struct(); end

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

Ns = getf(opts, 'Ns', [12 16 20 24]);
if ~isnumeric(Ns) || ~isvector(Ns) || any(~isfinite(Ns)) || ...
        any(Ns < 2) || any(Ns ~= round(Ns)) || numel(unique(Ns)) ~= numel(Ns)
    error('invzp:qcpGrid', 'opts.Ns must contain distinct integer grid sizes >= 2.');
end
Ns = reshape(Ns, 1, []);

dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', true);
dipole = getf(opts, 'dipole', 'bruteforce');
delta_u_eV = getf(opts, 'delta_u_eV', [0 0.025 0.05 0.1 0.2 0.5 1.0]);
scaled_c = getf(opts, 'scaled_c', [1 2 4 8 16]);
if any(~isfinite(delta_u_eV)) || any(delta_u_eV < 0)
    error('invzp:qcpGrid', 'opts.delta_u_eV must be finite and nonnegative.');
end
if any(~isfinite(scaled_c)) || any(scaled_c <= 0)
    error('invzp:qcpGrid', 'opts.scaled_c must be finite and positive.');
end

ion = invz_ion();
blank = struct('N', NaN, 'n_modes', NaN, 'J0', NaN, 'Jaa0', NaN, ...
    'Jmin', NaN, 'Jmax', NaN, 'Jbar', NaN, 'W', NaN, ...
    'control_2_over_W', NaN, 'gamma_gap', NaN, ...
    'edge_multiplicity', NaN, 'Jnext', NaN, 'dlevel', NaN, ...
    'delta_fixed_meV', [], 'S_fixed', [], ...
    'delta_scaled_meV', [], 'S_scaled', [], ...
    'S_edge', NaN, 'sigma_edge_identity', NaN, 'sigma_crit_helper', NaN, ...
    'sigma_identity_absdiff', NaN);
rows = repmat(blank, 1, numel(Ns));

for k = 1:numel(Ns)
    N = Ns(k);
    bz = struct('grid', [N N N], 'dpRng', dpRng, ...
                'cache', cache, 'dipole', dipole);
    [J, info, Jaa0] = invz_bz_couplings(ion, bz);
    J = J(:);

    Jmin = min(J);
    Jmax = max(J);
    W = Jmax-Jmin;
    tol_top = 128*eps(max(1, abs(Jmax)));
    below = J(J < Jmax-tol_top);
    if isempty(below)
        error('invzp:qcpGrid', ...
            'N=%d has no distinct coupling level below Jmax within the declared tolerance.', N);
    end
    Jnext = max(below);
    dlevel = Jmax-Jnext;
    J0 = info.Jcc0;
    if J0 < Jmax-tol_top
        error('invzp:qcpGrid', ...
            'N=%d has J0=%.17g below Jmax=%.17g; excluded-edge assumption failed.', ...
            N, J0, Jmax);
    end

    delta_fixed = reshape(delta_u_eV, 1, [])*1e-3; % micro-eV -> meV
    y_fixed = J0 + delta_fixed;
    S_fixed = arrayfun(@(y) mean(1./(J-y)), y_fixed);

    delta_scaled = reshape(scaled_c, 1, [])*dlevel;
    y_scaled = J0 + delta_scaled;
    S_scaled = arrayfun(@(y) mean(1./(J-y)), y_scaled);

    Sc = invz_sigma_crit(J0, J);
    Sc_id = -J0*S_fixed(1)-1;

    rows(k) = struct('N', N, 'n_modes', numel(J), 'J0', J0, 'Jaa0', Jaa0, ...
        'Jmin', Jmin, 'Jmax', Jmax, 'Jbar', mean(J), 'W', W, ...
        'control_2_over_W', 2/W, 'gamma_gap', J0-Jmax, ...
        'edge_multiplicity', nnz(abs(J-Jmax) <= tol_top), ...
        'Jnext', Jnext, 'dlevel', dlevel, ...
        'delta_fixed_meV', delta_fixed, 'S_fixed', S_fixed, ...
        'delta_scaled_meV', delta_scaled, 'S_scaled', S_scaled, ...
        'S_edge', S_fixed(1), 'sigma_edge_identity', Sc_id, ...
        'sigma_crit_helper', Sc, 'sigma_identity_absdiff', abs(Sc_id-Sc));
end

x = 1./double(Ns);
fit_fixed = repmat(struct('delta_meV', NaN, 'S_inf_linear', NaN, ...
    'S_inf_quadratic', NaN, 'intercept_spread', NaN), 1, numel(delta_u_eV));
for j = 1:numel(delta_u_eV)
    v = arrayfun(@(r) r.S_fixed(j), rows);
    p1 = polyfit(x, v, 1);
    p2 = polyfit(x, v, 2);
    fit_fixed(j) = struct('delta_meV', rows(1).delta_fixed_meV(j), ...
        'S_inf_linear', p1(end), 'S_inf_quadratic', p2(end), ...
        'intercept_spread', abs(p1(end)-p2(end)));
end

gap = [rows.gamma_gap];
pg = polyfit(log(double(Ns)), log(gap), 1);
fit = struct('variable', '1/N', 'fixed_delta', fit_fixed, ...
    'gamma_gap_power_exponent', -pg(1), ...
    'gamma_gap_power_prefactor', exp(pg(2)));

R = struct('schema', 'invzp_qcp_coupling_scan/v1', ...
    'created_utc', char(datetime('now', 'TimeZone', 'UTC', ...
                                'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')), ...
    'grid_route', 'production_legacy_absent_policy', ...
    'config', struct('Ns', Ns, 'dpRng', dpRng, 'cache', cache, ...
                     'dipole', dipole, 'delta_u_eV', delta_u_eV, ...
                     'scaled_c', scaled_c), ...
    'rows', rows, 'fit', fit);
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
