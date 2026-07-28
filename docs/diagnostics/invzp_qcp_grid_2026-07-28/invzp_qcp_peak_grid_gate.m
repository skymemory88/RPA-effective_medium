function R = invzp_qcp_peak_grid_gate(opts)
%INVZP_QCP_PEAK_GRID_GATE Near-QCP susceptibility at matched B-Bc offsets.
%
% Uses the solver-grade Bc(N) values retained by invzp_qcp_state_grid_gate,
% then evaluates the real-axis response at one common set of physical field
% offsets.  This separates the grid drift of Bc from the shape of the soft
% mode once each grid is aligned to its own mass root.
if nargin < 1, opts = struct(); end

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

Ns = reshape(getf(opts, 'Ns', [12 16 20 24]), 1, []);
T = getf(opts, 'T', 0.1);
offsets_T = reshape(getf(opts, 'offsets_T', ...
    [-0.08 -0.04 -0.02 -0.01 -0.005 0.005 0.01 0.02 0.04 0.08]), 1, []);
w_GHz = reshape(getf(opts, 'w_GHz', 0:0.002:2.0), [], 1);
eta_meV = getf(opts, 'eta_meV', 5e-5);
use_parallel = getf(opts, 'parallel', true);
dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', true);
dipole = getf(opts, 'dipole', 'bruteforce');

if ~isnumeric(Ns) || any(~isfinite(Ns)) || any(Ns < 2) || ...
        any(Ns ~= round(Ns)) || numel(unique(Ns)) ~= numel(Ns)
    error('invzp:qcpPeakGrid', 'opts.Ns must contain distinct integers >= 2.');
end
if ~(isscalar(T) && isfinite(T) && T > 0)
    error('invzp:qcpPeakGrid', 'opts.T must be a finite positive scalar.');
end
if any(~isfinite(offsets_T)) || any(offsets_T == 0) || any(diff(offsets_T) <= 0)
    error('invzp:qcpPeakGrid', ...
        'opts.offsets_T must be strictly increasing, finite, and exclude zero.');
end
if numel(w_GHz) < 3 || any(~isfinite(w_GHz)) || w_GHz(1) < 0 || ...
        any(diff(w_GHz) <= 0) || ...
        max(abs(diff(w_GHz)-diff(w_GHz(1:2)))) > 64*eps(max(w_GHz))
    error('invzp:qcpPeakGrid', ...
        'opts.w_GHz must be a finite, nonnegative, uniformly spaced increasing grid.');
end
if ~(isscalar(eta_meV) && isfinite(eta_meV) && eta_meV > 0)
    error('invzp:qcpPeakGrid', 'opts.eta_meV must be a finite positive scalar.');
end

if isfield(opts, 'Bc_1z')
    Bc = reshape(opts.Bc_1z, 1, []);
else
    q = load(fullfile(here, 'state_grid_gate.mat'), 'R');
    if ~strcmp(q.R.schema, 'invzp_qcp_state_grid_gate/v1')
        error('invzp:qcpPeakGrid', 'state_grid_gate.mat has an unexpected schema.');
    end
    qNs = [q.R.rows.N];
    [tf, loc] = ismember(Ns, qNs);
    if any(~tf) || q.R.config.T ~= T
        error('invzp:qcpPeakGrid', ...
            'Requested grid/T is absent from the retained state gate; supply opts.Bc_1z.');
    end
    Bc = [q.R.rows(loc).Bc_1z];
end
if numel(Bc) ~= numel(Ns) || any(~isfinite(Bc))
    error('invzp:qcpPeakGrid', 'Bc_1z must supply one finite value per requested grid.');
end

C = invz_const();
GHz_per_meV = 1/C.Gh2mV;
w_meV = w_GHz/GHz_per_meV;
ion = invz_ion();

blank = struct('N', NaN, 'Bc_1z', NaN, 'fields_T', [], ...
    'offsets_T', [], 'phase_1z', [], 'phase_1z_reason', {{}}, ...
    'stability_1z', [], 'Epeak_meV', [], 'Epeak_GHz', [], ...
    'Epeak_rpa_GHz', [], 'w_GHz', [], 'chiz', [], 'chirpa', []);
rows = repmat(blank, 1, numel(Ns));

for k = 1:numel(Ns)
    N = Ns(k);
    bz = struct('grid', [N N N], 'dpRng', dpRng, ...
        'cache', cache, 'dipole', dipole);
    [J, info] = invz_bz_couplings(ion, bz);
    fields = Bc(k)+offsets_T;
    mopts = struct('Jnu', J, 'info', info, 'eta', eta_meV, ...
        'parallel', logical(use_parallel), 'verbose', false, ...
        'field_dir', [1 0 0], ...
        'solve_opts', struct('transverse_mf', 'legacy_x'));
    S = invz_spectra_map(ion, T, fields, w_meV, mopts);

    rows(k) = struct('N', N, 'Bc_1z', Bc(k), ...
        'fields_T', S.fields, 'offsets_T', offsets_T, ...
        'phase_1z', S.phase_1z, 'phase_1z_reason', {S.phase_1z_reason}, ...
        'stability_1z', S.stability_1z, ...
        'Epeak_meV', S.Epeak, 'Epeak_GHz', S.Epeak*GHz_per_meV, ...
        'Epeak_rpa_GHz', S.Epeak_rpa*GHz_per_meV, ...
        'w_GHz', S.w*GHz_per_meV, 'chiz', S.chiz, 'chirpa', S.chirpa);
end

peak = vertcat(rows.Epeak_GHz);
peak_hi = max(peak, [], 1, 'includemissing');
peak_lo = min(peak, [], 1, 'includemissing');
peak_mean = mean(peak, 1, 'includemissing');
expected_phase = ones(size(offsets_T));
expected_phase(offsets_T > 0) = 2;
phase = vertcat(rows.phase_1z);
R = struct('schema', 'invzp_qcp_peak_grid_gate/v1', ...
    'created_utc', char(datetime('now', 'TimeZone', 'UTC', ...
                                'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')), ...
    'grid_route', 'production_legacy_absent_policy', ...
    'alignment', 'B-Bc_1z(N)', ...
    'config', struct('Ns', Ns, 'T', T, 'Bc_1z', Bc, ...
        'offsets_T', offsets_T, 'w_GHz', w_GHz, 'eta_meV', eta_meV, ...
        'parallel', logical(use_parallel), 'dpRng', dpRng, ...
    'cache', cache, 'dipole', dipole), ...
    'rows', rows, ...
    'expected_phase_1z', expected_phase, ...
    'all_phase_matched', all(phase == expected_phase, 'all'), ...
    'all_peak_finite', all(isfinite(peak), 'all'), ...
    'peak_spread_GHz', peak_hi-peak_lo, ...
    'peak_mean_GHz', peak_mean, ...
    'peak_rel_spread', (peak_hi-peak_lo)./max(peak_mean, eps));
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
