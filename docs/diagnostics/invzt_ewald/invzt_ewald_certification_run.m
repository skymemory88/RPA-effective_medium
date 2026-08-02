function report = invzt_ewald_certification_run(output_file)
%INVZT_EWALD_CERTIFICATION_RUN  Execute the pre-registered tensor Ewald gates.
%
% The gate definitions and fixed tolerances live in
% docs/execution/invzt_ewald_certification.md. This runner does not alter the
% backend default; it evaluates explicit Ewald candidates and retains only
% compact metrics, never the production-grid coupling pages.

if nargin < 1, output_file = ''; end
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root);
addpath(fullfile(root, 'invz_common'));
addpath(fullfile(root, 'invz_tensor'));

ion = invz_ion();
candidate = invzt_ewald_defaults(ion);
alpha0 = candidate.alpha;
refined = struct('alpha', alpha0, 'r_cut', 6/alpha0, ...
    'g_cut', 12*alpha0, 'boundary', 'conducting_k0_omitted');

report = struct();
report.schema = 'invzt_ewald_certification/v1';
report.generated_utc = char(datetime('now', 'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
report.matlab_version = version;
report.ion_Vc = ion.Vc;
report.candidate = candidate;
report.joint_refinement = refined;
report.tolerances = struct('rel_tensor', 1e-8, 'abs_tensor_scale', 1e-8, ...
    'J_ref_meV', 0.006424435656, 'rel_coupling', 1e-8, ...
    'abs_coupling_Jref', 1e-8, 'solver_abs', 1e-7);

qfrozen = frozen_q_sample();
[lat0, t0] = timed_lattice(ion, qfrozen, candidate);
[lat1, t1] = timed_lattice(ion, qfrozen, refined);
report.frozen_sample = compare_lattices(lat0, lat1);
report.frozen_sample.nq = size(qfrozen, 1);
report.frozen_sample.seconds = [t0 t1];

grid = invzt_qgrid(16, 'halfopen');
[lat16, t16] = timed_lattice(ion, grid, candidate);
[lat16r, t16r] = timed_lattice(ion, grid, refined);
report.production_grid = compare_lattices(lat16, lat16r);
report.production_grid.nq = size(grid.qvec, 1);
report.production_grid.seconds = [t16 t16r];

[~, counts0] = invz_dipole_ewald([0 0 0], ion.a, ion.tau, candidate);
[~, counts1] = invz_dipole_ewald([0 0 0], ion.a, ion.tau, refined);
report.preflight = struct( ...
    'candidate_estimated_peak_bytes', counts0.preflight.estimated_peak_bytes, ...
    'refined_estimated_peak_bytes', counts1.preflight.estimated_peak_bytes, ...
    'candidate_recip_candidates', counts0.recip_candidates, ...
    'refined_recip_candidates', counts1.recip_candidates, ...
    'candidate_max_real_pair', max(counts0.real_pair(:)), ...
    'refined_max_real_pair', max(counts1.real_pair(:)));

tic;
pt0 = invzt_solve_point(ion, 0.1, [6 0 0], lat16, struct());
ts0 = toc;
tic;
pt1 = invzt_solve_point(ion, 0.1, [6 0 0], lat16r, struct());
ts1 = toc;
report.solver_anchor = struct( ...
    'temperature_K', 0.1, 'field_T', [6 0 0], ...
    'candidate_converged', logical(pt0.converged), ...
    'refined_converged', logical(pt1.converged), ...
    'candidate_iters', pt0.outer_iters, 'refined_iters', pt1.outer_iters, ...
    'seconds', [ts0 ts1], ...
    'Sigma0_abs_error', abs(pt0.Sigma0-pt1.Sigma0), ...
    'crit_abs_error', abs(pt0.crit-pt1.crit), ...
    'Sigma_max_abs_error', max_abs_diff(pt0.Sigma, pt1.Sigma), ...
    'K_max_abs_error', max_abs_diff(pt0.K, pt1.K), ...
    'G_max_abs_error', max_abs_diff(pt0.G, pt1.G), ...
    'candidate_Sigma0', pt0.Sigma0, 'refined_Sigma0', pt1.Sigma0, ...
    'candidate_crit', pt0.crit, 'refined_crit', pt1.crit);

dp = [10 20 30];
brute = repmat(struct('dpRng', [], 'seconds', [], 'Jcc0', [], 'Jaa0', [], ...
    'max_eig_abs_from_ewald', []), 1, numel(dp));
eig0 = sorted_page_eigs(lat0.Jt);
for k = 1:numel(dp)
    tic;
    latb = invzt_jq_tensor(ion, qfrozen, struct('dipole', 'bruteforce', ...
        'dpRng', dp(k), 'cache', false));
    tb = toc;
    brute(k).dpRng = dp(k);
    brute(k).seconds = tb;
    brute(k).Jcc0 = latb.info.Jcc0;
    brute(k).Jaa0 = latb.info.Jaa0;
    brute(k).max_eig_abs_from_ewald = max_abs_diff(sorted_page_eigs(latb.Jt), eig0);
end
report.bruteforce_diagnostic = brute;

report.gates = struct();
report.gates.G3_frozen = report.frozen_sample.tensor.pass ...
    && report.frozen_sample.eigenvalues.pass ...
    && report.frozen_sample.gamma.pass;
report.gates.G3_production_grid = report.production_grid.tensor.pass ...
    && report.production_grid.eigenvalues.pass ...
    && report.production_grid.gamma.pass;
s = report.solver_anchor;
report.gates.G5_solver_anchor = s.candidate_converged && s.refined_converged ...
    && s.Sigma0_abs_error <= 1e-7 && s.crit_abs_error <= 1e-7 ...
    && s.Sigma_max_abs_error <= 1e-7 && s.K_max_abs_error <= 1e-7 ...
    && s.G_max_abs_error <= 1e-7;
report.pass = report.gates.G3_frozen && report.gates.G3_production_grid ...
    && report.gates.G5_solver_anchor;

if ~isempty(output_file)
    outdir = fileparts(output_file);
    if ~isempty(outdir) && ~exist(outdir, 'dir'), mkdir(outdir); end
    save(output_file, 'report');
end

fprintf(['Ewald certification: frozen=%d grid16=%d solver6T=%d overall=%d\n' ...
    '  tensor ratios frozen/grid = %.3g / %.3g\n' ...
    '  coupling ratios frozen/grid = %.3g / %.3g\n' ...
    '  solver errors Sigma0/crit/Sigma/K/G = %.3g %.3g %.3g %.3g %.3g\n'], ...
    report.gates.G3_frozen, report.gates.G3_production_grid, ...
    report.gates.G5_solver_anchor, report.pass, ...
    report.frozen_sample.tensor.worst_ratio, ...
    report.production_grid.tensor.worst_ratio, ...
    report.frozen_sample.eigenvalues.worst_ratio, ...
    report.production_grid.eigenvalues.worst_ratio, ...
    s.Sigma0_abs_error, s.crit_abs_error, s.Sigma_max_abs_error, ...
    s.K_max_abs_error, s.G_max_abs_error);

if ~report.pass
    error('invzt:ewaldCertification', ...
        'One or more pre-registered Ewald certification gates failed.');
end
end

function [lat, seconds] = timed_lattice(ion, q, eopts)
tic;
lat = invzt_jq_tensor(ion, q, struct('dipole', 'ewald', ...
    'ewald', eopts, 'cache', false));
seconds = toc;
end

function out = compare_lattices(a, b)
out = struct();
out.tensor = metric_tensor_pages(a.Jt, b.Jt);
out.eigenvalues = metric_coupling(sorted_page_eigs(a.Jt), sorted_page_eigs(b.Jt));
out.gamma = metric_coupling([a.info.Jcc0; a.info.Jaa0], ...
    [b.info.Jcc0; b.info.Jaa0]);
out.Jcc0 = [a.info.Jcc0 b.info.Jcc0];
out.Jaa0 = [a.info.Jaa0 b.info.Jaa0];
end

function out = metric_tensor_pages(A, B)
assert(isequal(size(A), size(B)));
nq = size(A, 3);
worstMargin = -Inf;
worstRatio = 0;
maxError = 0;
for iq = 1:nq
    aa = A(:,:,iq); bb = B(:,:,iq);
    scale = max(max(abs(aa(:))), max(abs(bb(:))));
    err = abs(aa-bb);
    allowed = 1e-8*scale + 1e-8*max(abs(aa),abs(bb));
    worstMargin = max(worstMargin, max(err(:)-allowed(:)));
    worstRatio = max(worstRatio, max(err(:)./max(allowed(:),realmin)));
    maxError = max(maxError, max(err(:)));
end
out = struct('pass', worstMargin <= 0, 'worst_margin', worstMargin, ...
    'worst_ratio', worstRatio, 'max_abs_error', maxError);
end

function out = metric_coupling(A, B)
assert(isequal(size(A), size(B)));
Jref = 0.006424435656;
err = abs(A-B);
allowed = 1e-8*Jref + 1e-8*max(abs(A),abs(B));
out = struct('pass', max(err(:)-allowed(:)) <= 0, ...
    'worst_margin', max(err(:)-allowed(:)), ...
    'worst_ratio', max(err(:)./max(allowed(:),realmin)), ...
    'max_abs_error', max(err(:)));
end

function E = sorted_page_eigs(J)
nq = size(J, 3);
E = zeros(12, nq);
for iq = 1:nq
    page = (J(:,:,iq)+J(:,:,iq)')/2;
    E(:,iq) = sort(real(eig(page)));
end
end

function q = frozen_q_sample()
dq = 2^-40;
edge = [-0.5, 0.5-dq];
q = [0 0 0; 0.137 0.291 0.453; -0.311 0.173 -0.227; 0.25 0 0.1];
for axis = 1:3
    other = setdiff(1:3, axis, 'stable');
    for side = edge
        row = zeros(1,3); row(axis) = side; row(other) = [0.137 0.291];
        q(end+1,:) = row; %#ok<AGROW>
    end
end
for freeAxis = 1:3
    fixed = setdiff(1:3, freeAxis, 'stable');
    for i = 1:2
        for j = 1:2
            row = zeros(1,3); row(freeAxis) = 0.173;
            row(fixed) = [edge(i) edge(j)];
            q(end+1,:) = row; %#ok<AGROW>
        end
    end
end
[x,y,z] = ndgrid(edge, edge, edge);
q = [q; x(:) y(:) z(:)];
assert(isequal(size(q), [30 3]));
end

function d = max_abs_diff(a, b)
assert(isequal(size(a), size(b)));
d = max(abs(a(:)-b(:)));
end
