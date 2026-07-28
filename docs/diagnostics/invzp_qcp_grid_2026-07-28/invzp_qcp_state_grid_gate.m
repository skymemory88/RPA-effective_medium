function R = invzp_qcp_state_grid_gate(opts)
%INVZP_QCP_STATE_GRID_GATE State-only grid test of the ordered computability edge.
%
% This fixture deliberately does not evaluate real-axis susceptibility.  For
% each legacy production N^3 grid it:
%   1. obtains a solver-grade PM mass root Bc from converged invz_solve_point
%      evaluations;
%   2. runs the Jensen ordered state solver on one common physical-field mesh;
%   3. records both the total accepted count and the contiguous accepted suffix
%      immediately below Bc.
%
% Counts are meaningful only together with config.fields and config.field_step.
% The common mesh, rather than a per-grid normalized mesh, is what makes a
% proposed O(1/N) physical-field shrinkage falsifiable.
if nargin < 1, opts = struct(); end

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

Ns = reshape(getf(opts, 'Ns', [12 16 20 24]), 1, []);
T = getf(opts, 'T', 0.1);
field_step = getf(opts, 'field_step', 0.025);
fields = reshape(getf(opts, 'fields', 4.0:field_step:4.7), 1, []);
bc_bracket = reshape(getf(opts, 'bc_bracket', [4.55 4.85]), 1, []);
bc_tol = getf(opts, 'bc_tol', 2e-5);
use_parallel = getf(opts, 'parallel', true);
dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', true);
dipole = getf(opts, 'dipole', 'bruteforce');

if ~isnumeric(Ns) || ~isvector(Ns) || any(~isfinite(Ns)) || ...
        any(Ns < 2) || any(Ns ~= round(Ns)) || numel(unique(Ns)) ~= numel(Ns)
    error('invzp:qcpGrid', 'opts.Ns must contain distinct integer grid sizes >= 2.');
end
if ~(isscalar(T) && isfinite(T) && T > 0)
    error('invzp:qcpGrid', 'opts.T must be a finite positive scalar.');
end
if numel(fields) < 2 || any(~isfinite(fields)) || any(diff(fields) <= 0)
    error('invzp:qcpGrid', 'opts.fields must be a strictly increasing finite vector.');
end
if ~(isscalar(field_step) && isfinite(field_step) && field_step > 0) || ...
        max(abs(diff(fields)-field_step)) > 64*eps(max(abs(fields)))
    error('invzp:qcpGrid', ...
        'opts.fields must use the declared positive constant opts.field_step.');
end
if numel(bc_bracket) ~= 2 || any(~isfinite(bc_bracket)) || ...
        bc_bracket(1) >= bc_bracket(2)
    error('invzp:qcpGrid', 'opts.bc_bracket must be [Blo Bhi] with Blo < Bhi.');
end
if ~(isscalar(bc_tol) && isfinite(bc_tol) && bc_tol > 0)
    error('invzp:qcpGrid', 'opts.bc_tol must be a finite positive scalar.');
end

ion = invz_ion();
blank = struct('N', NaN, 'n_modes', NaN, ...
    'J0', NaN, 'Jmin', NaN, 'Jmax', NaN, 'W', NaN, ...
    'control_2_over_W', NaN, 'gamma_gap', NaN, 'dlevel', NaN, ...
    'Bc_1z', NaN, 'Bc_bracket_final', nan(1,2), ...
    'Bc_crit_final', nan(1,2), 'Bc_iterations', NaN, ...
    'fields', [], 'ordered', struct([]), ...
    'accepted_total_count', NaN, 'stable_total_count', NaN, ...
    'eligible_ordered_count', NaN, ...
    'accepted_qcp_contiguous_count', NaN, ...
    'accepted_qcp_contiguous_width_T', NaN, ...
    'accepted_qcp_lowest_field_T', NaN);
rows = repmat(blank, 1, numel(Ns));

if use_parallel
    pool = gcp('nocreate');
    if isempty(pool), parpool('local'); end
end

for k = 1:numel(Ns)
    N = Ns(k);
    bz = struct('grid', [N N N], 'dpRng', dpRng, ...
                'cache', cache, 'dipole', dipole);
    [J, info, Jaa0] = invz_bz_couplings(ion, bz);
    J = J(:);
    sopts = struct('hyp', true, 'J0eff', info.Jcc0, 'Jxx0', Jaa0, ...
        'transverse_mf', 'legacy_x', 'static_medium', 'resummed');

    [Bc, bc_final, crit_final, bc_iters] = ...
        mass_root(ion, T, J, sopts, bc_bracket, bc_tol);

    probes = repmat(blank_probe(), 1, numel(fields));
    if use_parallel
        parfor j = 1:numel(fields)
            soj = sopts;
            soj.ordered_mode = 'jensen';
            probes(j) = ordered_probe(ion, T, fields(j), J, soj);
        end
    else
        for j = 1:numel(fields)
            soj = sopts;
            soj.ordered_mode = 'jensen';
            probes(j) = ordered_probe(ion, T, fields(j), J, soj);
        end
    end

    accepted = [probes.accepted];
    stable = [probes.stable];
    eligible = fields < Bc;
    qcp_idx = find(eligible);
    n_suffix = 0;
    for j = fliplr(qcp_idx)
        if ~accepted(j), break; end
        n_suffix = n_suffix + 1;
    end
    if n_suffix == 0
        suffix_width = 0;
        suffix_low = NaN;
    else
        suffix_low = fields(qcp_idx(end)-n_suffix+1);
        suffix_width = Bc-suffix_low;
    end

    Jmax = max(J);
    below = J(J < Jmax-128*eps(max(1, abs(Jmax))));
    if isempty(below)
        error('invzp:qcpGrid', 'N=%d has no distinct level below Jmax.', N);
    end
    Jmin = min(J);
    rows(k) = struct('N', N, 'n_modes', numel(J), ...
        'J0', info.Jcc0, 'Jmin', Jmin, 'Jmax', Jmax, ...
        'W', Jmax-Jmin, 'control_2_over_W', 2/(Jmax-Jmin), ...
        'gamma_gap', info.Jcc0-Jmax, 'dlevel', Jmax-max(below), ...
        'Bc_1z', Bc, 'Bc_bracket_final', bc_final, ...
        'Bc_crit_final', crit_final, 'Bc_iterations', bc_iters, ...
        'fields', fields, 'ordered', probes, ...
        'accepted_total_count', nnz(accepted), ...
        'stable_total_count', nnz(stable), ...
        'eligible_ordered_count', nnz(eligible), ...
        'accepted_qcp_contiguous_count', n_suffix, ...
        'accepted_qcp_contiguous_width_T', suffix_width, ...
        'accepted_qcp_lowest_field_T', suffix_low);
end

R = struct('schema', 'invzp_qcp_state_grid_gate/v1', ...
    'created_utc', char(datetime('now', 'TimeZone', 'UTC', ...
                                'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX')), ...
    'grid_route', 'production_legacy_absent_policy', ...
    'acceptance_definition', ...
        'pt.is_ordered && pt.converged && finite(pt.Sigma0)', ...
    'qcp_suffix_definition', ...
        'highest contiguous accepted fields on the common mesh strictly below solver-grade Bc', ...
    'config', struct('Ns', Ns, 'T', T, 'fields', fields, ...
        'field_step', field_step, 'bc_bracket', bc_bracket, ...
        'bc_tol', bc_tol, 'parallel', logical(use_parallel), ...
        'dpRng', dpRng, 'cache', cache, 'dipole', dipole), ...
    'rows', rows);
end

function [Bc, bracket, crits, iters] = mass_root(ion, T, J, sopts, bracket, tol)
lo = bracket(1);
hi = bracket(2);
[flo, oklo] = mass_at(ion, T, lo, J, sopts);
[fhi, okhi] = mass_at(ion, T, hi, J, sopts);
if ~(oklo && okhi && flo <= 0 && fhi >= 0)
    error('invzp:qcpGrid', ...
        ['PM mass root is not bracketed by converged points: ' ...
         'crit(%.8g)=%.8g (ok=%d), crit(%.8g)=%.8g (ok=%d).'], ...
        lo, flo, oklo, hi, fhi, okhi);
end
iters = 0;
while hi-lo > tol
    iters = iters+1;
    mid = 0.5*(lo+hi);
    [fm, okm] = mass_at(ion, T, mid, J, sopts);
    if ~okm
        error('invzp:qcpGrid', ...
            'PM mass solve did not converge at B=%.17g inside a clean bracket.', mid);
    end
    if fm <= 0
        lo = mid;
        flo = fm;
    else
        hi = mid;
        fhi = fm;
    end
end
Bc = 0.5*(lo+hi);
bracket = [lo hi];
crits = [flo fhi];
end

function [crit, ok] = mass_at(ion, T, B, J, sopts)
try
    pt = invz_solve_point(ion, T, B, J, sopts);
    crit = pt.crit;
    ok = pt.converged && isfinite(crit) && isfinite(pt.Sigma0);
catch ME
    if any(strcmp(ME.identifier, {'invz:degenerateDoublet', ...
                                  'invz:mediumOutOfDomain'}))
        crit = NaN;
        ok = false;
    else
        rethrow(ME);
    end
end
end

function p = ordered_probe(ion, T, B, J, sopts)
p = blank_probe();
p.B = B;
try
    pt = invz_solve_point_ordered(ion, T, B, J, sopts);
    p.completed = true;
    p.is_ordered = logical(pt.is_ordered);
    p.converged = logical(pt.converged);
    p.accepted = p.is_ordered && p.converged && isfinite(pt.Sigma0);
    p.stable = isfield(pt, 'stable_1z') && logical(pt.stable_1z);
    p.Sigma0 = pt.Sigma0;
    p.hmf_status = pt.hmf_prof.status;
    if isfield(pt, 'hmf'), p.hstar = pt.hmf; end
    if isfield(pt, 'Dq_min'), p.Dq_min = pt.Dq_min; end
    if isfield(pt, 'D_uni'), p.D_uni = pt.D_uni; end
    if isfield(pt, 'final_resid'), p.final_resid = pt.final_resid; end
    if isfield(pt, 'crit_1z'), p.crit_1z = pt.crit_1z; end
    if isfield(pt, 'hmf_prof') && isstruct(pt.hmf_prof) && ...
            isfield(pt.hmf_prof, 'status_detail') && ...
            isstruct(pt.hmf_prof.status_detail) && ...
            isscalar(pt.hmf_prof.status_detail) && ...
            isfield(pt.hmf_prof.status_detail, 'binding_node') && ...
            isstruct(pt.hmf_prof.status_detail.binding_node) && ...
            isscalar(pt.hmf_prof.status_detail.binding_node)
        b = pt.hmf_prof.status_detail.binding_node;
        p.binding_node_id = b.node_id;
        p.binding_h = b.h;
        p.binding_term_reason = b.term_reason;
        p.binding_medium_status = b.medium_status;
        p.binding_Dq_min = b.Dq_min;
        p.binding_Dq_abs_min = b.Dq_abs_min;
        p.binding_resid = [b.resid_A b.resid_B b.resid_C b.resid_D b.resid_static];
    end
catch ME
    p.error_id = ME.identifier;
    p.error_message = ME.message;
end
end

function p = blank_probe()
p = struct('B', NaN, 'completed', false, 'is_ordered', false, ...
    'converged', false, 'accepted', false, 'stable', false, ...
    'hmf_status', 'not_attempted', 'hstar', NaN, 'Sigma0', NaN, ...
    'crit_1z', NaN, 'D_uni', NaN, 'Dq_min', NaN, 'final_resid', NaN, ...
    'binding_node_id', NaN, 'binding_h', NaN, ...
    'binding_term_reason', 'not_evaluated', ...
    'binding_medium_status', 'not_applicable', ...
    'binding_Dq_min', NaN, 'binding_Dq_abs_min', NaN, ...
    'binding_resid', nan(1,5), 'error_id', '', 'error_message', '');
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
