function outputs = emt_run_from_workspace(data, emt_params)
% EMT_RUN_FROM_WORKSPACE High-level EMT runner for workspace-style inputs.

if nargin < 2
    emt_params = struct();
end

defaults = emt_default_params();
params = emt_merge_params(defaults, emt_params);
inputs = emt_prepare_workspace_inputs(data);

fprintf('\n=== Real-frequency EMT (rewritten backbone) ===\n');
fprintf('Scan mode: %s | points: %d | omega: %d | q: %d\n', ...
    inputs.scan_mode, inputs.n_cVar, inputs.n_omega, inputs.n_q);

results = emt_solve_scan(inputs, params);
report = emt_validate_results(results, inputs);

fprintf('Converged: %d / %d points (%.1f%%)\n', ...
    results.info.n_converged, results.info.n_total, 100 * report.converged_fraction);
fprintf('Max residual: %.3e | Max closure residual: %.3e\n', ...
    report.max_residual, report.max_closure);

outputs = results;
outputs.validation = report;
outputs.inputs = inputs;

end
