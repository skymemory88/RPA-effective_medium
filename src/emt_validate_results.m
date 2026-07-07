function report = emt_validate_results(results, inputs)
% EMT_VALIDATE_RESULTS Run basic validation checks on EMT outputs.

report = struct();
report.all_converged = all(results.converged);
report.converged_fraction = sum(results.converged) / max(numel(results.converged), 1);
report.max_residual = max(results.info.residual);
report.max_closure = max(results.info.closure);

% Symmetry check on K and local chi.
max_k_asym = 0;
max_chi_asym = 0;
for ii = 1:size(results.K_emt, 4)
    for iw = 1:size(results.K_emt, 3)
        k = results.K_emt(:,:,iw,ii);
        c = results.chi_emt(:,:,iw,ii);
        max_k_asym = max(max_k_asym, norm(k - k', 'fro'));
        max_chi_asym = max(max_chi_asym, norm(c - c', 'fro'));
    end
end
report.max_k_asymmetry = max_k_asym;
report.max_chi_asymmetry = max_chi_asym;

% Optional non-interacting reference check if Jq is exactly zero.
if all(abs(inputs.Jq(:)) < 1e-14)
    ref = squeeze(mean(inputs.chi_seed, 5));
    err = max(abs(results.chi_emt(:) - ref(:)));
    report.non_interacting_error = err;
else
    report.non_interacting_error = NaN;
end

end
