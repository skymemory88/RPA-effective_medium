function tmf = invz_check_transverse_mf(sxtra, By)
%INVZ_CHECK_TRANSVERSE_MF Resolve transverse_mf and guard legacy_x against a b-axis field.
% Shared by the trusted-result entry points (invz_spectra_map, invz_spectra_qpath,
% invz_chi_tensor_ref): under 'legacy_x' (x-only mean field) a nonzero b-axis (y)
% field component is C4-inconsistent (17 ueV a/b asymmetry at 4 T) and errors
% invz:transverseMF; use 'vector_ab' for genuine in-plane rotation. By may be a
% scalar or a column of per-field y-components (e.g. BvecM(:,2)).
tmf = getf(sxtra, 'transverse_mf', 'legacy_x');
if strcmp(tmf, 'legacy_x') && any(abs(By(:)) > 0)
    error('invz:transverseMF', ['field has a b-axis (y) component but transverse_mf is ' ...
        '''legacy_x'' (x-only mean field; C4-inconsistent, 17 ueV a/b asymmetry at 4 T). ' ...
        'Set opts.solve_opts.transverse_mf = ''vector_ab'' (or ''none'' for bare diagnostics).']);
end
end
