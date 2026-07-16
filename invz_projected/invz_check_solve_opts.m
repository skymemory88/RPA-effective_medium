function invz_check_solve_opts(sxtra)
%INVZ_CHECK_SOLVE_OPTS Error if solve_opts sets a driver-owned field.
% Shared by invz_spectra_map/invz_spectra_qpath: J0eff/Jxx0/hyp are computed by
% the driver from the BZ coupling sum and must not be overridden via solve_opts.
if any(isfield(sxtra, {'J0eff', 'Jxx0', 'hyp'}))
    error('invz:solveOpts', 'solve_opts fields J0eff/Jxx0/hyp are reserved (driver-owned).');
end
end
