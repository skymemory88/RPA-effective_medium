function invz_check_solve_opts(sxtra)
%INVZ_CHECK_SOLVE_OPTS Error if solve_opts sets a driver-owned field.
% Shared by invz_spectra_map/invz_spectra_qpath: J0eff/Jxx0/hyp are computed by
% the driver from the BZ coupling sum and must not be overridden via solve_opts.
% ordered_mode is likewise driver-owned (stage-2 P1-6): the map's 1/z leg sets it
% internally (opts.ordered_1z); the auto/overlay leg must never see it.
if any(isfield(sxtra, {'J0eff', 'Jxx0', 'hyp', 'ordered_mode'}))
    error('invz:solveOpts', ...
          'solve_opts fields J0eff/Jxx0/hyp/ordered_mode are reserved (driver-owned).');
end
end
