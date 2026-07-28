function invz_check_solve_opts(sxtra)
%INVZ_CHECK_SOLVE_OPTS Error if solve_opts sets a driver-owned field.
% Shared by invz_spectra_map/invz_spectra_qpath: J0eff/Jxx0/hyp are computed by
% the driver from the BZ coupling sum and must not be overridden via solve_opts.
% ordered_mode is likewise driver-owned (stage-2 P1-6): the map's 1/z leg sets it
% internally (opts.ordered_1z); the auto/overlay leg must never see it.
% static_medium/ref_margin are driver-owned at THIS level too (spec SS4.2): the driver
% resolves the scheme ONCE with invz_check_static_medium and stamps it into every
% per-field solve, so a value smuggled in through solve_opts could leave the sweep's
% columns on two different truncation orders while S.static_medium still advertised one.
% hmf_seed is also driver-owned: invz_spectra_map's optional physical-field continuation
% updates it between columns; a caller-supplied fixed seed would have false provenance.
if any(isfield(sxtra, ...
        {'J0eff', 'Jxx0', 'hyp', 'ordered_mode', 'static_medium', 'ref_margin', 'hmf_seed'}))
    error('invz:solveOpts', ...
          ['solve_opts fields J0eff/Jxx0/hyp/ordered_mode/static_medium/ref_margin/hmf_seed ' ...
           'are reserved (driver-owned).']);
end
end
