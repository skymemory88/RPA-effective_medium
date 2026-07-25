function [value, completed, error_id] = invz_try_solver_call(fn)
%INVZ_TRY_SOLVER_CALL Outer dispatcher boundary: absorb only reviewed recoverable signals.
% Strict-medium domain outcomes are statuses and do not arrive here.
%
% CONTRACT: fn must be callable for ONE output. This function evaluates `value = fn()`, so a
% VOID fn makes MATLAB raise MATLAB:maxlhs before fn's own error can escape -- which would mask
% a genuine recoverable signal (invz:orderedPhase / invz:degenerateDoublet) as an unrecognised
% MATLAB error and rethrow it. That is the exact inversion this boundary exists to prevent, so a
% void fn is rejected up front as a wiring error rather than silently mishandled.
%
% The rejection catches NAMED void functions, for which nargout is 0. It cannot catch the
% anonymous void form `@() error(...)`: nargout reports -1 for any anonymous handle, so the two
% are indistinguishable before the call. Callers writing a throwing test double must therefore
% declare an output (`function v = thrower(id)`), not use a bare `@() error(...)`.
if ~isa(fn, 'function_handle')
    error('invz:solverCall', 'fn must be a function handle.');
end
% nargout() itself throws for a direct handle to an unresolvable name (a stale, mistyped or
% str2func-built one): handle CONSTRUCTION succeeds, resolution fails only when probed. Left
% unguarded, that raw MATLAB:narginout:* identifier would escape this boundary un-namespaced --
% the very leak the check below exists to prevent. Report it as the wiring error it is, keeping
% the original text so the cause stays diagnosable.
try
    n_out = nargout(fn);
catch probe_err
    error('invz:solverCall', ...
        'fn is not a resolvable function handle: %s (%s)', ...
        probe_err.message, probe_err.identifier);
end
if n_out == 0
    error('invz:solverCall', ['fn must return a value; a void fn surfaces as MATLAB:maxlhs ' ...
        'and masks its own identifier, so a recoverable signal would be rethrown as fatal.']);
end
try
    value = fn();
    completed = true;
    error_id = '';
catch err
    if ~invz_is_recoverable_solver_error(err.identifier)
        rethrow(err);
    end
    value = [];
    completed = false;
    error_id = err.identifier;       % exact category; never collapse to "nonconverged"
end
end
