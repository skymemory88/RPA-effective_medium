function tests = test_invz_is_recoverable_solver_error
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

% The whitelist is exactly the two existing branch/domain signals (spec §5.1).
function test_whitelisted_ids_are_recoverable(testCase)
verifyTrue(testCase, invz_is_recoverable_solver_error('invz:orderedPhase'));
verifyTrue(testCase, invz_is_recoverable_solver_error('invz:degenerateDoublet'));
end

% Wiring errors must NOT be recoverable -- this is the whole point. A masked column caused by a
% swallowed wiring error is the exact failure mode that hid the original defect for a stage.
function test_wiring_ids_are_not_recoverable(testCase)
for id = {'invz:staticMedium', 'invz:emtJnu', 'invz:hzFixed', 'invz:nodeSolveNode', ...
          'invz:residualNode', 'invz:residualState', 'invz:couplingMoments', ...
          'invz:transverseMF', 'invz:hmfOpts', 'invz:oddArgs'}
    verifyFalse(testCase, invz_is_recoverable_solver_error(id{1}), ...
        sprintf('%s must NOT be recoverable', id{1}));
end
end

% Unclassified ids rethrow by DEFAULT -- including unseen invz:* ones. Adding to the whitelist
% is a reviewed contract change, not a response to a failing run.
function test_unknown_invz_id_is_not_recoverable(testCase)
verifyFalse(testCase, invz_is_recoverable_solver_error('invz:somethingBrandNew'));
verifyFalse(testCase, invz_is_recoverable_solver_error('MATLAB:nomem'));
verifyFalse(testCase, invz_is_recoverable_solver_error(''));
end

function test_non_char_input_is_not_recoverable(testCase)
verifyFalse(testCase, invz_is_recoverable_solver_error(42));
verifyFalse(testCase, invz_is_recoverable_solver_error([]));
end

function test_string_class_accepted(testCase)
verifyTrue(testCase, invz_is_recoverable_solver_error("invz:orderedPhase"));
end

% The classifier must be usable directly on an MException identifier.
function test_works_on_a_caught_exception(testCase)
try
    error('invz:degenerateDoublet', 'synthetic');
catch err
    verifyTrue(testCase, invz_is_recoverable_solver_error(err.identifier));
end
try
    error('invz:staticMedium', 'synthetic');
catch err2
    verifyFalse(testCase, invz_is_recoverable_solver_error(err2.identifier));
end
end

% A thrower that DECLARES an output. This matters: invz_try_solver_call requests one output
% (`value = fn()`), and a VOID body such as a bare @() error(...) makes MATLAB raise
% MATLAB:maxlhs BEFORE the intended error escapes -- masking the very identifier under test.
% Declaring the output lets the body run, so err.identifier is the real one. Verified:
% nargout(@() error(...)) is -1, so no nargout test can rescue the anonymous void form.
function v = local_thrower(id)
error(id, 'synthetic');
v = [];   %#ok<UNRCH> -- never reached; present so the handle declares an output
end

function local_void_fn()
% a NAMED void function: nargout(@local_void_fn) == 0, which the wrapper rejects up front
end

function test_try_call_returns_exact_recoverable_category(testCase)
[v, completed, id] = invz_try_solver_call(@() local_thrower('invz:degenerateDoublet'));
verifyEmpty(testCase, v);
verifyFalse(testCase, completed);
verifyEqual(testCase, id, 'invz:degenerateDoublet');
end

function test_try_call_rethrows_fatal_and_returns_success(testCase)
verifyError(testCase, @() invz_try_solver_call(@() local_thrower('invz:staticMedium')), ...
    'invz:staticMedium');
[v, completed, id] = invz_try_solver_call(@() 42);
verifyEqual(testCase, v, 42);
verifyTrue(testCase, completed);
verifyEqual(testCase, id, '');
end

% The recoverable id must survive the wrapper EXACTLY, never collapse to a generic category.
function test_try_call_preserves_the_other_whitelisted_id(testCase)
[~, completed, id] = invz_try_solver_call(@() local_thrower('invz:orderedPhase'));
verifyFalse(testCase, completed);
verifyEqual(testCase, id, 'invz:orderedPhase');
end

% A void fn cannot deliver its own identifier through a one-output call, so it is rejected up
% front as a contract violation rather than silently masking a recoverable signal as
% MATLAB:maxlhs. Detectable only for NAMED void functions (nargout 0); the anonymous void form
% reports nargout -1 and is documented, not detected.
function test_void_fn_is_rejected_as_a_contract_violation(testCase)
verifyError(testCase, @() invz_try_solver_call(@local_void_fn), 'invz:solverCall');
end

function test_non_function_handle_is_rejected(testCase)
verifyError(testCase, @() invz_try_solver_call(42), 'invz:solverCall');
verifyError(testCase, @() invz_try_solver_call('not a handle'), 'invz:solverCall');
end

% A handle to a name that is not on the path CONSTRUCTS fine and only fails when nargout probes
% it. That probe must not leak a raw MATLAB:narginout:* identifier out of this boundary -- doing
% so would be the same un-namespaced escape the void-fn guard exists to prevent.
function test_unresolvable_named_handle_is_a_wiring_error(testCase)
verifyError(testCase, @() invz_try_solver_call(@invz_no_such_function_xyz123), ...
            'invz:solverCall');
verifyError(testCase, @() invz_try_solver_call(str2func('invz_no_such_function_xyz123')), ...
            'invz:solverCall');
end
