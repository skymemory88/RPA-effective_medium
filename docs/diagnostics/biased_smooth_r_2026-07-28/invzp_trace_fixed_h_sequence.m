function trace = invzp_trace_fixed_h_sequence(ctx,hValues,initialState,spec)
%INVZP_TRACE_FIXED_H_SEQUENCE Predictor-controlled natural-parameter trace.
%
%   HVALUES is a caller-frozen strictly monotone sequence. This helper
%   neither inserts, sorts, adapts, nor skips fields. Each promoted point
%   must pass Newton's independent A-D audit and a frozen full-state
%   predictor-tube gate. A failure stops the path and leaves explicit
%   not-run records for the rest of the submitted sequence.

invalidId = 'invzp:FixedHTrace:InvalidInput';
[hValues,spec] = validateInputs(ctx,hValues,initialState,spec,invalidId);
npoint = numel(hValues);
records = repmat(blankRecord(),1,npoint);
acceptedIndices = zeros(0,1);
stopped = false;
stopDetail = '';

for k = 1:npoint
    record = blankRecord();
    record.index = k;
    record.h = hValues(k);
    if stopped
        record.status = 'not_run_after_stop';
        record.reason = stopDetail;
        records(k) = record;
        continue
    end

    if k == 1
        predictor = normalizeState(initialState);
        predictorKind = 'initial_seed';
        record.seed_id = spec.initial_state_id;
    elseif isscalar(acceptedIndices)
        predictor = normalizeState(records(acceptedIndices(end)).returned_state);
        predictorKind = 'copy';
    else
        last = acceptedIndices(end);
        previous = acceptedIndices(end-1);
        ratio = (hValues(k)-hValues(last))/ ...
            (hValues(last)-hValues(previous));
        currentState = records(last).returned_state;
        previousState = records(previous).returned_state;
        predictor = struct( ...
            'Sigma',currentState.Sigma(:)+ratio* ...
                (currentState.Sigma(:)-previousState.Sigma(:)), ...
            'K0s',currentState.K0s+ratio* ...
                (currentState.K0s-previousState.K0s));
        predictorKind = 'secant';
    end
    record.predictor_kind = predictorKind;
    record.predictor_state = predictor;
    if any(~isfinite(predictor.Sigma),'all') || ~isfinite(predictor.K0s)
        record.status = 'predictor_nonfinite';
        record.reason = 'the declared predictor produced a nonfinite state';
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    end

    try
        [node,nodeMeta] = invz_ordered_make_node(ctx,hValues(k));
    catch exception
        record.status = 'node_construction_failed';
        record.reason = ['node_exception:',exception.identifier];
        record.node_meta = exceptionInfo(exception);
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    end
    record.node_meta = nodeMeta;
    if isempty(node)
        record.status = 'node_unavailable';
        record.reason = nodeMeta.status;
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    end

    try
        [state,info] = invz_ordered_node_newton(node,predictor,spec.newton);
    catch exception
        record.status = 'solver_failed';
        record.reason = ['solver_exception:',exception.identifier];
        record.returned_state = predictor;
        record.info = struct('accepted',false,'reason','solver_exception', ...
            'exception_identifier',exception.identifier, ...
            'exception_message',exception.message);
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    end
    record.returned_state = state;
    record.info = info;
    if ~validSolverOutput(state,info,numel(ctx.wn))
        record.status = 'solver_failed';
        record.reason = 'invalid_solver_output';
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    end
    record.fixed_h_accepted = info.accepted && info.audit.accepted;
    record.predictor_distance = stateDistance( ...
        state,predictor,spec.scale.sigma,spec.scale.k0);
    if ~record.fixed_h_accepted
        record.status = 'solver_failed';
        record.reason = char(string(info.reason));
        if isempty(record.reason), record.reason = 'fixed_h_audit_failed'; end
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    elseif ~isfinite(record.predictor_distance)
        record.status = 'solver_failed';
        record.reason = 'nonfinite_predictor_distance';
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    elseif record.predictor_distance > spec.predictor_tube_max
        record.status = 'predictor_tube_exceeded';
        record.reason = sprintf( ...
            'predictor distance %.17g exceeds %.17g', ...
            record.predictor_distance,spec.predictor_tube_max);
        records(k) = record;
        stopped = true;
        stopDetail = record.reason;
        continue
    end

    record.status = 'accepted';
    record.reason = 'accepted';
    record.path_accepted = true;
    records(k) = record;
    acceptedIndices(end+1,1) = k; %#ok<AGROW>
end

trace = struct( ...
    'schema','invzp_fixed_h_sequence/v1', ...
    'status','ok', ...
    'status_detail',sprintf('%d scheduled fields accepted',npoint), ...
    'initial_state_id',spec.initial_state_id, ...
    'h_values',hValues, ...
    'spec',spec, ...
    'records',records, ...
    'accepted_indices',acceptedIndices, ...
    'accepted',records(acceptedIndices), ...
    'handoff',struct('last_accepted',struct([]), ...
        'previous_accepted',struct([])));
if stopped
    trace.status = 'stopped';
    trace.status_detail = stopDetail;
end
if ~isempty(acceptedIndices)
    trace.handoff.last_accepted = records(acceptedIndices(end));
end
if numel(acceptedIndices) >= 2
    trace.handoff.previous_accepted = records(acceptedIndices(end-1));
end
end

function [hValues,spec] = validateInputs(ctx,hValues,initialState,spec,invalidId)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1') || ~isfield(ctx,'wn')
    error(invalidId,'ctx must come from invz_ordered_node_context.');
end
if ~isnumeric(hValues) || ~isreal(hValues) || ~isvector(hValues) || ...
        isempty(hValues) || any(~isfinite(hValues),'all')
    error(invalidId,'hValues must be a nonempty finite real vector.');
end
hValues = hValues(:).';
dh = diff(hValues);
if ~isempty(dh) && ~(all(dh > 0) || all(dh < 0))
    error(invalidId,'hValues must be strictly monotone in traversal order.');
end
validateState(initialState,numel(ctx.wn),invalidId);
if ~isstruct(spec) || ~isscalar(spec)
    error(invalidId,'spec must be a scalar struct.');
end
required = {'initial_state_id','newton','scale','predictor_tube_max'};
for k = 1:numel(required)
    if ~isfield(spec,required{k})
        error(invalidId,'spec.%s is required.',required{k});
    end
end
spec.initial_state_id = normalizeId(spec.initial_state_id,invalidId);
if ~isstruct(spec.newton) || ~isscalar(spec.newton)
    error(invalidId,'spec.newton must be a scalar struct.');
end
if ~isstruct(spec.scale) || ~isscalar(spec.scale) || ...
        ~isfield(spec.scale,'sigma') || ~isfield(spec.scale,'k0')
    error(invalidId,'spec.scale requires sigma and k0.');
end
sigmaScale = spec.scale.sigma;
nw = numel(ctx.wn);
if isscalar(sigmaScale), sigmaScale = repmat(sigmaScale,nw,1); end
if ~isnumeric(sigmaScale) || ~isreal(sigmaScale) || ...
        ~isvector(sigmaScale) || numel(sigmaScale) ~= nw || ...
        any(~isfinite(sigmaScale),'all') || any(sigmaScale <= 0)
    error(invalidId,'spec.scale.sigma must be a positive scalar or nw-vector.');
end
spec.scale.sigma = sigmaScale(:);
spec.scale.k0 = positive(spec.scale.k0,'scale.k0',invalidId);
spec.predictor_tube_max = positive( ...
    spec.predictor_tube_max,'predictor_tube_max',invalidId);
end

function validateState(state,nw,invalidId)
if ~isstruct(state) || ~isscalar(state) || ...
        ~isfield(state,'Sigma') || ~isfield(state,'K0s') || ...
        ~isnumeric(state.Sigma) || ~isreal(state.Sigma) || ...
        ~isvector(state.Sigma) || numel(state.Sigma) ~= nw || ...
        any(~isfinite(state.Sigma),'all') || ...
        ~isnumeric(state.K0s) || ~isreal(state.K0s) || ...
        ~isscalar(state.K0s) || ~isfinite(state.K0s)
    error(invalidId, ...
        'initialState requires finite real nw-vector Sigma and scalar K0s.');
end
end

function state = normalizeState(state)
state = struct('Sigma',state.Sigma(:),'K0s',state.K0s);
end

function distance = stateDistance(state,predictor,sigmaScale,k0Scale)
if ~validStateCoordinates(state,numel(sigmaScale))
    distance = NaN;
    return
end
dsigma = (state.Sigma(:)-predictor.Sigma(:))./sigmaScale;
dK0 = (state.K0s-predictor.K0s)/k0Scale;
distance = sqrt(mean(abs(dsigma).^2)+abs(dK0).^2);
end

function valid = validSolverOutput(state,info,nw)
baseValid = validStateCoordinates(state,nw) && ...
    isstruct(info) && isscalar(info) && ...
    isfield(info,'accepted') && islogical(info.accepted) && ...
    isscalar(info.accepted) && ...
    isfield(info,'audit') && ...
    isfield(info,'reason') && validText(info.reason);
if ~baseValid
    valid = false;
    return
end
auditValid = isstruct(info.audit) && isscalar(info.audit) && ...
    isfield(info.audit,'accepted') && islogical(info.audit.accepted) && ...
    isscalar(info.audit.accepted);
valid = auditValid || (~info.accepted && isempty(info.audit));
end

function valid = validStateCoordinates(state,nw)
valid = isstruct(state) && isscalar(state) && ...
    isfield(state,'Sigma') && isnumeric(state.Sigma) && ...
    isreal(state.Sigma) && isvector(state.Sigma) && ...
    numel(state.Sigma) == nw && all(isfinite(state.Sigma),'all') && ...
    isfield(state,'K0s') && isnumeric(state.K0s) && ...
    isreal(state.K0s) && isscalar(state.K0s) && isfinite(state.K0s);
end

function valid = validText(value)
valid = (ischar(value) && (isempty(value) || isrow(value))) || ...
    (isstring(value) && isscalar(value) && ~ismissing(value));
end

function info = exceptionInfo(exception)
info = struct('status','exception', ...
    'exception_identifier',exception.identifier, ...
    'exception_message',exception.message);
end

function record = blankRecord()
state = struct('Sigma',[],'K',[],'lam',[],'K0s',NaN);
record = struct( ...
    'index',NaN, ...
    'h',NaN, ...
    'seed_id',"", ...
    'status','not_run_after_stop', ...
    'predictor_kind','not_computed', ...
    'predictor_state',struct('Sigma',[],'K0s',NaN), ...
    'predictor_distance',NaN, ...
    'node_meta',struct(), ...
    'returned_state',state, ...
    'info',struct(), ...
    'fixed_h_accepted',false, ...
    'path_accepted',false, ...
    'reason','');
end

function id = normalizeId(value,invalidId)
if ischar(value)
    valid = isrow(value) && ~isempty(value);
elseif isstring(value)
    valid = isscalar(value) && ~ismissing(value) && strlength(value) > 0;
else
    valid = false;
end
if ~valid
    error(invalidId, ...
        'spec.initial_state_id must be a nonempty character row or string scalar.');
end
id = string(value);
end

function value = positive(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error(invalidId,'spec.%s must be a finite positive scalar.',name);
end
end
