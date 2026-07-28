function trace = invzp_trace_fixed_h_adaptive(ctx,handoff,hTarget,spec)
%INVZP_TRACE_FIXED_H_ADAPTIVE Bounded two-state continuation in physical h.
%
%   This driver adapts only the proposed h spacing. Every attempted point is
%   solved and audited by INVZP_TRACE_FIXED_H_SEQUENCE using its certified
%   two-state handoff. A point is promoted only if that trace accepts it,
%   its scaled state correction passes the frozen predictor tube, and its r
%   value passes a separate secant-prediction tube. Failure shrinks the
%   spacing and retries from the last accepted pair; it never substitutes
%   the failed root. Growth requires a caller-frozen clean-success streak.

invalidId = 'invzp:AdaptiveFixedH:InvalidInput';
[handoff,hTarget,spec] = validateInputs( ...
    ctx,handoff,hTarget,spec,invalidId);

previous = handoff.previous;
current = handoff.current;
direction = sign(hTarget-current.h);
step = min(spec.step_initial,abs(hTarget-current.h));
streak = 0;
attempts = repmat(blankAttempt(),1,0);
points = repmat(blankPoint(),1,0);
status = 'ok';
detail = 'target reached';

while direction*(hTarget-current.h) > 0
    if numel(attempts) >= spec.max_attempts
        status = 'budget_exhausted';
        detail = 'max_attempts reached';
        break
    elseif numel(points) >= spec.max_points
        status = 'budget_exhausted';
        detail = 'max_points reached';
        break
    end

    remaining = abs(hTarget-current.h);
    stepTrial = min(step,remaining);
    hTrial = current.h+direction*stepTrial;
    if hTrial == current.h
        attempt = blankAttempt();
        attempt.index = numel(attempts)+1;
        attempt.h = hTrial;
        attempt.step = stepTrial;
        attempt.status = 'nonrepresentable_step';
        attempt.reason = 'the proposed h step does not change the floating-point field';
        attempts(end+1) = attempt; %#ok<AGROW>
        status = 'step_collapse';
        detail = attempt.reason;
        break
    end
    ratio = (hTrial-current.h)/(current.h-previous.h);
    rPredictor = current.r+ratio*(current.r-previous.r);

    fixedSpec = spec.fixed;
    fixedSpec.handoff = struct( ...
        'previous_state',previous.state, ...
        'previous_h',previous.h, ...
        'current_h',current.h);
    one = invzp_trace_fixed_h_sequence( ...
        ctx,hTrial,current.state,fixedSpec);

    attempt = blankAttempt();
    attempt.index = numel(attempts)+1;
    attempt.h = hTrial;
    attempt.step = stepTrial;
    attempt.fixed_status = one.status;
    attempt.fixed_detail = one.status_detail;
    attempt.r_predictor = rPredictor;
    if isscalar(one.records)
        attempt.fixed_record = one.records;
    end

    accepted = strcmp(one.status,'ok') && isscalar(one.accepted);
    if accepted
        record = one.accepted(1);
        attempt.state_predictor_distance = record.predictor_distance;
        if ~isfield(record.info,'r') || ...
                ~isnumeric(record.info.r) || ~isreal(record.info.r) || ...
                ~isscalar(record.info.r) || ~isfinite(record.info.r)
            accepted = false;
            attempt.status = 'invalid_fixed_output';
            attempt.reason = 'accepted fixed-h record has no finite real scalar info.r';
        else
            attempt.r = record.info.r;
            attempt.r_predictor_distance = abs(attempt.r-rPredictor);
            accepted = attempt.r_predictor_distance <= ...
                spec.r_predictor_tube_max;
            if accepted
                attempt.status = 'accepted';
                attempt.reason = 'accepted';
                point = blankPoint();
                point.index = numel(points)+1;
                point.h = hTrial;
                point.step = stepTrial;
                point.r_predictor = rPredictor;
                point.r_predictor_distance = attempt.r_predictor_distance;
                point.record = record;
                points(end+1) = point; %#ok<AGROW>
                previous = current;
                current = struct( ...
                    'h',hTrial, ...
                    'state',record.returned_state, ...
                    'r',record.info.r);
                streak = streak+1;
                if streak >= spec.successes_before_grow && ...
                        record.predictor_distance <= spec.grow_state_distance && ...
                        attempt.r_predictor_distance <= spec.grow_r_distance
                    step = min(spec.step_max,step*spec.step_grow);
                    streak = 0;
                end
            else
                attempt.status = 'r_predictor_rejected';
                attempt.reason = sprintf( ...
                    'r predictor distance %.17g exceeds %.17g', ...
                    attempt.r_predictor_distance,spec.r_predictor_tube_max);
            end
        end
    else
        attempt.status = 'fixed_h_rejected';
        attempt.reason = one.status_detail;
    end
    attempts(end+1) = attempt; %#ok<AGROW>

    if ~accepted
        streak = 0;
        nextStep = stepTrial*spec.step_shrink;
        if nextStep < spec.step_min
            status = 'step_collapse';
            detail = attempt.reason;
            break
        end
        step = max(spec.step_min,nextStep);
    end
end

trace = struct( ...
    'schema','invzp_fixed_h_adaptive/v1', ...
    'status',status, ...
    'status_detail',detail, ...
    'target_h',hTarget, ...
    'spec',spec, ...
    'initial_handoff',handoff, ...
    'attempts',attempts, ...
    'points',points, ...
    'last_pair',struct('previous',previous,'current',current));
end

function [handoff,hTarget,spec] = validateInputs( ...
        ctx,handoff,hTarget,spec,invalidId)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1') || ~isfield(ctx,'wn')
    error(invalidId,'ctx must come from invz_ordered_node_context.');
end
if ~isstruct(handoff) || ~isscalar(handoff) || ...
        ~all(isfield(handoff,{'previous','current'}))
    error(invalidId,'handoff requires scalar previous and current fields.');
end
handoff.previous = validatePoint(handoff.previous,numel(ctx.wn), ...
    'handoff.previous',invalidId);
handoff.current = validatePoint(handoff.current,numel(ctx.wn), ...
    'handoff.current',invalidId);
if handoff.previous.h == handoff.current.h
    error(invalidId,'The two handoff fields must be distinct.');
end
validateFiniteScalar(hTarget,'hTarget',invalidId);
if (hTarget-handoff.current.h)*(handoff.current.h-handoff.previous.h) < 0
    error(invalidId,'hTarget must continue the oriented handoff direction.');
end
if ~isstruct(spec) || ~isscalar(spec)
    error(invalidId,'spec must be a scalar struct.');
end
required = {'fixed','step_initial','step_min','step_max', ...
    'step_grow','step_shrink','successes_before_grow', ...
    'grow_state_distance','grow_r_distance','r_predictor_tube_max', ...
    'max_attempts','max_points'};
for k = 1:numel(required)
    if ~isfield(spec,required{k})
        error(invalidId,'spec.%s is required.',required{k});
    end
end
if ~isstruct(spec.fixed) || ~isscalar(spec.fixed)
    error(invalidId,'spec.fixed must be a scalar fixed-h trace specification.');
end
positiveNames = {'step_initial','step_min','step_max','step_grow', ...
    'step_shrink','grow_state_distance','grow_r_distance', ...
    'r_predictor_tube_max'};
for k = 1:numel(positiveNames)
    positive(spec.(positiveNames{k}),positiveNames{k},invalidId);
end
if spec.step_min > spec.step_initial || spec.step_initial > spec.step_max
    error(invalidId,'Require step_min <= step_initial <= step_max.');
end
if spec.step_grow <= 1 || spec.step_shrink >= 1
    error(invalidId,'Require step_grow > 1 and step_shrink < 1.');
end
integerNames = {'successes_before_grow','max_attempts','max_points'};
for k = 1:numel(integerNames)
    positiveInteger(spec.(integerNames{k}),integerNames{k},invalidId);
end
end

function point = validatePoint(point,nw,name,invalidId)
if ~isstruct(point) || ~isscalar(point) || ...
        ~all(isfield(point,{'h','state','r'}))
    error(invalidId,'%s requires h, state, and r.',name);
end
validateFiniteScalar(point.h,[name '.h'],invalidId);
validateFiniteScalar(point.r,[name '.r'],invalidId);
state = point.state;
if ~isstruct(state) || ~isscalar(state) || ...
        ~isfield(state,'Sigma') || ~isfield(state,'K0s') || ...
        ~isnumeric(state.Sigma) || ~isreal(state.Sigma) || ...
        ~isvector(state.Sigma) || numel(state.Sigma) ~= nw || ...
        any(~isfinite(state.Sigma),'all') || ...
        ~isnumeric(state.K0s) || ~isreal(state.K0s) || ...
        ~isscalar(state.K0s) || ~isfinite(state.K0s)
    error(invalidId,'%s.state has invalid Sigma or K0s.',name);
end
point.state = struct('Sigma',state.Sigma(:),'K0s',state.K0s);
end

function attempt = blankAttempt()
attempt = struct( ...
    'index',NaN, ...
    'h',NaN, ...
    'step',NaN, ...
    'status','not_run', ...
    'reason','', ...
    'fixed_status','not_run', ...
    'fixed_detail','', ...
    'state_predictor_distance',NaN, ...
    'r',NaN, ...
    'r_predictor',NaN, ...
    'r_predictor_distance',NaN, ...
    'fixed_record',struct());
end

function point = blankPoint()
point = struct( ...
    'index',NaN, ...
    'h',NaN, ...
    'step',NaN, ...
    'r_predictor',NaN, ...
    'r_predictor_distance',NaN, ...
    'record',struct());
end

function positive(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error(invalidId,'spec.%s must be a finite positive scalar.',name);
end
end

function positiveInteger(value,name,invalidId)
positive(value,name,invalidId);
if value ~= floor(value)
    error(invalidId,'spec.%s must be an integer.',name);
end
end

function validateFiniteScalar(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ~isfinite(value)
    error(invalidId,'%s must be a finite real scalar.',name);
end
end
