function trace = invzp_trace_fixed_w_sequence( ...
        ctx,wValues,initialPoint,spec)
%INVZP_TRACE_FIXED_W_SEQUENCE Trace accepted roots in w=z-K0.
%
%   TRACE = INVZP_TRACE_FIXED_W_SEQUENCE(CTX,WVALUES,INITIALPOINT,SPEC)
%   corrects the unchanged ordered equations at a caller-frozen strictly
%   monotone sequence of target values of
%
%       w = z-K0.
%
%   The square corrector appends (w-w_target)/w_scale to the ordered
%   residual. It uses the analytic state gradient of w and the same
%   Richardson field column as the retained arclength problem. This helper
%   does not adapt targets, infer an endpoint, select a branch, or claim a
%   complete Jensen section.

invalidId = 'invzp:FixedWTrace:InvalidInput';
[wValues,initialPoint,cfg] = validateInputs( ...
    ctx,wValues,initialPoint,spec,invalidId);
arc = invzp_ordered_arclength_problem(ctx,cfg.problem_options);
npoint = numel(wValues);
ncoordinate = numel(ctx.wn)+2;
records = repmat(blankRecord(ncoordinate),1,npoint);
acceptedY = zeros(ncoordinate,0);
stopped = false;
stopReason = '';

for k = 1:npoint
    record = blankRecord(ncoordinate);
    record.index = k;
    record.target_w = wValues(k);
    if stopped
        record.reason = stopReason;
        records(k) = record;
        continue
    end

    if k == 1
        predictor = arc.pack(initialPoint.state,initialPoint.h);
        predictorKind = 'initial_point';
    elseif k == 2
        predictor = acceptedY(:,end);
        predictorKind = 'copy';
    else
        ratio = (wValues(k)-wValues(k-1))/ ...
            (wValues(k-1)-wValues(k-2));
        predictor = acceptedY(:,end)+ ...
            ratio*(acceptedY(:,end)-acceptedY(:,end-1));
        predictorKind = 'secant';
    end
    record.predictor_kind = predictorKind;
    record.predictor_y = predictor;

    [corrected,callback,corrector] = correctTarget( ...
        arc,predictor,wValues(k),cfg);
    record.corrected_y = corrected;
    record.corrector = corrector;
    if ~corrector.accepted
        record.status = 'corrector_failed';
        record.reason = corrector.reason;
        if isstruct(callback) && isscalar(callback)
            record.callback = callback;
        end
        records(k) = record;
        stopped = true;
        stopReason = record.reason;
        continue
    end

    correctionRms = norm(corrected-predictor)/sqrt(ncoordinate);
    record.predictor_correction_rms = correctionRms;
    if ~isfinite(correctionRms) || ...
            correctionRms > cfg.predictor_tube_rms_max
        record.status = 'predictor_tube_exceeded';
        record.reason = sprintf( ...
            'predictor correction %.17g exceeds %.17g', ...
            correctionRms,cfg.predictor_tube_rms_max);
        record.callback = callback;
        records(k) = record;
        stopped = true;
        stopReason = record.reason;
        continue
    end

    record.status = 'accepted';
    record.reason = 'accepted';
    record.path_accepted = true;
    record.h = callback.h;
    record.w = callback.w;
    record.state = callback.state;
    record.physics = callback.physics;
    record.event_margins = callback.event_margins;
    record.callback = callback;
    records(k) = record;
    acceptedY(:,end+1) = corrected; %#ok<AGROW>
end

acceptedCount = size(acceptedY,2);
trace = struct( ...
    'schema','invzp_fixed_w_sequence/v1', ...
    'status','ok', ...
    'status_detail',sprintf('%d target values accepted',npoint), ...
    'initial_point_id',cfg.initial_point_id, ...
    'w_values',wValues, ...
    'spec',cfg, ...
    'records',records, ...
    'accepted_indices',(1:acceptedCount).', ...
    'accepted_y',acceptedY, ...
    'complete_candidate_claim',false);
if stopped
    trace.status = 'stopped';
    trace.status_detail = stopReason;
end
end

function [y,callback,info] = correctTarget(arc,predictor,target,cfg)
y = predictor(:);
ncoordinate = numel(y);
callback = struct();
residualHistory = nan(1,cfg.max_iterations+1);
info = struct( ...
    'accepted',false,'reason','max_iterations', ...
    'iterations',0,'equation_evaluations',0, ...
    'residual_inf',Inf,'bordered_rcond',0, ...
    'launch_used',false,'launch_correction_rms',NaN, ...
    'residual_history',zeros(1,0), ...
    'audit',struct());

for iteration = 0:cfg.max_iterations
    if info.equation_evaluations >= cfg.max_evaluations
        info.reason = 'evaluation_budget';
        break
    end
    [F,A,callback] = fixedWEquations(arc,y,target,cfg.w_scale);
    info.equation_evaluations = info.equation_evaluations+1;
    residual = norm(F,Inf);
    residualHistory(iteration+1) = residual;
    info.iterations = iteration;
    info.residual_inf = residual;
    if all(isfinite(A),'all')
        info.bordered_rcond = rcond(A);
    end

    if any(~isfinite(F),'all') || any(~isfinite(A),'all') || ...
            ~isfinite(info.bordered_rcond)
        info.reason = 'nonfinite_system';
        break
    elseif info.bordered_rcond <= cfg.rcond_min
        % A small residual is not sufficient evidence that the fixed-w
        % section is locally regular. Apply the rank gate before accepting
        % a residual-passing initial/corrected point.
        info.reason = 'rank_gate';
        break
    elseif residual <= cfg.tol_residual
        eventReason = arc.event(callback);
        [auditAccepted,auditReason,auditPayload] = ...
            arc.audit(y,callback);
        info.audit = auditPayload;
        if ~isempty(eventReason)
            info.reason = ['event:',eventReason];
        elseif ~auditAccepted
            info.reason = auditReason;
        else
            info.accepted = true;
            info.reason = 'accepted';
        end
        break
    elseif iteration == cfg.max_iterations
        break
    end

    step = -(A\F);
    stepRms = norm(step)/sqrt(ncoordinate);
    if iteration == 0 && isfinite(stepRms) && ...
            stepRms <= cfg.launch_correction_rms_max
        if info.equation_evaluations >= cfg.max_evaluations
            info.reason = 'evaluation_budget';
            break
        end
        trial = y+step;
        Ft = fixedWEquations(arc,trial,target,cfg.w_scale);
        info.equation_evaluations = info.equation_evaluations+1;
        if all(isfinite(Ft),'all')
            y = trial;
            info.launch_used = true;
            info.launch_correction_rms = stepRms;
            continue
        end
    end

    acceptedTrial = false;
    damping = 1;
    for line = 1:cfg.max_linesearch
        if info.equation_evaluations >= cfg.max_evaluations
            info.reason = 'evaluation_budget';
            break
        end
        trial = y+damping*step;
        Ft = fixedWEquations(arc,trial,target,cfg.w_scale);
        info.equation_evaluations = info.equation_evaluations+1;
        if all(isfinite(Ft),'all') && norm(Ft,Inf) < residual
            y = trial;
            acceptedTrial = true;
            break
        end
        damping = damping/2;
    end
    if ~acceptedTrial
        if ~strcmp(info.reason,'evaluation_budget')
            info.reason = 'line_search';
        end
        break
    end
end
info.residual_history = residualHistory(isfinite(residualHistory));
end

function [F,A,record] = fixedWEquations(arc,y,target,wScale)
[R,J,record] = arc.equations(y);
if ~isfield(record,'w') || ~isfield(record,'w_gradient_scaled')
    F = nan(size(J,2),1);
    A = nan(size(J,2));
    return
end
F = [R;(record.w-target)/wScale];
A = [J;record.w_gradient_scaled/wScale];
end

function [wValues,point,cfg] = validateInputs( ...
        ctx,wValues,point,spec,invalidId)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1') || ...
        ~isfield(ctx,'wn')
    error(invalidId,'ctx must come from invz_ordered_node_context.');
end
if ~isa(wValues,'double') || ~isreal(wValues) || ...
        ~isvector(wValues) || isempty(wValues) || ...
        any(~isfinite(wValues),'all')
    error(invalidId,'wValues must be a nonempty finite real double vector.');
end
wValues = wValues(:).';
if numel(wValues) > 1 && ...
        ~(all(diff(wValues) > 0) || all(diff(wValues) < 0))
    error(invalidId,'wValues must be strictly monotone.');
end
if ~isstruct(point) || ~isscalar(point) || ...
        ~all(isfield(point,{'state','h'})) || ...
        ~isstruct(point.state) || ~isscalar(point.state) || ...
        ~all(isfield(point.state,{'Sigma','K0s'})) || ...
        ~isnumeric(point.state.Sigma) || ~isreal(point.state.Sigma) || ...
        ~isvector(point.state.Sigma) || ...
        numel(point.state.Sigma) ~= numel(ctx.wn) || ...
        any(~isfinite(point.state.Sigma),'all') || ...
        ~isnumeric(point.state.K0s) || ~isreal(point.state.K0s) || ...
        ~isscalar(point.state.K0s) || ~isfinite(point.state.K0s) || ...
        ~isnumeric(point.h) || ~isreal(point.h) || ...
        ~isscalar(point.h) || ~isfinite(point.h)
    error(invalidId, ...
        'initialPoint requires a finite state (Sigma,K0s) and scalar h.');
end
point.state = struct( ...
    'Sigma',point.state.Sigma(:),'K0s',point.state.K0s);

if ~isstruct(spec) || ~isscalar(spec)
    error(invalidId,'spec must be a scalar struct.');
end
required = {'schema','version','initial_point_id','problem_options', ...
    'w_scale','tol_residual','rcond_min','max_iterations', ...
    'max_linesearch','max_evaluations','launch_correction_rms_max', ...
    'predictor_tube_rms_max'};
for k = 1:numel(required)
    if ~isfield(spec,required{k})
        error(invalidId,'spec.%s is required.',required{k});
    end
end
if stringScalar(spec.schema,'spec.schema',invalidId) ~= ...
        "invzp_fixed_w_sequence_spec/v1"
    error(invalidId, ...
        'spec.schema must be invzp_fixed_w_sequence_spec/v1.');
end
cfg = struct( ...
    'schema','invzp_fixed_w_sequence_spec/v1', ...
    'version',stringScalar(spec.version,'spec.version',invalidId), ...
    'initial_point_id',stringScalar( ...
        spec.initial_point_id,'spec.initial_point_id',invalidId), ...
    'problem_options',spec.problem_options, ...
    'w_scale',positive(spec.w_scale,'spec.w_scale',invalidId), ...
    'tol_residual',positive( ...
        spec.tol_residual,'spec.tol_residual',invalidId), ...
    'rcond_min',positive(spec.rcond_min,'spec.rcond_min',invalidId), ...
    'max_iterations',positiveInteger( ...
        spec.max_iterations,'spec.max_iterations',invalidId), ...
    'max_linesearch',positiveInteger( ...
        spec.max_linesearch,'spec.max_linesearch',invalidId), ...
    'max_evaluations',positiveInteger( ...
        spec.max_evaluations,'spec.max_evaluations',invalidId), ...
    'launch_correction_rms_max',positive( ...
        spec.launch_correction_rms_max, ...
        'spec.launch_correction_rms_max',invalidId), ...
    'predictor_tube_rms_max',positive( ...
        spec.predictor_tube_rms_max, ...
        'spec.predictor_tube_rms_max',invalidId));
if ~isstruct(cfg.problem_options) || ~isscalar(cfg.problem_options)
    error(invalidId,'spec.problem_options must be a scalar struct.');
end
end

function value = positive(value,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error(invalidId,'%s must be a finite positive scalar.',name);
end
end

function value = positiveInteger(value,name,invalidId)
value = positive(value,name,invalidId);
if value ~= floor(value)
    error(invalidId,'%s must be a positive integer.',name);
end
end

function value = stringScalar(value,name,invalidId)
if ischar(value) && isrow(value)
    value = string(value);
elseif isstring(value) && isscalar(value) && ~ismissing(value)
    value = string(value);
else
    error(invalidId,'%s must be a nonmissing text scalar.',name);
end
if strlength(value) == 0
    error(invalidId,'%s must not be empty.',name);
end
end

function record = blankRecord(ncoordinate)
record = struct( ...
    'index',NaN,'target_w',NaN,'status','not_run_after_stop', ...
    'reason','','predictor_kind','not_computed', ...
    'predictor_y',nan(ncoordinate,1), ...
    'corrected_y',nan(ncoordinate,1), ...
    'predictor_correction_rms',NaN, ...
    'path_accepted',false,'h',NaN,'w',NaN, ...
    'state',struct('Sigma',[],'K0s',NaN), ...
    'physics',struct(),'event_margins',struct(), ...
    'callback',struct(),'corrector',struct());
end
