function trace = invzp_trace_pseudo_arclength(problem, y0, opts)
%INVZP_TRACE_PSEUDO_ARCLENGTH Trace one regular residual curve through folds.
%
%   TRACE = INVZP_TRACE_PSEUDO_ARCLENGTH(PROBLEM,Y0,OPTS) follows R(Y)=0,
%   where R has one fewer component than the scaled coordinate vector
%   Y=[SCALED_STATE;SCALED_PARAMETER]. The continuation parameter is
%   required to be the last coordinate.
%   PROBLEM.equations(Y) must return [R,J,DIAG], with J=dR/dY.
%   PROBLEM.event(DIAG) optionally returns an empty reason for a traversable
%   point or a nonempty reason for a hard event. PROBLEM.audit(Y,DIAG)
%   optionally returns [ACCEPTED,REASON,PAYLOAD] at corrected roots.
%
%   This diagnostic primitive traces one connected curve. It does not find
%   initial roots, enumerate disconnected components, split the curve into
%   single-valued sections, or select a physical path.

if nargin < 3 || isempty(opts), opts = struct(); end
invalidId = 'invzp:ArcLength:InvalidInput';
validateProblem(problem,invalidId);
if ~isnumeric(y0) || ~isreal(y0) || ~isvector(y0) || ...
        numel(y0) < 2 || any(~isfinite(y0),'all')
    fail(invalidId,'y0 must be a finite real vector with at least two entries.');
end
y0 = y0(:);
nvar = numel(y0);
nres = nvar-1;

cfg = resolveOptions(opts,nvar,invalidId);
[R0,J0,diag0] = evaluate(problem,y0,nres,nvar,invalidId);
event0 = equationEvent(problem,R0,J0,diag0);
if ~isempty(event0)
    fail(invalidId,'The initial root lies on event "%s".',event0);
end
if norm(R0,Inf) > cfg.tol_residual
    fail(invalidId,'y0 residual %.6g exceeds opts.tol_residual %.6g.', ...
        norm(R0,Inf),cfg.tol_residual);
end
[accepted0,auditReason0,payload0] = auditPoint(problem,y0,diag0);
if ~accepted0
    fail(invalidId,'The initial root failed its audit: %s.',auditReason0);
end

if isempty(cfg.initial_tangent)
    Ju = J0(:,1:nres);
    Jh = J0(:,end);
    if rcond(Ju) <= cfg.rcond_min
        fail(invalidId, ...
            'Initial fixed-parameter Jacobian is singular; supply opts.initial_tangent.');
    end
    tangent = [-(Ju\Jh)*cfg.initial_direction; cfg.initial_direction];
else
    tangent = cfg.initial_tangent;
    nullResidual = norm(J0*tangent,Inf);
    if nullResidual > cfg.tangent_residual_tol
        fail(invalidId, ...
            'opts.initial_tangent has null residual %.6g above %.6g.', ...
            nullResidual,cfg.tangent_residual_tol);
    end
end
tangent = tangent/norm(tangent);
if tangent(end)*cfg.initial_direction < 0
    tangent = -tangent;
end

trace = struct( ...
    'schema','invzp_pseudo_arclength/v1', ...
    'status','running', ...
    'status_detail','', ...
    'y',y0, ...
    'tangent',tangent, ...
    'records',initialRecord(y0,R0,diag0,payload0,tangent), ...
    'attempts',repmat(blankAttempt(),1,0), ...
    'options',cfg);

ds = cfg.step_initial;
for istep = 1:cfg.max_steps
    acceptedStep = false;
    lastReason = 'step_retry_exhausted';
    iattempt = 0;
    attemptBudgetHit = false;
    while ds >= cfg.step_min
        if iattempt >= cfg.max_step_attempts
            attemptBudgetHit = true;
            break
        end
        iattempt = iattempt+1;
        yPred = trace.y(:,end)+ds*trace.tangent(:,end);
        [ok,yNew,Rnew,Jnew,diagNew,corr] = corrector( ...
            problem,yPred,trace.tangent(:,end),cfg,nres,nvar,invalidId);
        if ~ok
            lastReason = corr.reason;
            trace.attempts(end+1) = attemptRecord( ...
                istep,iattempt,ds,yPred,yNew,lastReason,corr,NaN);
            ds = ds*cfg.step_shrink;
            continue
        end

        correctionRatio = norm(yNew-yPred)/ds;
        if correctionRatio > cfg.correction_ratio_max
            lastReason = 'branch_correction_too_large';
            trace.attempts(end+1) = attemptRecord( ...
                istep,iattempt,ds,yPred,yNew,lastReason,corr,correctionRatio);
            ds = ds*cfg.step_shrink;
            continue
        end

        [accepted,auditReason,payload] = auditPoint(problem,yNew,diagNew);
        if ~accepted
            lastReason = ['audit:',auditReason];
            trace.attempts(end+1) = attemptRecord( ...
                istep,iattempt,ds,yPred,yNew,lastReason,corr,correctionRatio);
            ds = ds*cfg.step_shrink;
            continue
        end

        B = [Jnew; trace.tangent(:,end).'];
        borderedRcond = rcond(B);
        if ~isfinite(borderedRcond) || borderedRcond <= cfg.rcond_min
            lastReason = 'augmented_rank_loss';
            trace.attempts(end+1) = attemptRecord( ...
                istep,iattempt,ds,yPred,yNew,lastReason,corr,correctionRatio);
            ds = ds*cfg.step_shrink;
            continue
        end
        tangentNew = B\[zeros(nres,1);1];
        tangentNew = tangentNew/norm(tangentNew);
        if dot(tangentNew,trace.tangent(:,end)) < 0
            tangentNew = -tangentNew;
        end
        tangentOverlap = dot(tangentNew,trace.tangent(:,end));
        if tangentOverlap < cfg.tangent_overlap_min
            lastReason = 'tangent_overlap';
            trace.attempts(end+1) = attemptRecord( ...
                istep,iattempt,ds,yPred,yNew,lastReason,corr,correctionRatio);
            ds = ds*cfg.step_shrink;
            continue
        end

        acceptedStep = true;
        trace.attempts(end+1) = attemptRecord( ...
            istep,iattempt,ds,yPred,yNew,'accepted',corr,correctionRatio);
        trace.y(:,end+1) = yNew;
        trace.tangent(:,end+1) = tangentNew;
        trace.records(end+1) = makeRecord( ...
            istep,yNew,Rnew,diagNew,payload,ds,corr.iterations, ...
            corr.constraint_abs,correctionRatio,tangentOverlap,borderedRcond);
        trace.records(end).tangent_parameter = tangentNew(end);

        if corr.iterations <= cfg.target_iterations && ...
                correctionRatio <= cfg.grow_correction_ratio
            ds = min(cfg.step_max,ds*cfg.step_grow);
        elseif corr.iterations > cfg.target_iterations+2
            ds = max(cfg.step_min,ds*cfg.step_shrink);
        end
        break
    end

    if ~acceptedStep
        attemptBudgetHit = attemptBudgetHit || ...
            iattempt >= cfg.max_step_attempts;
        if attemptBudgetHit
            trace.status = 'budget_exhausted';
            trace.status_detail = ['step_attempt_budget:',lastReason];
        elseif startsWith(lastReason,'event:')
            trace.status = 'event';
            trace.status_detail = lastReason;
        else
            trace.status = 'step_collapse';
            trace.status_detail = lastReason;
        end
        return
    end
end

trace.status = 'max_steps';
trace.status_detail = 'the configured accepted-step budget was reached';
end

function [ok,y,R,J,diag,corr] = corrector( ...
        problem,yPred,tangent,cfg,nres,nvar,invalidId)
y = yPred;
R = nan(nres,1);
J = nan(nres,nvar);
diag = struct();
corr = struct('reason','max_corrector','iterations',0, ...
    'constraint_abs',NaN,'equation_evaluations',0);
for iter = 1:cfg.max_corrector
    if corr.equation_evaluations >= cfg.max_evaluations_per_attempt
        corr.reason = 'evaluation_budget';
        ok = false;
        return
    end
    corr.equation_evaluations = corr.equation_evaluations+1;
    [R,J,diag] = evaluate(problem,y,nres,nvar,invalidId);
    event = equationEvent(problem,R,J,diag);
    constraint = tangent.'*(y-yPred);
    F = [R;constraint];
    corr.iterations = iter;
    corr.constraint_abs = abs(constraint);
    if isempty(event) && norm(R,Inf) <= cfg.tol_residual && ...
            abs(constraint) <= cfg.tol_constraint
        corr.reason = 'accepted';
        ok = true;
        return
    end
    if ~isempty(event)
        corr.reason = ['event:',event];
        ok = false;
        return
    end

    B = [J;tangent.'];
    borderedRcond = rcond(B);
    if ~isfinite(borderedRcond) || borderedRcond <= cfg.rcond_min
        corr.reason = 'augmented_rank_loss';
        ok = false;
        return
    end
    delta = -(B\F);
    if any(~isfinite(delta))
        corr.reason = 'nonfinite_step';
        ok = false;
        return
    end

    baseNorm = norm(F,2);
    alpha = 1;
    acceptedLine = false;
    for ils = 1:cfg.max_linesearch
        if corr.equation_evaluations >= cfg.max_evaluations_per_attempt
            corr.reason = 'evaluation_budget';
            ok = false;
            return
        end
        corr.equation_evaluations = corr.equation_evaluations+1;
        trial = y+alpha*delta;
        [Rt,Jt,diagt] = evaluate(problem,trial,nres,nvar,invalidId);
        eventt = equationEvent(problem,Rt,Jt,diagt);
        Ft = [Rt;tangent.'*(trial-yPred)];
        if isempty(eventt) && norm(Ft,2) < baseNorm
            y = trial;
            acceptedLine = true;
            break
        end
        alpha = alpha/2;
    end
    if ~acceptedLine
        corr.reason = 'line_search_failed';
        ok = false;
        return
    end
end
ok = false;
end

function [R,J,diag] = evaluate(problem,y,nres,nvar,invalidId)
try
    [R,J,diag] = problem.equations(y);
catch exception
    wrapped = MException('invzp:ArcLength:EquationFailure', ...
        'problem.equations failed: %s',exception.message);
    wrapped = addCause(wrapped,exception);
    throw(wrapped);
end
if ~isnumeric(R) || ~isreal(R) || ~isequal(size(R),[nres,1])
    fail(invalidId,'problem.equations must return a [%d x 1] real R.',nres);
end
if ~isnumeric(J) || ~isreal(J) || ~isequal(size(J),[nres,nvar])
    fail(invalidId, ...
        'problem.equations must return a [%d x %d] real J.',nres,nvar);
end
if ~isstruct(diag) || ~isscalar(diag)
    fail(invalidId,'problem.equations DIAG output must be a scalar struct.');
end
end

function reason = equationEvent(problem,R,J,diag)
reason = eventReason(problem,diag);
if ~isempty(reason)
    return
elseif any(~isfinite(R),'all')
    reason = 'nonfinite_residual';
elseif any(~isfinite(J),'all')
    reason = 'nonfinite_jacobian';
end
end

function reason = eventReason(problem,diag)
reason = '';
if ~isfield(problem,'event') || isempty(problem.event)
    return
end
reason = problem.event(diag);
if isstring(reason)
    if ~isscalar(reason) || ismissing(reason)
        error('invzp:ArcLength:InvalidInput', ...
            'problem.event must return a nonmissing string scalar or character row.');
    end
    reason = char(reason);
elseif ~(ischar(reason) && (isempty(reason) || isrow(reason)))
    error('invzp:ArcLength:InvalidInput', ...
        'problem.event must return a string scalar or character row.');
end
end

function [accepted,reason,payload] = auditPoint(problem,y,diag)
if ~isfield(problem,'audit') || isempty(problem.audit)
    accepted = true;
    reason = '';
    payload = struct();
    return
end
[accepted,reason,payload] = problem.audit(y,diag);
if ~islogical(accepted) || ~isscalar(accepted)
    error('invzp:ArcLength:InvalidInput', ...
        'problem.audit ACCEPTED must be a scalar logical.');
end
if isstring(reason)
    if ~isscalar(reason) || ismissing(reason), reason = "invalid"; end
    reason = char(reason);
end
if ~(ischar(reason) && (isempty(reason) || isrow(reason)))
    error('invzp:ArcLength:InvalidInput', ...
        'problem.audit REASON must be a string scalar or character row.');
end
end

function rec = initialRecord(y,R,diag,payload,tangent)
rec = makeRecord(0,y,R,diag,payload,0,0,0,0,NaN,NaN);
rec.tangent_parameter = tangent(end);
end

function rec = makeRecord(index,y,R,diag,payload,ds,iterations, ...
        constraintAbs,correctionRatio,tangentOverlap,borderedRcond)
rec = struct( ...
    'index',index, ...
    'parameter',y(end), ...
    'residual_inf',norm(R,Inf), ...
    'step',ds, ...
    'corrector_iterations',iterations, ...
    'constraint_abs',constraintAbs, ...
    'correction_ratio',correctionRatio, ...
    'tangent_overlap',tangentOverlap, ...
    'tangent_parameter',NaN, ...
    'bordered_rcond',borderedRcond, ...
    'diag',diag, ...
    'audit',payload);
end

function rec = blankAttempt()
rec = struct( ...
    'step_index',0, ...
    'attempt_index',0, ...
    'step',NaN, ...
    'predicted_parameter',NaN, ...
    'corrected_parameter',NaN, ...
    'reason','', ...
    'corrector_iterations',0, ...
    'equation_evaluations',0, ...
    'constraint_abs',NaN, ...
    'correction_ratio',NaN);
end

function rec = attemptRecord(stepIndex,attemptIndex,step,yPred,yCorrected, ...
        reason,corr,correctionRatio)
rec = blankAttempt();
rec.step_index = stepIndex;
rec.attempt_index = attemptIndex;
rec.step = step;
rec.predicted_parameter = yPred(end);
if ~isempty(yCorrected), rec.corrected_parameter = yCorrected(end); end
rec.reason = reason;
rec.corrector_iterations = corr.iterations;
rec.equation_evaluations = corr.equation_evaluations;
rec.constraint_abs = corr.constraint_abs;
rec.correction_ratio = correctionRatio;
end

function cfg = resolveOptions(opts,nvar,invalidId)
if ~isstruct(opts) || ~isscalar(opts)
    fail(invalidId,'opts must be a scalar struct.');
end
cfg = struct( ...
    'step_initial',getf(opts,'step_initial',1e-2), ...
    'step_min',getf(opts,'step_min',1e-8), ...
    'step_max',getf(opts,'step_max',5e-2), ...
    'step_grow',getf(opts,'step_grow',1.25), ...
    'step_shrink',getf(opts,'step_shrink',0.5), ...
    'max_steps',getf(opts,'max_steps',100), ...
    'max_corrector',getf(opts,'max_corrector',10), ...
    'max_linesearch',getf(opts,'max_linesearch',10), ...
    'max_step_attempts',getf(opts,'max_step_attempts',8), ...
    'max_evaluations_per_attempt', ...
        getf(opts,'max_evaluations_per_attempt',16), ...
    'target_iterations',getf(opts,'target_iterations',4), ...
    'tol_residual',getf(opts,'tol_residual',1e-10), ...
    'tol_constraint',getf(opts,'tol_constraint',1e-10), ...
    'rcond_min',getf(opts,'rcond_min',1e-12), ...
    'correction_ratio_max',getf(opts,'correction_ratio_max',0.5), ...
    'grow_correction_ratio',getf(opts,'grow_correction_ratio',0.1), ...
    'tangent_overlap_min',getf(opts,'tangent_overlap_min',0.5), ...
    'tangent_residual_tol',getf(opts,'tangent_residual_tol',1e-8), ...
    'initial_direction',getf(opts,'initial_direction',1), ...
    'initial_tangent',getf(opts,'initial_tangent',[]));

positive = {'step_initial','step_min','step_max','step_grow', ...
    'step_shrink','tol_residual','tol_constraint','rcond_min', ...
    'correction_ratio_max','grow_correction_ratio','tangent_residual_tol'};
for k = 1:numel(positive)
    name = positive{k};
    value = cfg.(name);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value <= 0
        fail(invalidId,'opts.%s must be a finite positive scalar.',name);
    end
end
if cfg.step_min > cfg.step_initial || cfg.step_initial > cfg.step_max
    fail(invalidId,'Require step_min <= step_initial <= step_max.');
end
if cfg.step_grow <= 1 || cfg.step_shrink >= 1
    fail(invalidId,'Require step_grow > 1 and step_shrink < 1.');
end
if cfg.grow_correction_ratio > cfg.correction_ratio_max
    fail(invalidId, ...
        'grow_correction_ratio cannot exceed correction_ratio_max.');
end
integers = {'max_steps','max_corrector','max_linesearch', ...
    'max_step_attempts','max_evaluations_per_attempt','target_iterations'};
for k = 1:numel(integers)
    name = integers{k};
    value = cfg.(name);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
            ~isfinite(value) || value < 1 || value ~= floor(value)
        fail(invalidId,'opts.%s must be a positive integer.',name);
    end
end
if ~isnumeric(cfg.initial_direction) || ~isreal(cfg.initial_direction) || ...
        ~isscalar(cfg.initial_direction) || abs(cfg.initial_direction) ~= 1
    fail(invalidId,'opts.initial_direction must be +1 or -1.');
end
if ~isnumeric(cfg.tangent_overlap_min) || ~isreal(cfg.tangent_overlap_min) || ...
        ~isscalar(cfg.tangent_overlap_min) || ~isfinite(cfg.tangent_overlap_min) || ...
        cfg.tangent_overlap_min < -1 || cfg.tangent_overlap_min > 1
    fail(invalidId,'opts.tangent_overlap_min must lie in [-1,1].');
end
if ~isempty(cfg.initial_tangent)
    if ~isnumeric(cfg.initial_tangent) || ~isreal(cfg.initial_tangent) || ...
            ~isvector(cfg.initial_tangent) || numel(cfg.initial_tangent) ~= nvar || ...
            any(~isfinite(cfg.initial_tangent),'all') || norm(cfg.initial_tangent) == 0
        fail(invalidId, ...
            'opts.initial_tangent must be an optional finite nonzero nvar-vector.');
    end
    cfg.initial_tangent = cfg.initial_tangent(:);
end
end

function validateProblem(problem,invalidId)
if ~isstruct(problem) || ~isscalar(problem) || ...
        ~isfield(problem,'equations') || ~isa(problem.equations,'function_handle')
    fail(invalidId,'problem.equations must be a function handle in a scalar struct.');
end
for name = {'event','audit'}
    if isfield(problem,name{1}) && ~isempty(problem.(name{1})) && ...
            ~isa(problem.(name{1}),'function_handle')
        fail(invalidId,'problem.%s must be a function handle.',name{1});
    end
end
end

function value = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end

function fail(identifier,message,varargin)
error(identifier,message,varargin{:});
end
