function problem = invzp_ordered_squared_field_problem(ctx,hReference,opts)
%INVZP_ORDERED_SQUARED_FIELD_PROBLEM Positive-h endpoint coordinate.
%
%   Y=[Sigma(:)./sigma_scale; K0/K0_scale; q], where
%
%       q = (h/hReference)^2,  h >= 0.
%
%   For q>0 this is an exact one-to-one reparameterization of the original
%   fixed-h residual. It neither averages +h/-h nor alters an equation.
%   Every q>0 Jacobian uses centred Richardson differences with half-width
%   min(opts.q_fd_step,q/2), so every stencil node remains strictly
%   positive. A fourth-order forward derivative is retained only to report
%   q=0 diagnostics. An otherwise valid q=0 point is rejected as an
%   unresolved endpoint event: a finite stencil cannot prove that the
%   residual has no sqrt(q) term.
%   opts.q_domain bounds accepted central trace points only; derivative
%   stencils may cross that reporting bound and must remain node-valid.
%   opts.static_polish (default false) exposes one optional one-shot K0
%   Newton proposal for machine-resolution stalls. It is eligible only
%   after all Sigma residuals pass, is capped by
%   opts.static_polish_max_ulps in physical K0, and never changes an
%   acceptance tolerance. The tracer independently recomputes the complete
%   proposal and applies its ordinary constraint, event, A--D, rank, and
%   tangent gates.

if nargin < 3 || isempty(opts), opts = struct(); end
invalidId = 'invzp:OrderedSquaredField:InvalidInput';
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1')
    error(invalidId,'ctx must come from invz_ordered_node_context.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error(invalidId,'opts must be a scalar struct.');
end
hReference = positive(hReference,'hReference',invalidId);

nw = numel(ctx.wn);
sigmaScale = getf(opts,'sigma_scale',1);
if isscalar(sigmaScale), sigmaScale = repmat(sigmaScale,nw,1); end
if ~isnumeric(sigmaScale) || ~isreal(sigmaScale) || ...
        ~isvector(sigmaScale) || numel(sigmaScale) ~= nw || ...
        any(~isfinite(sigmaScale),'all') || any(sigmaScale <= 0)
    error(invalidId, ...
        'opts.sigma_scale must be a positive scalar or nw-vector.');
end
sigmaScale = sigmaScale(:);
K0Scale = positive(getf(opts,'K0_scale',ctx.Jscale), ...
    'K0_scale',invalidId);
qFdStep = positive(getf(opts,'q_fd_step',0.1), ...
    'q_fd_step',invalidId);
qDriftMax = positive(getf(opts,'q_jacobian_drift_max',1e-2), ...
    'q_jacobian_drift_max',invalidId);
staticPolish = logicalScalar(getf(opts,'static_polish',false), ...
    'static_polish',invalidId);
staticPolishMaxUlps = positiveInteger( ...
    getf(opts,'static_polish_max_ulps',4096), ...
    'static_polish_max_ulps',invalidId);
tolOuter = positive(getf(opts,'tol_outer',1e-8), ...
    'tol_outer',invalidId);
poleMarginMin = positive(getf(opts,'pole_margin_min',1e-10), ...
    'pole_margin_min',invalidId);
meanMarginMin = positive(getf(opts,'mean_margin_min',1e-10), ...
    'mean_margin_min',invalidId);
qDomain = getf(opts,'q_domain',[0 Inf]);
if ~isnumeric(qDomain) || ~isreal(qDomain) || ...
        ~isequal(size(qDomain),[1 2]) || isnan(qDomain(1)) || ...
        isnan(qDomain(2)) || ~isfinite(qDomain(1)) || qDomain(1) < 0 || ...
        qDomain(1) >= qDomain(2)
    error(invalidId, ...
        'opts.q_domain must be an increasing nonnegative [lower upper] row.');
end

stateScale = [sigmaScale;K0Scale];
cachedY = [];
cachedR = [];
cachedJ = [];
cachedRecord = [];
cfg = struct( ...
    'schema','invzp_ordered_squared_field_problem/v1', ...
    'h_reference',hReference, ...
    'sigma_scale',sigmaScale,'K0_scale',K0Scale, ...
    'q_fd_step',qFdStep, ...
    'q_jacobian_drift_max',qDriftMax, ...
    'static_polish',staticPolish, ...
    'static_polish_max_ulps',staticPolishMaxUlps, ...
    'tol_outer',tolOuter, ...
    'pole_margin_min',poleMarginMin, ...
    'mean_margin_min',meanMarginMin, ...
    'q_domain',qDomain);
polishCallback = [];
if staticPolish
    polishCallback = @polishStaticClosure;
end
problem = struct( ...
    'equations',@equations, ...
    'event',@event, ...
    'audit',@audit, ...
    'polish',polishCallback, ...
    'pack',@pack, ...
    'unpack',@unpack, ...
    'config',cfg);

    function [R,Jscaled,record] = equations(y)
        y = y(:);
        if ~isempty(cachedY) && isequal(y,cachedY)
            R = cachedR;
            Jscaled = cachedJ;
            record = cachedRecord;
            return
        end
        [u,q] = unpack(y);
        [R,physics,state,Ju,node,local] = equationAtQ(u,q,true);
        if isempty(node)
            Jscaled = nan(nw+1,nw+2);
            record = invalidRecord(q,u,local,'q_center_domain');
            updateCache(y,R,Jscaled,record);
            return
        end
        [Jq,reference,derivativeStatus,scheme,usedStep] = ...
            qDerivative(u,q,R);
        derivativeDrift = norm(Jq-reference,Inf)/max(1,norm(Jq,Inf));
        Jscaled = [Ju.*stateScale.',Jq];
        margins = struct( ...
            'q_lower',q-qDomain(1), ...
            'q_upper',qDomain(2)-q, ...
            'pole',physics.pole_margin-poleMarginMin, ...
            'mean',physics.mean_margin-meanMarginMin, ...
            'derivative_drift',qDriftMax-derivativeDrift);
        record = struct( ...
            'q',q,'h',hReference*sqrt(q),'u',u, ...
            'node',node,'state',state,'local',local, ...
            'physics',physics,'event_margins',margins, ...
            'q_fd_step',usedStep, ...
            'q_derivative_scheme',scheme, ...
            'q_jacobian_drift',derivativeDrift, ...
            'q_derivative_status',derivativeStatus);
        updateCache(y,R,Jscaled,record);
    end

    function [column,reference,status,scheme,step] = ...
            qDerivative(u,q,Rcenter)
        if q > 0
            step = min(qFdStep,q/2);
            scheme = 'adaptive_centered_richardson';
            half = step/2;
            stencil = [q-step,q-half,q+half,q+step];
            if ~(isfinite(step) && step > 0 && isfinite(half) && half > 0) || ...
                    any(~isfinite(stencil),'all') || ...
                    ~(stencil(1) < stencil(2) && stencil(2) < q && ...
                    q < stencil(3) && stencil(3) < stencil(4))
                column = nan(nw+1,1);
                reference = column;
                status = 'q_fd_resolution';
                return
            end
            [Rplus,statusPlus] = residualAtQ(u,q+step);
            [Rminus,statusMinus] = residualAtQ(u,q-step);
            coarse = (Rplus-Rminus)/(2*step);
            [RplusHalf,statusPlusHalf] = residualAtQ(u,q+half);
            [RminusHalf,statusMinusHalf] = residualAtQ(u,q-half);
            fine = (RplusHalf-RminusHalf)/(2*half);
            column = (4*fine-coarse)/3;
            reference = fine;
            status = firstFailure( ...
                {statusPlus,statusMinus,statusPlusHalf,statusMinusHalf});
            return
        end

        step = qFdStep;
        scheme = 'forward_fourth_order';
        values = zeros(nw+1,5);
        values(:,1) = Rcenter;
        statuses = repmat({''},1,4);
        for offset = 1:4
            [values(:,offset+1),statuses{offset}] = ...
                residualAtQ(u,q+offset*step);
        end
        column = (-25*values(:,1)+48*values(:,2)- ...
            36*values(:,3)+16*values(:,4)-3*values(:,5))/ ...
            (12*step);
        reference = (-3*values(:,1)+4*values(:,2)-values(:,3))/ ...
            (2*step);
        status = firstFailure(statuses);
    end

    function [R,status] = residualAtQ(u,q)
        [R,~,~,~,node,local] = equationAtQ(u,q,false);
        if isempty(node)
            status = local.status;
        else
            status = 'ok';
        end
    end

    function [R,physics,state,Ju,node,local] = equationAtQ(u,q,needJacobian)
        physics = invalidPhysics();
        state = blankState();
        Ju = nan(nw+1);
        if ~(isnumeric(q) && isreal(q) && isscalar(q) && ...
                isfinite(q) && q >= 0)
            R = nan(nw+1,1);
            node = [];
            local = struct('status','q_out_of_domain','q',q,'h',NaN);
            return
        end
        h = hReference*sqrt(q);
        [node,local] = invz_ordered_make_node(ctx,h);
        local.q = q;
        if isempty(node)
            R = nan(nw+1,1);
            return
        end
        if needJacobian
            [R,physics,state,Ju] = invz_ordered_node_equations(node,u);
        else
            [R,physics,state] = invz_ordered_node_equations(node,u);
        end
    end

    function reason = event(record)
        reason = '';
        m = record.event_margins;
        p = record.physics;
        if ~strcmp(record.local.status,'ok')
            reason = record.local.status;
        elseif ~strcmp(record.q_derivative_status,'ok')
            reason = record.q_derivative_status;
        elseif m.q_lower < 0
            reason = 'q_below_domain';
        elseif m.q_upper < 0
            reason = 'q_above_domain';
        elseif record.q == 0
            reason = 'q_zero_endpoint_unresolved';
        elseif ~isfinite(m.derivative_drift) || m.derivative_drift < 0
            reason = 'q_jacobian_drift';
        elseif ~isfinite(p.r)
            reason = 'unbounded_integrand';
        elseif ~isfinite(p.Gtil0)
            reason = 'nonfinite_Gtil0';
        elseif ~isfinite(m.pole) || m.pole <= 0
            reason = 'pole_event';
        elseif ~isfinite(m.mean) || m.mean <= 0
            reason = 'mean_event';
        end
    end

    function [accepted,reason,payload] = audit(~,record)
        eventReason = event(record);
        if isempty(eventReason)
            auditRecord = invz_ordered_residual( ...
                record.node,record.state,struct('tol_outer',tolOuter));
            accepted = auditRecord.accepted;
            if accepted
                reason = '';
            else
                reason = 'A-D audit rejected the corrected root';
            end
        else
            auditRecord = struct([]);
            accepted = false;
            reason = ['event:',eventReason];
        end
        payload = struct( ...
            'q',record.q,'h',record.h,'state',record.state, ...
            'local',record.local,'physics',record.physics, ...
            'event_margins',record.event_margins, ...
            'q_fd_step',record.q_fd_step, ...
            'q_derivative_scheme',record.q_derivative_scheme, ...
            'q_jacobian_drift',record.q_jacobian_drift, ...
            'audit',auditRecord);
    end

    function y = pack(state,h)
        if ~isstruct(state) || ~isscalar(state) || ...
                ~isfield(state,'Sigma') || ~isfield(state,'K0s') || ...
                ~isnumeric(state.Sigma) || ~isreal(state.Sigma) || ...
                ~isvector(state.Sigma) || numel(state.Sigma) ~= nw || ...
                any(~isfinite(state.Sigma),'all') || ...
                ~isnumeric(state.K0s) || ~isreal(state.K0s) || ...
                ~isscalar(state.K0s) || ~isfinite(state.K0s)
            error(invalidId, ...
                'pack requires a scalar state with nw-vector Sigma and K0s.');
        end
        if ~isnumeric(h) || ~isreal(h) || ~isscalar(h) || ...
                ~isfinite(h) || h < 0
            error(invalidId,'pack requires a finite nonnegative h.');
        end
        q = (h/hReference)^2;
        if ~isfinite(q) || (h > 0 && q == 0)
            error(invalidId, ...
                'h/hReference must square to a finite positive q when h is positive.');
        end
        y = [state.Sigma(:)./sigmaScale; ...
            state.K0s/K0Scale;q];
    end

    function [candidate,info] = polishStaticClosure( ...
            y,R,Jscaled,record,tolerance,budget)
        candidate = y(:);
        info = struct('applied',false,'reason','not_ready', ...
            'equation_evaluations',0);
        k0Index = nw+1;
        if isempty(record.node) || numel(R) ~= nw+1 || ...
                ~isequal(size(Jscaled),[nw+1,nw+2]) || ...
                norm(R(1:nw),Inf) > tolerance || ...
                abs(R(end)) <= tolerance
            return
        elseif budget < 1
            info.reason = 'no_budget';
            return
        end
        derivative = Jscaled(end,k0Index);
        if ~(isfinite(derivative) && derivative ~= 0)
            info.reason = 'invalid_static_derivative';
            return
        end
        x0 = candidate(k0Index);
        coordinateUlp = eps(abs(x0));
        physicalUlp = eps(abs(record.u(end)));
        step = -R(end)/derivative;
        if ~(isfinite(step) && step ~= 0)
            info.reason = 'outside_ulp_envelope';
            return
        end

        x1 = x0+step;
        if x1 == x0
            x1 = x0+sign(step)*coordinateUlp;
        end
        K1 = x1*K0Scale;
        if ~isfinite(K1) || ...
                abs(K1-record.u(end))/physicalUlp > staticPolishMaxUlps
            info.reason = 'outside_ulp_envelope';
            return
        end
        [R1,valid] = atPhysicalK0(K1);
        info.equation_evaluations = info.equation_evaluations+1;
        if valid && norm(R1,Inf) <= tolerance
            accept(x1);
        elseif valid
            info.reason = 'static_K0_candidate_rejected';
        else
            info.reason = 'invalid_static_trial';
        end

        function [Rt,valid] = atPhysicalK0(K0)
            ut = record.u;
            ut(end) = K0;
            Rt = invz_ordered_node_equations(record.node,ut);
            valid = isnumeric(Rt) && isreal(Rt) && ...
                isequal(size(Rt),[nw+1,1]) && all(isfinite(Rt),'all');
        end

        function accept(x)
            candidate(k0Index) = x;
            info.applied = true;
            info.reason = 'accepted_static_K0_polish';
        end
    end

    function [u,q] = unpack(y)
        if ~isnumeric(y) || ~isreal(y) || ~isvector(y) || ...
                numel(y) ~= nw+2 || any(~isfinite(y),'all')
            error(invalidId, ...
                'y must be a finite scaled vector with nw+2 entries.');
        end
        y = y(:);
        u = [y(1:nw).*sigmaScale;y(nw+1)*K0Scale];
        q = y(end);
    end

    function record = invalidRecord(q,u,local,derivativeStatus)
        if isfinite(q) && q >= 0
            h = hReference*sqrt(q);
        else
            h = NaN;
        end
        margins = struct('q_lower',NaN,'q_upper',NaN, ...
            'pole',NaN,'mean',NaN,'derivative_drift',NaN);
        record = struct( ...
            'q',q,'h',h,'u',u,'node',[],'state',blankState(), ...
            'local',local,'physics',invalidPhysics(), ...
            'event_margins',margins,'q_fd_step',qFdStep, ...
            'q_derivative_scheme','not_evaluated', ...
            'q_jacobian_drift',NaN, ...
            'q_derivative_status',derivativeStatus);
    end

    function updateCache(y,R,Jscaled,record)
        cachedY = y;
        cachedR = R;
        cachedJ = Jscaled;
        cachedRecord = record;
    end
end

function status = firstFailure(statuses)
status = 'ok';
for k = 1:numel(statuses)
    if ~strcmp(statuses{k},'ok')
        status = statuses{k};
        return
    end
end
end

function physics = invalidPhysics()
physics = struct('r',NaN,'Gtil0',NaN, ...
    'pole_margin',NaN,'mean_margin',NaN);
end

function state = blankState()
state = struct('Sigma',[],'K',[],'lam',[],'K0s',NaN);
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

function value = logicalScalar(value,name,invalidId)
if ~islogical(value) || ~isscalar(value)
    error(invalidId,'%s must be a scalar logical.',name);
end
end

function value = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
