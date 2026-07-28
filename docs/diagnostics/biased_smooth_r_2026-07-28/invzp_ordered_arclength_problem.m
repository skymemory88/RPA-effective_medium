function problem = invzp_ordered_arclength_problem(ctx, opts)
%INVZP_ORDERED_ARCLENGTH_PROBLEM Adapt resummed invZ nodes to scaled arclength.
%
%   Y=[Sigma(:)./sigma_scale; K0/K0_scale; h/h_scale].
%   The h-column of the residual Jacobian is obtained by centred,
%   Richardson-refined node reconstruction at fixed physical [Sigma;K0].

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1')
    error('invzp:OrderedArc:InvalidInput', ...
        'ctx must come from invz_ordered_node_context.');
end
if ~isstruct(opts) || ~isscalar(opts)
    error('invzp:OrderedArc:InvalidInput','opts must be a scalar struct.');
end

nw = numel(ctx.wn);
sigmaScale = getf(opts,'sigma_scale',1);
if isscalar(sigmaScale), sigmaScale = repmat(sigmaScale,nw,1); end
if ~isnumeric(sigmaScale) || ~isreal(sigmaScale) || ...
        ~isvector(sigmaScale) || numel(sigmaScale) ~= nw || ...
        any(~isfinite(sigmaScale),'all') || any(sigmaScale <= 0)
    error('invzp:OrderedArc:InvalidInput', ...
        'opts.sigma_scale must be a positive scalar or nw-vector.');
end
sigmaScale = sigmaScale(:);
K0Scale = positive(getf(opts,'K0_scale',ctx.Jscale),'K0_scale');
hScale = positive(getf(opts,'h_scale',ctx.Jscale),'h_scale');
hFdRelative = positive(getf(opts,'h_fd_relative',1e-4),'h_fd_relative');
hFdAbsolute = getf(opts,'h_fd_absolute',0);
if ~isnumeric(hFdAbsolute) || ~isreal(hFdAbsolute) || ...
        ~isscalar(hFdAbsolute) || ~isfinite(hFdAbsolute) || hFdAbsolute < 0
    error('invzp:OrderedArc:InvalidInput', ...
        'opts.h_fd_absolute must be a finite nonnegative scalar.');
end
tolOuter = positive(getf(opts,'tol_outer',1e-8),'tol_outer');
poleMarginMin = positive(getf(opts,'pole_margin_min',1e-10),'pole_margin_min');
meanMarginMin = positive(getf(opts,'mean_margin_min',1e-10),'mean_margin_min');
hDomain = getf(opts,'h_domain',[-Inf Inf]);
if ~isnumeric(hDomain) || ~isreal(hDomain) || ~isequal(size(hDomain),[1 2]) || ...
        isnan(hDomain(1)) || isnan(hDomain(2)) || hDomain(1) >= hDomain(2)
    error('invzp:OrderedArc:InvalidInput', ...
        'opts.h_domain must be an increasing real [lower upper] row.');
end

stateScale = [sigmaScale;K0Scale];
cfg = struct('schema','invzp_ordered_arclength_problem/v1', ...
    'sigma_scale',sigmaScale,'K0_scale',K0Scale,'h_scale',hScale, ...
    'h_fd_relative',hFdRelative,'h_fd_absolute',hFdAbsolute, ...
    'tol_outer',tolOuter,'pole_margin_min',poleMarginMin, ...
    'mean_margin_min',meanMarginMin,'h_domain',hDomain);
problem = struct( ...
    'equations',@equations, ...
    'event',@event, ...
    'audit',@audit, ...
    'pack',@pack, ...
    'unpack',@unpack, ...
    'config',cfg);

    function [R,Jscaled,record] = equations(y)
        [u,h] = unpack(y);
        [node,local] = invz_ordered_make_node(ctx,h);
        if isempty(node)
            R = nan(nw+1,1);
            Jscaled = nan(nw+1,nw+2);
            record = invalidRecord(h,u,local,'center_local_domain');
            return
        end
        [R,physics,state,Ju] = invz_ordered_node_equations(node,u);

        dh = max(hFdAbsolute,hFdRelative*max(hScale,abs(h)));
        [RhCoarse,statusCoarse] = hDerivative(u,h,dh);
        [RhHalf,statusHalf] = hDerivative(u,h,dh/2);
        Rh = (4*RhHalf-RhCoarse)/3;
        derivativeDrift = norm(Rh-RhHalf,Inf)/max(1,norm(Rh,Inf));
        Jscaled = [Ju.*stateScale.',Rh*hScale];
        derivativeStatus = statusCoarse;
        if strcmp(derivativeStatus,'ok'), derivativeStatus = statusHalf; end

        margins = struct( ...
            'h_lower',h-hDomain(1), ...
            'h_upper',hDomain(2)-h, ...
            'pole',physics.pole_margin-poleMarginMin, ...
            'mean',physics.mean_margin-meanMarginMin);
        record = struct('h',h,'u',u,'node',node,'state',state, ...
            'local',local,'physics',physics,'event_margins',margins, ...
            'h_fd_step',dh,'h_jacobian_scaled_drift',derivativeDrift, ...
            'h_derivative_status',derivativeStatus);
    end

    function [Rh,status] = hDerivative(u,h,dh)
        [nodePlus,metaPlus] = invz_ordered_make_node(ctx,h+dh);
        [nodeMinus,metaMinus] = invz_ordered_make_node(ctx,h-dh);
        if isempty(nodePlus) || isempty(nodeMinus)
            Rh = nan(nw+1,1);
            status = sprintf('centred_local_domain:%s/%s', ...
                metaMinus.status,metaPlus.status);
            return
        end
        Rplus = invz_ordered_node_equations(nodePlus,u);
        Rminus = invz_ordered_node_equations(nodeMinus,u);
        Rh = (Rplus-Rminus)/(2*dh);
        status = 'ok';
    end

    function reason = event(record)
        reason = '';
        p = record.physics;
        m = record.event_margins;
        if ~strcmp(record.local.status,'ok')
            reason = record.local.status;
        elseif ~strcmp(record.h_derivative_status,'ok')
            reason = record.h_derivative_status;
        elseif m.h_lower < 0
            reason = 'h_below_domain';
        elseif m.h_upper < 0
            reason = 'h_above_domain';
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
        auditRecord = invz_ordered_residual( ...
            record.node,record.state,struct('tol_outer',tolOuter));
        accepted = logical(auditRecord.accepted);
        if accepted
            reason = '';
        else
            reason = 'A-D audit rejected the corrected root';
        end
        payload = struct('h',record.h,'state',record.state, ...
            'local',record.local,'physics',record.physics, ...
            'event_margins',record.event_margins, ...
            'h_fd_step',record.h_fd_step, ...
            'h_jacobian_scaled_drift',record.h_jacobian_scaled_drift, ...
            'audit',auditRecord);
    end

    function record = invalidRecord(h,u,local,derivativeStatus)
        physics = struct('r',NaN,'Gtil0',NaN,'pole_margin',NaN, ...
            'mean_margin',NaN);
        margins = struct('h_lower',h-hDomain(1), ...
            'h_upper',hDomain(2)-h,'pole',NaN,'mean',NaN);
        state = struct('Sigma',[],'K',[],'lam',[],'K0s',NaN);
        record = struct('h',h,'u',u,'node',[],'state',state, ...
            'local',local,'physics',physics,'event_margins',margins, ...
            'h_fd_step',NaN,'h_jacobian_scaled_drift',NaN, ...
            'h_derivative_status',derivativeStatus);
    end

    function y = pack(state,h)
        if ~isstruct(state) || ~isfield(state,'Sigma') || ~isfield(state,'K0s')
            error('invzp:OrderedArc:InvalidInput', ...
                'pack requires a state with Sigma and K0s.');
        end
        if numel(state.Sigma) ~= nw
            error('invzp:OrderedArc:InvalidInput', ...
                'state.Sigma does not match the context Matsubara grid.');
        end
        y = [state.Sigma(:)./sigmaScale;state.K0s/K0Scale;h/hScale];
    end

    function [u,h] = unpack(y)
        if ~isnumeric(y) || ~isreal(y) || ~isvector(y) || ...
                numel(y) ~= nw+2 || any(~isfinite(y),'all')
            error('invzp:OrderedArc:InvalidInput', ...
                'y must be a finite scaled vector with nw+2 entries.');
        end
        y = y(:);
        u = [y(1:nw).*sigmaScale;y(nw+1)*K0Scale];
        h = y(end)*hScale;
    end
end

function value = positive(value,name)
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || value <= 0
    error('invzp:OrderedArc:InvalidInput', ...
        'opts.%s must be a finite positive scalar.',name);
end
end

function value = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
