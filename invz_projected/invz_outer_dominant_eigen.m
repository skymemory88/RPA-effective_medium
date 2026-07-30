function out = invz_outer_dominant_eigen(mapfun, Sigma, opts)
%INVZ_OUTER_DOMINANT_EIGEN Matrix-free dominant eigenvalue of dF/dSigma.
% Uses central finite-difference Jacobian-vector products and power iteration.
% This is a local diagnostic, not an outer solver. If either perturbed state
% leaves the admissible map domain, the result reports domain_boundary.
if nargin < 3, opts = struct(); end
if ~isa(mapfun,'function_handle')
    error('invz:outerEigen','mapfun must be a function handle.');
end
Sigma = Sigma(:);
n = numel(Sigma);
maxit = getf(opts,'maxit',30);
tol = getf(opts,'tol',1e-5);
h = getf(opts,'fd_step',1e-6*max(1,norm(Sigma,inf)));
if ~(isscalar(h) && isfinite(h) && h > 0)
    error('invz:outerEigen','fd_step must be a finite positive scalar.');
end
if isfield(opts,'seed') && ~isempty(opts.seed)
    v = opts.seed(:);
    if numel(v) ~= n || any(~isfinite(v)) || norm(v) == 0
        error('invz:outerEigen','seed must be a finite nonzero vector matching Sigma.');
    end
else
    v = cos((1:n)'*sqrt(2));
end
v = v/norm(v);

base = mapfun(Sigma);
out = struct('status','base_map_undefined','converged',false,'lambda',NaN, ...
    'eigen_residual',NaN,'vector',nan(n,1),'fd_step',h,'iterations',0, ...
    'lambda_history',nan(maxit,1),'residual_history',nan(maxit,1), ...
    'map_calls',1,'base_status',string(base.status), ...
    'boundary_status_plus',"",'boundary_status_minus',"");
if ~(isfield(base,'defined') && base.defined)
    return;
end

for it = 1:maxit
    [jv,ok,sp,sm] = jacobian_vector(v);
    if ~ok
        out.status = 'domain_boundary';
        out.boundary_status_plus = sp;
        out.boundary_status_minus = sm;
        out.iterations = it;
        trim_history();
        return;
    end
    scale = norm(jv);
    if ~(isfinite(scale) && scale > 0)
        out.status = 'zero_or_nonfinite_jacobian_action';
        out.iterations = it;
        trim_history();
        return;
    end
    lambda = real(v'*jv);
    eigres = norm(jv-lambda*v)/scale;
    out.lambda_history(it) = lambda;
    out.residual_history(it) = eigres;
    out.iterations = it;
    v = jv/scale;
    if eigres < tol && (it == 1 || ...
            abs(lambda-out.lambda_history(it-1)) <= tol*max(1,abs(lambda)))
        out.converged = true;
        break;
    end
end
out.lambda = out.lambda_history(out.iterations);
out.eigen_residual = out.residual_history(out.iterations);
out.vector = v;
if out.converged, out.status = 'ok';
else,             out.status = 'power_iteration_not_converged'; end
trim_history();

    function [jv,ok,sp,sm] = jacobian_vector(direction)
        plus = mapfun(Sigma+h*direction);
        minus = mapfun(Sigma-h*direction);
        out.map_calls = out.map_calls+2;
        sp = string(plus.status);
        sm = string(minus.status);
        ok = isfield(plus,'defined') && plus.defined && ...
             isfield(minus,'defined') && minus.defined;
        if ok
            jv = (plus.Sigma_map(:)-minus.Sigma_map(:))/(2*h);
        else
            jv = nan(n,1);
        end
    end

    function trim_history()
        out.lambda_history = out.lambda_history(1:out.iterations);
        out.residual_history = out.residual_history(1:out.iterations);
    end
end
