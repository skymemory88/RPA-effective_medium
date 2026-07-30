function out = invz_outer_arnoldi_diagnostic(mapfun, Sigma, opts)
%INVZ_OUTER_ARNOLDI_DIAGNOSTIC Matrix-free leading modes of dF/dSigma.
% Uses central finite-difference Jacobian-vector products and MATLAB's
% nonsymmetric Arnoldi eigensolver. Unlike power iteration, this diagnostic
% can resolve a leading complex-conjugate pair. It is not an outer solver.
%
% OUT.status is domain_boundary if any finite-difference perturbation leaves
% the admissible map domain. OUT.converged additionally requires ARPACK's
% convergence flag to be zero and every returned eigenpair to pass the
% requested residual tolerance.
if nargin < 3, opts = struct(); end
if ~isa(mapfun,'function_handle')
    error('invz:outerArnoldi','mapfun must be a function handle.');
end
Sigma = Sigma(:);
n = numel(Sigma);
if n < 3 || any(~isfinite(Sigma)) || ~isreal(Sigma)
    error('invz:outerArnoldi','Sigma must be a finite real vector of length >=3.');
end

k = getf(opts,'n_eigs',min(6,n-2));
tol = getf(opts,'tol',1e-5);
maxit = getf(opts,'maxit',100);
h = getf(opts,'fd_step',1e-6*max(1,norm(Sigma,inf)));
subspace = getf(opts,'subspace_dimension',min(n,max(20,2*k+1)));
mode_rel_floor = getf(opts,'mode_relative_floor',1e-6);
mode_abs_floor = getf(opts,'mode_absolute_floor',1e-12);
if ~(isscalar(k) && k == round(k) && k >= 1 && k <= n-2)
    error('invz:outerArnoldi','n_eigs must be an integer in [1,numel(Sigma)-2].');
end
if ~(isscalar(h) && isfinite(h) && h > 0)
    error('invz:outerArnoldi','fd_step must be a finite positive scalar.');
end
if ~(isscalar(subspace) && subspace == round(subspace) && ...
        subspace > k+1 && subspace <= n)
    error('invz:outerArnoldi', ...
        'subspace_dimension must be an integer in (n_eigs+1,numel(Sigma)].');
end
if ~(isscalar(mode_rel_floor) && isfinite(mode_rel_floor) && ...
        mode_rel_floor >= 0 && mode_rel_floor < 1 && ...
        isscalar(mode_abs_floor) && isfinite(mode_abs_floor) && ...
        mode_abs_floor >= 0)
    error('invz:outerArnoldi', ...
        'mode residual floors must be finite, nonnegative, and relative_floor < 1.');
end

base = mapfun(Sigma);
out = struct('status','base_map_undefined','converged',false, ...
    'eigenvalues',complex(nan(k,1)),'eigen_residuals',nan(k,1), ...
    'vectors',complex(nan(n,k)),'spectral_radius',NaN,'eigs_flag',NaN, ...
    'active_mode_mask',false(k,1),'max_active_residual',NaN, ...
    'mode_relative_floor',mode_rel_floor,'mode_absolute_floor',mode_abs_floor, ...
    'fd_step',h,'n_eigs',k,'subspace_dimension',subspace,'map_calls',1, ...
    'base_status',string(base.status),'boundary_status_plus',"", ...
    'boundary_status_minus',"",'error_message',"");
if ~(isfield(base,'defined') && base.defined)
    return;
end

boundary_plus = "";
boundary_minus = "";
eopts = struct('issym',false,'isreal',true,'tol',tol, ...
    'maxit',maxit,'p',subspace,'disp',0, ...
    'v0',cos((1:n)'*sqrt(2)));
try
    [V,D,flag] = eigs(@apply_columns,n,k,'largestabs',eopts);
    lambda = diag(D);
    [~,order] = sort(abs(lambda),'descend');
    lambda = lambda(order);
    V = V(:,order);
    eigres = nan(k,1);
    for j = 1:k
        Av = apply_columns(V(:,j));
        eigres(j) = norm(Av-lambda(j)*V(:,j)) / max(norm(Av),eps);
    end
catch err
    if strcmp(err.identifier,'invz:outerArnoldiBoundary')
        out.status = 'domain_boundary';
        out.boundary_status_plus = boundary_plus;
        out.boundary_status_minus = boundary_minus;
    else
        out.status = 'arnoldi_error';
        out.error_message = string(err.message);
    end
    return;
end

out.eigenvalues = lambda;
out.eigen_residuals = eigres;
out.vectors = V;
out.spectral_radius = max(abs(lambda));
out.eigs_flag = flag;
mode_floor = max(mode_abs_floor,mode_rel_floor*out.spectral_radius);
active = abs(lambda) >= mode_floor;
out.active_mode_mask = active;
if any(active)
    out.max_active_residual = max(eigres(active));
end
out.converged = flag == 0 && any(active) && ...
    all(isfinite(eigres(active))) && all(eigres(active) <= tol);
if out.converged
    out.status = 'ok';
elseif flag ~= 0
    out.status = 'arnoldi_not_converged';
elseif ~any(active)
    out.status = 'no_resolved_nonzero_modes';
else
    out.status = 'eigenpair_residual_failed';
end

    function Y = apply_columns(X)
        Y = complex(zeros(size(X)));
        for col = 1:size(X,2)
            xr = real(X(:,col));
            xi = imag(X(:,col));
            Y(:,col) = apply_real(xr);
            if any(xi)
                Y(:,col) = Y(:,col) + 1i*apply_real(xi);
            end
        end
    end

    function y = apply_real(v)
        nv = norm(v);
        if nv == 0
            y = zeros(n,1);
            return;
        end
        delta = h/nv;
        plus = mapfun(Sigma+delta*v);
        minus = mapfun(Sigma-delta*v);
        out.map_calls = out.map_calls+2;
        okplus = isfield(plus,'defined') && plus.defined;
        okminus = isfield(minus,'defined') && minus.defined;
        if ~(okplus && okminus)
            boundary_plus = string(plus.status);
            boundary_minus = string(minus.status);
            error('invz:outerArnoldiBoundary', ...
                'finite-difference perturbation left the admissible map domain.');
        end
        y = (plus.Sigma_map(:)-minus.Sigma_map(:))/(2*delta);
    end
end
