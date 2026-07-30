function out = invz_outer_picard_diagnostic(mapfun, Sigma0, opts)
%INVZ_OUTER_PICARD_DIAGNOSTIC Domain-gated Picard existence probe.
% Diagnostic only. Every iterate must have a uniquely defined admissible map;
% otherwise the probe stops without exporting a state.
if nargin < 3, opts = struct(); end
if ~isa(mapfun,'function_handle')
    error('invz:outerPicard','mapfun must be a function handle.');
end
mix = getf(opts,'mix',1);
tol = getf(opts,'tol',1e-8);
maxit = getf(opts,'maxit',200);
if ~(isscalar(mix) && isfinite(mix) && mix > 0 && mix <= 1)
    error('invz:outerPicard','mix must lie in (0,1].');
end
Sigma = Sigma0(:);
residual_history = nan(maxit,1);
map_status = strings(maxit,1);
last = [];
last_admissible_map = [];
last_admissible_Sigma = nan(size(Sigma));
status = 'max_iterations';
for it = 1:maxit
    last = mapfun(Sigma);
    map_status(it) = string(last.status);
    if ~(isfield(last,'defined') && last.defined)
        status = 'left_admissible_domain';
        break;
    end
    last_admissible_map = last;
    last_admissible_Sigma = Sigma;
    residual_history(it) = last.residual_norm;
    if last.residual_norm <= tol
        status = 'ok';
        break;
    end
    Sigma = Sigma + mix*(last.Sigma_map(:)-Sigma);
end
converged = strcmp(status,'ok');
out = struct('status',status,'converged',converged,'Sigma',nan(size(Sigma)), ...
    'iterations',it,'residual_history',residual_history(1:it), ...
    'map_status',map_status(1:it),'last_map',last, ...
    'last_admissible_map',last_admissible_map, ...
    'last_admissible_Sigma',last_admissible_Sigma,'mix',mix,'tol',tol);
if converged
    out.Sigma = Sigma;
end
end
