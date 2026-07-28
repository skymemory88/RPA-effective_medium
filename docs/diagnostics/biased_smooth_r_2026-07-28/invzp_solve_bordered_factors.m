function [solution, info] = invzp_solve_bordered_factors( ...
        factors, stateScale, hColumn, tangent, rhs, opts)
%INVZP_SOLVE_BORDERED_FACTORS Exact O(n) solve for the ordered arc border.
%
%   Solves
%
%     [J_u*diag(stateScale)   hColumn] [du] = rhs
%     [tangent(1:end-1)'      tangent(end)] [dh]
%
%   using the exact diagonal-plus-rank-two frequency block in FACTORS and
%   a 2-by-2 Schur complement for [K0,h]. No rank truncation is performed.

if nargin < 6 || isempty(opts), opts = struct(); end
invalidId = 'invzp:BorderedFactors:InvalidInput';
validateFactors(factors,invalidId);
nw = numel(factors.sigma_diagonal);
nstate = nw+1;
nvar = nstate+1;
validateVector(stateScale,nstate,'stateScale',invalidId);
if any(stateScale <= 0)
    error(invalidId,'stateScale entries must be positive.');
end
validateVector(hColumn,nstate,'hColumn',invalidId);
validateVector(tangent,nvar,'tangent',invalidId);
validateVector(rhs,nvar,'rhs',invalidId);
stateScale = stateScale(:);
hColumn = hColumn(:);
tangent = tangent(:);
rhs = rhs(:);

rcondMin = getf(opts,'rcond_min',1e-12);
residualTol = getf(opts,'residual_tol',1e-10);
if ~isnumeric(rcondMin) || ~isreal(rcondMin) || ~isscalar(rcondMin) || ...
        ~isfinite(rcondMin) || rcondMin <= 0
    error(invalidId,'opts.rcond_min must be a finite positive scalar.');
end
if ~isnumeric(residualTol) || ~isreal(residualTol) || ~isscalar(residualTol) || ...
        ~isfinite(residualTol) || residualTol <= 0
    error(invalidId,'opts.residual_tol must be a finite positive scalar.');
end

sigmaScale = stateScale(1:nw);
K0Scale = stateScale(end);
d = factors.sigma_diagonal.*sigmaScale;
U = factors.sigma_left;
V = factors.sigma_right.*sigmaScale;
diagonalMargin = min(abs(d))/max(abs(d));
info = blankInfo(diagonalMargin);
solution = nan(nvar,1);
if any(~isfinite(d)) || any(d == 0) || ...
        ~isfinite(diagonalMargin) || diagonalMargin <= rcondMin
    info.reason = 'frequency_diagonal_rank_loss';
    return
end

DinvU = U./d;
woodbury = eye(size(U,2))+V.'*DinvU;
info.woodbury_rcond = rcond(woodbury);
if ~isfinite(info.woodbury_rcond) || info.woodbury_rcond <= rcondMin
    info.reason = 'woodbury_rank_loss';
    return
end

topBorder = [factors.K0_column*K0Scale,hColumn(1:nw)];
bottomLeft = [factors.static_sigma_row.*sigmaScale.'; ...
    tangent(1:nw).'];
bottomRight = [ ...
    factors.static_K0*K0Scale,hColumn(end); ...
    tangent(nw+1),tangent(end)];

applyInverse = @(Q) applyAInverse(Q,d,DinvU,V,woodbury);
AinvBorder = applyInverse(topBorder);
AinvRhs = applyInverse(rhs(1:nw));
schur = bottomRight-bottomLeft*AinvBorder;
info.schur_rcond = rcond(schur);
info.rank_margin = min([diagonalMargin,info.woodbury_rcond,info.schur_rcond]);
if ~isfinite(info.schur_rcond) || info.schur_rcond <= rcondMin
    info.reason = 'border_schur_rank_loss';
    return
end

borderSolution = schur\(rhs(nw+1:end)-bottomLeft*AinvRhs);
frequencySolution = AinvRhs-AinvBorder*borderSolution;
solution = [frequencySolution;borderSolution];

topCheck = d.*frequencySolution+ ...
    U*(V.'*frequencySolution)+topBorder*borderSolution;
bottomCheck = bottomLeft*frequencySolution+bottomRight*borderSolution;
residual = [topCheck;bottomCheck]-rhs;
info.residual_inf = norm(residual,Inf);
info.relative_residual = info.residual_inf/max(1,norm(rhs,Inf));
if ~all(isfinite(solution)) || ~isfinite(info.relative_residual)
    info.reason = 'nonfinite_solution';
elseif info.relative_residual > residualTol
    info.reason = 'factor_solve_residual';
else
    info.accepted = true;
    info.reason = 'accepted';
end
end

function X = applyAInverse(Q,d,DinvU,V,woodbury)
DinvQ = Q./d;
X = DinvQ-DinvU*(woodbury\(V.'*DinvQ));
end

function info = blankInfo(diagonalMargin)
info = struct('accepted',false,'reason','not_solved', ...
    'diagonal_margin',diagonalMargin,'woodbury_rcond',NaN, ...
    'schur_rcond',NaN,'rank_margin',NaN, ...
    'residual_inf',NaN,'relative_residual',NaN);
end

function validateFactors(factors,invalidId)
required = {'schema','sigma_diagonal','sigma_left','sigma_right', ...
    'K0_column','static_sigma_row','static_K0'};
if ~isstruct(factors) || ~isscalar(factors)
    error(invalidId,'factors must be a scalar struct.');
end
for k = 1:numel(required)
    if ~isfield(factors,required{k})
        error(invalidId,'factors.%s is required.',required{k});
    end
end
nw = numel(factors.sigma_diagonal);
if ~strcmp(factors.schema,'invzp_ordered_node_jacobian_factors/v1') || ...
        ~isequal(size(factors.sigma_left),[nw 2]) || ...
        ~isequal(size(factors.sigma_right),[nw 2]) || ...
        numel(factors.K0_column) ~= nw || ...
        numel(factors.static_sigma_row) ~= nw || ...
        ~isscalar(factors.static_K0)
    error(invalidId,'factors has an invalid schema or shape.');
end
numeric = [factors.sigma_diagonal(:);factors.sigma_left(:); ...
    factors.sigma_right(:);factors.K0_column(:); ...
    factors.static_sigma_row(:);factors.static_K0];
if ~isreal(numeric) || any(~isfinite(numeric))
    error(invalidId,'factors must contain finite real values.');
end
end

function validateVector(value,n,name,invalidId)
if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        numel(value) ~= n || any(~isfinite(value),'all')
    error(invalidId,'%s must be a finite real %d-vector.',name,n);
end
end

function value = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
