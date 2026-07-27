function r = invz_finite_max_abs(a, b)
%INVZ_FINITE_MAX_ABS Fail-closed max(abs(a-b)) for solver/checker residuals.
% Inputs must have identical shape. Any NaN or Inf in either operand makes the scalar residual
% NaN, so an isfinite(residual) acceptance gate fails closed. This avoids MATLAB max's default
% omission of isolated NaNs. Empty, shape-matched inputs have residual zero.
if ~isequal(size(a), size(b))
    error('invz:residualShape', ...
        'residual operands must have identical size: got %s and %s.', ...
        mat2str(size(a)), mat2str(size(b)));
end
if any(~isfinite(a(:))) || any(~isfinite(b(:)))
    r = NaN;
elseif isempty(a)
    r = 0;
else
    r = max(abs(a(:) - b(:)));
end
end
