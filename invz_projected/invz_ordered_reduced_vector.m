function residual = invz_ordered_reduced_vector(y, ctx, opts)
%INVZ_ORDERED_REDUCED_VECTOR Numeric residual adapter for nonlinear solvers.
% Returns NaNs when the diagnostic reduced residual is outside its defined
% physical domain; callers must use a globalization strategy that treats such
% evaluations as rejected trials rather than roots.
if nargin < 3, opts = struct(); end
out = invz_ordered_reduced_residual(y,ctx,opts);
residual = out.residual;
end
