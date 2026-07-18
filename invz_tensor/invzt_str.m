function s = invzt_str(x)
%INVZT_STR Compact string form of x for error messages.
%   s = INVZT_STR(x) is char(x) when x is a char row or a scalar string, and
%   mat2str(x) otherwise (numeric/logical/array). Shared error-message helper
%   for the invz_tensor drivers -- replaces the per-file local_str /
%   local_conv_str copies.
if ischar(x) || (isstring(x) && isscalar(x))
    s = char(x);
else
    s = mat2str(x);
end
end
