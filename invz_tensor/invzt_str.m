function s = invzt_str(x)
%INVZT_STR Compact string form of x for error messages -- never throws.
%   s = INVZT_STR(x) is char(x) when x is a char row or a scalar string,
%   mat2str(x) for numeric/logical/string arrays, and a class/size
%   placeholder (e.g. '<1x1 struct>') for anything mat2str cannot format --
%   so an error MESSAGE can always be built, whatever malformed value
%   triggered it (mat2str itself throws on structs/cells/objects, which
%   would mask the intended error identifier; 2026-07-18 T-cut review F2).
%   Shared error-message helper for the invz_tensor drivers -- replaces the
%   per-file local_str / local_conv_str copies.
if ischar(x) || (isstring(x) && isscalar(x))
    s = char(x);
elseif isnumeric(x) || islogical(x) || isstring(x)
    s = mat2str(x);
else
    s = sprintf('<%s %s>', strjoin(cellstr(string(size(x))), 'x'), class(x));
end
end
