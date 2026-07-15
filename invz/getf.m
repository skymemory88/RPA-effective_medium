function v = getf(s, f, d)
%GETF Value of struct field s.f, or default d when the field is absent.
% Shared option-struct accessor for the invz drivers (replaces per-file copies).
if isfield(s, f), v = s.(f); else, v = d; end
end
