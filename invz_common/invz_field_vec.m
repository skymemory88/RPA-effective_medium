function B = invz_field_vec(B)
%INVZ_FIELD_VEC Normalize a field argument to a 1x3 row [Bx By Bz] in Tesla.
% Scalar b -> [b 0 0] (the historical transverse-along-a convention); any real,
% finite, numeric 3-element vector (row or column) -> 1x3 row, values unchanged.
% Anything else (NaN/Inf, complex, empty, wrong length, non-numeric) errors with
% identifier 'invz:fieldVec'. Idempotent, so boundaries may normalize freely.
if ~isnumeric(B) || ~isreal(B) || ~all(isfinite(B(:)))
    error('invz:fieldVec', 'Field must be real, finite, and numeric.');
end
if isscalar(B)
    B = [B 0 0];
elseif numel(B) == 3
    B = reshape(B, 1, 3);
else
    error('invz:fieldVec', ...
        'Field must be a scalar or a 3-element vector; got %d element(s).', numel(B));
end
end
