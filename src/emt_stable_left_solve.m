function X = emt_stable_left_solve(A, B, reg_epsilon, rcond_tol)
% EMT_STABLE_LEFT_SOLVE Robust linear solve with adaptive diagonal regularization.

rc = rcond(A);
if ~isfinite(rc) || rc < rcond_tol
    A = A + reg_epsilon * eye(size(A));
end

X = A \ B;
if any(~isfinite(X(:)))
    A = A + 10 * reg_epsilon * eye(size(A));
    X = A \ B;
end

if any(~isfinite(X(:)))
    X = zeros(size(B), 'like', B);
end

end
