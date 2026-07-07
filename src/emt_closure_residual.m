function max_res = emt_closure_residual(K, G_local, J_q, params)
% EMT_CLOSURE_RESIDUAL Compute closure residual max over frequencies.

I3 = eye(size(G_local, 1));
n_q = size(J_q, 3);
n_omega = size(G_local, 3);
max_res = 0;

for iw = 1:n_omega
    G_iw = G_local(:,:,iw);
    G_rhs = G_iw;
    if strcmpi(params.frequency_domain, 'real')
        G_rhs = G_rhs + 1i * params.eta * I3;
    end

    sum_jg = zeros(size(G_iw), 'like', G_iw);
    for iq = 1:n_q
        D = I3 + (J_q(:,:,iq) - K(:,:,iw)) * G_rhs;
        G_q = emt_stable_left_solve(D, G_rhs, params.reg_epsilon, params.rcond_tol);
        sum_jg = sum_jg + J_q(:,:,iq) * G_q;
    end
    sum_jg = sum_jg / n_q;

    denom = norm(K(:,:,iw) * G_iw, 'fro') + eps;
    res = norm(sum_jg - K(:,:,iw) * G_iw, 'fro') / denom;
    max_res = max(max_res, res);
end

end
