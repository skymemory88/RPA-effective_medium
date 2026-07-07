function K0 = emt_initial_k_from_seed(G0, J_q, params)
% EMT_INITIAL_K_FROM_SEED Build initial K guess from one closure pass.

if strcmpi(params.init_k, 'zero')
    K0 = zeros(size(G0), 'like', G0);
    return;
end

I3 = eye(size(G0, 1));
n_omega = size(G0, 3);
n_q = size(J_q, 3);
K0 = zeros(size(G0), 'like', G0);

for iw = 1:n_omega
    G_iw = G0(:,:,iw);
    G_rhs = G_iw;
    if strcmpi(params.frequency_domain, 'real')
        G_rhs = G_rhs + 1i * params.eta * I3;
    end

    sum_jg = zeros(size(G_iw), 'like', G_iw);
    for iq = 1:n_q
        D = I3 + J_q(:,:,iq) * G_rhs;
        G_q = emt_stable_left_solve(D, G_rhs, params.reg_epsilon, params.rcond_tol);
        sum_jg = sum_jg + J_q(:,:,iq) * G_q;
    end
    sum_jg = sum_jg / n_q;

    rc = rcond(G_iw);
    if isfinite(rc) && rc > params.rcond_tol
        K_tmp = sum_jg / G_iw;
        if all(isfinite(K_tmp(:)))
            K0(:,:,iw) = 0.5 * (K_tmp + K_tmp');
        end
    end
end

end
