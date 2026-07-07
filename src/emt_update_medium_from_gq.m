function K_new = emt_update_medium_from_gq(G_q, G_local, J_q, K_old, params)
% EMT_UPDATE_MEDIUM_FROM_GQ Update K from q-resolved Green function closure.

n_omega = size(G_local, 3);
n_q = size(J_q, 3);
K_new = zeros(size(K_old), 'like', K_old);

for iw = 1:n_omega
    sum_jg = zeros(size(G_local, 1), size(G_local, 2), 'like', G_local);
    for iq = 1:n_q
        sum_jg = sum_jg + J_q(:,:,iq) * G_q(:,:,iq,iw);
    end
    sum_jg = sum_jg / n_q;

    G_iw = G_local(:,:,iw);
    rc = rcond(G_iw);
    if isfinite(rc) && rc > params.rcond_tol
        K_tmp = sum_jg / G_iw;
        if all(isfinite(K_tmp(:)))
            K_new(:,:,iw) = 0.5 * (K_tmp + K_tmp');
        else
            K_new(:,:,iw) = K_old(:,:,iw);
        end
    else
        K_new(:,:,iw) = K_old(:,:,iw);
    end
end

end
