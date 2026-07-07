function [G_q, G_avg, track_state] = emt_compute_gq_and_average(G_local, K, J_q, params)
% EMT_COMPUTE_GQ_AND_AVERAGE Compute q-resolved Green functions and average.

I3 = eye(size(G_local, 1));
n_q = size(J_q, 3);
n_omega = size(G_local, 3);

[G_embed, track_state] = emt_apply_tracka(G_local, K, params);

G_q = zeros(size(G_local, 1), size(G_local, 2), n_q, n_omega, 'like', G_local);
G_avg = zeros(size(G_local), 'like', G_local);

for iw = 1:n_omega
    G_iw = G_embed(:,:,iw);
    G_rhs = build_rhs(G_iw, I3, params);

    for iq = 1:n_q
        D = I3 + (J_q(:,:,iq) - K(:,:,iw)) * G_rhs;
        G_q(:,:,iq,iw) = emt_stable_left_solve( ...
            D, G_rhs, params.reg_epsilon, params.rcond_tol);
    end

    G_avg(:,:,iw) = mean(G_q(:,:,:,iw), 3);
end

end

function G_rhs = build_rhs(G_iw, I3, params)
    if strcmpi(params.frequency_domain, 'real')
        G_rhs = G_iw + 1i * params.eta * I3;
    elseif strcmpi(params.frequency_domain, 'matsubara')
        G_rhs = G_iw;
    else
        error('emt_compute_gq_and_average:badDomain', ...
            'Unknown frequency_domain: %s', params.frequency_domain);
    end
end
