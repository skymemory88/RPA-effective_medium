function ex_all = exchange_batch(qvec, Jex, a, tau)
% exchange_batch Vectorized exchange tensor over all q points.
% ex_all: [3,3,Nions,Nions,Nq], matching exchange() per-q output.

% Convert tau to Cartesian
tau_cart = tau * a;

% Reciprocal lattice and q in Cartesian
vol = sum(a(1,:).*cross(a(2,:),a(3,:)));
b = [2*pi*cross(a(2,:),a(3,:))/vol;
     2*pi*cross(a(3,:),a(1,:))/vol;
     2*pi*cross(a(1,:),a(2,:))/vol];
qcart = qvec * b; % [Nq x 3]

% Lattice offsets (fixed small shell)
[x,y,z] = meshgrid(-1:1,-1:1,-1:1);
hkl = [x(:) y(:) z(:)];

Nions = size(tau,1);
Nq = size(qvec,1);
ex_all = zeros(3,3,Nions,Nions,Nq);

% Precompute squared lattice norm cutoff similar to exchange.m
a_norm = max(vecnorm(a,1));
rr_max = (a_norm).^2;

for nt = 1:Nions
    for mt = 1:nt % avoid double counting
        r0 = hkl * a;           % base lattice vectors
        r = r0;
        r(:,1) = r(:,1) - tau_cart(nt,1) + tau_cart(mt,1);
        r(:,2) = r(:,2) - tau_cart(nt,2) + tau_cart(mt,2);
        r(:,3) = r(:,3) - tau_cart(nt,3) + tau_cart(mt,3);
        rr = sum(r.^2,2);
        mask = ~(rr < 0.01 | rr > rr_max);
        r = r(mask,:);
        rr = rr(mask);
        if isempty(r)
            continue;
        end
        % Nearest-neighbor mask as in original exchange()
        nn_mask = (rr <= 14);

        % Phase matrix: [Nq x Pr]
        E = exp(-1i * (qcart * r.'));
        % Weighted sum per q
        w = double(nn_mask) * Jex; % [Pr x 1]
        S = E * w;                 % [Nq x 1]

        % Fill diagonal spin components; off-diagonals are zero by construction
        for c = 1:3
            ex_all(c,c,nt,mt,:) = reshape(S, 1,1,1,1,[]);
            ex_all(c,c,mt,nt,:) = conj(reshape(S, 1,1,1,1,[]));
        end
    end
end
end

