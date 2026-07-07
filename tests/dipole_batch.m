function dip_all = dipole_batch(qvec, N, a, tau, blockSize)
% dipole_batch Vectorized dipole tensor over all q points (batched in r).
%
% dip_all: [3,3,Nions,Nions,Nq], same convention as MF_dipole per-q outputs stacked.
% Implements the same lattice sum as MF_dipole.m but evaluates all q points
% together using block matrix multiplies to reduce overhead.
%
% Inputs
%   qvec      : [Nq x 3] q-points in Miller indices [h k l]
%   N         : scalar, dipole summation range (MF_dipole's N)
%   a         : [3x3] direct lattice vectors (rows)
%   tau       : [Nions x 3] fractional basis positions
%   blockSize : optional, number of r-vectors per block (default 8192)
%
% Notes
% - Matches MF_dipole cutoff: keep r with r2 <= (N*min(||a_i||))^2 and r2 >= 0.01
% - Uses qcart = qvec*b with reciprocal basis b built from a, same as MF_dipole.

if nargin < 5 || isempty(blockSize)
    blockSize = 8192; % tune based on available memory
end

% Reciprocal basis and q in Cartesian
vol = sum(a(1,:).*cross(a(2,:),a(3,:)));
b = [2*pi*cross(a(2,:),a(3,:))/vol;
     2*pi*cross(a(3,:),a(1,:))/vol;
     2*pi*cross(a(1,:),a(2,:))/vol];
qcart = qvec * b; % [Nq x 3]

% Convert tau to Cartesian coordinates
tau_cart = tau * a;

% Build base lattice vectors r0 = hkl * a
[x,y,z] = meshgrid(-N:N, -N:N, -N:N);
hkl = [x(:) y(:) z(:)];
r0 = hkl * a; % [Pr x 3]

% Sphere cutoff radius (MF_dipole)
R = N * min(vecnorm(a,2,2));
R2 = R^2;

Nq = size(qcart,1);
Nions = size(tau,1);
dip_all = zeros(3,3,Nions,Nions,Nq);

for nt = 1:Nions
    for mt = 1:nt
        % Shift r by tau differences (Cartesian)
        r = r0;
        r(:,1) = r(:,1) - tau_cart(nt,1) + tau_cart(mt,1);
        r(:,2) = r(:,2) - tau_cart(nt,2) + tau_cart(mt,2);
        r(:,3) = r(:,3) - tau_cart(nt,3) + tau_cart(mt,3);

        % Apply MF_dipole cutoff
        r2 = sum(r.^2, 2);
        keep = (r2 <= R2) & (r2 >= 0.01);
        r = r(keep, :);
        r2 = r2(keep);

        if isempty(r)
            continue;
        end

        rn = sqrt(r2);
        r3 = r2 .* rn;
        r5 = r2 .* r3;

        % Accumulators for each tensor component across q
        S = zeros(3,3,Nq);

        % Process r in manageable blocks to bound memory
        P = size(r,1);
        for s = 1:blockSize:P
            e = min(P, s+blockSize-1);
            rb = r(s:e, :);           % [B x 3]
            r3b = r3(s:e, :);         % [B x 1]
            r5b = r5(s:e, :);         % [B x 1]

            % Phase matrix for this block: [Nq x B]
            E = exp(-1i * (qcart * rb.'));

            % For each tensor component (n,m), accumulate sum(E * w_nm)
            for n = 1:3
                for m = 1:3
                    w = 3.0 * (rb(:,n) .* rb(:,m)) ./ r5b - (n==m) ./ r3b;
                    Snm = E * w; % [Nq x 1]
                    S(n,m,:) = S(n,m,:) + reshape(Snm, 1,1,[]);
                end
            end
        end

        % Store (with the negative sign per MF_dipole) and symmetrize
        for n = 1:3
            for m = 1:3
                dip_all(n,m,nt,mt,:) = reshape(-S(n,m,:), 1,1,1,1,[]);
                dip_all(n,m,mt,nt,:) = conj(reshape(-S(n,m,:), 1,1,1,1,[]));
            end
        end
    end
end

end

