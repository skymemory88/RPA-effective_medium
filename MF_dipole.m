
function [dip, nN] = MF_dipole(q, N, a, tau)
% This function performs a brute force summation of
% the q-dependent dipole coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% N is the number of unit cells that should be summed in each direction.

tau = tau*a; % Convert tau to a
vol = sum(a(1,:).*cross(a(2,:),a(3,:))); % Unit cell volume
% Reciprocal lattice unit vectors
b = [2*pi*cross(a(2,:),a(3,:))/vol
     2*pi*cross(a(3,:),a(1,:))/vol
     2*pi*cross(a(1,:),a(2,:))/vol];
q = q*b; % Convert q from Miller indicies to reciprocal aangstroms
% % Length of q
% qq=sqrt(sum(q.*q));

% For N moments in the unit cell, there will be N 
% coupling parameters J_ij. Many of these will be symmetry related
% (e.g. J_ij=J_ji), so we just calculate J_1j. 
% The result is a (3x3xN) matrix, where the first two dimensions
% are the x,y and z components. The last dimension holds the coupling
% between different ions in the unit cell.

[x,y,z] = meshgrid(-N:N,-N:N,-N:N);
hkl = [x(:) y(:) z(:)];
% hkl = [z(:) x(:) y(:)]; % use z x y to get nicer order in the list - but of no importance

r0 = hkl*a;
nN = size(tau,1);
dip = zeros(3,3,size(tau,1),size(tau,1));
for nt = 1:size(tau,1)
    for mt = 1:nt
        r = r0;
        r(:,1) = r(:,1) - tau(nt,1) + tau(mt,1);
        r(:,2) = r(:,2) - tau(nt,2) + tau(mt,2);
        r(:,3) = r(:,3) - tau(nt,3) + tau(mt,3);
        r2 = vecnorm(r,2,2).^2;
%         cut_off = find(rr<0.01); % singularities
        cut_off = find(r2>(N * min(vecnorm(a,2,2)))^2 | r2<0.01); % define interaction range
        r(cut_off,:) = []; % clear Spins outside of the sphere
        r2(cut_off) = []; % and the central spin itself
        nN(nt) = length(r2);
        r3 = r2 .* vecnorm(r,2,2);
        r5 = r2 .* r3;
        exp_qr = exp(-1i*q*r');
        for nn = 1:3
            for mm = 1:3
                dip(nn,mm,nt,mt) = -exp_qr*(3*r(:,nn).*r(:,mm)./r5 - eq(nn,mm)./r3);
            end
        end
        dip(:,:,mt,nt) = conj(dip(:,:,nt,mt)); % J_ij = J_ji*
    end
end
dip = squeeze(dip);
return
end