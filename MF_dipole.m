
function [dip, nN, geom] = MF_dipole(q, N, a, tau, geom)
% This function performs a brute force summation of
% the q-dependent dipole coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% N is the number of unit cells that should be summed in each direction.
%
% The q-INDEPENDENT real-space lattice geometry (reciprocal basis, the
% meshgrid/r0, and per sublattice-pair the shifted cutoff r together with the
% 9 tensor-factor columns) is built once and returned in the struct `geom`.
% Pass it back in as the 5th argument to skip the rebuild across a q-sweep:
%   [~,~,geom] = MF_dipole(q1,N,a,tau);
%   dip        = MF_dipole(q2,N,a,tau,geom);
% The 4-argument call MF_dipole(q,N,a,tau) is unchanged and bit-identical.

if nargin < 5 || isempty(geom)
    geom = MF_dipole_geom(N, a, tau);
end

% --- q-dependent evaluation (mirrors the original arithmetic exactly) ---
q = q*geom.b; % Convert q from Miller indicies to reciprocal aangstroms
ntau = geom.ntau;
dip = zeros(3,3,ntau,ntau);
for nt = 1:ntau
    for mt = 1:nt
        r  = geom.r{nt,mt};
        Tf = geom.Tf{nt,mt};
        exp_qr = exp(-1i*q*r');
        for nn = 1:3
            for mm = 1:3
                dip(nn,mm,nt,mt) = -exp_qr*Tf(:,nn,mm);
            end
        end
        dip(:,:,mt,nt) = conj(dip(:,:,nt,mt)); % J_ij = J_ji*
    end
end
dip = squeeze(dip);
nN = geom.nN;
return
end

function geom = MF_dipole_geom(N, a, tau)
% Build the q-independent lattice geometry used by MF_dipole.
tau = tau*a; % Convert tau to a
vol = sum(a(1,:).*cross(a(2,:),a(3,:))); % Unit cell volume
% Reciprocal lattice unit vectors
b = [2*pi*cross(a(2,:),a(3,:))/vol
     2*pi*cross(a(3,:),a(1,:))/vol
     2*pi*cross(a(1,:),a(2,:))/vol];
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
ntau = size(tau,1);
nN = ntau;
r_cell  = cell(ntau,ntau);
Tf_cell = cell(ntau,ntau);
for nt = 1:ntau
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
        Tf = zeros(numel(r2),3,3);
        for nn = 1:3
            for mm = 1:3
                Tf(:,nn,mm) = 3*r(:,nn).*r(:,mm)./r5 - eq(nn,mm)./r3;
            end
        end
        r_cell{nt,mt}  = r;
        Tf_cell{nt,mt} = Tf;
    end
end
geom.b    = b;
geom.ntau = ntau;
geom.r    = r_cell;
geom.Tf   = Tf_cell;
geom.nN   = nN;
end
