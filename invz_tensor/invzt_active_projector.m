function P = invzt_active_projector(X, rtol)
%INVZT_ACTIVE_PROJECTOR Frequency-consistent union-of-ranges active subspace.
%   P = INVZT_ACTIVE_PROJECTOR(X, rtol) returns an orthonormal basis P [3 x r]
%   for the FREQUENCY-CONSISTENT active subspace of the [3,3,nz] tensor stack X
%   (one Cartesian 3x3 block per frequency): the UNION over frequencies of the
%   per-frequency active channels.
%
%   G3 = sum_k Ck^H*Ck (Ck = Hermitized X(:,:,k)) is a 3x3 PSD Gram whose null
%   space is exactly the INTERSECTION of the per-frequency null spaces, so
%   range(G3) is the UNION of the active subspaces. eig gives an orthonormal
%   basis; the columns are rank-revealed on the magnitude scale sqrt(max(eig,0))
%   relative to the largest active channel, at the relative tolerance rtol.
%
%   Using ONE common P (not a per-frequency rank reveal) keeps the downstream
%   subspace solves continuous across frequencies when a channel skims the rank
%   threshold. Shared by INVZT_EMT_MATRIX (on ctil) and INVZT_SIGMA_TENSOR (on
%   G0) -- see either for the physics context; each keeps its own r = size(P,2).
%
%   See also INVZT_EMT_MATRIX, INVZT_SIGMA_TENSOR.
nz = size(X, 3);
G3 = zeros(3);
for k = 1:nz
    Ck = (X(:,:,k) + X(:,:,k)')/2;              % Hermitized channel content
    G3 = G3 + Ck'*Ck;
end
G3 = (G3 + G3')/2;
[V, dv] = eig(G3, 'vector');
dv = real(dv);
mag = sqrt(max(dv, 0));                          % per-eigvec magnitude scale
smax = max(mag);  if smax <= 0, smax = 1; end    % all-zero guard (P empty -> caller 0)
active = mag > rtol*smax;
P = V(:, active);
end
