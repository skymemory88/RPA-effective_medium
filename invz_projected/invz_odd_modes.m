function [Jnu, Juni] = invz_odd_modes(Vcc, dJ)
%INVZ_ODD_MODES Sorted cc mode spectrum of the (ODD-corrected) 4x4 sublattice blocks.
% Jnu(iq,:) = sort(real(eig(M))), M = Hermitized Vcc(:,:,iq) [+ dJ(:,:,iq)].
% dJ may be [] (no correction). Second output Juni(iq) = real(v'*M*v), v = ones(4,1)/2
% (uniform FM-mode projection), computed only when requested.
% VALUES-ONLY by design (eig jobz 'N'): the eigenvector sites in
% invz_solve_point's optional retarded path deliberately uses
% separate eig(M,'vector') calls -- do not fold them in here (last-bit identity).
nq   = size(Vcc, 3);
Jnu  = zeros(nq, 4);
addJ = ~isempty(dJ);
uni  = nargout > 1;
if uni
    Juni = zeros(nq, 1);
    v = ones(4,1)/2;                 % uniform (all-sublattices-in-phase) FM mode
end
for iq = 1:nq
    if addJ
        M = Vcc(:,:,iq) + dJ(:,:,iq);
    else
        M = Vcc(:,:,iq);
    end
    M = (M + M')/2;                  % both terms Hermitian; cleans rounding only
    Jnu(iq,:) = sort(real(eig(M))).';
    if uni
        Juni(iq) = real(v.'*M*v);    % uniform FM-mode coupling (physical dispersion)
    end
end
end
