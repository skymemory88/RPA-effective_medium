function [Gcc, diag4] = invzt_gcc_lattice(chi0, lat)
%INVZT_GCC_LATTICE Weighted BZ + sublattice average of the site-diagonal cc RPA response.
%   [Gcc, diag4] = INVZT_GCC_LATTICE(chi0, lat) runs INVZT_CHI_RPA once per
%   frequency slice of chi0 [3,3,nz] over the FULL q-grid carried in lat.Jt
%   [12,12,nq], then BZ-averages (quadrature weights lat.w) the site-diagonal cc
%   entries X(3*(s-1)+3, 3*(s-1)+3, :) for each sublattice s = 1..4:
%       diag4(s,k) = sum_q lat.w(q) * real(X_k(3*(s-1)+3, 3*(s-1)+3, q))
%       Gcc(k)     = mean_s diag4(s,k)                 (sublattice average)
%
%   chi0 [3,3,nz] : single-ion susceptibility tensor, one page per frequency.
%   lat           : LOCKED lattice struct from INVZT_JQ_TENSOR (fields Jt, w, ...).
%
%   Gcc   [1,nz] : BZ- and sublattice-averaged site-diagonal cc response.
%   diag4 [4,nz] : per-sublattice BZ-averaged site-diagonal cc response.
%
%   See also INVZT_CHI_RPA, INVZT_JQ_TENSOR.
nz = size(chi0, 3);
w = lat.w(:).';
diag4 = zeros(4, nz);
for k = 1:nz
    X = invzt_chi_rpa(chi0(:,:,k), lat.Jt);
    for s = 1:4
        diag4(s,k) = sum(w .* squeeze(real(X(3*(s-1)+3, 3*(s-1)+3, :))).');
    end
end
Gcc = mean(diag4, 1);
end
