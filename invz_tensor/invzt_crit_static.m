function [crit, cmass, arank] = invzt_crit_static(ctil0, JtGamma, rank_tol)
%INVZT_CRIT_STATIC Static Gamma-point criticality of a renormalized local 3x3.
%   [crit, cmass, arank] = INVZT_CRIT_STATIC(ctil0, JtGamma, rank_tol) returns
%   the min real eigenvalue of (M+M')/2, M = I - S*JtGamma*S, with S the
%   rank-clipped PSD square root of C12 = kron(eye(4), herm(ctil0)) (Hermitian
%   eigendecomposition, NOT sqrtm). crit shares the zeros of I - C12*JtGamma on
%   the active subspace; crit > 0 in a locally stable phase (PM above Bc, and
%   the converged FM state away from Bc). cmass = clipped eigenvalue mass,
%   arank = active rank. Extracted verbatim from INVZT_SOLVE_POINT local_crit
%   (2026-07-19) so the ordered solver shares one implementation.
if nargin < 3, rank_tol = 1e-12; end
C12 = kron(eye(4), (ctil0 + ctil0')/2);
[U, D] = eig((C12 + C12')/2);
d = real(diag(D));
clip = d < rank_tol;
cmass = sum(abs(d(clip)));
arank = sum(~clip);
d(clip) = 0;
S = U * diag(sqrt(max(d, 0))) * U';
M = eye(size(S,1)) - S*JtGamma*S;
crit = min(real(eig((M + M')/2)));
end
