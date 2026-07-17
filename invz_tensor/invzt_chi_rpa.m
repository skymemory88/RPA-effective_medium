function X = invzt_chi_rpa(chi0, Jt)
%INVZT_CHI_RPA Page-wise 12x12 tensor RPA propagator at one frequency.
%   X = INVZT_CHI_RPA(chi0, Jt) takes the single-ion Cartesian susceptibility
%   chi0 [3,3] at ONE frequency (real or complex) and a [12,12,nq] coupling
%   tensor Jt (see INVZT_JQ_TENSOR; index i = 3*(s-1)+mu, s = 1..4 sublattice,
%   mu = 1(a),2(b),3(c)) and returns the page-wise lattice RPA response
%   X [12,12,nq]:
%       X0        = kron(eye(4), chi0)              [12,12] (uniform local response)
%       X(:,:,iq) = (eye(12) - X0*Jt(:,:,iq)) \ X0
%   i.e. the RPA resummation X = chi0*(I - J*chi0)^-1 rewritten as a left
%   division (I - X0*J)\X0 (equal by the resolvent identity (I-AB)^-1*A =
%   A*(I-BA)^-1) -- solved via mldivide rather than forming an explicit
%   inverse.
%
%   chi0 [3,3]      : single-ion susceptibility tensor at one frequency.
%   Jt   [12,12,nq] : coupling tensor, one page per q (or a single [12,12]
%                     matrix, treated as nq = 1).
%
%   See also INVZT_JQ_TENSOR, INVZT_GCC_LATTICE.
X0 = kron(eye(4), chi0);
nq = size(Jt, 3);
X = zeros(12, 12, nq);
for iq = 1:nq
    X(:,:,iq) = (eye(12) - X0*Jt(:,:,iq)) \ X0;
end
end
