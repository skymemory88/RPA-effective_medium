function out = invzf_local_bilocal_hessian_modes(labels,C2,C3,C4,beta)
%INVZF_LOCAL_BILOCAL_HESSIAN_MODES Exact local mode-pair Hessian oracle.
%
% LABELS are distinct nonnegative Matsubara representatives and must begin
% at zero.  With
%
%   Y_n=X_n/sqrt(beta),  A_n=Y_n*Y_-n,
%   mu=<Y_0>,             D_n=<A_n>-delta_n0*mu^2,
%
% introduce dimensionless representative-mode sources
%
%   W(j,q)=log <exp(j*Y_0+sum_n q_n*A_n)>.
%
% The double Legendre transform has canonical one-point source
% a=j+2*q_0*mu.  At q=0, its exact response block is
%
%   d(mu,D)/d(a,q) = [ C2(0), B'
%                      B,     K ],
%
%   B_n = C3(0,n,-n)/sqrt(beta),
%
%   K_nm = C4(n,-n,m,-m)/beta
%          +delta_nm*(1+delta_n0)*C2(n)^2.
%
% C2, C3, and C4 are connected local cumulants, not 1PI vertices.  The
% extra zero-mode pairing is exact: A_0=Y_0^2 has two crossed Gaussian
% pairings, whereas A_n with n>0 has one within the nonnegative
% representative basis.
%
% For the energy functional Gamma/beta in variables (m,D), mu=sqrt(beta)m.
% The full Hessian is beta^(-1)*L'*inv(R)*L with
% L=diag(sqrt(beta),1,...,1).  At fixed m, the bilocal Hessian is
%
%   beta^(-1)*inv(K-B*B'/C2(0)).
%
% This is a local finite-mode oracle.  It does not supply a Matsubara-tail
% completion, lattice trace, double counting, or stationary skeleton.

validateattributes(labels,{'numeric'}, ...
    {'real','vector','finite','integer','nonnegative'},mfilename,'labels');
validateattributes(C2,{'numeric'}, ...
    {'real','vector','finite','positive'},mfilename,'C2');
validateattributes(C3,{'numeric'}, ...
    {'real','vector','finite'},mfilename,'C3');
validateattributes(C4,{'numeric'}, ...
    {'real','2d','finite'},mfilename,'C4');
validateattributes(beta,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'beta');

labels = labels(:);
C2 = C2(:);
C3 = C3(:);
nmodes = numel(labels);
if nmodes == 0 || labels(1) ~= 0 || any(diff(labels) <= 0)
    error('invzf:bilocalModeLabels', ...
        'labels must be distinct, increasing, nonnegative, and begin at zero.');
end
if numel(C2) ~= nmodes || numel(C3) ~= nmodes || ...
        ~isequal(size(C4),[nmodes,nmodes])
    error('invzf:bilocalModeShape', ...
        'C2, C3, and C4 must match the number of representative labels.');
end

symmetryResidual = norm(C4-C4.',inf);
symmetryScale = max(1,norm(C4,inf));
symmetryTolerance = 4096*eps(symmetryScale);
pairingDiagonal = C2.^2;
pairingDiagonal(1) = 2*C2(1)^2;

out = struct( ...
    'schema','invzf_local_bilocal_hessian_modes/v1', ...
    'status','invalid_connected_c4_symmetry', ...
    'labels',labels,'C2',C2,'C3',C3,'C4',C4,'beta',beta, ...
    'c4_symmetry_residual',symmetryResidual, ...
    'c4_symmetry_tolerance',symmetryTolerance, ...
    'pairing_diagonal',pairingDiagonal, ...
    'one_bilinear_response',nan(nmodes,1), ...
    'bilinear_covariance',nan(nmodes), ...
    'normalized_response',nan(nmodes+1), ...
    'normalized_response_eigenvalues',nan(nmodes+1,1), ...
    'fixed_m_schur',nan(nmodes), ...
    'fixed_m_schur_eigenvalues',nan(nmodes,1), ...
    'dimensionless_hessian',nan(nmodes+1), ...
    'energy_hessian_m_D',nan(nmodes+1), ...
    'energy_bilocal_hessian_fixed_m',nan(nmodes));

if symmetryResidual > symmetryTolerance
    return
end

C4s = (C4+C4.')/2;
B = C3/sqrt(beta);
K = C4s/beta+diag(pairingDiagonal);
response = [C2(1),B.';B,K];
response = (response+response.')/2;
schur = K-(B*B.')/C2(1);
schur = (schur+schur.')/2;
responseEigenvalues = eig(response);
schurEigenvalues = eig(schur);

out.status = 'invalid_local_bilinear_covariance';
out.C4 = C4s;
out.one_bilinear_response = B;
out.bilinear_covariance = K;
out.normalized_response = response;
out.normalized_response_eigenvalues = responseEigenvalues;
out.fixed_m_schur = schur;
out.fixed_m_schur_eigenvalues = schurEigenvalues;
if ~(all(isfinite(responseEigenvalues)) && ...
        all(isfinite(schurEigenvalues)) && ...
        min(responseEigenvalues) > 0 && min(schurEigenvalues) > 0)
    return
end

dimensionlessHessian = response\eye(nmodes+1);
dimensionlessHessian = ...
    (dimensionlessHessian+dimensionlessHessian.')/2;
scaleTransform = diag([sqrt(beta);ones(nmodes,1)]);
energyHessian = ...
    (scaleTransform.'*dimensionlessHessian*scaleTransform)/beta;
energyHessian = (energyHessian+energyHessian.')/2;
bilocalHessian = (schur\eye(nmodes))/beta;
bilocalHessian = (bilocalHessian+bilocalHessian.')/2;

out.status = 'ok';
out.dimensionless_hessian = dimensionlessHessian;
out.energy_hessian_m_D = energyHessian;
out.energy_bilocal_hessian_fixed_m = bilocalHessian;
end
