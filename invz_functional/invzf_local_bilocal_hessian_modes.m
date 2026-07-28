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
% At fixed source h, the partial bilocal Legendre curvature is instead
% beta^(-1)*inv(K).  A common functional that varies h and m independently
% must use this fixed-source block; the Schur inverse applies only after
% eliminating the one-point coordinate at fixed m.
%
% The oracle also grades the finite-cutoff three-site mixed-chain
% coefficient.  The representative weights are w_0=1 and w_n=2 for n>0.
% The Gaussian trace has local curvature inv(Kg)/beta, with
%
%   Kg=diag((1+delta_n0)*C2(n)^2).
%
% Replacing that curvature by inv(K)/beta and re-stationarizing the leading
% central-site return C2(n)^3*(a^2+b^2) must change the mixed coefficient by
%
%   -sum_nm w_n*w_m*C4_nm*C2_n*C2_m/(4*beta^2).
%
% This is the exact finite-cutoff connected-C4 residual.  It fixes the
% trace/local double-counting normalization without constructing a lattice
% solver.
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
frequencyWeights = 2*ones(nmodes,1);
frequencyWeights(1) = 1;

out = struct( ...
    'schema','invzf_local_bilocal_hessian_modes/v2', ...
    'status','invalid_connected_c4_symmetry', ...
    'labels',labels,'C2',C2,'C3',C3,'C4',C4,'beta',beta, ...
    'c4_symmetry_residual',symmetryResidual, ...
    'c4_symmetry_tolerance',symmetryTolerance, ...
    'frequency_weights',frequencyWeights, ...
    'pairing_diagonal',pairingDiagonal, ...
    'one_bilinear_response',nan(nmodes,1), ...
    'bilinear_covariance',nan(nmodes), ...
    'bilinear_covariance_eigenvalues',nan(nmodes,1), ...
    'normalized_response',nan(nmodes+1), ...
    'normalized_response_eigenvalues',nan(nmodes+1,1), ...
    'fixed_m_schur',nan(nmodes), ...
    'fixed_m_schur_eigenvalues',nan(nmodes,1), ...
    'dimensionless_hessian',nan(nmodes+1), ...
    'energy_hessian_m_D',nan(nmodes+1), ...
    'energy_bilocal_hessian_gaussian',nan(nmodes), ...
    'energy_bilocal_hessian_fixed_source',nan(nmodes), ...
    'energy_bilocal_hessian_fixed_m',nan(nmodes), ...
    'chain_ring_coefficient',NaN, ...
    'chain_exact_c4_residual',NaN, ...
    'chain_exact_coefficient',NaN, ...
    'chain_fixed_source_stationary_residual',NaN, ...
    'chain_fixed_source_stationary_coefficient',NaN, ...
    'chain_fixed_source_identity_error',NaN, ...
    'chain_fixed_m_stationary_residual',NaN, ...
    'chain_fixed_m_stationary_coefficient',NaN, ...
    'chain_fixed_m_vs_exact_error',NaN);

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
bilinearEigenvalues = eig(K);
schurEigenvalues = eig(schur);

out.status = 'invalid_local_bilinear_covariance';
out.C4 = C4s;
out.one_bilinear_response = B;
out.bilinear_covariance = K;
out.bilinear_covariance_eigenvalues = bilinearEigenvalues;
out.normalized_response = response;
out.normalized_response_eigenvalues = responseEigenvalues;
out.fixed_m_schur = schur;
out.fixed_m_schur_eigenvalues = schurEigenvalues;
if ~(all(isfinite(responseEigenvalues)) && ...
        all(isfinite(bilinearEigenvalues)) && ...
        all(isfinite(schurEigenvalues)) && ...
        min(responseEigenvalues) > 0 && ...
        min(bilinearEigenvalues) > 0 && min(schurEigenvalues) > 0)
    return
end

dimensionlessHessian = response\eye(nmodes+1);
dimensionlessHessian = ...
    (dimensionlessHessian+dimensionlessHessian.')/2;
scaleTransform = diag([sqrt(beta);ones(nmodes,1)]);
energyHessian = ...
    (scaleTransform.'*dimensionlessHessian*scaleTransform)/beta;
energyHessian = (energyHessian+energyHessian.')/2;
bilocalHessianFixedSource = (K\eye(nmodes))/beta;
bilocalHessianFixedSource = ...
    (bilocalHessianFixedSource+bilocalHessianFixedSource.')/2;
bilocalHessianFixedMean = (schur\eye(nmodes))/beta;
bilocalHessianFixedMean = ...
    (bilocalHessianFixedMean+bilocalHessianFixedMean.')/2;
bilocalHessianGaussian = diag(1./(beta*pairingDiagonal));

% At leading order the central local return of a three-site open chain is
% v*(a^2+b^2), v_n=C2(n)^3.  The following is the change in the stationary
% mixed coefficient when the Gaussian curvature is replaced by the exact
% local curvature.  Keeping the Gaussian term explicit makes the
% double-counting cancellation checkable.
leadingReturn = C2.^3;
gaussianForce = bilocalHessianGaussian*leadingReturn;
gaussianStationaryTerm = ...
    leadingReturn.'*bilocalHessianGaussian*leadingReturn;
fixedSourceResidual = ...
    -gaussianForce.'*(bilocalHessianFixedSource\gaussianForce) ...
    +gaussianStationaryTerm;
fixedMeanResidual = ...
    -gaussianForce.'*(bilocalHessianFixedMean\gaussianForce) ...
    +gaussianStationaryTerm;
weightedC2 = frequencyWeights.*C2;
exactResidual = ...
    -(weightedC2.'*C4s*weightedC2)/(4*beta^2);
ringCoefficient = ...
    -sum(frequencyWeights.*C2.^4)/(2*beta);

out.status = 'ok';
out.dimensionless_hessian = dimensionlessHessian;
out.energy_hessian_m_D = energyHessian;
out.energy_bilocal_hessian_gaussian = bilocalHessianGaussian;
out.energy_bilocal_hessian_fixed_source = bilocalHessianFixedSource;
out.energy_bilocal_hessian_fixed_m = bilocalHessianFixedMean;
out.chain_ring_coefficient = ringCoefficient;
out.chain_exact_c4_residual = exactResidual;
out.chain_exact_coefficient = ringCoefficient+exactResidual;
out.chain_fixed_source_stationary_residual = fixedSourceResidual;
out.chain_fixed_source_stationary_coefficient = ...
    ringCoefficient+fixedSourceResidual;
out.chain_fixed_source_identity_error = fixedSourceResidual-exactResidual;
out.chain_fixed_m_stationary_residual = fixedMeanResidual;
out.chain_fixed_m_stationary_coefficient = ...
    ringCoefficient+fixedMeanResidual;
out.chain_fixed_m_vs_exact_error = fixedMeanResidual-exactResidual;
end
