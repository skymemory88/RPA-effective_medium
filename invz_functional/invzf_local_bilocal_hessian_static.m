function out = invzf_local_bilocal_hessian_static(C2,C3,C4,beta)
%INVZF_LOCAL_BILOCAL_HESSIAN_STATIC Exact Q=0 static local Hessian oracle.
%
% This is an isolated oracle for the local kernel missing from the rejected
% Gaussian varied-covariance ansatz.  Let X0=int_0^beta X(tau) dtau and
% Y=X0/sqrt(beta).  Start from
%
%   W(j,K)=log <exp(j*Y+K*Y^2/2)>,
%
% and write q=K/2, mu=<Y>, and D=<Y^2>-mu^2=C2(0).  For
%
%   Gamma=-W+j*mu+q*(D+mu^2),
%
% the canonical one-form is
%
%   dGamma=(j+2*q*mu)*dmu+q*dD.
%
% Thus a=j+2*q*mu is conjugate to mu and q is conjugate to D.  At q=0,
% the exact response d(mu,D)/d(a,q) is
%
%   R = [ C2,             C3/sqrt(beta)
%         C3/sqrt(beta),  2*C2^2+C4/beta ].
%
% C2, C3, and C4 are connected zero-mode cumulants; in particular,
% C4=C4(0,0,0,0)=d^2 C2(0)/dh^2.  C4 is not a 1PI gamma4 vertex.
% The C2^2 term is the pair-disconnected part of the covariance of two
% bilinears; omitting it or replacing the block inverse by a 1PI gamma4 is
% the structural error exposed by the three-site mixed-bond oracle.
% For the energy functional Gamma/beta in variables (m,D),
%
%   H = beta^(-1)*L'*inv(R)*L,  L=diag(sqrt(beta),1).
%
% In particular, the exact bilocal curvature at fixed m is
%
%   H_DD = 1/(C4+2*beta*C2^2-C3^2/C2).
%
% At fixed local source h instead, the curvature is
%
%   H_DD|h = 1/(C4+2*beta*C2^2).
%
% The distinction is material away from Z2 symmetry.  A common functional
% that keeps h and m as independent variables must use the fixed-h
% curvature; the Schur curvature applies only after eliminating the
% one-point coordinate at fixed m.
%
% OUT reports both curvatures and their quadratic coefficients.  It also
% reports the Gaussian coefficient 1/(4*beta*C2^2) and the corresponding
% corrections.  This does not derive the
% finite-frequency/finite-Q bilocal kernel, lattice trace, double counting,
% or stationary functional, and performs no production dispatch.

validateattributes(C2,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'C2');
validateattributes(C3,{'numeric'}, ...
    {'real','scalar','finite'},mfilename,'C3');
validateattributes(C4,{'numeric'}, ...
    {'real','scalar','finite'},mfilename,'C4');
validateattributes(beta,{'numeric'}, ...
    {'real','scalar','finite','positive'},mfilename,'beta');

response = [C2,C3/sqrt(beta); ...
            C3/sqrt(beta),2*C2^2+C4/beta];
response = (response+response.')/2;
responseEigenvalues = eig(response);
schur = C4+2*beta*C2^2-C3^2/C2;
fixedSourceDenominator = C4+2*beta*C2^2;

out = struct( ...
    'schema','invzf_local_bilocal_hessian_static/v2', ...
    'status','invalid_local_bilinear_covariance', ...
    'C2',C2,'C3',C3,'C4',C4,'beta',beta, ...
    'normalized_response',response, ...
    'normalized_response_eigenvalues',responseEigenvalues, ...
    'schur_fixed_m',schur, ...
    'denominator_fixed_source',fixedSourceDenominator, ...
    'dimensionless_hessian',nan(2), ...
    'energy_hessian_m_D',nan(2), ...
    'fixed_m_D_curvature',NaN, ...
    'fixed_m_D_quadratic_coefficient',NaN, ...
    'fixed_source_D_curvature',NaN, ...
    'fixed_source_D_quadratic_coefficient',NaN, ...
    'gaussian_D_curvature',1/(2*beta*C2^2), ...
    'gaussian_D_quadratic_coefficient',1/(4*beta*C2^2), ...
    'fixed_m_bilocal_quadratic_correction',NaN, ...
    'fixed_source_bilocal_quadratic_correction',NaN, ...
    'static_ring_coefficient',NaN, ...
    'static_exact_c4_coefficient',NaN, ...
    'fixed_source_stationary_coefficient',NaN, ...
    'fixed_source_identity_error',NaN, ...
    'fixed_m_stationary_coefficient',NaN, ...
    'fixed_m_vs_exact_error',NaN, ...
    'zero_source_static_ring_coefficient',NaN, ...
    'zero_source_static_exact_coefficient',NaN, ...
    'zero_source_stationary_coefficient',NaN, ...
    'zero_source_identity_error',NaN);

if ~(all(isfinite(response),'all') && all(isfinite(responseEigenvalues)) && ...
        min(responseEigenvalues) > 0 && isfinite(schur) && schur > 0 && ...
        isfinite(fixedSourceDenominator) && fixedSourceDenominator > 0)
    return
end

dimensionlessHessian = response\eye(2);
dimensionlessHessian = (dimensionlessHessian+dimensionlessHessian.')/2;
scaleTransform = diag([sqrt(beta),1]);
energyHessian = (scaleTransform.'*dimensionlessHessian*scaleTransform)/beta;
energyHessian = (energyHessian+energyHessian.')/2;
fixedCurvature = 1/schur;
fixedQuadratic = fixedCurvature/2;
fixedSourceCurvature = 1/fixedSourceDenominator;
fixedSourceQuadratic = fixedSourceCurvature/2;
gaussianQuadratic = out.gaussian_D_quadratic_coefficient;
fixedMeanCorrection = fixedQuadratic-gaussianQuadratic;
fixedSourceCorrection = fixedSourceQuadratic-gaussianQuadratic;

out.status = 'ok';
out.dimensionless_hessian = dimensionlessHessian;
out.energy_hessian_m_D = energyHessian;
out.fixed_m_D_curvature = fixedCurvature;
out.fixed_m_D_quadratic_coefficient = fixedQuadratic;
out.fixed_source_D_curvature = fixedSourceCurvature;
out.fixed_source_D_quadratic_coefficient = fixedSourceQuadratic;
out.fixed_m_bilocal_quadratic_correction = fixedMeanCorrection;
out.fixed_source_bilocal_quadratic_correction = fixedSourceCorrection;

% This fixed-source restricted-static identity is a necessary weak-coupling
% consistency gate.  The correction replaces, rather than supplements, the
% Gaussian local curvature.  The fixed-m correction differs by a reducible
% C3^2 term away from symmetry and therefore fails this fixed-source gate.
ring = -C2^4/(2*beta);
exact = ring-C4*C2^2/(4*beta^2);
fixedSourceStationaryCorrection = 2*fixedSourceCorrection*C2^6 / ...
    (1+4*beta*fixedSourceCorrection*C2^2);
fixedSourceStationary = ring+fixedSourceStationaryCorrection;
fixedMeanStationaryCorrection = 2*fixedMeanCorrection*C2^6 / ...
    (1+4*beta*fixedMeanCorrection*C2^2);
fixedMeanStationary = ring+fixedMeanStationaryCorrection;
out.static_ring_coefficient = ring;
out.static_exact_c4_coefficient = exact;
out.fixed_source_stationary_coefficient = fixedSourceStationary;
out.fixed_source_identity_error = fixedSourceStationary-exact;
out.fixed_m_stationary_coefficient = fixedMeanStationary;
out.fixed_m_vs_exact_error = fixedMeanStationary-exact;

% Preserve the separately named exact-symmetry gate.  Use exact equality
% deliberately: a merely small C3 at nonzero source must not be silently
% rounded into the Z2 oracle.
if C3 == 0
    out.zero_source_static_ring_coefficient = ring;
    out.zero_source_static_exact_coefficient = exact;
    out.zero_source_stationary_coefficient = fixedSourceStationary;
    out.zero_source_identity_error = fixedSourceStationary-exact;
end
end
