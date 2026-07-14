function sig = invz_sigma_ordered(tl, lam, K, g, beta)
%INVZ_SIGMA_ORDERED Ordered-phase (m~=0) Jensen self-energy, HTML eqs 37-38 (J 2.26-2.27).
%
%   Sigma(iwn) = (alpha - alpha_m) + [ gamma(iwn) - (2 m^2 / M^2) gamma(0) ] g(iwn)   (eq 37)
%
% with alpha, gamma the paramagnetic pieces (eqs 27, 23) and the elastic correction
%
%   alpha_m = (m^2/n01^2){ lambda_2 - g0*lambda_1 + (4/g0) lambda_3
%                          - (1-n01^2)(1 + 1/2 beta Delta n01) K(0) g0 }               (eq 38)
%
% Reduces to the paramagnetic invz_sigma exactly at m = 0 (alpha_m = 0, gamma(0) term drops).
% Requires lam = [lambda_1; lambda_2; lambda_3] (invz_lambdas(...,[1 2 3])). The elastic weight
% (1-n01^2) is exponentially small at low T, so alpha_m there is ~ m^2{lambda_2 - g0 lambda_1
% + (4/g0)lambda_3}. K(:) is on the bosonic n>=0 grid so K(1) is the omega_n = 0 slot.
prefM = tl.M2 / tl.n01^2;
K = K(:);  g = g(:);
K0 = K(1);  g0 = tl.g0;

alpha  = prefM * (lam(2) - 0.5*(g0 + beta*(1 - tl.n01^2))*lam(1));      % eq 27
gamma  = prefM * (lam(1) - (1 - tl.n01^2)*K);                          % eq 23, vector over iwn
gamma0 = prefM * (lam(1) - (1 - tl.n01^2)*K0);                         % gamma at omega_n = 0

alpha_m = (tl.m^2 / tl.n01^2) * ( lam(2) - g0*lam(1) + (4/g0)*lam(3) ...
            - (1 - tl.n01^2)*(1 + 0.5*beta*tl.Delta*tl.n01)*K0*g0 );   % eq 38

sig.alpha   = alpha;
sig.alpha_m = alpha_m;
sig.gamma   = gamma;
sig.Sigma   = (alpha - alpha_m) + (gamma - (2*tl.m^2/tl.M2)*gamma0) .* g;
end
